#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 09:13:30 2024

@author: vrath
"""

import os
import sys
from sys import exit as error
import time
from time import process_time
from datetime import datetime
import warnings
import inspect
import copy

import scipy
import scipy.linalg
import scipy.sparse.linalg
import scipy.special
import scipy.fftpack
import scipy.sparse

from numba import njit, prange
import numpy
import functools


import multiprocessing
import joblib

import aesys
import core1d
import util

warnings.simplefilter(action="ignore", category=FutureWarning)


rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")



def run_inv_flightline(data_dir= None,
                        data_file=None,
                        ctrl=None,
                        prior_file=None,
                        result_dir=None,
                        result_strng=None,
                        results_out=True,
                        out=False):
    """
    Wrapper for data_dict parallel set inversion

    Parameters
    ----------

    data_file : strng
        Input data_dict file. The default is None.
    result_file : strng
            site_dict output file. The default is None.
    prior_file : strng, optional
        Prior file. Only required if the read option for priors
        is set. The default is None.

    ctrl : dict
        Contains control variables for this run. The default is None.

    out : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    results_dict : dict
        Contains items which will be stored into resultfile.

    vr, Bloomsday 2024


    ctrl_dict ={
        "system":
            [AEM_system, FwdCall],
        "header":
            [titstrng, ""],
        "inversion":
            numpy.array([RunType, RegFun, Tau0, Tau1, Maxiter, ThreshRMS,
                      LinPars, SetPrior, Delta, RegShift], dtype=object),
        "covar":
            numpy.array([L0, Cm0, L1, Cm1], dtype=object),
        "uncert":
            [Uncert],

        "data":
            numpy.array([DataTrans, data_active, DatErr_add, DatErr_mult, Direction], dtype=object),
        "model":
            numpy.array([ParaTrans, mod_act, mod_apr, mod_var, mod_bnd], dtype=object),
                }


    """
    if data_file is None:
        error("run_inv: No data_dict file given! Exit.")

    name, ext = os.path.splitext(data_file)


    if result_strng is None:
        result_file = name+"_results.npz"
    else:
        result_file = result_dir+name+result_strng+"_results.npz"


    start = time.time()

    AEM_system = ctrl["system"][0]
    _ ,NN, _, _, _, = aesys.get_system_params(System=AEM_system)

    print("\n Reading file " + data_file)
    DataObs, Header, _ = aesys.read_aempy(File=data_dir+data_file,
                                   System=ctrl["system"][0], OutInfo=False)


    fl_name = DataObs[0, 0]
    ctrl["name"] = fl_name

    print("site_dict: ",ctrl.keys())
    ctrl_file = result_file.replace("_results", "_ctrl")
    numpy.savez_compressed(ctrl_file,**ctrl)


    site_x = DataObs[:, 1]
    site_y = DataObs[:, 2]
    site_gps = DataObs[:, 3]
    site_alt = DataObs[:, 4]
    site_dem = DataObs[:, 5]
    dat_obs =  DataObs[:, 6:6+NN[2]]

    [nsite,ndata] = numpy.shape(dat_obs)
    sites = numpy.arange(nsite)


    """
    Setup data-related parameter dict
    """

    data_act = ctrl["data"][1]
    data_err_add = ctrl["data"][2]
    data_err_mult = ctrl["data"][3]
    direction = ctrl["data"][4]

    dat_act = numpy.tile(data_act,(nsite,1))
    dat_err = numpy.zeros_like(dat_obs)
    for ii in sites:
        dat_err[ii, :], _ = set_errors(dat_obs[ii, :],
                                            daterr_add=data_err_add,
                                            daterr_mult=data_err_mult)
        # print("\n",ii)
        # print(dat_obs[ii, :])
        # print(dat_err[ii, :])




    """
    Setup model-related parameter dict
    """
    runtype = ctrl["inversion"][0].lower()
    setprior = ctrl["inversion"][7].lower()
    maxiter =  ctrl["inversion"][4]

    print("XXXXXXXXX", runtype)

    mod_act = ctrl["model"][1].copy()
    mod_apr = ctrl["model"][2].copy()
    mod_var = ctrl["model"][3].copy()
    mod_bnd = ctrl["model"][4].copy()

    # print(mod_apr.shape)
    site_prior = numpy.zeros((nsite,numpy.shape(mod_apr)[0]))

    if "read" in setprior:
        if name.lower() not in prior_file.lower():
            error("Halfspace file name does not match! Exit.")
        site_prior = load_prior(prior_file,
                                       m_ref=mod_apr,
                                       m_apr=mod_apr,
                                       m_act=mod_act)
    if ("set" in setprior.lower()) or ("upd" in setprior.lower()):
        for ii in sites:
                site_prior[ii, :] = mod_apr


    uncert = ctrl["uncert"][0]
    header = ctrl["header"][0]

# This is the main loop over sites in a flight line or within an area:

    """
    Loop over sites
    """
    sequence = numpy.arange(nsite)
    if "rev" in direction.lower():
        sites = sequence[::-1]
    else:
        sites = sequence

    # logsize = (2 + 7*maxiter)
    # site_log = numpy.full((len(sites),logsize), numpy.nan)
    mtmp = numpy.array([])
    icnt = -1
    for ii in sites:
    # for ii in [0, 1, 2]:
        icnt =icnt + 1
        print("\n Invert site #"+str(ii)+"/"+str(len(sites)))

        """
        Setup parameter dict
        """
        data_dict = dict([
            ("d_act", dat_act[ii,:]),
            ("d_obs", dat_obs[ii,:]),
            ("d_err", dat_err[ii,:]),
            ("alt", site_alt[ii])
            ])

        # print(ii)
        # print(dat_obs[ii,:])


        if "read" in setprior:
            mod_apr = site_prior[ii,:]
            mod_ini = mod_apr.copy()

        elif "upd" in setprior:

            if icnt == 0:
                mod_ini = site_prior[ii,:]
                mod_apr = mod_ini.copy()
            else:
                mod_ini = insert_mod(M=mod_apr, m=mtmp[0], m_act=mod_act)
                mod_apr = mod_ini.copy()

        elif "set" in setprior:
            mod_apr = site_prior[ii,:]
            mod_ini = mod_apr.copy()

        model = mod_ini.copy()

        # print("\n",ii)
        # print(mod_ini[ii, :])
        # print(ii, mod_apr)
        # print(mod_var[ii, :])

        model_dict = dict([
            ("m_act", mod_act),
            ("m_apr", mod_apr),
            ("m_var", mod_var),
            ("m_bnd", mod_bnd),
            ("m_ini", mod_ini)
            ])

        # print(ii, mod_apr)

        print("run_inv_flightline: runtype "+runtype+" started.")


        if ("tikh" in runtype) and ("opt" in runtype):
            site_dict = \
                run_tikh_opt(Ctrl=ctrl, Model=model_dict, Data=data_dict,
                                          OutInfo=out)

        elif ("tikh" in runtype) and ("occ" in runtype):
            site_dict = \
                run_tikh_occ(Ctrl=ctrl, Model=model_dict, Data=data_dict,
                                          OutInfo=out)

        elif ("map" in runtype):
            site_dict = \
                run_map(Ctrl=ctrl, Model=model_dict, Data=data_dict,
                                          OutInfo=out)


        else:
            error("run_inv_flightline: runtype "+runtype+" not implementd! Exit.")

#         Now store inversion results for this site:

        if out:
            print("site_dict: ",site_dict.keys())




        mtmp = site_dict["model"]
        dtmp = site_dict["data"]
        ctmp = site_dict["log"]

        # print(mtmp[0].shape)
        if icnt==0:
            site_num  = numpy.array([ii])
            site_conv = ctmp[1]
            site_nrms = ctmp[2]
            site_smap = ctmp[3]
            site_modl = mtmp[0]
            site_merr = mtmp[1]
            site_sens = mtmp[2]
            site_dobs = dtmp[0].reshape((1,-1))
            site_dcal = dtmp[1].reshape((1,-1))
            site_derr = dtmp[2].reshape((1,-1))
            # clog = numpy.hstack((ctmp[0], ctmp[1], ctmp[2], ctmp[3],
            #                   ctmp[4].ravel(),
            #                   ctmp[5].ravel(),
            #                   ctmp[6].ravel(),
            #                   ctmp[7].ravel()))
            # site_log[ii,0:len(clog)] = clog
            if uncert:
                jacd = site_dict["jacd"]
                site_jacd = jacd.reshape((1,numpy.size(jacd)))
                pcov = site_dict["cpost"]
                site_pcov = pcov.reshape((1,numpy.size(pcov)))
        else:
           site_num = numpy.vstack((site_num, ii))
           site_conv = numpy.vstack((site_conv, ctmp[1]))
           site_nrms = numpy.vstack((site_nrms, ctmp[2]))
           site_smap = numpy.vstack((site_smap, ctmp[3]))
           site_modl = numpy.vstack((site_modl, mtmp[0]))
           site_merr = numpy.vstack((site_merr, mtmp[1]))
           site_sens = numpy.vstack((site_sens, mtmp[2]))
           site_dobs = numpy.vstack((site_dobs, dtmp[0]))
           site_dcal = numpy.vstack((site_dcal, dtmp[1]))
           site_derr = numpy.vstack((site_derr, dtmp[2]))
           # clog = numpy.hstack((ctmp[0], ctmp[1], ctmp[2], ctmp[3],
           #                    ctmp[4].ravel(),
           #                    ctmp[5].ravel(),
           #                    ctmp[6].ravel(),
           #                    ctmp[7].ravel()))
           # site_log[ii,0:len(clog)] = clog

           if uncert:
               jacd = site_dict["jacd"]
               site_jacd = numpy.vstack((site_jacd,jacd.reshape((1,numpy.size(jacd)))))
               pcov = site_dict["cpost"]
               site_pcov = numpy.vstack((site_pcov, pcov.reshape((1,numpy.size(pcov)))))

# The _Ctrl_ paramter set as well as the results for data_dict set (flight line or area) are stored in _.npz_ files, which strings _"ctrl.npz"_ and _"results.npz"_ appended:

# +

    results_dict ={
        "fl_data" : result_file,
        "fl_name" : fl_name,
        "header" : header,
        # "site_log" :  site_log,
        "mod_ref" : mod_apr,
        "mod_act" : mod_act,
        "dat_act" : dat_act,
        "site_modl" : site_modl,
        "site_sens" : site_sens,
        "site_merr" : site_merr,
        "site_dobs" : site_dobs,
        "site_dcal" : site_dcal,
        "site_derr" : site_derr,
        "site_nrms" : site_nrms,
        "site_smap" : site_smap,
        "site_conv" : site_conv,
        "site_num" : site_num,
        "site_y" : site_y,
        "site_x" : site_x,
        "site_gps" : site_gps,
        "site_alt" : site_alt,
        "site_dem" : site_dem
        }

    if uncert:
        results_dict["site_jacd"] = site_jacd
        results_dict["site_pcov"] = site_pcov

    if out:
        print(list(results_dict.keys()))
        elapsed = (time.time() - start)
        print (" Used %7.4f sec for %6i sites" % (elapsed, ii+1))
        print (" Average %7.4f sec/site\n" % (elapsed/(ii+1)))


    if results_out:
        numpy.savez_compressed(result_file, **results_dict)
        print("\n\nResults stored to "+result_file)
    else:

        return results_dict

    # joblib.memory.clear()

# def run_tikh_ensemble(data_file=None,
#             prior_file=None,
#             result_strng=None,
#             ctrl=None,
#             results_out=True,
#             out=False):
#     """
#     Wrapper for data_dict parallel set inversion

#     Parameters
#     ----------

#     data_file : strng
#         Input data_dict file. The default is None.
#     result_file : strng
#             site_dict output file. The default is None.
#     prior_file : strng, optional
#         Prior file. Only required if the read option for priors
#         is set. The default is None.

#     ctrl : dict
#         Contains control variables for this run. The default is None.

#     out : TYPE, optional
#         DESCRIPTION. The default is True.

#     Returns
#     -------
#     results_dict : dict
#         Contains items which will be stored into resultfile.

#     vr, Bloomsday 2024


#     ctrl_dict ={
#         "system":
#             [AEM_system, FwdCall],
#         "header":
#             [titstrng, ""],
#         "inversion":
#             numpy.array([RunType, RegFun, Tau0, Tau1, Maxiter, ThreshRMS,
#                       LinPars, SetPrior, Delta, RegShift], dtype=object),
#         "covar":
#             numpy.array([L0, Cm0, L1, Cm1], dtype=object),
#         "uncert":
#             [Ensemble, Percentiles],

#         "data":
#             numpy.array([DataTrans, data_active, DatErr_add, DatErr_mult, Direction], dtype=object),
#         "model":
#             numpy.array([ParaTrans, mod_act, mod_apr, mod_var, mod_bnd], dtype=object),
#                 }


#     """
#     if data_file is None:
#         error("run_inv: No data_dict file given! Exit.")

#     name, ext = os.path.splitext(data_file)

#     if result_strng is None:
#         result_file = name+"_results.npz"
#     else:
#         result_file = name+result_strng+"_results.npz"


#     start = time.time()

#     AEM_system = ctrl["system"][0]
#     _ ,NN, _, _, _, = aesys.get_system_params(System=AEM_system)

#     print("\n Reading file " + data_file)
#     DataObs, Header, _ = aesys.read_aempy(File=data_file,
#                                    System=ctrl["system"][0], out=False)


#     fl_name = DataObs[0, 0]
#     ctrl["name"] = fl_name

#     print("site_dict: ",ctrl.keys())
#     ctrl_file = result_file.replace("_results", "_ctrl")
#     numpy.savez_compressed(ctrl_file,**ctrl)


#     site_x = DataObs[:, 1]
#     site_y = DataObs[:, 2]
#     site_gps = DataObs[:, 3]
#     site_alt = DataObs[:, 4]
#     site_dem = DataObs[:, 5]
#     dat_obs =  DataObs[:, 6:6+NN[2]]


#     """
#     Setup data-related parameter dict
#     """
#     [nsite,ndata] = numpy.shape(dat_obs)
#     sites = numpy.arange(nsite)

#     data_act = ctrl["data"][1]
#     data_err_add = ctrl["data"][2]
#     data_err_mult = ctrl["data"][3]

#     dat_act = numpy.tile(data_act,(nsite,1))
#     dat_err = numpy.zeros_like(dat_obs)
#     for ii in sites:
#         dat_err[ii, :], _ = set_errors(dat_obs[ii, :],
#                                             daterr_add=data_err_add,
#                                             daterr_mult=data_err_mult)
#         # print("\n",ii)
#         # print(dat_obs[ii, :])
#         # print(dat_err[ii, :])




#     """
#     Setup model-related parameter dict
#     """
#     runtype = ctrl["inversion"][0].lower()
#     setprior = ctrl["inversion"][7].lower()
#     maxiter =  ctrl["inversion"][4]

#     mod_act = ctrl["model"][1].copy()
#     mod_apr = ctrl["model"][2].copy()
#     mod_var = ctrl["model"][3].copy()
#     mod_bnd = ctrl["model"][4].copy()


#     site_prior = numpy.zeros((nsite,numpy.shape(mod_apr)[0]))

#     if "read" in setprior:
#         if name.lower() not in prior_file.lower():
#             error("Halfspace file name does not match! Exit.")
#         site_prior = load_prior(prior_file,
#                                        m_ref=mod_apr,
#                                        m_apr=mod_apr,
#                                        m_act=mod_act)
#     if "set" in setprior:
#         for ii in sites:
#                 site_prior[ii, :] = mod_apr


#     ensout = ctrl["uncert"][0]
#     percentiles = ctrl["uncert"][1]
#     header = ctrl["header"][0]

# # This is the main loop over sites in a flight line or within an area:

#     """
#     Loop over sites
#     """

#     # logsize = (2 + 7*maxiter)
#     # site_log = numpy.full((len(sites),logsize), numpy.nan)
#     mtmp = numpy.array([])
#     for ii in sites:
#         print("\n Invert site #"+str(ii)+"/"+str(len(sites)))

#         """
#         Setup parameter dict
#         """
#         data_dict = dict([
#             ("d_act", dat_act[ii,:]),
#             ("d_obs", dat_obs[ii,:]),
#             ("d_err", dat_err[ii,:]),
#             ("alt", site_alt[ii])
#             ])

#         # print(ii)
#         # print(dat_obs[ii,:])


#         if "read" in setprior:
#             mod_apr = site_prior[ii,:]
#             mod_ini = mod_apr.copy()

#         elif "upd" in setprior:

#             if ii == 0:
#                 mod_ini = site_prior[ii,:]
#                 mod_apr = mod_ini.copy()
#             else:
#                 mod_ini = mtmp[0].copy()
#                 mod_apr = mod_ini.copy()

#         elif "set" in setprior:
#             mod_apr = site_prior[ii,:]
#             mod_ini = mod_apr.copy()


#         # print("\n",ii)
#         # print(mod_ini[ii, :])
#         # print(mod_apr)
#         # print(mod_var[ii, :])

#         model_dict = dict([
#             ("m_act", mod_act),
#             ("m_apr", mod_apr),
#             ("m_var", mod_var),
#             ("m_bnd", mod_bnd),
#             ("m_ini", mod_ini)
#             ])


#         """
#         Call inversion algorithms
#         """
#         if "opt" in runtype:
#             ens_dict =\
#                 run_tikh_opt(Ctrl=ctrl, Model=model_dict, Data=data_dict,
#                                   OutInfo=out)

#         if "occ" in runtype:
#             ens_dict =\
#                 run_tikh_occ(Ctrl=ctrl, Model=model_dict, Data=data_dict,
#                                   OutInfo=out)

#         if "map" in runtype:
#             ens_dict =\
#                 run_map(Ctrl=ctrl, Model=model_dict, Data=data_dict,
#                                   OutInfo=out)

# #         Now store inversion results for this site:
#         if out:
#             print("ens_dict: ",ens_dict.keys())


#         M = ens_dict["model"]
#         D = ens_dict["data"]
#         C = ens_dict["log"]

#         if ii==0:
#             ens_num  = numpy.array([ii])
#             ens_nrms = C[2]
#             ens_modl = M[0]
#             ens_merr = M[1]
#             ens_sens = M[2]
#             ens_dobs = D[0].reshape((1,-1))
#             ens_dcal = D[1].reshape((1,-1))
#             ens_derr = D[2].reshape((1,-1))

#         else:
#            ens_num = numpy.vstack((ens_num, ii))
#            ens_nrms = numpy.vstack((ens_nrms, C[2]))

#            ens_modl = numpy.vstack((ens_modl, M[0]))
#            ens_merr = numpy.vstack((ens_merr, M[1]))
#            ens_sens = numpy.vstack((ens_sens, M[2]))
#            ens_dobs = numpy.vstack((ens_dobs, D[0]))
#            ens_dcal = numpy.vstack((ens_dcal, D[1]))
#            ens_derr = numpy.vstack((ens_derr, D[2]))


#     m_quants, m_mean, m_stdv, m_skew, m_kurt, m_mode = \
#         calc_stat_ens(ensemble=ens_modl, quantiles=percentiles, sum_stats=True)
#     stat_modl = numpy.vstack((m_quants, m_mean, m_stdv, m_skew, m_kurt, m_mode))

#     d_quants, d_mean, d_stdv, d_skew, d_kurt, d_mode = \
#         calc_stat_ens(ensemble=ens_modl, quantiles=percentiles, sum_stats=True)
#     stat_dcal = numpy.vstack((d_quants, d_mean, d_stdv, d_skew, d_kurt, d_mode))

#     r_quants, r_mean, r_stdv, r_skew, r_kurt, r_mode = \
#         calc_stat_ens(ensemble=ens_nrms, quantiles=percentiles, sum_stats=True)
#     stat_nrms= numpy.vstack((r_quants, r_mean, r_stdv, r_skew, r_kurt, r_mode))

#     mod_alt =  data_dict["alt"]

#     if not ensout:

#         ens_dict ={
#             "fl_data" : result_file,
#             "fl_name" : fl_name,
#             "header" : header,
#             "ens_log" :  ens_log,
#             "mod_ref" : mod_apr,
#             "mod_act" : mod_act,
#             "dat_act" : dat_act,
#             "ens_modl" : ens_modl,
#             "ens_sens" : ens_sens,
#             "ens_merr" : ens_merr,
#             "ens_dobs" : ens_dobs,
#             "ens_dcal" : ens_dcal,
#             "ens_derr" : ens_derr,
#             "ens_nrms" : ens_nrms,
#             "ens_smap" : ens_smap,
#             "ens_conv" : ens_conv,
#             "ens_num" : ens_num,
#             "ens_y" : ens_y,
#             "ens_x" : ens_x,
#             "ens_gps" : ens_gps,
#             "ens_alt" : mod_alt,
#             "ens_dem" : mod_dem
#             }

#         if out:
#             print("ens_dict: ",ens_dict.keys())
#     else:
#         print("not implemeted yet!!!!")

#     if out:
#         print(list(ens_dict_dict.keys()))
#         elapsed = (time.time() - start)
#         print (" Used %7.4f sec for %6i sites" % (elapsed, ii+1))
#         print (" Average %7.4f sec/site\n" % (elapsed/(ii+1)))

#     if ens_dict_out:
#         numpy.savez_compressed(result_file, **ens_dict)
#         print("\n\nens_dict stored to "+result_file)
#     else:

#         return ens_dict

def run_tikh_opt(Ctrl=None, Model=None, Data=None, OutInfo=False):
    """
    Tikhonov inversion with optimal tau, using several popular
    methods of determining the regularization parameter

    vr July 2022

    Hansen, P. C.
    Rank Deficient and Discrete Ill-Posed Problems
    SIAM, Philadelphia, 1998

    Vogel, C.
    Computational Methods for Inverse Problems
    SIAM, 2002

    Doicu, A.; Trautmann, T. & Schreier, F.
    Numerical regularization for atmospheric inverse problems
    Springer, 2010

    Hansen, Per Christian
    Discrete Inverse Problems: Insight and Algorithms
    SIAM, 2010

    Aster, R. C.; Borchers, B. & Thurber, C. H.
    Parameter Estimation and Inverse Problems
    Elsevier, 2019

    """

    """
    unpack control block:

    Ctrl_TikhOpt = dict([
        ("system", [AEM_system, FwdCall]),
        ("inversion",[InvType, RegFun, Tau0, Tau1, Maxiter,thresh, LinPars, SetPrior, Delta, ]),
        ("covar", [L0, Cm0, L1, Cm1]),
        ("transform", [ParaTrans, DataTrans]),
        ("uncert", True)
        ])
    """

    system, fwdcall = Ctrl["system"]
    invtype, regfun, tau0, tau1, maxiter, thresh, linepars, setprior, delta, gshift = Ctrl[
        "inversion"]
    L0, Cm0, L1, Cm1 = Ctrl["covar"]
    uncert = Ctrl["uncert"]
    profname = Ctrl["name"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    nreg = numpy.size(tau0)*numpy.size(tau1)

    if linepars == []:
        do_linesearch = False
    else:
        do_linesearch = True
        maxreduce = linepars[0]
        facreduce = linepars[1]

    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_act = Data["d_act"]  # .reshape(-1,1)
    d_obs = Data["d_obs"]  # .reshape(-1,1)
    d_err = Data["d_err"]  # .reshape(-1,1)
    alt = Data["alt"]

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    d_state = 0
    d_obs, d_err, dobs_state = transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)


    obs = extract_dat(d_obs, d_act)
    err = extract_dat(d_err, d_act)

    Wd = numpy.diagflat(1.0/err, 0)
    Wd = scipy.sparse.csr_matrix(Wd)
    Cdi = numpy.diagflat(1.0/err**2, 0)
    Cdi = scipy.sparse.csr_matrix(Cdi)

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    """
    unpack model block
    Model =\
    [model_act, model_prior, model_var, model_bounds, model_ini, ParaTrans]
    """
    m_act = Model["m_act"]
    m_apr = Model["m_apr"]
    m_var = Model["m_var"]
    m_bnd = Model["m_bnd"]
    m_ini = Model["m_ini"]

    m_err = numpy.sqrt(m_var)

    # print(m_apr)
    m_state = 0
    m_apr, _ = transform_parameter(
        m_vec=m_apr, m_trn=m_trn, m_state=m_state, mode="f")
    m_ini, _ = transform_parameter(
        m_vec=m_ini, m_trn=m_trn, m_state=m_state, mode="f")
    m_state = m_trn

    mpara = numpy.shape(numpy.flatnonzero(m_act))[0]

    merr = numpy.nan * numpy.ones_like(m_apr)
    sens = numpy.nan * numpy.ones_like(m_apr)

    model = m_ini.copy()
    model_old = model.copy()
    model_error = m_err
    model_error_old = model_error

    """
    start inversion loop
    """
    niter = -1
    dfit_iter = 1.0e30
    dfit_old = dfit_iter

    while (niter < maxiter):

        niter = niter + 1
        d_cal, dcal_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                  m_vec=model,
                                                  m_trn=m_trn, m_state=m_state,
                                                  d_trn=d_trn)

        nrmse_iter, smape_iter = calc_datafit(data_obs=d_obs,
                                                      data_cal=d_cal,
                                                      data_err=d_err,
                                                      data_act=d_act)
        # print(niter,"####",model, m_trn, m_state,d_trn)
        # print("dfit0", nrmse_iter, smape_iter)
        if niter == 0:
            conv_status = 1
            print(" ThreshVal =", thresh)
            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter
                dfit_old = nrmse_iter
                dfit_0 = nrmse_iter


            if "smp" in thresh[3]:
                dfit_iter = smape_iter
                dfit_old = smape_iter
                dfit_0 = smape_iter

            model_old = model.copy()
            dnorm_iter = numpy.array([calc_dnorm(data_obs=d_obs, data_cal=d_cal,
                                                         data_err=d_err, data_act=d_act)])
            # print(model)
            mnorm_iter = numpy.array([scipy.linalg.norm(model)])
            rvals_iter = numpy.array([0., 0.])
            dfits_iter = numpy.array([nrmse_iter, smape_iter])

            if OutInfo:
                print("Starting NRMSE      =  %7.3f,  SMAPE = %4.1f percent"
                      % (nrmse_iter, smape_iter))
        else:

            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter

            if "smp" in thresh[3]:
                dfit_iter = smape_iter

            if dfit_iter < dfit_old:
                if OutInfo == True:
                    print("Iteration %6i NRMSE =  %7.3f" % (niter, dfit_iter))

                dfit_change = numpy.abs((dfit_iter - dfit_old)/dfit_old)
                modl_change = scipy.linalg.norm((model - model_old)/model_old)

                if (dfit_iter < thresh[0]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i data fit =  %7.3f <= threshold = %7.3f percent" %
                              (niter, dfit_iter, thresh[0]))
                    break

                if (dfit_change < thresh[1]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i data fit change =  %7.3f <= threshold = %4.1f" %
                              (niter, dfit_change, thresh[1]))
                    break

                if (modl_change < thresh[2]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i model change =  %7.3f <= threshold = %4.1f" %
                              (niter, modl_change, thresh[2]))
                    break

                if "rms" in thresh[3]:
                    dfit_old = nrmse_iter
                if "smp" in thresh[3]:
                    dfit_old = smape_iter
                model_old = model.copy()

            else:

                if OutInfo == True:
                    print("Iteration %6i dfit =  %8.4f >= dfit_old = %8.4f." %
                          (niter, dfit_iter, dfit_old))

                dfit_iter = dfit_old
                model = model_old.copy()
                break

        Jac = calc_jac(fwdcall=fwdcall, alt=alt,
                               m_vec=model, m_act=m_act, m_trn=m_trn, m_state=m_state,
                               d_vec=d_obs, d_act=d_act, d_trn=d_trn,
                               delta=delta, scalejac=False, out=False)

        cal = extract_dat(d_cal, d_act)
        # print("cal:", cal)
        # print("obs:", obs)

        m_iter = extract_mod(model, m_act)
        m_apri = extract_mod(m_apr, m_act)

        diff_m = m_iter - m_apri

        Jd = Wd@Jac
        JJ = Jd.T@Jd

        sensi = calc_sensitivity(Jac=Jd)

        m_test = numpy.zeros((nreg, mpara))
        m_err = numpy.zeros((nreg, mpara))


        reg_choice = numpy.zeros((nreg, 1))

        dnorm = numpy.zeros((nreg, 1))
        mnorm = numpy.zeros((nreg, 1))
        tau = numpy.zeros((nreg, 2))

        itest = -1

        for t0 in tau0:
            Cmi0 = t0 * Cm0
            for t1 in tau1:
                Cmi1 = t1 * Cm1

                model_test = model.copy()
                itest = itest + 1
                tau[itest, :] = [t0, t1]

                A = JJ + Cmi0 + Cmi1
                r = Jac.T@Cdi@(cal - obs).T+Cmi0@diff_m.T+Cmi1@diff_m.T

                m_delta = scipy.linalg.solve(A, r)
                m_test[itest] = m_iter - m_delta
                model_test = insert_mod(M=model_test, m=m_test[itest],
                                                m_act=m_act)

                cali, d_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model_test,
                                                      m_trn=m_trn, m_state=m_state, d_trn=d_trn, d_act=d_act)

                r_test = Wd@(obs - cali[d_act != 0]).T


                dnorm[itest] = scipy.linalg.norm(r_test)
                mnorm[itest] = scipy.linalg.norm(m_test[itest])

                """
                 Cov a posteriori & Generalized Inverse
                """
                C = scipy.linalg.inv(A)
                m_err[itest] = numpy.sqrt(C.diagonal())
                G = C@Jd.T

                """
                Calc target paramter for tau optimization
                """

                if "ufc" in regfun.lower():
                    reg_choice[itest] = calc_ufc(
                        dnorm[itest], mnorm[itest])

                elif "upr" in regfun.lower():
                    M = Jd@G
                    reg_choice[itest] = calc_upr(
                        dnorm[itest], M, d_err)

                elif "gcv" in regfun.lower():
                    M = Jd@G
                    if iter==0:
                        reg_choice[itest] = calc_gcv(dnorm[itest], M)
                    else:
                        reg_choice[itest] = numpy.amin(0.5*rvals_iter[iter-1], calc_gcv(dnorm[itest], M))

                elif "mle" in regfun.lower():
                    M = Jd@G
                    reg_choice[itest] = calc_mle(dnorm[itest], M)

                elif "fix" in regfun.lower() or "lcc" in regfun.lower():
                    pass

                else:
                    error("Regularisation method "+regfun.lower() +
                          " not yet implemented! Exit.")

        if "fix" in regfun.lower():
            g_index = 0

        elif "lcc" in regfun.lower():
            g_index = calc_lc_corner(dnorm, mnorm)+gshift
            g_index = g_index.item()
            g_index = numpy.amax([g_index, 0])
            g_index = numpy.amin([g_index, numpy.shape(tau)[0]-1])
        elif any(s in regfun.lower() for s in ["gcv", "upr", "ufc", "mle"]):

            g_index = numpy.argmin(reg_choice, axis=0)+gshift
            g_index = g_index.item()
            g_index = numpy.amax([g_index, 0])
            g_index = numpy.amin([g_index, numpy.shape(tau)[0]-1])

        else:
            error(regfun.lower() + " not implemented!")



        reg_iter = [tau[g_index, 0].item(), tau[g_index, 1].item()]


        mdl = m_test[g_index, :]
        model = insert_mod(M=model, m=mdl, m_act=m_act)

        """
        Line Search
        """
        if do_linesearch:
            model, dfit_iter = run_linesearch(fwdcall, alt,
                                                      d_obs=d_obs, d_err=d_err, d_trn=d_trn, d_act=d_act, d_state=d_state,
                                                      model=model, m_delta=m_delta, m_act=m_act, m_trn=m_trn, m_state=m_state,
                                                      dfit=dfit_iter, mdfit=thresh[3],
                                                      facreduce=facreduce, maxreduce=maxreduce, out=True) #OutInfo)
            d_cal, _ = calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model,
                                                      m_trn=m_trn, m_state=m_state, d_trn=d_trn, d_act=d_act)

            nrmse_iter, smape_iter = calc_datafit(data_obs=d_obs,
                                                          data_cal=d_cal,
                                                          data_err=d_err,
                                                          data_act=d_act)
            dnorm_ii= scipy.linalg.norm(Wd@(obs - d_cal[d_act != 0]).T)
            mnorm_ii = scipy.linalg.norm(model)


        else:
            dnorm_ii = dnorm[g_index]
            mnorm_ii = mnorm[g_index]

        if niter > 0:
            dnorm_iter = numpy.append(dnorm_iter, dnorm_ii)
            mnorm_iter = numpy.append(mnorm_iter, mnorm_ii)
            rvals_iter = numpy.append(rvals_iter, reg_iter)
            dfits_iter = numpy.append(dfits_iter, [nrmse_iter, smape_iter])


    # if OutInfo==True:
    print(" %s initial model: dfit =  %7.3f"
          % (profname, dfit_0))
    print(" %s final model:  iter  %6i NRMSE =  %7.3f,  SMAPE = %4.1f percent, %s RegPars are %10.4g / %10.4g"
          % (profname, niter, nrmse_iter, smape_iter, regfun.lower(),  tau[g_index, 0],  tau[g_index, 1]))
    # print(numpy.shape(sensi))
    sens = sensi  # insert_mod(M=sens, m=sensi,m_act=m_act)
    modl, m_state = transform_parameter(m_vec=model, m_trn=m_trn, m_state=m_state,  mode="b")
    modl = extract_mod(modl, m_act)
    merr = m_err[g_index, :].flat
    # insert_mod(M=merr, m=tmp, m_act=m_act)
    # print(numpy.shape(modl),numpy.shape(merr),numpy.shape(sens))
    print()
    results = \
        dict([
            ("model", [modl, merr, sens, m_state, m_act]),
            ("data", [d_obs, d_cal, d_err, d_state, d_act]),
            ("log", [niter, conv_status, nrmse_iter, smape_iter,
             dnorm_iter, mnorm_iter, rvals_iter, dfits_iter]),
        ])

    # print(type(dnorm_iter))

    if uncert:
        A = JJ + Cm0.multiply(reg[0]) + Cm1.multiply(reg[1])
        # Cov a posteriori & Generalized Inverse
        C = scipy.linalg.inv(A)
        E = numpy.sqrt(C.diagonal())

        G = C@Jd.T

        # Resolution martices & spread
        Rm = G@Jd
        Nm = numpy.sum(Rm.diagonal())
        Sm = scipy.linalg.norm(numpy.identity(numpy.shape(Rm)[0])-Rm)

        Rd = Jd@G
        Nd = numpy.sum(Rd.diagonal())
        Sd = scipy.linalg.norm(numpy.identity(numpy.shape(Rd)[0])-Rd)

        uncpars =\
            dict([
                ("jacd", Jd),               # jacobian
                ("cpost", C),               # cov a-post
                ("merr", E),                # model error
                ("mresm", [Rm, Sm, Nm]),    # model resolution
                ("dresm", [Rd, Sd, Nd]),    # data resolution
                ("gi", G),                  # generalized inverse
            ])

        results.update(uncpars)

    return results

# @ray.remote


def run_tikh_occ(Ctrl=None, Model=None, Data=None, OutInfo=False):
    """
    Tikhonov/MAP inversion with tau, determined by cooling scheme (Occam)

    vr July 2022

    Constable, S. C.; Parker, R. L. & Constable, C. G.
    Occam's inversion: a practical algorithm for generating smooth
    models from electromagnetic sounding data
    Geophysics, 1987, 52

    DeGroot‐Hedlin, C. & Constable, S.
    Occam's inversion to generate smooth, two‐dimensional models from
    magnetotelluric data
    Geophysics, 1990, 55, 1613-1624


    """

    """
    unpack control block:

    Ctrl_Occam = dict([
        ("system", [AEM_system, FwdCall]),
        ("inversion",[InvType, TauSeq, Maxiter,thresh, LinPars, SetPrior, Delta, ]),
        ("covar", [L0, Cm0, L1, Cm1]),
        ("transform", [ParaTrans, DataTrans]),
        ("uncert", True)
        ])
    """

    system, fwdcall = Ctrl["system"]
    invtype, tauseq, tau, maxiter, thresh, linepars, setprior, delta = Ctrl["inversion"]

    L0, Cm0, L1, Cm1 = Ctrl["covar"]

    uncert = Ctrl["uncert"]
    profname = Ctrl["name"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    if linepars == []:
        do_linesearch = False
    else:
        do_linesearch = True
        maxreduce = linepars[0]
        facreduce = linepars[1]

    if len(tauseq) == 2:
        set_taustart = False
        taustart = tauseq[1]
        taufac = tauseq[0]
    if len(tauseq) == 1:
        set_taustart = True
        taufac = tauseq[0]
    if len(tauseq) == 0:
        set_taustart = True
        taufac = 0.66
    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_act = Data["d_act"]  # .reshape(-1,1)
    d_obs = Data["d_obs"]  # .reshape(-1,1)
    d_err = Data["d_err"]  # .reshape(-1,1)
    alt = Data["alt"]

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    d_state = 0
    d_obs, d_err, dobs_state = transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)

    obs = extract_dat(d_obs, d_act)
    err = extract_dat(d_err, d_act)
    Wd = numpy.diagflat(1.0/err, 0)
    Wd = scipy.sparse.csr_matrix(Wd)
    Cdi = numpy.diagflat(1.0/err**2, 0)
    Cdi = scipy.sparse.csr_matrix(Cdi)

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    """
    unpack model block
    Model =\
    [model_act, model_prior, model_var, model_bounds, model_ini, ParaTrans]
    """
    m_act = Model["m_act"]
    m_apr = Model["m_apr"]
    m_var = Model["m_var"]
    m_bnd = Model["m_bnd"]
    m_ini = Model["m_ini"]
    m_err = numpy.sqrt(m_var)

    m_state = 0
    m_apr, _ = transform_parameter(
        m_vec=m_apr, m_trn=m_trn, m_state=m_state, mode="f")
    m_ini, _ = transform_parameter(
        m_vec=m_ini, m_trn=m_trn, m_state=m_state, mode="f")
    m_state = m_trn

    mpara = numpy.shape(numpy.flatnonzero(m_act))[0]

    merr = numpy.nan * numpy.ones_like(m_apr)
    sens = numpy.nan * numpy.ones_like(m_apr)

    model = m_ini.copy()
    model_old = model.copy()
    model_error = m_err
    model_error_old = model_error

    """
    start inversion loop
    """

    niter = -1
    while (niter < maxiter):

        niter = niter + 1
        d_cal, dcal_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                  m_vec=model,
                                                  m_trn=m_trn, m_state=m_state,
                                                  d_trn=d_trn)

        nrmse_iter, smape_iter = calc_datafit(data_obs=d_obs,
                                                      data_cal=d_cal,
                                                      data_err=d_err,
                                                      data_act=d_act)

        if niter == 0:
            conv_status = 1

            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter
                dfit_old = nrmse_iter
                dfit_0 = nrmse_iter

            if "smp" in thresh[3]:
                dfit_iter = smape_iter
                dfit_old = smape_iter
                dfit_0 = smape_iter

            model_old = model.copy()
            dnorm_iter = numpy.array([calc_dnorm(data_obs=d_obs, data_cal=d_cal,
                                                         data_err=d_err, data_act=d_act)])
            mnorm_iter = numpy.array([scipy.linalg.norm(model)])
            rvals_iter = numpy.array([0., 0.])
            dfits_iter = numpy.array([nrmse_iter, smape_iter])

            if OutInfo:
                print("Starting NRMSE      =  %7.3f,  SMAPE = %4.1f percent"
                      % (nrmse_iter, smape_iter))
        else:

            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter

            if "smp" in thresh[3]:
                dfit_iter = smape_iter

            if dfit_iter < dfit_old:
                if OutInfo == True:
                    print("Iteration %6i NRMSE =  %7.3f" % (niter, dfit_iter))

                dfit_change = numpy.abs((dfit_iter - dfit_old)/dfit_old)
                modl_change = scipy.linalg.norm((model - model_old)/model_old)

                if (dfit_iter < thresh[0]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i data fit =  %7.3f <= threshold = %7.3f percent" %
                              (niter, dfit_iter, thresh[0]))
                    break

                if (dfit_change < thresh[1]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i data fit change =  %7.3f <= threshold = %4.1f" %
                              (niter, dfit_change, thresh[1]))
                    break

                if (modl_change < thresh[2]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i model change =  %7.3f <= threshold = %4.1f" %
                              (niter, modl_change, thresh[2]))
                    break

                if "rms" in thresh[3]:
                    dfit_old = nrmse_iter
                if "smp" in thresh[3]:
                    dfit_old = smape_iter
                model_old = model.copy()

            else:

                if OutInfo == True:
                    print("Iteration %6i dfit =  %8.4f >= dfit_old = %8.4f." %
                          (niter, dfit_iter, dfit_old))

                dfit_iter = dfit_old
                model = model_old.copy()
                break

        Jac = calc_jac(fwdcall=fwdcall, alt=alt,
                               m_vec=model, m_act=m_act, m_trn=m_trn, m_state=m_state,
                               d_vec=d_obs, d_act=d_act, d_trn=d_trn,
                               delta=delta, scalejac=False, out=False)

        cal = extract_dat(d_cal, d_act)

        m_iter = extract_mod(model, m_act)
        m_apri = extract_mod(m_apr, m_act)

        diff_m = m_iter - m_apri

        Jd = Wd@Jac
        JJ = Jd.T@Jd
        # JJT = Jac.T@Cdi@Jac

        sensi = calc_sensitivity(Jac=Jd)

        # Cmi0 = Cm0.multiply(tau)

        if niter == 0 and set_taustart:
            tau_start = calc_regstart(D=JJ, M=Cm1)
            tau_iter = tau_start
            print("initial tau1: ", tau_iter)
        else:
            tau_iter = tau_iter*taufac
            print("tau1: ", tau_iter)

        Cmi1 = Cm1.multiply(tau_iter)

        model_test = model.copy()

        A = JJ + Cmi1
        r = Jac.T@Cdi@(obs - cal).T+Cmi1@diff_m.T

        m_delta = scipy.linalg.solve(A, r)
        m_test = m_iter - m_delta
        model_test = insert_mod(M=model_test, m=m_test,
                                        m_act=m_act)

        cali, d_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model_test,
                                              m_trn=m_trn, m_state=m_state, d_trn=d_trn, d_act=d_act)

        r_test = Wd@(obs - cali).T

        dnorm = scipy.linalg.norm(r_test)
        mnorm = scipy.linalg.norm(m_test)

    reg = [tau_iter]
    mdl = m_test
    model = insert_mod(M=model, m=mdl, m_act=m_act)

    """
    Line Search
    """
    if do_linesearch:
        model, dfit_iter = run_linesearch(fwdcall, alt,
                                        d_obs=d_obs, d_err=d_err, d_trn=d_trn, d_act=d_act, d_state=d_state,
                                        model=model, m_delta=m_delta, m_act=m_act, m_trn=m_trn, m_state=m_state,
                                        dfit=dfit_iter, mdfit=thresh[3],
                                        facreduce=facreduce, maxreduce=maxreduce, out=OutInfo)
        nrmse_iter, smape_iter = calc_datafit(data_obs=d_obs,
                                                      data_cal=d_cal,
                                                      data_err=d_err,
                                                      data_act=d_act)
    if niter > 0:
        dnorm_iter = numpy.append(dnorm_iter, dnorm)
        mnorm_iter = numpy.append(mnorm_iter, mnorm)
        rvals_iter = numpy.append(rvals_iter, reg)
        dfits_iter = numpy.append(dfits_iter, [nrmse_iter, smape_iter])

    # if OutInfo==True:
    print(" %s initial model: dfit =  %7.3f"
          % (profname, dfit_0))
    print(" %s final model:  iter  %6i NRMSE =  %7.3f,  SMAPE = %4.1f percent, occ RegPar is %10.4g "
          % (profname, niter, nrmse_iter, smape_iter, tau_iter))
    # print(numpy.shape(sensi))
    sens = sensi  # insert_mod(M=sens, m=sensi,m_act=m_act)
    modl, m_state = transform_parameter(
        m_vec=model, m_trn=m_trn, m_state=m_state,  mode="b")
    modl = extract_mod(modl, m_act)
    merr = m_err.flat
    # insert_mod(M=merr, m=tmp, m_act=m_act)
    # print(numpy.shape(modl),numpy.shape(merr),numpy.shape(sens))
    print()
    results = \
        dict([
            ("model", [modl, merr, sens, m_state, m_act]),
            ("data", [d_obs, d_cal, d_err, d_state, d_act]),
            ("log", [niter, conv_status, nrmse_iter, smape_iter,
             dnorm_iter, mnorm_iter, rvals_iter, dfits_iter]),
        ])

    if uncert:

        A = JJ + Cm1.multiply(tau_iter)
        # Cov a posteriori & Generalized Inverse
        C = scipy.linalg.inv(A)
        E = numpy.sqrt(C.diagonal())
        G = C@Jd.T

        # Resolution martices & spread
        Rm = G@Jd
        Nm = numpy.sum(Rm.diagonal())
        Sm = scipy.linalg.norm(numpy.identity(numpy.shape(Rm)[0])-Rm)

        Rd = Jd@G
        Nd = numpy.sum(Rd.diagonal())
        Sd = scipy.linalg.norm(numpy.identity(numpy.shape(Rd)[0])-Rd)

        uncpars =\
            dict([
                ("jacd", Jd),               # jacobian
                ("cpost", C),               # cov a-post
                ("merr", E),                # model error
                ("mresm", [Rm, Sm, Nm]),    # model resolution
                ("dresm", [Rd, Sd, Nd]),    # data resolution
                ("gi", G),                  # generalized inverse
            ])

        results.update(uncpars)

    return results


# @ray.remote
def run_map(Ctrl=None, Model=None, Data=None, OutInfo=False):
    """
    Maximum aposteriori (MAP) inversion with optimal tau, using
    several popular methods of determining the regularization parameter

    vr July 2022

    Tarantola, A. & Valette, B.
    Inverse problem = quest for information
    J. Geophysics, 1982, 50, 159-170

    Tarantola, A. & Valette, B.
    Generalized nonlinear inverse problems solved
    using the least squares criterion
    Rev. Geophys. Space Phys., 1982, 20, 219-232

    Rodgers, C. D.
    Inverse Methods for atmospheric sounding
    World Scientific, 2000

    Siripunvaraporn, W. & Egbert, G.
    An efficient data-subspace inversion method for
    2-D magnetotelluric data
    Geophysics, 2000, 65, 791-803

    Siripunvaraporn, W.; Egbert, G. & Lenbury, Y.
    Three-dimensional magnetotelluric inversion: data-space method
    Phys. Earth Planet. Inter., 2005, 150
    doi:10.1016/j.pepi.2004.08.023

    Tarantola, A.
    Inverse Problem Theory and Methods for Model Parameter Estimation
    SIAM, 2005

    Aster, R. C.; Borchers, B. & Thurber, C. H.
    Parameter Estimation and Inverse Problems
    Elsevier, 2019

"""

    """
    unpack control block:

    Ctrl_MAP = dict([
        ("system", [AEM_system, FwdCall]),
        ("inversion",
         [InvType, InvSpace, RegFun, Tau1, Maxiter,ThreshRMS, LinPars, SetPrior, Delta, RegShift]),
        ("covar", [C, sC]),
        ("transform", [DataTrans, ParaTrans]),
        ("uncert", True)
       ])

    """

    system, fwdcall = Ctrl["system"]
    invtype, invspace, regfun, tau, maxiter, thresh, linepars, setprior, delta, gshift = Ctrl[
        "inversion"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    if "par" in invspace.lower():
        Cmi, CmiS = Ctrl["covar"]
    else:
        Cm, CmS = Ctrl["covar"]

    nreg = numpy.size(tau)

    if linepars == []:
        do_linesearch = False
    else:
        do_linesearch = True
        maxreduce = linepars[0]
        facreduce = linepars[1]

    uncert = Ctrl["uncert"]
    profname = Ctrl["name"]

    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_obs = Data["d_obs"]  # .reshape(-1,1)
    d_act = Data["d_act"]  # .reshape(-1,1)
    d_err = Data["d_err"]  # .reshape(-1,1)
    alt = Data["alt"]

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    d_state = 0
    d_obs, d_err, dobs_state = transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)

    obs = extract_dat(d_obs, d_act)
    err = extract_dat(d_err, d_act)

    Cd = numpy.diagflat(err**2, 0)
    Cd = scipy.sparse.csr_matrix(Cd)
    Sd = numpy.diagflat(err, 0)
    Cdi = numpy.diagflat(1.0/err**2, 0)
    Cdi = scipy.sparse.csr_matrix(Cdi)
    Sdi = numpy.diagflat(1.0/err, 0)

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    """
    unpack model block
    Model =\
    [model_act, model_prior, model_var, model_bounds, model_ini, ParaTrans]
    """
    m_act = Model["m_act"]
    m_apr = Model["m_apr"]
    m_var = Model["m_var"]
    m_bnd = Model["m_bnd"]
    m_ini = Model["m_ini"]
    m_err = numpy.sqrt(m_var)

    m_state = 0
    m_apr, _ = transform_parameter(
        m_vec=m_apr, m_trn=m_trn, m_state=m_state, mode="f")
    m_ini, _ = transform_parameter(
        m_vec=m_ini, m_trn=m_trn, m_state=m_state, mode="f")
    m_state = m_trn

    mpara = numpy.shape(numpy.flatnonzero(m_act))[0]

    merr = numpy.nan * numpy.ones_like(m_apr)
    sens = numpy.nan * numpy.ones_like(m_apr)

    model = m_ini.copy()
    model_old = model.copy()
    model_error = m_err
    model_error_old = model_error

    """
    start inversion loop
    """
    niter = -1

    while (niter < maxiter):

        niter = niter + 1
        d_cal, dcal_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                  m_vec=model,
                                                  m_trn=m_trn, m_state=m_state,
                                                  d_trn=d_trn)

        nrmse_iter, smape_iter = calc_datafit(data_obs=d_obs,
                                                      data_cal=d_cal,
                                                      data_err=d_err,
                                                      data_act=d_act)

        if niter == 0:
            conv_status = 1

            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter
                dfit_old = nrmse_iter
                dfit_0 = nrmse_iter

            if "smp" in thresh[3]:
                dfit_iter = smape_iter
                dfit_old = smape_iter
                dfit_0 = smape_iter

            model_old = model.copy()
            dnorm_iter = numpy.array([calc_dnorm(data_obs=d_obs, data_cal=d_cal,
                                                         data_err=d_err, data_act=d_act)])
            mnorm_iter = numpy.array([scipy.linalg.norm(model)])
            rvals_iter = numpy.array([0., 0.])
            dfits_iter = numpy.array([nrmse_iter, smape_iter])

            if OutInfo:
                print("Starting NRMSE      =  %7.3f,  SMAPE = %4.1f percent"
                      % (nrmse_iter, smape_iter))
        else:

            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter

            if "smp" in thresh[3]:
                dfit_iter = smape_iter

            if dfit_iter < dfit_old:
                if OutInfo == True:
                    print("Iteration %6i NRMSE =  %7.3f" % (niter, dfit_iter))

                dfit_change = numpy.abs((dfit_iter - dfit_old)/dfit_old)
                modl_change = scipy.linalg.norm((model - model_old)/model_old)

                if (dfit_iter < thresh[0]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i data fit =  %7.3f <= threshold = %7.3f percent" %
                              (niter, dfit_iter, thresh[0]))
                    break

                if (dfit_change < thresh[1]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i data fit change =  %7.3f <= threshold = %4.1f" %
                              (niter, dfit_change, thresh[1]))
                    break

                if (modl_change < thresh[2]):
                    conv_status = 1
                    if OutInfo == True:
                        print("Iteration %6i model change =  %7.3f <= threshold = %4.1f" %
                              (niter, modl_change, thresh[2]))
                    break

                if "rms" in thresh[3]:
                    dfit_old = nrmse_iter
                if "smp" in thresh[3]:
                    dfit_old = smape_iter
                model_old = model.copy()

            else:

                if OutInfo == True:
                    print("Iteration %6i dfit =  %8.4f >= dfit_old = %8.4f." %
                          (niter, dfit_iter, dfit_old))

                dfit_iter = dfit_old
                model = model_old.copy()
                break

        Jac = calc_jac(fwdcall=fwdcall, alt=alt,
                               m_vec=model, m_act=m_act, m_trn=m_trn, m_state=m_state,
                               d_vec=d_obs, d_act=d_act, d_trn=d_trn,
                               delta=delta, scalejac=False, out=False)

        cal = extract_dat(d_cal, d_act)

        m_iter = extract_mod(model, m_act)
        m_apri = extract_mod(m_apr, m_act)

        diff_m = m_iter - m_apri

        Jd = Sd@Jac
        sensi = calc_sensitivity(Jac=Jd)

        m_test = numpy.zeros((nreg, mpara))
        m_err = numpy.zeros((nreg, mpara))

        reg_choice = numpy.zeros((nreg, 1))
        dnorm = numpy.zeros((nreg, 1))
        mnorm = numpy.zeros((nreg, 1))
        # tau = numpy.zeros((nreg,1))

        itest = -1

        for tt in tau:

            model_test = model.copy()
            itest = itest + 1
            if "par" in invspace.lower():
                # print(numpy.shape(JJ), numpy.shape(Cmi), tt)
                JJ = Jac.T@Cdi@Jac
                A = JJ + tt*Cmi.todense()
                r = Jac.T@Cdi@((obs - cal).T+Jac@diff_m.T)
                m_delta = scipy.linalg.solve(A, r)
                # cov a posteriori
                C = scipy.linalg.inv(A)
            else:
                Ctmp = tt*Cm
                JJ = Jac@Ctmp@Jac.T
                A = JJ + Cd
                r = (obs - cal).T+Jac@diff_m.T
                m_delta = Ctmp*Jac.T@scipy.linalg.solve(A, r)
                # cov a posteriori
                C = Ctmp + Ctmp*Jac.T@scipy.linalg.inv(A)@Jac*Ctmp

            m_test[itest] = m_apri + m_delta
            model_test = insert_mod(M=model_test, m=m_test[itest],
                                            m_act=m_act)
            cali, d_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model_test,
                                                  m_trn=m_trn, m_state=m_state, d_trn=d_trn, d_act=d_act)

            # print(model_test)
            # print(obs)
            # print(cali)

            r_test = Sdi@(obs - cali).T

            dnorm[itest] = scipy.linalg.norm(r_test)
            mnorm[itest] = scipy.linalg.norm(m_test[itest])

            # print(numpy.shape(C))
            m_err[itest] = numpy.sqrt(C.diagonal())

            """
            Calc target function for tau optimization
            """

            if "gcv" in regfun.lower():
                M = Jd@C@Jd.T
                reg_choice[itest] = calc_gcv(dnorm[itest], M)

            elif "mle" in regfun.lower():
                M = Jd@C@Jd.T
                reg_choice[itest] = calc_mle(dnorm[itest], M)

            elif "fix" in regfun.lower() or "lcc" in regfun.lower():
                pass

            else:
                error("Regularisation method "+regfun.lower() +
                      " not yet implemented! Exit.")

        if "fix" in regfun.lower():
            g_index = 0

        elif "lcc" in regfun.lower():
            g_index = calc_lc_corner(dnorm, mnorm)+gshift
            g_index = g_index.item()
            g_index = numpy.amax([g_index, 0])
            g_index = numpy.amin([g_index, numpy.shape(tau)[0]-1])

        elif any(s in regfun.lower() for s in ["gcv", "upr", "ufc", "mle"]):
            g_index = numpy.argmin(reg_choice, axis=0)+gshift
            g_index = g_index.item()
            g_index = numpy.amax([g_index, 0])
            g_index = numpy.amin([g_index, numpy.shape(tau)[0]-1])

        else:
            error(regfun.lower() + " not iplemented!")

        reg = tau[g_index]
        mdl = m_test[g_index, :]
        model = insert_mod(M=model, m=mdl, m_act=m_act)

        """
        Line Search
        """
        if do_linesearch:
            model, dfit_iter = run_linesearch(fwdcall, alt,
                                                                              d_obs=d_obs, d_err=d_err, d_trn=d_trn, d_act=d_act, d_state=d_state,
                                                                              model=model, m_delta=m_delta, m_act=m_act, m_trn=m_trn, m_state=m_state,
                                                                              dfit=dfit_iter, mdfit=thresh[3],
                                                                              facreduce=facreduce, maxreduce=maxreduce, out=OutInfo)
            nrmse_iter, smape_iter = calc_datafit(data_obs=d_obs,
                                                          data_cal=d_cal,
                                                          data_err=d_err,
                                                          data_act=d_act)
        if niter > 0:
            dnorm_iter = numpy.append(dnorm_iter, dnorm[g_index])
            mnorm_iter = numpy.append(mnorm_iter, mnorm[g_index])
            rvals_iter = numpy.append(rvals_iter, reg)
            dfits_iter = numpy.append(dfits_iter, [nrmse_iter, smape_iter])

    # if OutInfo==True:
    print(" %s initial model: dfit =  %7.3f"
          % (profname, dfit_0))
    print(" %s final model:  iter  %6i NRMSE =  %7.3f,  SMAPE = %4.1f percent, %s RegPar is %10.4g "
          % (profname, niter, nrmse_iter, smape_iter, regfun.lower(),  tau[g_index]))
    # print(numpy.shape(sensi))
    sens = sensi  # insert_mod(M=sens, m=sensi,m_act=m_act)
    modl, m_state = transform_parameter(
        m_vec=model, m_trn=m_trn, m_state=m_state,  mode="b")
    modl = extract_mod(modl, m_act)
    merr = m_err[g_index, :].flat
    # insert_mod(M=merr, m=tmp, m_act=m_act)
    # print(numpy.shape(modl),numpy.shape(merr),numpy.shape(sens))

    results = \
        dict([
            ("model", [modl, merr, sens, m_state, m_act]),
            ("data", [d_obs, d_cal, d_err, d_state, d_act]),
            ("log", [niter, conv_status, nrmse_iter, smape_iter,
             dnorm_iter, mnorm_iter, rvals_iter, dfits_iter]),
        ])

    if uncert:
        # Cov a posteriori & Generalized Inverse
        # see: Sun, W. & Durlofsky, L. J.
        # A New Data-Space Inversion Procedure for Efficient Uncertainty
        # Quantification in Subsurface Flow Problems
        # Mathematical Geoscience, 2017, 49, 679-715
        if "par" in invspace.lower():
            JJ = Jac.T@Cdi@Jac
            A = JJ + tau[g_index]*Cmi
            C = scipy.linalg.inv(A)

        else:
            Ctmp = tau[g_index]*Cm
            JJ = Jac@Ctmp@Jac.T
            A = JJ + Cd
            C = Ctmp + Ctmp*Jac.T@scipy.linalg.inv(A)@Jac*Ctmp

        E = numpy.sqrt(C.diagonal())
        G = C@Jd.T

        # Resolution martices & spread
        Rm = G@Jd
        print(numpy.shape(Rm))
        Nm = numpy.sum(Rm.diagonal())
        Sm = scipy.linalg.norm(numpy.identity(numpy.shape(Rm)[0])-Rm)

        Rd = Jd@G
        print(numpy.shape(Rd))
        Nd = numpy.sum(Rd.diagonal())
        Sd = scipy.linalg.norm(numpy.identity(numpy.shape(Rd)[0])-Rd)

        print("Nm =", str(Nm), "Nd =", str(Nd),)

        uncpars =\
            dict([
                ("jacd", Jd),               # jacobian
                ("cpost", C),               # cov a-post
                ("merr", E),                # model error
                ("mresm", [Rm, Sm, Nm]),    # model resolution
                ("dresm", [Rd, Sd, Nd]),    # data resolution
                ("gi", G),                  # generalized inverse
            ])

        results.update(uncpars)

    return results

# @ray.remote


def run_jcn(Ctrl=None, Model=None, Data=None, OutInfo=False):
    """
    Calculate Jackknife variance:

    See:
        Efron, B. & Stein, C.
        The Jackknife Estimate of Variance
        Annals Stat., 1981, 9, 586-596

        Efron, B. The jackknife, the bootstrap and other resampling plans
        SIAM, 1982

        Efron, B. & Hastie, T.
        Computer Age Statistical Inference. Algorithms,
        Evidence, and Data Science
        Cambridge University Press, 2016

    vr October 2022

    """

    """
    unpack contol variables

    Ctrl_jcn = dict([
        ("system", [AEM_system, FwdCall]),
        ("inversion",[RegFun, Tau0, Tau1, Maxiter,thresh, LinPars, SetPrior, delta]),
        ("covar", [L0, Cm0, L1, Cm1]),
        ("transform", [ParaTrans, DataTrans]),
        ("uncert", True)
        ])

    """

    system, fwdcall = Ctrl["system"]
    invtype, regfun, tau0, tau1, maxiter, thresh, linepars, setprior, delta, regshift = Ctrl[
        "inversion"]
    L0, Cm0, L1, Cm1 = Ctrl["covar"]
    jcn_out = Ctrl["output"][0]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    """
    """
    d_act = Data["d_act"]
    jcn_num = numpy.sum(d_act)
    d_err = Data["d_act"]

    if "tikh" in invtype.lower():
        results =\
            run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)

    if "occ" in invtype.lower():
        results =\
            run_tikh_occ(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)

    if "map" in invtype.lower():
        results =\
            run_map(Ctrl=Ctrl, Model=Model, Data=Data,
                    OutInfo=OutInfo)

    jcn_results = results

    ref_mod = results["model"][0]
    ref_rms = results["log"][2]
    ref_act = d_act.reshape((-1, 1))

    """
    loop over leave-out-one samples
    """

    for isample in numpy.arange(jcn_num):

        """
        leave-out-one data set:
        """
        d_sample = copy.deepcopy(d_act)

        if d_act[isample] == 0:
            break
        else:
            d_sample[isample] = 0

        Data["d_act"] = d_sample

        if "tikh" in invtype.lower():

            results =\
                run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                             OutInfo=OutInfo)

        if "occ" in invtype.lower():
            results =\
                run_tikh_occ(Ctrl=Ctrl, Model=Model, Data=Data,
                             OutInfo=OutInfo)

        if "map" in invtype.lower():
            results =\
                run_map(Ctrl=Ctrl, Model=Model, Data=Data,
                        OutInfo=OutInfo)

        model = results["model"][0]
        nrmse = results["log"][2]
        d_calc = results["data"][0]

        if isample == 0:
            w_jcn = numpy.array([1./nrmse])
            m_jcn = model.reshape((-1, 1))
            a_jcn = d_act.reshape((-1, 1))
        else:
            w_jcn = numpy.append(w_jcn, numpy.array([1./nrmse]))
            m_jcn = numpy.append(m_jcn, model.reshape((-1, 1)), axis=1)
            a_jcn = d_act.reshape((-1, 1))

    w_jcn = w_jcn/numpy.sum(w_jcn)

    """
    jackknife average & variance
    """

    jcn_avg = numpy.mean(m_jcn, axis=1).reshape((-1, 1))
    fac = (jcn_num-1)/jcn_num
    jcn_var = fac * numpy.sum((m_jcn - jcn_avg)**2)

    jcn_med = numpy.median(m_jcn, axis=1).reshape((-1, 1))
    d = numpy.abs(m_jcn - jcn_med)
    jcn_mad = numpy.nanmedian(d)

    jcnfields = (
        ("ref_mod", ref_mod),
        ("ref_act", ref_act),
        ("ref_rms", ref_rms),
        ("jcn_avg", jcn_avg),
        ("jcn_var", jcn_var),
        ("jcn_med", jcn_med),
        ("jcn_mad", jcn_mad),
    )
    jcn_results.update(jcnfields)

    if "ens" in jcn_out:
        jcn_ens = numpy.vstack((m_jcn, w_jcn))
        jcn_results["jcn_ens"] = jcn_ens

    return jcn_results


def run_rto(Ctrl=None, Model=None, Data=None, OutInfo=False):
    """
    Run the randomize-then-optimize (RTO) algorithm:

        for i = 1 : nsamples do
            Draw perturbed data set: d_pert∼ N (d, Cd)
            Draw prior model: m̃ ∼ N (0, 1/mu (LT L)^−1 )
            Solve determistic problem  to get the model m_i
        end

    See:

    Bardsley, J. M.; Solonen, A.; Haario, H. & Laine, M.
        Randomize-Then-Optimize: a Method for Sampling from Posterior
        Distributions in Nonlinear Inverse Problems
        SIAM J. Sci. Comp., 2014, 36, A1895-A1910

    Blatter, D.; Morzfeld, M.; Key, K. & Constable, S.
        Uncertainty quantification for regularized inversion of electromagnetic
        geophysical data. Part I: Motivation and Theory
        Geophysical Journal International, doi:10.1093/gji/ggac241, 2022

    Blatter, D.; Morzfeld, M.; Key, K. & Constable, S.
        Uncertainty quantification for regularized inversion of electromagnetic
        geophysical data – Part II: application in 1-D and 2-D problems
        Geophysical Journal International, , doi:10.1093/gji/ggac242, 2022

    vr July 2022

    """

    """
    unpack control variables and run reference model

    Ctrl_RTO = dict([
        ("system", [AEM_system, FwdCall]),
        ("rto", [nsamples, Cdpost, mu])
        ("inversion",[InvType, RegFun, Tau0, Tau1, Maxiter,thresh, LinPars, SetPrior, delta]),
        ("covar", [L0, Cm0, L1, Cm1]),
        ("transform", [ParaTrans, DataTrans]),
        ("uncert", True)
        ])
    Ctrl["rto"] = [nsamples, Percentiles ]
    Ctrl["output"] = ["ens "]

    """
    system, fwdcall = Ctrl["system"]
    invtype, regfun, tau0, tau1, maxiter, thresh, linepars, setprior, delta, regshift = Ctrl[
        "inversion"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    L0, Cm0, L1, Cm1 = Ctrl["covar"]
    nsamples, Percentiles = Ctrl["rto"]

    rto_out = Ctrl["output"][0]

    """
    run reference model
    """
    if "opt" in invtype.lower():

        results =\
            run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)
    if "map" in invtype.lower():
        results =\
            run_map(Ctrl=Ctrl, Model=Model, Data=Data,
                    OutInfo=OutInfo)
    # print(results.keys())
    rto_results = results

    """
    unpack results:
            ("model", [modl, merr, sens, m_state, m_act]),
            ("data", [d_obs, d_cal, d_err, d_state, d_act]),
            ("log", [niter, conv_status, nrmse_iter, smape_iter,
                     dnorm_iter, mnorm_iter, rvals_iter, dfits_iter]),
    """
    d_obs = results["data"][0]  # .reshape(-1,1)
    d_err = results["data"][2]
    d_act = results["data"][4]

    """
    unpack model
    """
    m_act = Model["m_act"]
    m_bas = Model["m_apr"].copy()
    m_opt = results["model"][0]

    m_ref = insert_mod(M=m_bas, m_act=m_act, m=m_opt)
    c_ref = results["cpost"]

    """
    Draw perturbed data set: d  ̃ ∼ N (d, Cd)
    """
    d_ens = generate_data_ensemble(dref=d_obs, dact=d_act,
                                           nens=nsamples,
                                           perturb=["gauss", d_err],
                                           out=OutInfo)
    # print(d_ens)
    # print("@@@")
    # print(d_err)
    # print("@@@")
    """
    Draw prior model: m̃ ∼ N (0, 1 (LT L)−1 )
    """
    m_ens = generate_param_ensemble(mref=m_ref, mact=m_act,
                                            nens=nsamples,
                                            perturb=["gauss", c_ref,
                                                     numpy.array([])],
                                            out=OutInfo)
    """
    Solve inverse problem for rto_ens
    """
    rto_ens = m_ref.copy()

    for isample in numpy.arange(nsamples):
        Data["d_obs"] = d_ens[isample, :]
        # print(isample, numpy.shape(m_ref))
        Model["m_apr"] = insert_mod(M=m_ref.copy(), m_act=m_act,
                                            m=m_ens[isample, :])

        if "opt" in invtype.lower():

            results =\
                run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                             OutInfo=OutInfo)

        m = insert_mod(M=m_ref, m_act=m_act, m=results["model"][0])
        if isample == 0:
            rto_ens = m
        else:
            rto_ens = numpy.vstack((rto_ens, m))

    ne = numpy.shape(rto_ens)
    rto_avg = numpy.mean(rto_ens, axis=1)
    # rto_std = numpy.std(rto_ens, axis=1)
    rto_var = numpy.var(rto_ens, axis=1)
    rto_med = numpy.median(rto_ens, axis=1)
    # print(numpy.shape(rto_ens), numpy.shape(rto_med))
    # print(ne)
    mm = numpy.tile(rto_med, (ne[1], 1))
    # print(numpy.shape(mm))

    rto_mad = numpy.median(
        numpy.abs(rto_ens.T - numpy.tile(rto_med, (ne[1], 1))))

    rto_prc = numpy.percentile(rto_ens, Percentiles)

    rtofields = (
        ("rto_avg", rto_avg),
        ("rto_var", rto_var),
        ("rto_med", rto_med),
        ("rto_mad", rto_mad),
        ("rto_prc", rto_prc)
    )
    rto_results.update(rtofields)

    if "ens" in rto_out:
        rto_results["jcn_ens"] = rto_ens

    return rto_results


# @ray.remote
def run_EnKI(Ctrl=None, Model=None, Data=None, OutInfo=True):
    """
    Ensemble Kalman inversion
    M. A. Iglesias, K. J. H. Law, and A. M. Stuart,
        “The ensemble Kalman filter for inverse problems”,
        Inverse Problems, 29, 2013,
        doi:10. 1088/0266-5611/29/4/045001.

    C.-H. M. Tso, M. Iglesias, P. Wilkinson, O. Kuras, J. Chambers, A. Binley,
        “Efficient multiscale imaging of subsurface resistivity with
        uncertainty quantification using ensemble Kalman inversion”,
        Geophysical Journal International, 225, 2021
        doi: 10.1093/gji/ggab013.

    M. Iglesias, D. M. McGrath, M. V. Tretyakov, and Susan T Francis,xxx
        “Ensemble Kalman inversion for magnetic resonance elastography” ,
        Phys. Med. Biol., vol. 67, p. 235003, 2022, doi: 10.1088/1361-6560/ac9fa1.

    C.-H. M. Tso and M. Iglesias and A. Binley
        “Ensemble Kalman inversion of induced polarization data”
        Geophysical Journal International, 236, 2024,
        doi: 10.1093/gji/ggae012.


    M. Y. Matveev, A. Endruweit, A. C. Long, M. A. Iglesias, and M. V. Tretyakov,
        “Bayesian inversion algorithm for estimating local variations
        in permeability and porosity of reinforcements using experimental data”,
        Composites Part A: Applied Science and Manufacturing, 143, 106323, 2021,
        doi: 10.1016/j.compositesa.2021.106323.

    S. Lan, S. Li, and M. Pasha
        "Bayesian spatiotemporal modeling for inverse problems",
        Stat Comput 33, 89 (2023). https://doi.org/10.1007/s11222-023-10253-z


    Created on Jan 17, 2022

    @author: vrath

    """

    """
    unpack control block:

        Ctrl_ENK= dict([("aemsys", AEM_system),
                    ("system", [AEM_system, FwdCall]),
                    ("transform", [ParaTrans, DataTrans]),
                    ("eki", nsamples, Maxiter, Tau, Cm, Cd, percentiles, ens_out)
                    ])

    """
    system, fwdcall = Ctrl["system"]
    invtype, regfun, tau0, tau1, maxiter, thresh, linepars, setprior, delta, regshift = Ctrl[
        "inversion"]
    L0, Cm0, L1, Cm1 = Ctrl["covar"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    nsamples, percentiles = Ctrl["enk"]

    ens_out = Ctrl["output"][0]
    # print(Cm1)

    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_act = Data["d_act"]  # .reshape(-1,1)
    d_obs = Data["d_obs"]  # .reshape(-1,1)
    d_err = Data["d_err"]  # .reshape(-1,1)
    alt = Data["alt"]
    n_data = round(numpy.sum(d_act[d_act != 0]/d_act[d_act != 0]))

    print("N_samples=", nsamples, ",   N_data=", n_data)
    d_cal = 0. * numpy.ones((nsamples, n_data))
    d_res = 0. * numpy.ones((nsamples, n_data))

    d_state = 0
    d_obs, d_err, dobs_state = transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)

    Cd = numpy.diag(d_err**2, 0)
    CdiSq = numpy.diag(1./d_err)

    """
    unpack model block

    Model =\
        [model_act, model_prior, model_var, model_bounds, model_ini]
    """

    m_act = Model["m_act"]
    m_ref = Model["m_apr"]
    # print(" mref", m_ref[:36])
    m_state = 0
    m_ref, m_state = transform_parameter(m_vec=m_ref,
                                                 m_trn=m_trn, m_state=m_state, mode="f")

    # print(m_state, " mref", m_ref[:36])
    """
    Draw prior model: m_p∼ N (0, 1 (LT L)−1 )
    """
    m_ens = generate_param_ensemble(mref=m_ref, mact=m_act,
                                            nens=nsamples,
                                            perturb=["gauss", Cm1],
                                            out=OutInfo)
    # """
    # Draw perturbed data set: d  ̃ ∼ N (d, Cd)
    # """
    # d_ens = generate_data_ensemble(dref=d_obs, dact = d_act,
    #                                       nens=nsamples,
    #                                       perturb=["gauss" ,0.,Cd],
    #                                       out=OutInfo)

    """
    Solve inverse problem for ensemble

    """
    s = 0.
    for niter in numpy.arange(maxiter):
        """
        Prediction
        """

        for isample in numpy.arange(nsamples):
            # print("XXX ",numpy.shape(m_ref))
            # print("XXX ",type(m_ref))
            model = m_ref.copy()
            m_smp = m_ens[isample, :]
            # print(isample, " m_ens ", m_smp)
            # print("model", model[m_act!=0])
            print(m_state)
            m = insert_mod(model, m_smp, m_act)
            d_cal[isample, :], _ = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                         m_vec=m, m_trn=m_trn, m_state=m_state,
                                                         d_act=d_act, d_trn=d_trn, d_state=d_state, out=True)
            print(" dobs ", d_obs)
            print(" dcal ", d_cal[isample, :])
            d_res[isample, :] = CdiSq@(d_obs[:]-d_cal[isample, :])
            print(" dres", d_res[isample, :])
        """
        Analysis
        """
        a = numpy.sum(d_res**2)/(n_data*nsamples)
        s, alpha = update_reg(s=s, alpha_ast=a)
        print("s, a,  alpha:  ",  s, a, alpha)

        Gyy = calc_encovar(x=d_cal, y=d_cal, method=0, out=True)
        Gxy = calc_encovar(x=m_ens, y=d_cal, method=0, out=True)

        # print("Gyy ",Gyy)

        """
        Ensemble update
        """

        sqalph = numpy.sqrt(alpha)

        for isample in numpy.arange(nsamples):

            p = sqalph*Cd @ rng.standard_normal(numpy.shape(d_err))
            rhs = d_obs - d_cal[isample, :] + p
            sol = numpy.linalg.solve(Gyy + alpha*Cd, rhs)
            # print("pert  ",numpy.shape(p),numpy.shape(pert))
            # print("rhs  ",numpy.shape(rhs))
            # print(numpy.shape(Gxy))
            # print(numpy.shape(sol))

            m_ens[isample, :] = m_ens[isample, :] + Gxy@sol

        if s >= 1.:
            break

    mod_avg = numpy.mean(m_ens, axis=1)
    mod_std = numpy.std(m_ens, axis=1)
    mod_med = numpy.percentile(m_ens, 50.)
    mod_prc = numpy.percentile(m_ens, percentiles)

    ekiresults =\
        dict([
            ("avg", mod_avg),
            ("std", mod_std),
            ("med", mod_med),
            ("percentiles", mod_prc),
        ])
    if ens_out:
        ekiresults["ens"] = m_ens

    return ekiresults


def update_reg(s=None, alpha_ast=None):
    """
    Update regularisation parameter in EKI

    Parameters
    ----------
    s : float, optional
        paramter s_n in EKI. The default is 0..
    ens_norm : float, optional
        DESCRIPTION. The default is 0..

    Returns
    -------
    s : float
        paramter s_n+1 in EKI.
    alpha : float
        regularisation paramter updated.

    M. Iglesias, D. M. McGrath, M. V. Tretyakov, and Susan T Francis,
        “Ensemble Kalman inversion for magnetic resonance elastography” ,
        Phys. Med. Biol., vol. 67, p. 235003, 2022, doi: 10.1088/1361-6560/ac9fa1.


    """

    alphast = alpha_ast

    if (s+1./alphast) >= 1.:
        alpha = 1./(1-s)
        s = 1.
    else:
        alpha = alphast
        s = s + 1./alpha

    return s, alpha


#  @ray.remote
def run_nullspace(Ctrl=None, Model=None, Data=None, OutInfo=True):
    """
    nullspace shuttle

    Parameters
    ----------

    OutInfo : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    nssresults

    Referencea:

    Deal M, Nolet G (1996)
    Nullspace shuttles
    Geophys. J. Int.,124, 372-380


    Muñoz G, Rath V (2006)
    Beyond smooth inversion: the use of nullspace projection for
    the exploration of non-uniqueness in MT
    Geophys. J. Int.,164, 301-311

    """
    system, fwdcall = Ctrl["system"]
    invtype, regfun, tau0, tau1, maxiter, thresh, linepars, setprior, delta, regshift = Ctrl[
        "inversion"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    L0, Cm0, L1, Cm1 = Ctrl["covar"]
    nsamples, percentiles, k, randsvd = Ctrl["nss"]
    ens_out = Ctrl["output"][0]

    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_act = Data["d_act"]  # .reshape(-1,1)
    d_obs = Data["d_obs"]  # .reshape(-1,1)
    d_err = Data["d_err"]  # .reshape(-1,1)
    alt = Data["alt"]

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    d_state = 0
    d_obs, d_err, dobs_state = transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)
    """
    unpack model block

    Model =\
    [model_act, model_prior, model_var, model_bounds, model_ini]
    unpack model
    """
    m_act = Model["m_act"]
    m_bas = Model["m_apr"].copy()

    """
    run reference model
    """
    if "opt" in invtype.lower():
        results =\
            run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)
    if "occ" in invtype.lower():
        results =\
            run_tikh_occ(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)

    if "map" in invtype.lower():
        results =\
            run_map(Ctrl=Ctrl, Model=Model, Data=Data,
                    OutInfo=OutInfo)

    m_opt = results["model"][0]
    m_ref = insert_mod(M=m_bas, m_act=m_act, m=m_opt)
    c_ref = results["cpost"]

    """
    now calculate the SVD of the Jacobian
    """
    Jacd = results["jacd"]

    if randsvd:
        U, S, Vt = rsvd(
            Jacd, rank=k, n_oversamples=0, n_subspace_iters=2)
    else:
        U, S, Vt = scipy.linalg.svd(Jacd, full_matrices=False)

    """
    truncation
    """

    V = Vt.T
    V = V[:, :k]
    S = S[:k]
    U = U[:, :k]

    """
    check  how much of Jacobian is explained by k
    """
    D = U@scipy.sparse.diags(S[:])@Vt - Jacd
    x_op = numpy.random.default_rng().normal(size=numpy.shape(D)[1])
    n_op = numpy.linalg.norm(D@x_op)/numpy.linalg.norm(x_op)
    j_op = numpy.linalg.norm(Jacd@x_op)/numpy.linalg.norm(x_op)
    if OutInfo:
        print(" Op-norm J_k = "+str(n_op)+", explains "
              + str(100. - n_op*100./j_op)+"% of variations")

    """
    Draw prior model: m̃ ∼ N (0, 1 (LT L)−1 )
    """
    m_ens = generate_param_ensemble(mref=m_ref, mact=m_act,
                                            nens=nsamples,
                                            perturb=["gauss", c_ref,
                                                     numpy.array([])],
                                            out=OutInfo)

    m_prj = numpy.zeros_like(m_ens)
    for isample in numpy.arange(nsamples):
        m_prj[isample, :] = project_nullspace(
            U=U, m_test=m_ens[isample, :])

    nss_avg = numpy.mean(m_ens, axis=1)
    nss_std = numpy.std(m_ens, axis=1)
    nss_med = numpy.percentile(m_ens, 50.)
    nss_prc = numpy.percentile(m_ens, percentiles)

    nss_results = results
    nss_results["nss_avg"] = nss_avg
    nss_results["nss_std"] = nss_std
    nss_results["nss_med"] = nss_med
    nss_results["nss_percentiles"] = nss_prc
    if ens_out:
        nss_results["ens"] = m_ens
        nss_results["prj"] = m_prj
    return nss_results


def run_sample_pcovar(Ctrl=None, Model=None, Data=None, OutInfo=True):
    """
    Algorithm given by  Osypov (2013)

    Parameters
    ----------

    OutInfo : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    spcresults

    References:

    Osypov K, Yang Y, Fournier A, Ivanova N, Bachrach R,
    Can EY, You Y, Nichols D, Woodward M (2013)
    Model-uncertainty quantification in seismic tomography: method and applications
    Geophysical Prospecting, 61, pp. 1114–1134, 2013, doi: 10.1111/1365-2478.12058.


    """
    system, fwdcall = Ctrl["system"]
    invtype, regfun, tau0, tau1, maxiter, thresh, linepars, setprior, delta, regshift = Ctrl[
        "inversion"]

    d_trn = Ctrl["data"][0]
    m_trn = Ctrl["model"][0]

    L0, Cm0, L1, Cm1 = Ctrl["covar"]
    nsamples, perc, k, randsvd = Ctrl["nss"]
    nss_out = Ctrl["output"][0]

    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_act = Data["d_act"]  # .reshape(-1,1)
    d_obs = Data["d_obs"]  # .reshape(-1,1)
    d_err = Data["d_err"]  # .reshape(-1,1)
    alt = Data["alt"]

    d_cal = numpy.nan * numpy.ones_like(d_obs)

    d_state = 0
    d_obs, d_err, dobs_state = transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)
    """
    unpack model block

    Model =\
    [model_act, model_prior, model_var, model_bounds, model_ini]
    """

    m_act = Model["m_act"]
    m_ref = Model["m_apr"]

    m_state = 0
    m_ref, m_state = transform_parameter(
        m_vec=m_ref, m_trn=m_trn, m_state=m_state, mode="f")
    """
    run reference model
    """
    if "opt" in invtype.lower():
        results =\
            run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)
    if "occ" in invtype.lower():
        results =\
            run_tikh_occ(Ctrl=Ctrl, Model=Model, Data=Data,
                         OutInfo=OutInfo)

    if "map" in invtype.lower():
        results =\
            run_map(Ctrl=Ctrl, Model=Model, Data=Data,
                    OutInfo=OutInfo)

    nss_results = results

    """
    now calculate the SVD of the Jacobian
    """
    Jacd = results["jacd"]

    if randsvd:
        U, S, Vt = rsvd(
            Jacd.T, rank=k, n_oversamples=2, n_subspace_iters=2)
    else:
        U, S, Vt = scipy.linsvd(Jacd.T, full_matrices=False)

    """
    truncation
    """

    V = Vt.T
    V = V[:, :k]
    S = S[:k]
    U = U[:, :k]

    """
    chek  how mauch of Jacd is explained by k
    """
    D = U@scipy.sparse.diags(S[:])@Vt - Jacd.T
    x_op = numpy.random.default_rng().normal(size=numpy.shape(D)[1])
    n_op = numpy.linalg.norm(D@x_op)/numpy.linalg.norm(x_op)
    j_op = numpy.linalg.norm(Jacd.T@x_op)/numpy.linalg.norm(x_op)
    if OutInfo:
        print(" Op-norm J_k = "+str(n_op)+", explains "
              + str(100. - n_op*100./j_op)+"% of variations")

    print("This algorithm is not yet implemented! Exit.")
    return nss_results

# @ray.remote


def run_most_squares(OutInfo=True):
    """
    Most Squares iteration

    Parameters
    ----------

    OutInfo : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    msqresults

    Reference:

    Meju MA, Hutton VRS (1992)
    Iterative most-squares inversion: application to magnetotelluric data
    Geophys. J. Int., 108, 758-766

    Meju MA (2009)
    Regularized extremal bounds analysis (REBA): Anapproach to
    quantifying uncertainty in nonlineargeophysical inverse problems
    Geophyical Research Letters, 36, doi:10.1029/2008GL036407

    Kalscheuer T,  De Los Angeles Garcia Juanatey M, Meqbel N,
    Pedersen, LB (2010)
    Non-linear model error and resolution properties from two-dimensional
    single and joint inversions of direct current resistivity and
    radiomagnetotelluric data, Geophysical Journal International,
    182(3), doi:10.1111/j.1365-246X.2010.04686.x

    Ren Z  & Kalscheuer T (2019)
    Uncertainty and Resolution Analysis of 2D and 3D Inversion Models
    Computed from Geophysical Electromagnetic Data
    Surveys in Geophysics, 1-66

    """
    msqresults = dict([])
    error("This algorithm is not yet implemented! Exit.")
    return msqresults



def calc_fwdmodel(fwdcall=None,
                  alt=None,
                  m_vec=numpy.array([]),
                  m_trn=0,
                  m_state=0,
                  m_params=numpy.array([]),
                  d_act=numpy.array([]),
                  d_trn=0,
                  d_state=0,
                  out=True):
    """
    Calculate forward model.

    last changes: VR 6/22
    """
    if numpy.size(m_vec) == 0:
        error("calc_fwdmodel: No model given! Exit.")

    nlyr = get_nlyr(m_vec)
    m, _ = transform_parameter(
        m_vec=m_vec, m_trn=m_trn, m_state=m_state, mode="b")

    d_vec = eval(fwdcall)

    if d_trn == 0:
        return d_vec, d_state
    else:
        d_vec, _, d_state = transform_data(d_vec=d_vec, d_trn=d_trn)

    if numpy.size(d_act) == 0:
        return d_vec, d_state
    else:
        d_vec = extract_dat(D=d_vec, d_act=d_act)

    return d_vec, d_state


# def run_linesearch_nested(fwdcall, alt,
#                    d_obs=numpy.array([]), d_err=numpy.array([]),
#                    d_trn=0, d_act=numpy.array([]), d_state=0,
#                    model=numpy.array([]), m_delta=numpy.array([]),
#                    m_act=numpy.array([]), m_trn=1, m_state=0,
#                    dfit=[999999.9,999999.9], numdiv = 5,
#                    maxreduce=3, mdfit="rms", out=True):
#     """
#     Run simple line search with 0.>alpha<1.

#     Parameters
#     ----------
#     fwdcall : str
#         Forward call.
#     alt : TYPE
#         Altitude.
#     d_obs, d_cal, d_err : float
#         Date (obs, cal) and errors.
#     d_trn , d_act : integer
#         Data transformation.
#     model : float nd arrray.
#         Base model (full). The default is numpy.array([]).
#     m_delta : TYPE, optional
#         Model update.
#     m_act, m_trn : integer.
#           Model transform.
#     mdfit : string
#         Controls which metric is usedL "rms" or "smp"
#     dfit : float
#         Initial dfit.
#     facreduce : float, optional
#         Reduction factor. The default is 0.6666.
#     maxreduce : integer, optional
#         maximal number of reduction steps. The default is 6.
#     out : logical, optional
#         Extended output. The default is False.

#     Returns
#     -------
#     model: float
#     linfit: float
#     linrms: float
#     linsmp: float
#     """

#     # print(dfit)

#     liniter = 0
#     linfit = dfit
#     linfit_old = dfit

#     mbase = model.copy()
#     model_old =  model.copy()
#     mact = extract_mod(M=mbase, m_act=m_act)

#     test_points = numpy.linspace(0.05, 0.95, numdiv )
#     while (liniter < maxreduce) and (linfit <= dfit):
#         liniter = liniter + 1
#         mfull = mbase.copy()
#         linfit_iter = []
#         for ptest in test_points:
#             mtest = mact - ptest*m_delta.reshape(numpy.size(mact))
#             m = insert_mod(M=mfull, m=mtest, m_act=m_act)
#             d_calc, d_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
#                                             m_vec=m, m_trn=m_trn, m_state=m_state,
#                                             d_trn=d_trn, d_state=d_state)
#             linrms_iter, linsmp_iter = calc_datafit(data_cal=d_calc,
#                                                  data_obs=d_obs,
#                                                  data_err=d_err,
#                                                  data_act=d_act)
#             if "rms" in mdfit: linfit_iter.append(linrms_iter)
#             if "smp" in mdfit: linfit_iter.append(linsmp_iter)

#         valfit_min = min(linfit_iter)
#         indfit_min = linfit_iter.index(valfit_min)



#         newmin = max(indfit_min-1, 0)
#         newmax = min(indfit_min+1, numdiv)
#         test_points =numpy.linspace(test_points[newmin], test_points[newmax], numdiv)



#         if linfit_iter < linfit:
#             fact = fact * facreduce
#             linfit = linfit_iter
#             linfit_old = linfit_iter
#             model = m
#             model_old = m


#         else:
#             linfit = linfit_old
#             model = model_old
#             break

#         if out:
#             print("Linesearch: ", liniter,"/",linfit, "/",linfit_iter)

#     return model, linfit

def run_linesearch(fwdcall, alt,
                   d_obs=numpy.array([]), d_err=numpy.array([]),
                   d_trn=0, d_act=numpy.array([]), d_state=0,
                   model=numpy.array([]), m_delta=numpy.array([]),
                   m_act=numpy.array([]), m_trn=1, m_state=0,
                   dfit=999999.9,
                   facreduce=0.6666, maxreduce=6, mdfit="rms",
                   out=False):
    """
    Run simple line search with 0.>alpha<1.

    Parameters
    ----------
    fwdcall : str
        Forward call.
    alt : TYPE
        Altitude.
    d_obs, d_cal, d_err : float
        Date (obs, cal) and errors.
    d_trn , d_act : integer
        Data transformation.
    model : float nd arrray.
        Base model (full). The default is numpy.array([]).
    m_delta : TYPE, optional
        Model update.
    m_act, m_trn : integer.
          Model transform.
    mdfit : string
        Controls which metric is usedL "rms" or "smp"
    dfit : float
        Initial dfit.
    facreduce : float, optional
        Reduction factor. The default is 0.6666.
    maxreduce : integer, optional
        maximal number of reduction steps. The default is 6.
    out : logical, optional
        Extended output. The default is False.

    Returns
    -------
    model: float
    linfit: float
    linrms: float
    linsmp: float
    """

    # print(dfit)

    liniter = 0
    linfit = dfit
    linfit_old = dfit

    fact = facreduce
    mbase = model.copy()
    model_old =  model.copy()
    mact = extract_mod(M=mbase, m_act=m_act)
    while (liniter < maxreduce) and (linfit <= dfit):
        liniter = liniter + 1
        mfull = mbase.copy()
        mtest = mact - fact*m_delta.reshape(numpy.size(mact))
        m = insert_mod(M=mfull, m=mtest, m_act=m_act)
        d_calc, d_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                        m_vec=m, m_trn=m_trn, m_state=m_state,
                                        d_trn=d_trn, d_state=d_state)
        linrms_iter, linsmp_iter = calc_datafit(data_cal=d_calc,
                                             data_obs=d_obs,
                                             data_err=d_err,
                                             data_act=d_act)
        if "rms" in mdfit: linfit_iter =  linrms_iter
        if "smp" in mdfit: linfit_iter =  linsmp_iter
        # if out:
        #     print("Linesearch: ", liniter,"/",linfit, "/",linfit_iter)

        if linfit_iter < linfit:
            fact = fact * facreduce
            linfit = linfit_iter
            linfit_old = linfit_iter
            model = m
            model_old = m


        else:
            linfit = linfit_old
            model = model_old
            break

        if out:
            print("Linesearch: ", liniter,"/",linfit, "/",linfit_iter)

    return model, linfit


def perturb_random(vbase=numpy.array([]),
                   n_ens=128,
                   std=1., avg=0.,
                   covar=numpy.array([]),
                   base_add=False,
                   out=True):
    """
    Random perturbation of vector.

    Default draws perturbation from Gaussian.

    Created on Sept 26,2023

    @author: vrath

    """
    if vbase.size == 0:
        error("generate_ensemble: No Reference model given! Exit.")

    l = 1.
    if covar.ndim == 2:
         if scipy.sparse.issparse(covar):
             l = msqrt_sparse(covar)
         else:
             l = msqrt(covar)
        # l = scipy.linalg.cholesky(icovar)


    nd = numpy.shape(vbase)

    for iens in numpy.arange(n_ens):
        p_normal = std*rng.standard_normal(nd)

        pert = p_normal*l

        if base_add:
            pert = avg + pert

        if iens == 0:
            p_ens = pert

        else:
            p_ens = numpy.vstack((p_ens, pert))

    return p_ens


def generate_data_ensemble(dref=numpy.array([]),
                           dact=numpy.array([]),
                           dtrn=0,
                           nens=128,
                           perturb=["Gauss", 1.],
                           inchol=numpy.array([]),
                           incovar=numpy.array([]),
                           out=True):
    """
    Generate Data Ensemble for Kalman methods

    Default draws perturbation from Gaussian.

    Created on Jul 8, 2022

    @author: vrath

    """
    if dref.size == 0:
        error("generate_data_ensemble: No Reference model given! Exit.")

    ndat = numpy.shape(dref[dact != 0])[0]
    if ndat == 0:
        error("generate_data_ensemble: no active data given! exit.")


    dbas, _, _ = transform_data(d_vec= dref.copy(), d_trn=dtrn, d_state=0)
    dbas = dref[dact != 0]
    dshp = numpy.shape(dbas)



    if incovar.size == 0:
        for iens in numpy.arange(nens):

            if iens == 0:
                dens = dbas.copy() + perturb[1]*rng.standard_normal(dshp)
                # print(dshp)
                # print(perturb[1])
                # print(rng.standard_normal(dshp))
            else:
                dens = numpy.vstack((dens,
                                    dbas.copy() + perturb[1]*rng.standard_normal(dshp)))
            # print(dens)
    else:
        if inchol.size == 0:
            if scipy.sparse.issparse(incovar):
                l = msqrt_sparse(incovar)
            else:
                l = msqrt(incovar)
           # l = scipy.linalg.cholesky(incovar)
        else:
            l = inchol

        for iens in numpy.arange(nens):
            if iens == 0:
                dens = dbas.copy() + l @ rng.standard_normal(dshp)
            else:
                dens = numpy.vstack(
                    (dens, dbas.copy() + l @ rng.standard_normal(dshp)))

            # print(dens)

    return dens


def generate_param_ensemble(mref=numpy.array([]),
                            mact=numpy.array([]),
                            mtrn=1, mstate=1,
                            nens=128,
                            perturb=["Gauss", 1.],
                            out=True):
    """
    Generate model Ensemble for Kalman-type methods

    Default draws perturbation from Gaussian.

    Created on Jan 18, 2022

    @author: vrath

    """
    if mref.size == 0:
        error("generate_param_ensemble: no reference model given! exit.")

    mpar = numpy.shape(mref[mact != 0])[0]
    if mpar == 0:
        error("generate_param_ensemble: no active parameter given! exit.")


    meth = perturb[0]
    covm = perturb[1]
    chol = numpy.array([])
    if len(perturb)==3: chol = perturb[2]


    # print("mbas")
    # print(mref[mact != 0])
    mbas, mtrn = transform_parameter(m_vec= mref.copy(), mode="f",
                                     m_trn=mtrn, m_state=mstate)

    mwrk = mbas[mact != 0]
    mshp = numpy.shape(mwrk)

    # print(mshp)
    # print("mbas", mbas)
    # print(numpy.exp(mbas))

    if "gau" in meth.lower():
        print("generate_param_ensemble: normal perturbation assumed.")
        covm = perturb[1]
        chol = numpy.array([])
        if len(perturb)==3: chol = perturb[2]
        if chol.size != 0:
            for iens in numpy.arange(nens):
                p = perturb[2]*rng.standard_normal(mshp)
                mtmp = mwrk.copy() + p
                # print(mtmp[1:5])
                mtmp = insert_mod(M=mbas.copy(), m=mtmp, m_act=mact)[mact != 0]
                if iens == 0:
                    mens = mtmp
                    # print(numpy.shape(mens))
                else:
                    mens = numpy.vstack((mens, mtmp))
                    # print(numpy.shape(mens))
                # print(mtmp[0:5])

        else:
            # print(numpy.shape(covm))
            if scipy.sparse.issparse(covm):
                l = msqrt_sparse(covm)
            else:
                l = msqrt(covm)
               # l = scipy.linalg.cholesky(covm)

            # print(inchol.size)
            # print(" mwrk", mwrk[0:7])

            for iens in numpy.arange(nens):
                mtmp = mwrk.copy() + l @ rng.standard_normal(mshp)
                # print(mtmp[1:5])
                mtmp = insert_mod(M=mbas.copy(), m=mtmp, m_act=mact)[mact != 0]
                # mtmp, _ = transform_parameter(m_vec= mtmp, mode="b", m_trn=mtrn, m_state=1)
                if iens == 0:
                    mens = mtmp
                else:
                    mens = numpy.vstack((mens, mtmp))

                print(mtmp[0:7])

    else:
        print("generate_param_ensemble: uniform perturbation assumed. ")
        low, high = perturb[1]
        for iens in numpy.arange(nens):
            mtmp = mwrk.copy() + rng.uniform(low=low, high=high, size=mshp)
            # print(mtmp[1:5])
            mtmp = insert_mod(M=mbas.copy(), m=mtmp, m_act=mact)[mact != 0]
            # mtmp, _ = transform_parameter(m_vec= mtmp, mode="b", m_trn=mtrn, m_state=1)
            if iens == 0:
                mens = mtmp
            else:
                mens = numpy.vstack((mens, mtmp))

    return mens


def calc_encovar(x=numpy.array([]),
                y=numpy.array([]),
                cscale=numpy.array([]),
                method=0, out=True):
    """
    Calculate empirical covariance for Kalman gain

    Created on Jul 6, 2022

    @author: vrath


    """

    if (x.size == 0) and (y.size == 0):
        error("calc_encovar: No data given! Exit.")

    X = x - numpy.mean(x, axis=0)
    if (y.size == 0):
        Y = X
    else:
        Y = y - numpy.mean(y, axis=0)

    [N_e, N_x] = numpy.shape(X)
    [N_e, N_y] = numpy.shape(Y)

    if method == 0:
        # print(N_e, N_x, N_y)
        # naive version, library versions probably faster)
        C = numpy.zeros((N_x, N_y))
        for n in numpy.arange(N_e):
            # print("XT  ",X.T)
            # print("Y   ",Y)
            Cn = X.T@Y
            # print(Cn)
            C = C + Cn

        C = C/(N_e-1)

    else:
        # numpy version
        for n in numpy.arange(N_e):
            X = numpy.stack((X, Y), axis=0)
            # C = numpy.cov((X,Y))
            C = numpy.cov((X))

    if out:
        print("Ensemble covariance is "+str(numpy.shape(C)))

    return C


def calc_keg_update(m=numpy.array([]), r=numpy.array([]),
                    Cxy=numpy.array([]),
                    Cyy=numpy.array([]),
                    Cd=numpy.array([]),
                    out=True):
    """
    Calculate Kalman update

    see:
        C. Bobe
        Efficient probabilistic processing of frequency-
        domain electromagnetic data for subsurface modelling,
        PhD thesis, Ghent University, 2020.

        Nowak, W.
        Best unbiased ensemble linearization and the
        quasi-linear Kalman ensemble generator
        Water Resour. Res., 45, W04431, doi:10.1029/2008WR007328, 2009.

    Created on Jul 9, 2022

    @author: vrath

    """
    if (m.size == 0):
        error("calc_kupdate: No model given! Exit.")
    if (r.size == 0):
        error("calc_kupdate: No residuals given! Exit.")

    ms = numpy.shape(m)
    mr = numpy.shape(r)
    if mr[0] != ms[0]:
        error("calc_kupdate: Ensemble sizes don't match! Exit.")

    print(ms, mr)
    if (Cxy.size == 0) or (Cyy.size == 0):
        error("calc_keg_update: One or more covariances missing! Exit.")
    if (Cd.size == 0):
        error("calc_keg_update: Data covariance missing! Exit.")
    else:
        C = C = numpy.eye(mr[1])/Cd

    Ci = scipy.linalg.inv(Cyy+C)
    K = Cxy@Ci

    m_update = m.copy()
    for ii in numpy.arange(mr[0]):
        m_update[ii, :] = m[ii, :] + K@r[ii, :]

    return m_update


def calc_eki_update(m=numpy.array([]), r=numpy.array([]),
                    Cxx=numpy.array([]),
                    Cxy=numpy.array([]),
                    Cyy=numpy.array([]),
                    Cd=numpy.array([]),
                    out=True):
    """
    Calculate Ensemble Kalman update

    see:
        N. K. Chada, Y. Chen, and D. Sanz-Alonso
        Iterative ensemble Kalman methods: A unified perspective
        with some new variants
        Foundations of Data Science, 3, 331-369. –, 2021.

        S. Duffield and S. S. Singh
        Ensemble Kalman inversion for general likelihoods
        Statistics \& Probability Letters, 187
        doi:10.1016/j.spl.2022.109523, 2022

    Created on Jul 9, 2022

    @author: vrath

    """
    if (m.size == 0):
        error("calc_kupdate: No model given! Exit.")
    if (r.size == 0):
        error("calc_kupdate: No residuals given! Exit.")

    ms = numpy.shape(m)
    mr = numpy.shape(r)
    if mr[0] != ms[0]:
        error("calc_kupdate: Ensemble sizes don't match! Exit.")

    if (Cxy.size == 0) or (Cyy.size == 0):
        error("calc_keg_update: One or more covariances missing! Exit.")
    if (Cd.size == 0):
        error("calc_keg_update: Data covariance missing! Exit.")
    else:
        C = numpy.eye(mr[1])/Cd

    Ci = scipy.linalg.inv(Cyy+C)
    K = Cxy@Ci

    m_update = m.copy()
    for ii in numpy.arange(mr[0]):
        m_update[ii, :] = m[ii, :] + K@r[ii, :]

    return m_update


def get_nlyr(model_vec=numpy.array([])):

    if numpy.size(model_vec) == 0 or numpy.size(model_vec) % 7 != 0:
        error("get_nlyr: model_vec is 0 or not a multiple of 7! Exit.")

    nlyr = int(numpy.size(model_vec)/7)

    return nlyr


def calc_jac(
        fwdcall=None,
        alt=None,
        m_vec=numpy.array([]), m_act=numpy.array([]), m_trn=0, m_state=0, m_params=None,
        d_vec=numpy.array([]), d_act=numpy.array([]), d_trn=0,
        delta=[1.0e-4], scalejac=False, out=False):
    """
    Calculate Jacobian

    Parameters
    ----------
    fwdcall:string
        Call to forward modeling tol
    altitude:float
        Flight altitude.  Default is None.
    m_vec:float
        Model vector (7*nlayer). Default is None.
        It is expected in transformed parameters
    m_act:int, optional
        Ativity control vector (7*nlayer). Default is None.
    params:unknown, optional
        Additional parameters for different systems. Default is None.
    m_trn:int, optional
        controls transforms of model to different parametrizations.
        Default is 0.
        0 = only resitivity is transformed to log
        1 = all parameters are transformed to log
    delta:float, optional
        controls divided differences differentiation. Default is 1e-4.
        if adaptive, it ia a numpy array with
        [initial delta, min_delta, delta reduction, min differential change]
        e.g., [1.e-4, 1.e-8, 0.1, 1.e-8]
    adaptive: logical, optional
        switched to adaptive mode. Default is False
    scalejac: logical, optional
        sale jacobian with errors.

    Returns
    -------
    jacobian:float
        Jacobian matrix (nlayer*number of data)
    w_vec:float
        calculated data weights (i.e.inverse errors)
    d_vec:float
        calculated data at central value.

    Author:VR 06/22
    """

    if numpy.isscalar(delta):
        delta = [delta]

    jacobian = numpy.zeros((numpy.size(d_vec), numpy.shape(m_vec)[0]))

    data0, d_state = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                   m_vec=m_vec, m_trn=m_trn, m_state=m_state,
                                   d_trn=d_trn)

    adaptive = False
    if delta[0] < 0:
        adaptive = True
        delta = [1.e-2,  0.1, 1.e-8, 1.e-8]

    if adaptive:
        count = 0
        for ipert in numpy.arange(numpy.size(m_act)):
            if m_act[ipert] != 0:
                count = count + 1

                deliter = delta[0]
                deljacb = 1.e-3
                while (deljacb >= delta[2]) and (deliter >= delta[3]):

                    m_current = m_vec.copy()
                    m_current[ipert] = m_current[ipert] + deliter
                    data1, _ = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                             m_vec=m_current, m_trn=m_trn, m_state=m_state,
                                             d_trn=d_trn)
                    mcurrent = m_vec.copy()
                    mcurrent[ipert] = mcurrent[ipert] - deliter
                    data2, _ = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                             m_vec=m_current, m_trn=m_trn, m_state=m_state,
                                             d_trn=d_trn)
                    deljacb = numpy.linalg.norm(
                        ((data1-data0) - (data0-data2))/deliter)
                    # print(deljacb, ipert, deliter)
                    jacobian[:, ipert] = 0.5*(data1 - data2) / (deliter)
                    deliter = deliter*delta[1]

    else:

        for ipert in numpy.arange(numpy.size(m_act)):
            if m_act[ipert] != 0:

                m_current = m_vec.copy()
                m_current[ipert] = m_current[ipert] + delta[0]
                data1, _ = calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                         m_vec=m_current, m_trn=m_trn, m_state=m_state,
                                         d_trn=d_trn)
                jacobian[:, ipert] = (data1 - data0) / delta[0]

    if numpy.size(d_act) > 0 and numpy.size(m_act) > 0:
        jacobian = extract_jac(J=jacobian, m_act=m_act, d_act=d_act)

    return jacobian

# def calc_sensitivity(Jac=numpy.array([]),
#                      Type = "euclidean", UseSigma = False, small_val = 1.e-30, OutInfo = False):
#     """
#     Calculate sensitivities.
#     Expects that Jacobian is already scaled, i.e Jac = C^(-1/2)*J.

#     Several options exist for calculating sensiotivities, all of them
#     used in the literature.
#     Type:
#         "raw"     sensitivities summed along the data axis
#         "abs"     absolute sensitivities summed along the data axis
#                     (often called coverage)
#         "euc"     squared sensitivities summed along the data axis.
#         "cum"     cummulated sensitivities as proposed by
#                   Christiansen & Auken, 2012. Not usable for negative data.

#     Usesigma:
#         if true, sensitivities with respect to sigma  are calculated.

#     Christiansen, A. V. & Auken, E.
#     A global measure for depth of investigation
#     Geophysics, 2012, 77, WB171-WB177

#     from UBC:
#     def depth_of_investigation_christiansen_2012(self, std, thres_hold=0.8):
#         pred = self.survey._pred.copy()
#         delta_d = std * numpy.log(abs(self.survey.dobs))
#         J = self.getJ(self.model)
#         J_sum = abs(Utils.sdiag(1/delta_d/pred) * J).sum(axis=0)
#         S = numpy.cumsum(J_sum[::-1])[::-1]
#         active = S-thres_hold > 0.
#         doi = abs(self.survey.depth[active]).max()
#         return doi, active

#     T. Guenther
#         Inversion Methods and Resolution Analysis for the 2D/3D Reconstruction
#         of Resistivity Structures from DC Measurements
#         Fakultaet für Geowissenschaften, Geotechnik und Bergbau,
#         Technische Universitaet Bergakademie Freiberg, 2004.

#     author:VR 9/23

#     """

#     if numpy.size(Jac)==0:
#         error("calc_sensitivity: Jacobian size is 0! Exit.")

#     if UseSigma:
#         Jac = -Jac



#     if "raw" in  Type.lower():
#         S = Jac.sum(axis=0)
#         if OutInfo:
#             print("raw:", S)
#         # else:
#         #     print("raw sensitivities")
#         # smax = Jac.max(axis = 0)
#         # smin = Jac.max(axis = 0)

#     elif "cov" in Type.lower():
#         S = Jac.abs().sum(axis=0)
#         if OutInfo:
#             print("cov:", S)
#         # else:
#         #     print("coverage")

#     elif "euc" in Type.lower():
#         S = Jac.power(2).sum(axis=0)
#         if OutInfo:
#             print("euc:", S)
#         # else:
#         #     print("euclidean (default)")

#     elif "cum" in Type.lower():
#         S = Jac.abs().sum(axis=0)
#         # print(numpy.shape(S))
#         # S = numpy.sum(Jac,axis=0)

#         S = numpy.append(0.+1.e-10, numpy.cumsum(S[-1:0:-1]))
#         S = numpy.flipud(S)
#         if OutInfo:
#            print("cumulative:", S)
#         # else:
#         #    print("cumulative sensitivity")

#     else:
#         print("calc_sensitivity: Type "
#               +Type.lower()+" not implemented! Default assumed.")
#         S = Jac.power(2).sum(axis=0)

#         if OutInfo:
#             print("euc (default):", S)
#         # else:
#         #     print("euclidean (default)")

#         # S = S.reshape[-1,1]

#     S[numpy.where(numpy.abs(S)<small_val)]=small_val
#     print("calc: ", numpy.any(S==0))
#     # S=S.A1
#     S = numpy.asarray(S).ravel()
#     return S

def calc_sensitivity(Jac=numpy.array([]),
                     sens_type="euclidean", use_sigma=False, small_val = 1.e-30, OutInfo=False):
    """
    Calculate sensitivities.
    Expects that Jacobian is already sclaed, i.e Jac = C^(-1/2)*J.

    Several options exist for calculating sensiotivities, all of them
    used in the literature.
    type:
        "raw"     sensitivities summed along the data axis
        "abs"     absolute sensitivities summed along the data axis
                    (often called coverage)
        "euc"     squared sensitivities summed along the data axis.
        "cum"     cummulated sensitivities as proposed by
                  Christiansen & Auken, 2012. Not usable for negative data.

    Usesigma:
        if true, sensitivities with respect to sigma  are calculated.

    Christiansen, A. V. & Auken, E.
    A global measure for depth of investigation
    Geophysics, 2012, 77, WB171-WB177

    from UBC:
    def depth_of_investigation_christiansen_2012(self, std, thres_hold=0.8):
        pred = self.survey._pred.copy()
        delta_d = std * numpy.log(abs(self.survey.dobs))
        J = self.getJ(self.model)
        J_sum = abs(Utils.sdiag(1/delta_d/pred) * J).sum(axis=0)
        S = numpy.cumsum(J_sum[::-1])[::-1]
        active = S-thres_hold > 0.
        doi = abs(self.survey.depth[active]).max()
        return doi, active

    T. Guenther
    Inversion Methods and Resolution Analysis for the 2D/3D Reconstruction
    of Resistivity Structures from DC Measurements
    Fakultaet für Geowissenschaften, Geotechnik und Bergbau,
    Technische Universitaet Bergakademie Freiberg, 2004.

    author:VR 4/23

    """

    if numpy.size(Jac) == 0:
        error("calc_sensitivity: Jacobian size is 0! Exit.")

    if use_sigma:
        Jac = -Jac

    if "raw" in sens_type.lower():
        S = numpy.sum(Jac, axis=0)
        if OutInfo:
            print("raw:", S)


    elif "cov" in sens_type.lower():
        S = numpy.sum(numpy.abs(Jac), axis=0)
        if OutInfo:
            print("cov:", S)


    elif "euc" in sens_type.lower():
        S = numpy.sum(numpy.power(Jac, 2), axis=0)
        if OutInfo:
            print("euc:", S)

    elif "cum" in sens_type.lower():
        S = numpy.sum(numpy.abs(Jac), axis=0)
        S = numpy.append(0.+1.e-10, numpy.cumsum(S[-1:0:-1]))
        S = numpy.flipud(S)
        if OutInfo:
            print("cumulative:", S)

    else:
        print("calc_sensitivity: type "
              + type.lower()+" not implemented! Default assumed.")
        S = numpy.sum(numpy.power(Jac, 2), axis=0)

        if OutInfo:
            print("euc (default):", S)

    S[numpy.where(numpy.abs(S)<small_val)]=small_val

    S = numpy.asarray(S).ravel()
    return S

def transform_sensitivity(S=numpy.array([]), vol=numpy.array([]),
                          transform=["size","max", "sqrt"],
                          asinhpar=[0.], max_val=None, small_val= 1.e-30, OutInfo=False):
    """
    Transform sensitivities.

    Several options exist for transforming sensitivities, all of them
    used in the literature.

    Normalize options:
        "siz"       Normalize by the values optional array V ("volume"),
                    i.e in our case layer thickness.
        "max"       Normalize by maximum value.
        "sur"       Normalize by surface value.
        "sqr"       Take the square root. Only usefull for euc sensitivities.
        "log"       Take the logaritm. This should always be the
                    last value in Transform list

        "asinh"     asinh transform. WARNING: excludes log option,
                    and should be used only for raw sensitivities
                    (C. Scholl, Die Periodizitaet von Sendesignalen
                    bei Long-Offset Transient Electromagnetics,
                    Diploma Thesis, Universität zu Koeln, 2001).

    author:VR 4/23

    """

    if numpy.size(S)==0:
        error("transform_sensitivity: Sensitivity size is 0! Exit.")

    # ns = numpy.shape(S)


    for item in transform:


        if "sqr" in item.lower():
            S = numpy.sqrt(S)
            # print("S0s", numpy.shape(S))

        if "log" in item.lower():
            S = numpy.log10(S)

        if "asinh" in item.lower():
            maxval = numpy.amax(S)
            minval = numpy.amin(S)
            if maxval>0 and minval>0:
                print("transform_sensitivity: No negatives, switched to log transform!")
                S = numpy.log10(S)
            else:
                if len(asinhpar)==1:
                    scale = asinhpar[0]
                else:
                    scale = get_scale(S, method=asinhpar[0])

                    S = numpy.arcsinh(S/scale)

        if "siz" in item.lower():
             print("trans_sensitivity: Transformed by volumes/layer thickness.")
             if numpy.size(vol)==0:
                 error("Transform_sensitivity: no volumes given! Exit.")

             else:
                 S = S/vol

        if "max" in item.lower():
             print("trans_sensitivity: Transformed by maximum value.")
             if max_val is None:
                 maxval = numpy.amax(numpy.abs(S))
             else:
                 maxval = max_val
             print("maximum value: ", maxval)
             S = S/maxval
             # print("S0m", numpy.shape(S))


        S[numpy.where(numpy.abs(S)<small_val)]=small_val


    return S, maxval


def get_scale(d=numpy.array([]), f=0.1, method = "other", OutInfo = False):
    """
    Get optimal Scale for arcsin transformation.

    Parameters
    ----------
    d : float, required.
        Data vector.
    F : float, optional
        Weight for arcsinh transformation, default from Scholl & Edwards (2007)

    Returns
    -------
    S : float
        Scale value for arcsinh

    C. Scholl
        Die Periodizitaet von Sendesignalen bei Long-Offset Transient Electromagnetics
        Diploma Thesis, Institut für Geophysik und Meteorologie der Universität zu Koeln, 2001.


    """

    if numpy.size(d)==0:
        error("get_S: No data given! Exit.")

    if "s2007" in method.lower():
        scale = f * numpy.nanmax(numpy.abs(d))

    else:
        dmax = numpy.nanmax(numpy.abs(d))
        dmin = numpy.nanmin(numpy.abs(d))
        denom =f *(numpy.log(dmax)-numpy.log(dmin))
        scale = numpy.abs(dmax/denom)

    if OutInfo:
        print("Scale value S is "+str(scale)+", method "+method)

    return scale


def calc_resolution_matrices(J=numpy.array([]), G=numpy.array([]),
                             OutInfo=True):
    """
    Calculate coverages.

    T. Guenther
    Inversion Methods and Resolution Analysis for the 2D/3D Reconstruction
    of Resistivity Structures from DC Measurements
    Fakultaet für Geowissenschaften, Geotechnik und Bergbau,
    Technische Universitaet Bergakademie Freiberg, 2004.

    S. Friedel
    Resolution, stability and efficiency of resistivity tomography estimated
    from a generalized inverse approach
    Geophysical Journal International, vol. 153, Art. no. 2, 2003.

    author:VR  last change 3/23
    """
    if numpy.size(J) == 0 or numpy.size(G) == 0:
        error("calc_resolution_matrices: J or G not given. Exit")

    Rd = J@G
    Rm = G@J

    return Rm, Rd


def calc_model_resolution(J=numpy.array([]), G=numpy.array([]),
                          Spread=[],
                          alpha=1.E-6,
                          OutInfo=True):
    """
    Calculate coverages.

    T. Guenther
    Inversion Methods and Resolution Analysis for the 2D/3D Reconstruction
    of Resistivity Structures from DC Measurements
    Fakultaet für Geowissenschaften, Geotechnik und Bergbau,
    Technische Universitaet Bergakademie Freiberg, 2004.

    S. Friedel
    Resolution, stability and efficiency of resistivity tomography estimated
    from a generalized inverse approach
    Geophysical Journal International, vol. 153, Art. no. 2, 2003.


    Z. Adavi, R. Weber, and W. Rohm
    Pre-analysis of GNSS tomography solution using the
    concept of spread of model resolution matrix
    Journal of Geodesy, vol. 96, 2023, doi: 10.1007/s00190-022-01620-1.

    C. R. Miller and P.S. Routh
    Resolution analysis of geophysical images: Comparison between point
    spread function and region of data influence measures
    Geophysical Prospecting, 2007, 55, 835–852
    doi:10.1111/j.1365-2478.2007.00640.x



    author:VR  last change 4/23
    """

    if numpy.size(J) == 0 or numpy.size(G) == 0:
        error("calc_model_resolution: J or G not given. Exit")

    smallval = 100.*numpy.finfo(float).eps

    Rm = G@J

    N = numpy.shape(Rm)[0]

    if len(Spread) != 0:
        stype = Spread[0]
        if OutInfo:
            print("calc_model_resolution: Spread definition is "+stype.lower())

        if "men" in stype.lower() or "fro" in stype.lower():
            R = Rm - numpy.eye(N)
            S = numpy.linalg.norm(R)

        W = numpy.zeros((N, N))
        if "too" in stype.lower() or "back" in stype.lower():
            lam = 2
            if len(Spread) > 1:
                lam = Spread[1]
            for i1 in numpy.arange(N):
                for i2 in numpy.arange(N):
                    W[i1, i2] = numpy.abs(i1-i2)

            R = numpy.power((W * Rm), lam)
            S = numpy.sum(R, axis=0)

        if "mil" in stype.lower():
            R = Rm - numpy.eye(N)

            lam = 2
            if len(Spread) > 1:
                lam = Spread[1]

            for i1 in numpy.arange(N):
                for i2 in numpy.arange(N):
                    W[i1, i2] = 1.0 + numpy.power(numpy.abs(i1-i2), lam)

            S = numpy.sum(W * R**2, axis=0)
            D = alpha + numpy.sum(Rm, axis=0)
            S = S/D

        if "mic" in stype.lower():
            d = numpy.zeros_like(Rm)
            S = numpy.zeros(N)

            for i1 in numpy.arange(N):
                for i2 in numpy.arange(N):
                    d[i2] = numpy.abs(i1-i2)

                R = Rm[i1, :]
                D = d[:]
                f = 1./scipy.linalg.norm(R+smallval)
                S[i1] = numpy.log(f*numpy.sum(numpy.power((f*R), 2)*D))

        return Rm, S

    else:
        if OutInfo:
            print("calc_model_resolution: Spread not calculated. Exit")
        return Rm


def transform_data(d_vec=numpy.array([]), e_vec=numpy.array([]),
                   d_trn=0, d_state=0, get_weights=False, S=1000., OutInfo=False):
    """

    Transform data.

    Parameters
    ----------
    d_vec : float, optional
        Data vector. The default is None.
    d_trn : integer, optional
        Data transform. 1: log, 2: arcsinh.
        The default is 2.
    get_weights : boolean, optional
        Detemines whether error weights for inversion are generated.
        The default is False.
    e_vec : float, optional
        Errors, required if get_weights is True. The default is None.
    S : float, optional
        Scale value for arcsing transform. Required for arcsinh transform = 1000.
        The default is None.

    Returns
    -------
    d_trans : float
        Transformed data.
    w_trans : float, optional
        weights for inversion, if get_weights is True.


    C. Scholl
    Die Periodizitaet von Sendesignalen bei Long-Offset Transient Electromagnetics
    Diploma Thesis, Institut für Geophysik und Meteorologie der Universität zu Koeln, 2001.

    A. Haroon
    Development of Novel Time-domain Electromagnetic Methods for Offshore
    Groundwater Studies: A Data Application from Bat Yam, Israel
    Mathematisch-Naturwissenschaftliche Fakultaet der Universität zu Koeln, 2016.

    author: VR 2/22
    """

    if numpy.size(d_vec) == 0:
        error("transform_data: No data given! Exit.")

    if numpy.size(e_vec) == 0:
        e_vec = numpy.ones_like(d_vec)
        if OutInfo:
            print("transform_data: No errors given! Set = 1.")

    if d_state != 0:
        print("transform_data: no transform possible, d_trans !=0!")
        return d_vec, e_vec, d_state

    e_trans = numpy.zeros_like(d_vec)
    w_trans = numpy.zeros_like(d_vec)

    if numpy.abs(d_trn) == 0:
        d_trans = d_vec
        if get_weights:

            w_trans = w_trans + \
                numpy.abs(numpy.ones_like(e_vec) /
                          e_vec).reshape(numpy.shape(w_trans))
        else:
            e_trans = e_vec

    elif numpy.abs(d_trn) == 1:
        if OutInfo:
            print("Data transformed with ln function")

        d_vec[numpy.where(numpy.abs(d_vec) < 1.e-16)] = numpy.nan
        d_trans = numpy.log(d_vec)

        if get_weights:
            w_trans = w_trans + (d_vec/e_vec).reshape(numpy.shape(w_trans))
        else:
            e_trans = e_trans + (e_vec/d_vec).reshape(numpy.shape(e_trans))

    elif numpy.abs(d_trn) == 2:
        if OutInfo:
            print("Data transformed with arcsinh function")

        d_vec[numpy.where(numpy.abs(d_vec) < 1.e-16)] = numpy.nan

        if S == 0.:
            S = get_S(d_vec)

            if OutInfo:
                print("S value = "+str(S))
        d_trans = numpy.arcsinh(d_vec/S)

        if get_weights:
            w_trans = w_trans + (numpy.sqrt(d_vec**2 + S**2) /
                                 e_vec).reshape(numpy.shape(w_trans))
        else:
            e_trans = e_trans + \
                (e_vec/numpy.sqrt(d_vec**2 + S**2)).reshape(numpy.shape(e_trans))

    else:
        error(
            "Data transformation " + str(d_trn) + " not implemented! Exit."
        )

    d_state_new = d_trn

    if get_weights:
        return d_trans, w_trans, d_state_new

    else:

        return d_trans, e_trans, d_state_new


def get_S(d=numpy.array([]), F=0.1, method="other", OutInfo=False):
    """
    Get optimal Scale for arcsin transformation.

    Parameters
    ----------
    d : float, required.
        Data vector.
    F : float, optional
        Weight for arcsinh transformation, default from Scholl & Edwards (2007)

    Returns
    -------
    S : float
        Scale value for arcsinh

    C. Scholl
    Die Periodizitaet von Sendesignalen bei Long-Offset Transient Electromagnetics
    Diploma Thesis, Institut für Geophysik und Meteorologie der Universität zu Koeln, 2001.


    """

    if numpy.size(d) == 0:
        error("get_S: No data given! Exit.")

    if "s2007" in method.lower():
        S = F * numpy.nanmax(numpy.abs(d))

    else:
        dmax = numpy.nanmax(numpy.abs(d))
        dmin = numpy.nanmin(numpy.abs(d))
        denom = F * (numpy.log(dmax)-numpy.log(dmin))
        S = numpy.abs(dmax/denom)

    if OutInfo:
        print("Scale value S is "+str(S)+", method "+method)

    return S


def transform_parameter(m_vec=numpy.array([]),
                        mode="fwd", m_trn=0, m_state=0,
                        bounds=None, dp=False, deltap=numpy.array([1.e-5])):
    """
    m_trn model to different parametrizations.

    C. Scholl
    Die Periodizitaet von Sendesignalen bei Long-Offset Transient Electromagnetics
    Diploma Thesis, Institut für Geophysik und Meteorologie der Universität zu Koeln, 2001.

    A. Haroon
    Development of Novel Time-domain Electromagnetic Methods for Offshore
    Groundwater Studies: A Data Application from Bat Yam, Israel
    Mathematisch-Naturwissenschaftliche Fakultaet der Universität zu Koeln, 2016.

    author: VR 2/22
    """
    npar = 7

    if m_vec.size == 0:
        error("transform_para: transform possible, no model given! Exit.")

    mshape = numpy.shape(m_vec)
    # print(mshape)
    nlyr = mshape[0]//npar

    m_trans = m_vec.copy()
    # print(m_trans.shape)
    m_trans = numpy.reshape(m_trans, (npar, nlyr))

    if "b" in mode:
        dp = False
        if m_state == 0:
            print("transform_parameter: no transform possible, m status is " +
                  str(m_state)+"!")
            return m_vec, m_state

        elif numpy.abs(m_trn) == 1:
            m_trans[0, :] = numpy.exp(m_trans[0, :])

        elif numpy.abs(m_trn) == 2:
            m_trans[0, :] = numpy.exp(m_trans[0, :])
            m_trans[6, :-1] = numpy.exp(m_trans[6, :-2])

        elif numpy.abs(m_trn) == 3:
            # e^mu/(e^mu + 1), INVERSE normalized chargeability Ghorbani 2007
            m_trans[0, :] = numpy.exp(m_trans[0, :])
            m_trans[3, :] = numpy.exp(m_trans[3, :]) / \
                (numpy.exp(m_trans[3, :]) + 1.0e0)
            m_trans[4, :] = numpy.exp(m_trans[4, :])
            m_trans[6, :-1] = numpy.exp(m_trans[6, :-2])

        else:
            error(
                "Parameter transformation " +
                str(m_trn) + " not implemented! Exit."
            )

        m_state_new = 0

    elif "f" in mode:

        if m_state != 0:
            print("transform_parameter: no transform possible, m status is " +
                  str(m_state)+"!")
            return m_vec, m_state

        if dp:
            if deltap.size == 1:
                delt = deltap * numpy.ones_like(m_vec)
            else:
                delt = numpy.ones_like(m_vec)*deltap

        if numpy.abs(m_trn) == 0:
            pass

        elif numpy.abs(m_trn) == 1:
            m_trans[0, :] = numpy.log(m_trans[0, :])
        elif numpy.abs(m_trn) == 2:
            m_trans[0, :] = numpy.log(m_trans[0, :])
            m_trans[6, :-1] = numpy.log(m_trans[6, :-2])

        elif numpy.abs(m_trn) == 3:
            # e^mu/(1e^mu + 1), INVERSE normalized chargeability Ghorbani 2001
            m_trans[0, :] = numpy.log(m_trans[0, :])
            m_trans[3, :] = numpy.log(m_trans[3, :] / (1e0 - m_trans[3, :]))
            m_trans[4, :] = numpy.log(m_trans[4, :])
            m_trans[6, :-1] = numpy.log(m_trans[6, :-2])
        else:
            error(
                "Parameter transformation " +
                str(m_trn) + " not implemented! Exit."
            )

        m_state_new = m_trn

    if not bounds == None:
        m_trans = impose_bounds(m=m_trans, bounds=bounds, mode=mode)

    m_trans = numpy.reshape(m_trans, mshape)

    if dp:
        return m_trans, delt, m_state_new
    else:
        return m_trans, m_state_new


def impose_bounds(m=None, bounds=None, mode="fwd",
                  method = "kim",f=1.):
    """
    Generate barrier operator used for enforcing limits for parameters.

    Based on:
    H J Kim & Y H Kim, A unified transformation function for lower and upper
    bounding constraints on model parameters in electrical and electromagnetic
    inversion, J. Geophys. Eng. 8 (2011) 21–26

    f = 1.  results in  the scheme described in Kim (1999).
    f = 2.  results in  the scheme described by Commer (2008)
            or Abubakar (2008).

    The logistic sigmoid transform is suggested in SimPEG's documentation at:
        https://docs.simpeg.xyz/content/api/generated/SimPEG.maps.LogisticSigmoidMap.dot.html

    VR Apr 2024

    """
    lower = bounds[0]
    upper = bounds[1]

    if "kim" in method.lower():

        if mode[0:1] == "f":
            m_trans = 1./f*numpy.log((m-lower)/(upper-m))
        elif mode[0:1] == "b":
            ex = numpy.exp(f*m)
            m_trans = (lower + upper*ex)/(1. + ex)#


    elif "sig" in method.lower():
        fac = bounds[:,0] +(bounds[:,1]-bounds[:,0])
        if mode[0:1] == "f":
            m_trans = fac*scipy.special.expit(m)
        elif mode[0:1] == "b":
            m_trans = scipy.special.logit(m)/fac

    return m_trans


def init_obsdata(nD, value="active"):
    """
    Initializes observational data structure inversion as all active,
    all incative, or NaN.
    last changes: VR 06/22
    """
    d_state = 0

    onobs = numpy.ones(nD)
    d_obs = numpy.zeros(nD)

    if value.lower() == "active":
        d_act = onobs.astype("int8")
    elif value.lower() == "inactive":
        d_act = 0 * onobs.astype("int8")
    elif value.lower() == "nan":
        d_act = numpy.nan * onobs.astype("int8")
    else:
        error("d_act initial  " + value + " not implemeted !")

    return d_obs, d_act, d_state


def init_1dmod(nlyr, npar=7):
    """
    Initializes data structures for 1D inversion
    parameters as all inactive (m_act)
    priors (here:gaussian) to reasonable values for res/dz-only inversion
    bounds to practical +-infinite
    last changes: VR 04/23
    """
    epsi = 1.e-12
    m_state = 0

    # m_act determines the active paramteters, currently npar per layer
    m_act = numpy.zeros((npar*nlyr)).astype("int8")

    # m_prior is the default prior model
    m_prior = numpy.ones((npar*nlyr))
    m_prior[3*nlyr:4*nlyr] = epsi
    m_prior[4*nlyr:5*nlyr] = epsi

    # mvar is the default prior variance
    m_var = numpy.ones((npar*nlyr))

    # m_bounds are bounds on the parameters (currently no used)
    m_bounds = numpy.ones((npar*nlyr, 2))
    m_max = 1e30
    m_min = -1e30
    m_bounds[:, 0] = m_min
    m_bounds[:, 1] = m_max
    #   =  6. #prior_avg[m_act!=0] + 3*prior_std[m_act!=0]
    #   = -1. #prior_avg[m_act!=0] - 3*prior_std[m_act!=0]

    return m_act, m_prior, m_var, m_bounds, m_state


def generate_random_vec(m0, m_act, covar, seed=None):
    """
    Generate random models with prescribed covariance.

    VR Jan 2021
    """
    sizepar = numpy.shape(m0)

    if seed is None:
        rng = numpy.random.default_rng()
    else:
        rng = numpy.random.default_rng(seed)

    m_gauss = rng.standard_normal(sizepar[0])

    l = scipy.linalg.cholesky(covar)

    m = m0 + m_act * m_gauss * l

    return m


def calc_aniso_euclidean_norm(x1, x2, scale=None):
    """
    anisotropic scaling for radial basic function interpolation
    """

    if scale == None:
        fn = numpy.sqrt(((x1 - x2)**2).sum(axis=0))
    else:
        scale = scale/numpy.sum(scale)
        f = scale*(x1 - x2)**2
        fn = numpy.sqrt((f).sum(axis=0))

    return fn


def insert_dat(D=numpy.array([]), d=numpy.array([]), d_act=numpy.array([])):
    """
    Unfold data according to both, m_act and d_act.

    VR Jan 2021

    """
    if numpy.size(D) == 0:
        error("insert_dat: D not defined! Exit.")
    if numpy.size(d) == 0:
        error("insert_dat: d not defined! Exit.")
    if numpy.size(d_act) == 0:
        error("insert_dat: d_act not defined! Exit.")

    A = D.copy()

    if d.ndim == 1:
        A[d_act != 0] = d
    else:
        for ii in numpy.arange(A.shape[0]):
            A[ii, d_act.flat != 0] = d[ii, :]

    return A


def extract_dat(D=numpy.array([]), d_act=numpy.array([])):
    """
    extract data according to d_act .

    VR Jan 2021

    """
    if numpy.size(D) == 0:
        error("extract_dat: D not defined! Exit.")

    if numpy.size(d_act) == 0:
        error("extract_dat: d_act not defined! Exit.")

    if numpy.size(d_act) != numpy.size(D):
        error("extract_dat: sizes do not match! Exit.")

    tmp = D.copy()

    if D.ndim == 1:
        d = tmp[d_act.flat != 0]
    else:
        for ii in numpy.arange(D.shape[0]):
            if ii == 0:
                d = tmp[0, d_act.flat != 0]
            else:
                d = numpy.vstack((d, tmp[ii, d_act.flat != 0]))

    return d


def insert_jac(M=numpy.array([]),
               m=numpy.array([]),
               m_act=numpy.array([]),
               d_act=numpy.array([])):
    """
    Unfold Jacobian according to m_act and d_act.

    VR Jun 2022

    """
    if numpy.size(M) == 0:
        error("insert_dat: M not defined! Exit.")

    if m.size == 0:
        error("insert_dat: m not defined! Exit.")

    if numpy.size(m_act) == 0:
        error("insert_dat: m_act not defined! Exit.")

    if numpy.size(d_act) == 0:
        error("insert_jac: d_act not defined! Exit.")

    A = M.copy()
    A[d_act.flat != 0, m_act.flat != 0] = m

    return A


def extract_jac(J=numpy.array([]),
                d_act=numpy.array([]),
                m_act=numpy.array([])):
    """
    extract Jacobian according to m_act and d_act.

    VR Jan 2021

    """
    if J.size == 0:
        error("extract_jac: J not defined! Exit.")
    if numpy.size(m_act) == 0:
        error("extract_jac: m_act not defined! Exit.")
    if numpy.size(d_act) == 0:
        error("extract_jac: d_act not defined! Exit.")

    tmp = J.copy()
    tmp = tmp[d_act.flat != 0, :]
    A = tmp[:, m_act.flat != 0]

    return A


def insert_mod(M=numpy.array([]),
               m=numpy.array([]),
               m_act=numpy.array([])):
    """
    insert according to m_act.
    Allows for mdel ensemble .

    VR Jul 2022

    """
    if numpy.size(M) == 0:
        error("insert_dat: M not defined! Exit.")

    if numpy.size(m) == 0:
        error("insert_dat: m not defined! Exit.")

    if numpy.size(m_act) == 0:
        error("insert_dat: m_act not defined! Exit.")

    A = M.copy()
    # print(numpy.shape(A), numpy.shape(m), numpy.shape(m_act),)
    if M.ndim == 1:
        A[m_act.flat != 0] = m.flat
    else:
        for ii in numpy.arange(M.shape[0]):
            A[ii, m_act.flat != 0] = m[ii, :]

    return A


def extract_mod(M=numpy.array([]), m_act=numpy.array([])):
    """
    extract model according to m_act .
    Allows for mdel ensemble .

    VR Jul 2021

    """
    if numpy.size(M) == 0:
        error("extract_mod: M not defined! Exit.")
    if numpy.size(m_act) == 0:
        error("extract_mod: m_act not defined! Exit.")

    tmp = M.copy()
    if M.ndim == 1:
        m = tmp[m_act.flat != 0]
    else:
        for ii in numpy.arange(M.shape[0]):
            m[ii, :] = tmp[ii, m_act.flat != 0]

    return m


def full_cov(C=[numpy.array([])]):
    """
    setup full block-diagonal covar matrix
    """
    # CmiS = scipy.sparse.block_diag([CmiS for ii in range(7)])

    if len(C) == 0:
        error("ful_cov: C not defined! Exit.")

    if len(C) == 1:
        tmp = [C[0] for ii in numpy.arange(7)]
        block_cov = scipy.linalg.block_diag(*tmp)
    else:
        block_cov = scipy.linalg.block_diag(*C)

    return block_cov


def extract_cov(C=numpy.array([]), m_act=numpy.array([])):
    """
    extract covariance according to m_act.

    VR Jan 2021

    """
    if C.size == 0:
        error("extract_cov: C not defined! Exit.")
    if numpy.size(m_act) == 0:
        error("extract_cov: m_act not defined! Exit.")
    print(numpy.shape(C), numpy.shape(m_act))
    tmp = C.copy()
    if scipy.sparse.issparse(tmp):
        tmp = tmp.todense()
    tmp = tmp[m_act.flat != 0, :]
    A = tmp[:, m_act.flat != 0]
    A = scipy.sparse.csr_matrix(A)

    return A


def extract_wgt(W=numpy.array([]), m_act=numpy.array([])):
    """
    extract weights (covar^1/2) according to m_act.

    VR Jan 2021

    """
    if W.size == 0:
        error("extract_wgt: W not defined! Exit.")
    if numpy.size(m_act) == 0:
        error("extract_wgt: m_act not defined! Exit.")

    tmp = W.copy()
    # tmp = tmp.todense()
    tmp = tmp[m_act != 0, :]
    A = tmp[:, m_act != 0]
    A = scipy.sparse.csr_matrix(A)

    return A


def load_prior(prior_file=None,
               m_ref=numpy.array([]),
               m_act=numpy.array([]),
               m_apr=numpy.array([]),
               out=False):
    """
    Loads pre-calculated prior file.
    File should be outout of earlier inversion)

    Parameters
    ----------
    prior_file : string
        name of the prior file.
    m_ref, m_act : numpy.arrays
        reference model and activation flag for current run


    OutInfo : logical, optional
        determines outputs.

    Returns
    -------
    mod_apr: numpy.array
        prior model

    VR April 2023

    """
    if prior_file == None:
        error("No prior file given! Exit.")

    if not os.path.exists(prior_file):
        error("Prior file "+prior_file+" does not exist! Exit.")

    tmp = numpy.load(prior_file, allow_pickle=True)

    prior_ref = tmp["mod_ref"]
    prior_act = tmp["mod_act"]
    prior_sit = tmp["site_modl"]

    m_sit = m_ref[m_act!=0]
    print("prior current:")
    print(numpy.shape(m_ref))
    print(numpy.shape(m_act))
    print(numpy.shape(m_apr))
    print(" prior read:")
    print(numpy.shape(prior_ref))
    print(numpy.shape(prior_act))
    print(numpy.shape(prior_sit))


    nsit = numpy.shape(prior_sit)[0]
    npar = numpy.shape(m_apr)[0]
    prior_mod = numpy.zeros((nsit,npar))
    print(" prior out:")
    print(numpy.shape(prior_mod))


    if (numpy.shape(prior_ref)==numpy.shape(m_ref)) and (numpy.shape(prior_act)==numpy.shape(m_act)):

        # prior run with same base model

        nsit = numpy.arange(numpy.shape(prior_sit)[0])
        for isit in nsit:
            prior_mod[isit,:] = insert_mod(m_ref, prior_sit[isit,:], m_act)

    elif get_nlyr(prior_ref)==1:
        # halfspace-prior

        m_sit = m_ref[m_act!=0]
        m_sit = prior_sit*numpy.ones_like(m_sit)
        nsit = numpy.arange(numpy.shape(prior_sit)[0])
        for isit in nsit:
            prior_mod[isit,:] = insert_mod(m_ref, m_sit[isit,:], m_act)

    return prior_mod


def diffops(dz=None, der=False,
            otype="L0", variant=0, mtype="sparse", mform="csr", flip=False, out=True):
    """
    Generate differential operators L0-L2 potentially based on dz.

    See:
    P. C. Hansen
    Discrete inverse problems: Insight and algorithms, SIAM, 2010
    pp 173ff

    A. Doicu, T. Trautmann, and F. Schreier
    Numerical regularization for atmospheric inverse problems,
    Springer, 2010.

    J. Xu, R. Lanlan, F. Schreier, D. Efremenko, A. Doicu, and T. Trautmann
    Insight into Construction of Tikhonov-Type Regularization for
    Atmospheric Retrievals,
    Atmosphere, 11, doi:10.3390/atmos11101052, 2020.

    M. Donatelli and L. Reichel
    Square smoothing regularization matrices with accurate boundary conditions,
    J. Comp. Appl. Math., 272, doi:10.1016/j.cam.2013.08.015, 2014.

    C. Scholl, S. Helwig, B. Tezkan, M. Goldman, U. Kafri (2009)
    1-D multimodel joint inversion of TEM-data over multidimensional structures
    Geophys J Int  , Vol. 176(1), 81-94



    Last change :VR July 2023

    """
    layers = numpy.shape(dz)
    nlyr = layers[0]

    if otype == "L0":
        d = numpy.ones((1, nlyr))
        L = scipy.sparse.spdiags(d, [0], nlyr, nlyr, format=mform)

    elif otype == "L1":

        if nlyr==1:
            d = numpy.ones((1, nlyr))
            L = scipy.sparse.spdiags(d, [0], nlyr, nlyr, format=mform)

        else:

            if der:
                z = numpy.append(0.0, numpy.cumsum(dz))
                zc = 0.5 * (z[0:nlyr] + z[1:nlyr+1])
                # py.shape(zc))
                # zc = numpy.append(zc, zc[nlyr - 2] + dz[nlyr - 2])
                h = 1.0 / numpy.diff(zc)

            if variant == 0:
                d = numpy.zeros((2, nlyr))
                d[0, :] = -1.
                d[1, :] = 1.
                if der:
                    d[:, 1:] = d[:, 1:]*h[:]

                L = scipy.sparse.spdiags(
                    d, [0, -1], nlyr, nlyr-1, format=mform).transpose()
                # if der:
                #     L = scipy.sparse.spdiags(h, [0], nlyr-1, nlyr-1, format=mform)*L
                if out: print("L1 matrix variant 0 is "+str(numpy.shape(L)))

            elif variant == 1:
                alpha = 1.
                d = numpy.zeros((2, nlyr))
                d[0, 1:] = 1.
                d[1, 1:] = -1.
                d[0, 0] = alpha
                if der:
                    d[:, 1:] = d[:, 1:]*h[:]

                L = scipy.sparse.spdiags(
                    d, [0, 1], nlyr, nlyr, format=mform).transpose()

                if out: print("L1 matrix variant 1 is "+str(numpy.shape(L)))

            elif variant == 2:
                alpha = 1.e-6
                d = numpy.zeros((2, nlyr))
                d[0, 1:] = 1.
                d[1, 1:] = -1.
                d[0, 0] = alpha
                if der:
                    d[:, 1:] = d[:, 1:]*h[:]

                L = scipy.sparse.spdiags(
                    d, [0, 1], nlyr, nlyr, format=mform).transpose()

                if out: print("L1 matrix variant 2 is "+str(numpy.shape(L)))

            elif variant == 3:
                alpha = 1.e-6
                d = numpy.zeros((2, nlyr))
                d[0, 1:] = 1.
                d[1, 1:] = -1.
                d[0, 0] = alpha
                d[-1, -1] = 1.
                if der:
                    d[:, 1:] = d[:, 1:]*h[:]

                L = scipy.sparse.spdiags(
                    d, [0, 1], nlyr, nlyr, format=mform).transpose()

                if out: print("L1 matrix variant 3 is "+str(numpy.shape(L)))

            elif variant == -1:
                alpha = 1.
                d = numpy.zeros((2, nlyr))
                d[0, :-1] = 1.
                d[1, :-1] = -1.
                d[-1, -1] = alpha
                if der:
                    d[:, 1:] = d[:, 1:]*h[:]

                L = scipy.sparse.spdiags(
                    d, [0, 1], nlyr, nlyr, format=mform).transpose()

                if out: print("L1 matrix variant -1 is "+str(numpy.shape(L)))

            elif variant == -2:
                alpha = 1.e-6
                d = numpy.zeros((2, nlyr))
                d[0, :-1] = 1.
                d[1, :-1] = -1.
                d[-1, -1] = alpha
                if der:
                    d[:, 1:] = d[:, 1:]*h[:]

                L = scipy.sparse.spdiags(
                    d, [0, 1], nlyr, nlyr, format=mform).transpose()

                if out: print("L1 matrix variant -3 is "+str(numpy.shape(L)))

            else:
                error("DiffOperator variant " + variant + " not implemeted ! Exit.")

    else:
        error("DiffOperator " + otype + " not implemeted ! Exit.")

    if mtype == "dense":
        L = L.todense()

    # print("\nRegularisation weight matrix dims:")
    # print(numpy.shape(L))
    return L


def matern(dist=numpy.array([]), var=numpy.array([]), nu=1):
    """
        Compute the Matern kernel function. nu determines the
        "smoothness".

        nu > inf:           gaussian (inf should not be larger 100!)
        nu = 1/2:           exponential

    VR Jul 10, 2022

    """

    if nu <= 0.:
        error("matern: input argument nu must be positive! Exit.")
    if nu > 100:
        nu = 100.
        print("matern: input argument nu too large! Set to 100.")

    dd = numpy.clip(dist.ravel, 1.e-8, None)

    f1 = numpy.power(2., nu-1.)/scipy.special.gamma(nu)
    f2 = dd * numpy.sqrt(2.*nu)

    matern = f1 * numpy.power(f2, nu) * scipy.special.kv(nu, f2)

    return matern


def covar(dx, dy, dz,
          covtype=["exp", numpy.array([1.0, 1.0, 1.0])],
          var=numpy.array([]), sparse=False, thresh=0.05,
          dist=False, inverse=True, OutInfo=True):
    """
    Calculate exponential, matern, or gaussian 3D covariance.

    Uses correlation lengths Lx,Ly,Lz.

    """
    z = set_zcenters(dz)
    x = numpy.zeros_like(z)
    y = numpy.zeros_like(z)
    shapez = numpy.shape(z)
    npoints = shapez[0]
    Cov = numpy.zeros((npoints, npoints))

    print(covtype[0].lower())
    if "exp" in covtype[0].lower():
        L = covtype[1]

        if OutInfo:
            print(str(npoints)+" exponential covariance")
        if dist:
            for il in range(npoints):
                for jl in range(il, npoints):
                    wdist = -(numpy.abs(x[il] - x[jl]) / L[0]
                              + numpy.abs(y[il] - y[jl]) / L[1]
                              + numpy.abs(z[il] - z[jl]) / L[2])
                    Cov[jl, il] = numpy.exp(wdist)
                    Cov[il, jl] = Cov[jl, il]
        else:
            for il in range(npoints):
                for jl in range(il, npoints):
                    ll = numpy.abs(il - jl)
                    wdist = -(ll/L[0] + ll/L[1] + ll / L[2])
                    Cov[jl, il] = numpy.exp(wdist)
                    Cov[il, jl] = Cov[jl, il]

    if "gau" in covtype[0].lower():
        L = covtype[1]
        if OutInfo:
            print(str(npoints)+" gaussian covariance")

        if dist:
            for il in range(npoints):
                for jl in range(il, npoints):
                    wdist = -(numpy.power(numpy.abs(x[il] - x[jl])/L[0], 2)
                              + numpy.power(numpy.abs(y[il] - y[jl])/L[1], 2)
                              + numpy.power(numpy.abs(z[il] - z[jl])/L[2], 2))
                    Cov[jl, il] = numpy.exp(wdist)
                    Cov[il, jl] = Cov[jl, il]
        else:
            for il in range(npoints):
                for jl in range(il, npoints):
                    ll = numpy.abs(il - jl)
                    wdist = -(numpy.power(ll/L[0], 2)
                              + numpy.power(ll/L[1], 2)
                              + numpy.power(ll/L[2], 2))
                    Cov[jl, il] = numpy.exp(wdist)
                    Cov[il, jl] = Cov[jl, il]

        # enforce positive definiteness
        threshpd = 1.0e2 * numpy.finfo(float).eps
        eigval, eigvec = numpy.linalg.eig(Cov)
        q = numpy.matrix(eigvec)
        xdiag = numpy.matrix(numpy.diag(numpy.maximum(eigval, threshpd)))
        Cov = q * xdiag * q.T

    if "mat" in covtype[0].lower():

        L = covtype[1]
        nu = covtype[2]

        if OutInfo:
            print(str(npoints)+" matern covariance")

        if nu <= 0.:
            error("matern: input argument nu must be positive! Exit.")
        if nu > 100:
            nu = 100.
            print("matern: input argument nu too large! Set to 100.")

        f1 = numpy.power(2., nu-1.)/scipy.special.gamma(nu)
        f2 = numpy.sqrt(2.*nu)

        if dist:

            for il in range(npoints):
                for jl in range(il, npoints):
                    wdist = numpy.sqrt(
                        numpy.power(numpy.abs(x[il] - x[jl])/L[0], 2)
                        + numpy.power(numpy.abs(y[il] - y[jl])/L[1], 2)
                        + numpy.power(numpy.abs(z[il] - z[jl])/L[2], 2)
                    )

                    if wdist <= 1.e-8:
                        wdist = 1.e-8

                    Cov[jl, il] = f1 * \
                        numpy.power(f2*wdist, nu) * \
                        scipy.special.kv(nu, f2*wdist)
                    Cov[il, jl] = Cov[jl, il]
        else:
            for il in range(npoints):
                for jl in range(il, npoints):
                    ll = numpy.abs(il - jl)
                    wdist = numpy.sqrt(numpy.power(ll/L[0], 2)
                                       + numpy.power(ll/L[1], 2)
                                       + numpy.power(ll/L[2], 2))

                    if wdist <= 1.e-8:
                        wdist = 1.e-8
                    Cov[jl, il] = f1 * \
                        numpy.power(f2*wdist, nu) * \
                        scipy.special.kv(nu, f2*wdist)
                    Cov[il, jl] = Cov[jl, il]

    print(numpy.shape(Cov))

    if numpy.size(var) > 0:

        Cov = numpy.dot(Cov, numpy.diag(var))

    if inverse:
        C = scipy.linalg.inv(Cov)
        SqrtC = msqrt(C, "cholesky")
    else:
        C = Cov
        SqrtC = msqrt(C, "cholesky")

    if sparse:

        threshC = thresh*numpy.amax(numpy.diag(C))
        threshS = thresh*numpy.amax(numpy.diag(SqrtC))

        C[numpy.abs(C) < threshC] = 0.0
        C = scipy.sparse.csr_matrix(C)

        SqrtC[numpy.abs(SqrtC) < threshS] = 0.0
        SqrtC = scipy.sparse.csr_matrix(SqrtC)

    return C, SqrtC


def msreg(dz=None, m=None,
          seps=1.e-4, otype="MS", diffop=None, mtype="sparse", mform="csr"):
    """
    Calculate minimum support weights for MS or MGS.
    Literature:

    VR Jan 2021
    """
    if otype == "MS":
        mabs = numpy.abs(m)
        W = 1. / (mabs + seps*seps)

    elif otype == "MGS":

        Lm = diffop * m

        L2 = Lm@Lm
        W = 1./(L2 + seps*seps)

    elif otype == "MSG":

        Lm = diffop * m

        L2 = Lm@Lm
        W = 1./(L2 + seps*seps)

#    W = scipy.linalg.inv(M)

    if mtype == "dense":
        W = W.todense()

    return W


def calc_regstart(D=numpy.array([]), M=numpy.array([]), Fac=1., out=True):
    """
    Estimate starting regularisation parameter value for
    methods like Occam and similar cooling schemes.

    see:

    Zhdanov, M. S.
        Geophysical inverse theory and regularization problems
        Elsevier, 2002
    Zhdanov, M. S.
        Inverse theory and applications in geophysics
        Elsevier, 2015

    Parameters
    ----------
    D : float
        Data fit norm.
    M : float
       Model norm.

    Returns
    -------
    RegParBase: float
        Estimate of

    """
    if D.size == 0 or M.size == 0:
        error("cal_regstart: D or M are 0! Exit.")

    taustart = numpy.nan

    if numpy.size(D) > 1:
        if scipy.sparse.issparse(D):
            D0 = scipy.sparse.linalg.norm(D)
        else:
            D0 = scipy.linalg.norm(D)
    else:
        D0 = D

    if numpy.size(M) > 1:
        if scipy.sparse.issparse(M):
            M0 = scipy.sparse.linalg.norm(M)
        else:
            M0 = scipy.linalg.norm(M)
    else:
        M0 = M

    taustart = Fac * D0/M0
    if out:
        print("Initial Tau = ", taustart)

    return taustart


def calc_regstart_base(J=numpy.array([]), W=numpy.array([]), Fac=1., out=True):
    """
    Get starting tau for Occam/Cooling inversion

    Inputs:
    -----------
    Jacd:       scaled Jacobian
    W:          regularization matrix, diff_op or sqrt (C_d)^-1

    Returns:
    -----------
    taustart

    Reference:

    Chen, T. & Yang, D.
    Modeling and Inversion of Airborne and Semi-Airborne Transient
    Electromagnetic Data with Inexact Transmitter and Receiver Geometries
    Remote Sensing, 2022, 14, doi:10.3390/rs14040915


    vr feb 20, 2023
    """
    if J.size == 0 or W.size == 0:
        error("cal_regstart: J or D are 0! Exit.")

    taustart = numpy.nan

    # S1 = Jacd.T@Jacd
    t1 = numpy.amax(scipy.linalg.svd(J.T@J, compute_uv=False))
    # S2 = W.T@W
    t2 = numpy.amax(scipy.linalg.svd(W.T@W, compute_uv=False))

    taustart = Fac*t1/t2

    return taustart


def rademacher_sample(N=None):
    """
    Calculate N samples drawn from a Rademacher distribution.

    Parameters
    ----------
    N: int
        Number of samples drawn. The default is None.

    Returns:float
        array of N samples from Rademacher distribution

    Author:VR 12/20

    """
    tmp = numpy.random.default_rng().uniform(0.0, 1.0, N)
    v = numpy.ones_like(tmp)
    v[tmp < 0.5] = -1.0

    return v

def msqrt_sparse(M=numpy.array([]), smallval=1.e-12):
    """
    Calculate sparse Cholesky.

    Missing in scipy.

    Parameters
    ----------
    A : double
        Positive definite sparse matrix.

    Returns
    -------
    CholA: double
        Cholesky factor of A.

    VR Feb 2021

    """
    n =M.shape[0]
    MM = M.copy() + numpy.identity(n)*smallval

    LU = scipy.sparse.linalg.splu(
        MM, diag_pivot_thresh=0)  # sparse LU decomposition

    # check the matrix A is positive definite.
    if (LU.perm_r == numpy.arange(n)).all() and (LU.U.diagonal() > 0).all():
        SqrtM = LU.L.dot(scipy.sparse.diags(LU.U.diagonal() ** 0.5))

    else:
        error("The matrix is not positive definite")

    return SqrtM


def msqrt(M=numpy.array([]), method="cho", smallval=1.e-12):
    """
    Computes a matrix square-root (Choleky, or eig).

    Parameter:
    M: M is a positive Hermitian (or positive definite) matrix.

    Return:
    SqrtM, Mevals, Mevecs:
    Here, SqrtM is a matrix such that SqrtM * SqrtM.T = M.
    The vector Mevals contains the eigenvectors of M,
    and the matrix Mevecs the corresponding eigenvectors.

    Also Calculate sparse Cholesky, missing in scipy.

    Parameters
    ----------
    A : double
        Positive definite sparse matrix.
    smallval: double
        small value to guarantee positive definiteness in
        the case of numueerical noise.
    method: str
        eigenvalue or cholesky in case of dense input matrices

    Returns
    -------
    CholM: double
        Cholesky factor of A.

    Last change: VR Mar 2024


    """
    n = numpy.shape(M)[0]
    MM = M.copy() + numpy.identity(n)*smallval

    if "eig" in method.lower():
        # compute eigenvalues and eigenvectors
        Mevals, Mevecs = scipy.linalg.eigh(MM)
        Mevals = Mevals.clip(min=0.0)
        SqrtM = Mevecs * numpy.sqrt(Mevals)
        return SqrtM, Mevals, Mevecs

    if "cho" in method.lower():
        SqrtM = scipy.linalg.cholesky(MM)


    return SqrtM

def isspd(A):

    n = A.shape[0]

    AAT = A@A.T
    if numpy.allclose(AAT, numpy.identity(n), rtol = 1.e-8, atol=1.e-8):
        print("A is symmetric.")
    else:
        print("A is NOT symmetric.")

    spd = numpy.all(numpy.linalg.eigvals(A) > 1.e-12)


    return spd

def rsvd(A, rank=300,
         n_oversamples=300,
         n_subspace_iters=None,
         return_range=False):
    """
    =============================================================================
    Randomized SVD. See Halko, Martinsson, Tropp's 2011 SIAM paper:

    "Finding structure with randomness: Probabilistic algorithms for constructing
    approximate matrix decompositions"
    Author: Gregory Gundersen, Princeton, Jan 2019
    =============================================================================
    Randomized SVD (p. 227 of Halko et al).

    :param A:                (m x n) matrix.
    :param rank:             Desired rank approximation.
    :param n_oversamples:    Oversampling parameter for Gaussian random samples.
    :param n_subspace_iters: Number of power iterations.
    :param return_range:     If `True`, return basis for approximate range of A.
    :return:                 U, S, and Vt as in truncated SVD.
    """
    if n_oversamples is None:
        # This is the default used in the paper.
        n_samples = 2 * rank
    else:
        n_samples = rank + n_oversamples

    # Stage A.
    # print(' stage A')
    Q = find_range(A, n_samples, n_subspace_iters)

    # Stage B.
    # print(' stage B')
    B = Q.T @ A
    # print(numpy.shape(B))
    # print(' stage B before linalg')
    U_tilde, S, Vt = numpy.linalg.svd(B)
    # print(' stage B after linalg')
    U = Q @ U_tilde

    # Truncate.
    U, S, Vt = U[:, :rank], S[:rank], Vt[:rank, :]

    # This is useful for computing the actual error of our approximation.
    if return_range:
        return U, S, Vt, Q
    return U, S, Vt


# ------------------------------------------------------------------------------


def find_range(A, n_samples, n_subspace_iters=None):
    """Algorithm 4.1: Randomized range finder (p. 240 of Halko et al).

    Given a matrix A and a number of samples, computes an orthonormal matrix
    that approximates the range of A.

    :param A:                (m x n) matrix.
    :param n_samples:        Number of Gaussian random samples.
    :param n_subspace_iters: Number of subspace iterations.
    :return:                 Orthonormal basis for approximate range of A.
    """
    # print('here we are in range-finder')
    rng = numpy.random.default_rng()

    m, n = A.shape
    # print(A.shape)
    O = rng.normal(size=(n, n_samples))
    # print(O.shape)
    Y = A @ O

    if n_subspace_iters:
        return subspace_iter(A, Y, n_subspace_iters)
    else:
        return ortho_basis(Y)


# ------------------------------------------------------------------------------


def subspace_iter(A, Y0, n_iters):
    """Algorithm 4.4: Randomized subspace iteration (p. 244 of Halko et al).

    Uses a numerically stable subspace iteration algorithm to down-weight
    smaller singular values.

    :param A:       (m x n) matrix.
    :param Y0:      Initial approximate range of A.
    :param n_iters: Number of subspace iterations.
    :return:        Orthonormalized approximate range of A after power
                    iterations.
    """
    # print('herere we are in subspace-iter')
    Q = ortho_basis(Y0)
    for _ in range(n_iters):
        Z = ortho_basis(A.T @ Q)
        Q = ortho_basis(A @ Z)
    return Q


# ------------------------------------------------------------------------------


def ortho_basis(M):
    """Computes an orthonormal basis for a matrix.

    :param M: (m x n) matrix.
    :return:  An orthonormal basis for M.
    """
    # print('herere we are in ortho')
    Q, _ = numpy.linalg.qr(M)
    return Q


def project_nullspace(U=numpy.array([]), m_test=numpy.array([])):
    """
    Calculates nullspace projection of a vector

    Parameters
    ----------
    U : numpy array, float
         npar*npar matrix from SAVD oj Jacobian.
    m_test : numpy array, float
         npar*vector to be projected.

    Returns
    -------
    m: numpy array, float
        projected model

    """
    if numpy.size(U) == 0:
        error("project_nullspace: V not defined! Exit.")

    m_proj = m_test - U@(U.T@m_test)

    return m_proj

# def project_model(m=None, U=None, small=1.0e-14, out=True):
#     """
#     Project to Nullspace.

#     (see Munoz & Rath, 2006)
#     author: vrath
#     last changed: Sep 25, 2020
#     """
#     b = numpy.dot(U.T, m)
#     # print(m.shape)
#     # print(b.shape)
#     # print(U.shape)

#     mp = m - numpy.dot(U, b)

#     return mp


def calc_upr(dnorm=numpy.array([]), M=numpy.array([]), err=numpy.array([])):
    """
    Calculates UPRE.

    Parameters
    ----------
    dnorm=None              data residual norm
    M=None                J times generalized inverse, R-Matrix
    err=None              data error

    Returns
    -------
    upr_val

    see:

        Vogel, C. R.:
        Computational Methods for Inverse Problems,
        Society for Industrial Mathematics, 2002

        Lin, Youzuo, Wohlberg, Brendt
        Application of the UPRE Method to Optimal Parameter Selection
        for Large Scale Regularization Problems
        IEEE Southwest Symposium on Image Analysis and Interpretation
        Santa Fe, NM, USA, p. 89-92, 2008

    VR May 2022
    """
    if (numpy.size(dnorm) == 0) or (numpy.size(M) == 0) or (numpy.size(err) == 0):
        error("calc_upr: parameters missing! Exit.")

    nd = numpy.size(err)
    traceM = numpy.trace(M)
    sqerr = numpy.linalg.alg.norm(err, 2)
    upr_val = numpy.power(dnorm, 2)/nd + (2.*sqerr/nd)*traceM - sqerr

    return upr_val


def calc_mle(M=numpy.array([]),
             d=numpy.array([])):
    """
    Calculates MLE (Maximum Likelihood Estimate).

    Parameters
    ----------
    d = None              data
    M = None              Generalized inverse times J^T, R-Matrix

    Returns
    -------
    mle_val

    see:

        J. Xu, F. Schreier, A. Doicu, and T. Trautmann:
        Assessment of Tikhonov-type regularization methods for
        solving atmospheric inverse problems
        Journal of Quantitative Spectroscopy and Radiative Transfer
        184, 274–286, 2016, doi:10.1016/j.jqsrt.2016.08.003.

        A. Doicu, T. Trautmann, and F. Schreier
        Numerical regularization for atmospheric inverse problems.
        Berlin, DE: Springer, 2010.

        J. Xu, R. Lanlan, F. Schreier, D. Efremenko, A. Doicu, and T. Trautmann
        Insight into Construction of Tikhonov-Type Regularization for
        Atmospheric Retrievals
        Atmosphere, 11, p. 1052, 2020, doi:10.3390/atmos11101052.



    VR Apr 2022
    """
    if (numpy.size(M) == 0) or (numpy.size(d) == 0):
        error("calc_mle: parameters missing! Exit.")

    nd = numpy.size(d)
    m1 = (d*numpy.trace(numpy.eye(nd) - M)*d.T)
    m2 = numpy.power(numpy.det(numpy.eye(nd) - M), 1./nd)
    mle_val = m1/m2

    return mle_val


def calc_lc_corner(dnorm=numpy.array([]), mnorm=numpy.array([])):
    """
    Calculates corner of thhe L-curve.

    Parameters
    ----------
    dnorm                   data norm
    mnorm                   Generalized inverse times J^T

    Returns
    -------
    lcc_val                 value of gcv function)

    see:

        Per Christian Hansen:
        Discrete Inverse Problems: Insight and Algorithms
        SIAM, Philadelphia, 2010

        Per Christian Hansen:
        The L-Curve and its Use in the Numerical Treatment of Inverse Problems
        In: P. Johnston ,Computational Inverse Problems in Electrocardiology
        WIT Press, 2001
        119-142

        Per Christian Hansen:
        Rank Deficient and Discrete Ill-Posed Problems
        SIAM, Philadelphia, 1998

    VR June 2022
    """
    if (numpy.size(dnorm) == 0) or (numpy.size(mnorm) == 0):
        error("calc_lcc: parameters missing! Exit.")

    lcurvature = curvature(numpy.log(dnorm), numpy.log(mnorm))

    indexmax = numpy.argmax(lcurvature)

    return indexmax


def calc_gcv(dnorm=numpy.array([]), M=numpy.array([]), a=1.):
    """
    Calculates (modified) GCV.

    Parameters
    ----------
    dnorm=None              data norm
    M=None                  Generalized inverse times J^T
    a=1.                    Value of modificator parameter (optional)

    Returns
    -------
    gcv_val                 value of gcv function)

    see:

    Golub,G. H., Heath, M. T. & Wahba, G.:
    Generalized Cross-Validation as a Method for Choosing a Good Ridge Parameter
    Technometrics, 21, pp. 215–223, 1979.

    Farquharson C. G. & Oldenburg D. W.:
    A comparison of automatic techniques for estimating the regularization
    parameter in non-linear inverse problems
    Geophysical Journal International, 156, Art. no. 3, 2004.

    Bauer, F. & Lukas, M. A.:
    Comparing parameter choice methods for regularization of ill-posed problems
    Mathematics and Computers in Simulation, 2011
    doi:10.1016/j.matcom.2011.01.016

    VR May 2021
    """
    if (numpy.size(dnorm) == 0) or (numpy.size(M) == 0):
        error("calc_gcv: parameters missing! Exit.")

    nd = numpy.shape(M)[0]
    gcv_val = (nd * numpy.power(dnorm, 2)
               / numpy.power(numpy.trace(numpy.eye(nd) - a*M), 2))

    return gcv_val


def calc_ufc(dnorm=numpy.array([]),
             mnorm=numpy.array([])):
    """
    Calculates U-Curve.

    Parameters
    ----------
    dnorm=None              data norm
    mnorm=None              model norm

    Returns
    -------
    u_val                  U-curve value

    see:

    Krawczyk-Stando, D. & Rudnicki, M.:
    Regularization Parameter Selection in Discrete Ill–posed
    Problems — the Use of the U–curve
    Int. J. Appl. Math. Comput. Sci., 2007, 17, 157-164

    VR March 2023
    """
    if numpy.size(dnorm) == 0 or numpy.size(mnorm) == 0:
        error("No model or data norms given! Exit.")

    sqdnorm = numpy.power(dnorm, 2)
    sqmnorm = numpy.power(mnorm, 2)

    u_val = numpy.log(1./sqdnorm + 1./sqmnorm)

    return u_val


def curvature(x_data, y_data):
    """
    Calculates curvature for all interior points
    on a curve whose coordinates are provided
    Used for l-curve corner estimation.
    Input:
        - x_data: list of n x-coordinates
        - y_data: list of n y-coordinates
    Output:
        - curvature: list of n-2 curvature values

    originally written by Hunter Ratliff on 2019-02-03
    """
    curvature = []
    for i in range(1, len(x_data)-1):
        R = circumradius(x_data[i-1:i+2], y_data[i-1:i+2])
        if (R == 0):
            print("Failed: points are either collinear or not distinct")
            return 0
        curvature.append(1/R)
    return curvature


def circumradius(xvals, yvals):
    """
    Calculates the circumradius for three 2D points

    originally written by Hunter Ratliff on 2019-02-03
    """
    x1, x2, x3, y1, y2, y3 = xvals[0], xvals[1], xvals[2], yvals[0], yvals[1], yvals[2]
    den = 2.*((x2-x1)*(y3-y2)-(y2-y1)*(x3-x2))
    num = ((((x2-x1)**2) + ((y2-y1)**2))
           * (((x3-x2)**2)+((y3-y2)**2))
           * (((x1-x3)**2)+((y1-y3)**2)))**(0.5)
    if (den == 0.):
        print("Failed: points are either collinear or not distinct")
        return 0.
    R = abs(num/den)

    return R


def circumcenter(xvals, yvals):
    """
    Calculates the circumcenter for three 2D points

    originally written by Hunter Ratliff on 2019-02-03
    """
    x1, x2, x3, y1, y2, y3 = xvals[0], xvals[1], xvals[2], yvals[0], yvals[1], yvals[2]
    A = 0.5*((x2-x1)*(y3-y2)-(y2-y1)*(x3-x2))
    if (A == 0):
        print("Failed: points are either collinear or not distinct")
        return 0
    xnum = ((y3 - y1)*(y2 - y1)*(y3 - y2)) - \
        ((x2**2 - x1**2)*(y3 - y2)) + ((x3**2 - x2**2)*(y2 - y1))
    x = xnum/(-4*A)
    y = (-1*(x2 - x1)/(y2 - y1))*(x-0.5*(x1 + x2)) + 0.5*(y1 + y2)
    return x, y


def calc_dnorm(data_obs=numpy.array([]),
               data_cal=numpy.array([]),
               data_act=numpy.array([]),
               data_err=numpy.array([]),
               p=2, calc_res=False):
    """
    Calculate the p-norm of the residuals.

    VR Jan 2021

    """

    if (numpy.size(data_obs) == 0 or
        numpy.size(data_cal) == 0 or
            numpy.size(data_err) == 0):
        error("calc_dnorm: parameters missing! Exit.")

    if numpy.size(data_act) == 0:
        dat_obs = data_obs
        dat_cal = data_cal
        dat_err = data_err
    else:
        dat_obs = extract_dat(D=data_obs, d_act=data_act)
        dat_cal = extract_dat(D=data_cal, d_act=data_act)

    if numpy.ndim(data_err) == 1:
        dat_err = extract_dat(D=data_err, d_act=data_act)
        w = 1./dat_err
    else:
        w = 1./dat_err

    resid = w * (dat_obs - dat_cal)

    rnorm = numpy.linalg.norm(resid, p)

    if calc_res:
        return rnorm, resid
    else:
        return rnorm


def calc_datafit(data_obs=numpy.array([]),
             data_cal=numpy.array([]),
             data_err=numpy.array([]),
             data_act=numpy.array([]),
             smaptype=1):
    """
    Calculate the NRMS ans SRMS.

    See:
    (SRMS)
    G. A. Wilson and A. P. Raiche and F. Sugeng
    2.5D Inversion of airborne electromagnetic data
    Exploration Geophysics,37, 363-371, 2006

    VR  Jan 2021
        Jun 2022
        Aug 2024
    """

    if (numpy.size(data_obs) == 0 or
        numpy.size(data_cal) == 0 or
            numpy.size(data_err) == 0):
        error("calc_datafit: parameters missing! Exit.")

    if numpy.size(data_act) == 0:
        dat_obs = data_obs
        dat_cal = data_cal
        dat_err = data_err
    else:
        dat_obs = extract_dat(D=data_obs, d_act=data_act)
        dat_cal = extract_dat(D=data_cal, d_act=data_act)

    if numpy.ndim(data_err) == 1:
        dat_err = extract_dat(D=data_err, d_act=data_act)
        w = 1./dat_err
    else:
        w = 1./dat_err

    nd = numpy.size(dat_cal)
    rscal = w*(dat_obs - dat_cal)
    # normalized root mean square error
    nrmse = numpy.sqrt(numpy.sum(numpy.power(abs(rscal), 2)) / nd)

    # sum squared scaled symmetric error (Raiche2007)
    if smaptype==1:
        fac = 1.
        d0 = numpy.abs(dat_obs - dat_cal)
        d1 = (numpy.abs(dat_obs) + numpy.abs(dat_cal))/fac
        smape0 = numpy.sum(d0/d1)
        smape = 100.0 * smape0 / nd

    if smaptype==2:
        fac = 2.
        d0 = numpy.abs(dat_obs - dat_cal)
        d1 = (numpy.abs(dat_obs) + numpy.abs(dat_cal))/fac
        smape0 = numpy.sum(d0/d1)
        smape = 100.0 * smape0 / nd

    if smaptype==3:
        fac = 1.
        d0 = numpy.abs(dat_obs - dat_cal)
        d1 = (numpy.abs(dat_obs) + numpy.abs(dat_cal))/fac
        smape0 = numpy.sum(d0)/numpy.sum(d1)
        smape = 100.0 * smape0 / nd


    return nrmse, smape


def set_errors(data_obs=numpy.array([]),
               daterr_add=0., daterr_mult=0., perturb=False, pdf=["gauss", 0., 0.]):
    """
    Generate errors.

    Error model including multiplicative and additive noise
    following Brodie (2015) GALEISBSTDEM and Green & Lane (2003)


    A. Green and R. Lane,
    Estimating Noise Levels in AEM Data,ASEG 2003.
    R. Brodie,
    GALEISBSTDEM: A deterministic algorithm for 1D sample by sample
    inversion of time-domain AEM data — theoretical details,
    Geoscience Australia, 2015.


    VR Apr 2021

    """
    data_obs = data_obs.reshape(numpy.size(data_obs), 1)

    if daterr_add == 0. and daterr_mult == 0:
        error("set_errors: additive and multiplicative error is zero! Exit.")

    daterr_a = daterr_add * numpy.ones_like(data_obs)
    daterr_m = daterr_mult * numpy.ones_like(data_obs)

    data_err = \
        numpy.sqrt(numpy.power(daterr_m * data_obs, 2) +
                   numpy.power(daterr_a, 2))

    # print(numpy.shape(data_err))

    if perturb:
        data_perturb = perturb_data(data_obs=data_obs, pdf=["gauss", data_err])

    else:
        data_perturb = data_obs

    return data_err.flatten(), data_perturb.flatten()


def perturb_data(data_obs=numpy.array([]), pdf=["gauss", 0.,]):
    """
    Generate errors.

    Error model including multiplicative and additive noise
    following Brodie (2015) GALEISBSTDEM and Green & Lane (2003)


    A. Green and R. Lane,
    Estimating Noise Levels in AEM Data,ASEG 2003.
    R. Brodie,
    GALEISBSTDEM: A deterministic algorithm for 1D sample by sample
    inversion of time-domain AEM data — theoretical details,
    Geoscience Australia, 2015.


    VR Apr 2021

    """
    data_obs = data_obs.reshape(numpy.size(data_obs), 1)

    if "gau" in pdf[0].lower() or "nor" in pdf.lower():
        perturb = pdf[1] * numpy.random.standard_normal(numpy.shape(data_obs))
    if "uni" in pdf[0].lower():
        perturb = numpy.random.uniform(-pdf[1], +pdf[1], numpy.shape(data_obs))

    data_perturb = data_obs + perturb

    return data_perturb


def set_zcenters(dz):
    """
    Define cell centers.

    VR Jun 2022

    """
    nlyr = numpy.shape(dz)[0]

    dz = numpy.append(dz, dz[-1])
    z = numpy.append(0.0, numpy.cumsum(dz))
    # print(numpy.shape(z))
    zc = 0.5 * (z[0:nlyr] + z[1:nlyr+1])

    return zc


def set_znodes(dz):
    """
    Define vertical nodes (depths).

    VR Jan 2021

    """
    zn = numpy.append(0.0, numpy.cumsum(dz))

    return zn


def get_taustart(Jacd=None, W=None, out=True):
    """
    Get starting tau for Occam  or Cooling inversion

    Inputs:
    -----------
    Jacd:       scaled Jacobian
    W:          regularization matrix, diff_op or sqrt (C_d)^-1

    Returns:
    -----------
    taustart

    Reference:

    Chen, T. & Yang, D.
    Modeling and Inversion of Airborne and Semi-Airborne Transient
    Electromagnetic Data with Inexact Transmitter and Receiver Geometries
    Remote Sensing, 2022, 14, doi:10.3390/rs14040915


    vr feb 20, 2023
    """
    S1 = Jacd.T@Jacd
    t1 = numpy.amax(scipy.linalg.svd(S1, compute_uv=False))
    S2 = W.T@W
    t2 = numpy.amax(scipy.linalg.svd(S2, compute_uv=False))

    taustart = t1/t2

    return taustart


def impute_matrix_isvd(Y, k=None, tol=1E-3, maxiter=10):
    """
    Approximate SVD on data with missing values via expectation-maximization

    Inputs:
    -----------
    Y:          (nobs, ndim) data matrix, missing values denoted by NaN/Inf
    k:          number of singular values/vectors to find (default: k=ndim)
    tol:        convergence tolerance on change in trace norm
    maxiter:    maximum number of EM steps to perform (default: no limit)

    Returns:
    -----------
    Y_hat:      (nobs, ndim) reconstructed data matrix
    mu_hat:     (ndim,) estimated column means for reconstructed data
    U, s, Vt:   singular values and vectors (see numpy.linalg.svd and
                scipy.sparse.linalg.svds for details)


    vr oct 29, 2022
    """

    if k is None:
        svdmethod = functools.partial(numpy.linalg.svd, full_matrices=False)
    else:
        svdmethod = functools.partial(scipy.sparse.linsvds, k=k)

    if maxiter is None:
        maxiter = numpy.inf

    # initialize the missing values to their respective column means
    mu_hat = numpy.nanmean(Y, axis=0, keepdims=1)
    valid = numpy.isfinite(Y)
    Y_hat = numpy.where(valid, Y, mu_hat)

    halt = False
    ii = 1
    v_prev = 0

    while not halt:

        # SVD on filled-in data
        U, s, Vt = svdmethod(Y_hat - mu_hat)

        # impute missing values
        Y_hat[~valid] = (U.dot(numpy.diag(s)).dot(Vt) + mu_hat)[~valid]

        # update bias parameter
        mu_hat = Y_hat.mean(axis=0, keepdims=1)

        # test convergence using relative change in trace norm
        v = s.sum()
        if ii >= maxiter or ((v - v_prev) / v_prev) < tol:
            halt = True
        ii += 1
        v_prev = v

    return Y_hat, mu_hat, U, s, Vt


def calc_mad(datavec=numpy.array([]), median=None, Out=False):
    """
    Calculate Median Absolute Deviation (MAD)

    Parameters
    ----------
    d_vec : array, float
        Vector for which MAD is calculated..
    median : float, optional
        If None, median is calulated.

    Returns
    -------
    MAD

    """
    if numpy.size(datavec) == 0:
        error("find_nearest: No vector given! exit.")

    if median == None:
        median = numpy.nanmedian(datavec)

    d = numpy.abs(datavec - median)
    mad = numpy.nanmedian(d)

    if Out:
        print("MAD = "+str(mad))

    return mad, median


def calc_stat_ens(ensemble=numpy.array([]),
                  quantiles=[2.3, 15.9,],
                  qmode="perc",
                  sum_stats=False):
    """
    calculates basic statistics  for ensemble

    Parameters
    ----------
    ensemble : numpy.array
        Ensemble of models or data sets. The default is numpy.array([]).
    quantiles : list of floats, optional
        quantiles or percentiles. The default is [2.3, 15.9,] in percent
    inmode : strng, optional
        switch between percentiles and quantiles. The default is "perc".

    Returns
    -------
    if sum_stat = False:
        median, mupp, mlow : numpy.array of floats
            Median snd symmetric quantiles corresponding to quantiles.
    else:
        additional mean, standard dev, skew, kurtosis, mode

    Last change: vrath, May 18, 2023

    """
    if ensemble.size == 0:
        error("calc_stat_ens: ensemble not defined")

    if ("quant" in qmode.lower()) and numpy.any(quantiles > 1.):
        inmode = "perc"

    quants = numpy.median(ensemble, axis=0)

    if "quant" in qmode.lower():
        for per in numpy.arange(len(quantiles)):
            qc = [quantiles[per], 1.-quantiles[per]]
            quants = numpy.vstack(
                [quants, numpy.quantile(ensemble, qc, axis=0)])
    else:
        for per in numpy.arange(len(quantiles)):
            qc = [quantiles[per], 100.-quantiles[per]]
            quants = numpy.vstack(
                [quants, numpy.percentile(ensemble, qc, axis=0)])

    if sum_stats:
        mean = numpy.mean(ensemble, axis=0)
        stdv = numpy.std(ensemble, axis=0)
        mode = scipy.stats.mode(ensemble, axis=0)
        kurt = scipy.stats.kurtosis(ensemble, axis=0)
        skew = scipy.stats.skew(ensemble, axis=0)
        return quants, mean, stdv, skew, kurt, mode

    else:
        return quants


def calc_made(ensemble=numpy.array([]), median=None, Out=False):
    """
    Calculate Median Absolute Deviation (MAD)

    Parameters
    ----------
    ensemble : array, float
        nd array for which MAD is calculated along first axis.
    median : float, optional
        If None, median is calulated.

    Returns
    -------
    MAD

    """
    if numpy.size(ensemble) == 0:
        error("find_nearest: No vector uiven! exit.")

    if median == None:
        median = numpy.nanmedian(ensemble, axis=0)

    d = numpy.abs(ensemble - median)
    mad = numpy.nanmedian(d)

    if Out:
        print("MAD = "+str(mad))

    return mad, median


def dctn(x, norm="ortho"):
    """
    Discrete cosine transform (fwd)
    https://stackoverflow.com/questions/13904851/use-pythons-scipy-dct-ii-to-do-2d-or-nd-dct
    """
    import scipy.fftpack
    for i in range(x.ndim):
        x = scipy.fftpack.dct(x, axis=i, norm=norm)
    return x


def idctn(x, norm="ortho"):
    """
    Discrete cosine transform (inv)
    https://stackoverflow.com/questions/13904851/use-pythons-scipy-dct-ii-to-do-2d-or-nd-dct
    """
    import scipy.fftpack

    for i in range(x.ndim):
        x = scipy.fftpack.idct(x, axis=i, norm=norm)
    return x


def KLD(P=numpy.array([]), Q=numpy.array([]), epsilon=1.e-8):
    """
    Calculates Kullback-Leibler distance

    Parameters
    ----------
    P, Q: numpy.array
        pdfs
    epsilon : TYPE
        Epsilon is used here to avoid conditional code for
        checking that neither P nor Q is equal to 0.

    Returns
    -------

    distance: float
        KL distance


    """
    if P.size * Q.size == 0:
        error("KLD: P or Q not defined! Exit.")

    # You may want to instead make copies to avoid changing the np arrays.
    PP = P.copy()+epsilon
    QQ = Q.copy()+epsilon

    distance = numpy.sum(PP*numpy.log(PP/QQ))

    return distance


def sample_pcovar(cpsqrti=None, m=None, tst_sample=None,
                  nsamp=1, small=1.0e-14, out=True):
    """
    Sample Posterior Covariance.
    Algorithm given by  Osypov (2013)

    Parameters
    ----------

    Returns
    -------
    spc_sanple

    References:

    Osypov K, Yang Y, Fournier A, Ivanova N, Bachrach R,
    Can EY, You Y, Nichols D, Woodward M (2013)
    Model-uncertainty quantification in seismic tomography: method and applications
    Geophysical Prospecting, 61, pp. 1114–1134, 2013, doi: 10.1111/1365-2478.12058.


    """
    error("sample_pcovar: Not yet fully implemented! Exit.")

    if (cpsqrti == None) or (m == None):
        error("sample_pcovar: No covariance or ref model given! Exit.")

    if m.ndim(m) > 1:
        m = m.flatten(order="F")

    if tst_sample == None:
        print("sample_pcovar: "+str(nsamp)+" sample models will be generated!")
        if nsamp == 0:
            error("sample_pcovar: No number of samples given! Exit.")
        tst_sample = numpy.random.default_rng().normal(0., 1., (nsamp, len(m)))

    else:
        nsamp = numpy.shape(tst_sample)[0]

    spc_sample = numpy.zeros(nsamp, len(m))

    for isamp in numpy.arange(nsamp):
        spc_sample[isamp, :] = m + cpsqrti@tst_sample[isamp, :]

    return spc_sample

def init_layers(nlyr=26, start=1., end=10.,
               logspace=True, out=True):
    """
    Set layer Parameters

    Parameters
    ----------
    nlyr : integer, optional
        Number of layers. The default is 26.
    start : float, optional
        Thickness of first layer. The default is 1 m.
    end : float, optional
        Largest thickness. The default is 10 m.
    logspace : logical, optional
        If True logarithmically spaced dz are generated. The default is True.
    out : logical, optional
        print output. The default is True.

    Returns
    -------
    dz : numpy array
        Layer thicknesses (nlyr).
    z_node : numpy array
        Node coordinates (nlyr+1).
    z_cent : numpy array
       Layer center coordinates (nlyr+1).
    """

    if logspace:
        dz = numpy.logspace(numpy.log10(start), numpy.log10(end), nlyr)
    else:
        dz = numpy.linalg.space(start, end, nlyr)

    z_node = numpy.append(0.0, numpy.cumsum(dz))

    dz_tmp = numpy.append(dz, dz[-1])
    z_tmp = numpy.append(0.0, numpy.cumsum(dz_tmp))
    z_cent = 0.5 * (z_tmp[0:nlyr] + z_tmp[1:nlyr+1])

    if out:
        print("\n")
        if logspace:
            print(" Number of log-spaced layers:", nlyr)
        else:
            print(" Number of lin-spaced (possibly equidistant) layers:", nlyr)
        print("Number of nodes:", nlyr+1)
        print("Number of cell centers:", nlyr+1)

    return dz, z_node, z_cent

def set_prior(pval=numpy.array([]), flightline=None):
    """
    set prior according to different methods

    """
    if flightline is None:
        error ("set_prior: no flightine given! Exit.")


    dims =numpy.shape(flightline)
    prior = numpy.zeros(dims)
    dims =numpy.shape(flightline)

    if numpy.size()==1:
        prior = prior + pval[0]

    if numpy.shape(pval)==numpy.shape(flightline):
        prior = pval
    elif numpy.shape(pval)==numpy.shape(flightline[:, 0]):
        for ii in numpy.arange(dims[1]):
            prior[:, ii] = pval
    else:
        for ii in numpy.arange(dims[1]):
            prior[:, ii] = pval

    return prior

def merge_data_sets(infile_list=None, outfile_name="./dat_tmp.npz",
                    aem_system="aem05", dictout=False, out=False):
    """
    Merge data sets from file list
    Created on Sun Jan 22 12:18:58 2023
    @author: vrath
    """
    dateform="%m/%d/%Y, %H:%M:%S"
    _,NN, _, _, _, = aesys.get_system_params(System=aem_system)

    k = 0
    for infile in infile_list:
        k = k+1
        if out: print("\nData read from: %s" % infile)
        data_k, header, _ = aesys.read_aempy(File=infile, System=aem_system, OutInfo=False)
        if k == 1:
            merged_data = data_k
        else:
            merged_data = numpy.vstack((merged_data, data_k))
    if outfile_name is not None:
        header = "merged data set:"+"".join("Date " + datetime.now().strftime(dateform))
        aesys.write_aempy(File=outfile_name, Data=merged_data, System=aem_system,
                Header=header, OutInfo=out)

    if dictout:
        merged_data = {
                "fnam": merged_data[:,0],
                "x": merged_data[:,1],
                "y": merged_data[:,2],
                "gps": merged_data[:,3],
                "alt": merged_data[:,4],
                "dem": merged_data[:,5],
                "data": merged_data[:,6:6+NN[2]]
                        }

    return merged_data

def merge_model_sets(infile_list=None, outfile_name="./mod_tmp.npz",
                     qthresh=None, out=True):
    """
    Merge models from file list
    Created on Sun Jan 22 12:18:58 2023
    @author: vrath
    """
    dateform="%m/%d/%Y, %H:%M:%S"
    # _,NN, _, _, _, = aesys.get_system_params(System=aem_system)
    if outfile_name is not None:
        if not ".npz" in os.path.splitext(outfile_name)[1]:
            error("merge_model_sets: Only npz format implemented.! Exit.")

    k = 0
    all_sites = 0
    chs_sites = 0
    for infile in infile_list:
        k=k+1
        print("\nData read from: %s" % infile)
        results = numpy.load(infile, allow_pickle=True)

        # ctrl = numpy.load(file.replace("_results.npz","_ctrl.npz"))
        # ctrl = results["ctrl"]

        if k==1 and outfile_name is not None :
            ctrl_file_in = infile.replace("_results.npz","_ctrl.npz")
            Ctrl =  numpy.load(ctrl_file_in, allow_pickle=True)
            ctrl_file_out = outfile_name.replace("_results.npz", "_ctrl.npz")
            numpy.savez_compressed(file=ctrl_file_out, **Ctrl)



        print(list(results.keys()))

        mod_ref = results["mod_ref"]
        mod_act = results["mod_act"]

        site_mod = results["site_modl"]
        site_sns = results["site_sens"]
        site_rms = results["site_nrms"]
        site_smp = results["site_smap"]
        site_con = results["site_conv"]

        site_y = results["site_y"]
        site_x = results["site_x"]
        site_gps = results["site_gps"]
        site_alt = results["site_alt"]
        site_dem = results["site_dem"]
        site_cov = results["site_pcov"]

        # tmp = results["site_log"]
        # print(numpy.shape(tmp))
        # print(tmp)

        # site_log = numpy.array([tmp[1], tmp[2], tmp[3]])

        site_d = numpy.zeros_like(site_mod)
        site_z = numpy.zeros_like(site_mod)

        nlyr = get_nlyr(mod_ref)
        tmp = numpy.cumsum(mod_ref[6*nlyr:7*nlyr-1])
        tmp = numpy.insert(tmp, 0, 0.)
        tmp = numpy.append(tmp, tmp[-1])

        choose =[]
        for isite in numpy.arange(numpy.shape(site_mod)[0]):
            site_d[isite] = 0.5*(tmp[0:len(tmp)-1]+tmp[1:len(tmp)])
            site_z[isite] = site_dem[isite] - site_d[isite]

            if ("smp" in qthresh[0].lower()):
                if (site_smp[isite] < qthresh[1]):
                    choose.append(isite)
                else:
                    continue

            if ("rms" in qthresh[0].lower()):
                if (site_rms[isite] > qthresh[1]):
                    choose.append(isite)
                else:
                    continue
        if out:
            print("\n",infile+":")
            print("original sites:",numpy.shape(site_mod)[0] )
            print("chosen sites:  ",len(choose) )

        all_sites = all_sites + numpy.shape(site_mod)[0]
        chs_sites = chs_sites + len(choose)
        if k == 1:
            # merged_log = site_log

            merged_mod = site_mod[choose]
            merged_sns = site_sns[choose]
            merged_rms = site_rms[choose]
            merged_smp = site_smp[choose]
            merged_con = site_con[choose]
            merged_x = site_x[choose].reshape(-1,1)
            merged_y = site_y[choose].reshape(-1,1)
            merged_z = site_z[choose]
            merged_d = site_d[choose]
            merged_gps = site_gps[choose].reshape(-1,1)
            merged_alt = site_alt[choose].reshape(-1,1)
            merged_dem = site_dem[choose].reshape(-1,1)
            merged_cov = site_cov[choose].T
            print(numpy.shape(merged_cov), numpy.shape(site_cov), k)
        else:
            # merged_log = numpy.vstack((merged_log, site_log))
            merged_mod = numpy.vstack((merged_mod, site_mod[choose]))
            merged_sns = numpy.vstack((merged_sns, site_sns[choose]))
            merged_rms = numpy.vstack((merged_rms, site_rms[choose]))
            merged_smp = numpy.vstack((merged_smp, site_smp[choose]))
            merged_con = numpy.vstack((merged_con, site_con[choose]))
            merged_x = numpy.vstack((merged_x, site_x[choose].reshape(-1,1)))
            merged_y = numpy.vstack((merged_y, site_y[choose].reshape(-1,1)))
            merged_z = numpy.vstack((merged_z, site_z[choose]))
            merged_d = numpy.vstack((merged_d, site_d[choose]))
            merged_gps = numpy.vstack((merged_gps, site_gps[choose].reshape(-1,1)))
            merged_alt = numpy.vstack((merged_alt, site_alt[choose].reshape(-1,1)))
            merged_dem = numpy.vstack((merged_dem, site_dem[choose].reshape(-1,1)))
            # print(numpy.shape(merged_cov), numpy.shape(site_cov), k)
            merged_cov = numpy.hstack((merged_cov, site_cov[choose].T))




    header = "merged model set:"+"".join("Date " + datetime.now().strftime(dateform))
    merged_models = {
            "header": header,
            "mod_ref": mod_ref,
            "mod_act": mod_act,
            "cov": merged_cov,
            "mod": merged_mod,
            "sns": merged_sns,
            "rms": merged_rms,
            "smp": merged_smp,
            "con": merged_con,
            "x": merged_x,
            "y": merged_y,
            "z": merged_z,
            "d": merged_d,
            "gps": merged_gps,
            "alt": merged_alt,
            "dem": merged_dem,
                    }

    print(" \n\n From ",all_sites, ", ",chs_sites, "were merged!")

    if outfile_name is not None:
            numpy.savez_compressed(file=outfile_name, **merged_models)

    return merged_models
