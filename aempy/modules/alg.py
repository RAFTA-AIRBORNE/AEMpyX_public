# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 16:13:37 2022

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
import scipy.sparse

from numba import njit

import numpy
import numpy.random
import functools

import aesys
import util
import inverse
import alg


warnings.simplefilter(action="ignore", category=FutureWarning)


rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")


def run_tikh_flightline(data_dir= None,
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


    """
    Setup data-related parameter dict
    """
    [nsite,ndata] = numpy.shape(dat_obs)
    sites = numpy.arange(nsite)

    data_act = ctrl["data"][1]
    data_err_add = ctrl["data"][2]
    data_err_mult = ctrl["data"][3]
    direction = ctrl["data"][4]

    dat_act = numpy.tile(data_act,(nsite,1))
    dat_err = numpy.zeros_like(dat_obs)
    for ii in sites:
        dat_err[ii, :], _ = inverse.set_errors(dat_obs[ii, :],
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

    mod_act = ctrl["model"][1].copy()
    mod_apr = ctrl["model"][2].copy()
    mod_var = ctrl["model"][3].copy()
    mod_bnd = ctrl["model"][4].copy()

    # print(mod_apr.shape)
    site_prior = numpy.zeros((nsite,numpy.shape(mod_apr)[0]))

    if "read" in setprior:
        if name.lower() not in prior_file.lower():
            error("Halfspace file name does not match! Exit.")
        site_prior = inverse.load_prior(prior_file,
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
    sequence = range(nsite)
    if "rev" in direction.lower():
        sites = sequence[::-1]
    else:
        sites = sequence

    # logsize = (2 + 7*maxiter)
    # site_log = numpy.full((len(sites),logsize), numpy.nan)
    mtmp = numpy.array([])
    # for ii in sites:
    for ii in [0, 1, 2]:
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

            if ii == 0:
                mod_ini = site_prior[ii,:]
                mod_apr = mod_ini.copy()
            else:
                mod_ini = inverse.insert_mod(M=mod_apr, m=mtmp[0], m_act=mod_act)
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

        site_dict = \
            alg.run_tikh_opt(Ctrl=ctrl, Model=model_dict, Data=data_dict,
                                      OutInfo=out)

#         Now store inversion results for this site:

        if out:
            print("site_dict: ",site_dict.keys())




        mtmp = site_dict["model"]
        dtmp = site_dict["data"]
        ctmp = site_dict["log"]

        # print(mtmp[0].shape)
        if ii==0:
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
#         dat_err[ii, :], _ = inverse.set_errors(dat_obs[ii, :],
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
#         site_prior = inverse.load_prior(prior_file,
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
#                 alg.run_tikh_opt(Ctrl=ctrl, Model=model_dict, Data=data_dict,
#                                   OutInfo=out)

#         if "occ" in runtype:
#             ens_dict =\
#                 alg.run_tikh_occ(Ctrl=ctrl, Model=model_dict, Data=data_dict,
#                                   OutInfo=out)

#         if "map" in runtype:
#             ens_dict =\
#                 alg.run_map(Ctrl=ctrl, Model=model_dict, Data=data_dict,
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
#         inverse.calc_stat_ens(ensemble=ens_modl, quantiles=percentiles, sum_stats=True)
#     stat_modl = numpy.vstack((m_quants, m_mean, m_stdv, m_skew, m_kurt, m_mode))

#     d_quants, d_mean, d_stdv, d_skew, d_kurt, d_mode = \
#         inverse.calc_stat_ens(ensemble=ens_modl, quantiles=percentiles, sum_stats=True)
#     stat_dcal = numpy.vstack((d_quants, d_mean, d_stdv, d_skew, d_kurt, d_mode))

#     r_quants, r_mean, r_stdv, r_skew, r_kurt, r_mode = \
#         inverse.calc_stat_ens(ensemble=ens_nrms, quantiles=percentiles, sum_stats=True)
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
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)


    obs = inverse.extract_dat(d_obs, d_act)
    err = inverse.extract_dat(d_err, d_act)

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
    m_apr, _ = inverse.transform_parameter(
        m_vec=m_apr, m_trn=m_trn, m_state=m_state, mode="f")
    m_ini, _ = inverse.transform_parameter(
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
        d_cal, dcal_state = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                  m_vec=model,
                                                  m_trn=m_trn, m_state=m_state,
                                                  d_trn=d_trn)

        nrmse_iter, smape_iter = inverse.calc_datafit(data_obs=d_obs,
                                                      data_cal=d_cal,
                                                      data_err=d_err,
                                                      data_act=d_act)
        print(niter,"####",model, m_trn, m_state,d_trn)
        # print("dfit0", nrmse_iter, smape_iter)
        if niter == 0:
            conv_status = 1
            print("ThreshVal =", thresh)
            if "rms" in thresh[3]:
                dfit_iter = nrmse_iter
                dfit_old = nrmse_iter
                dfit_0 = nrmse_iter


            if "smp" in thresh[3]:
                dfit_iter = smape_iter
                dfit_old = smape_iter
                dfit_0 = smape_iter

            model_old = model.copy()
            dnorm_iter = numpy.array([inverse.calc_dnorm(data_obs=d_obs, data_cal=d_cal,
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

        Jac = inverse.calc_jac(fwdcall=fwdcall, alt=alt,
                               m_vec=model, m_act=m_act, m_trn=m_trn, m_state=m_state,
                               d_vec=d_obs, d_act=d_act, d_trn=d_trn,
                               delta=delta, scalejac=False, out=False)

        cal = inverse.extract_dat(d_cal, d_act)
        # print("cal:", cal)
        # print("obs:", obs)

        m_iter = inverse.extract_mod(model, m_act)
        m_apri = inverse.extract_mod(m_apr, m_act)

        diff_m = m_iter - m_apri

        Jd = Wd@Jac
        JJ = Jd.T@Jd

        sensi = inverse.calc_sensitivity(Jac=Jd)

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
                model_test = inverse.insert_mod(M=model_test, m=m_test[itest],
                                                m_act=m_act)

                cali, d_state = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model_test,
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
                Calc target function for tau optimization
                """

                if "ufc" in regfun.lower():
                    reg_choice[itest] = inverse.calc_ufc(
                        dnorm[itest], mnorm[itest])

                elif "upr" in regfun.lower():
                    M = Jd@G
                    reg_choice[itest] = inverse.calc_upr(
                        dnorm[itest], M, d_err)

                elif "gcv" in regfun.lower():
                    M = Jd@G
                    reg_choice[itest] = inverse.calc_gcv(dnorm[itest], M)

                elif "mle" in regfun.lower():
                    M = Jd@G
                    reg_choice[itest] = inverse.calc_gcv(dnorm[itest], M)

                elif "fix" in regfun.lower() or "lcc" in regfun.lower():
                    pass

                else:
                    error("Regularisation method "+regfun.lower() +
                          " not yet implemented! Exit.")

        if "fix" in regfun.lower():
            g_index = 0

        elif "lcc" in regfun.lower():
            g_index = inverse.calc_lc_corner(dnorm, mnorm)+gshift
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

        reg = [tau[g_index, 0].item(), tau[g_index, 1].item()]
        mdl = m_test[g_index, :]
        model = inverse.insert_mod(M=model, m=mdl, m_act=m_act)

        """
        Line Search
        """
        if do_linesearch:
            model, dfit_iter = inverse.run_linesearch(fwdcall, alt,
                                                      d_obs=d_obs, d_err=d_err, d_trn=d_trn, d_act=d_act, d_state=d_state,
                                                      model=model, m_delta=m_delta, m_act=m_act, m_trn=m_trn, m_state=m_state,
                                                      dfit=dfit_iter, mdfit=thresh[3],
                                                      facreduce=facreduce, maxreduce=maxreduce, out=OutInfo)
            d_cal, _ = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model,
                                                      m_trn=m_trn, m_state=m_state, d_trn=d_trn, d_act=d_act)

            nrmse_iter, smape_iter = inverse.calc_datafit(data_obs=d_obs,
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
            rvals_iter = numpy.append(rvals_iter, reg)
            dfits_iter = numpy.append(dfits_iter, [nrmse_iter, smape_iter])

    # if OutInfo==True:
    print(" %s initial model: dfit =  %7.3f"
          % (profname, dfit_0))
    print(" %s final model:  iter  %6i NRMSE =  %7.3f,  SMAPE = %4.1f percent, %s RegPars are %10.4g / %10.4g"
          % (profname, niter, nrmse_iter, smape_iter, regfun.lower(),  tau[g_index, 0],  tau[g_index, 1]))
    # print(numpy.shape(sensi))
    sens = sensi  # inverse.insert_mod(M=sens, m=sensi,m_act=m_act)
    modl, m_state = inverse.transform_parameter(
        m_vec=model, m_trn=m_trn, m_state=m_state,  mode="b")
    modl = inverse.extract_mod(modl, m_act)
    merr = m_err[g_index, :].flat
    # inverse.insert_mod(M=merr, m=tmp, m_act=m_act)
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
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)

    obs = inverse.extract_dat(d_obs, d_act)
    err = inverse.extract_dat(d_err, d_act)
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
    m_apr, _ = inverse.transform_parameter(
        m_vec=m_apr, m_trn=m_trn, m_state=m_state, mode="f")
    m_ini, _ = inverse.transform_parameter(
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
        d_cal, dcal_state = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                  m_vec=model,
                                                  m_trn=m_trn, m_state=m_state,
                                                  d_trn=d_trn)

        nrmse_iter, smape_iter = inverse.calc_datafit(data_obs=d_obs,
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
            dnorm_iter = numpy.array([inverse.calc_dnorm(data_obs=d_obs, data_cal=d_cal,
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

        Jac = inverse.calc_jac(fwdcall=fwdcall, alt=alt,
                               m_vec=model, m_act=m_act, m_trn=m_trn, m_state=m_state,
                               d_vec=d_obs, d_act=d_act, d_trn=d_trn,
                               delta=delta, scalejac=False, out=False)

        cal = inverse.extract_dat(d_cal, d_act)

        m_iter = inverse.extract_mod(model, m_act)
        m_apri = inverse.extract_mod(m_apr, m_act)

        diff_m = m_iter - m_apri

        Jd = Wd@Jac
        JJ = Jd.T@Jd
        # JJT = Jac.T@Cdi@Jac

        sensi = inverse.calc_sensitivity(Jac=Jd)

        # Cmi0 = Cm0.multiply(tau)

        if niter == 0 and set_taustart:
            tau_start = inverse.calc_regstart(D=JJ, M=Cm1)
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
        model_test = inverse.insert_mod(M=model_test, m=m_test,
                                        m_act=m_act)

        cali, d_state = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model_test,
                                              m_trn=m_trn, m_state=m_state, d_trn=d_trn, d_act=d_act)

        r_test = Wd@(obs - cali).T

        dnorm = scipy.linalg.norm(r_test)
        mnorm = scipy.linalg.norm(m_test)

    reg = [tau_iter]
    mdl = m_test
    model = inverse.insert_mod(M=model, m=mdl, m_act=m_act)

    """
    Line Search
    """
    if do_linesearch:
        model, dfit_iter = inverse.run_linesearch(fwdcall, alt,
                                                                          d_obs=d_obs, d_err=d_err, d_trn=d_trn, d_act=d_act, d_state=d_state,
                                                                          model=model, m_delta=m_delta, m_act=m_act, m_trn=m_trn, m_state=m_state,
                                                                          dfit=dfit_iter, mdfit=thresh[3],
                                                                          facreduce=facreduce, maxreduce=maxreduce, out=OutInfo)
        nrmse_iter, smape_iter = inverse.calc_datafit(data_obs=d_obs,
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
    sens = sensi  # inverse.insert_mod(M=sens, m=sensi,m_act=m_act)
    modl, m_state = inverse.transform_parameter(
        m_vec=model, m_trn=m_trn, m_state=m_state,  mode="b")
    modl = inverse.extract_mod(modl, m_act)
    merr = m_err.flat
    # inverse.insert_mod(M=merr, m=tmp, m_act=m_act)
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
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
                                                      e_vec=d_err,
                                                      d_trn=d_trn,
                                                      d_state=d_state)

    obs = inverse.extract_dat(d_obs, d_act)
    err = inverse.extract_dat(d_err, d_act)

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
    m_apr, _ = inverse.transform_parameter(
        m_vec=m_apr, m_trn=m_trn, m_state=m_state, mode="f")
    m_ini, _ = inverse.transform_parameter(
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
        d_cal, dcal_state = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt,
                                                  m_vec=model,
                                                  m_trn=m_trn, m_state=m_state,
                                                  d_trn=d_trn)

        nrmse_iter, smape_iter = inverse.calc_datafit(data_obs=d_obs,
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
            dnorm_iter = numpy.array([inverse.calc_dnorm(data_obs=d_obs, data_cal=d_cal,
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

        Jac = inverse.calc_jac(fwdcall=fwdcall, alt=alt,
                               m_vec=model, m_act=m_act, m_trn=m_trn, m_state=m_state,
                               d_vec=d_obs, d_act=d_act, d_trn=d_trn,
                               delta=delta, scalejac=False, out=False)

        cal = inverse.extract_dat(d_cal, d_act)

        m_iter = inverse.extract_mod(model, m_act)
        m_apri = inverse.extract_mod(m_apr, m_act)

        diff_m = m_iter - m_apri

        Jd = Sd@Jac
        sensi = inverse.calc_sensitivity(Jac=Jd)

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
            model_test = inverse.insert_mod(M=model_test, m=m_test[itest],
                                            m_act=m_act)
            cali, d_state = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt, m_vec=model_test,
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
                reg_choice[itest] = inverse.calc_gcv(dnorm[itest], M)

            elif "mle" in regfun.lower():
                M = Jd@C@Jd.T
                reg_choice[itest] = inverse.calc_gcv(dnorm[itest], M)

            elif "fix" in regfun.lower() or "lcc" in regfun.lower():
                pass

            else:
                error("Regularisation method "+regfun.lower() +
                      " not yet implemented! Exit.")

        if "fix" in regfun.lower():
            g_index = 0

        elif "lcc" in regfun.lower():
            g_index = inverse.calc_lc_corner(dnorm, mnorm)+gshift
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
        model = inverse.insert_mod(M=model, m=mdl, m_act=m_act)

        """
        Line Search
        """
        if do_linesearch:
            model, dfit_iter = inverse.run_linesearch(fwdcall, alt,
                                                                              d_obs=d_obs, d_err=d_err, d_trn=d_trn, d_act=d_act, d_state=d_state,
                                                                              model=model, m_delta=m_delta, m_act=m_act, m_trn=m_trn, m_state=m_state,
                                                                              dfit=dfit_iter, mdfit=thresh[3],
                                                                              facreduce=facreduce, maxreduce=maxreduce, out=OutInfo)
            nrmse_iter, smape_iter = inverse.calc_datafit(data_obs=d_obs,
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
    sens = sensi  # inverse.insert_mod(M=sens, m=sensi,m_act=m_act)
    modl, m_state = inverse.transform_parameter(
        m_vec=model, m_trn=m_trn, m_state=m_state,  mode="b")
    modl = inverse.extract_mod(modl, m_act)
    merr = m_err[g_index, :].flat
    # inverse.insert_mod(M=merr, m=tmp, m_act=m_act)
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

    m_ref = inverse.insert_mod(M=m_bas, m_act=m_act, m=m_opt)
    c_ref = results["cpost"]

    """
    Draw perturbed data set: d  ̃ ∼ N (d, Cd)
    """
    d_ens = inverse.generate_data_ensemble(dref=d_obs, dact=d_act,
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
    m_ens = inverse.generate_param_ensemble(mref=m_ref, mact=m_act,
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
        Model["m_apr"] = inverse.insert_mod(M=m_ref.copy(), m_act=m_act,
                                            m=m_ens[isample, :])

        if "opt" in invtype.lower():

            results =\
                run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                             OutInfo=OutInfo)

        m = inverse.insert_mod(M=m_ref, m_act=m_act, m=results["model"][0])
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
def run_EnK(Ctrl=None, Model=None, Data=None, OutInfo=True):
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

    M. Iglesias, D. M. McGrath, M. V. Tretyakov, and Susan T Francis,
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
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
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
    m_ref, m_state = inverse.transform_parameter(m_vec=m_ref,
                                                 m_trn=m_trn, m_state=m_state, mode="f")

    # print(m_state, " mref", m_ref[:36])
    """
    Draw prior model: m_p∼ N (0, 1 (LT L)−1 )
    """
    m_ens = inverse.generate_param_ensemble(mref=m_ref, mact=m_act,
                                            nens=nsamples,
                                            perturb=["gauss", Cm1],
                                            out=OutInfo)
    # """
    # Draw perturbed data set: d  ̃ ∼ N (d, Cd)
    # """
    # d_ens = inverse.generate_data_ensemble(dref=d_obs, dact = d_act,
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
            m = inverse.insert_mod(model, m_smp, m_act)
            d_cal[isample, :], _ = inverse.calc_fwdmodel(fwdcall=fwdcall, alt=alt,
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

        Gyy = inverse.calc_encovar(x=d_cal, y=d_cal, method=0, out=True)
        Gxy = inverse.calc_encovar(x=m_ens, y=d_cal, method=0, out=True)

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
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
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
    m_ref = inverse.insert_mod(M=m_bas, m_act=m_act, m=m_opt)
    c_ref = results["cpost"]

    """
    now calculate the SVD of the Jacobian
    """
    Jacd = results["jacd"]

    if randsvd:
        U, S, Vt = inverse.rsvd(
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
    chek  how mauch of Jacd is explained by k
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
    m_ens = inverse.generate_model_ensemble(mref=m_ref, mact=m_act,
                                            nens=nsamples,
                                            perturb=["gauss", c_ref,
                                                     numpy.array([])],
                                            out=OutInfo)

    m_prj = numpy.zeros_like(m_ens)
    for isample in numpy.arange(nsamples):
        m_prj[isample, :] = inverse.project_nullspace(
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
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
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
    m_ref, m_state = inverse.transform_parameter(
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
        U, S, Vt = inverse.rsvd(
            Jacd.T, rank=k, n_oversamples=0, n_subspace_iters=2)
    else:
        U, S, Vt = scipy.linalg.svd(Jacd.T, full_matrices=False)

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
