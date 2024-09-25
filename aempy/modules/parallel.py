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

# from numba import njit

import numpy
import numpy.random
import functools

import aesys
import util
import inverse
import core1d
import alg

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
            numpy.array([DataTrans, data_active, DatErr_add, DatErr_mult, ReverseDir], dtype=object),
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
    
    
    site_prior = numpy.zeros((nsite,numpy.shape(mod_apr)[0]))
    
    if "read" in setprior:
        if name.lower() not in prior_file.lower():
            error("Halfspace file name does not match! Exit.")
        site_prior = inverse.load_prior(prior_file,
                                       m_ref=mod_apr,
                                       m_apr=mod_apr,
                                       m_act=mod_act)
    if "set" in setprior:      
        for ii in sites:
                site_prior[ii, :] = mod_apr


    uncert = ctrl["uncert"][0]
    header = ctrl["header"][0]

# This is the main loop over sites in a flight line or within an area:

    """
    Loop over sites
    """

    # logsize = (2 + 7*maxiter)
    # site_log = numpy.full((len(sites),logsize), numpy.nan)
    mtmp = numpy.array([])
    for ii in sites:
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
                mod_ini = mtmp[0].copy()
                mod_apr = mod_ini.copy()

        elif "set" in setprior:
            mod_apr = site_prior[ii,:]
            mod_ini = mod_apr.copy()


        # print("\n",ii) 
        # print(mod_ini[ii, :])
        # print(mod_apr)
        # print(mod_var[ii, :])
        
        model_dict = dict([
            ("m_act", mod_act),
            ("m_apr", mod_apr),
            ("m_var", mod_var),
            ("m_bnd", mod_bnd),
            ("m_ini", mod_ini)
            ])



        site_dict = \
            alg.run_tikh_opt(Ctrl=ctrl, Model=model_dict, Data=data_dict,
                                      OutInfo=out)

#         Now store inversion results for this site:

        if out:
            print("site_dict: ",site_dict.keys())

        
        mtmp = site_dict["model"]
        dtmp = site_dict["data"]
        ctmp = site_dict["log"]
        
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
#             numpy.array([DataTrans, data_active, DatErr_add, DatErr_mult, ReverseDir], dtype=object),
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

