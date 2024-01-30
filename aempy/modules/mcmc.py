#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 18:42:34 2023

@author: vrath
"""
# import pymcmcstat 
# from pymcmcstat.MCMC import MCMC
# from pymcmcstat.ParallelMCMC import ParallelMCMC
import matplotlib.pyplot 
# import pymcmcstat
# import mcmcplot


# @ray.remote


def run_EMCEE(Ctrl=None, Model=None, Data=None, OutInfo=True):
    """
    

    Returns
    -------
    nix : TYPE
        DESCRIPTION.
        
    References:
        
        Daniel Foreman-Mackey and David W. Hogg and Dustin Lang and Jonathan Goodman (2013)
        "emcee": The {MCMC} Hammer
        Publications of the Astronomical Society of the Pacific, doi:10.1086/670067

    """
    import emcee

    nix = [0]
    
    return nix

def EMCEE_modelfun(xdata, parameter):
    """

    Parameters
    ----------
    xdata : TYPE
        DESCRIPTION.
    parameter : TYPE
        DESCRIPTION.

    Returns
    -------
    y : TYPE
        DESCRIPTION.

    """
    
    m = inverse.calc_fwdmodel(fwdcall=None,
              alt=None,
              m_vec = numpy.array([]),
              m_trn=0,
              m_state = 0,
              m_params=numpy.array([]),
              d_act=numpy.array([]),
              d_trn=0,
              d_state=0,
              OutInfo=True):
    m = parameter[0]
    b = parameter[1]
    nrow, ncol = xdata.shape
    y = np.zeros([nrow,1])
    y[:, 0] = m*xdata.reshape(nrow,) + b
    
    return y

def EMCEE_ssfun(parameter, data):
    """

    Parameters
    ----------
    parameter : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    xdata = data.xdata[0]
    ydata = data.ydata[0]
    # eval model
    ymodel = aempy_modelfun(xdata, parameter)
    # calc sos
    res = ymodel[:, 0] - ydata[:, 0]
    return (res**2).sum(axis = 0)


def run_dram(Ctrl=None, Model=None, Data=None, OutInfo=True):
    """
    

    Returns
    -------
    nix : TYPE
        DESCRIPTION.
        
    References:
        
        Haario H, Laine M, Mira A, Saksman E (2006)
        DRAM: Efficient adaptive MCMC 
        Statistics and Computing , 16, 339-354
        
        Miles P (2019)
        pymcmcstat: A Python Package for Bayesian Inference Using 
        Delayed Rejection Adaptive Metropolis 
        Journal of Open Source Software, 4(38), 1417
        https://doi.org/10.21105/joss.01417
        
        Miles P (2020)
        prmiles/pymcmcstat: v1.9.1 (v1.9.1)
        Zenodo, doi:10.5281/zenodo.4127047
        
       

    # define test model function
    def modelfun(xdata, parameter):
        m = parameter[0]
        b = parameter[1]
        nrow, ncol = xdata.shape
        y = np.zeros([nrow,1])
        y[:, 0] = m*xdata.reshape(nrow,) + b
        return y
    
    def ssfun(parameter, data):
        xdata = data.xdata[0]
        ydata = data.ydata[0]
        # eval model
        ymodel = modelfun(xdata, parameter)
        # calc sos
        res = ymodel[:, 0] - ydata[:, 0]
        return (res**2).sum(axis = 0)

 
        

    """
    import pymcmcstat 
    from pymcmcstat.MCMC import MCMC
    from pymcmcstat.ParallelMCMC import ParallelMCMC 
    import mcmcplot
    
    nix = [0]
    
    return nix

def dram_modelfun(xdata, parameter, **kwargs):
    """

    Parameters
    ----------
    xdata : TYPE
        DESCRIPTION.
    parameter : TYPE
        DESCRIPTION.

    Returns
    -------
    y : TYPE
        DESCRIPTION.

    """
    
    
          fwdcal = 
          alt=None,
          m_vec = numpy.array([]),
          m_trn=0,
          m_state = 0,
          m_params=numpy.array([]),
          d_act=numpy.array([]),
          d_trn=0,
          d_state=0,
          OutInfo=True)
    
    
    m = inverse.calc_fwdmodel(fwdcall, alt,
              m_vec, m_trn, m_state, m_params, 
              d_act, d_trn, d_state, OutInfo):
    
    return y

def dram_ssfun(parameter, data, **kwargs):
    """

    Parameters
    ----------
    parameter : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    xdata = data.xdata[0]
    ydata = data.ydata[0]
    # eval model
    ymodel = aempy_modelfun(xdata, parameter)
    # calc sos
    res = ymodel[:, 0] - ydata[:, 0]
    return (res**2).sum(axis = 0)


# def run_hmc(Ctrl=None, Model=None, Data=None, OutInfo=True):
#     """
    

#     Returns
#     -------
#     nix : TYPE
#         DESCRIPTION.
        
    
#         Fichtner A, Zunino A (2019)
#         Hamiltonian Nullspace Shuttles 
#         Geophysical Research Letters, 46
#         doi:10.1029/2018GL080931
    
#         Zunino A, Gebraad L, Ghirotto A,  Fichtner, A (2023)
#         HMCLab: a framework for solving diverse geophysical inverse problems 
#         using the Hamiltonian Monte Carlo method 
#         arXiV, 2023
    
#         Fichtner A (2021)
#         Lecture Notes on Inverse Theory 
#         Cambridge Open Engage Platform

#     """
    
#     import hmclab
#     from hmclab import HMC
    
    
#     target = aempy_hmc(temperature=20)
    
    
#     hmc().sample(
#         "bin_samples/tutorial_3_styblinski-tang.h5",
#         target,
#         proposals=20000,
#         stepsize=1.0,
#         overwrite_existing_file=True,
#         )

#     nix = [0]
        
#     return nix


# class hmc(hmclab.Distributions._AbstractDistribution):
    
#     def __init__(self, dimensions=None, temperature=20.0):

#         # Default to 2 dimensions if not specified
#         if dimensions is None:
#             self.dimensions = 2
#         else:
#             self.dimensions = dimensions

#         # Default to temperature = 10 if not specified
#         self.temperature = temperature

#     def misfit(self, m):

#         # Check the shape of the passed model
#         assert m.shape == (self.dimensions, 1)
#         calc_fwdmodel(fwdcall=None,
#                   alt=None,
#                   m_vec = numpy.array([]),
#                   m_trn=0,
#                   m_state = 0,
#                   m_params=numpy.array([]),
#                   d_act=numpy.array([]),
#                   d_trn=0,
#                   d_state=0,
#                   OutInfo=True):
#         return mf / self.temperature

#     def gradient(self, m):

#         assert m.shape == (self.dimensions, 1)

#         return 0.5 * (4 * (m**3) - 32 * m + 5) / self.temperature

#     def generate(self, m):
#         raise NotImplementedError()





# def run_dream(Ctrl=None, Model=None, Data=None, OutInfo=True):
#     """
    

#     Returns
#     -------
#     nix : TYPE
#         DESCRIPTION.
        
#     References:
    
#         Shockley, Erin M. and Vrugt, Jasper A. and Lopez, Carlos F. (2018)
#         PyDREAM: high-dimensional parameter inference for biological models in Python
#         Bioinformatics, 34, 695--697, doi:10.1093/bioinformatics/btx626

#         E Laloy and JA Vrugt (2012)
#         Prior exploration of hydrologic models using 
#         multiple-try DREAM(ZS) and high-performance computing
#         Water Resour. Res., 48, W01526
#         doi10.1029/2011WR010608
#     """
#     import pydream
    
#     return nix

# def dream_modelfun(xdata, parameter):
#     """

#     Parameters
#     ----------
#     xdata : TYPE
#         DESCRIPTION.
#     parameter : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     y : TYPE
#         DESCRIPTION.

#     """
    
#     m = inverse.calc_fwdmodel(fwdcall=None,
#               alt=None,
#               m_vec = numpy.array([]),
#               m_trn=0,
#               m_state = 0,
#               m_params=numpy.array([]),
#               d_act=numpy.array([]),
#               d_trn=0,
#               d_state=0,
#               OutInfo=True):
#     m = parameter[0]
#     b = parameter[1]
#     nrow, ncol = xdata.shape
#     y = np.zeros([nrow,1])
#     y[:, 0] = m*xdata.reshape(nrow,) + b
    
#     return y

# def dream_ssfun(parameter, data):
#     """

#     Parameters
#     ----------
#     parameter : TYPE
#         DESCRIPTION.
#     data : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     xdata = data.xdata[0]
#     ydata = data.ydata[0]
#     # eval model
#     ymodel = aempy_modelfun(xdata, parameter)
#     # calc sos
#     res = ymodel[:, 0] - ydata[:, 0]
#     return (res**2).sum(axis = 0)


def run_mhsimple(Ctrl=None, Model=None, Data=None, OutInfo=True):
    """
    run MCMC chains - straightforward Metropolis-Hastings.

    Returns
    -------
    nix : TYPE
        DESCRIPTION.
        
    References:

   
    N Metropolis, AW Rosenbluth, MN Rosenbluth, AH Teller, and E Teller (1953)
    Equation of state calculations by fast computing machines 
    J. Chem. Phys., 21, 1953, doi:10.1063/1.1699114.
  
    WK Hastings (1970)
    Monte Carlo sampling methods using Markov chains and their applications
    Biometrika, 57, 97–109    
    
    
   Adapted from:
   Kiyan, D. & Rath, V.
   Inverse Methods for Airborne Electromagnetic Data from the Tellus Surveys.
   The aempy Toolbox (2017)  
   Dublin Institute for Advanced Studies (DIAS) &
   Geological survey of Ireland, Dublin

   Kiyan, D.; Rath, V.; Muller, M. R.; Ture, M. D. & Hodgson, J. (2022)
   1-D inversion of frequency-domain airborne electromagnetic data
   using the open-source aempy toolbox
   Journal of Applied Geophysics, 198, 104562


        Ctrl_mh = dict([
           ("system", [AEM_system, FwdCall]),
           ("mh",[nsample,nburnin,
                   accp_check, accp_max, accp_min, stepini, stepfac]),
           ("covar", [Cm, Cd]),
           ("transform", [DataTrans, ParaTrans])
           ])

   vrath  Feb 10, 2023

   @author: vrath
       """
    system, fwdcall = Ctrl["system"]
    nsample, nburnin, accp_check, accp_max, accp_min, stepini, stepfac = Ctrl["mh"]
    cm, _ = Ctrl["covar"]
    d_trans, m_trans = Ctrl["transform"]

    """
    unpack data block
    Data = [data_act[ii,:], data_obs[ii,:], data_error[ii,:], site_alt[ii]]
    """
    d_act = Data["d_act"] #.reshape(-1,1)
    d_obs = Data["d_obs"] #.reshape(-1,1)
    d_err = Data["d_err"]# .reshape(-1,1)
    alt   = Data["alt"]

    d_cal = numpy.nan * numpy.ones_like(d_obs)
    
    d_state = 0
    d_obs, d_err, dobs_state = inverse.transform_data(d_vec=d_obs,
                                              e_vec = d_err,
                                              d_trn=d_trn,
                                              d_state = d_state)
    """
    unpack model block
    
    Model =\
    [model_act, model_prior, model_var, model_bounds, model_ini]
    """
    
    m_act = Model["m_act"]
    m_ref = Model["m_apr"]
    
    m_state = 0
    m_ref, m_state = inverse.transform_parameter(m_vec=m_ref, m_trn=m_trn, m_state=m_state, mode="f")


       scount = 0
       step = stepini
       if onlyactive==False:
           mpara = 7*nlyr-1
           mchain = np.zeros((nsample+1,mpara+ndata+1))
           mchain[scount,0:mpara]= m_old
       else:
           mchain = np.zeros((nsample+1,mpara+ndata+1))
           mchain[scount,0:mpara]= m_old[mactive==1]

       mchain[scount,mpara:mpara+ndata]= dc_old
       mchain[scount,mpara+ndata]= ss_old


       reject = 0
       accept = 0
       for n in range(nsample):

       #    print('\n',n)
           #draw random numbers from the prior model (prior_avg)
           m_g = np.random.default_rng().normal(sizepar[0])
       #    print(' m_g   ',m_g[0*nlyr:1*nlyr])
           m_c = np.matmul(LC,m_g)
           m_s = step*mactive*m_c
       #    print(' m_s   ', m_s[0*nlyr:1*nlyr])
       #    print(' m_c   ', m_c[1*nlyr:2*nlyr])
       #    print(' m_s   ', m_s[1*nlyr:2*nlyr])
           #m_sample = m_old + tep*prior_std*mactive*np.random.default_rng().normal(sizepar[0])
           m_sample = m_old + m_s
       #    print(' m_smp   ', m_sample[1*nlyr:2*nlyr])
           if np.all(m_sample < m_upper) and np.all(m_sample > m_lower):
               scount = scount+1
               m_new = m_sample
               fp_new=np.sum(np.power((m_new-m_old)/prior_std,2));
               #data_calc = core1d.aemfwd1d_aem05(mode,alt,nlyr,m_new)
               data_calc = eval('core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,m_new)')
               
               # data_calc = calc_fwdmodel(fwdcall=fwdcall, 
               #                           alt=alt,
               #    m_vec = m_new,
               #    m_trn=0,
               #    m_state = 0,
               #    m_params=numpy.array([]),
               #    d_act=numpy.array([]),
               #    d_trn=0,
               #    d_state=0,
               #    OutInfo=True):
               
               
               
               
               
               
               ss_new = np.sum(np.power(abs(data_obs - data_calc)/data_std,2))               
               dc_new = data_calc
           #fl_new=np.exp(-0.5*ss_new)
       #        a1 = -0.5*(ss_new-ss_old)
       #        a2 = -0.5*(fp_new-fp_old)
       #        a3 = a1+a2
       #        print(ss_new,ss_old)
               alpha = min(1,np.exp(-0.5*(ss_new-ss_old) -0.5*(fp_new-fp_old)));
           #print (alpha)
               p = np.random.rand(1,1)
               if p < alpha:
               # accept new candidate

                 accept = accept+1

                   if onlyactive:
                         mchain[scount,0:mpara]= m_new[mactive==1]
                   else:
                         mchain[scount,0:mpara]= m_new

                   mchain[scount,mpara:mpara+ndata]= dc_new
                   mchain[scount,mpara+ndata]= ss_new

                   m_old  = m_new
                   ss_old = ss_new;
                   fp_old = fp_new;
                   dc_old = dc_new
               else:
                   if onlyactive:
                         mchain[scount,0:mpara]= m_old[mactive==1]
                   else:
                         mchain[scount,0:mpara]= m_old

                   mchain[scount,mpara:mpara+ndata]= dc_old
                   mchain[scount,mpara+ndata]= ss_old
           else:
                   reject=reject+1

           if np.mod(n,accp_check)==0 and not n==0:
               accpp=100*accept/scount
               if accpp > accp_max: step = step/stepfac
               if accpp < accp_min: step = step*stepfac
               print (' percentage of samples accepted = %6.2f of %6i, with  %6i out of bounds, step set to %6.4f'\
               % (accpp, scount, reject, step))

       print (' final percentage of samples accepted = %6.2f ' % (100*accept/scount))
       #np.savez('Mh_chain',mchain)
       print(np.shape(mchain[nburnin:scount,1]))
       # np.savez('Mh_chain1',mchain=mchain)
       #np.savez_compressed('Mh_chain2',mchain=mchain[nburnin:scount,:])




       mod_avg = numpy.mean(ensemble, axis=1)
       mod_std = numpy.std(ensemble, axis=1)

       mod_prc = numpy.percentile(ensemble, Percentiles)

       mod_med, mod_mad = inverse.calc_made(ensemble)

       mhresults =\
           dict([
           ("avg", mod_avg),
           ("std", mod_std),
           ("med", mod_med),
           ("mad", mod_mad),
           ("percentiles", mod_prc),
           ])

       if ens_out:
           mhresults["ens"] = ensemble


       return  mhresults




def aempy_modelfun(xdata, parameter):
    """

    Parameters
    ----------
    xdata : TYPE
        DESCRIPTION.
    parameter : TYPE
        DESCRIPTION.

    Returns
    -------
    y : TYPE
        DESCRIPTION.

    """
    
    m = inverse.calc_fwdmodel(fwdcall=None,
              alt=None,
              m_vec = numpy.array([]),
              m_trn=0,
              m_state = 0,
              m_params=numpy.array([]),
              d_act=numpy.array([]),
              d_trn=0,
              d_state=0,
              OutInfo=True)
        
    m = parameter[0]
    b = parameter[1]
    nrow, ncol = xdata.shape
    y = np.zeros([nrow,1])
    y[:, 0] = m*xdata.reshape(nrow,) + b
    
    return y

def aempy_ssfun(parameter, data):
    """

    Parameters
    ----------
    parameter : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    xdata = data.xdata[0]
    ydata = data.ydata[0]
    # eval model
    ymodel = aempy_modelfun(xdata, parameter)
    # calc sos
    res = ymodel[:, 0] - ydata[:, 0]
    return (res**2).sum(axis = 0)
