#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
# ---


"""
Plot diverse uncertainty parameters

Created May 2023

@author: vrath
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings
import random
import functools
from cycler import cycler

import numpy
import scipy.interpolate
import scipy.linalg



import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.axis



AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import util
import aesys
import inverse
import mcmc


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")
print(pymcmcstat.__version__)

now = datetime.now()

"""
System related settings.
Data transformation is now allowed with three possible options:
DataTrans   = 0           raw data
            = 1           natural log of data
            = 2           asinh transformation
An error model is applied for the raw data, which is
mixed additive/multiplicative. in case of data transformation,
errors are also transformed.
"""
AEM_system = "genesis"
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, Pars, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.03
    data_active = numpy.ones(NN[2], dtype="int8")
    CompDict = Pars[3]
    CompLabl = list(CompDict.keys())
    print(CompLabl)
    # Pars[0] = numpy.round(Pars[0],1)/1000.

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, Pars, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + "good" hoizontals"
    CompDict = Pars[3]
    CompLabl = list(CompDict.keys())


"""
Input format is "npz"
"""
FileList = "search"  # "search", "read"# 
# SearchStrng = "*PLM3s_k3.npz"
# SearchStrng = "*1379*k[1,2,3,5].npz"
SearchStrng = ""


FileList = "set"  # "search", "read"

InDatDir =  AEMPYX_DATA + "/Projects/InvParTest/proc_delete_PLM3s/"
if not InDatDir.endswith("/"): InDatDir=InDatDir+"/"

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    # dat_files = []
    dat_files = [InDatDir+"A1_rect_StGormans_FL11379-0_proc_delete_PLM3s_k3.npz"]

    
    dat_files = [os.path.basename(f) for f in dat_files]  
else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)


ns = numpy.size(dat_files)
if ns ==0:
    error("No modfiles set!. Exit.")

print("Filelist:")
print(dat_files)

"""
Output format is ".npz"
"""

ResultsDir =  InDatDir + "results_dram"
if not ResultsDir.endswith("/"): ResultsDir=ResultsDir+"/"
print("Models written to dir: %s " % ResultsDir)


if not os.path.isdir(ResultsDir):
    print("File: %s does not exist, but will be created" % ResultsDir)
    os.mkdir(ResultsDir)


"""
Model definition
"""

SetPrior = "set"
ParaTrans = 1

Nlyr = 30
dzstart = 2.5
dzend = 10.
dz = numpy.logspace(numpy.log10(dzstart), numpy.log10(dzend), Nlyr)
z = numpy.append(0.0, numpy.cumsum(dz))


mod_act, mod_apr, mod_var, mod_bnd, m_state = inverse.init_1dmod(Nlyr)

mod_act[0*Nlyr:1*Nlyr] = 1
sizepar = numpy.shape(mod_act)
mpara = sizepar[0]

Guess_r = 100.0  # initial guess for resistivity in mod_apr
Guess_s = 0.3   # mod_std defines standard deviation of mod_apr
mod_apr[0*Nlyr:1*Nlyr] = Guess_r
mod_var[0*Nlyr:1*Nlyr] = numpy.power(Guess_s,2)
mod_apr[6*Nlyr:7*Nlyr-1] = dz[0:Nlyr - 1]
mod_var[6*Nlyr:7*Nlyr-1] = numpy.power(1.,2)


# mod_bnd = mumpy.array([])
max_val = 1.e+30
min_val = 1.e-30
# max_val = mod_apr[mod_act!=0] + 3*mod_std[mod_act!=0]
# mod_bnd[mod_act!=0, 1] = max_val
# min_val = mod_apr[mod_act!=0] - 3*mod_std[mod_act!=0]
# mod_bnd[mod_act!=0, 0] = min_val
mod_bnd[:,0] = min_val
mod_bnd[:,1] = max_val


if OutInfo:
    #   print \
    #   (" Parameter set for inverting: \n", mod_act)
    print(" Layer thicknesses: \n", dz)
    print(" Layer interface depths: \n", z)
    print(" Initial halfspace resistivity of %6.2f Ohm*m" % (Guess_r))
    print(" Log Standard error of %6.2f " % (Guess_s))
    if not (mod_bnd == None) or (numpy.size(mod_bnd) == 0):
        print(" Upper limits: \n", mod_bnd[:, 1])
        print(" Lower limits: \n", mod_bnd[:, 0])






"""
Define inversion type  optional additional parameters (e.g., Waveforms )
"""

RunType = "DRAM"  # "DREAM" "EMCEE"


# Sample = "random"
# Sample = "distance list"
Sample = "distance list"
if "rand" in Sample:
    NSamples = 1

elif "list" in Sample:
    if "pos" in Sample:
        Samplist = [100, 200]
    if "dis" in Sample:
        Distlist = [ 1500.]



"""
Define inversion type  optional additional parameters (e.g., Waveforms )
"""

RunType = "DRAM"  # "DREAM" "EMCEE"





fcount =0
for file in dat_files:

    start = time.time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    fileout = OutDatDir + name + outstrng

    numpy.savez_compressed(file=fileout+"_ctrl"+OutFileFmt,**Ctrl)

    print("\n Reading file " + filein)
    DataObs, Header, _ = aesys.read_aempy(File=filein,
                                   System=AEM_system, OutInfo=False)

    fl_name = DataObs[0, 0]
    site_x = DataObs[:, 1]
    site_y = DataObs[:, 2]
    site_gps = DataObs[:, 3]
    site_alt = DataObs[:, 4]
    site_dem = DataObs[:, 5]
    dat_obs =  DataObs[:, 6:6+NN[2]]
    [nsite,ndata] = numpy.shape(dat_obs)
    dat_act = numpy.tile(data_active,(nsite,1))
     
    """
    construct site_list
    """
    site_x = site_x - site_x[0]
    site_y = site_y - site_y[0]
    site_r = numpy.sqrt(numpy.power(site_x, 2.0) + numpy.power(site_y, 2.0))
    
    site_list = []
    if "rand" in Sample:
        site_list = random.sample(range(len(site_x)), NSamples)
 
    elif "list" in Sample:
        if "posi" in Sample:
            site_list = Samplist
        if "dist" in Sample:
            for nid in numpy.arange(len(Distlist)):
                nds = (numpy.abs(Distlist[nid] - site_r)).argmin()
                site_list.append(nds)
    else:
        site_list = numpy.arange(len(site_x))
                
                
    for isite in site_list:

     
    
    
    
    
    

    if "read" in SetPrior.lower():
        halfspace ="halfspace_results"
        file, filext0 = os.path.splitext(file)
        prior_file = file+halfspace+filext0
        mod_prior, var_prior = inverse.load_prior(prior_file)


    start = time.time()
    """
    Loop over sites
    """
    sequence = range(nsite)
    if ReverseDir:
        sites = sequence[::-1]
    else:
        sites = sequence


    logsize = (2 + 7*Maxiter)
    site_log = numpy.full((len(sites),logsize), numpy.nan)


    for isite in site_list:

        print("\n Invert site #"+str(isite)+"/"+str(len(sites)))

        """
        Setup model-related parameter dict
        """

        if "read" in SetPrior.lower():
            mod_apr = mod_prior[isite]
            mod_ini = mod_apr.copy()

        elif "upd" in SetPrior:

            if isite == 0:
                mod_ini = mod_apr.copy()
                model = mod_ini.copy()
            else:
                mod_ini = model.copy()
                model = mod_ini.copy()

        elif "set" in SetPrior:

                mod_ini = mod_apr.copy()
                model = mod_ini.copy()

        Model = dict([
            ("m_act", mod_act),
            ("m_apr", mod_apr),
            ("m_var", mod_var),
            ("m_bnd", mod_bnd),
            ("m_ini", mod_ini)
            ])

        """
        Setup data-related parameter dict
        """

        dat_err = numpy.zeros_like(dat_obs)
        dat_err[isite, :], _ = inverse.set_errors(dat_obs[isite, :],
                                                daterr_add=DatErr_add,
                                                daterr_mult=DatErr_mult)

        Data = dict([
            ("d_act", dat_act[isite,:]),
            ("d_obs", dat_obs[isite,:]),
            ("d_err", dat_err[isite,:]),
            ("alt", site_alt[isite])
            ])

        """
        Call inversion algorithms
        """
        if "dram" in RunType.lower():
            
            from pymcmcstat.MCMC import MCMC
            from pymcmcstat.ParallelMCMC import ParallelMCMC
            import pymcmcstat
            from mcmc import dram_modelfun as modelfun
            from mcmc import dram_ssfun as ssfun
            
            

            mcset = MCMC()
            # Add data
            mcset.data.add_data_set(x, y)
            datestr = datetime.now().strftime('%Y%m%d_%H%M%S')
            savedir = 'resources' + os.sep + str('{}_{}'.format(datestr, 'parallel_chains'))
            mcset.simulation_options.define_simulation_options(
                nsimu=5.0e3, updatesigma=True, method='dram',
                savedir=savedir, savesize=1000, save_to_json=True,
                verbosity=0, waitbar=False, save_lightly=True, save_to_bin=True)
            mcset.model_settings.define_model_settings(sos_function=ssfun)
            
            for ipar in 
            mcset.parameters.add_model_parameter(name='m',
                                                 theta0=2.,
                                                 minimum=-10,
                                                 maximum=200,
                                                 sample=1)
            mcset.parameters.add_model_parameter(name='b',
                                                 theta0=2.75,
                                                 minimum=-10,
                                                 maximum=100,
                                                 sample=1)
         
            
        else:
            
            error("Runtype "+RunType+" does not exist! Exit.")


        # setup parallel MCMC
        parMC = ParallelMCMC()
        initial_values = np.array([[2.5, 2.5], [1.8, 3.8], [2.05, 3.42]])
        parMC.setup_parallel_simulation(mcset=mcset,
                                        initial_values=initial_values,
                                        num_chain=10,
                                    num_cores=11)
