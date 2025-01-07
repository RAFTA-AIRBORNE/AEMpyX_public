#!/usr/bin/env python3
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


import os
import sys
from sys import exit as error
from datetime import datetime
# from time import process_time
# from random import randrange
import time
import warnings
# import inspect
import copy

import numpy
import scipy

# %logstart -o

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse
import alg

warnings.simplefilter(action="ignore", category=FutureWarning)

AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = False
now = datetime.now()


"""
System related settings.
Data transformation is allowed with three possible options:
DataTrans   = 0           raw data
            = 1           natural log of data
            = 2           asinh transformation
An error model is applied for the raw data, which is
mixed additive/multiplicative. in case of data transformation,
errors are also transformed.
"""
# AEM_system = "genesis"
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  75.
    DatErr_mult = 0.05
    data_active = numpy.ones(NN[2], dtype="int8")


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'


"""
input formats are ".npz",".nc4",".asc"
"""
Direction =  "normal"

InDatDir =  AEMPYX_DATA + "/Projects/StGormans/proc_delete_PLM3s/"
print("Data files read from dir:  %s" % InDatDir)
dat_files = [""]

ns = len(dat_files)
if ns ==0:
    error("No data files set!. Exit.")

SampleType="dist"           # "dist", "site", "rand"
if "dist" in SampleType.lower():
    SampleSites = [[500., 1000.],]

if "site" in SampleType.lower():
    SampleSites = [[ 300 ,  500],]

if "rand" in SampleType.lower():
        
    NumSites = [3, 2]  

"""
Output format is ".npz"
"""
OutFileFmt = ".npz"
OutDatDir =  AEMPYX_DATA + "/Projects/Compare/results/New/"
print("Models written to dir: %s " % OutDatDir)


if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


"""
Define inversion type  optional additional parameters (e.g., Waveforms )
"""
RunType = "emcee" # "dram",  "dream", "hmc", "mhsimple"


"""
Model and prior covariance definition
"""

Nlyr = 36
DzStart = 5.
DzEnd = 10.
dz = numpy.logspace(numpy.log10(DzStart), numpy.log10(DzEnd), Nlyr)
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

zc = inverse.set_zcenters(dz)
xc = numpy.zeros_like(zc)
yc = numpy.zeros_like(zc)
CorrL = numpy.array([30.0, 30.0, 30.0])

"""
This setup is a workaround, correct only for rho-only inversion
"""

mvar  = mod_var[0*Nlyr:1*Nlyr]
# inverse.extract_mod(mod_var, mod_act)

InvSpace = "dat"
Cm, CmS = inverse.covar(xc, yc, zc, covtype= ["exp", CorrL],
          var=mvar, sparse=False, thresh=0.05, inverse=False)
Cm=inverse.extract_cov(Cm, mod_act)
Cm = scipy.sparse.block_diag([Cm for Ci in range(7)])

"""
Setup Controls for different Algorithms
"""
if "emcee" in  RunType.lower():
    NumSample = 10000
    NumChains = 8
    
 
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]), 
        ("transform",
         [DataTrans, ParaTrans])
        ("inversion",
         numpy.array([NumSample, NumChains], dtype=object)),
        ("prior", 
         numpy.array([Cm], dtype=object)),
       ])

if "dram" in  RunType.lower():
    NumSample = 10000
    NumChains = 8
    
 
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]), 
        ("transform",
         [DataTrans, ParaTrans])
        ("inversion",
         numpy.array([NumSample, NumChains], dtype=object)),
        ("prior", 
         numpy.array([Cm], dtype=object)),
       ])

if "dream" in  RunType.lower():
    NumSample = 10000
    NumChains = 8
    
 
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]), 
        ("transform",
         [DataTrans, ParaTrans])
        ("inversion",
         numpy.array([NumSample, NumChains], dtype=object)),
        ("prior", 
         numpy.array([Cm], dtype=object)),
       ])

if "hmc" in  RunType.lower():
    NumSample = 10000
    NumChains = 8
    
 
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]), 
        ("transform",
         [DataTrans, ParaTrans])
        ("inversion",
         numpy.array([NumSample, NumChains], dtype=object)),
        ("prior", 
         numpy.array([Cm], dtype=object)),
       ])

if "mh" in  RunType.lower():
    NumSample = 10000
    NumChains = 8
    
 
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]), 
        ("transform",
         [DataTrans, ParaTrans])
        ("inversion",
         numpy.array([NumSample, NumChains], dtype=object)),
        ("prior", 
         numpy.array([Cm], dtype=object)),
       ])



if OutInfo:
    print(Ctrl.keys())



OutStrng = "_nlyr"+str(Nlyr)\
            +"_"+RunType\
            +"_results"
print("ID string: input file + %s " % OutStrng)


fcount = -1
for file in dat_files:

    start = time.time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    print("\n Reading file " + filein)

    fileout = OutDatDir + name + OutStrng
    numpy.savez_compressed(file=fileout.replace("_results","_ctrl")+OutFileFmt,**Ctrl)


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

    start = time.time()
    
      
    """
    construct site list
    """

    site_list = SampleSite[fcount]
   
    sitex = site_x - site_x[0]
    sitey = site_y - site_y[0]
    siter = numpy.sqrt(numpy.power(sitex, 2.0) + numpy.power(sitey, 2.0))
    
    sites = []
    
    if "site" in SampleType.lower():
        sites = site_list
        
    elif "dist" in SampleType.lower():
        site_x = site_x - site_x[0]
        site_y = site_y - site_y[0]
        site_r = numpy.sqrt(numpy.power(site_x, 2.0) + numpy.power(site_y, 2.0))
        for nid in numpy.arange(len(site_list)):
            nds = (numpy.abs(site_list[nid] - site_r)).argmin()
            sites.append(nds)
    
    elif "rand" in Sample:
        sites = random.sample(range(len(site_x)), NSamples)


    else:
        sites = numpy.arange(len(site_x))



    """
    Loop over site_list
    """

    for ii in sites:
        
        print("\n Invert site #"+str(ii)+"/"+str(len(site_list)))

        """
        Setup model-related parameter dict
        """

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
        dat_err[ii, :], _ = inverse.set_errors(dat_obs[ii, :],
                                                daterr_add=DatErr_add,
                                                daterr_mult=DatErr_mult)

        Data = dict([
            ("d_act", dat_act[ii,:]),
            ("d_obs", dat_obs[ii,:]),
            ("d_err", dat_err[ii,:]),
            ("alt", site_alt[ii])
            ])

        """
        Call inversion algorithms
        """

        if "emc" in RunType.lower():
            results =\
                alg.run_EMCEE(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)
        if "dram" in RunType.lower():
            results =\
                alg.run_DRAM(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)
                
        if "dream" in RunType.lower():
            results =\
                alg.run_DREAM(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)
        if "mh" in RunType.lower():
            results =\
                alg.run_MH(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

                
       if "hmc" in RunType.lower():
           results =\
               alg.run_HMC(Ctrl=Ctrl, Model=Model, Data=Data,
                                 OutInfo=OutInfo)




        """
        Store inversion Results
        """
        if OutInfo:
            print("Results: ",results.keys())

            
            if "ens" in Ctrl["output"]:
               site_rto_ens = numpy.vstack((site_rto_ens, model_ens))


    numpy.savez_compressed(
        file=Fileout+OutFileFmt,
        fl_data=file,
        fl_name=fl_name,
        header=titstrng,
        site_log =site_log,
        mod_ref=mod_apr,
        mod_act=mod_act,
        dat_act=dat_act,
        site_modl=site_modl,
        site_sens=site_sens,
        site_merr=site_merr,
        site_dobs=site_dobs,
        site_dcal=site_dcal,
        site_derr=site_derr,
        site_nrms=site_nrms,
        site_num=site_num,
        site_y=site_y,
        site_x=site_x,
        site_gps=site_gps,
        site_alt=site_alt,
        site_dem=site_dem,            
        site_jacd=site_jacd,
        site_pcov=site_pcov,
        site_rto_avg=site_rto_avg,
        site_rto_var=site_rto_var,
        site_rto_med=site_rto_med,
        site_rto_mad=site_rto_mad,
        site_rto_prc=site_rto_prc)
           
    if "ens" in Ctrl["output"]:
        util.add_object_npz(filein=Fileout+OutFileFmt,
                   xkeys=["site_rto_ens"], xobjects=[site_rto_ens])

    print("\n\nResults stored to "+Fileout)
    elapsed = (time.time() - start)
    print (" Used %7.4f sec for %6i site_list" % (elapsed, ii+1))
    print (" Average %7.4f sec/site\n" % (elapsed/(ii+1)))

print("\n\nAll done!")

