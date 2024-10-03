#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,ipynb
#     main_language: python
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (Spyder)
#     language: python3
#     name: python3
# ---


import os
import sys
from sys import exit as error
# from datetime import datetime
# from time import process_time
# from random import randrange
# import time
# import warnings
# import inspect
# import copy

import numpy
import scipy

import multiprocessing
# from numba import njit

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse
# import alg
# -

AEMPYX_DATA = os.environ["AEMPYX_DATA"]
rng = numpy.random.default_rng()
nan = numpy.nan  

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = False

Parallel = True
if Parallel:
    import parallel
    
    Njobs = 17

    if Njobs<0:
        Njobs=multiprocessing.cpu_count()
    else:
        Njobs=min(Njobs, multiprocessing.cpu_count())

    print(str(Njobs)+" processors will be used in parallel")

else:
    Njobs = 1
# -

# The following cell gives values to AEM-system related settings. 
#
# Data transformation is activated by the variable _DataTrans_. Currently 
# three possible options are allowed: _DataTrans = 0_: No transformation, 
# i.e., the raw data are used. _DataTrans = 1_: The natural log of data 
# is taken, only allowed for strictly positive values. _DataTrans = 2_: 
# If data scale logarithmically, an _asinh_ transformation (introduced by
# Scholl, 2000) is applied. It allows negatives, which may occur in TDEM, 
# when IP effects are present.
#        
# A general additive/multiplicative error model is applied on the raw data
# before transformation, and errors are also transformed.

# +
# AEM_system = "genesis"
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  70.
    DatErr_mult = 0. #03
    data_active = numpy.ones(NN[2], dtype="int8")
    # data_active[0] = 0   # real at 900Hz

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' horizontals'

ReverseDir = False


FileList = "search"  # "search", "read"
SearchStrng = "*FL*k2.npz"

AEMPYX_DATA  = "/home/vrath/Mohammednur/"
InDatDir =  AEMPYX_DATA + "/proc/"
if not InDatDir.endswith("/"): InDatDir=InDatDir+"/"
# +
"""
Output format is ".npz"
"""
OutFileFmt = ".npz"
OutResDir =   AEMPYX_DATA + "/results/"
if not OutResDir.endswith("/"): OutResDir=OutResDir+"/"
print("Models written to dir: %s " % OutResDir)
if not os.path.isdir(OutResDir):
    print("File: %s does not exist, but will be created" % OutResDir)
    os.mkdir(OutResDir)


if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    # dat_files = []
    dat_files = [InDatDir+"StGormans_FL11379-0_raw.npz"]
    # dat_files =  numpy.load(AEMPYX_DATA + "/Projects/Compare/BundoranSubsets.npz")["setC"]
    
    dat_files = [os.path.basename(f) for f in dat_files]  
else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              fullpath=False, out= True, sort=True)


ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")
if Njobs<=0:
    Njobs= min(Njobs, ns)

# +
"""
Define inversion type  optional additional parameters (e.g., Waveforms )
"""

RunType = "TikhOpt" # "TikhOcc",  "MAP_ParSpace", "MAP_DatSpace","Jack","DoI", "RTO""
Uncert = True
RegFun = "gcv" # "fix", "lcc", "gcv", "mle"
RegVal0 = 1.e-5
NTau0 = 1
Tau0min = numpy.log10(RegVal0)
Tau0max = numpy.log10(RegVal0)
Tau0 = numpy.logspace(Tau0min, Tau0max, NTau0)

if any(s in RegFun.lower() for s in ["gcv", "upr", "ufc", "mle", "lcc"]):
    RegVal1Min = 0.1
    RegVal1Max = 1000.
    NTau1 =64
    Tau1min = numpy.log10(RegVal1Min)
    Tau1max = numpy.log10(RegVal1Max)
else:
    RegVal1 =0.00001
    NTau1 =1
    Tau1min = numpy.log10(RegVal1)
    Tau1max = numpy.log10(RegVal1)

Tau1 = numpy.logspace(Tau1min, Tau1max, NTau1)
nreg = NTau0 * NTau1

# +
"""
Model definition
"""

SetPrior = "set"
ParaTrans = 1

Nlyr = 21
dzstart = 4.
dzend = 10.
dz = numpy.logspace(numpy.log10(dzstart), numpy.log10(dzend), Nlyr)
print(dz)
z = numpy.append(0.0, numpy.cumsum(dz))
print(z)


mod_act, mod_apr, mod_var, mod_bnd, m_state = inverse.init_1dmod(Nlyr)

mod_act[0*Nlyr:1*Nlyr] = 1
sizepar = numpy.shape(mod_act)
mpara = sizepar[0]

Guess_r = 100.0  # initial guess for resistivity in mod_apr
# Guess_r = 10.0    # low value for DoI estimate
# Guess_r = 1000.0  # high value for DoI estimate
Guess_s = 0.3   # mod_std (logscale) defines standard deviation of mod_apr
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
# -

# Setup controls for different slgorithms, here in particular prepare 
# differential operator base methods for regularization matrices

# +
if "tikhopt" in  RunType.lower():

    D0 = inverse.diffops(dz, der=False, mtype="sparse", otype="L0")
    L = [D0 for D in range(7)]
    L0 = scipy.sparse.block_diag(L)
    Cm0 = L0.T@L0
    Cm0 = inverse.extract_cov(Cm0, mod_act)

    D1 = inverse.diffops(dz, der=False, mtype="sparse", otype="L1")
    L = [D1 for D in range(7)]
    L1 = scipy.sparse.block_diag(L)
    Cm1 = L1.T@L1
    Cm1 = inverse.extract_cov(Cm1, mod_act)

    Maxiter = 10
    Maxreduce = 10
    Rfact = 0.66
    LinPars = [Maxreduce, Rfact]

    # ThreshFit = [0.9, 1.0e-2, 1.0e-2, "rms"]
    ThreshFit = [0.5, 1.0e-2, 1.0e-2, "rms"]
    Delta = [1.e-5]
    RegShift = -1


    ctrl_dict ={
        "system":
            [AEM_system, FwdCall],
        "header":
            [titstrng, ""],
        "inversion":
            numpy.array([RunType, RegFun, Tau0, Tau1, Maxiter, ThreshFit,
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

if OutInfo:
    print(ctrl_dict.keys())
# -

outstrng =  "_t"+RunType+"_q"+ThreshFit[3]+"_r"+RegFun+\
            "_a"+str(round(DatErr_add,0))+"_m"+str(round(DatErr_mult*100,0))
print("ID string: input file + %s " % outstrng)




if Parallel:
    import joblib
    # from joblib import Parallel, delayed, parallel_config
    joblib.Parallel(n_jobs=Njobs, verbose=100)(
        joblib.delayed(parallel.run_tikh_flightline)(ctrl=ctrl_dict, 
                                                     data_dir=InDatDir,
                                                     data_file=filin,
                                                     result_dir=OutResDir, 
                                                     result_strng=outstrng) for filin in dat_files)
else:
    for filin in dat_files:
        _ = parallel.run_tikh_flightline(ctrl=ctrl_dict, 
                                         data_dir=InDatDir,
                                         data_file=filin,
                                         result_dir=OutResDir, 
                                         result_strng=outstrng)

print("\n\nAll done!")