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
import scipy.stats

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
Data transformation is now allowed with three possible options:
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
    DatErr_add =  50.
    DatErr_mult = 0.03
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
input format is ".npz"
"""
InDatDir = AEMPYX_DATA + "/SynthData/data/"
FileList = "search"  # "search", "read"
SearchStrng = "*.npz"
# SearchStrng = "GEN*3LayerMod*Alt_120*.npz"

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = []

else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(dat_files)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")

"""
Output formats is ".npz"
"""
OutDatDir = AEMPYX_DATA + "/SynthData/results/"
print("Models written to dir: %s " % OutDatDir)


if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


"""
Define inversion type  optional additional parameters (e.g., Waveforms )
"""
RunType = "TikhOpt" # "TikhOcc",  "MAP_ParSpace", "MAP_DatSpace","Jack","DoI", "RTO""
Uncert = True

RegFun = "lcc" # "fix", "lcc", "gcv", "mle"
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
    RegVal1 =20.
    NTau1 =1
    Tau1min = numpy.log10(RegVal1)
    Tau1max = numpy.log10(RegVal1)

Tau1 = numpy.logspace(Tau1min, Tau1max, NTau1)
nreg = NTau0 * NTau1

EnsOut = True
Percentiles = [10., 20., 30., 40.] # linear
# Percentiles = [2.3, 15.9,]                   # 95/68

"""
Model definition
"""

SetPrior = "set"
ParaTrans = 1

Nlyr = 33
dzstart = 3.
dzend = 10.
dz = numpy.logspace(numpy.log10(dzstart), numpy.log10(dzend), Nlyr)
z = numpy.append(0.0, numpy.cumsum(dz))

zerolayer = numpy.zeros(Nlyr)

"""
Background model: default settings is rho only, - IP is nonexistent 
Neeeds to be adapted for reasonable IP
""" 
mod_act, mod_apr, mod_var, mod_bnd, m_state = inverse.init_1dmod(Nlyr)

mod_act[0*Nlyr:1*Nlyr] = 1
"""
For activating chargeability:
"""
# mod_act[2*Nlyr:3*Nlyr] = 1

"""
For activating Thicknesses (few layers only):
"""
# mod_act[6*Nlyr:7*Nlyr-1] = 1


sizepar = numpy.shape(mod_act)
mpara = sizepar[0]
"""
All activated parameter need reasonable priors
"""
Guess_rv = 100.0  # initial guess for resistivity in mod_apr
Guess_rs = 0.3    # std defines standard deviation 
mod_apr[0*Nlyr:1*Nlyr] = Guess_rv
mod_var[0*Nlyr:1*Nlyr] = numpy.power(Guess_rs,2)
Guess_mv = 0.5    # initial guess for chargeability in mod_apr
Guess_ms = 0.05   # std defines standard deviation 
mod_apr[2*Nlyr:3*Nlyr] = Guess_rv
mod_var[2*Nlyr:3*Nlyr] = numpy.power(Guess_rs,2)
"""
Thicknesses are kept constant (not activated)
"""
mod_apr[6*Nlyr:7*Nlyr-1] = dz[0:Nlyr - 1]
mod_var[6*Nlyr:7*Nlyr-1] = numpy.power(1.,2)


# mod_bnd = numpy.array([])
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
    print(" Initial halfspace resistivity of %6.2f Ohm*m" % (Guess_rv))
    print(" Log Standard error of %6.2f " % (Guess_rs))
    if not (mod_bnd == None) or (numpy.size(mod_bnd) == 0):
        print(" Upper limits: \n", mod_bnd[:, 1])
        print(" Lower limits: \n", mod_bnd[:, 0])

"""
Setup Controls for different Algorithms
"""
if "tikhopt" in  RunType.lower():
    """
    Prepare differential operator base methods for regularization matrices
    """
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
    Maxreduce = 5
    Rfact = 0.66
    LinPars = [Maxreduce, Rfact]


    ThreshRMS = [0.9, 1.0e-2, 1.0e-2]
    Delta = [1.e-5]
    RegShift = 2

    Ctrl = dict([
        ("system", [AEM_system, FwdCall]),
        ("name", ""),
        ("inversion",
         numpy.array([RunType, RegFun, Tau0, Tau1, Maxiter, ThreshRMS, 
                      LinPars, SetPrior, Delta, RegShift], dtype=object)),
        ("covar", 
         numpy.array([L0, Cm0, L1, Cm1], dtype=object)),
        ("transform",
         [DataTrans, ParaTrans]),
        ("uncert", 
         Uncert)
       ])

if "occ" in RunType.lower():
    """
    Prepare differential operator base methods for regularization matrices
    """
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
    Maxreduce = 5
    Rfact = 0.66
    LinPars = [Maxreduce, Rfact]

    Maxiter = 10
    Maxreduce = 5
    Rfact = 0.66
    ThreshRMS = [0.5, 1.0e-2, 1.0e-2]
    L = L1
    TauSeq = [0.5]
    Delta = [1.e-5]
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]),
        ("name", ""),
        ("inversion", 
         numpy.array([RunType, TauSeq, Tau0, Maxiter,ThreshRMS, 
                      LinPars, SetPrior, Delta],dtype=object)),
        ("covar", 
         numpy.array( [L0, Cm0, L1, Cm1], dtype=object)),
        ("transform",
         [DataTrans, ParaTrans]),
        ("uncert", 
         Uncert)
       ])

if "map" in  RunType.lower():

    """
    Prepare explicit covariances for MAP and related methods
    """

    zc = inverse.set_zcenters(dz)
    xc = numpy.zeros_like(zc)
    yc = numpy.zeros_like(zc)
    CorrL = numpy.array([30.0, 30.0, 30.0])

    """
    This setup is a workaround, correct only for rho-only inversion
    """

    mvar  = mod_var[0*Nlyr:1*Nlyr]
    # inverse.extract_mod(mod_var, mod_act)

    if "par"in RunType.lower():
        InvSpace = "par"
        Cmi, CmiS = inverse.covar(xc, yc, zc, covtype= ["exp", CorrL],
                  var=mvar, sparse=False, thresh=0.05, inverse=True)
        Cmi=inverse.extract_cov(Cmi, mod_act)
        Cmi = scipy.sparse.block_diag([Cmi for Cmi in range(7)])
        CmiS=inverse.extract_cov(CmiS, mod_act)
        CmiS = scipy.sparse.block_diag([CmiS for Cmis in range(7)])
        C, sC = Cmi, CmiS
    else:
        InvSpace = "dat"
        Cm, CmS = inverse.covar(xc, yc, zc, covtype= ["exp", CorrL],
                  var=mvar, sparse=False, thresh=0.05, inverse=False)
        Cm=inverse.extract_cov(Cm, mod_act)
        Cm = scipy.sparse.block_diag([Cm for Ci in range(7)])
        CmS=inverse.extract_cov(CmS, mod_act)
        CmS = scipy.sparse.block_diag([CmS for CmS in range(7)])
        C, sC = Cm, CmS

    Maxiter = 10
    Maxreduce = 5
    Rfact = 0.66
    ThreshRMS = [0.5, 1.0e-2, 1.0e-2]
    Delta = [1.e-5]
    TauSeq = [0.5]
    RegShift = 0
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]),
        ("name", ""),
        ("inversion",
         numpy.array([RunType, InvSpace, RegFun, Tau0, Tau1, Maxiter,ThreshRMS,
                      LinPars, SetPrior, Delta, RegShift], dtype=object)),
        ("covar",
         numpy.array([C, sC], dtype=object)),
        ("transform",
         [DataTrans, ParaTrans]),
        ("uncert",
         Uncert)
       ])




if OutInfo:
    print(Ctrl.keys())


OutStrng = "_nlyr"+str(Nlyr)\
            +"_"+RunType\
            +"_"+RegFun\
            +"_Prior"+str(int(Guess_rv))\
            +"_Err_a"+ str(int(DatErr_add))+"-m"+str(int(100*DatErr_mult))\
            +"_results"
print("ID string: input file + %s " % OutStrng)


fcount =0
for file in dat_files:

    start = time.time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    print("\n Reading file " + filein)


    fileout = (OutDatDir + name + OutStrng).replace("_results", "_ctrl")
    numpy.savez_compressed(file=fileout, **Ctrl)

    tmp = numpy.load(filein, allow_pickle=True)
    #  numpy.savez_compressed(
    # file=NPZFile, model = Model, data = Data, para = Para)
    imod_para = tmp["para"]
    imod_data = tmp["data"]
    imod_modl = tmp["model"]
    # print(numpy.shape(imod_data))

    imod_num = imod_data[:,0]
    imod_smp = imod_data[:,1]
    imod_alt = imod_data[:,2]


    [nsample,ndata] = numpy.shape(imod_data)


    if "read" in SetPrior.lower():
        halfspace ="halfspace_results"
        file, filext0 = os.path.splitext(file)
        prior_file = file+halfspace+filext0
        mod_prior, mod_var = inverse.load_prior(prior_file)


    dat_act = numpy.tile(data_active,(nsample,1))
    dat_obs = imod_data[:,3:]
    dat_err = numpy.zeros_like(dat_obs)
    # print(numpy.shape(dat_obs))

    """
    Loop over samples

    """
    sequence = range(nsample)
    samples = sequence


    logsize = (2 + 7*Maxiter)
    imod_log = numpy.full((len(samples),logsize), numpy.nan)

    start = time.time()
    for ii in samples:
        print("\n Invert sample #"+str(ii))

        """
        Setup model-related paramter
        """
        if ii==0:
            mod_true = imod_modl
            dat_true = dat_obs[ii, :]

        if "read" in SetPrior.lower():
            mod_apr = mod_prior[ii]
            mod_ini = mod_apr.copy()

        elif "upd" in SetPrior:

            if ii == 0:
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
        Setup data-related paramter
        """
        dat_err[ii, :], _ = inverse.set_errors(dat_obs[ii, :],
                                                daterr_add=DatErr_add,
                                                daterr_mult=DatErr_mult)

        Data = dict([
            ("d_act", dat_act[ii,:]),
            ("d_obs", dat_obs[ii,:]),
            ("d_err", dat_err[ii,:]),
            ("alt", imod_alt[ii])
            ])

        """
        Call inversion algorithms
        """
        if "opt" in RunType.lower():
            results =\
                alg.run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "occ" in RunType.lower():
            results =\
                alg.run_tikh_occ(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "map" in RunType.lower():
            results =\
                alg.run_map(Ctrl=Ctrl, Model=Model, Data=Data,

                                  OutInfo=OutInfo)
        """
        Store inversion Results
        """
        if OutInfo:
            print("Results: ",results.keys())


        M = results["model"]
        D = results["data"]
        C = results["log"]

        if ii==0:
            ens_num  = numpy.array([ii])
            ens_nrms = C[2]
            ens_modl = M[0]
            ens_merr = M[1]
            ens_sens = M[2]
            ens_dobs = D[0].reshape((1,-1))
            ens_dcal = D[1].reshape((1,-1))
            ens_derr = D[2].reshape((1,-1))
 
        else:
           ens_num = numpy.vstack((ens_num, ii))
           ens_nrms = numpy.vstack((ens_nrms, C[2]))
           
           ens_modl = numpy.vstack((ens_modl, M[0]))
           ens_merr = numpy.vstack((ens_merr, M[1]))
           ens_sens = numpy.vstack((ens_sens, M[2]))
           ens_dobs = numpy.vstack((ens_dobs, D[0]))
           ens_dcal = numpy.vstack((ens_dcal, D[1]))
           ens_derr = numpy.vstack((ens_derr, D[2]))
 

    m_quants, m_mean, m_stdv, m_skew, m_kurt, m_mode = \
        inverse.calc_stat_ens(ensemble=ens_modl, quantiles=Percentiles, sum_stats=True)
    stat_modl = numpy.vstack((m_quants, m_mean, m_stdv, m_skew, m_kurt, m_mode))

    d_quants, d_mean, d_stdv, d_skew, d_kurt, d_mode = \
        inverse.calc_stat_ens(ensemble=ens_modl, quantiles=Percentiles, sum_stats=True)
    stat_dcal = numpy.vstack((d_quants, d_mean, d_stdv, d_skew, d_kurt, d_mode))

    r_quants, r_mean, r_stdv, r_skew, r_kurt, r_mode = \
        inverse.calc_stat_ens(ensemble=ens_nrms, quantiles=Percentiles, sum_stats=True) 
    stat_nrms= numpy.vstack((r_quants, r_mean, r_stdv, r_skew, r_kurt, r_mode))

    mod_alt =  imod_alt[0]
    
    if EnsOut:
        fileout = OutDatDir + name + OutStrng +".npz"
        numpy.savez_compressed(
            file=fileout,
            header= titstrng,
            aem_system=AEM_system,
            mod_true=mod_true, dat_true=dat_true, mod_alt=mod_alt,
            mod_ref=mod_apr,
            mod_act=mod_act,
            dat_act=dat_act,            
            ens_modl=ens_modl,
            ens_merr=ens_merr,
            ens_dobs=ens_dobs,
            ens_dcal=ens_dcal,
            ens_derr=ens_derr,
            ens_nrms=ens_nrms,            
            stat_dcal=stat_dcal,
            stat_modl=stat_modl,
            stat_nrms=stat_nrms)
    else:
        fileout = OutDatDir + name + OutStrng +"_stat.npz"
        numpy.savez_compressed(
            file=fileout,
            header= titstrng,
            aem_system=AEM_system,
            mod_true=mod_true, dat_true=dat_true, mod_alt=mod_alt,
            mod_ref=mod_apr,
            mod_act=mod_act,
            dat_act=dat_act,            
            stat_dcal=stat_dcal,
            stat_modl=stat_modl,
            stat_nrms=stat_nrms)


print("\n\nResults stored to "+fileout)
elapsed = (time.time() - start)
print (" Used %7.4f sec for %6i samples" % (elapsed, ii+1))
print (" Average %7.4f sec/imod\n" % (elapsed/(ii+1)))

print("\n\nAll done!")
