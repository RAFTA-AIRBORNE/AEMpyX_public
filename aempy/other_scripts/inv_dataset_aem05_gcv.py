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
#       jupytext_version: 1.14.7
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
titstrng = util.print_title(version=version, fname=inspect.getfile(inspect.currentframe()), out=False)
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
    DatErr_add = 500.
    DatErr_mult = 0.02
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'


"""
input formats are ".npz",".nc4",".asc"
"""
Direction =  "normal"

FileList = "search"   # "set",  "read"
InDatDir =  AEMPYX_DATA + "/SGL2023/data/"
SearchStrng = "*k1*.npz"

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    # dat_files = []
    dat_files = numpy.load(AEMPYX_DATA + "/Projects/Compare/BundoranSubsets.npz")["setC"]
    dat_files = [os.path.basename(f) for f in dat_files]  
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
Output format is ".npz"
"""
OutFileFmt = ".npz"
OutResDir =  AEMPYX_DATA + "/SGL2023/results/"
print("Models written to dir: %s " % OutResDir)


if not os.path.isdir(OutResDir):
    print("File: %s does not exist, but will be created" % OutResDir)
    os.mkdir(OutResDir)


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
    RegVal1 =20.
    NTau1 =1
    Tau1min = numpy.log10(RegVal1)
    Tau1max = numpy.log10(RegVal1)

Tau1 = numpy.logspace(Tau1min, Tau1max, NTau1)
nreg = NTau0 * NTau1

"""
Model definition
"""

SetPrior = "set"
ParaTrans = 1

Nlyr = 36
DzStart = 5.
DzEnd = 10.
dz = numpy.logspace(numpy.log10(DzStart), numpy.log10(DzEnd), Nlyr)
z = numpy.append(0.0, numpy.cumsum(dz))


mod_act, mod_apr, mod_var, mod_bnd, m_state = inverse.init_1dmod(Nlyr)

mod_act[0*Nlyr:1*Nlyr] = 1
sizepar = numpy.shape(mod_act)
mpara = sizepar[0]

Guess_r = 500.0  # initial guess for resistivity in mod_apr
Guess_s = 0.3   # mod_std defines standard deviation of mod_apr
mod_apr[0*Nlyr:1*Nlyr] = Guess_r
mod_var[0*Nlyr:1*Nlyr] = numpy.power(Guess_s,2)
mod_apr[6*Nlyr:7*Nlyr-1] = dz[0:Nlyr - 1]
mod_var[6*Nlyr:7*Nlyr-1] = numpy.power(1.,2)


# mod_bnd = mumpy.array([])
max_val = 1.e+30
min_val = 1.e-30
# max_val = mod_apr[mod_act==1] + 3*mod_std[mod_act==1]
# mod_bnd[mod_act == 1, 1] = max_val
# min_val = mod_apr[mod_act==1] - 3*mod_std[mod_act==1]
# mod_bnd[mod_act == 1, 0] = min_val
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
    RegShift = 1
    Ctrl = dict([
        ("system", [AEM_system, FwdCall]),
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
                  var=mvar, sparse=True, thresh=0.05, inverse=True)
        Cmi=inverse.extract_cov(Cmi, mod_act)
        Cmi = scipy.sparse.block_diag([Cmi for Cmi in range(7)])
        CmiS=inverse.extract_cov(CmiS, mod_act)
        CmiS = scipy.sparse.block_diag([CmiS for Cmis in range(7)])
        C, sC = Cmi, CmiS
    else:
        InvSpace = "dat"
        Cm, CmS = inverse.covar(xc, yc, zc, covtype= ["exp", CorrL],
                  var=mvar, sparse=True, thresh=0.05, inverse=False)
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
    RegShift = 1
    Ctrl = dict([
        ("system",
         [AEM_system, FwdCall]),
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


if "rto" in RunType.lower():
    """
    Prepre parameters for rto (randomize-then-optimize) algorithm
    """
    NSamples = 100
    Percentiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.] # linear
    # Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68

    Ctrl["rto"] = [NSamples]
    Ctrl["output"] = numpy.array(["quant", Percentiles], dtype=object)  # ["ens "]

if "eki" in RunType.lower():
    """
    Prepare pramters for ensemble kalman inversion (eki)
    """
    NSamples = 100
    Percentiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.]  # linear
    # Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68

    Ctrl = dict([
        ("system", [AEM_system, FwdCall]),
        ("covar", numpy.array([C, sC], dtype=object)),
        ("transform", [DataTrans, ParaTrans]),
        ("eki", [NSamples]),
        ("output", numpy.array(["quant", Percentiles], dtype=object))
           ])

if "jac" in RunType.lower():
    """
    no new parameteres for jackknife estimates
    """
    Percentiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.] # linear
    # Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68
    Ctrl["output"] = numpy.array(["quant", Percentiles], dtype=object)  # ["ens "]
    pass

if "nul" in RunType.lower():
    """
    no new parametrres for nullspace projection
    """
    NSing = 3
    NSamples = 1000
    Percentiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.] # linear
    # Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68
    Ctrl["nullspace"] = [NSing, NSamples]
    Ctrl["output"] = numpy.array(["quant", Percentiles], dtype=object)  # ["ens "]
    pass



if OutInfo:
    print(Ctrl.keys())



outstrng = "_nlyr"+str(Nlyr)\
            +"_"+RunType\
            +"_"+RegFun\
            +"_Prior"+str(int(Guess_r))\
            +"_Err_a"+ str(int(DatErr_add))+"-m"+str(int(100*DatErr_mult))
print("ID string: input file + %s " % outstrng)


fcount =0
for file in dat_files:

    start = time.time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    fileout = OutResDir + name + outstrng

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
    if "reverse" in Direction.lower():
        sites = sequence[::-1]
    else:
        sites = sequence


    logsize = (2 + 7*Maxiter)
    site_log = numpy.full((len(sites),logsize), numpy.nan)

    for ii in sites:
        print("\n Invert site #"+str(ii)+"/"+str(len(sites)))
        
        Ctrl["name"] = str(fl_name)+"/"+str(ii)
        print("\n Data set: " + str(fl_name))


        """
        Setup model-related parameter dict
        """

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
        if "opt" in RunType.lower():
            Results =\
                alg.run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "occ" in RunType.lower():
            Results =\
                alg.run_tikh_occ(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "map" in RunType.lower():
            Results =\
                alg.run_map(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "rto" in RunType.lower():

            Results =\
                alg.run_rto(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "eki" in RunType.lower():
            Results =\
                alg.run_eki(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "jac" in RunType.lower():
            Results =\
                alg.run_jac(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        if "nul" in RunType.lower():
            Results =\
                alg.run_jac(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)


        """
        Store inversion Results
        """
        if OutInfo:
            print("Results: ",Results.keys())


        M = Results["model"]
        D = Results["data"]
        C = Results["log"]

        if ii==0:
            site_num  = numpy.array([ii])
            site_nrms = C[2]
            site_modl = M[0]
            site_merr = M[1]
            site_sens = M[2]
            site_dobs = D[0].reshape((1,-1))
            site_dcal = D[1].reshape((1,-1))
            site_derr = D[2].reshape((1,-1))
            # print("ii=0")
            cc = numpy.hstack((C[0], C[1],
                              C[2],
                              C[3].ravel(),
                              C[4].ravel(),
                              C[5].ravel(),
                              C[6].ravel()))
            site_log[ii,0:len(cc)] = cc
            if Uncert:
                jacd = Results["jacd"]
                site_jacd = jacd.reshape((1,numpy.size(jacd)))
                pcov = Results["cpost"]
                site_pcov = pcov.reshape((1,numpy.size(pcov)))
        else:
           site_num = numpy.vstack((site_num, ii))
           site_nrms = numpy.vstack((site_nrms, C[2]))
           site_modl = numpy.vstack((site_modl, M[0]))
           site_merr = numpy.vstack((site_merr, M[1]))
           site_sens = numpy.vstack((site_sens, M[2]))
           site_dobs = numpy.vstack((site_dobs, D[0]))
           site_dcal = numpy.vstack((site_dcal, D[1]))
           site_derr = numpy.vstack((site_derr, D[2]))
           cc = numpy.hstack((C[0], C[1],
                              C[2],
                              C[3].ravel(),
                              C[4].ravel(),
                              C[5].ravel(),
                              C[6].ravel()))
           site_log[ii,0:len(cc)] = cc

           if Uncert:
               jacd = Results["jacd"]
               site_jacd = numpy.vstack((site_jacd,jacd.reshape((1,numpy.size(jacd)))))
               pcov = Results["cpost"]
               site_pcov = numpy.vstack((site_pcov, pcov.reshape((1,numpy.size(pcov)))))



    numpy.savez_compressed(
        file=fileout+"_results.npz",
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
        site_dem=site_dem)

    if Uncert:
        numpy.savez_compressed(        
            file=fileout+"_results.npz",
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
            site_jacd= site_jacd,
            site_pcov= site_pcov)

    print("\n\nResults stored to "+fileout)
    elapsed = (time.time() - start)
    print (" Used %7.4f sec for %6i sites" % (elapsed, ii+1))
    print (" Average %7.4f sec/site\n" % (elapsed/(ii+1)))

print("\n\nAll done!")
