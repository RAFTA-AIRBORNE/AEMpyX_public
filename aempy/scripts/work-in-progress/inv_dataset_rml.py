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
import inspect
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
    DatErr_add = 300.
    DatErr_mult = 0.03
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'


"""
input formats are ".npz",".nc4",".asc"
"""


FileList = "search"  # "search", "read"
InDatDir =  AEMPYX_DATA + "/Projects/StGormans/proc_delete_PLM3s/"
SearchStrng = "A1_rect_StGormans_FL11379-0_proc_delete_PLM3s_k3.npz"

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
OutResDir =  AEMPYX_DATA + "/Projects/StGormans/results_rto/"
print("Models written to dir: %s " % OutResDir)


if not os.path.isdir(OutResDir):
    print("File: %s does not exist, but will be created" % OutResDir)
    os.mkdir(OutResDir)


"""
script offers several methods do choose sites:
Sample = 
"rand"                  choose Nsample random sites
"step"                  choose sites with Start/Sop/Step
"list suboption"        define list, with suboptions position and distance

Any other string will choose full data set.

"""

Sample = ["random"]   # 

if "rand" in Sample[0].lower():
    Num_samples = 10
    Sample.append(Num_samples)
    
elif "step" in Sample[0].lower():
    Start, Stop, Step = 0, -1, 10
    Sample.extend((Start, Stop, Step))

elif "list" in Sample[0].lower():
    
    if "pos" in Sample[0].lower():
        Samplist = [100, 200]
    else:
        Distlist = [ 1500.]


Direction =  "normal"

"""
Define inversion type  optional additional parameters (e.g., Waveforms )
"""
RunType = "TikhOpt RTO" # "TikhOcc",  "MAP_ParSpace", "MAP_DatSpace","Jack","DoI", "RTO""
Uncert = True
RegFun = "gcv" # "fix", "lcc", "gcv", "mle"
RegVal0 = 1.e-5
NTau0 = 1
Tau0min = numpy.log10(RegVal0)
Tau0max = numpy.log10(RegVal0)
Tau0 = numpy.logspace(Tau0min, Tau0max, NTau0)

if any(s in RegFun.lower() for s in ["gcv", "upr", "ufc", "mle", "lcc"]):
    RegVal1Min = 0.1
    RegVal1Max = 10000.
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
        ("system", [AEM_system, FwdCall]),  ("name", ""),
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

if "map" in  RunType.lower():

    """
    Prepare explicit covariances for MAP and related methods
    """

    zc = inverse.set_zcenters(dz)
    xc = numpy.zeros_like(zc)
    yc = numpy.zeros_like(zc)

    """
    This setup is a workaround, correct only for rho-only inversion
    """
    Tau = Tau1
    CorrL = 40.
    mvar  = mod_var[0*Nlyr:1*Nlyr]
    # inverse.extract_mod(mod_var, mod_act)
  
    if "par"in RunType.lower():
        InvSpace = "par"
        Cmi, CmiS = inverse.covar(xc, yc, zc, covtype= ["exp", CorrL],
                  var=mvar, sparse=False, thresh=0.05, inverse=True)

        
        Cmi = inverse.full_cov([Cmi])
        Cmi = inverse.extract_cov(Cmi, mod_act)
        C = scipy.sparse.csr_matrix(Cmi)
        
  
        CmiS = inverse.full_cov([CmiS])
        CmiS = inverse.extract_cov(CmiS, mod_act)        
        sC = scipy.sparse.csr_matrix(CmiS)

        print(numpy.shape(C),numpy.shape(sC))
    else:
        InvSpace = "dat"
        Cm, CmS = inverse.covar(xc, yc, zc, covtype= ["exp", CorrL],
                  var=mvar, sparse=False, thresh=0.05, inverse=False)

        Cm = inverse.full_cov([Cm])
        Cm = inverse.extract_cov(Cm, mod_act)
        C = scipy.sparse.csr_matrix(Cm)
        
        CmS = inverse.full_cov([CmS])
        CmS = inverse.extract_cov(CmS, mod_act)
        sC = scipy.sparse.csr_matrix(CmS)

        
    Maxiter = 10
    Maxreduce = 5
    Rfact = 0.66
    LinPars = [Maxreduce, Rfact]
        

    ThreshRMS = [0.5, 1.0e-2, 1.0e-2]
    Delta = [1.e-5]
    TauSeq = [0.5]
    RegShift = 1
    Ctrl = dict([
        ("system",
         [AEM_system, FwdCall]),
        ("name", ""),
        ("inversion",
         numpy.array([RunType, InvSpace, RegFun, Tau, Maxiter,ThreshRMS,
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
    Percentiles = numpy.array([10., 20., 30., 40., 50., 60., 70., 80., 90.]) # linear
    # Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68
    Ctrl["rto"] =  numpy.array([NSamples, Percentiles], dtype=object)
    Ctrl["output"] = ["ens "]

if OutInfo:
    print(Ctrl.keys())



OutStrng = "_RTO_nlyr"+str(Nlyr)\
            +"_"+RunType\
            +"_"+RegFun\
            +"_Prior"+str(int(Guess_r))\
            +"_results"
print("ID string: input file + %s " % OutStrng)


fcount =0
for file in dat_files:

    start = time.time()

    start = time.time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    print("\n Reading file " + filein)

    Fileout = OutResDir + name + OutStrng

    numpy.savez_compressed(file=Fileout.replace("_results","_ctrl"), **Ctrl)

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
    fl_orig = [site_x[0], site_y[0]]
    
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
        site_list = sequence[::-1]
    else:
        site_list = sequence
        
        
    site_x = site_x - site_x[0]
    site_y = site_y - site_y[0]
    site_r = numpy.sqrt(numpy.power(site_x, 2.0) + numpy.power(site_y, 2.0))
    
    site_list = numpy.arange(len(site_x)) 
                             
    if "rand" in Sample[0].lower() or "step" in Sample[0].lower():
        site_list = util.sample_list(site_list, method=Sample)
        

    elif "list" in Sample[0].lower():
        
        if "posi" in Sample[0].lower():
            site_list = Samplist
        if "dist" in Sample[0].lower():
            for nid in numpy.arange(len(Distlist)):
                nds = (numpy.abs(Distlist[nid] - site_r)).argmin()
                site_list.append(nds)
   




    logsize = (2 + 7*Maxiter)
    site_log = numpy.full((len(site_list),logsize), numpy.nan)

    for ii in site_list:

        print("\n Invert site #"+str(ii)+"/"+str(len(site_list)))
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

        if "rto" in RunType.lower():
            results =\
                alg.run_RML(Ctrl=Ctrl, Model=Model, Data=Data,
                                  OutInfo=OutInfo)

        """
        Store inversion Results
        """
        if OutInfo:
            print("Results: ",results.keys())
       
        M = results["model"]
        D = results["data"]
        C = results["log"]            

        jacd = results["jacd"]
        pcov = results["cpost"]

        rto_avg = results["rto_avg"]
        rto_var = results["rto_var"]
        rto_med = results["rto_med"]
        rto_mad = results["rto_mad"]
        rto_prc = results["rto_prc"]
        
        if "ens" in Ctrl["output"]:
           rto_ens = results["rto_ens"]
        
        if ii==0:
            site_num  = numpy.array([ii])
            site_nrms = C[2]
            site_modl = M[0]
            site_merr = M[1]
            site_sens = M[2]
            site_dobs = D[0].reshape((1,-1))
            site_dcal = D[1].reshape((1,-1))
            site_derr = D[2].reshape((1,-1))
            site_jacd = jacd.reshape((1,numpy.size(jacd)))
            site_pcov = pcov.reshape((1,numpy.size(pcov)))
            site_rto_avg = rto_avg
            site_rto_var = rto_var
            site_rto_med = rto_med            
            site_rto_mad = rto_mad
            site_rto_prc = rto_prc
            
            if "ens" in Ctrl["output"]:
                site_rto_ens = rto_ens
                
                
        else:
            site_num = numpy.vstack((site_num, ii))
            site_nrms = numpy.vstack((site_nrms, C[2]))
            site_modl = numpy.vstack((site_modl, M[0]))
            site_merr = numpy.vstack((site_merr, M[1]))
            site_sens = numpy.vstack((site_sens, M[2]))
            site_dobs = numpy.vstack((site_dobs, D[0]))
            site_dcal = numpy.vstack((site_dcal, D[1]))
            site_derr = numpy.vstack((site_derr, D[2]))               
            site_jacd = numpy.vstack((site_jacd, jacd.reshape((1,numpy.size(jacd)))))
            site_pcov = numpy.vstack((site_pcov, pcov.reshape((1,numpy.size(pcov)))))
            site_rto_avg = numpy.vstack((site_rto_avg, rto_avg))
            site_rto_var = numpy.vstack((site_rto_var, rto_var))
            site_rto_med = numpy.vstack((site_rto_med, rto_med))  
            site_rto_mad = numpy.vstack((site_rto_mad, rto_mad))
            site_rto_prc = numpy.vstack((site_rto_prc, rto_prc))
            
            if "ens" in Ctrl["output"]:
               site_rto_ens = numpy.vstack((site_rto_ens, rto_ens))


    numpy.savez_compressed(
        file=Fileout+OutFileFmt,
        fl_data=file,
        fl_name=fl_name,
        fl_orig = fl_orig,
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
    print (" Used %7.4f sec for %6i sites" % (elapsed, ii+1))
    print (" Average %7.4f sec/site\n" % (elapsed/(ii+1)))

print("\n\nAll done!")

