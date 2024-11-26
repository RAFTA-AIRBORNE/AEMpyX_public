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
from time import process_time
# from random import randrange
# import time
import warnings
# import inspect
import copy
import getpass

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
AEMPYX_DATA = "/home/vrath/work/A1_StGormans/"

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
input formats are ".npz",".nc4",".asc"
"""
InFileFmt = ".npz"
AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/A1_StGormans/"
FileList = "search"
SearchStrng = "*FL*k3*data.npz"
InDatDir =  AEMPYX_DATA + "/proc/"
if not InDatDir.endswith("/"): InDatDir=InDatDir+"/"


"""
Output format is ".npz"
"""
OutFileFmt = ".npz"
OutResDir =   AEMPYX_DATA + "/results_halfspace/"
if not os.path.isdir(OutResDir):
    print("File: %s does not exist, but will be created" % OutResDir)
    os.mkdir(OutResDir)
print("Models written to dir: %s " % OutResDir)


if "set" in FileList.lower():
    # InDatDir = AEMPYX_DATA + "/ERT_AEM_Profiles/data/"
    InDatDir ="/home/vrath/work/AEM_Data/Tellus/data/SYNTH/"
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = ["SYNTH_AEM05_1Layer_Resistor.asc",
                "SYNTH_AEM05_1Layer_Conductor.asc",]

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
Define inversion type  optional additional parameters (e.g., Waveforms )
"""
RunType = "TikhOpt" # "TikhOcc",  "MAP_ParSpace", "MAP_DatSpace","Jack","DoI", "RTO""
Uncert = True
Direction =  "normal"


RegFun = "fix" # "fix", "lcc", "gcv", "mle"
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

NSamples = 100
Percentiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.] # linear
# Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68

"""
Model definition
"""

SetPrior = "set"
ParaTrans = 1

Nlyr = 1
dzstart = 3.
dzend = 10.
dz = numpy.logspace(numpy.log10(dzstart), numpy.log10(dzend), Nlyr)
z = numpy.append(0.0, numpy.cumsum(dz))

zerolayer = numpy.zeros(Nlyr)
izerolayer = zerolayer.astype(int)
onelayer = numpy.ones(Nlyr)
ionelayer = onelayer.astype(int)

mod_act, mod_apr, mod_var, mod_bnd, m_state = inverse.init_1dmod(Nlyr)

mod_act[0*Nlyr:1*Nlyr] = ionelayer[0:Nlyr]
sizepar = numpy.shape(mod_act)
mpara = sizepar[0]

Guess_r = 300.0  # initial guess for resistivity in mod_apr
Guess_s = 0.3   # mod_std defines standard deviation of mod_apr
mod_apr[0*Nlyr:1*Nlyr] = Guess_r * onelayer[0:Nlyr]
mod_var[0*Nlyr:1*Nlyr] = numpy.power(Guess_s * onelayer[0:Nlyr],2)
mod_apr[6*Nlyr:7*Nlyr-1] = dz[0:Nlyr - 1]
mod_var[6*Nlyr:7*Nlyr-1] = numpy.power(1. * onelayer[0:Nlyr-1],2)


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

    ThreshFit = [0.9, 1.0e-2, 1.0e-2, "rms"]
    # ThreshFit = [5., 1.0e-2, 1.0e-2, "smp"]
    Delta = [1.e-5]
    RegShift = 0


    Ctrl ={
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
            numpy.array([DataTrans, data_active, DatErr_add, DatErr_mult, Direction], dtype=object),
        "model":
            numpy.array([ParaTrans, mod_act, mod_apr, mod_var, mod_bnd], dtype=object),
                }

if OutInfo:
    print(Ctrl.keys())


OutStrng = "_halfspace"
print("ID string: input file + %s " % OutStrng)


fcount =0
for file in dat_files:

    start = process_time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    print("\n Reading file " + filein)

    fileout = OutResDir + name + OutStrng
    numpy.savez_compressed(file=fileout+"_ctrl"+OutFileFmt,**Ctrl)


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


    if "read" in SetPrior.lower():
        error("Prior model read not yet implemented! Exit.")
        # file,filext0 = os.path.splitext(file)
        # tmp = numpy.load(file, allow_pickle=True)
        # mod_aprel = tmp["site_modl"]
        # mod_apr = numpy.mat(inverse.extract_mod(mod_apr,mod_act))
    elif "set" in SetPrior.lower() or "upd" in SetPrior.lower():
        mod_apr = mod_apr.copy()
        mod_ini = mod_apr.copy()
        error_ini = mod_var
    else:
        error("Prior model method "+SetPrior.lower()+" not yet implemented! Exit.")

    Ctrl["name"] = fl_name


    start = process_time()
    """
    Loop over sites
    """
    sequence = range(nsite)
    if "rev" in Direction.lower():
        sites = sequence[::-1]
    else:
        sites = sequence


    logsize = (2 + 7*Maxiter)
    site_log = numpy.full((len(sites),logsize), numpy.nan)

    for ii in sites:
        print("\n Invert site #"+str(ii))

        """
        Setup model-related paramter
        """


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
        dat_act = numpy.tile(data_active,(nsite,1))
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
        Results = inverse.run_tikh_opt(Ctrl=Ctrl, Model=Model, Data=Data,
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
                jd = Results["jacd"]
                site_jacd = jd.reshape((1,numpy.size(jd)))
                postcov = Results["cpost"]
                site_postcov = postcov.reshape((1,numpy.size(postcov)))
                mres = Results["mresm"][0]
                site_nump = numpy.sum(numpy.diag(mres))

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


    header=numpy.array(Header, dtype=object),
    fileout = OutResDir + name + OutStrng + OutFileFmt

    numpy.savez_compressed(
        file=fileout,
        fl_data=file,
        fl_name=fl_name,
        header=header,
        aem_system=AEM_system,
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


    print("\n\nResults stored to "+fileout)
    elapsed = (process_time() - start)
    print (" Used %7.4f sec for %6i sites" % (elapsed, ii+1))
    print (" Average %7.4f sec/site\n" % (elapsed/(ii+1)))

print("\n\nAll done!")
