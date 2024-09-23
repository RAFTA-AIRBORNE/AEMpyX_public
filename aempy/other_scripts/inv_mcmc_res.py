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
#       jupytext_version: 1.11.4
# ---
import sys
from sys import exit as error
from datetime import datetime
from time import process_time
import time
import os
import warnings

import scipy.linalg
import numpy


import matplotlib
import matplotlib.pyplot

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)



import prep
import util
import core1d
import inverse
import aesys
import viz

from version import versionstrg

warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")

version, _ = versionstrg()
now = datetime.now()

Strng = "AEM Metropolis-Hastings\n"+"AEMpyX Version "+version
print("\n\n"+Strng)
print("".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

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
AEM_system = "aem05"

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 50.
    DatErr_mult = 0.0
    alt = 60.
    data_active = numpy.ones(NN[2], dtype="int8")


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 100.
    DatErr_mult = 0.01
    alt = 90.
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'

StoreChain = True
StoreSum = True
OutName="SYNTH_"

Nsample = 800000
Nburnin = 200000
Accept_check = 10000
Accept_max = 60.0
Accept_min = 30.0
StepStart = 0.5
StepFactor = 0.6666
CorrLength = 10.0

# percentiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.] # linear
Percentiles = [2.3, 15.9, 50., 84.1,97.7]                   # 95/68



AEMPYX_DATA = os.environ["AEMPYX_DATA"]

# InDatDir = AEMPYX_DATA+"/MH_Test/"
InDatDir = AEMPYX_ROOT+"/aempy/data/mh/"
print("Data read from dir:  %s" % InDatDir)

InFileFmt = ".npz"
dat_files = ["SYNTH_AEM05_1Layer_Conductor.npz"]

ReadFilelist = False
if ReadFilelist:
    ListFile = ""

SearchFilelist = False
if SearchFilelist:
    SearchStrng = "*.asc"
    print("Search string: %s " % SearchStrng)

# dat_files = ["SYNTH_AEM053Layer_NewCase_All_Alt_Perturb_50.asc"]
"""
Output formats are ".npz"
"""
OutFileFmt = ".npz"
# OutDatDir = AEMPYX_DATA + "/Nearest/fwd_compare/models/"
OutDatDir = InDatDir+"/results/"
print("Results written to dir: %s " % OutDatDir)

if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


if ReadFilelist:

    dat_files = []
    with open(ListFile, "r") as file:
        for line in file:
            dat_files.append(line[:-1])

if SearchFilelist:
    dat_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InDatDir)

dat_files = sorted(dat_files)
ns = numpy.size(dat_files)


"""
Define Model
"""

nlyr = 12

m_active, prior_avg, prior_std, m_bounds, m_state = inverse.init_1dmod(nlyr)

dzstart = 8.    #5.
dzend = 16.

dz = numpy.logspace(numpy.log10(dzstart), numpy.log10(dzend), nlyr) + 0.5 * numpy.random.randn(nlyr)
zc = inverse.set_zcenters(dz)
zn = inverse.set_znodes(dz)


m_active[0 * nlyr : 1 * nlyr] = 1
sizepar = numpy.shape(m_active)  # sizepar = sizepar[:]
sizeact = numpy.shape(numpy.flatnonzero(m_active))
mpara = sizeact[0]

guess_rho = 100
prior_avg[0 * nlyr : 1 * nlyr] = guess_rho
prior_avg[6 * nlyr : 7 * nlyr - 1] = dz[0 : nlyr - 1]

guess_std  = 0.1
prior_std[0 * nlyr : 1 * nlyr] = guess_std
para_var =prior_std[0 * nlyr : 1 * nlyr]**2

mod_prior, m_state = inverse.transform_parameter(m_vec=prior_avg, m_trn=ParaTrans, mode="f")

m_upper = m_bounds[:,0]
m_upper[0 * nlyr : 1 * nlyr] = 6.0  # prior_avg[m_active==1] + 3*prior_std[m_active==1]
m_u = m_upper[m_active==1]
m_lower = m_bounds[:,1]
m_lower[0 * nlyr : 1 * nlyr] = -1.0  # prior_avg[m_active==1] - 3*prior_std[m_active==1]
m_l = m_lower[m_active==1]

xc = yc = numpy.zeros(numpy.size(zc))
Cm_prior, _, _ = inverse.covar(xc, yc, zc, covtype= "exp", L=[CorrLength, CorrLength, CorrLength],
                var=para_var, sparse=False, thresh=0.05,
                calc_inverse=False)
print(Cm_prior)
CholFac = scipy.lininverse.cholesky(Cm_prior)

print(" Parameter set for inverting: \n", m_active)
print(" Layer thicknesses: \n", dz)
print(" Layer interface depths: \n", zn)
print(" Initial homogeneous halfspace resistivity of %6.2f Ohm*m \n" % (guess_rho))
print(" Log10 Standard error of %6.2f \n " % (guess_std))

print(" Upper limits: \n", m_upper)
print(" Lower limits: \n", m_lower)
print(
    " Assuming exponential vertical covariance with correlation length %6.2f m \n"
    % (CorrLength)
)



fcount =0
for file in dat_files:

    start = time.time()

    fcount=fcount+1

    name, ext = os.path.splitext(file)
    filein = InDatDir+file
    print("\n Reading file " + filein)
    Data, Header = aesys.read_aempy(File=filein,
                                   System=AEM_system, OutInfo=False)
    Data = numpy.reshape(Data,(1, -1))
    Alt =  Data[:, 4]
    DatObs =  Data[:, 6:6+NN[2]]

    data_err, data_obs = inverse.set_errors(DatObs,
                                                DatErr_add, DatErr_mult)


    sizedat = numpy.shape(data_obs)
    data_std = data_err
    data_cov = data_std * data_std
    ndata = sizedat[0]




    print("\n  starting simulation \n")

    # Metropolis-Hastings algorithm with Nsample iteration

    m_initial = mod_prior.copy()
    fp = numpy.sum(numpy.power((m_initial - mod_prior) / prior_std, 2))
    fp_old = fp


    data_calc, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall,
                                          alt=Alt,
                                          m_vec = m_initial,
                                          m_trn=ParaTrans,
                                          m_state = m_state,
                                          d_trn=DataTrans,
                                          d_act = data_active,
                                          OutInfo=True)


    dc_old = data_calc
    ss_old = numpy.sum(numpy.power(abs(data_obs - data_calc) / data_std, 2))

    m_old = m_initial


    scount = 0
    step = StepStart
    numpara = numpy.sum(m_active==1)

    mchain = numpy.zeros((Nsample + 1, mpara + ndata + 1))
    mchain[scount, 0:mpara] = m_old[m_active == 1]

    mchain[scount, mpara : mpara + ndata] = dc_old
    mchain[scount, mpara + ndata] = ss_old


    reject = 0
    accept = 0
    for n in range(Nsample):

        #    print('\n',n)
        # draw random numbers from the prior model (prior_avg)
        m_g = numpy.random.randn(numpara)
        #    print(' m_g   ',m_g[0*nlyr:1*nlyr])
        m_c = CholFac@m_g
        m_s = m_old[m_active==1] + step*m_c
        # print(' m_s   ', m_s)

        if numpy.all(m_s < m_u) and numpy.all(m_s > m_l):
            m_sample = inverse.insert_mod(M=m_old,m=m_s, m_act=m_active)
            scount = scount + 1
            m_new = m_sample
            fp_new = numpy.sum(numpy.power((m_new - m_old) / prior_std, 2))

            data_calc, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall,
                                              alt=Alt,
                                              m_vec = m_new,
                                              m_trn=ParaTrans,
                                              m_state = m_state,
                                              d_trn=DataTrans,
                                              d_act = data_active,
                                              OutInfo=True)

            dc_new = data_calc
            ss_new = numpy.sum(numpy.power(abs(data_obs - data_calc) / data_std, 2))
            # fl_new=numpy.exp(-0.5*ss_new)
            # alpha = min(1, numpy.exp(-0.5 * (ss_new - ss_old) - 0.5 * (fp_new - fp_old)))
            alpha = min(1, numpy.exp(-0.5 * (ss_new - ss_old) ))
            # print (alpha)
            p = numpy.random.rand(1, 1)
            if p < alpha:
                # accept new candidate
                # print(n, " accepted")
                accept = accept + 1

                mchain[scount, 0:mpara] = m_new[m_active == 1]


                mchain[scount, mpara : mpara + ndata] = dc_new
                mchain[scount, mpara + ndata] = ss_new

                m_old = m_new
                ss_old = ss_new
                fp_old = fp_new
                dc_old = dc_new
            else:
                mchain[scount, 0:mpara] = m_old[m_active == 1]
                mchain[scount, mpara : mpara + ndata] = dc_old
                mchain[scount, mpara + ndata] = ss_old
        else:
            # print(n, " rejected")
            reject = reject + 1
            continue

        if numpy.mod(n, Accept_check) == 0 and not n == 0:
            accpp = 100 * accept / scount
            if accpp > Accept_max:
                step = step / StepFactor
            if accpp < Accept_min:
                step = step * StepFactor
            print(
                " percentage of samples accepted = %6.2f %% of %6i, with  %6i out of bounds, step set to %6.4f"
                % (accpp, scount, reject, step)
            )


    print(" final percentage of samples accepted = %6.2f " % (100 * accept / scount))
    # numpy.savez('Mh_chain',mchain)
    print(numpy.shape(mchain[Nburnin:scount, 1]))
    # numpy.savez('Mh_chain1',mchain=mchain)
    # numpy.savez_compressed('Mh_chain2',mchain=mchain[Nburnin:scount,:])

    n_prc = len(Percentiles)


    chaindat = mchain[Nburnin:scount, mpara : mpara + ndata]
    avg_dat = numpy.mean(chaindat,axis=0)
    std_dat = numpy.std(chaindat,axis=0)
    prc_dat = numpy.percentile(chaindat, Percentiles)

    chainres= mchain[Nburnin:scount, :-1]

    avg_res = numpy.zeros(nlyr)
    std_res = numpy.zeros(nlyr)
    prc_res = numpy.zeros((nlyr, n_prc))

    for ilyr in range(nlyr):
        chainpar = mchain[Nburnin:scount, ilyr]
        avg_res[ilyr] = numpy.exp(numpy.mean(chainpar))
        std_res[ilyr] = numpy.exp(numpy.std(chainpar))
        prc_res[ilyr, :] = numpy.exp(numpy.percentile(mchain[Nburnin:scount, ilyr],
                                            Percentiles))

        print(
            " Layer %3i with center at %4.1f m Resistivity:"
            % (ilyr, zc[ilyr])
        )
        quantz =prc_res[ilyr,:]
        print(
            "Percentiles: "+str(quantz)
        )

    filen, filext = os.path.splitext(file)

    if StoreChain:
        filout = OutDatDir+ filen + "_MHchain" + "_" + str(Nsample) + "_CorrLength" + str(CorrLength) + filext
        print("MCMC Results written to: "+filout)
        numpy.savez_compressed(
            filout,
            mchain=mchain[Nburnin:scount, :],
            nlyr=nlyr, lcorr=CorrLength,
            avg_dat=avg_dat, std_dat=std_dat, prc_dat=prc_dat,
            avg_res=avg_res, std_res=std_res, prc_res=prc_res
       )
    if StoreSum:
        filout = OutDatDir+ filen + "_MHSummary" + "_" + str(Nsample) + "_CorrLength" + str(CorrLength) + filext
        print("MCMC Results written to: "+filout)
        numpy.savez_compressed(
            filout,
            nlyr=nlyr, lcorr=CorrLength,
            avg_dat=avg_dat, std_dat=std_dat, prc_dat=prc_dat,
            avg_res=avg_res, std_res=std_res, prc_res=prc_res
       )




Plots = False
FilesOnly = False
PlotFormat = [".pdf", ".png", ".svg"]
PdfCatalog = True
if ".pdf" in PlotFormat:
    PdfCName = OutName+"_Catalog.pdf"
else:
    error(" No pdfs generated. No catalog possible!")
    PdfCatalog = False
if PdfCatalog:
    pdf_list = []

if Plots:

    # Determine graphical parameter.
    # print(plt.style.available)
    matplotlib.pyplot.style.use("seaborn-paper")
    matplotlib.rcParams["figure.dpi"] = 400
    matplotlib.rcParams["axes.linewidth"] = 0.5
    matplotlib.rcParams["savefig.facecolor"] = "none"
    Fontsize = 10
    Labelsize = Fontsize
    Titlesize = Fontsize+2
    Linewidth= 1
    Markersize = 4

    Colors = ["r", "g", "b", "m", "c", "y", "k","r", "g", "b", "m"]
    Lines  = ["-", "--", ":", "-.","--", ":", "-.","--", ":","-","-."]

    if FilesOnly:
        matplotlib.use("cairo")

    plotfile = filout = filen + "_MHchain" + "_" + str(Nsample) + "_CorrLength" + str(CorrLength)
    # viz.plot_mcmc_results(
    #     PlotFile = plotfile,
    #     PlotTitle = plotfile.replace("_"," "),
    #     PlotFormat = PlotFormat,
    #     DPI = 400,
    #     Depth = zn,
    #     Param = quantz,
    #     Labels=["95%", "68%", "median (50%)"],
    #     Colors=["r", "g", "b", "g","r"],
    #     Lines=["-", ":", ";"],
    #     Linewidth=[1.5],
    #     Markersize = 4,
    #     Fontsizes=[Fontsize, Labelsize, Titlesize],
    #     MLimits= [],
    #     SLimits= [],
    #     DLimits= [0., 120.],
    #     PlotStrng="",
    #     StrngPos=[0.05,0.05])
