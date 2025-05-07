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
import getpass
import functools
import inspect


from cycler import cycler

import numpy
import scipy.interpolate
import scipy.linalg

import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.axis
import matplotlib.backends.backend_pdf #  matplotlib.backends. backend_pdf.PdfPages


AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import util
import aesys
import inverse
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)
cm = 1/2.54

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=inspect.getfile(inspect.currentframe()), out=False)
print(titstrng+"\n\n")

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
AEM_system = "aem05"
# AEM_system = "genesis"
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


Plotlist = ["model", "sens", "jac", "cov", "cor", "respar", "resdat"]

# Plotlist = [ "respar"]
# Plotlist = [ "model"]

Plotlist = [item.lower() for item in Plotlist]


if "model" in Plotlist:
    print("Model will be plotted")
    err = "lsq"
    err = "msq"
    err = "jsq"
    Modelcolor = ["b", "r", "r", ]
    Modellines = ["-", ":", ":" ]
    Modelwidth = [ 1,  1,  1,]
    ModelLimits = [1., 10000.]
    DepthLimits = [0., 100.]

if "sens" in Plotlist:
    print("Sensitivities will be plotted")
    whichsens = ["raw","cov", "euc" , "cum"]
    print("   Sensitivity type is ", str(whichsens))

    Senscolor = ["b", "g", "r", "m", "y"]
    Senslines =  ["-", "-", "-", "-", "-"]
    Senswidth = [ 1, 1,  1, 1, 1.]
    SensLimits = [0.001, 2.]
    DepthLimits = [0., 100.]

if "respar" in Plotlist:
    print("Parameter resolution will be plotted")
    whichspread = "fro"   #, too"", "euc", "mic"
    print("Spread type is "+whichspread)
    PhysAxes = True
    NoHalfspace = True

if "resdat" in Plotlist:
    print("Data resolution will be plotted")
    whichspread = "fro"   #, too"", "euc", "mic"
    print("   Spread type is "+whichspread)
    PhysAxes = True

if "cov" in Plotlist:
    print("Posterior covariance matrix will be plotted")
    NoHalfspace = True

if "cot" in Plotlist:
    print("Parameter correlation matrix will be plotted")
    NoHalfspace = True

if "jac" in Plotlist:
    print("Jacobian matrix will be plotted")
    NoHalfspace = True


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
input format is "npz"
"""
# GENESIS
# InModDir =  "/home/vrath/DuyguPoster/TD_uncert/"
# SearchStrng ="*901*results.npz"
# PDFCName = "GENESIS_FL901_Uncert-Catalog.pdf"
# AEM05
# InModDir =  AEMPYX_DATA + "/Projects/InvParTest/proc_delete_PLM3s/results_diffop/"
# InModDir =  "/home/vrath/DuyguPoster/FD_uncert/"
# SearchStrng ="A1*36*results.npz"

AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/A1_StGormans/"
InModDir =  AEMPYX_DATA+"/results_parallel/"
FileList = "search"
SearchStrng ="*results.npz"
if not InModDir.endswith("/"): InModDir=InModDir+"/"
print("Models read from dir: %s " % InModDir)
# FileList = "set" #"search"

FileList ="search"
if "search" in FileList.lower():
    how = ["search", SearchStrng, InModDir]
    # how = ["read", FileList, InModDir]
    print("Method is ", how )
    mod_files = util.get_data_list(how, out= True, sort=True)

# FileList ="set"
# mod_files = ["A1*results.npz"]

ns = numpy.size(mod_files)
if ns ==0:
    error("No modfiles set!. Exit.")

print("Filelist:")
print(mod_files)


"""
Plot formats are "".png", ".pdf", or any other
format matplotlib allows.
"""
PlotFmt = [".pdf", ".png"] #".png", ".pdf",]
# PlotDir = AEMPYX_DATA+"/ClaraUncert/plots/"
# PlotDir = InModDir+"/Lvar1/"
PlotDir = InModDir+"/plots/"
print("Plots written to dir: %s " % PlotDir)

if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)

PDFCatalog = True
PDFCatName = "AEM05_F11379_Uncert-Catalog.pdf"
if ".pdf" in PlotFmt:
    pass
else:
    error(" No pdfs generated. No catalog possible!")
    PDFCatalog = False


"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""
FilesOnly = False
matplotlib.pyplot.style.use("seaborn-v0_8-paper")
matplotlib.rcParams["figure.dpi"] = 600
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
matplotlib.rcParams["savefig.transparent"] = True
matplotlib.rcParams["savefig.bbox"] = "tight"
Fontsize = 6
Labelsize = Fontsize
Titlesize = 6
Fontsizes = [Fontsize, Labelsize, Titlesize]


# Markersize = 4
FigSize = [8.5*cm, 8.5*cm]

"""
https://matplotlib.org/stable/tutorials/colors/colormaps.html
"""
ColorMapResMat="seismic"
ColorMapCovMat="seismic"
ColorMapCorMat="seismic"
ColorMapJacMat="jet"

Grey20 = (0.2, 0.2, 0.2)
Grey50 = (0.5, 0.5, 0.5)
# Lines = (cycler("linewidth", [1.])
#          * cycler("linestyle", ["-", "--", ":", "-."])
#          * cycler("color", ["r", "g", "b", "m"]))


if FilesOnly:
    matplotlib.use("cairo")

if PDFCatalog:
    pdf_list = []
    catalog =matplotlib.backends.backend_pdf.PdfPages(PDFCatName)


for filein in mod_files:
    start = process_time()

    modfile = InModDir + filein
    ctrfile = modfile.replace("_results.npz", "_ctrl.npz")
    fnam, ext = os.path.splitext(os.path.basename(modfile))

    print("\nResults read from: %s" % modfile)
    results = numpy.load(modfile, allow_pickle=True)

    print("\nCtrl read from: %s" % modfile)
    control = numpy.load(ctrfile, allow_pickle=True)
    Runtyp = control["inversion"][0]
    Regfun = control["inversion"][1]
    OptnsStrng = "Opts: "+Runtyp+"|"+Regfun

    fline = results["fl_name"]
    site_x = results["site_x"]
    site_y = results["site_y"]
    site_z = results["site_dem"]

    m_act = results["mod_act"]
    m_ref = results["mod_ref"]


    site_mod = results["site_modl"]
    site_err = results["site_merr"]
    site_sns = results["site_sens"]
    site_rms = results["site_nrms"]

    site_dact = results["dat_act"]
    site_dobs = results["site_dobs"]
    site_dcal = results["site_dcal"]
    site_derr = results["site_derr"]
    site_rms = results["site_nrms"]

    site_jac= results["site_jacd"]
    site_cov= results["site_pcov"]

    nlyr = inverse.get_nlyr(m_ref)
    dz = m_ref[6*nlyr:7*nlyr-1]

    zn = inverse.set_znodes(dz)
    zm = inverse.set_zcenters(dz)
    DepthN = zn
    DepthC = numpy.append(zm, 999.9)
    LayThk = numpy.append(dz, 9999.)

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
        print ("ISITE="+str(isite))

        # calculation
        # generalized inverse
        npar = numpy.sum(m_act)
        ndat = numpy.sum(site_dact[isite,:])
        cov = site_cov[isite,:].reshape((npar,npar))
        jac = site_jac[isite,:].reshape((ndat,npar))
        # print(cov.shape)
        # print(jac.shape)
        dcal = site_dcal[isite,:]
        dact = site_dact[isite,:]
        dcal  = inverse.extract_dat(D=dcal, d_act=dact)
        scal = numpy.diag(1./dcal)

        v = numpy.sqrt(1./numpy.diag(cov))
        cor = cov*numpy.outer(v,v)

        # sensitivities
        sens = []

        sens0 = inverse.calc_sensitivity(Jac=jac, UseSigma=True, Type = "raw") #[:-1]
        sens0 = inverse.transform_sensitivity(S=sens0, V=LayThk,
                                              Transform=[" val","max"])
                                              # Transform=[" val","max", "sqr"])
        if NoHalfspace:
            sens1 = sens0[:-1]
        sens.append(numpy.abs(sens0))

        sens1 = inverse.calc_sensitivity(Jac=jac, UseSigma=True, Type = "cov") #[:-1]
        sens1 = inverse.transform_sensitivity(S=sens1, V=LayThk,
                                              Transform=["max"])
                                              # Transform=[" val","max", "sqr"])
        if NoHalfspace:
            sens1 = sens1[:-1]
        sens.append(numpy.abs(sens1))

        sens2 = inverse.calc_sensitivity(Jac=jac, UseSigma=True, Type = "euc") #[:-1]
        sens2 = inverse.transform_sensitivity(S=sens2, V=LayThk,
                                              Transform=[" max", "sqr"])
                                              # Transform=[" val","max", "sqr"])
        if NoHalfspace:
            sens2 = sens2[:-1]
        sens.append(numpy.abs(sens2))

        sens3 = inverse.calc_sensitivity(Jac=scal@jac, UseSigma=True, Type = "cum")
        sens3 = inverse.transform_sensitivity(S=sens3, V=LayThk,
                                             Transform=["max"])
                                             # Transform=[" val","max", "sqr"])
        if NoHalfspace:
            sens3 = sens3[:-1]
        sens.append(numpy.abs(sens3))
        sens.pop(0)


        # parameter resolution matrix & spread(s)

        gi =  cov@jac.T

        rm = gi@jac
        nm = numpy.sum(rm.diagonal())

        _, mspread0 = inverse.calc_model_resolution(J=jac, G=gi,
                                                    Spread=["frob"])
        _, mspread1 = inverse.calc_model_resolution(J=jac, G=gi,
                                                    Spread=["toomey"])
        _, mspread2 = inverse.calc_model_resolution(J=jac, G=gi,
                                                    Spread=["miller"])


        rd =  jac@gi
        nd = numpy.sum(rd.diagonal())


        FlineStrng = "FL: "+str(fline)


        if "sens" in Plotlist:



            if "dist" in Sample:
                TitleStrng = FlineStrng+", site "+str(numpy.rint(site_r[isite]))+" - sens"
            else:
                TitleStrng = FlineStrng+", site "+str(isite)+" - sens"

            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_sens"
            #PlotFile = PlotDir+fnam+"_site"+str(isite)+"_sens_unscaled"

            fig, _ =viz.plot_depth_prof(
                    PlotFile = PlotFile,
                    PlotTitle = TitleStrng+"\n"+OptnsStrng,
                    PlotFormat = PlotFmt,
                    FigSize = FigSize,
                    Depth = [zn],
                    Partyp = "sens",
                    Params = [sens],
                    DLabel = "depth (m)",
                    PLabel = "sensitivity (-)",
                    Legend = ["coverage", "euclidean","cumulative"],    #  "cummulative"
                    XScale = "log",
                    PlotType = "steps",
                    Linecolor=Senscolor,
                    Linetypes=Senslines,
                    Linewidth=Senswidth,
                    Fillcolor = [Grey50],
                    Fontsizes = Fontsizes,
                    PLimits = SensLimits,
                    DLimits = DepthLimits,
                    PlotStrng="", #Formula, #"", #"Error: mult="+str(DatErr_mult)+" add="+str(DatErr_add),
                    StrngPos=[0.05,0.05])


            if PDFCatalog:
                pdf_list.append(PlotDir+PlotFile+".pdf")
                catalog.savefig(fig)

        if "model" in Plotlist:

            if "dist" in Sample:
                TitleStrng = FlineStrng+", model, site "+str(numpy.rint(site_r[isite]))+" m "
            else:
                TitleStrng = FlineStrng+ "model, site "+str(isite)


            model = site_mod[isite, :]
            error = site_err[isite, :]
            print(error)
            val = numpy.log(model)
            errm = numpy.exp(val-error)
            errp = numpy.exp(val+error)
            model = [model, errm, errp]


            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_model"
            fig, _ = viz.plot_depth_prof(
                    PlotFile = PlotFile,
                    PlotTitle = TitleStrng+"\n"+OptnsStrng,
                    PlotFormat = PlotFmt,
                    FigSize=FigSize,
                    XScale = "log",
                    PlotType = "steps filled",
                    Depth = [zn],
                    Params = [model],
                    Partyp = "model",
                    DLabel = "depth (m)",
                    PLabel = "resistivity (Ohm m)",
                    Legend = [],
                    Linecolor=Modelcolor,
                    Linetypes=Modellines,
                    Linewidth=Modelwidth,
                    Fillcolor = [Grey50],
                    Fontsizes=Fontsizes,
                    PLimits = ModelLimits,
                    DLimits = DepthLimits,
                    PlotStrng="nRMS = "+str(numpy.around(site_rms[isite][0],2)),
                    StrngPos=[0.05,0.05])

            if PDFCatalog:
                 pdf_list.append(PlotFile+".pdf")
                 catalog.savefig(fig)

        if "jac" in Plotlist:


            if "dist" in Sample:
                TitleStrng = FlineStrng+", jac, site "+str(numpy.rint(site_r[isite]))+" m"
            else:
                TitleStrng = FlineStrng+", jac, site "+str(isite)

            Matrix = jac.T
            if NoHalfspace:
                Matrix  = Matrix[:-1]

            yticks = numpy.arange(nlyr)[0:-1:5]

            if PhysAxes:
                yticklabels = numpy.rint(DepthC[yticks]).astype(int).astype(str)
                ylabel = "depth (m)"
            else:
                yticklabels = yticks.astype(str)
                ylabel =  "layer #"


            if "aem05" in AEM_system:
                xlabel = "data #"
                xticks = numpy.arange(8)
                xticklabels = xticks.astype(str)

                if PhysAxes:
                    xlabel = "frequency (kHz)"
                    pars = Pars[0]*1.e-3
                    vals = numpy.concatenate((pars, pars))
                    iticks = numpy.array([0, 2, 4 , 6])
                    xticks = xticks[iticks]
                    xticklabels = numpy.round(vals,2).astype(str)[iticks]



            if "genes" in AEM_system:
                xlabel =  "data #"
                xticks = numpy.arange(11)
                xticklabels = xticks.astype(str)


                if PhysAxes:
                    xlabel = "window center (1e-6 s)"
                    vals = Pars[0]*1000.
                    iticks = numpy.array([0, 2, 4, 6, 8, 10])
                    xticks = xticks[iticks]
                    xticklabels = numpy.round(vals,1).astype(str)[iticks]



            Aspect = "auto" #aspect
            AxLabels = [xlabel, ylabel]
            AxTicks  = [xticks, yticks]
            AxTickLabels =  [xticklabels, yticklabels]
            TickStr=["", ""]

            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_jac"
            viz.plot_matrix(
                  PlotFile=PlotFile,
                  PlotTitle=TitleStrng+"\n"+OptnsStrng,
                  PlotFormat=PlotFmt,
                  Matrix=Matrix,
                  FigSize=FigSize,
                  ColorMap=ColorMapJacMat,
                  TickStr=TickStr,
                  AxLabels=AxLabels,
                  AxTicks=AxTicks ,
                  AxTickLabels=AxTickLabels,
                  Fontsizes=Fontsizes,
                  Aspect =Aspect,
                  PlotStrng="",
                  StrngPos=[0.05,0.05])

            if PDFCatalog:
                  pdf_list.append(PlotFile+".pdf")
                  catalog.savefig(fig)

        if "cov" in Plotlist:

            if "dist" in Sample:
                TitleStrng = FlineStrng+", p-cov, site "+str(numpy.rint(site_r[isite]))+" m"
            else:
                TitleStrng = FlineStrng+", p-cov, site "+str(isite)

            Matrix = cov

            xticks = numpy.arange(nlyr)[0:-1:5]
            xticklabels = xticks.astype(str)
            yticks = xticks
            if PhysAxes:
                yticklabels = numpy.rint(DepthC[yticks]).astype(int).astype(str)
                AxLabels = [" layer #"," depth (m)"]
            else:
                yticklabels = yticks.astype(str)
                AxLabels =  ["layer #","layer #"]
            Aspect = "equal"
            AxTicks = [xticks, yticks]
            AxTickLabels =  [xticklabels, yticklabels]
            TickStr=["", ""]

            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_pcov"
            viz.plot_matrix(
                 PlotFile=PlotFile,
                 PlotTitle=TitleStrng+"\n"+OptnsStrng,
                 PlotFormat=PlotFmt,
                 Matrix=Matrix,
                 FigSize=FigSize,
                 ColorMap=ColorMapCovMat,
                 TickStr=TickStr,
                 AxLabels=AxLabels,
                 AxTicks=AxTicks ,
                 AxTickLabels=AxTickLabels,
                 Aspect =Aspect,
                 Fontsizes=Fontsizes,
                 PlotStrng="",
                 StrngPos=[0.05,0.05])

            if PDFCatalog:
                 pdf_list.append(PlotFile+".pdf")
                 catalog.savefig(fig)

        if "cor" in Plotlist:

            if "dist" in Sample:
                TitleStrng = FlineStrng+", p-cor, site "+str(numpy.rint(site_r[isite]))+" m"
            else:
                TitleStrng = FlineStrng+", p-cor, site "+str(isite)

            Matrix = cor

            xticks = numpy.arange(nlyr)[0:-1:5]
            xticklabels = xticks.astype(str)
            yticks = xticks
            if PhysAxes:
                yticklabels = numpy.rint(DepthC[yticks]).astype(int).astype(str)
                AxLabels = [" layer #"," depth (m)"]
            else:
                yticklabels = yticks.astype(str)
                AxLabels =  ["layer #","layer #"]

            AxTicks = [xticks, yticks]
            AxTickLabels =  [xticklabels, yticklabels]
            TickStr=["", ""]

            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_pcor"
            viz.plot_matrix(
                 PlotFile=PlotFile,
                 PlotTitle=TitleStrng+"\n"+OptnsStrng,
                 PlotFormat=PlotFmt,
                 Matrix=Matrix,
                 FigSize=FigSize,
                 ColorMap=ColorMapCorMat,
                 TickStr=TickStr,
                 AxLabels=AxLabels,
                 AxTicks=AxTicks ,
                 AxTickLabels=AxTickLabels,
                 Fontsizes=Fontsizes,
                 Aspect =Aspect,
                 PlotStrng="",
                 StrngPos=[0.05,0.05])

            if PDFCatalog:
                 pdf_list.append(PlotFile+".pdf")
                 catalog.savefig(fig)

        if "respar" in Plotlist:

            if "dist" in Sample:
                TitleStrng = FlineStrng+", p-res, site "+str(numpy.rint(site_r[isite]))+" m"
            else:
                TitleStrng = FlineStrng+", p-res, site "+str(isite)

            Matrix = rm
            if NoHalfspace:
                Matrix = Matrix[:-1,:-1]

            Np = numpy.sum(numpy.diag(rm))
            PlotStrng = " Npar = "+numpy.around(Np,1).astype(str)
            StrngPos=[0.05,0.1]

            xticks = numpy.arange(nlyr)[0:-1:5]
            xticklabels = xticks.astype(str)
            yticks = xticks
            if PhysAxes:
                yticklabels = numpy.rint(DepthC[yticks]).astype(int).astype(str)
                AxLabels = [" layer #"," depth (m)"]
            else:
                yticklabels = yticks.astype(str)
                AxLabels =  ["layer #","layer #"]


            Aspect = "equal"
            AxTicks = [xticks, yticks]
            AxTickLabels =  [xticklabels, yticklabels]
            TickStr=["", ""]

            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_pres"
            viz.plot_matrix(
                PlotFile=PlotFile,
                PlotTitle=TitleStrng+"\n"+OptnsStrng,
                PlotFormat=PlotFmt,
                Matrix=Matrix,
                FigSize=FigSize,
                ColorMap=ColorMapResMat,
                AxLabels=AxLabels,
                AxTicks=AxTicks,
                AxTickLabels=AxTickLabels,
                TickStr=TickStr,
                Fontsizes=Fontsizes,
                Aspect =Aspect,
                PlotStrng=PlotStrng,
                StrngPos=StrngPos)

            if PDFCatalog:
                pdf_list.append(PlotFile+".pdf")
                catalog.savefig(fig)

        if "resdat" in Plotlist:

            if "dist" in Sample:
                TitleStrng = FlineStrng+", d-res, site "+str(numpy.rint(site_r[isite]))+" m"
            else:
                TitleStrng = FlineStrng+", d-res, site "+str(isite)

            Matrix = rd

            Nd = numpy.sum(numpy.diag(rd))
            PlotStrng = " Ndat = "+numpy.around(Nd,1).astype(str)
            StrngPos=[0.05,0.1]

            if "aem05" in AEM_system:
                Axlabels =  ["data #","data #"]
                xticks = numpy.arange(8)
                xticklabels = xticks.astype(str)
                yticks =  xticks
                yticklabels = xticklabels

                if PhysAxes:
                    AxLabels =  ["data #"," frequency (kHz)"]
                    pars = Pars[0]*1.e-3
                    vals = numpy.concatenate((pars, pars))
                    iticks = numpy.array([0, 2, 4 , 6])
                    yticks = yticks[iticks]
                    yticklabels = numpy.round(vals,2).astype(str)[iticks]



            if "genes" in AEM_system:
                AxLabels =  ["data #","data #"]
                xticks = numpy.arange(11)
                xticklabels = xticks.astype(str)

                yticks = xticks
                yticklabels = xticklabels

                if PhysAxes:
                    AxLabels =  ["data #"," window center (1e-6 s)"]
                    vals = Pars[0]*1000.
                    iticks = numpy.array([0, 2, 4, 6, 8, 10])
                    yticks = yticks[iticks]
                    yticklabels = numpy.round(vals,1).astype(str)[iticks]

            Aspect = "equal"
            AxTicks = [xticks, yticks]
            AxTickLabels =  [xticklabels, yticklabels]
            TickStr=["", ""]

            PlotFile = PlotDir+fnam+"_site"+str(isite)+"_dres"
            viz.plot_matrix(
                PlotFile=PlotFile,
                PlotTitle=TitleStrng+"\n"+OptnsStrng,
                PlotFormat=PlotFmt,
                Matrix=Matrix,
                FigSize=FigSize,
                ColorMap=ColorMapResMat,
                AxLabels=AxLabels,
                AxTicks=AxTicks,
                AxTickLabels=AxTickLabels,
                TickStr=TickStr,
                Fontsizes=Fontsizes,
                Aspect =Aspect,
                PlotStrng=PlotStrng,
                StrngPos=StrngPos)


            if PDFCatalog:
                pdf_list.append(PlotDir+PlotFile+".pdf")
                catalog.savefig(fig)






if PDFCatalog:
    print(pdf_list)
    # viz.make_pdf_catalog(PDFList=pdf_list, FileName=PDFCatName)
    # print(str(len(pdf_list))+" collected to "+PlotDir+PDFCatName)
    d = catalog.infodict()
    d["Title"] =  PDFCatName
    d["Author"] = getpass.getuser()
    d["CreationDate"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    catalog.close()
