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
#       format_version: "1.5"
#       jupytext_version: 1.11.4
# ---

"""
Created on Tue Aug  3 17:03:39 2021

@author: vrath
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings

from cycler import cycler
import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.axis
import mpl_toolkits.axes_grid1

import scipy.interpolate


AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys
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
# AEM_system = "genesis"
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _,Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.03
    data_active = numpy.ones(NN[2], dtype="int8")
    CompDict=Misc[2]
    CompLabl = list(CompDict.keys())
    print(CompLabl)

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + "good" hoizontals"
    CompDict = Misc[2]
    CompLabl = list(CompDict.keys())







"""
input formats are "npz","nc4","ascii"
"""
InFilFmt = ".npz"
InDatDir = AEMPYX_DATA+"/Projects/StGormans/raw/"
# InDatDir = AEMPYX_DATA+"/Blocks/A1O/proc_delete/data_asc/"
print("Data read from dir: %s " % InDatDir)
FileList = "set" #"search"
# FileList = "search"
SearchStrng = "*.npz"
# print("Searchstring: %s \n" % SearchStrng)

"""
Output formats are "npz","nc4","ascii"
"""
PlotFmt = [".pdf"] #".png", ".pdf",]

PdfCatalog = True
PdfCName = "StGormans_Catalog_raw.pdf"
if ".pdf" in PlotFmt:
    pass
else:
    error(" No pdfs generated. No catalog possible!")
    PdfCatalog = False

PlotDir = InDatDir+"/plots/"
print("Plots written to dir: %s " % PlotDir)
PlotName = "StGormans_raw.pdf"
print("Plot filname: %s " % PlotName)

if "set" in FileList.lower():
    InDatDir =""
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = [AEMPYX_DATA+"/Projects/StGormans/A1_rect_StGormans_Full.npz"]

else:
    print("da")
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(dat_files)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""
FilesOnly = False
matplotlib.pyplot.style.use("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
Fontsize = 8
Labelsize = Fontsize
Titlesize = 8
Fontsizes = [Fontsize, Labelsize, Titlesize]
Linewidths= [0.6]
Markersize = 4

ncols = 11
Colors = matplotlib.pyplot.cm.jet(numpy.linspace(0,1,ncols))
Grey = 0.7
Lcycle =Lcycle = (cycler("linestyle", ["-", "--", ":", "-."])
          * cycler("color", ["r", "g", "b", "y"]))



if FilesOnly:
    matplotlib.use("cairo")


if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)


pdf_list = []
for filein in dat_files:
    start = process_time()
    file = InDatDir + filein
    print("\nData read from: %s" % file)
    Data, header, _ = aesys.read_aempy(File=file, System=AEM_system, OutInfo=False)



    pdf_list = []



    fig, ax = matplotlib.pyplot.subplots()
    fig.set_figwidth(FigWidth)


    if ScatterPlot:
        print("Scatter Plot")


        im = matplotlib.pyplot.scatter(E, N, c=D, s=Markersize**2, cmap=cmp)

        # ax = matplotlib.pyplot.gca()
        ax.set_aspect("equal")
        ax.xaxis.set_major_formatter(xformatter)
        ax.set_xlabel("Easting "+XYUnits, size=Fontsizes[1])
        ax.yaxis.set_major_formatter(yformatter)
        ax.set_ylabel("Northing "+XYUnits, size=Fontsizes[1])

        ax.tick_params(axis="x", labelsize=Fontsizes[1]-2, labelrotation=0.)#-45)
        ax.tick_params(axis="y", labelsize=Fontsizes[1]-2, labelrotation=0.)#-45)
        ax.grid(which="major", axis="both", visible=True,linewidth= Linewidths[0],linestyle=" --")
        ax.set_title(AEM_system.upper()+": "+ Comp)

        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = matplotlib.pyplot.colorbar(im, cax=cax)
        cb.ax.set_title(Unit)

    if ContourPlot or ImagePlot:
        print("Contour/Image Plot")

        xi= numpy.linspace(E_min,E_max,numIndexes)
        yi= numpy.linspace(N_min,N_max,numIndexes)
        XI, YI = numpy.meshgrid(xi, yi)
        Mesh = numpy.vstack((XI.ravel(), YI.ravel())).T
        Pnts = numpy.vstack((E.ravel(), N.ravel())).T
        Dats = D.ravel()

        DI = scipy.interpolate.griddata(Pnts, Dats, Mesh, method=Imethod)
        DI = numpy.reshape(DI,(len(xi), len(yi)))

        RBF = scipy.interpolate.RBFInterpolator(Pnts, Dats,
                                                kernel=Imethod, smoothing=0.)
        DI = RBF(mesh.flat)

        im= ax.pcolor(XI, YI, DI, cmap=cmp)

        ax.set_aspect("equal")
        ax.xaxis.set_major_formatter(xformatter)
        ax.set_xlabel("Easting "+XYUnits, size=Fontsizes[1])
        ax.yaxis.set_major_formatter(yformatter)
        ax.set_ylabel("Northing "+XYUnits, size=Fontsizes[1])

        ax.set_title(titl,fontsize=Fontsize)

        ax.grid(color="k", alpha=0.5, linestyle="dotted", linewidth=1.5)
        ax.tick_params(labelsize=Labelsize)

        cb = matplotlib.pyplot.colorbar(im)
        # cb.set_label(Unit) #, rotation=90.)# 270)
        cb.ax.set_title(Unit)

    plotfile = PlotDir+PlotName+"_"+AEM_system+"_"+ Comp


    for F in PlotFmt:
        print("Plot written to "+plotfile+F)
        matplotlib.pyplot.savefig(plotfile+F, dpi=600)

    if PdfCatalog:
        pdf_list.append(plotfile+".pdf")

    matplotlib.pyplot.show()
    matplotlib.pyplot.clf()

if PdfCatalog:
    viz.make_pdf_catalog(PdfList=pdf_list, FileName=PdfCName)
