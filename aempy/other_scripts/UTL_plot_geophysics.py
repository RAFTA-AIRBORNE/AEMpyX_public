#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---

"""
@author: vr May 2022
"""

# Import required modules

import os
import sys
from sys import exit as error
import csv
import warnings
from time import process_time
from datetime import datetime
import simplekml
import numpy
from cycler import cycler
import matplotlib
import matplotlib.pyplot

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)


import aesys
import util
import viz

from version import versionstrg


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]


version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Generate KML for boreholes & geophysics "+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


# Define the path to your data-files
SiteDir = AEMPYX_DATA+"/Intersection/geophysics/"
print(" site files read from: %s" % SiteDir)
# site_files = ["Verified_Boreholes_utm_NM_A1.npz",
#                 "Verified_Boreholes_utm_NM_A2.npz",
#                 "Verified_Boreholes_utm_NM_TB.npz"]

FlightLines = True

MinProfLength = 200.
site_files = ["Geophysics_2018_utm_overlaps.npz", "Geophysics_2020_utm_overlaps.npz"]

LinesMin = 10
SearchRadius=300.
DataDir = SiteDir+"/data/"
InpFmt = ".asc"
OutFmt = InpFmt

print(" data files read from: %s" % DataDir)
# Define the path for plots
PlotDir= DataDir

FilesOnly = False

PlotFmt = [".pdf",".png",]
PdfCatalog = True
if not ".pdf" in PlotFmt:
    error(" No pdfs generated. No catalog possible!")
    PdfCatalog = False


DataTrans = "asinh"
LogPlot = False
LogSym = False
LinThresh =100.
if LogPlot == False:
    LogSym = False
Logparams=[LogPlot, LogSym, LinThresh]
HLimits = [50., 200.]
PLimits = [0., 25.]

# ProfScale = 0.001  # m to km
ProfScale = 1. # 0.001  # m to km
ProfUnit  = "m" #

PlotThresh =20
PosDegrees = False
if PosDegrees:
    EPSG=32629
"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""

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

"""
For just plot.ing to files, choose the cairo backend (eps, pdf, png, jpg...).
If you need to see the plots directly (plot window, or jupyter), simply
comment out the following line. In this case matplotlib may run into
memory problems after a few hundreds of high-resolution plot..
"""
if FilesOnly:
    matplotlib.use("cairo")


# Determine which geographical info is added to the KML-tags:
# define empty list
sites = []
for f in site_files:
    sdata = numpy.load(SiteDir+f, allow_pickle=True)["Data"]
    print("reading "+SiteDir+f)
    pname=sdata[:,0]

    slat, slon = util.project_utm_to_latlon(sdata[:,1], sdata[:,2])
    num_site =-1
    for p in pname:
        num_site =num_site+1
        this_site = "/Geophysics_"+p+"/"

        if float(sdata[num_site,5])<= MinProfLength:
            print("\n\n\n*****"+str(sdata[num_site,3])
                  +"("+str(sdata[num_site,3])+") layout, "+str(sdata[num_site,5])
                  +", for "+p+" is < "+str(MinProfLength)+"m\n\n\n")
            continue

        indir = DataDir+"/Geophysics_"+p+"/"
        if not os.path.isdir(indir):
            error(" File: %s does not exist! Exit." % indir)

        plot_dir = PlotDir+this_site
        plot_file = "Geophysics_"+p
        if PdfCatalog:
            pdf_list = []


        AEM_system = "aem05"
        FwdCall, nD, _, _, _, = aesys.get_system_params(System=AEM_system)
        fd_file = AEM_system.upper()+"_Geophysics"+p+"_Datafile"+"_SearchRadius"+str(round(SearchRadius))+"m"
        print(fd_file)
        fdata, head = aesys.read_aempy(File=indir+fd_file+InpFmt, System=AEM_system, OutInfo=False)
        fs = numpy.shape(fdata)

        flat, flon = util.project_utm_to_latlon(fdata[:,1], fdata[:,2])

        fdata[:, 1] = fdata[:, 1] * ProfScale
        fdata[:, 2] = fdata[:, 2] * ProfScale

        flines = numpy.unique(fdata[:,0])

        for fline in flines:

            tmp = fdata[numpy.where(fdata[:, 0] == fline), :]
            ns = numpy.shape(tmp)
            tmp = numpy.reshape(tmp, (ns[1], ns[2]))
            print("OutInfo: "+str(numpy.shape(tmp)))

            if numpy.size(tmp)<=nD[0]*LinesMin:
                print("Not enough data! Not written")
                continue

            if FlightLines:
                flname=AEM_system.upper()+"_FL"+str(fline)
                aesys.write_aempy(File=indir+flname+OutFmt,
                                  Data=tmp, System=AEM_system,
                                  Header=head, OutInfo=False)

            if PdfCatalog:
                pdf_list.append(plot_dir+plot_file+"_"+str(fline)+".pdf")

            viz.plot_flightline_aem05(
                PlotName = plot_file+"_"+str(fline),
                PlotDir = plot_dir,
                PlotFmt=PlotFmt,
                DataObs=tmp,
                DataCal=[],
                XLimits =[],
                YLimits =[],
                HLimits = HLimits,
                PLimits = PLimits,
                ProfUnit=ProfUnit,
                Colors=Colors,
                Linewidths=Linewidths,
                Fontsizes=Fontsizes,
                Logparams=Logparams,
                PlotStrng=" - raw",
                PlotPLM = True)



        AEM_system = "genesis"
        FwdCall, nD, _, _, _, = aesys.get_system_params(System=AEM_system)

        td_file = AEM_system.upper()+"_Geophysics"+p+"_Datafile"+"_SearchRadius"+str(round(SearchRadius))+"m"
        print(td_file)
        tdata, _ = aesys.read_aempy(File=indir+td_file+InpFmt, System=AEM_system, OutInfo=False)
        ts = numpy.shape(fdata)

        tlat, tlon = util.project_utm_to_latlon(tdata[:,1], tdata[:,2])

        tdata[:, 1] = tdata[:, 1] * ProfScale
        tdata[:, 2] = tdata[:, 2] * ProfScale

        tlines = numpy.unique(tdata[:,0])

        for tline in tlines:

            tmp = tdata[numpy.where(tdata[:, 0] == tline), :]
            ns = numpy.shape(tmp)
            tmp = numpy.reshape(tmp, (ns[1], ns[2]))
            print("OutInfo: "+str(numpy.shape(tmp)))

            if numpy.size(tmp)<=nD[0]*LinesMin:
                print("Not enough data! Not written")
                continue

            if FlightLines:
                flname=AEM_system.upper()+"_FL"+str(tline)
                aesys.write_aempy(File=indir+flname+OutFmt,
                                  Data=tmp, System=AEM_system,
                                  Header=head, OutInfo=False)

            if PdfCatalog:
                pdf_list.append(plot_dir+plot_file+"_"+str(tline)+".pdf")

            viz.plot_flightline_genesis(PlotName = plot_file+"_"+str(tline),
               PlotDir = plot_dir,
               PlotFmt=PlotFmt,
               DataObs=tmp,
               DataCal=[],
               DataTrans = DataTrans,
               XLimits =[],
               YLimits =[],
               HLimits =[],
               ProfUnit=ProfUnit,
               Colors=Colors,
               Linewidths=Linewidths,
               Fontsizes=Fontsizes,
               Logparams=Logparams,
               PlotStrng=" - raw")
        tlines = numpy.unique(tdata[:,0])



        if PdfCatalog:
            viz.make_pdf_catalog(PdfList=pdf_list, FileName=indir+p+"_flines.pdf")
