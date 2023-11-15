#!/usr/bin/env python3

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
@author: vrath Feb 2021
@author: duygu June 2021
Lines 139-140 are corrected.
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

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)


from version import versionstrg
# import util
# import prep
import aesys
import inverse
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]


version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus Genesis  "
      +"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


"""
Define AEM system, and optioal additional parameters (e.g., Waveforms )
"""
AEM_system = "genesis"
_,NN, _, _, _, = aesys.get_system_params(AEM_system)
nD = NN[0]

if "gen" in AEM_system.lower():
    _,NN, _, _, _, = aesys.get_system_params(AEM_system)
    nD = NN[0]
    Transform = 1
    wct =  [
        r"0.009 ms", r"0.026 ms", r"0.052 ms", r"0.095 ms",
        r"0.156 ms", r"0.243 ms", r"0.365 ms", r"0.547 ms",
        r"0.833 ms", r"1.259 ms", r"1.858 ms"
          ]
    data_files = ["NM.asc"]


if "aem" in AEM_system.lower():
    _,NN, _, _, _, = aesys.get_system_params(AEM_system)
    nD = NN[0]
    Transform = 0
    freq = []

    data_files = ["A1.asc","A2.asc","A3.asc","A4.asc","A5.asc",
                  "A6.asc","A7.asc","A8.asc","A9.asc",
                  "TB.asc","WF.asc",]


"""
input formats are .npz, .nc4, "ascii"
"""
InFileFmt = "asc"
Value = 0.01
bins = "auto"
InDatDir = AEMPYX_DATA+"/Blocks/Area/"
OutDatDir = InDatDir

PlotIt = False
FilesOnly = True
PlotDir = InDatDir+"plots"

PlotFmt = [".pdf",".png",]
PdfCatalog = True
if not ".pdf" in PlotFmt:
    error(" No pdfs generated. No catalog possible!")
    PdfCatalog = False

PosDegrees = False
if PosDegrees:
    EPSG=32629

Transform = 0

# """
# Determine graphical parameter.
# => print(matplotlib.pyplot.style.available)
# """

# matplotlib.pyplot.style.use("seaborn-paper")
# matplotlib.rcParams["figure.dpi"] = 400
# matplotlib.rcParams["axes.linewidth"] = 0.5
# matplotlib.rcParams["savefig.facecolor"] = "none"
# Fontsize = 8
# Labelsize = Fontsize
# Titlesize = 8
# Fontsizes = [Fontsize, Labelsize, Titlesize]

# Linewidths= [0.6]
# Markersize = 4

# ncols = 11
# Colors = matplotlib.pyplot.cm.jet(numpy.linspace(0,1,ncols))
# Grey = 0.7
# Lcycle =Lcycle = (cycler("linestyle", ["-", "--", ":", "-."])
#           * cycler("color", ["r", "g", "b", "y"]))

# """
# For just plot.ing to files, choose the cairo backend (eps, pdf, png, jpg...).
# If you need to see the plots directly (plot window, or jupyter), simply
# comment out the following line. In this case matplotlib may run into
# memory problems after a few hundreds of high-resolution plot..
# """
# if FilesOnly:
#     matplotlib.use("cairo")

# if not os.path.isdir(PlotDir):
#     print("File: %s does not exist, but will be created" % PlotDir)
#     os.mkdir(PlotDir)

# if PdfCatalog:
#     pdf_list = []


ifl = -1
for file in data_files:
    ifl =ifl +1
    name, ext = os.path.splitext(file)
    PlotName =name
    filein =InDatDir + file
    print("data file: "+filein)
    start = process_time()
    Data, Header = aesys.read_aempy(File=filein, Format=InFileFmt,
                                   System=AEM_system, OutInfo=False)
    print("time used for reading = "+str(process_time()-start))

    if ifl == 0:
        DataAll = Data
    else:
        DataAll = numpy.vstack((DataAll,Data))

    # X = Data[:, 6:17].copy()
    # Z = Data[:, 17:28].copy()
    # nd =numpy.size(X)

    # Xmax = numpy.nanmax(X)
    # Xmin = numpy.nanmin(X)
    # #print(Xmax,Xmin)

    # Zmax = numpy.nanmax(Z)
    # Zmin = numpy.nanmin(Z)
    # #print(Zmax,Zmin)

    # ZData = numpy.zeros_like(Z)
    # XData = numpy.zeros_like(X)

    # Xlabel="X-component (fT)"
    # Zlabel="Z-component (fT)"
    # # print(Transform)

    # if Transform>0:

    #     if Transform==1:
    #         Xlabel="ln "+Xlabel
    #         Zlabel="ln "+Zlabel

    #     if Transform==2:
    #         Xlabel="arcsinh "+Xlabel
    #         Zlabel="arcsinh "+Zlabel
    #         print(Xlabel)

    #     SX = inverse.get_S(d=X.flat)
    #     SZ = inverse.get_S(d=Z.flat)
    #     S = max([SX, SZ])
    #     for ivec in numpy.arange(11):

    #         XData[:, ivec], _, d_state = inverse.transform_data(d_vec=X[:, ivec],
    #                                           d_trn=Transform, S=S)
    #         ZData[:, ivec], _, d_state = inverse.transform_data(d_vec=Z[:,ivec],
    #                                           d_trn=Transform, S=S)

    # else:

    #         XData = X
    #         ZData = Z

    # Xmax = numpy.nanmax(XData)
    # Xmin = numpy.nanmin(XData)
    # # print(Xmax,Xmin)

    # Zmax = numpy.nanmax(ZData)
    # Zmin = numpy.nanmin(ZData)
    # # print(Zmax,Zmin)


    # d = XData.flat

    # if PlotIt:
    #     matplotlib.pyplot.hist(d, bins, color = "r", density=True)
    #     matplotlib.pyplot.xlabel(Xlabel,
    #                                 fontsize = Labelsize)
    #     matplotlib.pyplot.ylabel("Probability",
    #                                 fontsize = Labelsize)
    #     matplotlib.pyplot.title("Histogram X",
    #                             fontsize = Titlesize)
    #     matplotlib.pyplot.grid(True)
    #     matplotlib.pyplot.xlim(Xmin, Xmax)

    #     P = PlotDir+PlotName+"_X"
    #     for F in PlotFmt:
    #         matplotlib.pyplot.savefig(P+F, dpi=400)
    #     if PdfCatalog:
    #         pdf_list.append(P+".pdf")

    #     matplotlib.pyplot.show()

    # d = ZData.flat
    # if PlotIt:
    #     matplotlib.pyplot.hist(d, bins, color = "g", density=True)
    #     matplotlib.pyplot.xlabel(Zlabel,
    #                                 fontsize = Labelsize)
    #     matplotlib.pyplot.ylabel("Probability",
    #                                 fontsize = Labelsize)
    #     matplotlib.pyplot.title("Histogram Z",
    #                             fontsize = Titlesize)
    #     matplotlib.pyplot.grid(True)
    #     matplotlib.pyplot.xlim(Zmin, Zmax)

    #     P = PlotDir+PlotName+"_Z"
    #     for F in PlotFmt:
    #         matplotlib.pyplot.savefig(P+F, dpi=400)
    #     if PdfCatalog:
    #         pdf_list.append(P+".pdf")
    #     matplotlib.pyplot.show()

    print ("\n\nTime Windows:")
    for icol in numpy.arange(11):
            d = XData[:,icol]
            ds = numpy.size(d)
            dcount = (d.flat<Value).sum()
            print ("number of values less than "+str(S*Value)+" in X data: "+str(dcount/ds)
                    +"  =  "+str(100*dcount/ds)+"% ("+wct[icol]+")")

    #         if PlotIt:
    #             matplotlib.pyplot.hist(d, bins, color = "r", density=True)
    #             matplotlib.pyplot.xlabel(Xlabel,
    #                                     fontsize = Labelsize)
    #             matplotlib.pyplot.ylabel("Probability",
    #                                     fontsize = Labelsize)
    #             matplotlib.pyplot.title("Histogram X (window center = "+wct[icol]+")",
    #                                     fontsize = Titlesize)
    #             matplotlib.pyplot.grid(True)
    #             matplotlib.pyplot.xlim(Xmin, Xmax)

    #             P = PlotDir+PlotName+"_X"+wct[icol].split()[0]
    #             for F in PlotFmt:
    #                 matplotlib.pyplot.savefig(P+F, dpi=400)
    #             if PdfCatalog:
    #                 pdf_list.append(P+".pdf")

    #             matplotlib.pyplot.show()

    # for icol in numpy.arange(11):
    #         d = ZData[:,icol]
    #         ds = numpy.size(d)
    #         dcount = (d.flat<Value).sum()
    #         print ("number of values less than "+str(S*Value)+" in Z data: "+str(dcount/ds)
    #                 +"  =  "+str(100*dcount/ds)+"% ("+wct[icol]+")")

    #         if PlotIt:
    #             matplotlib.pyplot.hist(d, bins, color = "g", density=True)
    #             matplotlib.pyplot.xlabel(Zlabel)
    #             matplotlib.pyplot.ylabel("Probability")
    #             matplotlib.pyplot.title("Histogram X (window center = "+wct[icol]+")",
    #                                     fontsize = Titlesize)
    #             matplotlib.pyplot.grid(True)
    #             matplotlib.pyplot.xlim(Zmin, Zmax)
    #             P = PlotDir+PlotName+"_Z"+wct[icol].split()[0]
    #             for F in PlotFmt:
    #                 matplotlib.pyplot.savefig(P+F, dpi=400)
    #             if PdfCatalog:
    #                 pdf_list.append(P+".pdf")


    #             matplotlib.pyplot.show()

    # if PlotIt:
    #     if PdfCatalog:
    #         viz.make_pdf_catalog(PdfList=pdf_list, FileName="NM_Histograms.pdf")


    # print ("\n\nSummary:")
    # Xcount = (XData.flat<Value).sum()
    # print ("number of values less than "+str(S*Value)+" in X data: "+str(Xcount/nd)
    #         +"  =  "+str(100*Xcount/nd)+"%")

    # Zcount = (ZData.flat<Value).sum()
    # print ("number of values less than "+str(S*Value)+" in Z data: "+str(Zcount/nd)
    #         +"  =  "+str(100*Zcount/nd)+"%")

    # Acount = 2*nd
    # Dcount =Xcount + Zcount
    # print ("number of values less than "+str(S*Value)+" in all data: "+str(Dcount/Acount)
    #         +"  =  "+str(100*Dcount/Acount)+"%")

    # print("time used = "+str(process_time()-start))
Fileout = OutDatDir + AEM_system.upper()+"_AllData.npz"
numpy.savez_compressed(file=Fileout, Data = DataAll)
print("Results stored to "+Fileout)
