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
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (Spyder)
#     language: python3
#     name: python3
# ---

# +
#!/usr/bin/env python3
# -


# This script is used for the visualization of AEM data along a flightline.


# +
import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime

import warnings
import getpass

from cycler import cycler

import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.backends.backend_pdf #  matplotlib.backends. backend_pdf.PdfPages

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import prep
import aesys
import viz
# -

version, _ = versionstrg()
# script = "VIZ_data_flightline.py"
script = __file__  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")
Header = titstrng


OutInfo = False
AEMPYX_DATA = os.environ["AEMPYX_DATA"]


# The following cell gives values to AEM-system related settings.
#
# Data transformation is activated by the variable DataTrans. Currently three possible options are allowed: _DataTrans = 0_: No transformation, i.e., the raw data are used. _DataTrans = 1_: The natural log of data is taken, only allowed for strictly positive values. _DataTrans = 2_: If data scale logarithmically, an asinh transformation (introduced by Scholl, 2000) is applied. It allows negatives, which may occur in TDEM, when IP effects are present.
#
# A general additive/multiplicative error model is applied on the raw data before transformation, and errors are also transformed.
#


# +
# AEM_system = "genesis"
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    _ ,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.05
    data_active = numpy.ones(NN[2], dtype="int8")

if "genes" in AEM_system.lower():
    _ , NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + "good" hoizontals"
# -


# +
InFileFmt = ".npz"
AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/SmallTest/"
FileList = "search"


# un/comment according to which data  you want to plot
# # raw data
# InDatDir =  AEMPYX_DATA + "/raw/"
# PlotDir =  InDatDir + "/plots/"
# SearchStrng = "*FL*.npz"
# PlotStrng = " - raw"
# PDFCatName = PlotDir+"StGormans_raw.pdf"
# +


# processed data
InDatDir =  AEMPYX_DATA + "/data/"
PlotDir =  InDatDir + "/plots/"
SearchStrng = "*FL*nan*.npz" # if no interpolation was chosen
#SearchStrng = "*FL*.npz" # else
PlotStrng = " - proc"



# +

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = []

else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(dat_files)
# -

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")

if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)

# The next block determines the graphical output. if _PDFCatalog_ is set, a catalogue
# including all generated figures, named _PDFCatName_ is generated. This option is only
# available if ".pdf" is included in the output file format list (_PlotFmt_).

FilesOnly = False    # for headless plotting.
PlotFmt = [".png", ".pdf"]

PDFCatName = PlotDir+"StGormans_processed.pdf"
PDFCatalog = True
if ".pdf" in PlotFmt:
    pass
else:
    print(" No pdf files generated. No catalog possible!")
    PdfCatalog = False

if "aem05" in AEM_system.lower():
    IncludePlots = ["alt", "qdata", "idata",]
    # IncludePlots = ["qdata", "idata",]
    QLimits = []
    ILimits = []
    PlotSize = [18., 6.]
    PLimits = [0., 1.]
    HLimits = [30., 90.] #[40., 140.]
    LogPlot = False
    LogSym = False
    LinThresh =10.
    if LogPlot == False:
        LogSym = False
    Logparams=[LogPlot, LogSym, LinThresh]

if "genes" in AEM_system.lower():
    IncludePlots = ["alt", "xdata", "zdata",]
    # IncludePlots = ["xdata", "zdata",]
    PlotSize = [18., 6.]
    DataTrans = "asinh"
    XLimits = [3.5, 12.]
    ZLimits = [6., 14.]

    LogPlot = False
    LogSym = False
    LinThresh =100.
    if LogPlot == False:
        LogSym = False
    Logparams=[LogPlot, LogSym, LinThresh]


    HLimits = [80., 320.]
    PLimits = [0., 25.]

# +
PosDegrees = False
if PosDegrees:
    EPSG=32629
PlotThresh =10

ProfType = "distance"
if "dist" in ProfType.lower():
    ProfLabel = "profile distance (m) "
    ProfScale = 1. # 0.001  # m to km
else:
    ProfLabel = "site # (-)"
    ProfScale = 1. # 0.001  # m to km
# -


# This block sets graphical parameters related to the \textit{matplotlib}.
# package. A list of available plotting styles can be found on matplotlib"s
# website at https://matplotlib.org/stable/users/explain/customizing.htm, or
# entering the python command
# _print(matplotlib.pyplot.style.available)} in an appropriate_
# window.
#

matplotlib.pyplot.style.use("seaborn-v0_8-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
matplotlib.rcParams["savefig.transparent"] = True
matplotlib.rcParams["savefig.bbox"] = "tight"
Fontsize = 8
Labelsize = Fontsize
Titlesize = 8
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidths= [0.6]
Markersize = 4

ncols = 11
Colors = matplotlib.pyplot.cm.jet(numpy.linspace(0,1,ncols))
Grey = 0.7
Lcycle = (cycler("linestyle", ["-", "--", ":", "-."])
          * cycler("color", ["r", "g", "b", "y"]))

# For just plotting to files ("headless plotting"), choose the
# cairo backend (eps, pdf, png, jpg...).

if FilesOnly:
    matplotlib.use("cairo")

if PDFCatalog:
    pdf_list = []
    catalog =matplotlib.backends.backend_pdf.PdfPages(PDFCatName)

# Depending on the region of interest, the number of plots may be quite large.

ifl = 0
for file in dat_files:

    ifl = ifl+1

    name, ext = os.path.splitext(file)

    filein = InDatDir+file
    Data, Header, _ = aesys.read_aempy(File=filein,
                                   System=AEM_system, OutInfo=False)
    sD = numpy.shape(Data)
    print("flightline "+name+"  #"
          +str(ifl)+" of "
          +str(numpy.size(dat_files)) +" has shape: "+str(sD))
    print(sD, nD )
    if numpy.shape(Data)[0]<=nD:
        print("Not enough data! Not plotted")
        continue

    anynan = numpy.argwhere(numpy.isnan(Data))
    nnans = numpy.shape(anynan)[0]
    for ii in anynan:
        Data[ii[0],3:] = numpy.nan

    # if numpy.shape(Data)[0]-nnans < PlotThresh:
    #     print("Not enough data! Not plotted")
    #     continue


    if PDFCatalog:
        pdf_list.append(PlotDir+name+".pdf")

    fline = Data[:, 0]
    Data[:, 1] = Data[:, 1] * ProfScale
    Data[:, 2] = Data[:, 2] * ProfScale

    if "aem05" in AEM_system.lower():
        fig, _ = viz.plot_flightline_aem05(
                PlotName = name,
                PlotDir = PlotDir,
                PlotFmt=PlotFmt,
                IncludePlots=IncludePlots,
                PlotSize = PlotSize,
                DataObs=Data,
                DataCal=[],
                QLimits =[],
                ILimits =[],
                DLimits = [],
                HLimits = HLimits,
                PLimits = PLimits,
                ProfLabel=ProfLabel,
                Linecolor=Colors,
                Linewidth=Linewidths,
                Fontsizes=Fontsizes,
                Logparams=Logparams,
                PlotStrng=PlotStrng,
                PlotPLM = True)

    if "genes" in AEM_system.lower():
        fig, _ = viz.plot_flightline_genesis(
                PlotName = name,
                PlotDir = PlotDir,
                PlotFmt=PlotFmt,
                IncludePlots=IncludePlots,
                PlotSize = PlotSize,
                DataObs=Data,
                DataCal=[],
                DataTrans = DataTrans,
                DLimits = [],
                XLimits =XLimits,
                ZLimits =ZLimits,
                HLimits =[],
                ProfLabel=ProfLabel,
                Linecolor=Colors,
                Linewidth=Linewidths,
                Fontsizes=Fontsizes,
                Logparams=Logparams,
                PlotStrng=PlotStrng)

    if PDFCatalog:
        pdf_list.append(PlotDir+name+".pdf")
        catalog.savefig(fig)



if PDFCatalog:
    print(pdf_list)
    # viz.make_pdf_catalog(PDFList=pdf_list, FileName=PDFCatName)
    d = catalog.infodict()
    d["Title"] =  name
    d["Author"] = getpass.getuser()
    d["CreationDate"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    catalog.close()
