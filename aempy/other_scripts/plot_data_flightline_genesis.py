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
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)


from version import versionstrg
import util
import prep
import aesys
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=inspect.getfile(inspect.currentframe()), out=False)
print(titstrng+"\n\n")

OutInfo = True
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

AEM_system = "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  75.
    DatErr_mult = 0.05
    data_active = numpy.ones(NN[2], dtype="int8")


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'


"""
input formats are .npz, .nc4, 'ascii'
"""

InStrng = ""

FileList = "search"  # "search", "read"
InDatDir =  AEMPYX_DATA + "/Projects/Compare_systems/data_reduced/"
PlotDir  =  AEMPYX_DATA + "/Projects/Compare_systems/plots/"
PlotStrng = " - data reduced"
SearchStrng = "CGG*reduced.npz"
PdfCat = "Testlines_CGG_reduced.pdf"

CheckNaN = False

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = []
    f = numpy.load(AEMPYX_DATA+"/Projects/Compare_systems/BundoranSubsets.npz")

    dat_files = f["setC"][0]
    dat_files = [os.path.basename(d) for d in dat_files]
else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(dat_files)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")

print(dat_files)


FilesOnly = False
PlotFmt = [".pdf", ".png"]
PdfCatalog = True
if ".pdf" in PlotFmt:
    PdfCStr = SearchStrng.replace("*","").replace("_","")

else:
    print(" No pdfs generated. No catalog possible!")
    PdfCatalog = False


IncludePlots = ["alt", "xdata", "zdata",]
IncludePlots = ["xdata", "zdata",]
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

PlotThresh =20


ProfType = "distance"
if "dist" in ProfType.lower():
    ProfLabel = "profile distance (m) "
    ProfScale = 1. # 0.001  # m to km
else:
    ProfLabel = "site # (-)"
    ProfScale = 1. # 0.001  # m to km


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
Lcycle =Lcycle = (cycler("linestyle", ["-", "--", ":", "-."])
          * cycler("color", ["r", "g", "b", "y"]))

"""
For just plot.ing to files, choose the cairo backend (eps, pdf, png, jpg...).
If you need to see the plots directly (plot window, or jupyter), simply
comment out the following line. In this case matplotlib may run into
memory problems after a few hundreds of high-resolution plot..
"""

if FilesOnly:
    print("using cairo")
    matplotlib.use("cairo")


if PdfCatalog:
    pdf_list = []

ifl = 0
for file in dat_files:

    ifl = ifl+1
    name, ext = os.path.splitext(file)

    filein =InDatDir + file
    Data, Header, _ = aesys.read_aempy(File=filein,
                                    System=AEM_system, OutInfo=False)
    sD = numpy.shape(Data)
    print("flightline "+name+"  #"
            +str(ifl)+" of "
            +str(numpy.size(dat_files)) +" has shape: "+str(sD))

    sD = numpy.shape(Data)
    if CheckNaN:
        anynan = numpy.argwhere(numpy.isnan(Data))
        nnans = numpy.shape(anynan)[0]
        for ii in anynan:
            Data[ii[0],3:] = numpy.nan

        if numpy.shape(Data)[0]-nnans < PlotThresh:
            print("Not enough data! No plot.")
            continue

    if numpy.size(Data)<=nD:
        print("Not enough data! No plot.")
        continue

    if PdfCatalog:
        pdf_list.append(PlotDir+name+".pdf")

    fline = Data[:, 0]

    Data[:, 1] = Data[:, 1] * ProfScale
    Data[:, 2] = Data[:, 2] * ProfScale

    # not usually used
    # altitude = Data[:, 4]
    # alt_file = name+"_altitude.txt"
    # numpy.savetxt(alt_file, altitude)

    viz.plot_flightline_genesis(
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
        Colors=Colors,
        Linewidths=Linewidths,
        Fontsizes=Fontsizes,
        Logparams=Logparams,
        PlotStrng=PlotStrng)

if PdfCatalog:
    viz.make_pdf_catalog(PdfList=pdf_list, FileName=PdfCat)
