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
#       jupytext_version: 1.16.2
# ---


import os
import os.path

import sys
from sys import exit as error

from datetime import datetime
import warnings
import getpass

import numpy

import matplotlib.collections
import matplotlib.patches
import matplotlib.colors
import matplotlib.pyplot
import matplotlib
import matplotlib.cm

import matplotlib.backends.backend_pdf #  matplotlib.backends. backend_pdf.PdfPages

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
# mypath = ["/home/vrath/AEMpyX/aempy/modules/", "/home/vrath/AEMpyX/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import viz
import viz
import inverse


warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
cm = 1/2.54  # centimeters in inches

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")


now = datetime.now()
AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/Synth/"
InModDir = AEMPYX_DATA+"/results/"
print("Data/models read from dir:  %s" % InModDir)

FileList = "search"  # "search", "read"

if "search" in FileList.lower():

    SearchStrng = "*results.npz"
    print("Search flightline ID string: %s " % SearchStrng)
    data_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InModDir)
    data_files = sorted(data_files)


if "set" in FileList.lower():
   data_files =[]

PlotDir = AEMPYX_DATA+"/plots/"
print("Plots written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)

FilesOnly = False

PlotFormat = [".pdf", ".png"]

PDFCatalog = False
PDFCatName = "AEM05_EnsemblePlots.pdf"
if ".pdf" in PlotFormat:
    pass
else:
    print(" No pdfs generated. No catalog possible!")
    PDFCatalog = False

PlotTrue = True #True

"""
Placement of plots
"""
Horiz = True

"""
Parameter for data plot
"""

AEM_system = "aem05"
# AEM_system = "genesis"

FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)

if "aem05" in AEM_system.lower():

    DataLimits = [0., 2500.]
    FreqLimits = []


if "genes" in AEM_system.lower():
    # asinh trans (if negatives)
    DataLimits = []
    TimeLimits = []


"""
Parameter for model plots
"""
ModLimits = [3., 3000.]
DepthLimits = [0., 100.]

# Percentiles = [10., 20., 30., 40.] # linear
Percentiles = [2.3, 15.9 ]                   # 95/68



PlotTitle = "Aem05: 3-layer Model"
PlotSize = [8., 8.]

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
see:
MatplotlibDeprecationWarning: The seaborn styles shipped by Matplotlib
are deprecated since 3.6, as they no longer correspond to the styles s
hipped by seaborn. However, they will remain available as
'seaborn-v0_8-<style>'. Alternatively, directly use the seaborn API instead.

"""
matplotlib.pyplot.style.use("seaborn-v0_8-paper") # ("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
# matplotlib.rcParams["text.usetex"] = True

Fontsize = 8
Labelsize = Fontsize
Titlesize = 12
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidth = [1., 1., 0.75]
Linetypes = ["-", ":", "--", ";"]
Linecolors = ["k", "r", "g", "b", "c", "m"]

Markers = ["o"]
Markersize = [5]


ncols = len(Percentiles)+3

Greys = [0.1, 0.2, 0.3, 0.4, 0.5,.6, 0.7, 0.8, 0.9]

"""
see:
https://matplotlib.org/stable/gallery/color/colormap_reference.html
"""
cmap = matplotlib.cm.gist_rainbow
mycmap = matplotlib.pyplot.get_cmap(cmap)
"""
For just plotming to files, choose the cairo backend (eps, pdf, ,png, jpg...).
If you need to see the plot directly (plot window, or jupyter), simply
comment out the following line. In this case matplotlib may run into
memory problems after a few hundreds of high-resolution plot.
"""
if FilesOnly:
   matplotlib.use("cairo")


ns = numpy.size(data_files)
if PDFCatalog:
    pdf_list = []
    catalog =matplotlib.backends.backend_pdf.PdfPages(PDFCatName)

for file in data_files:

    filename, filext0 = os.path.splitext(file)

    PlotFile = filename

    results = numpy.load(InModDir+file)

    print("Data loaded from file ", InModDir+file)


    m_act    = results["mod_act"]
    m_ref    = results["mod_ref"]
    m_ens    = results["ens_modl"]

    print("minimum: ",numpy.amin(m_ens))
    print("maximum: ",numpy.amax(m_ens))

    nlyr     = inverse.get_nlyr(m_ref)
    dz       = m_ref[6*nlyr:7*nlyr-1]
    z_ens = inverse.set_znodes(dz)



    d_act    = results["dat_act"]
    d_ens    = results["ens_dcal"]
    m_alt    = results["mod_alt"]

    r_ens    = results["stat_nrms"]

    control = numpy.load((InModDir+file).replace("_results","_ctrl"), allow_pickle=True)
    method =  control["inversion"][1]


    nlyr = inverse.get_nlyr(m_ref)

    ens_modl = \
        inverse.calc_stat_ens(ensemble=m_ens, quantiles=Percentiles, sum_stats=True)
    ens_dcal = \
        inverse.calc_stat_ens(ensemble=d_ens, quantiles=Percentiles, sum_stats=True)
    ens_nrms = \
        inverse.calc_stat_ens(ensemble=r_ens, quantiles=Percentiles, sum_stats=True)


    if PlotTrue:
        m_true = results["mod_true"]
        d_true = results["dat_true"]
        print("read: ",m_true)
        l_true = inverse.get_nlyr(m_true)
        z_true = inverse.set_znodes(m_true[6* l_true:7* l_true-1])
        z_true = numpy.append(z_true, 10000.)

        m_true = m_true[0*l_true:1*l_true]
        m_true = numpy.append(m_true, m_true[-1])
        print("calc: ", m_true)
        print("calc: ", z_true)



    nplots = 2
    if Horiz:
        horz = nplots
        vert = 1
    else:
        horz = 1
        vert = nplots

    fig, ax = matplotlib.pyplot.subplots(1,nplots,
                                      figsize=(horz*PlotSize[0]*cm, vert*PlotSize[0]*cm),
                                      gridspec_kw={
                                          "height_ratios": [1.],
                                          "width_ratios": [1., 1.]})
    fig.suptitle(PlotTitle+" ("+method+")", fontsize=Fontsizes[2])


    ax[0] = viz.plot_model_ensemble(
            ThisAxis = ax[0],
            PlotType = "percentiles", # lines, percentiles. iso
            System  = AEM_system,
            ModEns = m_ens,
            Depth = z_ens,
            Percentiles=[2.5, 16.],
            Fillcolor=["0.8", "0.4"],
            Alphas = [0.3 , 0.6],
            Labels=[],
            Linecolor=Linecolors,
            Linetype=Linetypes,
            Linewidth=Linewidth,
            Markers = ["v"],
            Markersize =[4],
            Fontsizes=Fontsizes,
            XLimits=ModLimits,
            YLimits= DepthLimits,
            Legend=False)

    if PlotTrue:
        # print(m_true)
        # print(z_true)
        ax[0] = viz.plot_model(
                ThisAxis = ax[0],
                System  = AEM_system,
                Model = m_true,
                Depth = z_true,
                Labels=["true model"],
                Linecolor=["k"],
                Linetype=["--"],
                Linewidth=[1],
                Markers = ["v"],
                Markersize =[4],
                Fontsizes=Fontsizes,
                XLimits= ModLimits,
                YLimits= DepthLimits,
                Legend=True)


    ax[1] = viz.plot_data_ensemble(
            ThisAxis = ax[1],
            PlotType = "percentiles", # lines, percentiles. iso
            System  = AEM_system,
            DatEns = d_ens,
            Percentiles=[2.5, 16.],
            Fillcolor=["0.8", "0.4"],
            Alphas = [0.3 , 0.6],
            Labels=[],
            Linecolor=Linecolors,
            Linetype=Linetypes,
            Linewidth=Linewidth,
            Markers = [""],
            Markersize =[4],
            Fontsizes=Fontsizes,
            XLimits= FreqLimits,
            YLimits= DataLimits,
            Legend=False)

    if PlotTrue:

        ax[1] = viz.plot_data(
                ThisAxis = ax[1],
                System  = AEM_system,
                Data = d_true,
                Errs = [],
                Linecolor=["k","r","g","b"],
                Linetype=[""],
                Linewidth=Linewidth,
                Markers = ["s", "o"],
                Markersize = Markersize,
                Fontsizes=Fontsizes,
                XLimits= FreqLimits,
                YLimits= DataLimits,
                Legend=True)

    for F in PlotFormat:
        matplotlib.pyplot.savefig(PlotDir+PlotFile+F)

    if PDFCatalog:
        pdf_list.append(PlotDir+PlotFile+".pdf")
        catalog.savefig(fig)

if PDFCatalog:
    print(pdf_list)
    # viz.make_pdf_catalog(PDFList=pdf_list, FileName=PDFCatName)
    d = catalog.infodict()
    d["Title"] =  PDFCatName
    d["Author"] = getpass.getuser()
    d["CreationDate"] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    catalog.close()
