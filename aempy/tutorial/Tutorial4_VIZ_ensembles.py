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
Show several 1d block models as (stitched) section.

"""
import os
import os.path

import sys
from sys import exit as error

from datetime import datetime
import warnings

import numpy

import matplotlib.collections
import matplotlib.patches
import matplotlib.colors
import matplotlib.pyplot
import matplotlib
import matplotlib.cm

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
# mypath = ["/home/vrath/AEMpyX/aempy/modules/", "/home/vrath/AEMpyX/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import viz
import eviz
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
"""
input format is npz
"""

InFileFmt = ".npz"
InModDir ="/home/vrath/work/AEM_Data//Projects/StGormans/results_jcn"
if not InModDir.endswith("/"): InModDir = InModDir+"/"
print("Data/models read from dir:  %s" % InModDir)


FileList = "set"  # "search", "read"

if "search" in FileList.lower():

    SearchStrng = ""
    print("Search flightline ID string: %s " % SearchStrng)
    data_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InModDir)
    data_files = sorted(data_files)


if "set" in FileList.lower():
   data_files =[
    "A1_rect_StGormans_FL11379-0_proc_delete_PLM3s_k3_nlyr36_TikhOpt-JCN_gcv_Prior100_results",
   ]

PlotDir = InModDir+"/plots/"
print("Plots written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)

FilesOnly = False

PlotFmt = [".pdf", "png"]
PdfCatalog = False

if ".pdf" in PlotFmt:
    PdfCatName = "JCNPlots.pdf"
else:
    print(" No pdfs generated. No catalog possible!")
    PdfCatalog = False


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

nD = NN[0]

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

"""
Placement and size of plots
"""   
Nplots = 2
Horiz = True

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

ifl = 0
pdf_list = []
for file in data_files:

    FileName, filext0 = os.path.splitext(file)

    title=FileName

    """
    numpy.savez_compressed(
        file=Fileout+".npz",
        fl_data=file,
        fl_name=fl_name,
        fl_orig=fl_orig,
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
        site_jacd= site_jacd,
        site_pcov= site_pcov,
        site_jcn_avg=site_jcn_avg,
        site_jcn_var=site_jcn_var,
        site_jcn_med=site_jcn_med,
        site_jcn_mad=site_jcn_mad)
           
        if "ens" in Ctrl["output"]:
            util.add_object_npz(filein=Fileout+".npz",
                       xkeys=["site_jcn_ens"], xobjects=[site_jcn_ens])
            
        """
        
    results = numpy.load(InModDir+file)
    
    fl_name = results["fl_name"]
    fl_orig = results["fl_orig"]
    
    site_num = results["site_num"]
    site_x   = results["site_x"] - fl_orig[0]
    site_y   = results["site_y"] - fl_orig[1]
          
   
   
    site_alt = results["site_alt"]
    
    num_sites = len(site_num)
    
    pdf_list = []
    for isit in numpy.arange(num_sites):
        
        
        ensemble = results["site_jcn_ens"][isit]
        
        
        
        pos = site_x[isit]**2 + site_y[isit]**2 
        PlotTitle = FileName+"  site "+str(isit)+" at position "+str(numpy.around(pos),0)+" m"

        
        if Horiz: 
            horz = Nplots
            vert = 1
        else:
            horz = 1
            vert = Nplots
            
        fig, ax = matplotlib.pyplot.subplots(1,nplots,
                                          figsize=(horz*PlotSize[0]*cm, vert*PlotSize[0]*cm),
                                          gridspec_kw={
                                              "height_ratios": [1.],
                                              "width_ratios": [1., 1.]})
        fig.suptitle(PlotTitle+" ("+method+")", fontsize=Fontsizes[2])
     
        
        ax[0] = eviz.plot_model_ensemble(
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
            ax[0] = eviz.plot_model(
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
        
        
        ax[1] = eviz.plot_data_ensemble(
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


        
        
        


    for F in PlotFmt:
      matplotlib.pyplot.savefig(PlotDir+FileName+"_model"+F, dpi=400)

    if matplotlib.get_backend()!="cairo":
        matplotlib.pyplot.show()
    matplotlib.pyplot.LogrhoLimitsf()


    if PdfCatalog:
        pdf_list.append(PlotDir+FileName+".pdf")
