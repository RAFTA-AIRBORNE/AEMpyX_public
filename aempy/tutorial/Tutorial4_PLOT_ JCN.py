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
#       jupytext_version: 1.14.7
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

PlotFmt = [".pdf"]
PdfCatalog = False
if ".pdf" in PlotFmt:
    PdfCStr = "JCN_Plots.pdf"
else:
    print(" No pdfs generated. No catalog possible!")
    PdfCatalog = False



AEM_system = "aem05"
# AEM_system = "genesis"

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    """
    Parameter for data plot
    """
    QLimits = [0., 2500.]
    ILimits = [0., 2500.]

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    """
    Parameter for data plot
    """
    XLimits = [-8., 12.]
    ZLimits = [-8., 12.]


nD = NN[0]


PlotType =  0      # model
# PlotType =  1      #  model + data

if PlotType==0:
    PlotSize = [8., 8.]
    numplots = 1
else:
    PlotSize = [8., 8.]
    numplots = 2




"""
Parameter for model plots
"""
MinLogrho =  1.
MaxLogrho =  4.
LogrhoLimits = [MinLogrho, MaxLogrho]
LogrhoTicks = [-1, 0, 1, 2, 3, 4]

# Percentiles = [10., 20., 30., 40.] # linear
Percentiles = [2.3, 15.9 ]                   # 95/68
PlotPerc = True
FillPerc = True


"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""

matplotlib.pyplot.style.use("seaborn-v0_8-paper") # ("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
# matplotlib.rcParams["text.usetex"] = True

Fontsize = 8
Labelsize = Fontsize
Titlesize = 4
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidth= 1.5
Markersize = 4

ncols = len(Percentiles)+1
Colors = matplotlib.pyplot.cm.jet(numpy.linspace(0,1,ncols))
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
        file=Fileout,
        fl_data=file,
        fl_name=fl_name,
        header=Header,
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
        site_site_y,
        site_site_x,
        site_gps=site_gps,
        site_alt=site_alt,
        site_dem=site_dem)


    """
    tmp = numpy.load(InModDir+file)

    m_act    = tmp["mod_act"]
    m_ref    = tmp["mod_ref"]
    m_ens    = tmp["site_modl"]


    d_act    = tmp["dat_act"]
    d_ens    = tmp["site_dobs"]


    nlyr = inverse.get_nlyr(m_ref)


    site_alt = tmp["site_alt"]


        # uncpars =\
        #     dict([
        #     ("jacd", Jd),               # jacobian
        #     ("cpost", C),               # cov a-post
        #     ("merr", E),                # model error
        #     ("mresm", [Rm, Sm, Nm]),    # model resolution
        #     ("dresm", [Rd, Sd, Nd]),    # data resolution
        #     ("gi", G),                  # generalized inverse
        #     ])
        # results.update(uncpars)


    viz.plot_ensemble(
                    PlotFile = PlotDir+FileName,
                    PlotTitle =file,
                    PlotType = 2,
                    PlotOrder= "horizontal",
                    PlotSize = [8.],
                    PlotFormat = ["png",],
                    System  = "aem05",
                    Mod0 = [],
                    Depth0 = [],
                    ModEns = [],
                    Depth = [],
                    DatEns = [],
                    Percentiles=[2.5, 16.],
                    Fillcolor=["0.8", "0.4"],
                    Alphas = [0.3 , 0.6],
                    Labels=[],
                    Linecolor=["k", "r", "g", "b", "y", "m"],
                    Linetypes=["-", ":", ";"],
                    Linewidth=[1., 1.5, 2.],
                    # Marker = ["v"],                    Markersize =[4],
                    Fontsizes=[10,10,12],
                    MLimits= [],
                    ZLimits= [],
                    DLimits=[],
                    XLimits = [],
                    PlotStrng="",
                    PlotTrue = True,
                    Save = True)

        
        
        
    """       
        
            PlotFile = PlotDir+FileName,
            PlotTitle = file,
            PlotFormat = ["png",],
            Mod0 = [],
            Depth0 = [],
            ModEns = [],
            Depth = [],
            Percentiles=[2.5, 16.],
            Alphas = [0.3 , 0.6],
            Labels=[],
            Colors=["k", "r", "g", "b", "y", "m"],
            Lines=["-", ":", ";"],
            Linewidths=[1., 1.5, 2.],
            Markersize = 4,
            Fontsizes=[10,10,12],
            MLimits= [],
            SLimits= [],
            DLimits= [],
            PlotStrng="")

    # models = numpy.shape(m_ens)
    # smpls = models[0]
    # param = models[1]

    # if PlotType in [0, 1]:
    #     thk = m_ref[6 * nlyr : 7 * nlyr - 1].reshape(1,-1)
    #     thk = numpy.vstack((thk.T))
    #     dpt= numpy.hstack((0.0, numpy.cumsum(thk)))
    #     res = numpy.log10(m_ens[:,:])

    #     m_avg = numpy.nanmean(res, axis=0)
    #     m_std = numpy.nanstd(res, axis=0)

    #     m_prc =numpy.percentile(res,50., axis=0)
    #     m_lab = ["50 % (median)"]
    #     for pp in Percentiles:
    #         pc = [pp, 100.-pp]
    #         m_prc = numpy.vstack([m_prc, numpy.percentile(res,pc, axis=0)])
    #         m_lab.append(str(pc[0])+"/"+str(pc[1])+" %")



    # if PlotType in [1]:
    #     d_avg = numpy.nanmean(d_ens, axis=0)
    #     d_std = numpy.nanstd(d_ens, axis=0)
    #     d_true   = tmp["site_dobs"][0]
    #     d_prc =numpy.percentile(res,50., axis=0)
    #     d_lab = ["50 % (median)"]
    #     for pp in Percentiles:
    #         pc = [pp, 100.-pp]
    #         d_prc = numpy.vstack(m_prc, numpy.percentile(res,pc, axis=0))
    #         d_lab.append(str(pc[0])+"/"+str(pc[1])+" %")

    # print (numpy.shape(m_prc))



    # print(PlotType)

    # nplots = numplots
    # if PlotType == 0:
    #     fig, ax = matplotlib.pyplot.subplots(1, nplots,
    #                                       figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
    #                                       gridspec_kw={"height_ratios": [1]},
    #                                       squeeze=False)

    # if PlotType == 1:

    #     fig, ax = matplotlib.pyplot.subplots(1, nplots,
    #                                       figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
    #                                       gridspec_kw={"height_ratios": [1, 1]},
    #                                       squeeze=False)

    # fig.suptitle(title, x=0.5, y=0.9, fontsize=Titlesize, va="baseline") #, fontweight="bold")


    # ax = ax.flatten()
    # ii = - 1
    # if PlotType in [0, 1]:

    #     ii = ii+1


    #     if PlotPerc:

    #         ax[ii].step(m_prc[0], dpt, linewidth= Linewidth, color= 'red', label=m_lab[0])
    #         for jj in range(1, len(m_prc)):
    #             ax[ii].step(m_prc[jj], dpt, linewidth= Linewidth, color= 'red', label=m_lab[0])
    #         # ax[ii].plot(site_r[:-1], site_rms[:-1], "r", linewidth=Linewidth)
    #         # ax[ii].plot(site_r[:-1], avg[:-1], "b:", linewidth=Linewidth)
    #         # ax[ii].plot(site_r[:-1], med[:-1], "g:", linewidth=Linewidth)
    #         # # ax[ii][0].set_title(title, fontsize=Fontsize+1)
    #         # ax[ii].legend(m_lab, fontsize=Labelsize, loc="best")
    #         ax[ii].set_ylabel("nRMS ", fontsize=Fontsize-1)
    #         ax[ii].set_ylim(rms_limits)
    #         ax[ii].grid(True)
    #         ax[ii].tick_params(labelsize=Labelsize)


    #         legend =ax[ii].legend(
    #             [" 0.9 kHz", "3 kHz", "12 kHz", "24.5 kHz"],
    #             fontsize=Labelsize-2, loc="best", ncol=2)
    #         legend.set_title("Frequency", prop={"size":Labelsize})

    # #     # if LogPlot:
    #     #     if LogSym:
    #     #         ax[ii].set_yscale("symlog", linthresh=LinThresh)
    #     #     else:
    #     #         ax[ii].set_yscale("log")
    #     # else:
    #     #     ax[ii].set_yscale("linear")
    #     if ii == nplots:
    #         ax[ii].set_xlabel("Profile distance "+ProfUnit, fontsize=Labelsize)
    #     else:
    #         ax[ii].tick_params(labelbottom=False)

    #     ii = ii+1

    #     if Quantile95:
    #         ax[ii].fill_between(
    #             site_r, I_obs[:, 0]-I_err95[:, 0], I_obs[:, 0]+I_err95[:, 0],
    #             color="r",alpha=alpha_err95)
    #         ax[ii].fill_between(
    #             site_r, I_obs[:, 1]-I_err95[:, 1], I_obs[:, 1]+I_err95[:, 1],
    #             color="g",alpha=alpha_err95)
    #         ax[ii].fill_between(
    #           site_r, I_obs[:, 2]-I_err95[:, 2], I_obs[:, 2]+I_err95[:, 2],
    #           color="b",alpha=alpha_err95)
    #         ax[ii].fill_between(
    #               site_r, I_obs[:, 3]-I_err95[:, 3], I_obs[:, 3]+I_err95[:, 3],
    #               color="m",alpha=alpha_err95)
    #     if Quantile68:
    #         ax[ii].fill_between(
    #             site_r, I_obs[:, 0]-I_err68[:, 0], I_obs[:, 0]+I_err68[:, 0],
    #             color="r",alpha=alpha_err68)
    #         ax[ii].fill_between(
    #             site_r, I_obs[:, 1]-I_err68[:, 1], I_obs[:, 1]+I_err68[:, 1],
    #             color="g",alpha=alpha_err68)
    #         ax[ii].fill_between(
    #           site_r, I_obs[:, 2]-I_err68[:, 2], I_obs[:, 2]+I_err68[:, 2],
    #           color="b",alpha=alpha_err68)
    #         ax[ii].fill_between(
    #               site_r, I_obs[:, 3]-I_err68[:, 3], I_obs[:, 3]+I_err68[:, 3],
    #               color="m",alpha=alpha_err68)


    #     ax[ii].plot(
    #          site_r, I_obs[:, 0],
    #          color="r", linewidth=Linewidth, linestyle="-")
    #     ax[ii].plot(
    #         site_r, I_obs[:, 1],
    #         color="g", linewidth=Linewidth, linestyle="-")
    #     ax[ii].plot(
    #         site_r, I_obs[:, 2],
    #         color="b", linewidth=Linewidth, linestyle="-")
    #     ax[ii].plot(
    #         site_r, I_obs[:, 3],
    #         color="m", linewidth=Linewidth, linestyle="-")
    #     ax[ii].plot(
    #         site_r, I_cal[:, 0],
    #         color="r", linewidth=Linewidth, linestyle=":")
    #     ax[ii].plot(
    #         site_r, I_cal[:, 1],
    #         color="g", linewidth=Linewidth, linestyle=":")
    #     ax[ii].plot(
    #         site_r, I_cal[:, 2],
    #         color="b", linewidth=Linewidth, linestyle=":")
    #     ax[ii].plot(
    #         site_r, I_cal[:, 3],
    #         color="m", linewidth=Linewidth, linestyle=":")


    #     # ax[ii].errorbar(
    #     #     site_r, I_obs[:, 0], yerr=I_err95[:,0],
    #     #     color="r",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
    #     # ax[ii].errorbar(
    #     #     site_r, I_obs[:, 1], yerr=I_err95[:,1],
    #     #     color="g",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
    #     # ax[ii].errorbar(
    #     #     site_r, I_obs[:, 2], yerr=I_err95[:,2],
    #     #     color="b",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
    #     # ax[ii].errorbar(
    #     #     site_r, I_obs[:, 3], yerr=I_err95[:,3],
    #     #     color="k",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)


    #     ax[ii].set_title("In-Phase", fontsize=Fontsize, y=1.0, pad=-(Fontsize+2))
    #     ax[ii].set_ylim([I_min, I_max])
    #     if numpy.size(ILimits)>0:
    #         ax[ii].set_ylim(ILimits)
    #     ax[ii].set_ylabel("(ppm)", fontsize=Fontsize)



    #     ax[ii].grid(True)

    #     legend =ax[ii].legend(
    #         [" 0.9 kHz", "3 kHz", "12 kHz", "24.5 kHz"],
    #         fontsize=Labelsize-2, loc="best", ncol=2)
    #     legend.set_title("Frequency", prop={"size":Labelsize})

    #     # if LogPlot:
    #     #     if LogSym:
    #     #         ax[ii].set_yscale("symlog", linthresh=LinThresh)
    #     #     else:
    #     #         ax[ii].set_yscale("log")
    #     # else:
    #     #     ax[ii].set_yscale("linear")

    #     norm = matplotlib.colors.Normalize(vmin=LogrhoLimits[0], vmax=LogrhoLimits[1], LogrhoLimitsip=False)
    #     cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=mycmap), cax=ax[ii]
    #                       ticks=LogrhoTicks ,
    #                       orientation="horizontal", aspect=40, pad=0.25, shrink=0.65)
    #     cb.set_label(r"log10($\Omega$ m)", size=Fontsize)
    #     cb.ax.tick_params(labelsize=Fontsize)
    """

    for F in PlotFmt:
      matplotlib.pyplot.savefig(PlotDir+FileName+"_model"+F, dpi=400)

    if matplotlib.get_backend()!="cairo":
        matplotlib.pyplot.show()
    matplotlib.pyplot.LogrhoLimitsf()


    if PdfCatalog:
        pdf_list.append(PlotDir+FileName+".pdf")
