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
Show several 1d block models as (stitched) section.

"""
import os
import sys
from sys import exit as error
from datetime import datetime
import warnings

import numpy

import matplotlib.collections
import matplotlib.patches
import matplotlib.colors as col
import matplotlib.pyplot
import matplotlib

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
# mypath = ["/home/vrath/AEMpyX/aempy/modules/", "/home/vrath/AEMpyX/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)

from version import versionstrg


import util
import viz
import inverse


warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Plot  Inversion results"+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

cm = 1/2.54  # centimeters in inches

OutInfo = True
now = datetime.now()

"""
input formats are "npz","nc4","asc"
"""
InFileFmt = ".npz"
InModDir = AEMPYX_DATA + "/Nearest/fwd_compare/models/"
print("Data read from dir:  %s" % InModDir)

Search = False
if Search:
    SearchStrng = "*k2*gcv*"
    print("Search flightline ID string: %s " % SearchStrng)
    data_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InModDir+"/models/")
    data_files = sorted(data_files)
else:
    SearchStrng = ""
    data_files =["A2_NM_intersection_FL21514-0_delete_PLM10_TikhOpt_lcc_Results.npz",
                 "A2_NM_intersection_FL21514-0_delete_PLM10_TikhOpt_gcv_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_TikhOpt_lcc_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_TikhOpt_gcv_Results.npz",
                 "A2_NM_intersection_FL21514-0_delete_PLM10_k3_TikhOpt_lcc_Results.npz",
                 "A2_NM_intersection_FL21514-0_delete_PLM10_k3_TikhOpt_gcv_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_k3_TikhOpt_fix_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_k3_TikhOpt_lcc_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_k3_TikhOpt_gcv_Results.npz",
                 "A2_NM_intersection_FL21514-0_delete_PLM10_k1_TikhOpt_fix_Results.npz",
                 "A2_NM_intersection_FL21514-0_delete_PLM10_k1_TikhOpt_lcc_Results.npz",
                 "A2_NM_intersection_FL21514-0_delete_PLM10_k1_TikhOpt_gcv_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_k1_TikhOpt_fix_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_k1_TikhOpt_lcc_Results.npz",
                 "A2_NM_intersection_FL21513-0_delete_PLM10_k1_TikhOpt_gcv_Results.npz"]


PlotDir = AEMPYX_DATA + "/Nearest/fwd_compare/plots/"
print("Plots written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)

FilesOnly = False

PlotFmt = [".pdf"]
PdfCatalog = True
if ".pdf" in PlotFmt:
    PdfCStr = SearchStrng.replace("*","") #.replace("_","")
else:
    print(" No pdfs generated. No catalog possible!")
    PdfCatalog = False


PlotType =  0      # rms, datafit, model
# PlotType =  1      #  model + rms
# PlotType =  2      #  model
# PlotType =  3      #  datafit
PlotSize = [25., 5. ]


"""
Parameter for nRMS plot
"""
rms_limits  =[0., 4.]


"""
Parameter for model plot
"""
min_lrho =  1.
max_lrho =  4.
sl = [min_lrho, max_lrho]

plot_min = -10.0
plot_max = 130.0
blank = 10

low_sens = True
if low_sens:
    lowsens = -4.
    alpha_sens = 0.3

high_rms = True
if high_rms:
    highrms = 3.
    alpha_rms= 0.2

high_err = False
if high_err:
    higherr = 0.5
    alpha_err = 0.3

max_doi = True
if max_doi:
    maxdoi = 100.
    alpha_doi = 0.0


plot_adapt = True
if plot_adapt:
   if not max_doi:
       maxdoi = 100.
else:
    plot_min = -10.0
    plot_max = 130.0

topo_use_average = False #
scale_sens = ["size", "max"]
strng_sens = " sensitivity scaling = "+str(scale_sens)+", thresh ="+str(lowsens)
cb_ticks = [-1, 0, 1, 2, 3, 4]

"""
Parameter for data fit plot
"""


"""
General Parameter for all plots
"""
poslatlon = True
if poslatlon:
    EPSG=32629


ProfScale = 1. # 0.001  # m to km
ProfUnit  = "(m)" #

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""

matplotlib.pyplot.style.use("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
# matplotlib.rcParams["text.usetex"] = True

Fontsize = 8
Labelsize = Fontsize
Titlesize = 8
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidth= 1
Markersize = 4

ncols = 11
Colors = matplotlib.pyplot.cm.jet(numpy.linspace(0,1,ncols))
Grey = 0.7

"""
see:
https://matplotlib.org/stable/gallery/color/colormap_reference.html
"""
mycmap = matplotlib.cm.gist_rainbow

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

    nplots = 0

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
        site_y=site_y,
        site_x=site_x,
        site_gps=site_gps,
        site_alt=site_alt,
        site_dem=site_dem)


    """
    tmp = numpy.load(InModDir+file)

    m_active    = tmp["mod_act"]
    mod_ref     = tmp["mod_ref"]
    site_model  = tmp["site_modl"]
    site_error  = tmp["site_merr"]
    site_sens   = tmp["site_sens"]

    site_dobs   = tmp["site_dobs"]
    site_dcal   = tmp["site_dcal"]
    site_err    = tmp["site_derr"]

    nlyr = inverse.get_nlyr(mod_ref)

    site_rms = tmp["site_nrms"]

    site_x = tmp["site_x"] * ProfScale
    site_y = tmp["site_y"] * ProfScale

    site_alt = tmp["site_alt"]
    site_gps = tmp["site_gps"]
    site_dem = tmp["site_dem"]

    if topo_use_average:
        site_tref = numpy.mean(site_dem)
    else:
        site_tref = 0
    site_topo = site_dem - site_tref
    max_topo = numpy.amax(site_topo)

    models = numpy.shape(site_model)
    sites = models[0]
    param = models[1]

    beg_pos = [site_x[0],  site_y[0],  site_topo[0] ]
    end_pos = [site_x[-1], site_y[-1], site_topo[-1]]
    beg_strng  = f"Start:\nX: {beg_pos[0]:.0f} \nY: {beg_pos[1]:.0f} (EPSG="+str(EPSG)+")"
    end_strng  = f"End:\nX: {end_pos[0]:.0f} \nY: {end_pos[1]:.0f}"

    if poslatlon:
        beg_pos = util.project_utm_to_latlon(beg_pos[0], beg_pos[1], utm_zone=EPSG)
        end_pos = util.project_utm_to_latlon(end_pos[0], end_pos[1], utm_zone=EPSG)
        beg_strng  = f"Start\nLat: {beg_pos[0]:.5f} \nLon: {beg_pos[1]:.5f}"
        end_strng  = f"End\nLat: {end_pos[0]:.5f} \nLon: {end_pos[1]:.5f}"


    site_x = site_x - site_x[0]
    site_y = site_y - site_y[0]
    site_r = numpy.sqrt(numpy.power(site_x, 2.0) + numpy.power(site_y, 2.0))


    if PlotType in [0, 1, 2]:


        dxmed2 = numpy.median(numpy.diff(site_r)) / 2.0

        thk = mod_ref[6 * nlyr : 7 * nlyr - 1].reshape(1,-1)
        thk = numpy.vstack((thk.T, thk[-1,0]))
        z0 = numpy.hstack((0.0, numpy.cumsum(thk)))
        zm = 0.5 * (z0[0:nlyr] + z0[1 : nlyr + 1])
        size_lay = numpy.repeat(thk.T, sites, axis=0)
        site_val = numpy.log10(site_model[:,:])

        if "size" in scale_sens:
            site_sens = site_sens/size_lay

        if "max" in scale_sens:
            #site_sens = numpy.sqrt(site_sens)
            scale = numpy.amax(numpy.abs(site_sens.flat))
            site_sens = site_sens/scale

        site_sens = numpy.log10(site_sens[:,:])


        alpha = numpy.ones_like(site_val)
        site_dpth = numpy.ones_like(site_topo)
        patches = []
        for nmod  in  numpy.arange(sites):

            zmp = site_topo[nmod] - zm
            z0p = site_topo[nmod] - z0
            site_dpth = numpy.amin(z0p)

            if high_rms:
                if site_rms[nmod] > highrms:
                    alpha[nmod, :] = alpha_rms

            if low_sens:
                 for il in numpy.arange(nlyr):
                     if site_sens[nmod, il] < lowsens:
                             alpha[nmod, il] = alpha_sens

            if high_err:
                for il in numpy.arange(nlyr):
                    if numpy.abs(site_err[nmod,il]) > higherr:
                            alpha[nmod, il] = alpha_err

            if max_doi:
                for il in numpy.arange(nlyr):
                    if z0p[0]-z0p[il] > maxdoi:
                            alpha[nmod, il] = alpha_doi


            for j in numpy.arange(nlyr):
                rect = matplotlib.patches.Rectangle((site_r[nmod] - dxmed2, z0p[j]), dxmed2 * 2, thk[j])
                patches.append(rect)

            site_val[nmod, -1] = numpy.nan


        p_coll = matplotlib.collections.PatchCollection(patches,
                                                        cmap=matplotlib.pyplot.get_cmap(mycmap),
                                                        alpha=alpha, linewidths=0, edgecolor=None)
        p_coll.set_clim(sl)
        p_coll.set_array(site_val.ravel())

        max_topo = numpy.amax(site_topo)
        min_topo = numpy.amin(site_topo)

        if plot_adapt:
            min_dpth = min_topo - maxdoi
            plotmax = max_topo + blank
            plotmin = min_dpth - blank
        else:
            plotmax = plot_max
            plotmin = plot_min

    if PlotType in [0, 3]:
        #
        print("coming soon!")

    if PlotType in [0, 1]:

        avg_rms = round(numpy.nanmean(site_rms), 2)
        med_rms = round(numpy.nanmedian(site_rms),2)
        med = med_rms*numpy.ones_like(site_rms)
        avg = avg_rms*numpy.ones_like(site_rms)

    print(PlotType)
    """
    PlotType =  0           rms, datafit, model
    PlotType =  1           model + rms
    PlotType =  2           model
    PlotType =  3           datafit
    """
    # fig =  matplotlib.pyplot.figure(figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm))

    fig = matplotlib.pyplot.figure()
    fig.suptitle(title, fontsize=Titlesize, fontweight="bold")


    if PlotType == 0:
        nplots = 4
        figs = matplotlib.gridspec.GridSpec(nplots, 1,  height_ratios=[1, 2, 2, 3, 0.2], sharex=True)
        # fig = matplotlib.pyplot.subplots(nplots, 1,
        #                                   figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
        #                                   sharex=True,
        #                                   gridspec_kw={"height_ratios": [1, 2, 2, 3, 0.2]})
    if PlotType == 1:
        nplots = 2

    if PlotType == 2:
        nplots = 1

    if PlotType == 3:
        nplots = 2



    if PlotType in [0, 1]:

        ax = fig.add_subplot(figs[0])

        ax.plot(site_r[:-1], site_rms[:-1], "r", linewidth=Linewidth)
        ax.plot(site_r[:-1], avg[:-1], "b:", linewidth=Linewidth)
        ax.plot(site_r[:-1], med[:-1], "g:", linewidth=Linewidth)
        # ax[0].set_title(title, fontsize=Fontsize+1)
        ax.legend([" nRMS ",
                      "nRMS average=" +str(avg_rms),
                      "nRMS median="+str(med_rms)],
                     fontsize=Labelsize, loc="best")
        ax.set_ylabel("nRMS ", fontsize=Fontsize-1)
        ax.set_ylim(rms_limits)
        ax.grid(True)
        ax.tick_params(labelsize=Labelsize)

        # if PlotPLM:
        #     ax2=ax.twinx()   # make a plot with different y-axis using second axis object
        #     ax2.plot(prof_dist[:],data_plm[:],color="blue",linewidth=Linewidths[0])
        #     ax2.set_ylabel("plm (nT)", fontsize=Lsize)
        #     ax2.legend([" powerline monitor"], fontsize=Lsize, loc="upper right")
        #     if PLimits:
        #         ax2.set_ylim(PLimits)

    if PlotType in [0, 1, 2]:
        ax = fig.add_subplot(figs[2])

        ax.add_collection(p)

        ax.set_xlim((min(site_r) - dxmed2, max(site_r) + dxmed2))
        ax.set_ylabel("depth (m asl)", fontsize=Fontsize)
        ax.yaxis.set_label_position("left")
        ax.set_xlabel(" profile distance (m)", fontsize=Fontsize)
        ax.tick_params(labelsize=Fontsize)

        ax.set_ylim((plotmax, plotmin))
        ax.invert_yaxis()
        ax.grid(True)

        if len(strng_sens) > 0:
            ax.text(0.5, 0.1, strng_sens,
                       verticalalignment="bottom", horizontalalignment="center",
                       transform=ax.transAxes,
                       fontsize=Fontsize-2)

        ax.text(0.0, -0.3, beg_strng,
                   verticalalignment="top", horizontalalignment="left",
                   transform=ax.transAxes,
                   color="black", fontsize=Fontsize)
        ax.text(0.9, -0.3, end_strng,
                   verticalalignment="top", horizontalalignment="left",
                   transform=ax.transAxes,
                   color="black", fontsize=Fontsize)

        ax = fig.add_subplot(figs[nplots-1])
        cb = fig.colorbar(p, ticks=cb_ticks , orientation="horizontal", aspect=60, pad=0.4)
        cb.set_label(r"log10($\Omega$ m)", size=Fontsize)
        cb.ax.tick_params(labelsize=Fontsize)

    for F in PlotFmt:
      matplotlib.pyplot.savefig(PlotDir+FileName+"_model"+F, dpi=400)

    if matplotlib.get_backend()!="cairo":
        matplotlib.pyplot.show()
    matplotlib.pyplot.clf()


    if PdfCatalog:
        pdf_list.append(PlotDir+FileName+".pdf")
