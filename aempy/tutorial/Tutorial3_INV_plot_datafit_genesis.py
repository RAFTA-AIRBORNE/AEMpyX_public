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


"""
Show several 1d block models as (stitched) section.

"""
import os
import sys
from sys import exit as error
from datetime import datetime
import warnings

import numpy

# import matplotlib.collections
# import matplotlib.patches
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import matplotlib.cm

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
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

AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
cm = 1/2.54  # centimeters in inches

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = True
now = datetime.now()

"""
input formats are "npz","nc4","asc"
"""
InFileFmt = ".npz"
InModDir = AEMPYX_DATA + "/Aug2023/results/"
print("Data read from dir:  %s" % InModDir)

FileList = "search"  # "search", "read"
SearchStrng = "CGG*results.npz"

if "search" in FileList.lower():

    print("Search flightline ID string: %s " % SearchStrng)
    res_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InModDir, fullpath=False)
    res_files = sorted(res_files)

if "set" in FileList.lower():
    res_files =[]

print (res_files)

PlotType =  1

# PlotType =  0      # rms, datafit, model
# PlotType =  1      #  model + rms
# PlotType =  2      #  model
# PlotType =  3      #  datafit

# PlotDir = AEMPYX_DATA + "/Projects/Compare/plots/C36/"
PlotDir = AEMPYX_DATA + "/Aug2023/plots/"
print("Plots written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)


FilesOnly = True
PlotFmt = [".pdf"]
PDFCatalog = True
PDFStrng = ".pdf"
if ".pdf" in PlotFmt:
    PDFCatName = PlotDir+"GENESIS_Type"+str(PlotType)+PDFStrng
else:
    print(" No pdfs generated. No catalog possible!")
    PDFCatalog = False



PlotSize = [25., 5. ]

"""
Parameter for data fit plot
"""
Quantile95 = False
alpha_err95 = 0.2
Quantile68 = True
alpha_err68 = 0.2
# QLimits = [0., 2000.]
# ILimits = [0., 2000.]

QLimits = []
ILimits = []
"""
Parameter for nRMS plot
"""
rms_limits  =[0., 4.]


"""
Parameter for model plot
"""
min_lrho =  1.
max_lrho =  4.
cl = [min_lrho, max_lrho]
cb_ticks = [-1, 0, 1, 2, 3, 4]


blank = 10

low_sens = True
if low_sens:        
    lowsens = -2.
    # sens_pars = ["euc","size", "sqrt","max", "log"]
    sens_pars = ["euc","sqrt","max", "log"]
    
    strng_sens = " sensitivity scaling = "+str(sens_pars)+", thresh ="+str(lowsens)
    alpha_sens = 0.05
else:
    lowsens = -6.
    strng_sens = ""
    alpha_sens = 1.


high_rms = True
if high_rms:
    highrms = 3.
    alpha_rms= 0.05

high_err = False
if high_err:
    higherr = 0.5
    alpha_err95 = 0.3


max_doi = False
if max_doi:
    maxdoi = 150.
    alpha_doi = 0.0


plot_adapt = False
if plot_adapt:
   if not max_doi:
       maxdoi = 150.
else:
    plot_min = -200.0
    plot_max = 60.0


topo_use_average = False #


"""
General Parameter for all plots
"""
poslatlon = True
if poslatlon:
    EPSG=32629

ProfType = "distance"
if "dist" in ProfType.lower():
    ProfLabel = "distance (m)"
    ProfScale = 1. # 0.001  # m to km
else:
    ProfLabel= "site #"
    ProfScale = 1. # 0.001  # m to km

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""

matplotlib.pyplot.style.use("seaborn-v0_8-paper") # ("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
# matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.transparent"] = True
matplotlib.rcParams["savefig.bbox"] = "tight" 

Fontsize = 8
Labelsize = Fontsize
Titlesize = 10
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidth= 1.5
Markersize = 4

ncols = 4
Colors = matplotlib.pyplot.cm.jet(numpy.linspace(0,1,ncols))
Grey = 0.7

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


ns = numpy.size(res_files)

ifl = 0
pdf_list = []
for file in res_files:
    
    nplots = 0

    FileName, filext0 = os.path.splitext(file)

    title=FileName.replace("_"," ")


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

    m_active    = tmp["mod_act"]
    mod_ref     = tmp["mod_ref"]
    site_model  = tmp["site_modl"]
    site_error  = tmp["site_merr"]
    site_sens   = tmp["site_sens"]
    jac          = tmp["site_jacd"]



    d_active    = tmp["dat_act"]
    site_dobs   = tmp["site_dobs"]
    site_dcal   = tmp["site_dcal"]
    site_derr    = tmp["site_derr"]

    nlyr = inverse.get_nlyr(mod_ref)
    nsite, ndata = numpy.shape(site_dobs)

    site_rms = tmp["site_nrms"]

    site_x = tmp["site_x"] * ProfScale
    site_y = tmp["site_y"] * ProfScale

    site_alt = tmp["site_alt"]
    site_gps = tmp["site_gps"]
    site_dem = tmp["site_dem"]
    
    
    models = numpy.shape(site_model)
    sites = models[0]
    param = models[1]
    
    nlyr = inverse.get_nlyr(mod_ref)
    dz = mod_ref[6*nlyr:7*nlyr-1] 
    dz = numpy.append(dz, 1.*dz[-1]) 
    z0 = numpy.hstack((0.0, numpy.cumsum(dz)))        
    zm = 0.5 * (z0[0:nlyr] + z0[1 : nlyr + 1])
    thk = dz
    thk[-1] = 1.e6

    for isite in numpy.arange(sites):
        ndata = numpy.sum(d_active[isite,:])
        jaci = jac[isite,:].reshape(ndata,param)
        sens = inverse.calc_sensitivity(Jac=jaci, 
                                        UseSigma=False, 
                                        Type =sens_pars[0])
        site_sens[isite,:] = inverse.transform_sensitivity(S=sens, 
                                                           V=thk, 
                                                           Transform=sens_pars[1:])
        # if isite  in [ 3, 10, 70]: print(site_sens[isite,:])

    if topo_use_average:
        site_tref = numpy.mean(site_dem)
    else:
        site_tref = 0
    site_topo = site_dem - site_tref
    max_topo = numpy.amax(site_topo)
    
    avg_rms = round(numpy.nanmean(site_rms), 2)
    med_rms = round(numpy.nanmedian(site_rms),2)
    med = med_rms*numpy.ones_like(site_rms)
    avg = avg_rms*numpy.ones_like(site_rms)

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

    if "dist" in ProfType.lower():
        site_x = site_x - site_x[0]
        site_y = site_y - site_y[0]
        site_r = numpy.sqrt(numpy.power(site_x, 2.0) + numpy.power(site_y, 2.0))
    else:
        site_r = numpy.arange(numpy.shape(site_x)[0])

    if PlotType in [0, 1, 2]:

        dxmed2 = numpy.median(numpy.diff(site_r)) / 2.0

        thk = mod_ref[6 * nlyr : 7 * nlyr - 1].reshape(1,-1)
        thk = numpy.vstack((thk.T, thk[-1,0]))
        z0 = numpy.hstack((0.0, numpy.cumsum(thk)))
        zm = 0.5 * (z0[0:nlyr] + z0[1 : nlyr + 1])
        size_lay = numpy.repeat(thk.T, sites, axis=0)
        site_val = numpy.log10(site_model[:,:])
        
        # site_sens = numpy.log10(site_sens[:,:])


        alpha = numpy.ones_like(site_val)
        site_dpth = numpy.ones_like(site_topo)
        site_node = numpy.ones((sites, numpy.size(z0)))

        for nmod  in  numpy.arange(sites):
            # if nmod  in [ 3, 10, 70]: print(" nmod ", site_sens[nmod,:])
            site_val[nmod, -1] = numpy.nan
            zmp = site_topo[nmod] - zm
            z0p = site_topo[nmod] - z0
            site_dpth = numpy.amin(z0p)
            site_node [nmod, :] = z0p

            if high_rms:
                if site_rms[nmod] > highrms:
                    alpha[nmod, :] = alpha_rms

            if low_sens:
                 for il in numpy.arange(nlyr):
                     if site_sens[nmod, il] < lowsens:
                             alpha[nmod, il] = alpha_sens

            if high_err:
                for il in numpy.arange(nlyr):
                    if numpy.abs(site_derr[nmod,il]) > higherr:
                            alpha[nmod, il] = alpha_err95

            if max_doi:
                for il in numpy.arange(nlyr):
                    if z0p[0]-z0p[il] > maxdoi:
                            alpha[nmod, il] = alpha_doi


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

        I_obs = site_dobs[:, 0:4]
        Q_obs = site_dobs[:, 4:8]
        I_cal = site_dcal[:, 0:4]
        Q_cal = site_dcal[:, 4:8]
        #  quantiles
        I_err68 = site_derr[:, 0:4]
        Q_err68 = site_derr[:, 4:8]
        I_err95 = site_derr[:, 0:4]*2.
        Q_err95 = site_derr[:, 4:8]*2.


        I_min, I_max = numpy.nanmin(I_obs), numpy.nanmax(I_obs)
        Q_min, Q_max = numpy.nanmin(Q_obs), numpy.nanmax(Q_obs)

        # print("Imin,max = "+str(I_min)+",  "+str(I_max))
        # print("Qmin,max = "+str(Q_min)+",  "+str(Q_max))



    print(PlotType)
    # fig =  matplotlib.pyplot.figure(figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm))

    if PlotType == 0:
        nplots = 4
        fig, ax = matplotlib.pyplot.subplots(4, 1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [2, 2, 2, 6]})

    if PlotType == 1:
        nplots = 2
        fig, ax = matplotlib.pyplot.subplots(nplots, 1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [1,2]})
    if PlotType == 2:
        nplots = 1
        fig, ax = matplotlib.pyplot.subplots(nplots, 1,
                                          figsize=(PlotSize[0]*cm, 1.5*nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [3]})
    if PlotType == 3:
        nplots = 2
        fig, ax = matplotlib.pyplot.subplots(nplots, 1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [1, 1]})

    fig.suptitle(title, x=0.5, y=0.9, fontsize=Titlesize, va="baseline") #, fontweight="bold")

    ii = - 1
    if PlotType in [0, 1]:

        ii = ii+1

        ax[ii].plot(site_r[:-1], site_rms[:-1], "r", linewidth=Linewidth)
        ax[ii].plot(site_r[:-1], avg[:-1], "b:", linewidth=Linewidth)
        ax[ii].plot(site_r[:-1], med[:-1], "g:", linewidth=Linewidth)
        # ax[ii][0].set_title(title, fontsize=Fontsize+1)
        ax[ii].legend([" nRMS ",
                      "nRMS average=" +str(avg_rms),
                      "nRMS median="+str(med_rms)],
                     fontsize=Labelsize, loc="best")
        ax[ii].set_ylabel("nRMS " , fontsize=Fontsize-1)
        ax[ii].set_ylim(rms_limits)
        ax[ii].grid(True)
        ax[ii].tick_params(labelsize=Labelsize)

        # if PlotPLM:
        #     ax[ii]=ax.twinx()   # make a plot with different y-axis using second axis object
        #     ax[ii].plot(prof_dist[:],data_plm[:],color="blue",linewidth=Linewidths[0])
        #     ax[ii].set_ylabel("plm (nT)", fontsize=Labelsize)
        #     ax[ii].legend([" powerline monitor"], fontsize=Labelsize, loc="upper right")
        #     if PLimits:


    if PlotType in [0, 3]:

        ii = ii+1
        if Quantile95:
            ax[ii].fill_between(
                site_r, Q_obs[:, 0]-Q_err95[:, 0], Q_obs[:, 0]+Q_err95[:, 0],
                color="r",alpha=alpha_err95, label ="_nolegend_")
            ax[ii].fill_between(
                site_r, Q_obs[:, 1]-Q_err95[:, 1], Q_obs[:, 1]+Q_err95[:, 1],
                color="g",alpha=alpha_err95, label ="_nolegend_" )
            ax[ii].fill_between(
              site_r, Q_obs[:, 2]-Q_err95[:, 2], Q_obs[:, 2]+Q_err95[:, 2],
              color="b",alpha=alpha_err95, label ="_nolegend_")
            ax[ii].fill_between(
                  site_r, Q_obs[:, 3]-Q_err95[:, 3], Q_obs[:, 3]+Q_err95[:, 3],
                  color="m",alpha=alpha_err95,label ="_nolegend_")
        if Quantile68:
            ax[ii].fill_between(
                site_r, Q_obs[:, 0]-Q_err68[:, 0], Q_obs[:, 0]+Q_err68[:, 0],
                color="r",alpha=alpha_err68,label ="_nolegend_")
            ax[ii].fill_between(
                site_r, Q_obs[:, 1]-Q_err68[:, 1], Q_obs[:, 1]+Q_err68[:, 1],
                color="g",alpha=alpha_err68,label ="_nolegend_")
            ax[ii].fill_between(
              site_r, Q_obs[:, 2]-Q_err68[:, 2], Q_obs[:, 2]+Q_err68[:, 2],
              color="b",alpha=alpha_err68,label ="_nolegend_")
            ax[ii].fill_between(
                  site_r, Q_obs[:, 3]-Q_err68[:, 3], Q_obs[:, 3]+Q_err68[:, 3],
                  color="m",alpha=alpha_err68,label ="_nolegend_")

        ax[ii].plot(
             site_r, Q_obs[:, 0],
             color="r", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, Q_obs[:, 1],
            color="g", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, Q_obs[:, 2],
            color="b", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, Q_obs[:, 3],
            color="m", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, Q_cal[:, 0],
            color="r", linewidth=Linewidth, linestyle="-")
        ax[ii].plot(
            site_r, Q_cal[:, 1],
            color="g", linewidth=Linewidth, linestyle="-")
        ax[ii].plot(
            site_r, Q_cal[:, 2],
            color="b", linewidth=Linewidth, linestyle="-")
        ax[ii].plot(
            site_r, Q_cal[:, 3],
            color="m", linewidth=Linewidth, linestyle="-")

        # ax[ii].errorbar(
        #     site_r, Q_obs[:, 0], yerr=Q_err95[:, 0],
        #     color="r",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, Q_obs[:, 1], yerr=Q_err95[:, 1],
        #     color="g",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, Q_obs[:, 2], yerr=Q_err95[:, 2],
        #     color="b",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, Q_obs[:, 3], yerr=Q_err95[:, 3],
        #     color="k",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)

        ax[ii].set_title("Quadrature", fontsize=Fontsize, y=1.0, pad=-(Fontsize+2))
        ax[ii].set_ylim([Q_min, Q_max])
        if numpy.size(QLimits)>0:
            ax[ii].set_ylim(QLimits)

        ax[ii].set_ylabel("(ppm)", fontsize=Labelsize)

        ax[ii].tick_params(labelsize=Labelsize)
        ax[ii].grid(True)

        legend =ax[ii].legend(
            [" 0.9 kHz", "3 kHz", "12 kHz", "24.5 kHz"],
            fontsize=Labelsize-2, loc="best", ncol=2)
        legend.set_title("Frequency", prop={"size":Labelsize})

        # if LogPlot:
        #     if LogSym:
        #         ax[ii].set_yscale("symlog", linthresh=LinThresh)
        #     else:
        #         ax[ii].set_yscale("log")
        # else:
        #     ax[ii].set_yscale("linear")
        if ii == nplots:
            ax[ii].set_xlabel(ProfLabel, fontsize=Labelsize)
        else:
            ax[ii].tick_params(labelbottom=False)

        ii = ii+1

        if Quantile95:
            ax[ii].fill_between(
                site_r, I_obs[:, 0]-I_err95[:, 0], I_obs[:, 0]+I_err95[:, 0],
                color="r",alpha=alpha_err95,label ="_nolegend_")
            ax[ii].fill_between(
                site_r, I_obs[:, 1]-I_err95[:, 1], I_obs[:, 1]+I_err95[:, 1],
                color="g",alpha=alpha_err95,label ="_nolegend_")
            ax[ii].fill_between(
              site_r, I_obs[:, 2]-I_err95[:, 2], I_obs[:, 2]+I_err95[:, 2],
              color="b",alpha=alpha_err95,label ="_nolegend_")
            ax[ii].fill_between(
                  site_r, I_obs[:, 3]-I_err95[:, 3], I_obs[:, 3]+I_err95[:, 3],
                  color="m",alpha=alpha_err95,label ="_nolegend_")
        if Quantile68:
            ax[ii].fill_between(
                site_r, I_obs[:, 0]-I_err68[:, 0], I_obs[:, 0]+I_err68[:, 0],
                color="r",alpha=alpha_err68,label ="_nolegend_")
            ax[ii].fill_between(
                site_r, I_obs[:, 1]-I_err68[:, 1], I_obs[:, 1]+I_err68[:, 1],
                color="g",alpha=alpha_err68,label ="_nolegend_")
            ax[ii].fill_between(
              site_r, I_obs[:, 2]-I_err68[:, 2], I_obs[:, 2]+I_err68[:, 2],
              color="b",alpha=alpha_err68,label ="_nolegend_")
            ax[ii].fill_between(
                  site_r, I_obs[:, 3]-I_err68[:, 3], I_obs[:, 3]+I_err68[:, 3],
                  color="m",alpha=alpha_err68,label ="_nolegend_")


        ax[ii].plot(
             site_r, I_obs[:, 0],
             color="r", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, I_obs[:, 1],
            color="g", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, I_obs[:, 2],
            color="b", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, I_obs[:, 3],
            color="m", linewidth=Linewidth, linestyle=":")
        ax[ii].plot(
            site_r, I_cal[:, 0],
            color="r", linewidth=Linewidth, linestyle="-")
        ax[ii].plot(
            site_r, I_cal[:, 1],
            color="g", linewidth=Linewidth, linestyle="-")
        ax[ii].plot(
            site_r, I_cal[:, 2],
            color="b", linewidth=Linewidth, linestyle="-")
        ax[ii].plot(
            site_r, I_cal[:, 3],
            color="m", linewidth=Linewidth, linestyle="-")


        # ax[ii].errorbar(
        #     site_r, I_obs[:, 0], yerr=I_err95[:,0],
        #     color="r",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, I_obs[:, 1], yerr=I_err95[:,1],
        #     color="g",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, I_obs[:, 2], yerr=I_err95[:,2],
        #     color="b",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, I_obs[:, 3], yerr=I_err95[:,3],
        #     color="k",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)


        ax[ii].set_title("In-Phase", fontsize=Fontsize, y=1.0, pad=-(Fontsize+2))
        ax[ii].set_ylim([I_min, I_max])
        if numpy.size(ILimits)>0:
            ax[ii].set_ylim(ILimits)
        ax[ii].set_ylabel("(ppm)", fontsize=Fontsize)

        if ii == nplots:
            ax[ii].set_xlabel(ProfLabel, fontsize=Labelsize)
        else:
            ax[ii].tick_params(labelbottom=False)



        ax[ii].grid(True)

        legend =ax[ii].legend(
            [" 0.9 kHz", "3 kHz", "12 kHz", "24.5 kHz"],
            fontsize=Labelsize-2, loc="best", ncol=2)
        legend.set_title("Frequency", prop={"size":Labelsize})

        # if LogPlot:
        #     if LogSym:
        #         ax[ii].set_yscale("symlog", linthresh=LinThresh)
        #     else:
        #         ax[ii].set_yscale("log")
        # else:
        #     ax[ii].set_yscale("linear")

        if ii == nplots:
            ax[ii].set_xlabel(ProfLabel, fontsize=Labelsize)
        else:
            ax[ii].tick_params(labelbottom=False)

    if PlotType in [0, 1, 2]:

        ii = nplots-1

        if ii == 0:
            axii = ax
        else:
            axii = ax[ii]

        for isit in numpy.arange(sites):
            x = [site_r[isit] - dxmed2, site_r[isit] + dxmed2]
            y = site_node[isit,:]
            c = site_val[isit,:].reshape(-1,1)
            a = alpha[isit,:]
            axii.pcolormesh(x, y, c,
                            alpha = a,
                            shading='flat',
                            vmin=cl[0], vmax=cl[1],
                            cmap=mycmap)


        axii.set_xlim((min(site_r) - dxmed2, max(site_r) + dxmed2))
        axii.set_ylabel("depth (m asl)", fontsize=Fontsize)
        axii.yaxis.set_label_position("left")
        axii.set_xlabel(ProfLabel, fontsize=Fontsize)
        axii.tick_params(labelsize=Fontsize)
        axii.set_ylim((plotmax, plotmin))
        axii.invert_yaxis()
        axii.grid(True)
        
        # if PlotType not in [0, 1]:
        #     strng_rms = "nRMS: avg=" +str(avg_rms)+ "  med="+str(med_rms)
        #     axii.text(0.05, 0.8, strng_rms,
        #                verticalalignment="bottom", horizontalalignment="left",
        #                transform=axii.transAxes,
        #                fontsize=Fontsize-2)

        if len(strng_sens) > 0:
            axii.text(0.05, 0.1, strng_sens,
                       verticalalignment="bottom", horizontalalignment="left",
                       transform=axii.transAxes,
                       fontsize=Fontsize-2)
        axii.text(0.95, 0.1, "nRMS average=" +str(avg_rms)+", median="+str(med_rms),
                  verticalalignment="bottom", horizontalalignment="right",
                    transform=axii.transAxes,
                    color="black", fontsize=Fontsize-2)
        axii.text(0.0, -0.3, beg_strng,
                   verticalalignment="top", horizontalalignment="left",
                   transform=axii.transAxes,
                   color="black", fontsize=Fontsize)
        axii.text(0.9, -0.3, end_strng,
                   verticalalignment="top", horizontalalignment="left",
                   transform=axii.transAxes,
                   color="black", fontsize=Fontsize)

        norm = matplotlib.colors.Normalize(vmin=cl[0], vmax=cl[1], clip=False)
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=mycmap), 
                          ticks=cb_ticks , ax=matplotlib.pyplot.gca(),
                          orientation="horizontal", aspect=40, pad=0.25, shrink=0.65)
        cb.set_label(r"log10($\Omega$ m)", size=Fontsize)
        cb.ax.tick_params(labelsize=Fontsize)

    for F in PlotFmt:
      matplotlib.pyplot.savefig(PlotDir+FileName+"_PType"+str(PlotType)+F)

    if matplotlib.get_backend()!="cairo":
        matplotlib.pyplot.show()
        
    matplotlib.pyplot.clf()

    if PDFCatalog:
        pdf_list.append(PlotDir+FileName+"_PType"+str(PlotType)+".pdf")


if PDFCatalog:
    viz.make_pdf_catalog(PDFList=pdf_list, FileName=PDFCatName)
