#!/usr/bin/env python3

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
import matplotlib.colors
import matplotlib.cm

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)


from version import versionstrg
import inverse
import util
import prep
import aesys
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

cm = 1/2.54  # centimeters in inches
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

now = datetime.now()

"""
input formats are .npz, .nc4, 'ascii'
"""

InStrng = ""
PlotStrng = " - data "+InStrng
InModDir = AEMPYX_DATA + "/Projects/Compare_systems/results/"
print("Data read from dir:  %s" % InModDir)

FileList = "search"  # "search", "read
SearchStrng = "CGG*Prior3*1000*.npz"


if "set" in FileList.lower():
    InModDir = AEMPYX_DATA + ""
    print("Data files read from dir:  %s" % InModDir)
    res_files = []

else:
    # how = ["search", SearchStrng, InModDir]
    # how = ["read", FileList, InModDir]
    res_files = util.get_data_list(how=["search", SearchStrng, InModDir],
                              out= True, sort=True)

ns = numpy.size(res_files)
if ns ==0:
    error("No files set!. Exit.")



PlotType = 1       # rms, datafit, model
# PlotType =  0      # rms, datafit, model
# PlotType =  1      #  model + rms
# PlotType =  2      #  model
# PlotType =  3      #  datafit

PlotDir = AEMPYX_DATA + "Projects/Compare_systems/plots/"
print("Plots written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)

FilesOnly = False
PlotFmt = [".pdf", ".png"]
PdfCatalog = True
if ".pdf" in PlotFmt:
    PdfCatName = "GENESIS_Type"+str(PlotType)+"_Catalog.pdf"
else:
    print(" No pdfs generated. No catalog possible!")
    PdfCatalog = False



PlotSize = [25., 5. ]

"""
Parameter for nRMS plot
"""
rms_limits  =[0., 5.]

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
    sens_pars = ["euc","size", "sqrt","max"]
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
    maxdoi = 100.
    alpha_doi = 0.0


plot_adapt = True
if plot_adapt:
   if not max_doi:
       maxdoi = 200.
else:
    plot_min = -100.0
    plot_max = 60.0
topo_use_average = False #



"""
Parameter for data fit plot
"""
Quantile95 = False
alpha_err95 = 0.2
Quantile68 = True
alpha_err68 = 0.2

DataTrans = "asinh"
XLimits = [0., 16.]
ZLimits = [0., 16.]


LogPlot = False
LogSym = False
LinThresh =100.
if LogPlot == False:
    LogSym = False
Logparams=[LogPlot, LogSym, LinThresh]

"""
General Parameter for all plots
"""
PosLatLon = True
if PosLatLon:
    EPSG=32629

ProfType = "distance"
if "dist" in ProfType.lower():
    ProfLabel = "profile distance "
    ProfScale = 1. # 0.001  # m to km
    ProfUnit  = "(m)" #
else:
	ProfLabel = "site # "
	ProfScale = 1. # 0.001  # m to km
    ProfUnit  = "(#)" #

PlotThresh =20

PosLatLon = True
if PosLatLon:
    EPSG=32629

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""

matplotlib.pyplot.style.use("seaborn-paper") # ("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
# matplotlib.rcParams["text.usetex"] = True

Fontsize = 8
Labelsize = Fontsize
Titlesize = 10
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidth= 1.5

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


ns = numpy.size(res_files)

ifl = 0
pdf_list = []
for file in res_files:

    nplots = 0

    FileName, filext0 = os.path.splitext(file)

    title=FileName.replace("_"," ").replace("-"," ")



    """
    numpy.savez_compressed(
        file=Fileout,
        fl_data=file,
        fl_name=fl_name,
        header=Header,
        ctrl = Ctrl,
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
        site_nump=site_nump,
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
    site_rms = tmp["site_nrms"]
    jac          = tmp["site_jacd"]
    
    d_active    = tmp["dat_act"]
    site_dobs   = tmp["site_dobs"]
    site_dcal   = tmp["site_dcal"]
    site_derr   = tmp["site_derr"]


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
        print(site_sens[isite,:])
    if topo_use_average:
        site_tref = numpy.mean(site_dem)
    else:
        site_tref = 0
    site_topo = site_dem - site_tref
    max_topo = numpy.amax(site_topo)

    beg_pos = [site_x[0],  site_y[0],  site_topo[0] ]
    end_pos = [site_x[-1], site_y[-1], site_topo[-1]]
    beg_strng  = f"Start:\nX: {beg_pos[0]:.0f} \nY: {beg_pos[1]:.0f} (EPSG="+str(EPSG)+")"
    end_strng  = f"End:\nX: {end_pos[0]:.0f} \nY: {end_pos[1]:.0f}"

    if PosLatLon:
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

        size_lay = numpy.repeat(thk.T, sites, axis=0)
        site_val = numpy.log10(site_model[:,:])
        site_val[:, -1] = numpy.nan
        

        site_sens = numpy.log10(site_sens[:,:])


        alpha = numpy.ones_like(site_val)
        site_dpth = numpy.ones_like(site_topo)
        site_node = numpy.ones((sites, numpy.size(z0)))

        for nmod  in  numpy.arange(sites):
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

    if PlotType in [0, 1]:
        #
        avg_rms = round(numpy.nanmean(site_rms), 2)
        med_rms = round(numpy.nanmedian(site_rms),2)
        med = med_rms*numpy.ones_like(site_rms)
        avg = avg_rms*numpy.ones_like(site_rms)


    if PlotType in [0, 3]:
  
       

        X_obs = site_dobs[:, 0:11]
        Z_obs = site_dobs[:, 11:22]
        X_cal = site_dcal[:, 0:11]
        Z_cal = site_dcal[:, 11:22]
        dunit = "(fT)"

        # 98 % quantiles
        X_err68 = site_derr[:, 0:11]
        Z_err68 = site_derr[:, 11:22]

        X_err95 = site_derr[:, 0:11]*2.
        Z_err95 = site_derr[:, 11:22]*2.

        if Quantile95:
            QZerr =Z_err95
            QXerr =Z_err95
            alpha_err=alpha_err95
        else:
            QZerr =Z_err68
            QXerr =Z_err68
            alpha_err=alpha_err68

        if "asinh"in DataTrans.lower():

            LogPlot = False
            dunit = "(-)"
            Xc_obs = numpy.arcsinh(X_obs)
            Xplus = numpy.arcsinh(X_obs+QXerr)
            Xmins = numpy.arcsinh(X_obs-QXerr)
            Xc_cal = numpy.arcsinh(X_cal)

            Zc_obs = numpy.arcsinh(Z_obs)
            Zplus = numpy.arcsinh(Z_obs+QZerr)
            Zmins = numpy.arcsinh(Z_obs-QZerr)
            Zc_cal = numpy.arcsinh(Z_cal)



        X_min, X_max = numpy.nanmin(Xc_obs), numpy.nanmax(Xc_obs)
        Z_min, Z_max = numpy.nanmin(Zc_obs), numpy.nanmax(Zc_obs)

        print("Xmin,max = "+str(X_min)+",  "+str(X_max))
        print("Zmin,max = "+str(Z_min)+",  "+str(Z_max))


    print(PlotType)
    # fig =  matplotlib.pyplot.figure(figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm))

    if PlotType == 0:
        nplots = 4
        fig, ax = matplotlib.pyplot.subplots(4, 1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [1, 2, 2, 4]})

    if PlotType == 1:
        nplots = 2
        fig, ax = matplotlib.pyplot.subplots(nplots, 1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [1,3]})
    if PlotType == 2:
        nplots = 1
        fig, ax = matplotlib.pyplot.subplots(nplots, 1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm),
                                          sharex=True,
                                          gridspec_kw={"height_ratios": [3]})
    if PlotType == 3:
        nplots = 1
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
        ax[ii].legend(["nRMS ",
                      "nRMS average=" +str(avg_rms),
                      "nRMS median="+str(med_rms)],
                     fontsize=Labelsize, loc="best")
        ax[ii].set_ylabel("nRMS ", fontsize=Fontsize-1)
        ax[ii].set_ylim(rms_limits)
        ax[ii].grid(True)
        ax[ii].tick_params(labelsize=Labelsize)

    if PlotType in [0, 3]:

        ii = ii+1
        xtit = "Z"
        if "asinh"in DataTrans.lower():
            xtit = "Z (asinh)"

        for ichan in numpy.arange(ncols):


            ax[ii].fill_between(site_r, Zmins[:,ichan], Zplus[:,ichan],
                                color=Colors[ichan],alpha=alpha_err, label ="_nolegend_")
            ax[ii].plot(
                  site_r, Z_obs[:, ichan],
                  color=Colors[ichan], linewidth=Linewidth, linestyle=":")
            ax[ii].plot(
                site_r, Zc_cal[:,ichan],
                color=Colors[ichan], linewidth=Linewidth, linestyle="-")

      # ax[ii].errorbar(
        #     site_r, Z_obs[:, 0], yerr=Z_err95[:, 0],
        #     color=Colors[ichan],linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, Z_obs[:, 1], yerr=Z_err95[:, 1],
        #     color="g",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, Z_obs[:, 2], yerr=Z_err95[:, 2],
        #     color="b",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, Z_obs[:, 3], yerr=Z_err95[:, 3],
        #     color="k",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)

        ax[ii].set_title("Z (asinh)", fontsize=Fontsize, y=1.0, pad=-(Fontsize+2))
        ax[ii].set_ylim([Z_min, Z_max])

        ax[ii].set_ylabel(dunit, fontsize=Labelsize)
        ax[ii].tick_params(labelsize=Labelsize)
        ax[ii].grid(True)
        leg1 = ax[ii].legend(
            [r"0.009 ms", r"0.026 ms", r"0.052 ms", r"0.095 ms",
              r"0.156 ms", r"0.243 ms", r"0.365 ms", r"0.547 ms",
              r"0.833 ms", r"1.259 ms", r"1.858 ms"],
            fontsize=Labelsize-3, loc="best",ncol=3)
        leg1.set_title("Window (center)", prop={"size":Labelsize-1})

        if ii == nplots:
            ax[ii].set_xlabel(ProfLabel+ProfUnit, fontsize=Labelsize)
        else:
            ax[ii].tick_params(labelbottom=False)

        ii = ii+1
        xtit = "X"
        if "asinh"in DataTrans.lower():
            xtit = "X (asinh)"

        for ichan in numpy.arange(ncols):


            ax[ii].fill_between(site_r, Xmins[:,ichan], Xplus[:,ichan],
                                color=Colors[ichan],alpha=alpha_err, label ="_nolegend_")

            ax[ii].plot(
                  site_r, Xc_obs[:, ichan],
                  color=Colors[ichan], linewidth=Linewidth, linestyle=":")
            ax[ii].plot(
                site_r, Xc_cal[:,ichan],
                color=Colors[ichan], linewidth=Linewidth, linestyle="-")

        # ax[ii].errorbar(
        #     site_r, X_obs[:, 0], yerr=X_err95[:, 0],
        #     color=Colors[ichan],linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, X_obs[:, 1], yerr=X_err95[:, 1],
        #     color="g",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, X_obs[:, 2], yerr=X_err95[:, 2],
        #     color="b",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)
        # ax[ii].errorbar(
        #     site_r, X_obs[:, 3], yerr=X_err95[:, 3],
        #     color="k",linewidth=Linewidth, linestyle="-", elinewidth=Linewidth-1)

        ax[ii].set_title(xtit, fontsize=Fontsize, y=1.0, pad=-(Fontsize+2))
        ax[ii].set_ylim([Z_min, Z_max])


        ax[ii].set_ylabel(dunit, fontsize=Labelsize)

        ax[ii].tick_params(labelsize=Labelsize)
        ax[ii].grid(True)
        leg1 = ax[ii].legend(
            [r"0.009 ms", r"0.026 ms", r"0.052 ms", r"0.095 ms",
              r"0.156 ms", r"0.243 ms", r"0.365 ms", r"0.547 ms",
              r"0.833 ms", r"1.259 ms", r"1.858 ms"],
            fontsize=Labelsize-3, loc="best",ncol=3)
        leg1.set_title("Window (center)", prop={"size":Labelsize-1})

        if ii == nplots:
            ax[ii].set_xlabel(ProfLabel+ProfUnit, fontsize=Labelsize)
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
        axii.set_xlabel(ProfLabel+ProfUnit, fontsize=Fontsize)
        axii.tick_params(labelsize=Fontsize)

        axii.set_ylim((plotmax, plotmin))
        axii.invert_yaxis()
        axii.grid(True)

        if len(strng_sens) > 0:
            axii.text(0.5, 0.1, strng_sens,
                        verticalalignment="bottom", horizontalalignment="center",
                        transform=axii.transAxes,
                        fontsize=Fontsize-2)

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
                          ticks=cb_ticks ,
                          orientation="horizontal", aspect=40, pad=0.25, shrink=0.65)
        cb.set_label(r"log10($\Omega$ m)", size=Fontsize)
        cb.ax.tick_params(labelsize=Fontsize)

    for F in PlotFmt:
      matplotlib.pyplot.savefig(PlotDir+FileName+"_PType"+str(PlotType)+F, dpi=400)

    if matplotlib.get_backend()!="cairo":
        matplotlib.pyplot.show()
    matplotlib.pyplot.clf()

    if PdfCatalog:
        pdf_list.append(PlotDir+FileName+"_PType"+str(PlotType)+".pdf")


if PdfCatalog:
    viz.make_pdf_catalog(PdfList=pdf_list, FileName=PdfCatName)
