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
#       format_version: '1.5'
#       jupytext_version: 1.15.2
# ---

"""
Created on Tue Aug  3 17:03:39 2021

@author: vrath
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings

import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.axis
import mpl_toolkits.axes_grid1


import scipy.interpolate
import scipy.spatial
import skgstat
import shapely
# import rasterio
# from rasterio import features
# import affine

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys
import viz
import inverse



AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
script = "VIZ_models_area.py"
# script = __file__  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")

# for tutorial only...
AEMPYX_DATA = AEMPYX_ROOT+"/data/"

now = datetime.now()
cm = 1/2.54
OutInfo = True

# """
# System related settings.
# Data transformation is now allowed with three possible options:
# DataTrans   = 0           raw data
#             = 1           natural log of data
#             = 2           asinh transformation
# An error model is applied for the raw data, which is
# mixed additive/multiplicative. in case of data transformation,
# errors are also transformed.
# """
# AEM_system = "genesis"
AEM_system = "aem05"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1


InModDir = AEMPYX_DATA+"/aem05_limerick/merged/"

# FileList = "search"
# SearchStrng = "A9*k3*.npz"
# print("Searchstring: %s \n" % SearchStrng)

FileList = "set"
ListName = ""
SearchStrng = ""

if "set" in FileList.lower():
    # mod_files = [InModDir+"LimShale_proc_dec5_mean_merged_results.npz"]
    mod_files = [InModDir+"LimShale_proc_k3_dec5_mean_merged_results.npz"]
    # mod_files = [InModDir+"LimShale_proc_k2_dec5_mean_merged_results.npz"]
if "read" in FileList.lower():
    print("File names read from : "+ListName)
    how = ["read", ListName, InModDir]
    mod_files = util.get_data_list(how=how,
                              out= True, sort=True)

    mod_files = numpy.loadtxt("A9-7.dat", dtype=str)

if "search" in FileList.lower():
    print("Searchstring is : "+SearchStrng)
    how = ["search", SearchStrng, InModDir]
    mod_files = util.get_data_list(how=how,
                              out=True,
                              fullpath=True,
                              sort=True)

ns = numpy.size(mod_files)
if ns ==0:
    error("No files set!. Exit.")

MergeModels = False
ModelMergeFile = InModDir+"MUN_k3_data_merged.npz"

"""
Output formats are "npz","nc4","ascii"
"""
PlotFmt = [".pdf", ".png"] #".png", ".pdf",]

PDFCatalog = True
# PDFCName = "LimShale_proc_dec5_mean.pdf"
# PDFCName = "LimShale_proc_k3_dec5_mean.pdf"
PDFCName = "LimShale_proc_K2_dec5_mean.pdf"
if ".pdf" in PlotFmt:
    pass
else:
    error(" No pdfs generated. No catalog possible!")
    PDFCatalog = False

PlotDir = InModDir+"/plots/"
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)
print("Plots written to dir: %s " % PlotDir)
PlotName = "LimShale"
print("Plot filname: %s " % PlotName)


ImageType = "image"
ImageType = "contour"
# ImageType = "scatter"

Layers = [1, 5, 10, 15, 20, 25]
Prop = "rho"
Unit = r"log10 $\Omega$m"
Limits = [-1., 4.]
# Steps  = [-1., -0.5, 0., 0.5, 1., 1.5, 2., 2.5,  3., 3.5, 4.,]
Steps = numpy.arange(-1., 4.01, 0.2)
print("\nLayer parameters:")
LayList = []
for il in numpy.arange(len(Layers)):
    label = "Layer"+str(Layers[il])
    LayList.append([label, Layers[il], Prop, Unit, Limits, Steps])
    print(LayList[il])




"""
Kernel functions for RBF:
    The radial basis function, based on the radius, r,
    given by the norm (default is Euclidean distance); the default is ‘multiquadric’:
        ‘linear’ : -r
        ‘thin_plate_spline’ : r**2 * log(r)
        ‘cubic’ : r**3
        ‘quintic’ : -r**5

If a callable, then it must take 2 arguments (self, r). The epsilon parameter
will be available as self.epsilon. Other keyword arguments passed
in will be available as well.


Methods for griddata:
        'nearest'       data point closest to the point of interpolation
        'linear'        tessellate the input point set to N-D simplices
                        and interpolate linearly on each simplex
        'cubic'         return the value determined from a piecewise cubic,
                        continuously differentiable (C1), and approximately
                        curvature-minimizing polynomial surface.
"""
if ("image" in ImageType.lower()) or ("contour"in ImageType.lower()):
    step = 1

    InterpMethod = ["griddata","linear"]
    # InterpMethod = ["griddata", "cubic"]
    # InterpMethod = ["rbf", "linear", 0.0]
    # InterpMethod = ["rbf", "thin_plate_spline", 0.0]
    # InterpMethod = ["rbf", "cubic", 0.01]

    # InterpMethod = ["krig", "linear", 0.5, 340.]
    S = 500.
    numIndexes = [121, 141]
    smooth = 0.
    Levels = []
    MaskPoly = False
    MaskDist = True

    if MaskDist:
        DistMask = 100.

    if MaskPoly:
        PolyDir = AEMPYX_DATA+"/Blocks/polygons/"
        PolyFiles = [PolyDir+"A9_2019_utm.npz"]
        Polygon= numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]

if ("scatter" in ImageType.lower()):
    Decim=25
    step = min(1,abs(Decim))
    Markersize = 5
    Pixelsize = False
    if Pixelsize:
        step=10
        sfac = 10.



# XYFact = 1.
# XYUnits = "(m)"
# xformatter = matplotlib.ticker.FormatStrFormatter("%7f")
# yformatter = matplotlib.ticker.FormatStrFormatter("%6f")

XYUnits = "(km)"
XYFact = 0.001
xformatter = matplotlib.ticker.FormatStrFormatter("%.2f")
yformatter = matplotlib.ticker.FormatStrFormatter("%.2f")

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""
FilesOnly = False
matplotlib.pyplot.style.use("seaborn-v0_8-paper")
matplotlib.rcParams["text.usetex"] = False
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
matplotlib.rcParams["savefig.bbox"]= "tight"

Fontsize = 7
Labelsize = Fontsize
Titlesize = 8
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidths= [0.5]
FigWidth = 16.

"""
Determine colormap.
=> https://matplotlib.org/stable/gallery/color/colormap_reference.html
"""

Cmap ="viridis"
Cmap = "hsv"
# Cmap ="magma"
Cmap = "jet_r"
# Cmap = "seismic"
# Cmap = "Spectral"
cmp = matplotlib.colormaps[Cmap]

if FilesOnly:
    matplotlib.use("cairo")


if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)


if MergeModels:
    Models = util.merge_model_sets(infile_list=mod_files,
                                   outfile_name=ModelMergeFile,
                                   dictout= True, out=False)
    mod_files = [ModelMergeFile]

for filein in mod_files:
    start = process_time()
    print("\nNModels read from: %s" % filein)

    Models = numpy.load(filein, allow_pickle=True)

    E = Models["x"][::step]*XYFact
    E_min = numpy.amin(E)
    E_max = numpy.amax(E)
    N = Models["y"][::step]*XYFact
    N_min = numpy.amin(N)
    N_max = numpy.amax(N)
    Z = Models["d"][::step,:]

    DX = Models["mod"][::step]


    if ParaTrans==1:
       DX = numpy.log10(DX)


    if ("image" in ImageType.lower()) or ("contour"in ImageType.lower()):
        xi= numpy.linspace(E_min,E_max,numIndexes[0])
        yi= numpy.linspace(N_min,N_max,numIndexes[1])
        dx = numpy.around(numpy.diff(xi)[0]/XYFact, decimals=0)
        dy = numpy.around(numpy.diff(yi)[0]/XYFact, decimals=0)
        print("Interpolation mesh, dx = "+ str(dx)+" m, dy ="+ str(dy)+" m")
        XI, YI = numpy.meshgrid(xi, yi, indexing="ij" )
        Pnts = numpy.stack([ E.ravel(),  N.ravel()], -1)
        Mesh = numpy.stack([XI.ravel(), YI.ravel()], -1)

        if MaskDist:
            D_tree=scipy.spatial.KDTree(Pnts, leafsize=10,
                                        compact_nodes=True,
                                        copy_data=True,
                                        balanced_tree=True,
                                        boxsize=None)
            mindist, _ = D_tree.query(Mesh, k=1)
            blankdist = mindist>=DistMask*XYFact


        if MaskPoly:
            XIF = XI.flatten().reshape(-1,1)/XYFact
            YIF = YI.flatten().reshape(-1,1)/XYFact
            blankpoly=[]
            for ipnt in numpy.arange(numpy.size(XIF)):
                outside = not util.point_inside_polygon(XIF[ipnt], YIF[ipnt],
                                                        Polygon)
                blankpoly.append(outside)


    pdf_list = []
    for nc in numpy.arange(len(LayList)):

        layl = LayList[nc][0]
        layr = LayList[nc][1]
        layp = LayList[nc][2]
        lunt = LayList[nc][3]
        plim = LayList[nc][4]
        cstp = LayList[nc][5]


        dstr = str(numpy.round(Z[0,layr-1], decimals=0))
        if "scatter"in ImageType.lower():
            titl = layl+" ("+dstr+" m): trn="+str(ParaTrans)
        else:
            titl = layl+" ("+dstr+" m): "+str(ParaTrans)+"/"+InterpMethod[0]+"/"+InterpMethod[1]

        print("Plotting  "+titl)

        D  = DX[:, layr-1].copy()
        D_min = numpy.amin(D)
        D_max = numpy.amax(D)
        print("Models, read   min="+str( D_min)+"   max="+str( D_max))

        Unit = lunt




        fig, ax = matplotlib.pyplot.subplots()
        fig.set_figwidth(FigWidth)


        if ("scatter" in ImageType.lower()):
            print("Scatter Plot")


            D[D<=plim[0]]=plim[0]
            D[D>=plim[1]]=plim[1]

            if Pixelsize:
                Markersize =(72./fig.dpi)
                im = matplotlib.pyplot.scatter(E, N, color='black', marker='.', lw=0, s=(sfac*72./fig.dpi)**2)
            else:
                im = matplotlib.pyplot.scatter(E, N, c=D, s=Markersize**2, cmap=cmp)

            # ax = matplotlib.pyplot.gca()
            ax.set_aspect("equal")
            ax.xaxis.set_major_formatter(xformatter)
            ax.set_xlabel("Easting "+XYUnits, size=Fontsizes[1])
            ax.yaxis.set_major_formatter(yformatter)
            ax.set_ylabel("Northing "+XYUnits, size=Fontsizes[1])

            ax.tick_params(axis="x", labelsize=Fontsizes[1]-2, labelrotation=0.)#-45)
            ax.tick_params(axis="y", labelsize=Fontsizes[1]-2, labelrotation=0.)#-45)
            ax.grid(which="major", axis="both", visible=True,linewidth= Linewidths[0],linestyle="--")
            ax.set_title(AEM_system.upper()+": "+ titl)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            cb = matplotlib.pyplot.colorbar(im, cax=cax, extend="both")
            cb.ax.tick_params(labelsize=Fontsizes[1])
            cb.ax.set_title(lunt, fontsize=Fontsizes[1])

        if ("image" in ImageType.lower()) or ("contour"in ImageType.lower()):

            Dats = D.flatten()

            if "grid" in InterpMethod[0].lower():
                DI = scipy.interpolate.griddata(Pnts, Dats, Mesh,
                                                method=InterpMethod[1].lower())
                DI = numpy.reshape(DI,(len(xi), len(yi)))

            elif "rbf" in InterpMethod[0].lower():
                # RBF = scipy.interpolate.Rbf(E, N, D,
                #                             function=InterpMethod[1].lower(), smooth=InterpMethod[2])
                # DI  = RBF(XI, YI)
                Pnts = numpy.stack([ E.ravel(),  N.ravel()], -1)
                Mesh = numpy.stack([XI.ravel(), YI.ravel()], -1)
                Dats = D.ravel()
                RBF = scipy.interpolate.RBFInterpolator(
                            Pnts, Dats,
                            kernel=InterpMethod[1], smoothing=InterpMethod[2])
                DI = RBF(Mesh)
                DI = numpy.reshape(DI,(len(xi), len(yi)))


            elif "krig" in InterpMethod[0].lower():
                error("Kriging estimation not yet implemented! Exit.")

            D[D<=plim[0]]=plim[0]
            D[D>=plim[1]]=plim[1]

            if MaskPoly:
                DIF = DI.flatten().reshape(-1,1)
                DIF[blankpoly] = numpy.nan
                DI = numpy.reshape(DIF,(len(xi), len(yi)))

            if MaskDist:
                DIF = DI.flatten().reshape(-1,1)
                DIF[blankdist] = numpy.nan
                DI = numpy.reshape(DIF,(len(xi), len(yi)))


            D_min = numpy.nanmin(DI)
            D_max = numpy.nanmax(DI)
            print("Models, interpolated   min="+str( D_min)+"   max="+str( D_max))


            if len(plim)==0:
                if ("image" in ImageType.lower()):
                    im = ax.pcolor(XI, YI, DI, cmap=cmp)
                if ("contour" in ImageType.lower()):
                    im = ax.contourf(XI, YI, DI, cmap=cmp, levels=cstp)
            else:
                if ("image" in ImageType.lower()):
                    valmin, valmax = plim
                    im = ax.pcolor(XI, YI, DI,
                                   cmap=cmp,
                                   vmin=valmin, vmax=valmax)
                if ("contour" in ImageType.lower()):
                    valmin, valmax =  plim
                    im = ax.contourf(XI, YI, DI,
                                     cmap=cmp,
                                     vmin=valmin, vmax=valmax,
                                     levels=cstp)

            ax.set_aspect("equal")
            ax.xaxis.set_major_formatter(xformatter)
            ax.set_xlabel("Easting "+XYUnits, size=Fontsizes[1])
            ax.yaxis.set_major_formatter(yformatter)
            ax.set_ylabel("Northing "+XYUnits, size=Fontsizes[1])

            ax.set_title(titl,fontsize=Fontsize)

            ax.grid(color="k", alpha=0.5, linestyle="dotted", linewidth=1.5)
            ax.tick_params(labelsize=Labelsize)


            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)

            cb = matplotlib.pyplot.colorbar(im, cax=cax, extend="both")
            cb.ax.tick_params(labelsize=Fontsizes[1])
            cb.ax.set_title(lunt, fontsize=Fontsizes[1])



        if "scatter" in ImageType.lower():
            plotfile = PlotDir+PlotName+"_"+AEM_system\
                +"_"+ImageType\
                +"_"+ layl
        else:
            plotfile = PlotDir+PlotName+"_"+AEM_system\
                +"_"+ImageType\
                +"_"+ layl\
            +"_"+InterpMethod[0].lower()\
            +"_"+InterpMethod[1].lower()

        for F in PlotFmt:
            print("Plot written to "+plotfile+F)
            matplotlib.pyplot.savefig(plotfile+F,
                                      dpi=600,
                                      bbox_inches="tight",
                                      backend= "cairo",
                                      transparent=True)


        if PDFCatalog:
            pdf_list.append(plotfile+".pdf")

        matplotlib.pyplot.show()
        matplotlib.pyplot.clf()


if PDFCatalog:
    viz.make_pdf_catalog(PDFList=pdf_list, FileName=PlotDir+PDFCName)