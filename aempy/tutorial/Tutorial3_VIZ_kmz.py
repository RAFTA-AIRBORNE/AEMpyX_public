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
#       jupytext_version: 1.16.2
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
import skgstat
import shapely

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys
import viz
import inverse

warnings.simplefilter(action="ignore", category=FutureWarning)
cm = 1/2.54

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

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
# AEM_system = "genesis"
AEM_system = "aem05"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.03
    data_active = numpy.ones(NN[2], dtype="int8")
    CompDict=Misc[2]
    CompLabl = list(CompDict.keys())
    print(CompLabl)

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 2
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + "good" hoizontals"
    CompDict =Misc[2]
    CompLabl = list(CompDict.keys())




"""
input formats are "npz","nc4","ascii"
"""
InFilFmt = ".npz"
InDatDir = AEMPYX_DATA+"/Blocks/A5/raw/"
# # InDatDir = AEMPYX_DATA+"/Blocks/A1O/proc_delete/data_asc/"
print("Data read from dir: %s " % InDatDir)
FileList = "set" #"search"
SearchStrng = ""


PlotDir = InDatDir+"/plots/"
print("Plots written to dir: %s " % PlotDir)
PlotName = "A5"
print("Plot filname: %s " % PlotName)

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = [InDatDir+"A5_Merged.npz"]

else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, full=True, sort=True)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")



MergeData = True
DataMergeFile = InDatDir+"A5_Data_Merged.npz"
MergeModels = False
ModelMergeFile = InDatDir+"A5_Models_Merged.npz"


PlotModel = False
NLyr = 25
Layer = 10
if PlotModel:
    nlyr =NLyr
    numbers = numpy.arange(1,nlyr+1)
    keys = ["R"+str(ii) for ii in numbers]
    Compdict = dict(zip(keys, numbers))

PlotRMS = False
if PlotRMS:
    pass


ImageType = "image"
# ImageType = "contour"
# ImageType = "scatter"

Pixels = 2048 * 10
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
        "nearest"       data point closest to the point of interpolation
        "linear"        tessellate the input point set to N-D simplices
                        and interpolate linearly on each simplex
        "cubic"         return the value determined from a piecewise cubic,
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
    numIndexes = 1001
    smooth = 0.
    Levels = []
    MaskNeg = False
    MaskPoly = True
    if MaskPoly:
        PolyDir = AEMPYX_DATA+"/Blocks/polygons/"
        PolyFiles = [PolyDir+"A5_2019_utm.npz"]
        Polygon= numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]

if ("scatter" in ImageType.lower()):
    Decim=25
    step = min(1,abs(Decim))
    Markersize = 5
    MaskNeg = True
    Pixelsize = False
    if Pixelsize:
        step=10
        sfac = 10.

# CompList=[
#     ["Z3", []], #[0., 2000., 100.]],
#     ["Z6", []],#[0., 2000., 100.]],
#     ["Z9", []], #[0., 2000., 100.]],
#     ["H3", []], #[0., 2000., 100.]],
#     ["H6", []],#[0., 2000., 100.]],
#     ["H9", []], #[0., 2000., 100.]],
#     ["ALT", [80., 160., 20.], 240.]     # ALTthresh = 70.
          # ]

CompList=[
    # ["P1", [0., 2000., 100.]],
    # ["Q1", [0., 2000., 100.]],
    # ["P2", [0., 2000., 100.]],
    # ["Q2", [0., 2000., 100.]],
    ["P3", [0., 2000., 100.]],
    ["Q3", [0., 2000., 100.]],
    # ["P4", [0., 2000., 100.]],
    # ["Q4", [0., 2000., 100.]],
    #["PLM", [], 0.2],      # PLMthresh = 0.25
    # ["ALT", [40., 120., 20.], 300.]     # ALTthresh = 70.
          ]




# XYUnits = "(m)"
# xformatter = matplotlib.ticker.FormatStrFormatter("%7f")
# yformatter = matplotlib.ticker.FormatStrFormatter("%6f")
# XYFact = 1.
XYUnits = "(km)"
XYFact = 0.001
xformatter = matplotlib.ticker.FormatStrFormatter("%.2f")
yformatter = matplotlib.ticker.FormatStrFormatter("%.2f")
"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
"""
FilesOnly = False
matplotlib.pyplot.style.use("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
matplotlib.rcParams["savefig.bbox"]= "tight"

Fontsize = 8
Labelsize = Fontsize
Titlesize = 8
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidths= [0.5]
FigWidth = 16.
Pixels = 2048

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


if MergeData:
    Data = util.merge_data_sets(infile_list=dat_files,
                                outfile_name=DataMergeFile,
                                aem_system="aem05", out= True)
    dat_files = [DataMergeFile]

if MergeModels:
    Data = util.merge_model_sets(infile_list=dat_files,
                                outfile_name=ModelMergeFile,
                                aem_system="aem05", out= True)
    dat_files = [ModelMergeFile]

for filein in dat_files:
    start = process_time()
    print("\nData read from: %s" % filein)
    Data, header, _ = aesys.read_aempy(File=filein, System=AEM_system, OutInfo=False)


    E = Data[:,1][::step]*XYFact
    E_min = numpy.amin(E)
    E_max = numpy.amax(E)
    N = Data[:,2][::step]*XYFact
    N_min = numpy.amin(N)
    N_max = numpy.amax(N)
    Z = Data[:,5][::step]

    Lat, Lon = util.project_utm_to_latlon(E, N)
    LatMax = numpy.amax(Lat)
    LatMin = numpy.amin(Lat)
    LonMax = numpy.amax(Lon)
    LonMin = numpy.amin(Lon)


    # for Comp in CompList:
    for nc in numpy.arange(len(CompList)):

        Comp = CompList[nc][0]

        comp = CompDict[Comp][0]
        indx = CompLabl.index(Comp)

        if "scatter"in ImageType.lower():
            titl = CompLabl[indx]+CompDict[Comp][2]+"  trans="+str(DataTrans)
        else:
            titl = CompLabl[indx]+CompDict[Comp][2]+"  trans="+str(DataTrans)+" / "+InterpMethod[0]+" / "+InterpMethod[1]

        print("Plotting component "+titl)
        D = Data[:,comp][::step]

        if ("Z" in Comp) or ("H" in Comp):
            Unit = "ppm"
            if DataTrans ==1:
                Unit = "-"
                D = numpy.log10(D)
            if DataTrans ==2:
                Unit = "-"
                if not numpy.isfinite(S):
                   S = inverse.get_S(D)
                D= numpy.arcsinh(D/S)
                print("Scaling Value S for arcsinh:"+str(S))

        if ("P" in Comp) or ("Q" in Comp):
            Unit = "ppm"
            if DataTrans ==1:
                Unit = "-"
                D = numpy.log10(D)
            if DataTrans ==2:
                Unit = "-"
                if not numpy.isfinite(S):
                    S = inverse.get_S(D)
                D= numpy.arcsinh(D/S)
                print("Scaling Value S for arcsinh:"+str(S))

        if ("PL" in Comp):
            Unit = "(-)"
            PLMthresh= CompList[nc][2]
            titl = titl+" / thresh = "+str(PLMthresh)+" m"

        if ("A" in Comp):
            Unit = "m"
            ALTthresh= CompList[nc][2]
            titl = titl+" / thresh = "+str(ALTthresh)+" m"

        if ("R" in Comp):
            Unit = "Ohm.m"



        fig1, ax = viz.gearth_fig(llcrnrlon=LonMin,
                     llcrnrlat=LatMin,
                     urcrnrlon=LonMax,
                     urcrnrlat=LatMax,
                     pixels=Pixels)

        if ("scatter" in ImageType.lower()):
            print("Scatter Plot")

            plotfile_base = PlotDir+PlotName+"_"+AEM_system\
                +"_"+ImageType\
                +"_"+ Comp

            if Pixelsize:
                Markersize =(72./fig1.dpi)

            if ("PL" in Comp):
                D[D<=PLMthresh]=numpy.nan
            if ("A" in Comp):
                D[D>=ALTthresh]=numpy.nan

            if MaskNeg:
                D[D<=0.]= numpy.nan

            if Pixelsize:
                im = matplotlib.pyplot.scatter(E, N, color="black", marker=".", lw=0, s=(sfac*72./fig1.dpi)**2)
            else:
                im = matplotlib.pyplot.scatter(E, N, c=D, s=Markersize**2, cmap=cmp)

            ax.set_axis_off()

        if ("image" in ImageType.lower()) or ("contour"in ImageType.lower()):
            print("Contour/Image Plot")



            xi= numpy.linspace(E_min,E_max,numIndexes)
            yi= numpy.linspace(N_min,N_max,numIndexes)
            dx = numpy.around(numpy.diff(xi)[0]/XYFact, decimals=0)
            dy = numpy.around(numpy.diff(yi)[0]/XYFact, decimals=0)
            print("Interpolation mesh, dx = "+ str(dx)+" m, dy ="+ str(dy)+" m")

            XI, YI = numpy.meshgrid(xi, yi)
            Mesh = numpy.vstack((XI.flatten(), YI.flatten())).T
            Pnts = numpy.vstack((E.flatten(), N.flatten())).T
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
                RBF = scipy.interpolate.RBFInterpolator(Pnts, Dats,
                                                     kernel=InterpMethod[1], smoothing=InterpMethod[2])
                DI = RBF(Mesh)
                DI = numpy.reshape(DI,(len(xi), len(yi)))


            elif "krig" in InterpMethod[0].lower():
                error("Kriging estimation not yet implemented! Exit.")


            if ("PL" in Comp):
                DI[DI<=PLMthresh]=numpy.nan
            if ("A" in Comp):
                DI[DI>=ALTthresh]=numpy.nan

            if MaskNeg:
                DI[DI<=0.]= numpy.nan

            if MaskPoly:
                DIF = DI.flatten().reshape(-1,1)
                XIF = XI.flatten().reshape(-1,1)/XYFact
                YIF = YI.flatten().reshape(-1,1)/XYFact
                for ipnt in numpy.arange(numpy.size(DIF)):
                    outside = not util.point_inside_polygon(XIF[ipnt], YIF[ipnt], Polygon)

                    if outside:
                         DIF[ipnt] = numpy.nan
                DI = numpy.reshape(DIF,(len(xi), len(yi)))


            if len(CompList[nc][1])==0:
                if ("image" in ImageType.lower()):
                    im = ax.pcolor(XI, YI, DI, cmap=cmp)
                if ("contour" in ImageType.lower()):
                    im = ax.contourf(XI, YI, DI, cmap=cmp, levels=Levels)
            else:
                if ("image" in ImageType.lower()):
                    valmin, valmax, _ = CompList[nc][1]
                    im = ax.pcolor(XI, YI, DI,
                                   cmap=cmp,
                                   vmin=valmin, vmax=valmax)
                if ("contour" in ImageType.lower()):
                    valmin, valmax, valdel = CompList[nc][1]
                    levels = numpy.arange(valmin, valmax+valdel, valdel)
                    im = ax.contourf(XI, YI, DI,
                                     cmap=cmp,
                                     vmin=valmin, vmax=valmax,
                                     levels=levels)
            ax.set_axis_off()

            plotfile_base = PlotDir+PlotName+"_"+AEM_system\
                +"_"+ImageType\
                +"_"+ Comp\
                +"_"+InterpMethod[0].lower()\
                +"_"+InterpMethod[1].lower()


        matplotlib.pyplot.show()
        matplotlib.pyplot.clf()


        f1 =plotfile_base+"_image.png"
        fig1.savefig(f1, transparent=False, format="png")


        fig2 = matplotlib.pyplot.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
        ax = fig2.add_axes([0.0, 0.05, 0.2, 0.9])
        cb = fig2.colorbar(cmp, cax=ax)
        cb.set_label("Mean Dynamic Topography [m]", rotation=-90, color="k", labelpad=20)
        f2 = plotfile_base+"_legend.png"
        fig2.savefig(f2, transparent=False, format="png")  # Change transparent to True if your colorbar is not on space :)

        viz.make_kml(llcrnrlon=LonMin,
                     llcrnrlat=LatMin,
                     urcrnrlon=LonMax,
                     urcrnrlat=LatMax,
                     figs=[f1], colorbar=f2,
                     kmzfile=plotfile_base+"kmz", name="Mean Dynamic Topography and velocity")
