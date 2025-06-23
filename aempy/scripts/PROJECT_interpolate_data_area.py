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
#   kernelspec:
#     display_name: Python 3 (Spyder)
#     language: python3
#     name: python3
# ---

# +
#!/usr/bin/env python3
# -

# This script plots data over an spatial area.

# +
import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings
import getpass
import inspect

import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import matplotlib.axis
import mpl_toolkits.axes_grid1

import scipy.interpolate
import scipy.spatial
import shapely

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

# +
OutInfo = True
cm = 1/2.54
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
# script = "Tutorial1_VIZ_data_area.py"
script = inspect.getfile(inspect.currentframe())  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")
Header = titstrng
# -

# The following cell gives values to AEM-system related settings.
#
# Data transformation is activated by the variable _DataTrans_. Currently
# three possible options are allowed: _DataTrans = 0_: No transformation,
# i.e., the raw data are used. _DataTrans = 1_: The natural log of data
# is taken, only allowed for strictly positive values. _DataTrans = 2_:
# If data scale logarithmically, an _asinh_ transformation (introduced by
# Scholl, 2000) is applied. It allows negatives, which may occur in TDEM,
# when IP effects are present.
#
# A general additive/multiplicative error model is applied on the raw data
# before transformation, and errors are also transformed.

# +
# AEM_system = "genesis"
AEM_system = "aem05"
if "aem05" in AEM_system.lower():
    _, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.05
    data_active = numpy.ones(NN[2], dtype="int8")

    CompDict = Misc[3]
    CompLabl = list(CompDict.keys())
    print(CompLabl)

if "genes" in AEM_system.lower():
    _, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
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

# +

InFileFmt = ".npz"


##############################################################################
A2 Somaye
##############################################################################
# raw data
InDatDir =  AEMPYX_DATA + "/AllDataBlocks/Archive/A2/raw/"
SearchStrng = "*Full_data*.npz"
FileList = "search"
ProfileName = ""


print("Data read from dir: %s " % InDatDir)
print("Plots written to dir: %s " % PlotDir)
print("Plot filname: %s " % PlotName)

# FileList = "read"
ListName = ""
# -

if "set" in FileList.lower():
    dat_files = []

if "read" in FileList.lower():
    print("File names read from : "+ListName)
    how = ["read", ListName, InDatDir]
    dat_files = util.get_data_list(how=how,
                              out= True, sort=True)

    dat_files = numpy.loadtxt("A9-7.dat", dtype=str)

if "search" in FileList.lower():
    print("Searchstring is : "+SearchStrng)
    how = ["search", SearchStrng, InDatDir]
    dat_files = util.get_data_list(how=how, fullpath=True,
                              out= True, sort=True)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")

MergeData = False
DataMergeFile = InDatDir+ProfileName+"_merged.npz"


#
# Kernel functions for RBF:
#     The radial basis function, based on the radius, r,
#     given by the norm (default is Euclidean distance); the default is ‘multiquadric’:
#         ‘linear’ : -r
#         ‘thin_plate_spline’ : r**2 * log(r)
#         ‘cubic’ : r**3
#         ‘quintic’ : -r**5
#
# If a callable, then it must take 2 arguments (self, r). The epsilon parameter
# will be available as self.epsilon. Other keyword arguments passed
# in will be available as well.
#
#
# Methods for griddata:
#         'nearest'       data point closest to the point of interpolation
#         'linear'        tessellate the input point set to N-D simplices
#                         and interpolate linearly on each simplex
#         'cubic'         return the value determined from a piecewise cubic,
#                         continuously differentiable (C1), and approximately
#                         curvature-minimizing polynomial surface.
#

if ("image" in ImageType.lower()) or ("contour"in ImageType.lower()):
    step = 1

    InterpMethod = ["griddata","linear"]
    # InterpMethod = ["griddata", "cubic"]
    # InterpMethod = ["rbf", "linear", 0.0]
    # InterpMethod = ["rbf", "thin_plate_spline", 0.0]
    # InterpMethod = ["rbf", "cubic", 0.01]

    # InterpMethod = ["krig", "linear", 0.5, 340.]


    numIndexes = [121, 141]
    smooth = 0.
    Levels = []
    MaskNeg = True
    MaskPoly = False
    MaskDist = True

    if MaskDist:
        DistMask = 100.*XYFact

    if MaskPoly:
        PolyDir = AEMPYX_DATA+"/Blocks/polygons/"
        PolyFiles = [PolyDir+"A5_2019_utm.npz"]
        Polygon= numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]


# The following cell determines the settings for individual components. Each sublist associated to a componet contains the name, followed by a list of parameters determining the data limits, and a step determining the color bar, or the isolines. Further paramers, as e.g. the threshhold for the PLM, may be added.


CompList=[
    ["P1", [0., 2000., 200.]],
    ["Q1", [0., 2000., 200.]],
    ["P2", [0., 2000., 200.]],
    ["Q2", [0., 2000., 200.]],
    ["P3", [0., 2000., 200.]],
    ["Q3", [0., 2000., 200.]],
    ["P4", [0., 2000., 200.]],
    ["Q4", [0., 2000., 200.]],
    ["PLM", [], 0.5],      # PLMthresh = 0.25
    ["ALT", [40., 120., 20.], 300.]   # ALTthresh = 120
]


if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)


if MergeData:
    Data = inverse.merge_data_sets(infile_list=dat_files,
                                outfile_name=DataMergeFile,
                                aem_system="aem05", out= True)
    dat_files = [DataMergeFile]

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

    xi= numpy.linspace(E_min,E_max,numIndexes[0])
    yi= numpy.linspace(N_min,N_max,numIndexes[1])
    dx = numpy.around(numpy.diff(xi)[0]/XYFact, decimals=0)
    dy = numpy.around(numpy.diff(yi)[0]/XYFact, decimals=0)
    print("Interpolation mesh, dx = "+ str(dx)+" m, dy ="+ str(dy)+" m")
    XI, YI = numpy.meshgrid(xi, yi, indexing="ij")
    Pnts = numpy.stack([ E.ravel(),  N.ravel()], -1)
    Mesh = numpy.stack([XI.ravel(), YI.ravel()], -1)

    if MaskDist:
        D_tree=scipy.spatial.KDTree(Pnts, leafsize=10,
                                    compact_nodes=True,
                                    copy_data=True,
                                    balanced_tree=True,
                                    boxsize=None)
        mindist, _ = D_tree.query(Mesh, k=1)
        blankdist = mindist>=DistMask


    if MaskPoly:
        XIF = XI.flatten().reshape(-1,1)/XYFact
        YIF = YI.flatten().reshape(-1,1)/XYFact
        blankpoly=[]
        for ipnt in numpy.arange(numpy.size(XIF)):
            outside = not util.point_inside_polygon(XIF[ipnt], YIF[ipnt],
                                                    Polygon)
            blankpoly.append(outside)

    for nc in numpy.arange(len(CompList)):

        Comp = CompList[nc][0]

        comp = CompDict[Comp][0]
        indx = CompLabl.index(Comp)

        else:
            titl = CompLabl[indx]+CompDict[Comp][2]+": "\
                +PlotStrng+"/"+str(DataTrans)\
                +"/"+InterpMethod[0]+"/"+InterpMethod[1]

        print("Plotting component "+titl)
        D = Data[:,comp][::step]

        D_min = numpy.amin(D)
        D_max = numpy.amax(D)
        print("Data, read   min="+str( D_min)+"   max="+str( D_max))

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
            titl = titl # +"\nthresh = "+str(PLMthresh)

        if ("A" in Comp):
            Unit = "m"
            ALTthresh= CompList[nc][2]
            titl = titl # +"\nthresh = "+str(ALTthresh)+" m"

        if ("R" in Comp):
            Unit = "Ohm.m"






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


        if ("PL" in Comp):
            DI[DI<=PLMthresh]=numpy.nan
        if ("A" in Comp):
            DI[DI>=ALTthresh]=numpy.nan

        if MaskNeg:
            DI[DI<=0.]= numpy.nan

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
        print("Data, interpolated   min="+str( D_min)+"   max="+str( D_max))


