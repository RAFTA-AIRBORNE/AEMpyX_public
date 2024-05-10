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
import scipy.spatial
import skgstat
import shapely

# from evtk import *
import evtk

AEMPYX_pnts_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_pnts_ROOT+"/aempy/modules/", AEMPYX_pnts_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys
import inverse

warnings.simplefilter(action="ignore", category=FutureWarning)
cm = 1/2.54

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

now = datetime.now()

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
    FwdCall,YY, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = YY[0]
    ParaTrans = 1


if "genes" in AEM_system.lower():
    FwdCall, YY, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = YY[0]
    ParaTrans = 1


InModDir = AEMPYX_DATA+"/Blocks/A9/results/"

FileList = "search"
SearchStrng = "A9*k3*.npz"
print("Searchstring: %s \n" % SearchStrng)

FileList = "set"
ListName = ""

MergeModels = False
ModelMergeFile = InModDir+"MUN_k3_models_merged.npz"

if "set" in FileList.lower():
    mod_files = [InModDir+"MUN_k3_models_merged.npz"]

if "read" in FileList.lower():
    print("File names read from : "+ListYame)
    how = ["read", ListYame, InModDir]
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
    error("Yo files set!. X_pntsxit.")



VTKDir = InModDir+"/vtk/"
print("Plots written to dir: %s " % VTKDir)
VTKYame = "MUN_k3_Pnts.vtk"
print("Plot filname: %s " % VTKYame)



# Layers = [5, 10, 15, 20, 25]
# Prop = "rho"
# Unit = r"log10 $\Omega$m"
# Limits = [0., 4.]
# Steps  = [0., 0.5, 1., 1.5, 2., 2.5,  3., 3.5, 4.,]
# Steps = numpy.arange(0.4, 3.6001, 0.2)
# print("\nLayer parameters:")
# LayList = []
# for il in numpy.arange(len(Layers)):
#     label = "Layer "+str(Layers[il])
#     LayList.append([label, Layers[il], Prop, Unit, Limits, Steps])
#     print(LayList[il])

Step = 1
Interp =  False
InterpMethod = ["griddata","linear"]


MaskDist = True
if MaskDist:
    DistMask = 100.

MaskPoly = False
if MaskPoly:
    PolyDir = AEMPYX_DATA+"/Blocks/polygons/"
    PolyFiles = [PolyDir+"A9_2019_utm.npz"]
    Polygon= numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]



# XYFact = 1.
# XYUnits = "(m)"
XYUnits = "(km)"
XYFact = 0.001
Ztype = "depth"  # Ztype = "z"

if not os.path.isdir(VTKDir):
    print("File: %s does not exist, but will be created" % VTKDir)
    os.mkdir(VTKDir)


if MergeModels:
    Models = util.merge_model_sets(infile_list=mod_files,
                                   outfile_name=ModelMergeFile,
                                   dictout= True, out=False)
    mod_files = [ModelMergeFile]

for filein in mod_files:
    start = process_time()
    print("\nModels read from: %s" % filein)

    Models = numpy.load(filein, allow_pickle=True)

    X = Models["x"][::Step]*XYFact
    X_min = numpy.amin(X)
    X_max = numpy.amax(X)
    Y = Models["y"][::Step]*XYFact
    Y_min = numpy.amin(Y)
    Y_max = numpy.amax(Y)
    if "dep" in Ztype:
        Z = Models["d"][::Step,:]
    else:
        Z = Models["z"][::Step,:]



    D = Models["mod"][::Step]
    if ParaTrans==1:
       D = numpy.log10(D)
    S = Models["sns"][::Step]

    D_pnts = D
    S_pnts = S
    X_pnts = numpy.zeros_like(D)
    Y_pnts = numpy.zeros_like(D)
    Z_pnts = numpy.zeros_like(D)
    [sites, layers] = numpy.shape(D)
    ipnts = -1
    for isit in numpy.arange(sites):
        for ilay in numpy.arange(layers):
            X_pnts[ipnts] = X_pnts[isit]
            Y_pnts[ipnts] = Y_pnts[isit]
            Z_pnts[ipnts] = Z_pnts[isit,ilay]

    print(numpy.shape(X_pnts))
    print(numpy.shape(Y_pnts))
    print(numpy.shape(Z_pnts))
    print(numpy.shape(D_pnts))
    print(numpy.shape(S_pnts))
    X_pnts = X_pnts.flatten()
    Y_pnts = Y_pnts.flatten()
    Z_pnts = Z_pnts.flatten()
    D_pnts = D_pnts.flatten()
    S_pnts = S_pnts.flatten()

    # if Interp:
    #     xi= numpy.linspace(X_pnts_min,X_pnts_max,numIndexes)
    #     yi= numpy.linspace(Y_min,Y_max,numIndexes)
    #     dx = numpy.around(numpy.diff(xi)[0]/X_pntsYFact, decimals=0)
    #     dy = numpy.around(numpy.diff(yi)[0]/X_pntsYFact, decimals=0)
    #     print("Interpolation mesh, dx = "+ str(dx)+" m, dy ="+ str(dy)+" m")

    #     X_pntsI, YI = numpy.meshgrid(xi, yi)
    #     Pnts = numpy.stack([ X_pnts.ravel(),  Y.ravel()], -1)
    #     Mesh = numpy.stack([X_pntsI.ravel(), YI.ravel()], -1)
    #     Dats = D.flatten()

    #     if "grid" in InterpMethod[0].lower():
    #         DI = scipy.interpolate.griddata(Pnts, Dats, Mesh,
    #                                         method=InterpMethod[1].lower())
    #         DI = numpy.reshape(DI,(len(xi), len(yi)))

    #     elif "rbf" in InterpMethod[0].lower():
    #         # RBF = scipy.interpolate.Rbf(X_pnts, Y, D,
    #         #                             function=InterpMethod[1].lower(), smooth=InterpMethod[2])
    #         # DI  = RBF(X_pntsI, YI)
    #         Pnts = numpy.stack([ X_pnts.ravel(),  Y.ravel()], -1)
    #         Mesh = numpy.stack([X_pntsI.ravel(), YI.ravel()], -1)
    #         Dats = D.ravel()
    #         RBF = scipy.interpolate.RBFInterpolator(
    #                     Pnts, Dats,
    #                     kernel=InterpMethod[1], smoothing=InterpMethod[2])
    #         DI = RBF(Mesh)
    #         DI = numpy.reshape(DI,(len(xi), len(yi)))


    #     elif "krig" in InterpMethod[0].lower():
    #         error("Kriging estimation not yet implemented! X_pntsxit.")

    #     if MaskDist:
    #         D_tree=scipy.spatial.KDTree(Pnts, leafsize=10,
    #                                     compact_nodes=True,
    #                                     copy_data=True,
    #                                     balanced_tree=True,
    #                                     boxsize=Yone)
    #         mindist, _ = D_tree.query(Mesh, k=1)
    #         blankdist = mindist>=DistMask*X_pntsYFact


    #     if MaskPoly:
    #         X_pntsIF = X_pntsI.flatten().reshape(-1,1)/X_pntsYFact
    #         YIF = YI.flatten().reshape(-1,1)/X_pntsYFact
    #         blankpoly=[]
    #         for ipnt in numpy.arange(numpy.size(X_pntsIF)):
    #             outside = not util.point_inside_polygon(X_pntsIF[ipnt], YIF[ipnt],
    #                                                     Polygon)
    #             blankpoly.append(outside)

        # X_pnts = X_pntsI
        # Y = YI

    if MaskPoly:
        Xmask = X_pnts.flatten().reshape(-1,1)/XYFact
        Ymask = Y_pnts.flatten().reshape(-1,1)/XYFact
        blankpoly=[]
        for ipnt in numpy.arange(numpy.size(Xmask)):
            outside = not util.point_inside_polygon(Xmask[ipnt], Ymask[ipnt],
                                                    Polygon)
            blankpoly.append(outside)

        nD = numpy.shape(D_pnts)
        # D_pnts = D_pnts.flatten().reshape(-1,1)
        D_pnts[blankpoly] = numpy.nan
        D_pnts = numpy.reshape(D_pnts,nD)
        # S_pnts =  S_pnts.flatten().reshape(-1,1)
        S_pnts[blankpoly] = numpy.nan
        S_pnts = numpy.reshape(S_pnts,nD)

        X_pnts[blankpoly] = numpy.nan
        Y_pnts[blankpoly] = numpy.nan

    """
    Now store to VTK

    """
    Z_pnts = -Z_pnts

    D_pnts[numpy.logical_or(D_pnts<-1., D_pnts<4.)] = numpy.nan

    print("To VTK:")
    print(numpy.shape(X_pnts))
    print(numpy.shape(Y_pnts))
    print(numpy.shape(Z_pnts))
    print(numpy.shape(D_pnts))
    print(numpy.shape(S_pnts))

    evtk.hl.pointsToVTK(VTKDir+VTKYame, X_pnts, Y_pnts, Z_pnts,
                        data = {"res" : D_pnts,"sens" : D_pnts})

    # evtk.hl.pointsToVTK(VTKDir+VTKYame, X_pnts, Y, Z, data = {"temp" : temp, "pressure" : pressure})
