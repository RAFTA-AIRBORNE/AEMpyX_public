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
#       jupytext_version: 1.11.4
# ---

"""
Created on Mon Jul 26 10:09:33 2021

@author: vrath
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings

import numpy
import shapely.geometry as shg

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0, pth)

from version import versionstrg
import util

warnings.simplefilter(action="ignore", category=FutureWarning)


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus Read boreholes "
      + "".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

InputDir = AEMPYX_DATA+"/Intersection/boreholes/"
OutputDir = InputDir
PolyDir = AEMPYX_DATA+"/Intersection/polygons/"

print("Boreholes read from: %s" % InputDir)
print("Boreholes writen to: %s" % OutputDir)

if not os.path.isdir(OutputDir):
    print(" File: %s does not exist, but will be created" % OutputDir)
    os.mkdir(OutputDir)


BoreHoles = ["Verified_Boreholes.XYZ"]

print("\nReading from "+InputDir+BoreHoles[0])

name, ext = os.path.splitext(BoreHoles[0])

columns = 6
ncol=numpy.arange(columns)
Data = []
iline = 0
with open(InputDir+BoreHoles[0]) as fd:
    for line in fd:
        iline = iline + 1
        if (line[0].lower().startswith("#")
            or line[0].lower().startswith("/")
            or "line" in line[:24].lower()):
            continue
        t = line.split()
        # print(t)
        tmp = [t[ii] for ii in ncol]
        Data.append(tmp)

Data = numpy.asarray(Data, dtype=object)
Data = numpy.column_stack((Data[:,0],Data[:,1:columns].astype(numpy.float64)))

Data[:,1], Data[:,2] = util.project_itm_to_utm(Data[:,1], Data[:,2], utm_zone=32629)

Ddims=numpy.shape(Data)

# print("Writing to "+name+"_utm.txt")
# numpy.savetxt(OutputDir+name+"_utm.txt", Data, fmt="%20s"+"%20.6f"*5)

print("Writing to "+name+"_utm.npz")
numpy.savez_compressed(OutputDir+name+"_utm.npz", Data=Data, allow_pickle=True)

DataSelect = "Intersection"

InPoly= "A1_2019_utm.npz"
PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]
Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)
Boreholes = []
for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Boreholes.append(Data[row, :])
if not Boreholes:
    print("No Boreholes found in A1-NM intersection.")
else:
    Bdims=numpy.shape(Boreholes)
    print(str(Bdims[0])+" Boreholes found in A1-NM intersection.")
    print("Writing to "+name+"_utm_NM_A1.txt")
    numpy.savetxt(OutputDir+name+"_utm_NM_A1.txt", Boreholes, fmt="%20s"+"%20.6f"*5)
    print("Writing to "+name+"_utm_NM_A1.npz")
    numpy.savez_compressed(OutputDir+name+"_utm_NM_A1.npz", Data=Boreholes, allow_pickle=True)

InPoly = "A2_2019_utm.npz"
PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]
Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)
Boreholes = []
for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Boreholes.append(Data[row, :])
if not Boreholes:
    print("No Boreholes found in A2-NM intersection.")
else:
    Bdims=numpy.shape(Boreholes)
    print(str(Bdims[0])+" Boreholes found in A2-NM intersection.")
    print("Writing to "+name+"_utm_NM_A2.txt")
    numpy.savetxt(OutputDir+name+"_utm_NM_A2.txt", Boreholes, fmt="%20s"+"%20.6f"*5)
    print("Writing to "+name+"_utm_NM_A2.npz")
    numpy.savez_compressed(OutputDir+name+"_utm_NM_A2.npz", Data=Boreholes, allow_pickle=True)


InPoly = "TB_2019_utm.npz"
PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]
Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)
Boreholes = []
for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Boreholes.append(Data[row, :])

if not Boreholes:
    print("No Boreholes found in TB-NM intersection.")
else:
    Bdims=numpy.shape(Boreholes)
    print(str(Bdims[0])+" Boreholes found in "+name)
    print("Writing to "+name+"_utm_NM_TB.txt")
    numpy.savetxt(OutputDir+name+"_utm_NM_TB.txt", Boreholes, fmt="%20s"+"%20.6f"*5)
    print("Writing to "+name+"_utm_NM_TB.npz")
    numpy.savez_compressed(OutputDir+name+"_utm_NM_TB.npz", Data=Boreholes, allow_pickle=True)
