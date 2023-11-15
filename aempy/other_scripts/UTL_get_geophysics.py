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
import csv

import numpy
import shapely.geometry as shg

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

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
print("Tellus read geophysics "
      + "".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

InputDir = AEMPYX_DATA+"/Intersection/geophysics/"
OutputDir = InputDir
PolyDir = AEMPYX_DATA+"/Intersection/polygons/"

print("Geophysics items read from: %s" % InputDir)
print("Geophysics items written to: %s" % OutputDir)

if not os.path.isdir(OutputDir):
    print(" File: %s does not exist, but will be created" % OutputDir)
    os.mkdir(OutputDir)

Fmt ="%20s "+"%20.6f "*2+"%15s "*3

# InputFiles = ["Geophysics_2018.csv"]
# columns = [2,3,4,6,5,7]

InputFiles= ["Geophysics_2020.csv"]
columns = [1,3,2,4,5]


print("\nReading from "+InputDir+InputFiles[0])

name, ext = os.path.splitext(InputFiles[0])

Data = []
for file in InputFiles:
    with open(InputDir+file) as fd:
          rows = csv.reader(fd, delimiter=',', )
          h = next(rows)
          if h != None:
              for row in rows:
                 tmp = [row[ii] for ii in columns]
                 # tmp[3] ="ERI"
                 # tmp.insert(4,"NONE")
                 Data.append(tmp)

Data = numpy.asarray(Data, dtype=object)
Data[:,1], Data[:,2] = util.project_itm_to_utm(Data[:,1], Data[:,2], utm_zone=32629)
Ddims=numpy.shape(Data)

print("Writing to "+name+"_utm.npz")
numpy.savez_compressed(OutputDir+name+"_utm.npz", Data=Data, allow_pickle=True)

DataSelect = "Intersection"

AllItems = numpy.empty([1,6])


InPoly= "A1_2019_utm.npz"
PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]
Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)

Items = []
for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Items.append(Data[row, :])
Boreholes = []

for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Boreholes.append(Data[row, :])


if not Items:
    print("No Items found in A1-NM intersection.")
else:

    AllItems =numpy.append( AllItems, Items, axis=0)
    Bdims=numpy.shape(Items)
    print(str(Bdims[0])+" Items found in A1-NM intersection.")
    print("Writing to "+name+"_utm_NM_A1.txt")
    numpy.savetxt(OutputDir+name+"_utm_NM_A1.txt", Items, fmt=Fmt)
    print("Writing to "+name+"_utm_NM_A1.npz")
    numpy.savez_compressed(OutputDir+name+"_utm_NM_A1.npz", Data=Items, allow_pickle=True)

InPoly = "A2_2019_utm.npz"
PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]
Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)
Items = []
for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Items.append(Data[row, :])
if not Items:
    print("No Items found in A2-NM intersection.")
else:
    AllItems =numpy.append( AllItems, Items, axis=0)
    Bdims=numpy.shape(Items)
    print(str(Bdims[0])+" Items found in A2-NM intersection.")
    print("Writing to "+name+"_utm_NM_A2.txt")
    numpy.savetxt(OutputDir+name+"_utm_NM_A2.txt", Items, fmt=Fmt)
    print("Writing to "+name+"_utm_NM_A2.npz")
    numpy.savez_compressed(OutputDir+name+"_utm_NM_A2.npz", Data=Items, allow_pickle=True)


InPoly = "TB_2019_utm.npz"
PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]
Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)
Items = []
for row in numpy.arange(Ddims[0] - 1):
        if util.point_inside_polygon(Data[row, 1], Data[row, 2], Polygon):
            Items.append(Data[row, :])
if not Items:
    print("No Items found in TB-NM intersection.")
else:
    AllItems =numpy.append( AllItems, Items, axis=0)
    Bdims=numpy.shape(Items)
    print(str(Bdims[0])+" Items found in "+name)
    print("Writing to "+name+"_utm_NM_TB.txt")
    numpy.savetxt(OutputDir+name+"_utm_NM_TB.txt", Items, fmt=Fmt)
    print("Writing to "+name+"_utm_NM_TB.npz")
    numpy.savez_compressed(OutputDir+name+"_utm_NM_TB.npz", Data=Items, allow_pickle=True)

AllItems = AllItems[1:,:]
Bdims=numpy.shape(AllItems)
print(str(Bdims[0])+" Items found in "+name)
# print("Writing to "+name+"_utm_overlaps.txt")
# numpy.savetxt(OutputDir+name+"_utm_overlaps.txt", AllItems, fmt=Fmt)
print("Writing to "+name+"_utm_overlaps.npz")
numpy.savez_compressed(OutputDir+name+"_utm_overlaps.npz", Data=AllItems, allow_pickle=True)
