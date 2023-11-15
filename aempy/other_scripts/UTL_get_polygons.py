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
print("Tellus Read Polygons "
      + "".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

InputDir = AEMPYX_DATA+"/RAFTA/polygons/"
OutputDir = InputDir


print("Polygons read from: %s" % InputDir)
print("Polygons writen to: %s" % OutputDir)

if not os.path.isdir(OutputDir):
    print(" File: %s does not exist, but will be created" % OutputDir)
    os.mkdir(OutputDir)


PlyFiles = []

files = os.listdir(InputDir)
for entry in files:
    # print(entry)
    if entry.endswith(".ply"):
        PlyFiles.append(entry)

for file in PlyFiles:
    print("\nReading from "+file)
    name, ext = os.path.splitext(file)
    Poly = numpy.loadtxt(InputDir+file, skiprows=6)

    itm_e = Poly[:, 0]
    itm_n = Poly[:, 1]

    PolyITM = numpy.stack((itm_e, itm_n), axis=1)
    SPolyITM = shg.Polygon(PolyITM)

    utm_e, utm_n = util.project_itm_to_utm(itm_e, itm_n, utm_zone=32629)
    PolyUTM = numpy.stack((utm_e, utm_n), axis=1)
    SPolyUTM = shg.Polygon(PolyUTM)

    lat, lon = util.project_itm_to_latlon(itm_e, itm_n)
    PolyLatLon = numpy.stack((lat, lon), axis=1)
    SPolyLatLon = shg.Polygon(PolyLatLon)

    print("Writing to "+name+"_itm.txt")
    numpy.savetxt(OutputDir+name+"_itm.txt", PolyITM, fmt="%20.6f")
    print("Writing to "+name+"_utm.txt")
    numpy.savetxt(OutputDir+name+"_utm.txt", PolyUTM, fmt="%20.6f")
    print("Writing to "+name+"_latlon.txt")
    numpy.savetxt(OutputDir+name+"_latlon.txt", PolyLatLon, fmt="%14.7f")

    print("Writing Shapely Objects to "+name+"_itm.npz")
    numpy.savez_compressed(OutputDir+name+"_itm.npz", Poly=[SPolyITM], allow_pickle=True)
    print("Writing Shapely Objects to "+name+"_utm.npz")
    numpy.savez_compressed(OutputDir+name+"_utm.npz", Poly=[SPolyUTM], allow_pickle=True)
    print("Writing Shapely Objects to "+name+"_latlon.npz")
    numpy.savez_compressed(OutputDir+name+"_latlon.npz", Poly=[SPolyLatLon], allow_pickle=True)



