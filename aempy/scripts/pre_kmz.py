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

# ---

"""
@author: vr May 2022
"""

# Import required modules

import os
import sys
from sys import exit as error
import csv
import warnings
from time import process_time
from datetime import datetime
import simplekml
import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)


import aesys
import util
import viz

from version import versionstrg

warnings.simplefilter(action="ignore", category=FutureWarning)

AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = False
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
# AEM_systemtem = "genesis"
AEM_system = "aem05"

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]


# Define the path to your data-files
# DataDir =  AEMPYX_ROOT + "/work/data/raw/nan/"
# print(" data files read from: %s" % DataDir)
# PlotDir  =  AEMPYX_ROOT + "/work/data/raw/plots/"
# print(" plots read from: %s" % PlotDir)

DataDir =  AEMPYX_ROOT + "/work/data//proc_delete_PLM3s//nan/"
print(" data files read from: %s" % DataDir)
PlotDir  =  AEMPYX_ROOT + "/work/data/proc_delete_PLM3s/plots/"
print(" plots read from: %s" % PlotDir)


SearchStrng = "*.npz"
data_files = util.get_filelist(searchstr=[SearchStrng], searchpath=DataDir, fullpath=False)
data_files = sorted(data_files)
ns = numpy.size(data_files)

KMLDir = DataDir
KLMFile = KMLDir+"Limerick_shale_proc"

MarkStartPoints = True
MarkEndPoints = False
MarkCenterPoints = False
AddImages = True
ImageWidth= 600
plots_fmt = ".png"

MarkEvery = 50


# Determine what is added to the KML-tags:


kml = False
kmz = True

# Define the path for saving  kml files

icon_dir = AEMPYX_ROOT+"/aempy/share/icons/"

line_icon =  icon_dir + "star.png"
line_iscale = 1.5
line_icolor = simplekml.Color.yellow
line_tscale = 1.  # sc
line_tcolor = simplekml.Color.yellow

data_icon =  icon_dir + "square.png"
data_iscale = 0.8
data_icolor = simplekml.Color.red
data_tscale = 1.
data_tcolor = simplekml.Color.yellow
# simplekml.Color.rgb(0, 0, 255)
# "ffff0000"


# Determine which geographical info is added to the KML-tags:
# define empty list
kml = simplekml.Kml(open=1)
line_iref = kml.addfile(line_icon)
data_iref = kml.addfile(data_icon)


if (not os.path.isdir(DataDir)) or (not os.path.isdir(PlotDir)):
    error(" File: %s or %s does not exist! Exit." % (DataDir, PlotDir))

for f in data_files:
    print(f)

    file = DataDir+f
    name,  ext = os.path.splitext(f)
    Data, _, _ = aesys.read_aempy(File=file, System=AEM_system, OutInfo=False)

    data = Data
    nd = numpy.shape(data)[0]
    lat, lon = util.project_utm_to_latlon(data[:,1], data[:,2])
    line = str(round(data[0,0],2)).replace(".","-")

    folder_line = kml.newfolder(name="FL"+line)

    for idt in numpy.arange(nd):

        if numpy.mod(idt, MarkEvery) == 0:
            d = folder_line.newpoint()
            d.coords = [(lon[idt], lat[idt])]
            d.style.iconstyle.icon.href = data_iref
            d.style.iconstyle.scale = data_iscale
            d.style.iconstyle.color = data_icolor
            d.description = AEM_system.upper()+"\nFlightline: "+str(data[idt,0])


    if AddImages:
        d_plot = PlotDir+name+ plots_fmt
        if os.path.exists(d_plot)==True:
            src= kml.addfile(d_plot)
            imstring ='<img width="'+str(ImageWidth)+'" align="left" src="' + src + '"/>'
            # imstring = '<img width="1200" align="left" src="' + src + '"/>'
            d.description = (imstring)
        else:
            print(d_plot+ " does not exist!")

    if MarkStartPoints:
        d = folder_line.newpoint(name="S:"+str(data[0,0]))
        d.style.labelstyle.color = data_tcolor
        d.style.labelstyle.scale = data_tscale
        d.style.iconstyle.icon.href = data_iref
        d.style.iconstyle.scale = data_iscale*1.5
        d.style.iconstyle.color = line_icolor
        d.coords = [(lon[0], lat[0])]
        d.description = (imstring)
    if MarkEndPoints:
        d = folder_line.newpoint(name="E:"+str(data[0,0]))
        d.style.labelstyle.color = data_tcolor
        d.style.labelstyle.scale = data_tscale
        d.style.iconstyle.icon.href = data_iref
        d.style.iconstyle.scale = data_iscale*1.5
        d.coords = [(lon[nd-1], lat[nd-1])]
        d.description = (imstring)
    if MarkCenterPoints:
        d = folder_line.newpoint(name=str(data[round(nd/2),0]))
        d.coords = [(lon[round(nd/2)], lat[round(nd/2)])]
        d.style.labelstyle.color = data_tcolor
        d.style.labelstyle.scale = data_tscale
        d.style.iconstyle.icon.href = data_iref
        d.style.iconstyle.scale = data_iscale*1.5
        d.style.iconstyle.color = line_icolor
        d.description = (imstring)


    # Compressed kmz file:

kml.savekmz(KLMFile + ".kmz")
