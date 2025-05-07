#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 20:18:33 2022

@author: vrath
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings
import csv
import simplekml

import numpy
import shapely.geometry as shg

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0, pth)

from version import versionstrg
import util
import aesys

warnings.simplefilter(action="ignore", category=FutureWarning)


AEMPYX_DATA = os.environ["AEMPYX_DATA"]


version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=inspect.getfile(inspect.currentframe()), out=False)
print(titstrng+"\n\n")

# AEM_system = "genesis"
AEM_system = "aem05"

FileList = "search"  # "search", "read"
InDatDir =  AEMPYX_DATA + "/Projects/Compare_systems/data_reduced/"
SearchStrng = "SGL*proc_reduced.npz"

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = []

else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(dat_files)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")

KMLDir = AEMPYX_DATA + "/Projects/Compare_systems/"

if not os.path.isdir(KMLDir):
    print(" File: %s does not exist, but will be created" %KMLDir)
    os.mkdir(KMLDir)


sitestep = 5

#  Determine what is added to the KML-tags:
kml = simplekml.Kml(open=1)

# Define the path for saving  kml files
kml_dir =KMLDir
kml_file = kml_dir+"BundoranSGL_Data"
writekml = False
writekmz = True
icon_dir = AEMPYX_ROOT+"/aempy/share/icons/"

aem_icon =  icon_dir + "circle.png"
# aem_icon =  icon_dir + "square.png"
aem_icolor = simplekml.Color.blue
aem_tcolor = simplekml.Color.blue
aem_iscale = 1.0
aem_tscale = 1.0
aem_iref = kml.addfile(aem_icon)

annotation = []

# EPSG_in = 32629
# pname ="FL13490"
# lat, lon= util.project_utm_to_latlon(tdsite[:,1], tdsite[:,2],
#                                 utm_zone=32629)
for file in dat_files:
    
    
    name, ext = os.path.splitext(file)
    
    filein = InDatDir+file
    print("\n Reading file " + filein)
    Data, Header, _ = aesys.read_aempy(File=filein,
                                   System=AEM_system, OutInfo=False)

    fl_name = Data[0, 0]
    site_x = Data[:, 1]
    site_y = Data[:, 2]
    site_gps = Data[:, 3]
    nsite = len(site_x)
    
    lat, lon= util.project_utm_to_latlon(site_x, site_y,
                                    utm_zone=32629)
    folder = kml.newfolder(name=str(fl_name))
    
    for ii in numpy.arange(nsite)[0:-1:sitestep]:
        # pname = "FL"+str(fl_name)
        pname = str(ii)
        if "end" in annotation:
            if ii==len(lat)-1:
                site = folder.newpoint(name=pname)
        elif "start" in annotation:
            if ii==0:
                site = folder.newpoint(name=pname)            
        elif "cent" in annotation:
            if ii==nsite//2:
                site = folder.newpoint(name=pname)            
        else:
                site = folder.newpoint()
    
        site.coords = [(lon[ii], lat[ii], 0.)]
        
        site.style.labelstyle.color = aem_tcolor
        site.style.labelstyle.scale = aem_tscale
        site.style.iconstyle.icon.href = aem_iref
        site.style.iconstyle.scale = aem_iscale
        site.style.iconstyle.color = aem_icolor
        site.description = "" #pname+str(ii)

    
    # Compressed kmz file:
    kml.savekmz(kml_file + ".kmz")
