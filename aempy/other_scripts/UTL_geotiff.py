#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 10:39:12 2022

@author: vrath
"""
import os
import sys


import numpy
from osgeo import ogr, gdal
import subprocess

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)

import util

SearchStrng = ".tif"
InDatDir ="/home/vrath/Limerick2022/data_aem/A5_EM_INV_MODELS_OHMM_GRIDS/GeoTiff/"

# uleft_y = 5832233.00
# uleft_x = 497771.00

# lright_y = 5792752.00
# lright_x = 590986.00
upperlat =  52.65
upperlon =  -8.65+360.
lowerlat = 52.45
lowerlon =   -8.15+360.

# uleft_x, uleft_y =util.project_itm_to_utm(uleft_x, uleft_y)
# lright_x, lright_y = util.project_itm_to_utm(uleft_x, uleft_y)

window = (upperlon, lowerlat, lowerlon, upperlat)   #-te xmin ymin xmax ymax
print(window)

data_files = util.get_filelist(searchpath=InDatDir)
data_files = sorted(data_files)


for file in data_files:
    if file.endswith(SearchStrng):
        fileinp = InDatDir+file
        fileout = InDatDir+file.replace("_IDW", "_IDW_croppy")
        gdal.Warp(fileout, fileinp, projWin = window, projWinSRS = "epsg:4326")
# kwargs = {'format': 'GTiff', 'geoloc': True}
# ds = gdal.Warp('C:/test/MYD09.A2011093.0410.006.2015217030905.tif', 'C:/test/tel.vrt', **kwargs)
