#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:sphinx,ipynb
#     text_representation:
#       extension: .py
#       format_name: sphinx
#       format_version: '1.1'
#       jupytext_version: 1.16.2
# ---


# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: title,-all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

"""

This script present a work flow preparing data for further
preprocessing with PREP_Process.

Current options are :
    - change projection for spatial coordinates
    - split data into flightOutLines
    - choose data subset by rectangle (polygon in preparation)
    - write as ASCII file and NETCDF
    (- project data to profile in preparation)

@author: vrath nov 2020
"""
import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings
# import pickle

import numpy


AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys

warnings.simplefilter(action="ignore", category=FutureWarning)


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

now = datetime.now()
Header = titstrng


OutInfo = True


OutFileFmt = ".npz"#".asc" #".npz"
FileList = "search"  # "search", "read"
LinesOut = True
LinesMin = 30
CheckNaN = True
MergeOut = True
SearchStr = ".xyz"
"""
Change data projection

"""
SetProj = False
if SetProj:
    ProjInp = ""
    ProjOut = ""

"""
Put all flightlines in same direction
"""
TellusAng = 345.
Spread = 5.
CorrectDirection = True


"""
Data selection

"""
""

RectCorners = []
PolyFiles = []
DataSelect = ""

###############################################################################
# Full Blocks
# ##############################################################################
AEM_system = "aem05"
_, NN, _, _, _, = aesys.get_system_params(AEM_system)
nD = NN[0]

# InDatDir = AEMPYX_DATA+"/Blocks/A1/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A1/raw/"
# InSurvey = "A1"
# InPoly = "A1_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A2/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A2/raw/"
# InSurvey = "A2"
# InPoly = "A2_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A3/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A3/raw/"
# InSurvey = "A3"
# InPoly = "A3_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A4/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A4/raw/"
# InSurvey = "A4"
# InPoly = "A4_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A5/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A5/raw/"
# InSurvey = "A5"
# InPoly = "A5_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A6/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A6/raw/"
# InSurvey = "A6"
# InPoly = "A6_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A7/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A7/raw/"
# InSurvey = "A7"
# InPoly = "A7_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A8/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A8/raw/"
# InSurvey = "A8"
# # InPoly = "A7_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/A9/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/A9/raw/"
# InSurvey = "A9"
# # InPoly = "A7_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/TB/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/TB/raw/"
# InSurvey = "TB"
# InPoly = "TB_2019_utm.npz"
# OutStrng = InSurvey

# InDatDir = AEMPYX_DATA+"/Blocks/WF/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/WF/raw/"
# InSurvey = "WF"
# InPoly = "WF_2019_utm.npz"
# OutStrng = InSurvey

###############################################################################
# StGormans
# ##############################################################################

# DataSelect = "Rectangle"   # "Polygon", "Intersection", "Union"
# InDatDir = AEMPYX_DATA+"/Blocks/A1/orig/"
# OutDatDir = AEMPYX_DATA+"/Projects/StGormans/raw/"
# RectCorners = [638968.67, 5922331.93,  641519.17, 5924940.46]  # StGormans
# InSurvey = "A1"
# OutStrng = InSurvey+"_rect_StGormans"
###############################################################################
# LoughGur
# ##############################################################################
# DataSelect = "Rectangle"   # "Polygon", "Intersection", "Union"
# InDatDir = AEMPYX_DATA+"/Blocks/A5/orig/"
# OutDatDir = AEMPYX_DATA+"/Projects/LoughGur/raw/"
# RectCorners = [529600., 5816800., 534200., 5820250.] # StGormans
# InSurvey = "A5"
# OutStrng = InSurvey+"_rect_LoughGur"

###############################################################################
# # Limerick
# ###############################################################################
# DataSelect = "Rectangle"   # "Polygon", "Intersection", "Union"
# InDatDir = AEMPYX_ROOT+"/work/data/"
# OutDatDir = AEMPYX_DATA+"/Projects/LoughGur/raw/"
# RectCorners = [486000., 5815000., 498000., 5828000.] # StGormans
# InSurvey = "A5"
# OutStrng = InSurvey+"_rect_shale"

###############################################################################
# Munster
# ##############################################################################
# ##############################################################################
# Munster
# ##############################################################################

DataSelect = "Rectangle"   # "Polygon", "Intersection", "Union"
InDatDir = AEMPYX_DATA+"/Blocks/A9/orig/"
OutDatDir = AEMPYX_DATA+"/Projects/Munster/raw/"
# RectCorners = [516860.94, 5786658.92,   536925.58, 5768259.04 ]# Munster
RectCorners = [516000., 5768000.,   541000., 5788000. ]# Munster
InSurvey = "A9"
OutStrng = InSurvey+"_rect_Munster"


###############################################################################
# CGG NM
# ##############################################################################
# AEM_system = "genesis"
# _, NN, _, _, _, = aesys.get_system_params(AEM_system)
# nD = NN[0]
# InDatDir = AEMPYX_DATA+"/Blocks/NM/orig/"
# OutDatDir = AEMPYX_DATA+"/Blocks/NM/raw/"
# InSurvey = "NM"
# # InPoly = "NM_2019_utm.npz"
# OutStrng = InSurvey

###############################################################################
# Overlap Area
# ##############################################################################
# InSurvey = "A1"
# InPoly = "A1_2019_utm.npz"
# InDatDir = AEMPYX_DATA+"/Blocks/A1/orig/"
# OutStrng = InSurvey+"_NM_intersection"

# InSurvey = "A2"
# InPoly = "A2_2019_utm.npz"
# InDatDir = AEMPYX_DATA+"/Blocks/A2/orig/"
# OutStrng = InSurvey+"_NM_intersection"

# InSurvey = "TB"
# InPoly = "TB_2019_utm.npz"
# InDatDir = AEMPYX_DATA+"/Blocks/TB/orig/"
# OutStrng = InSurvey+"_NM_intersection"

# DataSelect = "Intersection"   # "Polygon", "Intersection", "Union"
# OutDatDir = AEMPYX_DATA+"/RAFTA/Intersection/raw/"
# PolyDir = AEMPYX_DATA+"/RAFTA/Intersection/polygons/"
# PolyFiles = [PolyDir+InPoly,PolyDir+"TNM_2019_utm.npz"]


"""
Output formats are "npz","nc4","asc"
"""



print("Data read from dir:  %s" % InDatDir)
print("Data written to dir: %s" % OutDatDir)
print("Flightline ID string: %s \n" % OutStrng)

if "search" in FileList.lower():
    DataSet = []

    files = os.listdir(InDatDir)
    for entry in files:
        # print(entry)
        if SearchStr in entry.lower():
            DataSet.append(entry)

if "set" in FileList.lower():
    DataSet = []

DataSet = sorted(DataSet)
ns = numpy.size(DataSet)

if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


dcount=0
for dset in DataSet:
    dcount=dcount+1
    start = process_time()
    file = InDatDir + dset
    print("\nRaw data read from: %s" % file)
    Datar = aesys.read_survey_data(DatFile=file, Survey=InSurvey, OutInfo=True)
    if dcount == 1:
        Data = Datar
    else:
        Data = numpy.vstack((Data, Datar))
    print("Read time taken = ", process_time() - start, "s \n")


if SetProj:
    start = process_time()
    itm_e = Data[:,1]
    itm_n = Data[:,2]
    utm_e, utm_n = util.project_itm_to_utm(itm_e, itm_n) #, utm_zone=32629)
    Data[:,1] = utm_e
    Data[:,2] = utm_n
    print("ITM Transformed to UTM")
    print("Projection time taken = ", process_time() - start, "s \n")


"""
Data subsets based on rectangle, polygons or operators on polygons
"""
start = process_time()
print("In: "+str(numpy.shape(Data)))
print("Data select is "+DataSelect+" \n")

if "rec" in DataSelect.lower():
    head = aesys.grow_header(Header, "Subset: " + str(RectCorners))
    start = process_time()
    Rect = util.extract_data_rect(Data, RectCorners)
    if Rect.size != 0:
        Data = Rect
    else:
        error("No data found in rectangle!\n")

if "pol" in DataSelect:
    head = aesys.grow_header(Header,
            " | Subset: " + str(DataSet)+": "+PolyFiles[0])
    Polygon = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
    start = process_time()
    Poly= util.extract_data_poly(Data, Polygon, method="shp")
    if Poly.size != 0:
        Data = Poly
    else:
        error("No data found in polygon!\n")

if ("uni" in DataSelect.lower()) or ("int" in DataSelect.lower()):
    for polyfile in PolyFiles:
        head = aesys.grow_header(Header,
                " | Subset: "+DataSelect.lower()[0:3]+": "+polyfile)
    Polygon1 = numpy.load(PolyFiles[0], allow_pickle=True)["Poly"][0]
    Polygon2 = numpy.load(PolyFiles[1], allow_pickle=True)["Poly"][0]
    start = process_time()
    Polygon= util.modify_polygon([Polygon1, Polygon2], Operator=DataSelect)
    Poly= util.extract_data_poly(Data, Polygon, method="shp")
    if Poly.size != 0:
        Data = Poly
    else:
        error("No data found in polygons!\n")

print("Data select time taken = ", process_time() - start, "s \n")
print("Out: "+str(numpy.shape(Data)))

if MergeOut:
    head = aesys.grow_header(Header,"All Lines")
    f = OutDatDir + OutStrng+"_Full"+OutFileFmt
    aesys.write_aempy(File=f, Data=Data, System=AEM_system,
                    Header=head, OutInfo=OutInfo)
    print("All data written to File: " + f )
    print("Header written: ")
    print(head)
    print("time taken = ", process_time() - start, "s \n")


if LinesOut:
    bad_files = 0
    startlines = process_time()
    Lines = sorted(numpy.unique(Data[:, 0]))
    print(">Flight lines in data set:")
    print(Lines)
    for s in Lines:
        tmp = Data[numpy.where(Data[:, 0] == s), :]
        ns = numpy.shape(tmp)
        tmp = numpy.reshape(tmp, (ns[1], ns[2]))
        print("OutInfo: "+str(numpy.shape(tmp)))

        if numpy.size(tmp)<=nD*LinesMin:
            print("Not enough data! Not written")
            continue

        if CheckNaN:
                nn = numpy.count_nonzero(numpy.isnan(tmp))
                print (str(nn)+" NaNs in Data Block")
                if nn >0:
                    bad_files = bad_files+1
                    print("Too many NaNs = "+str(nn)+" in block, not written")
                    continue


        if CorrectDirection:
            AngLimits = [TellusAng-5., TellusAng+5. ]
            nd =numpy.shape(tmp)[0]
            spoint = [tmp[round(nd*0.3),1], tmp[round(nd*0.3),2]]
            epoint = [tmp[round(nd*0.6),1], tmp[round(nd*0.6),2]]
            ang, _ = util.get_direction_angle(spoint, epoint)
            if (ang < TellusAng-Spread) or (ang > TellusAng+Spread):
                tmp = numpy.flipud(tmp)
                print(" Angle = "+str(round(ang,1))
                    +" not in interval "
                    +str(round(AngLimits[0],1))+" - "
                    +str(round(AngLimits[1],1)))
                print("Flightline direction has been reversed.")
                chdir = ", direction has been reversed"
            else:
                print("Flightline direction is approx. 345 degrees")
                chdir = ""

        head = aesys.grow_header(Header, "Flightline " + str(s))

        f = OutDatDir + OutStrng + "_FL" + str(s).replace(".", "-")+OutFileFmt
        aesys.write_aempy(File=f, Data=tmp, System=AEM_system,
                        Header=head, OutInfo=OutInfo)
        print("Flight line written to File: " + f)
        print("Header written: ")
        print(head)
        print("time taken = ", process_time() - start, "s \n")

    print("Flight line data, time taken = ",
          process_time() - startlines, "s \n")
