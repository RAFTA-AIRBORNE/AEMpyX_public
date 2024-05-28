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

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys


warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

AEM_system = "aem05"
FwdCall,nD,_,_,_, = aesys.get_system_params(System=AEM_system)

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus "
      +"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

Header = "AEMpyX Version "+version
Header = aesys.grow_header(
    Header,
    ["\nTellus ",
     "\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S"))],
)


ReadFilelist = False
LinesOut = True
MergeOut = True
"""
input formats are 'npz','nc4','ascii'
"""
InFileFmt = "asc"
InDatDir = AEMPYX_DATA+"/Blocks/A1O/proc_delete/data_asc/"
# InDatDir = AEMPYX_DATA+"/Blocks/A1O/proc_delete/data_asc/"
print("Data read from dir: %s " % InDatDir)
SearchStrng = "*k3*.asc"
print("Old flightline ID string: %s \n" % SearchStrng)

"""
Output formats are 'npz','nc4','ascii'
"""
OutFileFmt = "asc"
OutDatDir = AEMPYX_DATA+"/RAFTA/stgormans/data_delete/"
print("Data written to dir: %s " % OutDatDir)
OutName = "StGormans_delete_k3"
print("New filname string: %s " % OutName)

if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


Operation = "rectangle"

"""
# Operation = "Rectangle"
# Operation = "Polygon"
# Operation = "Intersect"
# Operation = "Union"
# Operation = "Nearest"
# Operation = "Radius"
# Operation = "ChangeProject"
"""

if "rect" in Operation.lower():
    # RectCorners = [515449.00, 522564.00, 5820187.00, 5827231.00]  # Limerick
    RectCorners = [638968.67, 641519.17, 5922331.93, 5924940.46]  # StGormans
    # RectCorners = []
    OutName = OutName+"_rect"

if "poly" in Operation.lower():
    PolyFiles = ["A1_2019_utm.txt"]
    OutName = OutName+"_poly"

if "union" in Operation.lower() or "intersect" in Operation.lower():
    PolyFiles = ["A1_2019_utm.txt","TNM_2019_utm.txt"]
    OutName = OutName+"_pol_"+Operation.lower()[0:3]

if "nearest" in Operation.lower():
    error("option "+Operation.lower()+"not implemented!. Exit.")
    Points = []
    OutName = OutName+"_near"


if "radius" in Operation.lower():
    error("option "+Operation.lower()+"not implemented!. Exit.")
    Points = []
    SearchRadius = 300
    OutName = OutName+"_rad_"+str(SearchRadius)+"m"

if "project"in Operation.lower():
    error("option "+Operation.lower()+"not implemented!. Exit.")
    InpProj = "latlon"
    OutProj = "utm"



if ReadFilelist:
    dat_files = []
    with open("Filelist.txt", 'r') as file:
        for line in file:
            dat_files.append(line[:-1])
else:
    dat_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InDatDir)

dat_files = sorted(dat_files)
ns = numpy.size(dat_files)

if ns ==0:
    error("No files corresponding to searchstring <"+SearchStrng+"> found!. Exit.")


Data = numpy.empty((1, nD[0]), dtype=float)

for filein in dat_files:
    start = process_time()
    file = InDatDir + filein
    print("\nRaw data read from: %s" % file)
    Datar, header = aesys.read_aempy(File=file, System=AEM_system, OutInfo=False)
    # print(numpy.shape(Datar))
    # print(numpy.shape(Data))
    Data = numpy.vstack((Data, Datar))
    print("Read time taken = ", process_time() - start, "s \n")

"""
Block: Select data subset

Data subsets based on rectangle, polygons or operators on polygons
"""

if ("rect" in Operation.lower()):
    Header = aesys.grow_header(Header, "\nSubset: " + str(RectCorners))
    start = process_time()
    Rect = util.extract_datata_rect(Data, RectCorners)
    if Rect.size != 0:
        Data = Rect
    else:
        print("No data found in rectangle!\n")

    print("Rectangle time taken = ", process_time() - start, "s \n")

if "poly" in Operation.lower():
    Header = aesys.grow_header(Header,
            "\nSubset: " +PolyFiles[0])
    PP = numpy.loadtxt(PolyFiles[0])
    start = process_time()
    Poly= util.extract_datata_poly(Data, PP)
    if Poly.size != 0:
        Data = Poly
    else:
        print("No data found in polygon!\n")

    print("time taken = ", process_time() - start, "s \n")

if "uni" in Operation.lower() or "int" in Operation.lower():
    for polyfile in PolyFiles:
        Header = aesys.grow_header(Header,
                "\nSubset: " +Operation.lower()[0:3]+polyfile)

    P1 = numpy.loadtxt(polyfile[0])
    P2 = numpy.loadtxt(polyfile[1])

    start = process_time()
    Poly= util.extract_datata_poper(Data, P1, P2, Operator=Operation)
    if Poly.size != 0:
        Data = Poly
    else:
        print("No data found in polygons!\n")

    print("time taken = ", process_time() - start, "s \n")


if MergeOut:
    head = aesys.grow_header(Header, "\nAll Lines")
    f = OutDatDir + OutName + "."+OutFileFmt
    aesys.write_aempy(File=f, Format=OutFileFmt, Data=Data, System=AEM_system,
                    Header=head, OutInfo=OutInfo)
    print("All Data written to File: " + f)
    print("Header written: ")
    print(head)
    print("time taken = ", process_time() - start, "s \n")



if LinesOut:
    bad_files = 0
    startlines = process_time()
    Lines = numpy.unique(Data[:, 0])
    for s in Lines:
        tmp = Data[numpy.where(Data[:, 0] == s), :]
        ns = numpy.shape(tmp)
        tmp = numpy.reshape(tmp, (ns[1], ns[2]))

        if numpy.size(tmp)<=nD[0]:
            print("Not enough data! Not written")
            continue

        nn = numpy.count_nonzero(numpy.isnan(tmp))
        print (str(nn)+" NaNs in Data Block")
        if nn >0:
            bad_files = bad_files+1
            print("Too many NaNs = "+str(nn)+" in block, not written")
            continue

        head = aesys.grow_header(Header, "\nFlightline " + str(s))
        f = OutDatDir+OutName+"_FL"+str(s).replace(".", "-")+"."+OutFileFmt
        aesys.write_aempy(File=f, Format=OutFileFmt, Data=tmp, System=AEM_system,
                        Header=head, OutInfo=OutInfo)
        print("Flight line written to File: "+f)
        print("Header written: ")
        print(head)
        print("time taken = ", process_time() - start, "s \n")

    print("Flight line data, time taken = ",
          process_time() - startlines, "s \n")
