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
Created on Sun Feb 28 17:12:28 2021

@author: vrath
"""

import os
import sys
import warnings
from time import process_time
from datetime import datetime

import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.append(pth)

import util
import aesys
from version import versionstrg


warnings.simplefilter(action="ignore", category=FutureWarning)

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus"+"\n"
      +"".join("Date "+now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

AEMPYX_DATA = os.environ["AEMPYX_DATA"]

InFileFmt = "asc"
InDatDir = AEMPYX_DATA+"_work/data/A5/"
Name = "A5_All"
print("Data read from : %s " % InDatDir+Name)


# XLimits=[515449.00, 522564.000]
# dX=300
# YLimits=[5820187.00, 5827231.00]
# dY=300
# AddStr = "_LIM""

XLimits = []
YLimits = []
AddStr = ""

Header = ("AEMpyX Version "+version
          +"Tellus"+"\n"
          +"".join("Date "+now.strftime("%m/%d/%Y, %H:%M:%S"))
          +"\n"+Name)

start = process_time()
Data, _ = aesys.read_aempy_asc(InDatDir + Name + ".asc")
Points = Data[:, 1:3]
print("data read, time taken = ", process_time() - start, "s \n")
# print(numpy.shape(Points))

if XLimits == []:
    dX = 300
    XLimits = [numpy.min(Points[:, 0]), numpy.max(Points[:, 0])]
if YLimits == []:
    dY = 300
    YLimits = [numpy.min(Points[:, 1]), numpy.max(Points[:, 1])]

Header = aesys.grow_header(Header, "\nUTM"
                         +"\nX: "+str(XLimits[0])+","
                         +str(XLimits[1])+",    dX = "+str(dX)
                         +"\nY: "+str(YLimits[0])+","
                         +str(YLimits[1])+",    dX = "+str(dY)
                         )

start = process_time()
p = util.gen_searchgrid(Points,
                       XLimits=XLimits, dX=dX, YLimits=YLimits, dY=dY)
print("grid ready,  time taken = ", process_time() - start, "s \n")

Fname = InDatDir+Name+AddStr+'_searchgrid.npz'
numpy.savez_compressed(Fname,
                    Header=Header,
                    dX=dX, XLimits=XLimits, dY=dY, YLimits=YLimits,
                    Searchgrid=p)
print("Search mesh written to: %s " % Fname)
