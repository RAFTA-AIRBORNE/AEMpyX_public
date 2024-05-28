#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 16:39:41 2022

@author: vrath
"""
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
#       format_version: "1.5"
#       jupytext_version: 1.11.4
# ---

"""
Created on Tue Sep  6 10:57:01 2016

@author: vrath

edited by dkiyan - Sep 30
edited by vrath  - May 7, 2021

"""
import time
import sys
from sys import exit as error
import os
import warnings
from time import process_time
from datetime import datetime
import math

import numpy
from cycler import cycler
import matplotlib
import matplotlib.pyplot

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import util
import core1d
import inverse
import aesys
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]


print(" AEMpyX search for nearest flight lines")

nan = numpy.nan  # float("NaN")
rng = numpy.random.default_rng()


version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("AEMpyX search for nearest flight lines"+"\n"
      +"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

"""
System related settings.

"""
AEM_system1 = "genesis"
FwdCall1,NN1, _, _, _, = aesys.get_system_params(System=AEM_system1)
nD1 = NN1[0]

AEM_system2 = "aem05"
FwdCall2,NN2, _, _, _, = aesys.get_system_params(System=AEM_system2)
nD2 = NN2[0]

"""
input formats are "npz","nc4","asc"
"""
InFileFmt = ".asc"
InDatDir = AEMPYX_DATA + "/Nearest/"

Block="TB"
SearchStrng1 = "NM_"+Block+"*.asc"
SearchStrng2 = Block +"_NM"+"*.asc"
Fileout = InDatDir+"/TD_"+Block+"_Nearest_FDLines.txt"

start = time.time()


print("Data read from dir:  %s" % InDatDir+"/raw/")
print("Search flightline ID string: %s " % SearchStrng1)
data_files1 = util.get_filelist(searchstr=[SearchStrng1], searchpath=InDatDir+"/raw/")
data_files1 = sorted(data_files1)

ifl = -1
for file in data_files1:

    ifl = ifl+1
    name, ext = os.path.splitext(file)
    filein =InDatDir +"/raw/" + file
    Data, Header = aesys.read_aempy(File=filein, Format=InFileFmt,
                                   System=AEM_system1, OutInfo=False)

    sD = numpy.shape(Data)
    print("flightline "+name+"  #"
          +str(ifl)+" of "
          +str(numpy.size(data_files1)) +" has shape: "+str(sD))

    nc = round(sD[0]/2)
    if ifl ==0:
        fl_numc = Data[0,0]
        pcenter = [Data[nc,1],Data[nc,2]]
    else:
        fl_numc = numpy.vstack((fl_numc, Data[0,0]))
        pcenter =numpy.vstack((pcenter, [Data[nc,1],Data[nc,2]]))

npoints = numpy.shape(pcenter)[0]


print("Data read from dir:  %s" % InDatDir)
print("Search flightline ID string: %s " % SearchStrng2)
data_files2 = util.get_filelist(searchstr=[SearchStrng2], searchpath=InDatDir+"/raw/")
data_files2 = sorted(data_files2)

ifl = -1
for file in data_files2:

    ifl = ifl+1
    name, ext = os.path.splitext(file)
    filein =InDatDir +"/raw/" + file
    Data, Header = aesys.read_aempy(File=filein, Format=InFileFmt,
                                   System=AEM_system2, OutInfo=False)

    sD = numpy.shape(Data)
    print("flightline "+name+"  #"
          +str(ifl)+" of "
          +str(numpy.size(data_files1)) +" has shape: "+str(sD))

    n1 = 0
    n2 = sD[0]-1
    if ifl ==0:
        fl_num = Data[0,0]
        p1 = [Data[n1,1],Data[n1,2]]
        p2 = [Data[n2,1],Data[n2,2]]

    else:
        fl_num = numpy.vstack((fl_num, Data[0,0]))
        p1 =numpy.vstack((p1, [Data[n1,1],Data[n1,2]]))
        p2 =numpy.vstack((p2, [Data[n2,1],Data[n2,2]]))

nsegments = numpy.shape(p1)[0]

print("\n\n")
print(" Number of center points: "+str(npoints))
print(" Number of lines: "+str(nsegments))
"""
Distance from p0 perpendicular to a line drawn between p1 and p2.
with: p1=(x1,y1), p2=(x2,y2), and p0=(x0,y0)
d = numpy.abs(numpy.lininverse.norm(numpy.cross(p2-p1, p1-p0))/norm(p2-p1))

https://stackoverflow.com/questions/39840030/distance-between-point-and-a-line-from-two-points
"""

icc = -1
for pp in numpy.arange(npoints):
    icc=icc+1
    p  = pcenter[pp]
    for ll in numpy.arange(nsegments):
        pl1 = p1[ll]
        pl2 = p2[ll]
        if ll==0:
            d = numpy.lininverse.norm(numpy.cross(pl2-pl1, pl1-p)/numpy.lininverse.norm(pl2-pl1))
            dist = d
        else:
            d = numpy.lininverse.norm(numpy.cross(pl2-pl1, pl1-p)/numpy.lininverse.norm(pl2-pl1))
            dist = numpy.hstack((dist,d))

    nmin = numpy.argmin(dist)
    nearest_line = data_files2[nmin]

    if icc==0:
        nearest = numpy.array([data_files1[pp], nearest_line, Block ], dtype=object)
    else:
        nearest =numpy.vstack((nearest, [data_files1[pp], nearest_line, Block ]))


print(" Results stored to "+Fileout)
numpy.savetxt(Fileout,nearest, fmt="%s    %s    %s")
elapsed = (time.time() - start)
print (" Used %7.4f sec for %6i nearest sites \n" % (elapsed, npoints))
