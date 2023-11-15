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
import pickle

import numpy
import matplotlib
import shapely.geometry as shg

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0, pth)

from version import versionstrg
import util
import aesys


warnings.simplefilter(action="ignore", category=FutureWarning)


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus find sites near geophysics sites "
      + "".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")



GeophysDir = AEMPYX_DATA+"/Intersection/geophysics/"
DataOutDir = GeophysDir+"/data/"

DataInDir  = AEMPYX_DATA+"/Intersection/area/"

# Area = "Box"
SearchBox = 500;

Area = "Radius"
SearchRadius = 500.

DistExp= 2
InFileFmt = "asc"
OutFileFmt =InFileFmt



DataInFiles = []
ReadFileList = False


###############################################################################
#  aem05
###############################################################################

SearchStrng = "*_NM*"
AEM_system = "aem05"
FwdCall,nD ,_,_,_, = aesys.get_system_params(System=AEM_system)

# ##############################################################################
#   genesis
# ##############################################################################

# SearchStrng = "NM*"
# AEM_system = "genesis"
# FwdCall,nD ,_,_,_, = aesys.get_system_params(System=AEM_system)


GeophysFiles = ["Geophysics_2018_utm_overlaps.npz","Geophysics_2020_utm_overlaps.npz"]


print("Geophysics items read from: %s" % GeophysDir)
print("Data read from to: %s" % DataInDir)
print("Data written to: %s" % DataOutDir)
print("Searchradius is: "+str(SearchRadius)+" m")

if not os.path.isdir(DataOutDir):
    print(" File: %s does not exist, but will be created" % DataOutDir)
    os.mkdir(DataOutDir)


dat_files =DataInFiles
if ReadFileList:
    dat_files = []
    with open(DataInFiles, "r") as file:
        for line in file:
            dat_files.append(line[:-1])
    dat_files = sorted(dat_files)
    ns = numpy.size(dat_files)
    if ns == 0:
        error("No files corresponding found in "+DataInFiles+"!. Exit.")
    else:
        print(str(ns)+" data files found in "+DataInFiles)
else:
    dat_files = util.get_filelist(searchstr=[SearchStrng], searchpath=DataInDir)
    dat_files = sorted(dat_files)
    ns = numpy.size(dat_files)
    if ns ==0:
        error("No files corresponding to searchstring <"+SearchStrng+"> found!. Exit.")
    else:
        print(str(ns)+" data files found in "+DataInDir)

for datafile in dat_files:

    data, header = aesys.read_aempy(File=DataInDir+datafile, Format=InFileFmt,
                                   System=AEM_system, OutInfo=False)
    now = datetime.now()
    header = aesys.grow_header(header, ["\nProcessed " + now.strftime("%m/%d/%Y, %H:%M:%S")])

    mindist_all = []
    nearest = numpy.empty_like(data[0,:])

    # name_in, ext = os.path.splitext(datafile)
    # outdir = DataOutDir+name_in+"-"+str(int(SearchRadius))
    # if not os.path.isdir(outdir):
    #     print("Directory %s does not exist, but will be created" % outdir)
    #     os.mkdir(outdir)

    for geopfile in GeophysFiles:
        name_in, ext = os.path.splitext(geopfile)
        # outdir = DataOutDir+name_in+"-"+str(int(SearchRadius))
        # if not os.path.isdir(outdir):
        #     print("Directory %s does not exist, but will be created" % outdir)
        #     os.mkdir(outdir)

        geophys_list =numpy.load(GeophysDir+geopfile, allow_pickle=True)["Data"]

        for line in geophys_list:

            outdata = []
            distance = []
            outnear = []
            numdata = 0

            outdir = DataOutDir+"/Geophysics_"+line[0]+"/"
            if not os.path.isdir(outdir):
                print(" File: %s does not exist, but will be created" % outdir)
                os.mkdir(outdir)

            outfile_data = outdir+AEM_system.upper()+"_Geophysics"+line[0]+"_Datafile"+"_SearchRadius"+str(round(SearchRadius))+"m"
            outfile_best = outdir+AEM_system.upper()+"_Geophysics"+line[0]+"_best"

            h = aesys.grow_header(header, ["\n"+str(line)])

            x0 = line[1]
            y0 = line[2]
            for site in data:
                x1 = site[1]
                y1 = site[2]
                if "radius" in Area.lower():
                    r = numpy.sqrt((x1-x0)**2+(y1-y0)**2)

                    if r <= SearchRadius:
                        site[nD[0]-1]=r
                        distance.append(r)
                        outdata.append(site)

                if "box" in Area.lower():
                    if (x1 >= x0-SearchBox and y1 >= y0-SearchBox and
                        x1 <= x0+SearchBox and y1 <=y0+SearchBox):
                        site[nD[0]-1]=r
                        outdata.append(site)
                        distance.append(r)


            outdata = numpy.asarray(outdata, dtype=float)
            distance =  numpy.asarray(distance, dtype=float)
            numdata = numpy.shape(outdata)[0]

            if numdata ==0:
                print("No data for "+line[0]+" found!")
                continue
            else:
                flines=numpy.unique(outdata[:,0])
                print(numpy.unique(outdata[:,0]))
                print(str(numdata)+" nearest data written to: %s" % outfile_data)
                f = outfile_data

                aesys.write_aempy(File=f, Data=outdata, Format=OutFileFmt,
                                   System=AEM_system, Header=header, OutInfo=False)

                near = []
                mindist = numpy.min(distance)

                for ll in numpy.arange(numpy.shape(outdata)[0]):
                    if numpy.isclose(distance[ll], mindist, atol =0.1):
                       outnear.append(outdata[ll,:])

                outbest = numpy.asarray(outnear, dtype=float)


                outave = numpy.nanmean(outdata,axis=0).reshape(1,nD[0])
                outave[0,0]=0.
                outbest = numpy.append(outbest, outave,axis=0)

                outmed = numpy.nanmedian(outdata,axis=0).reshape(1,nD[0])
                outmed[0,0]=0.
                outbest = numpy.append(outbest,outmed,axis=0)

                # Inverse distance
                wi = 1./distance**DistExp
                wd = numpy.nansum(wi)
                w=(wi/wd).reshape(numpy.shape(outdata)[0],1)
                outidw =numpy.nansum(w*outdata,axis=0).reshape(1,nD[0])
                outidw[0,0]=0.
                outidw[0,1]=x0
                outidw[0,2]=y0
                outbest = numpy.append(outbest,outidw,axis=0)
                f = outfile_best #+OutFileFmt
                aesys.write_aempy(File=f, Data=outbest, Format=OutFileFmt,
                                  System=AEM_system, Header=header, OutInfo=True)
                print("Best site data written to: %s" % outfile_best)


                mindist_all.append(mindist)




    # n, bins, patches = matplotlib.pyplot.hist(mindist_all,20, density=False, facecolor="g", alpha=0.75)

    # matplotlib.pyplot.xlabel("distance from borehole (m)")
    # matplotlib.pyplot.ylabel("n")
    # matplotlib.pyplot.title("Nearest points to borehole: "+datafile)
    # matplotlib.pyplot.grid(True)
    # matplotlib.pyplot.show()
