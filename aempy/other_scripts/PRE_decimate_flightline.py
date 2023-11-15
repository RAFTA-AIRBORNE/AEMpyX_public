#!/usr/bin/env python3
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
Created on 2020/11/11

This script present a work flow for prrocessing data for Inversions with
INV.py scripts

@author: vrath Oct 2020
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import inverse
import util
import prep
import aesys
import viz

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

now = datetime.now()

AEM_system = "genesis"
#AEM_system = "aem05"

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 50.
    DatErr_mult = 0.03
    alt = 60.
    data_active = numpy.ones(NN[2], dtype="int8")


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 100.
    DatErr_mult = 0.01
    alt = 90.
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontal

action = "decimate" # "
action = "split"
if "dec" in action:
    Window = 3
    Method = ["mean", Window]

if "spl" in action:
    Step = 3
    Method = [Step]


"""
input formats are 'npz','nc4','ascii'
"""

FileList = "search"  # "set", "read"
InFileFmt = ".npz"
#InDatDirprint(" New flightline ID string: %s " % OutStrng) = AEMPYX_DATA+"/Blocks/A1/raw/"
InDatDir = AEMPYX_DATA+"/Projects/TestLines/GENESIS/raw/"
InStrng = "CGGT"

SearchStrng = "CGGT_FL*.npz"
print("Data read from dir: %s " % InDatDir)
print("Old flightline ID string: %s \n" % InStrng)


if "search" in FileList.lower():
    flightlines = util.get_filelist(searchstr=[SearchStrng], searchpath=InDatDir)
    ns = numpy.size(flightlines)

if "set" in FileList.lower():
	flightlines = []

OutFileFmt = ".npz"
OutDatDir = InDatDir
print("Data written to dir: %s " % OutDatDir)

if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)

ns = numpy.size(flightlines)
flightlines = sorted(flightlines)

start = process_time()
nsites = 0
nfiles = 0
for filename in flightlines:
    nfiles = nfiles+1
    name, ext = os.path.splitext(filename)
    filein = InDatDir + filename
    print("\n Preprocessing file " + filein)
    Data, Header, _ = aesys.read_aempy(File=filein,
                                   System=AEM_system, OutInfo=False)
    now = datetime.now()
    Header = aesys.grow_header(Header, titstrng)

    fline = Data[:, 0]
    D = Data[:, :]
    sizedat = numpy.shape(D)
    nvars = sizedat[1]
    last = nvars - 1
    print("Flightline Data block on input has shape: ", numpy.shape(D))



    if "dec" in action:
        if Window > 1:

            print("\n Proc action: " + action)
            print(" method: ", Method[:])
            Header = aesys.grow_header(Header, "DECIMATE = " + str(Method))
            D, blkhead = prep.reduce_data(D, Method=Method, System = AEM_system)
            print(" data block now has shape: ", numpy.shape(D))

            OutStrng = InStrng+"_dec" + str(Window)
            newname = name.replace(InStrng, OutStrng)+"_spl"+str(start)
            filout = OutDatDir+newname+OutFileFmt
            head = Header
            aesys.write_aempy(File=filout, Data=D,
                            Header=head, OutInfo=False)
            print("Decimated data  written to File: " + filout)
            print("Info:")
            print(head)
        else:
            error("Window is 1, no decimation! Exit.")

    if "spl" in action:
        # error("Splitting not yet implemented! Exit.")
        if Step > 1:
            for start in numpy.arange(Step):
                Dsplit=D[start:-1:Step]
                print("\n Proc action: " + action)
                print(" method: ", Method[:])
                Header = aesys.grow_header(Header, "SPLIT = " + str(Method)+" "+str(start))
                print(" data block now has shape: ", numpy.shape(Dsplit))
                OutStrng = InStrng+"_spl" + str(Step)
                newname = name.replace(InStrng, OutStrng)+"_spl"+str(start)
                filout = OutDatDir+newname+OutFileFmt
                head = Header
                aesys.write_aempy(File=filout, Data=Dsplit,
                                Header=head, OutInfo=False)

                print("Decimated data  written to File: " + filout)
                print("Info:")
                print(head)

        else:
            error("Step is 1, no decimation! Exit.")

elapsed = process_time() - start
print(" Used %7.4f sec for %8i files\n" % (elapsed, nfiles))

#
