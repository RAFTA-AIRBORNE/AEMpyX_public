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
#       jupytext_version: 1.15.2
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

# action = "decimate"
action = "split"
# action = "cut"

if "dec" in action.lower():
    Window = 3
    Method = ["mean", Window]

if "spl" in action.lower():
    Step = 3
    Method = [Step]

if "cut" in action.lower():
    Interval = [1000., 2000.]
    interval = sorted(Interval)
    


"""
input formats are 'npz','nc4','ascii'
"""

FileList = "search"  # "set", "read"
# InDatDir = AEMPYX_DATA+"/Projects/TestLines/GENESIS/raw/"
InDatDir = "/home/vrath/AEM_Data/Blocks/A1/proc_delete_PLM10/k3split/"
SearchStrng = "*.npz"

print("Data read from dir: %s " % InDatDir)

if "set" in FileList.lower():
    flightlines = []
else:
    flightlines = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(flightlines)
    if ns ==0:
        error("No files corresponding to searchstring <"+SearchStrng+"> found!. Exit.")


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

    if "cut" in action:
        if len(flightlines)==1:
            profile = prep.get_profile(Data[:, 1], Data[:, 1])
            interval[0] = numpy.amax([profile[0],interval[0]])
            interval[1] = numpy.amin([profile[-1],interval[1]])
            intvl = numpy.where(numpy.logical_and((profile>=interval[0]),
                                                  (profile<=interval[1])))
            Dcut = D[intvl]
            print("\n Proc action: " + action)
            print(" interval ", interval)
            Header = aesys.grow_header(Header, filename +", CUT = "+str(interval))
            print(" data block now has shape: ", numpy.shape(Dcut))
            newname = name+"_cut_"+str(interval[0])+"-"+str(interval[1])
            filout = OutDatDir+newname+OutFileFmt
            head = Header
            aesys.write_aempy(File=filout, Data=Dcut,
                            Header=head, OutInfo=False)

            print("Decimated data  written to File: " + filout)
            print("Info:")
            print(head)
        else:
            error(" cut option should only be used for single flight line! Exit.")

    if "dec" in action:
        if Window > 1:
            print("\n Proc action: " + action)
            print(" method: ", Method[:])
            Header = aesys.grow_header(Header, "DECIMATE = "+str(Method))
            D, blkhead = prep.reduce_data(D, Method=Method, System = AEM_system)
            print(" data block now has shape: ", numpy.shape(D))

            newname = name + "_dec" + str(Window)
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

                newname = name+"_spl" + str(Step)+"-"+str(start)
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
