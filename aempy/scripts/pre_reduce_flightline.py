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

import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import inverse
import util
import prep
import aesys


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

now = datetime.now()

# AEM_system = "genesis"
AEM_system = "aem05"

# 
Action = "decimate"
# Action = "split"
# Action = "cut"

if "dec" in Action.lower():
    # Window = 5
    # Method = ["median", Window]

    Window = 5
    Method = ["mean", Window]

if "spl" in Action.lower():
    Step = 5
    Method = [Step]

if "cut" in Action.lower():
    Interval = [1000., 2000.]
    interval = sorted(Interval)
    


"""
input formats are 'npz','nc4','ascii'
"""

OutFileFmt = ".npz"
FileList = "search"  # "search", "read"


InDatDir = AEMPYX_ROOT+"/work/data/proc/"
OutDatDir = AEMPYX_ROOT+"/work/data/dec/"

SearchStrng = "*_k2.npz"


print("Data read from dir:  %s" % InDatDir)
print("Data written to dir: %s" % OutDatDir)
print("SearchStrng is %s\n" % SearchStrng)


if not InDatDir.endswith("/"): InDatDir=InDatDir+"/"
if not OutDatDir.endswith("/"): InDatDir=OutDatDir+"/"
if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = []
    
    
    dat_files = [os.path.basename(f) for f in dat_files]  
else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True, fullpath=False)
ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")



start = process_time()
nsites = 0
nfiles = 0
for filename in dat_files:
    nfiles = nfiles+1
    name, ext = os.path.splitext(filename)
    filein = InDatDir+filename
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

    if "cut" in Action:
        if len(dat_files)==1:
            profile = prep.get_profile(Data[:, 1], Data[:, 1])
            interval[0] = numpy.amax([profile[0],interval[0]])
            interval[1] = numpy.amin([profile[-1],interval[1]])
            intvl = numpy.where(numpy.logical_and((profile>=interval[0]),
                                                  (profile<=interval[1])))
            Dcut = D[intvl]
            print("\n Proc Action: " + Action)
            print(" interval ", interval)
            Header = aesys.grow_header(Header, filename +", CUT = "+str(interval))
            print(" data block now has shape: ", numpy.shape(Dcut))
            newname = name+"_cut_"+str(interval[0])+"-"+str(interval[1])
            filout = OutDatDir+newname+OutFileFmt
            head = Header
            aesys.write_aempy(File=filout, Data=Dcut,
                            Header=head, OutInfo=False)

            print("Reduced data  written to File: " + filout)
            print("Info:")
            print(head)
        else:
            error(" cut option should only be used for single flight line! Exit.")

    if "dec" in Action:
        if Window > 1:
            print("\n Proc Action: " + Action)
            print(" method: ", Method[:])
            Header = aesys.grow_header(Header, "DECIMATE = "+str(Method))
            D, blkhead = prep.reduce_data(D, Method=Method, System = AEM_system)
            print(" data block now has shape: ", numpy.shape(D))

            newname = name+"_dec"+str(Window)+"_"+Method[0]
            filout = OutDatDir+newname+OutFileFmt
            head = Header
            aesys.write_aempy(File=filout, Data=D,
                            Header=head, OutInfo=False)
            print("Decimated data  written to File: " + filout)
            print("Info:")
            print(head)
        else:
            error("Window is 1, no decimation! Exit.")

    if "spl" in Action:
        # error("Splitting not yet implemented! Exit.")
        if Step > 1:
            for start in numpy.arange(Step):
                Dsplit=D[start:-1:Step]

                print("\n Proc Action: " + Action)
                print(" method: ", Method[:])
                Header = aesys.grow_header(Header, "SPLIT = " + str(Method)+" "+str(start))
                print(" data block now has shape: ", numpy.shape(Dsplit))

                newname = name+"_spl"+str(Step)+"-"+str(start)
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
