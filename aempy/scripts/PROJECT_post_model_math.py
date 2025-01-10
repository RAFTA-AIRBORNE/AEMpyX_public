#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 10:26:24 2024

@author: vrath
"""

import os
import sys
from sys import exit as error
from datetime import datetime
from time import process_time, time
# import inspect
# import copy
from datetime import datetime
import warnings
import getpass
import getpass

import numpy
import scipy

# import multiprocessing
# from numba import njit

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse
# -

AEMPYX_DATA = os.environ["AEMPYX_DATA"]
rng = numpy.random.default_rng()
nan = numpy.nan
dateform="%m/%d/%Y, %H:%M:%S"

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = False


# AEM_system = "genesis"
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 75. #50.
    DatErr_mult = 0.0
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0] = 0   # real at 900Hz
    data_active[4] = 0   # imag at 900Hz

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.0
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' horizontals'


Method = "twoway"  # "mean", "weighted_mean"


"""
Input
"""
##############################################################################
# StGormans
##############################################################################
AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/A1_StGormans/"
# InResDir =  AEMPYX_DATA + "/proc/"
InResDir =  AEMPYX_DATA + "/results_rect/"
if not InResDir.endswith("/"): InResDir=InResDir+"/"
print("Data read from dir: %s " % InResDir)
# +

"""
Output"
"""
# FileList = "set"
FileList = "search"  # "search", "read"
SearchStrng = "*FL*_normal*results.npz"
# OutResDir =   AEMPYX_DATA + "/results_parallel/"
OutResDir =   AEMPYX_DATA + "/results_rect/"
if not OutResDir.endswith("/"): OutResDir=OutResDir+"/"
print("Models written to dir: %s " % OutResDir)
if not os.path.isdir(OutResDir):
    print("File: %s does not exist, but will be created" % OutResDir)
    os.mkdir(OutResDir)



if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InResDir)
    # res_files = []
    res_files = [InResDir+"StGormans_FL11379-0_raw.npz"]
    # res_files =  numpy.load(AEMPYX_DATA + "/Projects/Compare/BundoranSubsets.npz")["setC"]
    res_files = [os.path.basename(f) for f in res_files]
else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    res_files = util.get_data_list(how=["search", SearchStrng, InResDir],
                              fullpath=False, out= True, sort=True)

ns = numpy.size(res_files)
if ns ==0:
    error("No files set!. Exit.")


if  "two" in Method.lower():

    for f in res_files:

        filein_normal = f
        print("\n")
        print("Normal models read from: %s" % InResDir+filein_normal)
        f1  = numpy.load(InResDir+filein_normal, allow_pickle=True)
        m = f1["site_modl"]
        if ParaTrans==1:
           m = numpy.log(m)
        m_normal = m.copy()

        filein_revers = f.replace("_normal", "_reverse")
        print("Reverse models read from: %s" % InResDir+filein_normal)
        f2  = numpy.load(InResDir+filein_revers, allow_pickle=True)
        m = f2["site_modl"]
        if ParaTrans==1:
           m = numpy.log(m)
        m_revers = m.copy()

        m_avg = 0.5*(m_normal + m_revers)
        m_dif = numpy.abs(m_normal + m_revers)
        if ParaTrans==1:
           m = numpy.exp(m)

        d = dict(f1)
        d["site_modl"] = m_avg
        d["site_diff"] = m_dif
        d["header"] = "Averaged model set:"+"".join("Date " + datetime.now().strftime(dateform))
        fileout = filein_normal.replace("_normal", "_average")
        print("Averaged models saved to: %s" % OutResDir+fileout)
        numpy.savez_compressed(file=OutResDir+fileout, **d)


# if "mean" in Method.lower:

#  for f in res_files:


#      filein = InResDir+f
#      print("\nNModels read from: %s" % filein)

#      m  = numpy.load(filein, allow_pickle=True)
#      model = m["mod"]
#      if ParaTrans==1:
#         m = numpy.log10(m)
#      if "weight" in  Method.lower():
#         rms = m["nrms"]
