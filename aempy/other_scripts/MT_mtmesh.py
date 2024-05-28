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
Show several 1d block models as (stitched) section.

"""
import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings

import numpy

# import matplotlib.collections
# import matplotlib.patches
# import matplotlib
# import matplotlib.colors
# import matplotlib.pyplot
# import matplotlib.cm

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)

from version import versionstrg


import util
# import viz
import inverse
from mt import mt1dfwd


warnings.simplefilter(action="ignore", category=FutureWarning)

AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
cm = 1/2.54  # centimeters in inches

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = True
now = datetime.now()

MTDataType = ["rhophas", "imped"]
Freqs = [3.e4, 1.e4, 3.e3, 1.e3, 3.e2]


"""
input formats are "npz","nc4","asc"
"""
InFileFmt = ".npz"
InModDir = AEMPYX_DATA +"/Projects/Munster/results/"

print("Data read from dir:  %s" % InModDir)

FileList = "search"  # "search", "read"

if "search" in FileList.lower():
    SearchStrng = "*Results.npz"
    print("Search flightline ID string: %s " % SearchStrng)
    data_files = util.get_filelist(searchstr=[SearchStrng],
                                   searchpath=InModDir,
                                   fullpath=False)
    data_files = sorted(data_files)

if "set" in FileList.lower():
    data_files =["NM_A1_intersection_FL13490-0_cnlyr30_TikhOpt_gcv_Results.npz"]

ns = numpy.size(data_files)

print(data_files)

MTOutDir = AEMPYX_DATA +"/Projects/Munster/mt_from_aem/"
print("mt data written to dir: %s " % MTOutDir)
if not os.path.isdir(MTOutDir):
    print("File: %s does not exist, but will be created" % MTOutDir)
    os.mkdir(MTOutDir)

allsites = 0
for filein in data_files:
    start = process_time()
    print("\nData read from: %s" % filein)

    name, ext =os.path.splitext(filein)

    aem = numpy.load(InModDir+filein)
    easting  = aem["site_x"]
    northing = aem["site_y"]
    site_nrms = aem["site_nrms"]
    aemmodel = aem["site_modl"]
    refmodel = aem["mod_ref"]

    nlayer = inverse.get_nlyr(refmodel)
    thklayer = refmodel[6*nlayer:7*nlayer-1]

    [nsite, nvals] = numpy.shape(aemmodel)
    allsites = allsites + nsite
    for isite in numpy.arange(nsite):
        rholayer = aemmodel[isite,:]
        imp, rhoa, phas = mt1dfwd(Freqs, rholayer, thklayer,
                             inmod="res", outdat="both")
        if isite ==0:
            site_imp = imp
            site_rhoa = rhoa
            site_phas = phas
        else:
            site_imp = numpy.vstack((site_imp,imp))
            site_rhoa = numpy.vstack((site_rhoa, rhoa))
            site_phas = numpy.vstack((site_phas, phas))


    fileout  = MTOutDir+name.replace("Results","MT")+ext
    print(fileout)
    numpy.savez_compressed(fileout,
                           site_x=easting,
                           site_y=northing,
                           aem_model=aemmodel,
                           thklayer=thklayer,
                           freq=Freqs,
                           site_imp=site_imp,
                           site_rhoa=site_rhoa,
                           site_phas=site_phas,
                           site_nrms= site_nrms
                          )


elapsed = (process_time() - start)
print (" Used %7.4f sec for %6i files (%7i sites)" % (elapsed, ns, allsites))
print (" Average %7.4f sec/site\n" % (elapsed/(allsites)))

print("\n\nAll done!")
