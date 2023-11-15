#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadat_filter: -all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: "1.5"
#       jupytext_version: 1.11.4
# ---

    # ! /usr/bin/python
# ---
# jupyter:
#   jupytext:
#     cell_metadat_filter: comment_questions,-all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: "1.5"
#       jupytext_version: 1.11.3
# ---

import os
import sys
from sys import exit as error
from datetime import datetime
from time import process_time
from random import randrange
import time
import warnings

import numpy
import scipy

# %logstart -o

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse

warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Inversion Domain-of-Investigation"+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


OutInfo = True
now = datetime.now()

Method = "Oldenburg1999"
Method = "Variance"

"""
input and output format is ".npz",
"""
InModDir ="/home/vrath/work/Clara/work/AEM_FD/" # AEMPYX_DATA + "/Blocks/StGormans/proc_delete/"
# InDatDir = AEMPYX_DATA + "/Intersection/geophysics/data/Geophysics_R12A/"
# InDatDir = AEMPYX_DATA + "/Tests/"
print("Models read from dir:  %s" % InModDir)
OutModDir = InModDir
print("Models written to dir: %s " % OutModDir)
if not os.path.isdir(OutModDir):
    print("File: %s does not exist, but will be created" % OutModDir)
    os.mkdir(OutModDir)

FileList = "set"  # "search", "read"

if "set" in FileList.lower():
	mod_files = [
	"A1_NM_intersection_FL11126-0_delete_PLM10_k3_prof1_nlyr32_TikhOpt_gcv_Prior30_Results.npz",
	"A1_NM_intersection_FL11126-0_delete_PLM10_k3_prof1_nlyr32_TikhOpt_gcv_Prior3000_Results.npz"
				  ]


if "search" in FileList.lower():
    mod_files = []
    with open("Filelist.txt", "r") as file:
        for line in file:
            mod_files.append(line[:-1])


if "search" in FileList.lower():
    SearchStrng = "StGormans_*k5*.asc"
    print("Search string: %s " % SearchStrng)
    mod_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InModDir)
	mod_files = sorted(mod_files)
	ns = numpy.size(mod_files)
	if ns ==0:
		error("No files found!. Exit.")

Fileout = OutModDir +"AEM05_TikhOpt_P30-3000_DoI.npz"
Header = "DoI from FD models"

start = time.time()

if "olden" in Method.lower():
    icount = -1
    for file in mod_files:
        icount = icount+1

        m = numpy.load(OutModDir +file)
        m_ref = m["mod_ref"]
        m_act = m["mod_act"]
        nsit = len(m["site_num"])
        nlyr = inverse.get_nlyr(m_ref)
        prior = inverse.extract_mod(M=m_ref, m_act=m_act)


        if icount == 0:
            modls = m["site_modl"]
            merrs = m["site_merr"]
            mrefs = prior
            modshp0 = numpy.shape(m["site_modl"])
            # print( numpy.shape(modls))
        else:
            modls = numpy.vstack((modls, m["site_modl"]))
            merrs = numpy.vstack((merrs, m["site_merr"]))
            mrefs = numpy.vstack((mrefs, prior))
            modshp = numpy.shape(m["site_modl"])
            if modshp != modshp0:
                error("model dimenssions do not fit! Exit.")
            # print( numpy.shape(modls))


    modls = numpy.log10(modls.reshape(ns, nsit, nlyr))
    mrefs = numpy.log10(mrefs.reshape(ns, nlyr))

    # print("modls ", numpy.shape(modls))
    m_avg = numpy.sum(modls, axis=0).reshape(nsit,nlyr)
    m_dif = numpy.diff(modls, axis=0).reshape(nsit,nlyr)
    r_dif = numpy.diff(mrefs, axis=0)
    # print("avg ", numpy.shape(m_avg))
    # print("mdif", numpy.shape(m_dif))
    # print("rdif", numpy.shape(r_dif))
    doi = numpy.abs(m_dif)/numpy.abs(r_dif)
    # print("doi", numpy.shape(doi))

    print("\n\nMethod:  "+Method)
    print("DoI min = "+str(numpy.amin(doi)))
    print("    max = "+str(numpy.amax(doi)))
    numpy.savez_compressed(
       file=Fileout,
       mod_files=mod_files,
       header=Header,
    )

if "var" in Method.lower():
    icount = -1
    for file in mod_files:
        icount = icount+1

        m = numpy.load(OutModDir +file)
        m_ref = m["mod_ref"]
        m_act = m["mod_act"]
        nsit = len(m["site_num"])
        nlyr = inverse.get_nlyr(m_ref)
        prior = inverse.extract_mod(M=m_ref, m_act=m_act)


        if icount == 0:
            modls = m["site_modl"]
            merrs = m["site_merr"]
            mrefs = prior
            modshp0 = numpy.shape(m["site_modl"])
            # print( numpy.shape(modls))
        else:
            modls = numpy.vstack((modls, m["site_modl"]))
            merrs = numpy.vstack((merrs, m["site_merr"]))
            mrefs = numpy.vstack((mrefs, prior))
            modshp = numpy.shape(m["site_modl"])
            if modshp != modshp0:
                error("model dimenssions do not fit! Exit.")
            # print( numpy.shape(modls))

    modls = numpy.log10(modls.reshape(ns, nsit, nlyr))
    mrefs = numpy.log10(mrefs.reshape(ns, nlyr))

    m_avg = numpy.mean(modls, axis=0).reshape(nsit,nlyr)
    m_ano = numpy.zeros_like(modls)
    for nn in numpy.arange(ns):
        m_ano[nn,:,:] = modls[nn,:,:] - m_avg

    m_var = numpy.var(m_ano,axis = 0)
    doi = numpy.sqrt(m_var)
    print("\n\nMethod:  "+Method)
    print("DoI min = "+str(numpy.amin(doi)))
    print("    max = "+str(numpy.amax(doi)))


print("Results stored to "+Fileout)
elapsed = (time.time() - start)
print (" Used %7.4f sec \n" % (elapsed))

print("\nAll done!")
