#!/usr/bin/env python3
# ---
"""
"""
import os
import sys
from sys import exit as error
from datetime import datetime
from time import process_time
from random import randrange
import time
import warnings
from cycler import cycler

import numpy
import scipy

import matplotlib
import matplotlib.pyplot
import seaborn


AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Plot model-like inversion output"+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


OutInfo = True
now = datetime.now()

"""
input format is ".npz",
"""
InModDir ="/home/vrath/work/Clara/work/AEM_FD/" # AEMPYX_DATA + "/Blocks/StGormans/proc_delete/"
# InDatDir = AEMPYX_DATA + "/Intersection/geophysics/data/Geophysics_R12A/"
# InDatDir = AEMPYX_DATA + "/Tests/"
print("Models read from dir:  %s" % InModDir)
PlotDir = InModDir
print("Models written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)


mod_files = [
"A1_NM_intersection_FL11126-0_delete_PLM10_k3_prof1_nlyr32_TikhOpt_gcv_Prior30_Results.npz",
"A1_NM_intersection_FL11126-0_delete_PLM10_k3_prof1_nlyr32_TikhOpt_gcv_Prior3000_Results.npz"
              ]

ReadFilelist = False
if ReadFilelist:
    mod_files = []
    with open("Filelist.txt", "r") as file:
        for line in file:
            mod_files.append(line[:-1])

SearchFilelist = False

if SearchFilelist:
    SearchStrng = ""
    print("Search string: %s " % SearchStrng)
    mod_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InModDir)


mod_files = sorted(mod_files)
ns = numpy.size(mod_files)

if ns ==0:
    error("No files found!. Exit.")

Fileout = PlotDir +"AEM05_TikhOpt_P30-3000_DoI.npz"
Header = "DoI from FD models"

start = time.time()


icount = -1
for file in mod_files:
    icount = icount+1

    m = numpy.load(InModDir +file)
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
m_doi = numpy.abs(m_dif)/numpy.abs(r_dif)
# print("doi", numpy.shape(doi))

print("DoI min = "+str(numpy.amin(m_doi)))
print("    man = "+str(numpy.amax(m_doi)))

numpy.savez_compressed(
       file=Fileout,
       mod_files=mod_files,
       header=Header,
       mod_ref=m_ref,
       mod_ref=m_ref,
       mod_doi=m_doi,
    )
"""
    numpy.savez_compressed(
        file=Fileout,
        fl_data=file,
        fl_name=fl_name,
        header=Header,
        ctrl = **Ctrl,
        mod_ref=mod_apr,
        mod_act=mod_act,
        dat_act=dat_act,
        site_modl=site_modl,
        site_sens=site_sens,
        site_merr=site_merr,
        site_dobs=site_dobs,
        site_dcal=site_dcal,
        site_derr=site_derr,
        site_nrms=site_nrms,
        site_nump=site_nump,
        site_num=site_num,
        site_y=site_y,
        site_x=site_x,
        site_gps=site_gps,
        site_alt=site_alt,
        site_dem=site_dem)
"""
print("Results stored to "+Fileout)
elapsed = (time.time() - start)
print (" Used %7.4f sec \n" % (elapsed))

print("\nAll done!")
