#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadat_filter: -all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
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
# from time import process_time
# from random import randrange
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

# import aesys
import util
import inverse

warnings.simplefilter(action="ignore", category=FutureWarning)

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")

version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = True
now = datetime.now()

Method = "Oldenburg1999"
Method = "Variance"

# InDatDir =  AEMPYX_DATA + "/Projects/InvParTest/proc_delete_PLM3s/"
InModDir =  AEMPYX_ROOT + "/aempy/data/AEM05/results_doi/"
if not InModDir.endswith("/"): InModDir=InModDir+"/"


# Method = "Oldenburg1999"
# FileList = "set" 
# SearchStrng = ""

Method = "Variance"
FileList = "search"  
SearchStrng = "*results.npz"


# InDatDir =  AEMPYX_DATA + "/Projects/InvParTest/proc_delete_PLM3s/"
InDatDir =  AEMPYX_ROOT + "/aempy/data/AEM05/"
if not InDatDir.endswith("/"): InDatDir=InDatDir+"/"


if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    # dat_files = []
    mod_files = ["StGormans_FL11379-0_k3_nlyr36_TikhOpt_gcv_Prior10_Err_a75-m5_results.npz",
                 "StGormans_FL11379-0_k3_nlyr36_TikhOpt_gcv_Prior1000_Err_a75-m5_results.npz"]  
    # numpy.load(AEMPYX_DATA + "/Projects/Compare/BundoranSubsets.npz")["setC"]
    
else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    mod_files = util.get_data_list(how=["search", SearchStrng, InModDir],
                              out= True, sort=True)

    mod_files = [os.path.basename(f) for f in mod_files]  

ns = numpy.size(mod_files)
if ns ==0: error("No files set!. Exit.")

"""
Output format is ".npz"
"""
OutFileFmt = ".npz"
OutModDir =  InModDir
print("Results written to dir: %s " % OutModDir)

if not os.path.isdir(OutModDir):
    print("File: %s does not exist, but will be created" % OutModDir)
    os.mkdir(OutModDir)


Fileout = OutModDir +"AEM05_TikhOpt_P30-3000_DoI_"+Method+".npz"
Header = "DoI from FD models"

start = time.time()

results = numpy.load(InModDir + mod_files[0])

numpy.savez_compressed(Fileout,
        fline = results["fl_name"],
        m_act = results["mod_act"],
        m_ref = results["mod_ref"],
        site_x = results["site_x"],
        site_y = results["site_y"],
        site_z = results["site_dem"],
        site_dact = results["dat_act"],
        site_dobs = results["site_dobs"],
        site_derr = results["site_derr"],
        doimethod = Method.lower(),
        )


if "old" in Method.lower():
    icount = -1
    for file in mod_files:
        icount = icount+1

        m = numpy.load(OutModDir+file)
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
    v_doi = numpy.abs(m_dif)/numpy.abs(r_dif)
    # print("doi", numpy.shape(doi))

    print("\n\nMethod:  "+Method)
    print("DoI min = "+str(numpy.amin(v_doi)))
    print("    max = "+str(numpy.amax(v_doi)))
    
    util.add_object_npz(filein=Fileout, 
               xkey = ["doi_files"], xobject=[mod_files])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_mavg"], xobject=[m_avg])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_mdif"], xobject=[m_dif])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_rdif"], xobject=[r_dif])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_doi"], xobject=[v_doi])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_models"], xobject=[modls])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_mrefs"], xobject=[mrefs])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_merrs"], xobject=[merrs])    
    util.add_object_npz(filein=Fileout, 
               xkey = ["doi_meth"], xobject=[Method])

    

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
                error("model dimensions do not fit! Exit.")
            # print( numpy.shape(modls))

    modls = numpy.log10(modls.reshape(ns, nsit, nlyr))
    mrefs = numpy.log10(mrefs.reshape(ns, nlyr))

    m_avg = numpy.mean(modls, axis=0).reshape(nsit,nlyr)
    m_ano = numpy.zeros_like(modls)
    for nn in numpy.arange(ns):
        m_ano[nn,:,:] = modls[nn,:,:] - m_avg

    m_var = numpy.var(m_ano,axis = 0)
    v_doi = numpy.sqrt(m_var)
    
    print("\n\nMethod:  "+Method)
    print("DoI min = "+str(numpy.amin(v_doi)))
    print("    max = "+str(numpy.amax(v_doi)))
    
    util.add_object_npz(filein=Fileout, 
               xkey = ["doi_files"], xobject=[mod_files])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_mavg"], xobject=[m_avg])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_doi"], xobject=[v_doi])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_ano"], xobject=[m_ano])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_var"], xobject=[m_var])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_models"], xobject=[modls])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_mrefs"], xobject=[mrefs])
    util.add_object_npz(filein=Fileout, 
               xkey = ["site_doi_merrs"], xobject=[merrs])
    util.add_object_npz(filein=Fileout, 
               xkey = ["doi_meth"], xobject=[Method])


print("Results stored to "+Fileout)
elapsed = (time.time() - start)
print (" Used %7.4f sec \n" % (elapsed))

print("\nAll done!")
