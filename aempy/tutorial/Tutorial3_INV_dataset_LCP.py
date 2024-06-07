#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,ipynb
#     main_language: python
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (Spyder)
#     language: python3
#     name: python3
# ---
# -*- coding: utf-8 -*-


"""
This script realizes the Laterally Correlation Procedure approach
of Christensen (2009). 

References:
           
    N. B. Christensen & R. J. Tølbøll (2009)
    “A lateral model parameter correlation procedure for one-dimensional inverse modelling""
    Geophysical Prospecting, 57, 919–929, doi: 10.1111/j.1365-2478.2008.00756.x.
  
 
Created on Tue Aug  3 17:03:39 2021

vrath 10/23

"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime


import numpy
import scipy
import scipy.interpolate
import scipy.spatial
import scipy.linalg

import shapely

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys
import viz
import inverse
import post

#warnings.simplefilter(action="ignore", category=FutureWarning)
cm = 1/2.54
OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
script = "Tutorial3_INV_dataset_LCP.py"
# script = __file__  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")
Header = titstrng


now = datetime.now()

"""
System related settings.
Data transformation is now allowed with three possible options:
DataTrans   = 0           raw data
            = 1           natural log of data
            = 2           asinh transformation
An error model is applied for the raw data, which is
mixed additive/multiplicative. in case of data transformation,
errors are also transformed.
"""
# AEM_system = "genesis"
AEM_system = "aem05"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.03
    data_active = numpy.ones(NN[2], dtype="int8")
    CompDict=Misc[2]


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 2
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical +             print(numpy.shape(rect))"good" hoizontals"
    CompDict =Misc[2]

"""
input formats is "npz"
"""
InFilFmt = ".npz"
XYFact = 1. 


TileSize = 2000.
TileOverlap = 0.5
TileMinSites = 3


LayerWise = True
CovarThresh = 500.
Scale = 0.5*CovarThresh

ReCalc = "fwd"   # "inverse"

MergeModels = True

MergeFile = "Limerick_shale_dec5_merged.npz"
SearchStrng = "*delete_dec5*mean*results.npz"

#MergeFile = "Limerick_shale_k2_dec5_merged.npz"
#SearchStrng = "*k2_dec5*mean*results.npz"


AEMPYX_DATA =  AEMPYX_ROOT + "/data/"

InModDir = AEMPYX_DATA+"/aem05_limerick/dec/results/"
print("Data read from dir: %s " % InModDir)
FileList = "search" #"search"

OutModDir =  AEMPYX_DATA+"/aem05_limerick/merged/"


if "set" in FileList.lower():
    mod_files = []

if "read" in FileList.lower():
    ListName=""
    print("File names read from : "+ListName)
    how = ["read", ListName, InModDir]
    mod_files = util.get_data_list(how=how,
                              out= True, sort=True)

    mod_files = numpy.loadtxt("A9-7.dat", dtype=str)

if "search" in FileList.lower():
    print("Searchstring is : "+SearchStrng)
    how = ["search", SearchStrng, InModDir]
    mod_files = util.get_data_list(how=how,fullpath=True,
                              out= True, sort=True)

ns = numpy.size(mod_files)
if ns ==0:
    error("No files set!. Exit.")




print(mod_files[0])
print("Data read from dir: %s " % InModDir)

#   workaround!!!!!    
corrfile = MergeFile

if MergeModels:
    _ = util.merge_model_sets(infile_list=mod_files,
                                   outfile_name=MergeFile,
                                   dictout=True, out=False)
    mod_files = [corrfile]

"""
read  data set
"""

for filein in mod_files:

    print("\nMerged models read from: %s" % filein)

    models = numpy.load(filein, allow_pickle=True)
   

    e = models["x"]*XYFact
    e_min = numpy.amin(e)
    e_max = numpy.amax(e)
    n = models["y"]*XYFact
    n_min = numpy.amin(n)
    n_max = numpy.amax(n)
    
    d = models["d"]
    m = models["mod"]
    c = models["cov"]
    r = models["rms"]

    # being developed: good_index = post.mod_qc(model=m, model_error=numpy.diag(), data_fit=r)
    
    
    dims= numpy.shape(d)
    m = numpy.reshape(m, (dims[0], dims[1]))
    c = numpy.reshape(c, (dims[0], dims[1]*dims[1]))
    

    if ParaTrans==1:
       m = numpy.log(m)

    mod_cor = m.copy()
    

    """
    Step 1: calculate the laterally correlated moidle set
    setup overlapping tiles
    
    """
    dxtiles = TileSize
    xtiles = numpy.arange(numpy.floor(e_min), numpy.ceil(e_max), TileOverlap*dxtiles) 
    nx = len(xtiles)
    
    dytiles = TileSize
    ytiles = numpy.arange(numpy.floor(n_min), numpy.ceil(n_max), TileOverlap*dytiles) 
    ny = len(ytiles)
    
    
    start = process_time()
    total = start
    
    
    itile = -1
    ntile = nx*ny
       
    
    for ii in numpy.arange(nx):
        for jj in numpy.arange(ny):
            itile = itile+1
            ll=  [xtiles[ii],ytiles[jj]]
            ur = [xtiles[ii]+dxtiles,ytiles[jj]+dytiles]
            print("\n\n Tile",itile,"of", ntile)
            print("Rect lower left  (m): "+str(ll[0])+", "+str(ll[1]))
            print("Rect upper right (m_): "+str(ur[0])+", "+str(ur[1]))
            
            # rect = numpy.array([])
            inside = numpy.where((e>ll[0]) & (e<ur[0]) & (n>ll[1]) & (n<ur[1]))
         
            e_tile = e[inside]
            n_tile = n[inside]
            d_tile = d[inside[0],:]
            m_tile = m[inside[0],:]
            c_tile = c[inside[0],:]
            
            nsit, nlyr = numpy.shape(m_tile)
            print("Tile",itile,"contains", nsit, "sites with", nlyr, "layers.")
            
            if nsit > TileMinSites:
            
                c_tile = numpy.reshape(c_tile, (nsit, nlyr, nlyr))
                
                if LayerWise:
                    points = numpy.stack(      
                        [ e_tile.ravel(order="F"),   
                          n_tile.ravel(order="F")
                        ], -1)
                else:
                    points = numpy.stack(      
                        [ e_tile.ravel(order="F"),   
                          n_tile.ravel(order="F"),
                          d_tile.ravel(order="F")
                        ], -1)
                    
                dists  = scipy.spatial.distance.squareform(
                    scipy.spatial.distance.pdist(points, metric="euclidean"))
                cov_s = numpy.linalg.inv(numpy.exp(-dists/Scale))
            
                cov_i = c_tile.copy()       
                          # for isit  in numpy.arange(nsit):               
                #    cov_i[isit,:, :] = scipy.linalg.inv(c_tile[isit,:,:])
                if LayerWise:
                    for ilyr in numpy.arange(nlyr):
                        par_e = m_tile[:, ilyr]
                        
                        cov_e = numpy.diag(1./cov_i[:, ilyr, ilyr])
                        cov_c = numpy.linalg.inv(cov_e + cov_s)
                        par_c = cov_c@cov_e@par_e
                        
                        m_tile[:, ilyr] = par_c

                                                
                else:
                    par_e = m_tile
                    cov_e =  numpy.diag(1./numpy.diagonal(cov_i, axis1=1, axis2=2))
                    par_c = numpy.linalg.solve(cov_e + cov_s, cov_e@par_e )
                    m_tile = par_c
                    
               
            else:
                
                print("Not enough sites in tile.")
                
                

            mod_cor[inside[0],:] = m_tile
          
            
            print("Tile",itile,", norm of differences:", 
                  numpy.linalg.norm(mod_cor[inside[0]]-m[inside[0]])/(nsit*nlyr)) 
                
            
    elapsed = process_time()
    print("\n\n")
    print("Time used for correlating models:", elapsed-start,"s")
    print("Time per Tile:", (elapsed-start)/ntile)
    
    
    dateform="%m/%d/%Y, %H:%M:%S"
    header = "corr model set:"+"".join("Date " + datetime.now().strftime(dateform))
    
    models_dict = dict(models)
    models_dict["header"] = header 
    models_dict["mod_cor"] = numpy.exp(mod_cor)

    
    numpy.savez_compressed(corrfile, **models_dict)

    print(list(models_dict.keys()))
                 
    # """
    # Step 3:
    # Run either forward models to check data fit or re-run inversion with 
    # correlated models as prior. 
        
    # """
    # start  = process_time()    
    
    
    
    # print("\n\n")
    # print("Time used for recalulating models:", elapsed-start,"s")

    # elapsed = process_time()
    
    elapsed = process_time()
    print("\n\n")
    print("Total time:", elapsed-total,"s")

    
