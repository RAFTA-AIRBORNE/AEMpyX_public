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
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (Spyder)
#     language: python3
#     name: python3
# ---

# +
#!/usr/bin/env python3
# -

'''
PROJECT_interpolate_data_area.py - AEMpyX spatial data interpolation.

Provenance
----------
AEMpyX project.

@authors: Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)
'''
# This script plots data over an spatial area.

# +
import os
import sys

from time import process_time
import inspect

import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker

import scipy.interpolate
import scipy.spatial

AEMPYX_ROOT = os.environ['AEMPYX_ROOT']
mypath = [os.path.join(AEMPYX_ROOT, 'aempy/modules/')]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys

# +
OutInfo = True
cm = 1/2.54
AEMPYX_DATA = os.environ['AEMPYX_DATA']

version, _ = versionstrg()
# script = 'Tutorial1_VIZ_data_area.py'
script = inspect.getfile(inspect.currentframe())  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+'\n\n')
Header = titstrng
# -

# The following cell gives values to AEM-system related settings.
#
# Data transformation is activated by the variable _DataTrans_. Currently
# three possible options are allowed: _DataTrans = 0_: No transformation,
# i.e., the raw data are used. _DataTrans = 1_: The natural log of data
# is taken, only allowed for strictly positive values. _DataTrans = 2_:
# If data scale logarithmically, an _asinh_ transformation (introduced by
# Scholl, 2000) is applied. It allows negatives, which may occur in TDEM,
# when IP effects are present.
#
# A general additive/multiplicative error model is applied on the raw data
# before transformation, and errors are also transformed.

# +
# AEM_system = 'genesis'
AEM_system = 'aem05'
if 'aem05' in AEM_system.lower():
    _, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = [0]
    DatErr_add =  50.
    DatErr_mult = 0.05
    data_active = numpy.ones(NN[2], dtype='int8')
    CompDict = Misc[3]
    CompLabl = list(CompDict.keys())
    print(CompLabl)

if 'genes' in AEM_system.lower():
    _, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = [2, numpy.nan]
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype='int8')
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'
    CompDict =Misc[2]
    CompLabl = list(CompDict.keys())

# +

InFileFmt = '.npz'
XYFact = 1.e-3

##############################################################################
# A2 Somaye
##############################################################################
# raw data
InDatDir =  AEMPYX_DATA + '/AllDataBlocks/Archive/A2/raw'
DataFile= 'A2_full_data*.npz'
DataName = 'Profile_13_AEM'
DataSet = 'prof'

if 'mesh'  in DataSet:
    numIndexes = [121, 181]
    
elif 'prof'  in DataSet:
    
    # prof 13
    ProfFileName = AEMPYX_DATA + '/AllDataBlocks/Archive/A2/'+'19-BOL-13_CDP_Coords_1_utm29N.csv'
    NumIndexes = [2232]
    # prof 14
    # ProfFileName = AEMPYX_DATA + '/AllDataBlocks/Archive/A2/'+'19-BOL-14_CDP_Coords_1_utm29N.csv'
    # NumIndexes = [2016]
    
    with open(ProfFileName, 'r') as f:
        prof_data = f.readlines()
      
    tmp = prof_data[0].split(',')  
    prof_start = (float(tmp[0]), float(tmp[1]))
    tmp = prof_data[-1].split(',') 
    prof_end = (float(tmp[0]), float(tmp[1]))
    print ('profile 13 start', prof_start, '(utm 29 N(U)')
    print ('profile 13 end', prof_end, '(utm 29 N(U)')
    # print ('profile 13 start', prof_start, '(utm 29 N(U)')
    # print ('profile 13 end', prof_end, '(utm 29 N(U)')
else: 
    sys.exit('Data type',DataSet, 'not defined! Exit.')
    
    
OutDatDir = InDatDir

print('Data read from dir: %s ' % InDatDir)
print('Plot filname: %s ' % DataName)


dat_files = [DataFile]

ns = numpy.size(dat_files)
if ns ==0:
    sys.exit('No files set!. Exit.')


step = 1

# The following cell determines the settings for individual components. Each sublist associated to a componet contains the name, followed by a list of parameters determining the data limits, and a step determining the color bar, or the isolines. Further paramers, as e.g. the threshhold for the PLM, may be added.

CompList=['P1', 'Q1', 'P2', 'Q2', 'P3', 'Q3', 'P4', 'Q4', 'PLM', 'ALT']  

if not os.path.isdir(InDatDir):
    print('File: %s does not exist, but will be created' % OutDatDir)
    os.mkdir(OutDatDir)


for filein in dat_files:
    start = process_time()
    print('\nData read from: %s' % filein)
    Data, header, _ = aesys.read_aempy(File=filein, System=AEM_system, OutInfo=False)


    E = Data[:,1][::step]*XYFact
    E_min = numpy.amin(E)
    E_max = numpy.amax(E)
    N = Data[:,2][::step]*XYFact
    N_min = numpy.amin(N)
    N_max = numpy.amax(N)
    Z = Data[:,5][::step] 

    if 'mesh' in DataSet:
        xi= numpy.linspace(E_min,E_max,numIndexes[0])
        yi= numpy.linspace(N_min,N_max,numIndexes[1])
        dx = numpy.around(numpy.diff(xi)[0]/XYFact, decimals=0)
        dy = numpy.around(numpy.diff(yi)[0]/XYFact, decimals=0)
        print('Interpolation mesh, dx = '+ str(dx)+' m, dy ='+ str(dy)+' m')
        
        XI, YI = numpy.meshgrid(xi, yi, indexing='ij')
        
        
        Pnts = numpy.stack([ E.ravel(),  N.ravel()], -1)
        Mesh = numpy.stack([XI.ravel(), YI.ravel()], -1)
    

        D_tree=scipy.spatial.KDTree(Pnts, leafsize=10,
                                    compact_nodes=True,
                                    copy_data=True,
                                    balanced_tree=True,
                                    boxsize=None,
                                    workers=-1)
        mindist, index = D_tree.query(Mesh, k=1)
        DI = Dats[index].reshape(XI.shape)


    
    elif 'prof'  in DataSet:
        xi= numpy.linspace(E_min,E_max,NumIndexes[0])
        yi= numpy.linspace(N_min,N_max,NumIndexes[0])

    for nc in numpy.arange(len(CompList)-2):

        Comp = CompList[nc]
        comp = CompDict[Comp][0]

        D = Data[:,comp][::step]

        D_min = numpy.amin(D)
        D_max = numpy.amax(D)
        print('Data, read   min='+str( D_min)+'   max='+str( D_max))

        Dats = D.flatten()


        D_min = numpy.nanmin(DI)
        D_max = numpy.nanmax(DI)
        print('Data, interpolated   min='+str( D_min)+'   max='+str( D_max))


