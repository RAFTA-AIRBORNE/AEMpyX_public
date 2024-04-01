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
#       format_version: '1.5'
#       jupytext_version: 1.15.2
# ---

"""
Created on Tue Sep  6 10:57:01 2016

@author: vrath

edited by dkiyan - Sep 30
edited by vrath  - May 7, 2021

"""
import time
import sys
from sys import exit as error
import os
import warnings
from time import process_time
from datetime import datetime

import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import util
#import core1d_par as core1d
import core1d
import inverse
import aesys


warnings.simplefilter(action="ignore", category=FutureWarning)


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = False
now = datetime.now()


OutDir  = AEMPYX_DATA+"/SynthData/data/"

if not os.path.isdir(OutDir):
    print("File: %s does not exist, but will be created" % OutDir)
    os.mkdir(OutDir)


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
AEM_system = "aem05"
# AEM_system = "genesis"

print("AEM system: " + AEM_system + "\n \n")

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 50.
    DatErr_mult = 0.03
    alt = 60.
    DataActive = numpy.ones((1,NN[2]))

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 100.
    DatErr_mult = 0.01
    alt = 90.
    DataActive = numpy.ones((1,NN[2]))

nD = NN[0]

"""
Define models:
These are loops over different parameters, in this case for a 3-Layer case.
Should be adapted according to your needs.
"""

#Alt = [60., 120.]
Alt = [60]

Nsamples = 1000
# NSamples = 1
Perturb = True
SplitData= True

"""
Set up base model
"""

nlyr = 3
Model_active, Model_base, model_var, m_bounds, m_state = inverse.init_1dmod(nlyr)

"""
Background model: default settings is rho only, - IP is nonexistent 
Adapted for reasonable IP values
""" 
Model_base[0*nlyr:1*nlyr] =[100., 100., 100.]   #rho
Model_base[6*nlyr:7*nlyr-1] =[30.,30.]          #layers 

Model_base[3*nlyr:4*nlyr] =[0.,  0.5, 0.]      #chargeability
Model_base[4*nlyr:5*nlyr] =[0.,  0.5, 0.]      #exponent
Model_base[5*nlyr:6*nlyr] =[0., 100., 0.]      #frequency



"""
Currently, one parameter  and altitude can be varied within a loop. 
"""

"""
rho for layer 1 (starting from 0!)
"""

# FWDBaseName = "AEM05_Rho1"
# VarPar = [ 10., 100.,1000.]
# VarInd = 0 * nlyr+1

"""
thickness of layer 1 (starting from 0!)
"""
# FWDBaseName = "AEM05_Thk1"
# VarPar = [10., 30., 50.] 
# VarInd = 6*nlyr+1
"""
chargeability of layer 1 (starting from 0!)
"""
FWDBaseName = "AEM05_m1"
VarPar = [0.0001, 0.2, 0.4, 0.6, 0.8] 
VarInd = 3*nlyr+1 

"""
Generate Data

"""


mod_num = -1
for par in numpy.arange(len(VarPar)):

        mod_num += 1
        
        m_i = Model_base.copy()

        if VarInd==numpy.size(m_i):
            alt = Alt[par]
            p_i = numpy.array([mod_num, VarInd, Alt[par], DataTrans, DatErr_add, DatErr_mult])


        else:
            m_i[VarInd] = VarPar[par]
            p_i = numpy.array([mod_num, VarInd, VarPar[par], DataTrans, DatErr_add, DatErr_mult])


       

        d_state = 0
        m_state = 0

        m_current, m_state = inverse.transform_parameter(m_vec=m_i, m_trn=ParaTrans, m_state=m_state, mode="f")
        d_ref, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall, alt=Alt,
                                          m_vec = m_current, m_trn=ParaTrans, m_state=m_state,
                                          d_trn=0, d_state=d_state, d_act = DataActive )

        if mod_num==0:
            Model = m_i
            Data = numpy.insert(d_ref,0,[mod_num, -1, alt])
            Para = p_i
            print(numpy.shape(Para))
        else:
            Model = numpy.vstack((Model, m_i))
            Data =  numpy.vstack((Data, numpy.insert(d_ref,0,[mod_num, -1, alt])))
            Para =  numpy.vstack((Para, p_i))
        # print(mod_num, numpy.shape(Model))
       
        for ismp in numpy.arange(Nsamples):
            _, data_obs = inverse.set_errors(d_ref, DatErr_add, DatErr_mult, perturb=Perturb)
            data_obs =numpy.insert(data_obs,0,[mod_num, ismp, alt])
            Data =  numpy.vstack((Data, data_obs))

if SplitData:
    for imod in numpy.arange(mod_num+1):

        p_s = Para[imod]
        m_s = Model[imod]
        d_s = Data[numpy.isin(Data[:,0],imod)]

        SplitStrng = "_model"+str(imod)+"_"+str(Nsamples)+"samples"


        NPZSplit=OutDir+FWDBaseName+SplitStrng+".npz"
        print("Results written to "+NPZSplit)
        numpy.savez_compressed(file=NPZSplit, model=m_s, data=d_s, para=p_s)
else:
    print(numpy.shape(Data))
    NPZFile = OutDir+FWDBaseName+".npz"
    print("\n\nResults written to "+NPZFile)
    numpy.savez_compressed(
        file=NPZFile, model=Model, data=Data, para=Para)
