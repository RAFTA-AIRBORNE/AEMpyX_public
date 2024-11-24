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
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
# #!/usr/bin/env python3
# -

# This script allows you to do forward modelling, with several options on the output. The purpose for including this is multifold: (1) It is useful to see the response for a given model which may be hypothetical, to see what might be inverted for. (2) A series of models for parameter studies is possible. (3) a set of (perturbed) responses can be generated, which in turn may be fed into one of the inversion algorithms.

# +
import time
import sys
from sys import exit as error
import os
import warnings
from time import process_time
from datetime import datetime

import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]


for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import core1d
import inverse
import aesys
# -


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

# +
rng = numpy.random.default_rng()
nan = numpy.nan

version, _ = versionstrg()
#script = "Tutorial0_FWD_synth.py"
script = __file__  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")
Header = titstrng
# -

OutInfo = False
AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/Synth/"
OutDir  = AEMPYX_DATA+"/data/"
if not os.path.isdir(OutDir):
    print("File: %s does not exist, but will be created" % OutDir)
    os.mkdir(OutDir)


# The following cell gives values to AEM-system related settings.
#
# Data transformation is activated by the variable DataTrans. Currently three possible options are allowed: _DataTrans = 0_: No transformation, i.e., the raw data are used. _DataTrans = 1_: The natural log of data is taken, only allowed for strictly positive values. _DataTrans = 2_: If data scale logarithmically, an asinh transformation (introduced by Scholl, 2000) is applied. It allows negatives, which may occur in TDEM, when IP effects are present.
#
# A general additive/multiplicative error model is applied on the raw data before transformation, and errors are also transformed.

# +
AEM_system = "aem05"
# AEM_system = "genesis"
print("AEM system: " + AEM_system + "\n \n")

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 50.
    DatErr_mult = 0.00
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
# -

# In case an ensemble of model responses is desired, e.g. for future inversions, the resukting output can be controlled here.

Nsamples = 1000
# NSamples = 1
PerturbDat = True

SplitData= True

# To initialize loops over different parameters,
# first a reference model must be set up, with reasonable values for all parameters not within the loop. Default settings is rho only, no IP. Currently, one parameter and altitude can be varied within a loop.
# The following should be adapted according to the user's needs.

nlyr = 3
Model_active, Model_base, model_var, m_bounds, m_state = inverse.init_1dmod(nlyr)

Model_base[0*nlyr:1*nlyr] =[100., 100., 100.]   #rho
Model_base[6*nlyr:7*nlyr-1] =[30.,30.]          #layers

# +
# Adapted for reasonable IP values

# Model_base[3*nlyr:4*nlyr] =[0.,  0.5, 0.]      #chargeability
# Model_base[4*nlyr:5*nlyr] =[0.,  0.5, 0.]      #exponent
# Model_base[5*nlyr:6*nlyr] =[0., 100., 0.]      #frequency
# -

# rho for layer 1 (starting from 0!)
FWDBaseName = "AEM05_Rho1"
VarPar = [ 10., 100.,1000.]
VarInd = 0 * nlyr+1

# +

# thickness of layer 1 (starting from 0!)
# FWDBaseName = "AEM05_Thk1"
# VarPar = [10., 30., 50.]
# VarInd = 6*nlyr+1

# chargeability of layer 1 (starting from 0!)
# FWDBaseName = "AEM05_m1"
# VarPar = [0.0001, 0.2, 0.4, 0.6, 0.8]
# VarInd = 3*nlyr+1

#Alt = [60., 120.]
Alt = [60]
# -

# Now generate the response data:


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
            _, data_obs = inverse.set_errors(d_ref, DatErr_add, DatErr_mult, perturb=PerturbDat)
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
