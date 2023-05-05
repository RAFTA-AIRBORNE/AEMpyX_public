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
# AEM_system = "aem05"
AEM_system = "genesis"

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
Alt = [90., 180.]


IdString = AEM_system.upper()+"_5Layer"
nlyr = 5
Model_active, Model_base, model_var, m_bounds, m_state = inverse.init_1dmod(nlyr)
Model_active[0 * nlyr:1 * nlyr] = numpy.ones((nlyr)).astype(int)

# IdString = AEM_system.upper()+"_5Layer"
# Thick0 = [25.]
# Thick1 = [10., 20.]
# Thick2 = [10., 20.]
# Thick3 = [10., 20.]
# Rho0 = [100.]
# Rho1 = [10., 1000]
# Rho2 = [100.]
# Rho3 = [10., 1000]
# Rhob = [100.]

IdString = AEM_system.upper()+"_3Layer"
Thick0 = [25., 50., 100.]
Thick1 = [10., 30.]
Thick2 = [1.]
Thick3 = [1.]
Rho0 = [ 100.]
Rho1 = [ 10., 1000]
Rho2 = [100.]
Rho3 = [100.]
Rhob = [100.]

"""
Generate Data

"""
Nsamples = 1000
Perturb = True
SplitData= True

mod_num = -1
for alt in Alt:
    for thk0 in Thick0:
        for thk1 in Thick1:
            for thk2 in Thick2:
                for thk3 in Thick3:
                   for rho0 in Rho0:
                        for rho1 in Rho1:
                            for rho2 in Rho2:
                                for rho3 in Rho3:
                                    for rhob in Rhob:

                                        mod_num += 1

                                        p_i = [mod_num, alt,
                                               thk0, thk1, thk2, thk3,
                                               rho0, rho1,rho2, rho3, rhob,
                                               DataTrans, DatErr_add, DatErr_mult]

                                        description =    "{0:2d} ".format(p_i[0])\
                                                        +"{0:.0f} ".format(p_i[1])\
                                                        +"{0:.0f} ".format(p_i[2])\
                                                        +"{0:.0f} ".format(p_i[3])\
                                                        +"{0:.0f} ".format(p_i[4])\
                                                        +"{0:.0f} ".format(p_i[5])\
                                                        +"{0:.0f} ".format(p_i[6])\
                                                        +"{0:.0f} ".format(p_i[7])\
                                                        +"{0:.0f} ".format(p_i[8])\
                                                        +"{0:.0f} ".format(p_i[9])\
                                                        +"{0:.0f} ".format(p_i[10])\
                                                        +"{0:2d} ".format(p_i[11])\
                                                        +"{0:.0f}/{1:.2f} ".format(p_i[12],p_i[13])

                                        print("\n\n\nModel description:")
                                        print(description)


                                        m_i = Model_base.copy()

                                        m_i[0 * nlyr:1 * nlyr] = [rho0, rho1, rho2, rho3, rhob]
                                        m_i[6 * nlyr:7 * nlyr - 1] = [thk0, thk1, thk2, thk1]

                                        d_state = 0
                                        m_state = 0

                                        m_current, m_state = inverse.transform_parameter(m_vec=m_i, m_trn=ParaTrans, m_state=m_state, mode="f")
                                        d_ref, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall, alt=alt,
                                                                          m_vec = m_current, m_trn=ParaTrans, m_state=m_state,
                                                                          d_trn=0, d_state=d_state, d_act = DataActive )

                                        if mod_num==0:
                                            Model = m_i
                                            Info = [description]
                                            Data = numpy.insert(d_ref,0,[mod_num, -1, alt])
                                            Para = p_i
                                        else:
                                            Model = numpy.vstack((Model, m_i))
                                            Info.append(description)
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
        i_s = Info[imod]
        d_s = Data[numpy.isin(Data[:,0],imod)]


        SplitStrng = "Mod"+"_"\
                    +str(int(p_s[0]))+"_"\
                    +"Alt"+"_"\
                    +str(int(p_s[1]))+"_"\
                    +"Thick"+"_"\
                    +str(int(p_s[2]))+"_"\
                    +str(int(p_s[3]))+"_"\
                    +str(int(p_s[4]))+"_"\
                    +str(int(p_s[5]))+"_"\
                    +"Res"+"_"\
                    +str(int(p_s[6]))+"_"\
                    +str(int(p_s[7]))+"_"\
                    +str(int(p_s[8]))+"_"\
                    +str(int(p_s[9]))+"_"\
                    +str(int(p_s[10]))

        NPZSplit=OutDir+IdString+SplitStrng+".npz"
        print("Results written to "+NPZSplit)
        numpy.savez_compressed(file=NPZSplit, model=m_s, data=d_s, para=p_s, info=i_s)

print(numpy.shape(Data))
NPZFile = OutDir+IdString+".npz"
print("\n\nResults written to "+NPZFile)
numpy.savez_compressed(
    file=NPZFile, model=Model, data=Data, para=Para, info=Info)
