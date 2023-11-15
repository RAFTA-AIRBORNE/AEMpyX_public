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
from cycler import cycler
import matplotlib
import matplotlib.pyplot

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
import viz

warnings.simplefilter(action="ignore", category=FutureWarning)


OutInfo = True

AEMPYX_DATA = os.environ["AEMPYX_DATA"]


print(" AEMpyX FDEM forward modeling ")

rng = numpy.random.default_rng()
nan = numpy.nan  # float("NaN")
version, _ = versionstrg()
titstrng = util.print_title(version=version, fname=__file__, out=False)
print(titstrng+"\n\n")

OutInfo = False
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
AEM_system = "aem05"  # "genesis"
if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.03
    data_active = numpy.ones(NN[2], dtype="int8")


if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype="int8")
    data_active[0:11]=0  # only vertical component
    # data_active[10:11]=0  # Vertical + 'good' hoizontals'


nD = NN[0]


print("AEM system: " + AEM_system + "\n \n")

OutDir  = AEMPYX_DATA+"/SYNTH/data/"

if not os.path.isdir(OutDir):
    print("File: %s does not exist, but will be created" % OutDir)
    os.mkdir(OutDir)


ASCout = False
NPZout = True



"""
These are loops over different parameters, in this case for a 3-Layer case.
Should be adapted according to your needs.
"""

depth = [25., 50., 100.]
thick = [10., 30.]
rho0 = [ 100.]
rho1 = [ 10., 1000.]
rhob = [100.]
altitudes  = [60.,120.,180.,240.]


Nsamples = 1000
Perturb = True



"""
Define model
"""
nlyr = 3
Model_active, Model_base, model_var, m_bounds, m_state = inverse.init_1dmod(nlyr)
Model_active[0 * nlyr:1 * nlyr] = numpy.ones((nlyr))


"""
Generate Data

"""


ncount = 0
for di in depth:
    for ti in thick:
        for r1i in rho1:
            for r0i in rho0:
                for rbi in rhob:
                    for alt in altitudes:
                    # IdString = AEM_system.upper()+"_1Layer_Resistor"
                        IdString = AEM_system.upper()\
                                +"_Rho_"+str(int(r0i))+"_"+str(int(r1i))+"_"+str(int(rbi))\
                                +"_Thk_"+str(int(di))+"_"+str(int(ti))+"_Alt_"+str(int(alt))
                        Header = IdString+" | "+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S"))

                        m_i = Model_base.copy()

                        m_i[0 * nlyr:1 * nlyr]       = [r0i, r1i, rbi]
                        m_i[6 * nlyr:7 * nlyr - 1]   = [di, ti]

                        # m_i = Model_base.copy()
                        # m_i[0:nlyr] = r0i
                        # layer = (dl>=di) & (dl<=di+ti)
                        # for il in range(nlyr):
                        #      if layer[il]:
                        #          m_i[il] =r1i

                        d_state = 0
                        m_state = 0

                        m_current, m_state = inverse.transform_parameter(m_vec=m_i, m_trn=ParaTrans, m_state=m_state, mode="f")
                        d_current, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall, alt=alt,
                                                          m_vec = m_current, m_trn=ParaTrans, m_state=m_state,
                                                          d_trn=0, d_state=d_state, d_act=data_active )
                        p_current = [di, ti, r0i, r1i, DataTrans, DatErr_add, DatErr_mult]



                        for ismp in numpy.arange(Nsamples):
                            ncount=ncount+1
                            data_err, data_obs = inverse.set_errors(d_current, DatErr_add, DatErr_mult, perturb=Perturb)

                            p_current = [di, ti, r0i, r1i, rbi, DataTrans, DatErr_add, DatErr_mult]
                            description = "Parameters: "\
                                       +"{0:.0f} ".format(p_current[0])\
                                       +"{0:.0f} ".format(p_current[1])\
                                       +"{0:.0f} ".format(p_current[2])\
                                       +"{0:.1f}/{1:.1f} ".format(p_current[3],p_current[4])\
                                       +"{0:.0f}".format(p_current[5])\
                                       +"{0:.0f}/{1:.2f} ".format(p_current[6],p_current[7])\
                                       +"{0:2d}".format(ncount)
                            print("\n\n\nModel description:")
                            print(description)

                            if "aem" in AEM_system.lower():
                                line = numpy.concatenate(([99999.99], [float(ncount)],[0.], [0.],[alt],[0.], data_obs[:], [0.], [0.],[0.]))
                            else:
                                line = numpy.concatenate(([99999.99], [float(ncount)],[0.], [0.],[alt],[0.], data_obs[:], [0.], [0.]))

                            if ismp == 0:
                                Model = m_i
                                Data  = line
                                Info = [description]
                            else:
                                Model = numpy.vstack((Model, m_i))
                                Data =  numpy.vstack((Data, line))
                                Info.append(description)



                        if NPZout:
                            NPZFile = OutDir+IdString+".npz"
                            print("\n\nResults written to "+NPZFile)
                            aesys.write_aempy(File=NPZFile, Data=Data, System=AEM_system, Header=Header, OutInfo=OutInfo)
