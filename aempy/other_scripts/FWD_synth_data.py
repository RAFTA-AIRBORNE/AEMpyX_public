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
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

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

nan = numpy.nan  # float("NaN")
rng = numpy.random.default_rng()


version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("FWD"+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

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

FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
nD = NN[0]

if "aem05" in AEM_system.lower():
    FwdCall,NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 30.
    DatErr_mult = 0.02
    alt = 60.
    DataActive = numpy.ones((1,NN[2]))

if "genes" in AEM_system.lower():
    FwdCall, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    ParaTrans = 1
    DataTrans=0
    DatErr_add = 50.
    DatErr_mult = 0.01
    alt = 120.
    DataActive = numpy.ones((1,NN[2]))

nD = NN[0]


print("AEM system: " + AEM_system + "\n \n")

# IdString = AEM_system.upper()+"_1Layer_Resistor"
IdString = AEM_system.upper()+"_1Layer_Conductor"

OutDir  = AEMPYX_DATA
OutName =  OutDir+"/SYNTH_"

if not os.path.isdir(OutDir):
    print("File: %s does not exist, but will be created" % OutDir)
    os.mkdir(OutDir)


ASCout = False
NPZout = True


Plots= False
PlotFormat = [".pdf",]
PdfCatalog = True
if ".pdf" in PlotFormat:
    PdfCName = OutName+IdString+"_Catalog.pdf"
else:
    error(" No pdfs generated. No catalog possible!")
    PdfCatalog = False
if PdfCatalog:
    pdf_list = []

"""
For just plotting to files, choose the cairo backend (eps, pdf, png, jpg...).
If you need to see the plots directly (plot window, or jupyter), simply
comment out the following line. In this case matplotlib may run into
memory problems after a few hundreds of high-resolution plot..
"""
FilesOnly = False

Header = AEM_system+" Synthetics | "+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S"))

"""
These are loops over different parameters, in this case for a 3-Layer case.
Should be adapted according to your needs.
"""
# depth = [30., 60.]
# thick = [10., 20., 40.]
# rho0 = [ 100.]
# rho1 = [ 10., 1000.]
# rhob = [100.]

depth = [25.]
thick = [25.]
rho0 = [ 100.]
rho1 = [ 10.]
rhob = [100.]
altitudes = [90.]

Nsamples = 1000
Perturb = True

"""
Define model
"""
nlyr = 3
Model_active, Model_base, model_var, m_bounds, m_state = inverse.init_1dmod(nlyr)
onlyr = numpy.ones((nlyr)).astype(int)
Model_active[0 * nlyr:1 * nlyr] = numpy.ones((nlyr)).astype(int)


if Plots:
    # Determine graphical parameter.
    # print(plt.style.available)
    matplotlib.pyplot.style.use("seaborn-paper")
    matplotlib.rcParams["figure.dpi"] = 400
    matplotlib.rcParams["axes.linewidth"] = 0.5
    matplotlib.rcParams["savefig.facecolor"] = "none"
    Fontsize = 10
    Labelsize = Fontsize
    Titlesize = Fontsize+2
    Linewidth= 1
    Markersize = 4

    Colors = ["r", "g", "b", "m", "c", "y", "k","r", "g", "b", "m"]
    Lines  = ["-", "--", ":", "-.","--", ":", "-.","--", ":","-","-."]

if FilesOnly:
    matplotlib.use("cairo")

"""
Generate Data

"""


ncount = 0
for dtest in depth:
    for ttest in thick:
        for r1test in rho1:
            for r0test in rho0:
                for rbtest in rhob:

                    m_test = Model_base.copy()

                    m_test[0 * nlyr:1 * nlyr]       = [r0test, r1test, rbtest]
                    m_test[6 * nlyr:7 * nlyr - 1]   = [dtest, ttest]

                    # m_test = Model_base.copy()
                    # m_test[0:nlyr] = r0test
                    # layer = (dl>=dtest) & (dl<=dtest+ttest)
                    # for il in range(nlyr):
                    #      if layer[il]:
                    #          m_test[il] =r1test

                    d_state = 0
                    m_state = 0

                    m_current, m_state = inverse.transform_parameter(m_vec=m_test, m_trn=ParaTrans, m_state=m_state, mode="f")
                    d_current, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall, alt=alt,
                                                      m_vec = m_current, m_trn=ParaTrans, m_state=m_state,
                                                      d_trn=0, d_state=d_state, d_act = DataActive )
                    p_current = [dtest, ttest, r0test, r1test, DataTrans, DatErr_add, DatErr_mult]



                    for ismp in numpy.arange(Nsamples):
                        ncount=ncount+1
                        data_err, data_obs = inverse.set_errors(d_current, DatErr_add, DatErr_mult, perturb=Perturb)

                        p_current = [dtest, ttest, r0test, r1test, rbtest, DataTrans, DatErr_add, DatErr_mult]
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

                        if ncount == 1:
                            Model = m_test
                            Data  = line
                            Info = [description]
                        else:
                            Model = numpy.vstack((Model, m_test))
                            Data =  numpy.vstack((Data, line))
                            Info.append(description)
                        if Plots:
                            PlotFile = OutName+IdString+"_"+str(ncount)

                            if "gen" in AEM_system.lower():
                                viz.plot_site_genesis(
                                        PlotFile = PlotFile,
                                        PlotFormat = PlotFormat,
                                        PlotTitle = description,
                                        Data=numpy.reshape(data_obs, (1,-1)),
                                        Errors=numpy.reshape(data_err, (1,-1)),
                                        DataTrans = DataTrans,
                                        Labels=[],
                                        YLimits = [],
                                        Fontsizes=[Fontsize, Labelsize, Titlesize],
                                        Linewidth= Linewidth,
                                        Markersize = Markersize,
                                        Colors=Colors,
                                        Lines= Lines,
                                        PlotStrng="Error: mult="+str(DatErr_mult)+" add="+str(DatErr_add),
                                        StrngPos = [0.05, 0.95])


                            elif "aem" in AEM_system.lower():
                                viz.plot_site_aem05(
                                        PlotFile = PlotFile,
                                        PlotFormat = PlotFormat,
                                        PlotTitle = description,
                                        Data=numpy.reshape(data_obs, (1,-1)),
                                        Errors=numpy.reshape(data_err, (1,-1)),
                                        Labels=[],
                                        YLimits = [],
                                        Fontsizes=[Fontsize, Labelsize, Titlesize],
                                        Linewidth= Linewidth,
                                        Markersize = Markersize,
                                        Colors=Colors,
                                        Lines= Lines,
                                        PlotStrng="Error: mult="+str(DatErr_mult)+" add="+str(DatErr_add),
                                        StrngPos = [0.05, 0.9])

                            if PdfCatalog:
                                pdf_list.append(OutName+IdString+"_"+str(ncount)+".pdf")


if NPZout:
    NPZFile = OutName+IdString+".npz"
    print("\n\nResults written to "+NPZFile)
    aesys.write_aempy(File=NPZFile, Data=Data, System=AEM_system, Header=Header, OutInfo=OutInfo)


# if PdfCatalog:
#     print("\n\nPDF catalogue: "+PdfCName)
#     viz.make_pdf_catalog(PdfList=pdf_list, FileName=PdfCName)
