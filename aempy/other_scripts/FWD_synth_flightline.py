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
Show several 1d block models as (stitched) section.

"""
import os
import sys
from sys import exit as error
from datetime import datetime
import warnings

import numpy
import scipy

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.colors as col
import matplotlib.pyplot
import matplotlib

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
# mypath = ["/home/vrath/AEMpyX/aempy/modules/", "/home/vrath/AEMpyX/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)

from version import versionstrg


import util
import prep
import aesys
import core1d
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
print("Forward modelling of full flightline"+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
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
AEM_system_out = "genesis"
FwdCall_out,NN_out, _, _, _, = aesys.get_system_params(System=AEM_system_out)
DataTrans_out = 0

AEM_system_fd = "aem05"
FwdCall_fd,NN_fd, _, _, _, = aesys.get_system_params(System=AEM_system_fd)
dat_act = numpy.ones(NN_out[2], dtype="int8")

AEM_system_td = "aem05"
FwdCall_td,NN_td, _, _, _, = aesys.get_system_params(System=AEM_system_td)
dat_act = numpy.ones(NN_out[2], dtype="int8")


OutInfo = True
now = datetime.now()


InModDir = AEMPYX_DATA + "/Nearest/fwd_compare/models/"
print("Input model read from dir:  %s" % InModDir)
mod_files = [
            "A2_NM_intersection_FL21513-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21514-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",
            "A1_NM_intersection_FL11126-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21506-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21520-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21524-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21618-0_delete_PLM10_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k5_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k4_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k3_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k2_nlyr32_TikhOpt_gcv_Results.npz",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k1_nlyr32_TikhOpt_gcv_Results.npz",]

InDatDir = AEMPYX_DATA + "/Nearest/fwd_compare/data/"
print("Input data read from dir:  %s" % InDatDir)
dat_files_fd = [
            "A2_NM_intersection_FL21513-0_delete_PLM10.asc",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k5.asc",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k4.asc",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k3.asc",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k2.asc",
            "A2_NM_intersection_FL21513-0_delete_PLM10_k1.asc",
            "A2_NM_intersection_FL21514-0_delete_PLM10.asc",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k5.asc",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k4.asc",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k3.asc",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k2.asc",
            "A2_NM_intersection_FL21514-0_delete_PLM10_k1.asc",
            "A1_NM_intersection_FL11126-0_delete_PLM10.asc",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k5.asc",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k4.asc",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k3.asc",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k2.asc",
            "A1_NM_intersection_FL11126-0_delete_PLM10_k1.asc",
            "A2_NM_intersection_FL21506-0_delete_PLM10.asc",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k5.asc",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k4.asc",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k3.asc",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k2.asc",
            "A2_NM_intersection_FL21506-0_delete_PLM10_k1.asc",
            "A2_NM_intersection_FL21520-0_delete_PLM10.asc",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k5.asc",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k4.asc",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k3.asc",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k2.asc",
            "A2_NM_intersection_FL21520-0_delete_PLM10_k1.asc",
            "A2_NM_intersection_FL21524-0_delete_PLM10.asc",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k5.asc",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k4.asc",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k3.asc",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k2.asc",
            "A2_NM_intersection_FL21524-0_delete_PLM10_k1.asc",
            "A2_NM_intersection_FL21618-0_delete_PLM10.asc",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k5.asc",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k4.asc",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k3.asc",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k2.asc",
            "A2_NM_intersection_FL21618-0_delete_PLM10_k1.asc",]

dat_files_td = [
            "NM_A2_intersection_FL10420-0.asc",
            "NM_A2_intersection_FL10420-0.asc",
            "NM_A2_intersection_FL10420-0.asc",
            "NM_A2_intersection_FL10420-0.asc",
            "NM_A2_intersection_FL10420-0.asc",
            "NM_A2_intersection_FL10420-0.asc",
            "NM_A2_intersection_FL10430-0.asc",
            "NM_A2_intersection_FL10430-0.asc",
            "NM_A2_intersection_FL10430-0.asc",
            "NM_A2_intersection_FL10430-0.asc",
            "NM_A2_intersection_FL10430-0.asc",
            "NM_A2_intersection_FL10430-0.asc",
            "NM_A1_intersection_FL13490-0.asc",
            "NM_A1_intersection_FL13490-0.asc",
            "NM_A1_intersection_FL13490-0.asc",
            "NM_A1_intersection_FL13490-0.asc",
            "NM_A1_intersection_FL13490-0.asc",
            "NM_A1_intersection_FL13490-0.asc",
            "NM_A2_intersection_FL10350-0.asc",
            "NM_A2_intersection_FL10350-0.asc",
            "NM_A2_intersection_FL10350-0.asc",
            "NM_A2_intersection_FL10350-0.asc",
            "NM_A2_intersection_FL10350-0.asc",
            "NM_A2_intersection_FL10350-0.asc",
            "NM_A2_intersection_FL10491-0.asc",
            "NM_A2_intersection_FL10491-0.asc",
            "NM_A2_intersection_FL10491-0.asc",
            "NM_A2_intersection_FL10491-0.asc",
            "NM_A2_intersection_FL10491-0.asc",
            "NM_A2_intersection_FL10491-0.asc",
            "NM_A2_intersection_FL10530-0.asc",
            "NM_A2_intersection_FL10530-0.asc",
            "NM_A2_intersection_FL10530-0.asc",
            "NM_A2_intersection_FL10530-0.asc",
            "NM_A2_intersection_FL10530-0.asc",
            "NM_A2_intersection_FL10530-0.asc",
            "NM_A2_intersection_FL11470-0.asc",
            "NM_A2_intersection_FL11470-0.asc",
            "NM_A2_intersection_FL11470-0.asc",
            "NM_A2_intersection_FL11470-0.asc",
            "NM_A2_intersection_FL11470-0.asc",
            "NM_A2_intersection_FL11470-0.asc",]


ns = numpy.size(dat_files_fd)


for ifl in numpy.arange(ns):

    print("\n")
    mod_file = InModDir+mod_files[ifl]
    print("Model read from: "+mod_file)
    fd_file = InDatDir+dat_files_fd[ifl]
    print("Data read from: "+fd_file)


    fn, fext = os.path.splitext(fd_file)

    dat_fd, head = aesys.read_aempy(File=fd_file,
                                   System=AEM_system_fd, OutInfo=False)




    mod = numpy.load(mod_file)
    mod_sit = mod["site_modl"]
    mod_ref = mod["mod_ref"]
    mod_act = mod["mod_act"]

    site0_x = mod["site_x"]
    site0_y = mod["site_y"]
    site0_alt = mod["site_alt"]

    d_act=mod["dat_act"]


    alt_file = InDatDir+dat_files_td[ifl]

    site1, _= aesys.read_aempy(File=alt_file,
                                   System=AEM_system_td, OutInfo=False)
    site1_x = site1[:,1]
    site1_y = site1[:,2]
    site1_alt = site1[:,4]

    site_alt = numpy.full_like(site0_alt, numpy.nan)
    for site0 in numpy.arange(numpy.size(site0_x)):
        x0 = site0_x[site0]
        y0 = site0_y[site0]
        dist = numpy.full_like(site1_alt, numpy.nan)

        for site1 in numpy.arange(numpy.size(site1_x)):
            x1 = site1_x[site1]
            y1 = site1_y[site1]
            dist[site1] = numpy.sqrt((x1-x0)**2+(y1-y0)**2)

        minpos = numpy.argmin(dist)

        site_alt[site0]=site1_alt[minpos]


    nlyr = inverse.get_nlyr(mod_ref)
    sites, param = numpy.shape(mod_sit)

    for site in numpy.arange(sites):
        """
        calculate forward model
        """

        d_state = 0
        m_state = 0
        m_sit =mod_sit[site,:]
        #!!!!!! only workaround
        alt = site_alt[site]



        m_vec = inverse.insert_mod(M=mod_ref, m = m_sit, m_act=mod_act)
        d_current, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall_out, alt=alt,
                                          m_vec = m_vec, d_trn=DataTrans_out)
        if site == 0:
            d_new = d_current
        else:
            d_new = numpy.vstack(( d_new, d_current))

    dum = numpy.zeros((numpy.shape(dat_fd)[0],2))
    # D  = numpy.hstack((dat_fd[:,0:6], d_new, dum))
    D  = numpy.concatenate((dat_fd[:,0:3], site_alt.reshape(-1,1), dat_fd[:,5:7], d_new, dum), axis =1)
    dat_new = fn+"_genesis_225Hz_data "+fext
    print("Data written to: "+dat_new)
    aesys.write_aempy(File=dat_new, Data=D, System=AEM_system_out,
                    Header=head, OutInfo=False)
