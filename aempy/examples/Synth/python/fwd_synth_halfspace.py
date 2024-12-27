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


import matplotlib.collections
import matplotlib.patches
import matplotlib.colors
import matplotlib.pyplot
import matplotlib
import matplotlib.cm


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
# -


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

# +
rng = numpy.random.default_rng()
nan = numpy.nan
cm = 1/2.54


version, _ = versionstrg()
#script = "Tutorial0_FWD_synth.py"
script = __file__  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")
Header = titstrng
# -

OutInfo = False
AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/Synth/"
PlotDir = AEMPYX_DATA+"/plots/"
print("Plots written to dir: %s " % PlotDir)
if not os.path.isdir(PlotDir):
    print("File: %s does not exist, but will be created" % PlotDir)
    os.mkdir(PlotDir)
# The following cell gives values to AEM-system related settings.
#
# Data transformation is activated by the variable DataTrans. Currently three possible options are allowed: _DataTrans = 0_: No transformation, i.e., the raw data are used. _DataTrans = 1_: The natural log of data is taken, only allowed for strictly positive values. _DataTrans = 2_: If data scale logarithmically, an asinh transformation (introduced by Scholl, 2000) is applied. It allows negatives, which may occur in TDEM, when IP effects are present.
#
# A general additive/multiplicative error model is applied on the raw data before transformation, and errors are also transformed.

# +
AEM_system = "aem05"
# AEM_system = "genesis"
print("AEM system: " + AEM_system + "\n \n")

FwdCall,NN, _, _, Pars, = aesys.get_system_params(System=AEM_system)
ParaTrans = 0
DataTrans = 0
DatErr_add = 50.
DatErr_mult = 0.00
alt = 60.
DataActive = numpy.ones((1,NN[2]))

nD = NN[0]
# -

# To initialize loops over different parameters,
# first a reference model must be set up, with reasonable values for all parameters not within the loop. Default settings is rho only, no IP. Currently, one parameter and altitude can be varied within a loop.
# The following should be adapted according to the user's needs.

nlyr = 1
Model_active, Model_base, model_var, m_bounds, m_state = inverse.init_1dmod(nlyr)

Model_base[0*nlyr:1*nlyr] =[100.]           #rho
Model_base[6*nlyr:7*nlyr-1] =[]          #layers

# +
# Adapted for reasonable IP values

# Model_base[3*nlyr:4*nlyr] =[0.,  0.5, 0.]      #chargeability
# Model_base[4*nlyr:5*nlyr] =[0.,  0.5, 0.]      #exponent
# Model_base[5*nlyr:6*nlyr] =[0., 100., 0.]      #frequency
# -

# rho for layer 1 (starting from 0!)
FWDBaseName = "AEM05_RhoDepend_30m"
VarPar = [ 1., 3.,  10., 30., 100.,300., 1000.,3000., 10000.]
VarInd = 0

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
Alt = [30]
# -

# Now generate the response data:

data = []
mod_num = -1
for par in numpy.arange(len(VarPar)):

        mod_num += 1

        m_i = Model_base.copy()

        if VarInd==numpy.size(m_i):
            alt = Alt[par]
            # p_i = numpy.array([mod_num, VarInd, Alt[par], DataTrans, DatErr_add, DatErr_mult])


        else:
            m_i[VarInd] = VarPar[par]
            # p_i = numpy.array([mod_num, VarInd, VarPar[par], DataTrans, DatErr_add, DatErr_mult])


        d_state = 0
        m_state = 0

        m_current = m_i
        # m_current, m_state = inverse.transform_parameter(m_vec=m_i, m_trn=ParaTrans, m_state=m_state, mode="f")
        print("\n",par)
        print("modl:",m_current)
        d_ref, d_state = inverse.calc_fwdmodel(fwdcall=FwdCall, alt=Alt,
                                          m_vec = m_current, m_trn=ParaTrans, m_state=m_state,
                                          d_trn=0, d_state=d_state, d_act = DataActive )

        print("data:", d_ref)

        d = numpy.insert(d_ref, 0, VarPar[par])
        if par==0:
            data = d
        else:
            data = numpy.vstack((data, d))
        # print(numpy.shape(data))



PlotTitle = "AEM05: Halfspaces at "+str(int(alt))+" m"
PlotSize = [8., 8.]
FilesOnly = False
PlotFormat = [".pdf", ".png"]

"""
Determine graphical parameter.
=> print(matplotlib.pyplot.style.available)
see:
MatplotlibDeprecationWarning: The seaborn styles shipped by Matplotlib
are deprecated since 3.6, as they no longer correspond to the styles s
hipped by seaborn. However, they will remain available as
'seaborn-v0_8-<style>'. Alternatively, directly use the seaborn API instead.

"""
matplotlib.pyplot.style.use("seaborn-v0_8-paper") # ("seaborn-paper")
matplotlib.rcParams["figure.dpi"] = 400
matplotlib.rcParams["axes.linewidth"] = 0.5
matplotlib.rcParams["savefig.facecolor"] = "none"
# matplotlib.rcParams["text.usetex"] = True

Fontsize = 8
Labelsize = Fontsize
Titlesize = 12
Fontsizes = [Fontsize, Labelsize, Titlesize]

Linewidth = [1., 1., 0.75]
Linetypes = ["-", ":", "--", ";"]
Linecolors = ["k", "r", "g", "b", "c", "m"]

Markers = ["o"]
Markersize = [3]


ncols = len(VarPar)+3

Greys = [0.1, 0.2, 0.3, 0.4, 0.5,.6, 0.7, 0.8, 0.9]

"""
see:
https://matplotlib.org/stable/gallery/color/colormap_reference.html
"""
cmap = matplotlib.cm.gist_rainbow
mycmap = matplotlib.pyplot.get_cmap(cmap)
"""
For just plotming to files, choose the cairo backend (eps, pdf, ,png, jpg...).
If you need to see the plot directly (plot window, or jupyter), simply
comment out the following line. In this case matplotlib may run into
memory problems after a few hundreds of high-resolution plot.
"""
if FilesOnly:
   matplotlib.use("cairo")


nplots = 2
fig, ax = matplotlib.pyplot.subplots(1,nplots,
                           figsize=(2*PlotSize[0]*cm, PlotSize[0]*cm),
                           gridspec_kw={
                               "height_ratios": [1.],
                               "width_ratios": [1., 1.]})

fig.suptitle(PlotTitle, fontsize=Fontsizes[2])



freqs = Pars[0]/1000.


nsmp,ndat = numpy.shape(data)

pure =  data[:,1:]

qens =   pure[:,0:4]
iens =   pure[:,4:8]

# print(numpy.shape(data), numpy.shape(qens), numpy.shape(iens) )

nens = numpy.shape(qens)


for var in numpy.arange(nens[0]):

    ax[0].plot(freqs, qens[var,:],
                linewidth=Linewidth[0],marker=Markers[0], markersize=Markersize[0],
                label=str(VarPar[var]))

ax[0].set_xscale("log")
ax[0].set_xlabel("freqency (kHz)",fontsize=Fontsizes[0])
ax[0].set_ylabel("quadrature (ppm)",fontsize=Fontsizes[0])

ax[0].xaxis.set_label_position("bottom")
ax[0].xaxis.set_ticks_position("both")
ax[0].tick_params(labelsize=Fontsizes[1])
ax[0].grid(True)
ax[0].grid("major", "both", linestyle=":", lw=0.3)
ax[0].legend(fontsize=Fontsizes[0]-2, loc="best", title="Ohm.m", title_fontsize=Fontsizes[0]-2)


for  var in numpy.arange(nens[0]):
    ax[1].plot(freqs,iens[var,:],
                    linewidth=Linewidth[0], marker=Markers[0], markersize=Markersize[0],
                label=str(VarPar[var]))
ax[1].set_xscale("log")
ax[1].set_xlabel("freqency (kHz)",fontsize=Fontsizes[0])
ax[1].set_ylabel("in-phase (ppm)",fontsize=Fontsizes[0])

ax[1].xaxis.set_label_position("bottom")
ax[1].xaxis.set_ticks_position("both")
ax[1].tick_params(labelsize=Fontsizes[1])
ax[1].grid(True)
ax[1].grid("major", "both", linestyle=":", lw=0.3)
ax[1].legend(fontsize=Fontsizes[0]-2, loc="best", title="Ohm.m", title_fontsize=Fontsizes[0]-2)
fig.tight_layout()
for F in PlotFormat:
    matplotlib.pyplot.savefig(PlotDir+FWDBaseName+F)
