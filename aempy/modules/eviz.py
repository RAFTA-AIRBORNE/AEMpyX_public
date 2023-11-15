#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 10:49:39 2023

@author: vrath
"""

import os
import sys
from sys import exit as error
from datetime import datetime
# from time import process_time
# from random import randrange
import time
import warnings
# import inspect
import copy
import numpy

import matplotlib
import matplotlib.pyplot

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse

def plot_model_ensemble(
        ThisAxis = None, 
        PlotFile = None,        
        PlotFormat = ["png",],
        PlotTitle = None,
        PlotType = ["lines"], # lines, Percentiles. iso
        PlotSize = [8.],
        System  = "aem05",
        # ModTrue = [],
        # DepthTrue = [],
        # DatTrue = [],
        ModEns = [],
        Depth = [],
        # DatEns = [],
        Percentiles=[2.5, 16.],
        Quantiles = [.025, .16],
        Fillcolor=["0.8", "0.4"],
        Alphas = [0.3 , 0.6],
        Labels=[],
        Linecolor=["k", "r", "g", "b", "y", "m"],
        Linetype=["-", ":", ";"],
        Linewidth=[1., 1.5, 2.],
        Markers = ["v"],
        Markersize =[4],
        Fontsizes=[10,10,12],
        XLimits= [],
        YLimits=[],
        PlotStrng="",
        Invalid=1.e30,
        Maxlines=30):

    """
    
    """
    cm = 1/2.54  # centimeters to inches
    
    ax = ThisAxis
    
    if numpy.size(ModEns)==0:
       error("No parameter ensemble given!! Exit.")
       # nopar=True

    if ThisAxis==None:
        nplots = 1
        fig, ax = matplotlib.pyplot.subplots(1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[0]*cm),
                                          gridspec_kw={"height_ratios": [1]})
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])
    

    avgval = numpy.mean(ModEns, axis=0)     

    if "per" in PlotType:
        nper = numpy.size(Percentiles)
        medval=numpy.percentile(ModEns, 50., axis=0)
  
        for p in numpy.arange(nper):
                mlow = numpy.percentile(ModEns,      Percentiles[p], axis=0)
                mupp = numpy.percentile(ModEns, 100.-Percentiles[p], axis=0)
                plabel = None
                ax.fill_betweenx(Depth, mlow, mupp, step="post",
                                   linewidth=Linewidth[0],
                                   color= Fillcolor[p],
                                   alpha= Alphas[p],
                                   label=plabel)

        ax.step(medval, Depth,
            linewidth=Linewidth[0], color= Linecolor[1],
            label="medval")

        for p in numpy.arange(nper):
            plabel = "p="\
                +str(100.-Percentiles[p]-Percentiles[p])+" %"
            mlow = numpy.percentile(ModEns,      Percentiles[p], axis=0)
            mupp = numpy.percentile(ModEns, 100.-Percentiles[p], axis=0)
            ax.step(mlow, Depth, where="pre", 
                    linewidth=Linewidth[0], color= Linecolor[p+2])
            ax.step(mupp, Depth, where="pre", 
                    linewidth=Linewidth[0], color= Linecolor[p+2],
                    label=plabel)
    
    if "qua" in PlotType:
        nper = numpy.size(Quantiles)
        medval=numpy.quantile(ModEns, 0.5, axis=0)
  
        for p in numpy.arange(nper):
                mlow = numpy.quantile(ModEns,      Quantiles[p], axis=0)
                mupp = numpy.quantile(ModEns, 1.-Quantiles[p], axis=0)
                plabel = None
                ax.fill_betweenx(Depth, mlow, mupp, step="post",
                                   linewidth=Linewidth[0],
                                   color= Fillcolor[p],
                                   alpha= Alphas[p],
                                   label=plabel)

        ax.step(medval, Depth,
            linewidth=Linewidth[0], color= Linecolor[1],
            label="medval")

        for p in numpy.arange(nper):
            plabel = "q="\
                +str(1.-Quantiles[p]-Quantiles[p])
            mlow = numpy.quantile(ModEns,      Quantiles[p], axis=0)
            mupp = numpy.quantile(ModEns, 1.-Quantiles[p], axis=0)
            ax.step(mlow, Depth, where="pre", 
                    linewidth=Linewidth[0], color= Linecolor[p+2])
            ax.step(mupp, Depth, where="pre", 
                    linewidth=Linewidth[0], color= Linecolor[p+2],
                    label=plabel)
    
    elif "lin" in PlotType:

        nens = numpy.shape(ModEns)
        if nens[0] > Maxlines:
            error("plot_model_ensemble: too many lines! Exit.")
            
        medval=numpy.percentile(ModEns, 50., axis=0)

        for ne in numpy.arange(nens):
            plabel = None

            ax.step(ModEns[1], 
                          where="pre",
                          linewidth=Linewidth[0],
                          color= Fillcolor[ne],
                          alpha= Alphas[ne],
                          label=plabel)

        # ax.step(medval, Depth,step="pre",
        #     linewidth=Linewidth[0]+2, color= Linecolor[1],
        #     label="medval")

    ax.set_xscale("log")
    ax.set_xlim(MLimits)
    ax.set_xlabel("resistivity ($\Omega$m)",fontsize=Fontsizes[0])
    ax.set_ylim(ZLimits)
    ax.set_ylabel("depth (m)",fontsize=Fontsizes[0])
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(labelsize=Fontsizes[1])
    ax.legend(fontsize=Fontsizes[0]-1, loc="best")
    ax.grid(True)
    ax.invert_yaxis()
    ax.grid("major", "both", linestyle=":", lw=0.3)

    if ThisAxis==None:        
        for F in PlotFormat:
            matplotlib.pyplot.savefig(PlotFile+F)
            matplotlib.pyplot.show()
            matplotlib.pyplot.clf()
        
    return ax


def plot_data_ensemble(
        ThisAxis = None, 
        PlotFile = None,        
        PlotFormat = ["png",],
        PlotTitle = None,
        PlotType = ["lines"], # lines, Percentiles. iso
        PlotSize = [8.],
        System  = "aem05",
        DatTrue = [],
        DatEns = [],
        Percentiles=[2.5, 16.],        
        Quantiles = [.025, .16],
        Fillcolor=["0.8", "0.4"],
        Alphas = [0.3 , 0.6],
        Labels=[],
        Linecolor=["k", "r", "g", "b", "y", "m"],
        Linetype=["-", ":", ";"],
        Linewidth=[1., 1.5, 2.],
        Markers = ["v"],
        Markersize =[4],
        Fontsizes=[10,10,12],
        DataTrans=0,
        YLimits=[],
        YLabel = " ppm (-)",
        XLimits=[],        
        XLabel = "frequency (kHz)",
        Invalid=1.e30,
        Maxlines=30):

#     """

    cm = 1/2.54  # centimeters to inches
 
    if numpy.size(DatEns)==0:
        error("No data ensemble given!! Exit.")
        # nodat=True
   
    ax = ThisAxis

    if ThisAxis==None:
        nplots = 1
        fig, ax = matplotlib.pyplot.subplots(1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[0]*cm),
                                          gridspec_kw={"height_ratios": [1]})
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])

    _, NN, _, _, Pars =  aesys.get_system_params(System)

    FAxis = Pars[0]
    AUnit = Pars[1]
    DUnit = Pars[2]
    
    nper = numpy.size(Percentiles)
    
    if "aem05" in System.lower(): 
        XLabel = "frequency (kHz)"
        YLabel = "Q/I (ppm)"
        FAxis = FAxis/1000.
        
        nsmp,ndat = numpy.shape(DatEns)
        
        Qens =   DatEns[:,0:4]        
        Iens =   DatEns[:,4:8]
        
        if "per" in PlotType:
            
            medQens = numpy.percentile(Qens, 50., axis=0)
            for p in numpy.arange(nper):
                
                dlow = numpy.percentile(Qens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Qens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)    
   
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
 
            
            ax.plot(FAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
                  
            for p in numpy.arange(nper):
                 
                medIens = numpy.percentile(Qens, 50., axis=0)   
                dlow = numpy.percentile(Iens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Iens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)      

    
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2],
                        label=plabel)
                
            ax.plot(FAxis, medIens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")            
                
        elif "qua" in PlotType:
            
            medQens = numpy.percentile(Qens, 50., axis=0)
            for p in numpy.arange(nper):
                
                dlow = numpy.percentile(Qens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Qens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)    
   
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
 
            
            ax.plot(FAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
                  
            medIens = numpy.percentile(Iens, 50., axis=0)                   
            for p in numpy.arange(nper):
                 

                dlow = numpy.percentile(Iens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Iens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)      

    
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2],
                        label=plabel)
                
            ax.plot(FAxis, medIens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="I, median")             


        elif "lin" in PlotType:
            
            nens = numpy.shape(Qens)
            if nens[0] > Maxlines:
                error("plot_data_ensemble: too many lines! Exit.")       
            
            plabel = None
            medQens = numpy.percentile(Qens, 50., axis=0) 
            ax.plot(FAxis, Qens,
                        linewidth=Linewidth[0], 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            ax.plot(FAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
            
            medIens = numpy.percentile(Iens, 50., axis=0)  
            ax.plot(FAxis, Iens,
                        linewidth=Linewidth[0], 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            ax.plot(FAxis, medIens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="I, median")           
                
        ax.set_xscale("log")
        if len(XLimits) !=0:
            ax.set_xlim(XLimits)
        ax.set_xlabel(XLabel,fontsize=Fontsizes[0])
        
        if len(YLimits) !=0:
            ax.set_ylim(YLimits)
        ax.set_ylabel(YLabel,fontsize=Fontsizes[0])
        
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("both")
        ax.tick_params(labelsize=Fontsizes[1])
        ax.legend(fontsize=Fontsizes[0]-1, loc="best")
        ax.grid(True)
        ax.invert_yaxis()
        ax.grid("major", "both", linestyle=":", lw=0.3)
        
    elif "gen" in System.lower():
        XLabel = "time (ms)"
       
        
        nsmp,ndat = numpy.shape(DatEns)
        
        H =   DatEns[:,0:11]        
        Z =   DatEns[:,11:22]
        
        if DataTrans==2:
           S = numpy.amin(DatEns)        
           H = H, _, _ = inverse.transform_data(d_vec=H, d_trn=DataTrans, scale=S)
           Z = Z, _, _ = inverse.transform_data(d_vec=Z, d_trn=DataTrans, scale=S)
           
           YLabel = "atanh H/Z (-)"
           ax.set_yscale("linear")
           
        elif DataTrans==1:
           YLabel = "H/Z (ppm)"              
           ax.set_yscale("symlog")
        else:
           YLabel = "H/Z (ppm)"   
           ax.set_yscale("log", nonpositive="clip")

        if "per" in PlotType:
            
            medQens = numpy.percentile(H, 50., axis=0)
            for p in numpy.arange(nper):
                
                dlow = numpy.percentile(H,      Percentiles[p], axis=0)
                dupp = numpy.percentile(H, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)    
   
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
 
            
            ax.plot(FAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
                  
            for p in numpy.arange(nper):
                 
                medIens = numpy.percentile(Qens, 50., axis=0)   
                dlow = numpy.percentile(Iens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Iens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)      

    
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2],
                        label=plabel)
                
            ax.plot(FAxis, medIens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")            
                
        elif "qua" in PlotType:
            
            medQens = numpy.percentile(Qens, 50., axis=0)
            for p in numpy.arange(nper):
                
                dlow = numpy.percentile(Qens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Qens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)    
   
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
 
            
            ax.plot(FAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
                  
            medIens = numpy.percentile(Iens, 50., axis=0)                   
            for p in numpy.arange(nper):
                 

                dlow = numpy.percentile(Iens,      Percentiles[p], axis=0)
                dupp = numpy.percentile(Iens, 100.-Percentiles[p], axis=0)
        
                plabel = None
                ax.fill_between(FAxis, dlow, dupp, 
                                        linewidth=Linewidth[0],
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=plabel)      

    
                plabel = "p="\
                    +str(100.-Percentiles[p]-Percentiles[p])+" %"
    
                ax.plot(FAxis, dlow,
                        linewidth=Linewidth[0], color= Linecolor[p+2])
                ax.plot(FAxis, dupp,
                        linewidth=Linewidth[0], color= Linecolor[p+2],
                        label=plabel)
                
            ax.plot(FAxis, medIens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="I, median")             


        elif "lin" in PlotType:
            
            nens = numpy.shape(Qens)
            if nens[0] > Maxlines:
                error("plot_data_ensemble: too many lines! Exit.")       
            
            plabel = None
            medQens = numpy.percentile(Qens, 50., axis=0) 
            ax.plot(FAxis, Qens,
                        linewidth=Linewidth[0], 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            ax.plot(FAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
            
            medIens = numpy.percentile(Iens, 50., axis=0)  
            ax.plot(FAxis, Iens,
                        linewidth=Linewidth[0], 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            ax.plot(FAxis, medIens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="I, median")           
 


        
        ax.set_xscale("log")
        if len(XLimits) !=0:
            ax.set_xlim(XLimits)
        ax.set_xlabel(XLabel,fontsize=Fontsizes[0])
         
        if len(YLimits) !=0:
            ax.set_ylim(YLimits)
        ax.set_ylabel(YLabel,fontsize=Fontsizes[0])


        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("both")
        ax.tick_params(labelsize=Fontsizes[1])
        ax.legend(fontsize=Fontsizes[0]-1, loc="best")
        ax.grid(True)
        ax.invert_yaxis()
        ax.grid("major", "both", linestyle=":", lw=0.3)


    if ThisAxis==None:
        for F in PlotFormat:
            matplotlib.pyplot.savefig(PlotFile+F)      
            matplotlib.pyplot.show()
            matplotlib.pyplot.clf()
            
    return ax
