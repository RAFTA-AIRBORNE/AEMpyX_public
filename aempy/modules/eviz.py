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
        PlotFile = "",        
        PlotFormat = ["png",],
        PlotTitle = "",
        PlotType = ["lines"], # lines, Percentiles. iso
        PlotSize = [8.],
        System  = "aem05",
        ModEns = [],
        Depth = [],
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
        Maxlines=50,
        Median=True,
        Legend=True):

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
        Legend = True
    
     

    if "per" in PlotType.lower():
        np = numpy.size(Percentiles)
        medval=numpy.percentile(ModEns, 50., axis=0)
  
        for p in numpy.arange(np):
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

        for p in numpy.arange(np):
            plabel = "p="\
                +str(100.-Percentiles[p]-Percentiles[p])+" %"
            mlow = numpy.percentile(ModEns,      Percentiles[p], axis=0)
            mupp = numpy.percentile(ModEns, 100.-Percentiles[p], axis=0)
            ax.step(mlow, Depth, where="pre", 
                    linewidth=Linewidth[0]/2, color= Linecolor[p+2])
            ax.step(mupp, Depth, where="pre", 
                    linewidth=Linewidth[0]/2, color= Linecolor[p+2],
                    label=plabel)
    
    if "qua" in PlotType.lower():
        np = numpy.size(Quantiles)
        medval=numpy.quantile(ModEns, 0.5, axis=0)
  
        for p in numpy.arange(np):
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

        for p in numpy.arange(np):
            plabel = "q="\
                +str(1.-Quantiles[p]-Quantiles[p])
            mlow = numpy.quantile(ModEns,      Quantiles[p], axis=0)
            mupp = numpy.quantile(ModEns, 1.-Quantiles[p], axis=0)
            ax.step(mlow, Depth, where="pre", 
                    linewidth=Linewidth[0]/2, color= Linecolor[p+2])
            ax.step(mupp, Depth, where="pre", 
                    linewidth=Linewidth[0]/2, color= Linecolor[p+2],
                    label=plabel)
    
    elif "lin" in PlotType.lower():

        nens = numpy.shape(ModEns)
        if nens[0] > Maxlines:
            error("plot_model_ensemble: too many lines! Exit.")
            
        medval=numpy.percentile(ModEns, 50., axis=0)

        for ne in numpy.arange(nens):
            plabel = None

            ax.step(ModEns[1], Depth,
                          where="pre",
                          linewidth=Linewidth[0]/2,
                          color= Fillcolor[ne],
                          alpha= Alphas[ne],
                          label=plabel)

        # ax.step(medval, Depth,step="pre",
        #     linewidth=Linewidth[0]+2, color= Linecolor[1],
        #     label="medval")

    ax.set_xscale("log")
    ax.set_xlim(XLimits)
    ax.set_xlabel("resistivity ($\Omega$m)",fontsize=Fontsizes[0])
    ax.set_ylim(YLimits)
    ax.set_ylabel("depth (m)",fontsize=Fontsizes[0])
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(labelsize=Fontsizes[1])

    ax.grid(True)
    ax.invert_yaxis()
    ax.grid("major", "both", linestyle=":", lw=0.3)

    if Legend:
        ax.legend(fontsize=Fontsizes[0], loc="best")

    if ThisAxis==None:        
        for F in PlotFormat:
            matplotlib.pyplot.savefig(PlotFile+F)
        # matplotlib.pyplot.show()
        # matplotlib.pyplot.clf()      
        
    return ax


def plot_data_ensemble(
        ThisAxis = None, 
        PlotFile = "",        
        PlotFormat = [".png",],
        PlotTitle = "",
        PlotType = ["lines"], # lines, Percentiles. iso
        PlotSize = [8.],
        System  = "aem05",
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
        Maxlines=30,
        Legend=True,
        Median=True):

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
        Legend = True

    _, NN, _, _, Pars =  aesys.get_system_params(System)

    
    np = numpy.size(Percentiles)
    
    if "aem05" in System.lower(): 
        XLabel = "frequency (kHz)"
        YLabel = "Q/I (ppm)"
        XAxis = Pars[0]/1000.
        
        
        nsmp,ndat = numpy.shape(DatEns)
        
        Qens =   DatEns[:,0:4]        
        Iens =   DatEns[:,4:8]
        
        if ("per" in PlotType.lower()) or ("qua" in PlotType.lower()):
            

            for p in numpy.arange(np):
                
                
                
                if "per" in PlotType.lower():
                    medQens = numpy.percentile(Qens, 50., axis=0)
                    plabel = "p="\
                        +str(100.-Percentiles[p]-Percentiles[p])+" %"
                   
                    dlow = numpy.percentile(Qens,      Percentiles[p], axis=0)
                    dupp = numpy.percentile(Qens, 100.-Percentiles[p], axis=0)
                
                else:
                    medQens = numpy.percentile(Qens, 0.5, axis=0)
                    plabel = "p="\
                        +str(1.-Quantiles[p]-Quantiles[p])+" %"
                    
                    dlow = numpy.quantile(Qens,    Quantiles[p], axis=0)
                    dupp = numpy.quantile(Qens, 1.-Quantiles[p], axis=0)

                
                
                ax.fill_between(XAxis, dlow, dupp, 
                                        linewidth=Linewidth[0]/2,
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=None)    
       
                ax.plot(XAxis, dlow,
                        linewidth=Linewidth[0]/2, color= Linecolor[p+2])
                ax.plot(XAxis, dupp,
                        linewidth=Linewidth[0]/2, color= Linecolor[p+2])
     
            
            if Median:
                ax.plot(XAxis, medQens,
                        linewidth=Linewidth[0], color= Linecolor[0], linestyle=Linetype[2],
                        label="medval")
                  

            for p in numpy.arange(np):
                 
                if "per" in PlotType.lower():
                    medIens = numpy.percentile(Iens, 50., axis=0)
                    plabel = "p="\
                        +str(100.-Percentiles[p]-Percentiles[p])+" %"
                   
                    dlow = numpy.percentile(Iens,      Percentiles[p], axis=0)
                    dupp = numpy.percentile(Iens, 100.-Percentiles[p], axis=0)
                
                else:
                    medIens = numpy.percentile(Iens, 0.5, axis=0)
                    plabel = "p="\
                        +str(1.-Quantiles[p]-Quantiles[p])+" %"
                    
                    dlow = numpy.quantile(Iens,    Quantiles[p], axis=0)
                    dupp = numpy.quantile(Iens, 1.-Quantiles[p], axis=0)

                ax.fill_between(XAxis, dlow, dupp, 
                                        linewidth=Linewidth[0]/2,
                                        color= Fillcolor[p],
                                        alpha= Alphas[p],
                                        label=None)      

    
    
                ax.plot(XAxis, dlow,
                        linewidth=Linewidth[0]/2, color= Linecolor[p+2])
                ax.plot(XAxis, dupp,
                        linewidth=Linewidth[0]/2, color= Linecolor[p+2],
                        label=plabel)
            if Median:           
                ax.plot(XAxis, medIens,
                        linewidth=Linewidth[0], color= Linecolor[0], linestyle=Linetype[2],
                        label=None)            
                    
           


        elif "lin" in PlotType.lower():
            
            nens = numpy.shape(Qens)
            if nens[0] > Maxlines:
                error("plot_data_ensemble: too many lines! Exit.")       
            
            plabel = None
            medQens = numpy.percentile(Qens, 50., axis=0) 
            ax.plot(XAxis, Qens,
                        linewidth=Linewidth[0]/2, 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            
            if Median: 
                ax.plot(XAxis, medQens,
                    linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                    label="Q, median")
            
            medIens = numpy.percentile(Iens, 50., axis=0)  
            ax.plot(XAxis, Iens,
                        linewidth=Linewidth[0]/2, 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            if Median:
                ax.plot(XAxis, medIens,
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
        
        ax.grid(True)
        # ax.invert_yaxis()
        ax.grid("major", "both", linestyle=":", lw=0.3)
        
        if Legend:
            ax.legend(fontsize=Fontsizes[0], loc="best")
        
        
    elif "gen" in System.lower():
        XLabel = "time (ms)"
        YLabel = "ppm (-)"
        XAxis = Pars[0]
        
        nsmp,ndat = numpy.shape(DatEns)
        
        Hens =   DatEns[:,0:11]        
        Zens =   DatEns[:,11:22]
        
        if DataTrans==2:
           S = numpy.amin(DatEns)        
           Hens = Hens, _, _ = inverse.transform_data(d_vec=Hens, d_trn=DataTrans, scale=S)
           Zens = Zens, _, _ = inverse.transform_data(d_vec=Zens, d_trn=DataTrans, scale=S)
           
           YLabel = "atanh H/Z (-)"
           ax.set_yscale("linear")
           
        elif DataTrans==1:
           YLabel = "H/Z (ppm)"              
           ax.set_yscale("symlog")
        else:
           YLabel = "H/Z (ppm)"   
           ax.set_yscale("log", nonpositive="clip")

        
        if ("per" in PlotType.lower()) or ("qua" in PlotType.lower()):
            


                for p in numpy.arange(np):
                    
                    if "per" in PlotType.lower():
                        medHens = numpy.percentile(Hens, 50., axis=0)
                        plabel = "p="\
                            +str(100.-Percentiles[p]-Percentiles[p])+" %"
                       
                        dlow = numpy.percentile(Hens,      Percentiles[p], axis=0)
                        dupp = numpy.percentile(Hens, 100.-Percentiles[p], axis=0)
                    
                    else:
                        medHens = numpy.percentile(Hens, 0.5, axis=0)
                        plabel = "p="\
                            +str(1.-Quantiles[p]-Quantiles[p])
                        
                        dlow = numpy.quantile(Hens,    Quantiles[p], axis=0)
                        dupp = numpy.quantile(Hens, 1.-Quantiles[p], axis=0)
                
            
                    ax.fill_between(XAxis, dlow, dupp, 
                                            linewidth=Linewidth[0],
                                            color= Fillcolor[p],
                                            alpha= Alphas[p],
                                            label=None)    
       
        
                    ax.plot(XAxis, dlow,
                            linewidth=Linewidth[0]/2, color= Linecolor[p+2])
                    ax.plot(XAxis, dupp,
                            linewidth=Linewidth[0]/2, color= Linecolor[p+2])
     
                if Median:
                    ax.plot(XAxis, medHens,
                            linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                            label="medval")
                          

                for p in numpy.arange(np):
                     
                    if "per" in PlotType.lower():
                        medZens = numpy.percentile(Zens, 50., axis=0)
                        plabel = "p="\
                            +str(100.-Percentiles[p]-Percentiles[p])+" %"
                       
                        dlow = numpy.percentile(Zens,      Percentiles[p], axis=0)
                        dupp = numpy.percentile(Zens, 100.-Percentiles[p], axis=0)
                    
                    else:
                        medZens = numpy.percentile(Zens, 0.5, axis=0)
                        plabel = "p="\
                            +str(1.-Quantiles[p]-Quantiles[p])
                        
                        dlow = numpy.quantile(Zens,    Quantiles[p], axis=0)
                        dupp = numpy.quantile(Zens, 1.-Quantiles[p], axis=0)
            
                    plabel = None
                    ax.fill_between(XAxis, dlow, dupp, 
                                            linewidth=Linewidth[0]/2,
                                            color= Fillcolor[p],
                                            alpha= Alphas[p],
                                            label=None)      
    
        
        
                    ax.plot(XAxis, dlow,
                            linewidth=Linewidth[0]/2, color= Linecolor[p+2])
                    ax.plot(XAxis, dupp,
                            linewidth=Linewidth[0]/2, color= Linecolor[p+2],
                            label=plabel)

                if Median:
                    ax.plot(XAxis, medIens,
                            linewidth=Linewidth[0]/2, color= Linecolor[2], linestyle=Linetype[2],
                            label=None)            
       


        elif "lin" in PlotType.lower():
            
            nens = numpy.shape(Hens)
            if nens[0] > Maxlines:
                error("plot_data_ensemble: too many lines! Exit.")       
            
            plabel = None
            medHens = numpy.percentile(Hens, 50., axis=0) 
            ax.plot(XAxis, Hens,
                        linewidth=Linewidth[0]/2, 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            if Median:
                ax.plot(XAxis, medHens,
                        linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                        label="H, median")
                
            medZens = numpy.percentile(Zens, 50., axis=0)  
            ax.plot(XAxis, Zens,
                        linewidth=Linewidth[0]/2, 
                        color= Linecolor[0], alpha=0.5,
                        label=plabel)
            if Median:
                ax.plot(XAxis, medZens,
                        linewidth=Linewidth[0], color= Linecolor[2], linestyle=Linetype[2],
                        )           
  

        
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
 
        ax.grid(True)
        ax.grid("major", "both", linestyle=":", lw=0.3)
        
        if Legend:
            ax.legend(fontsize=Fontsizes[0], loc="best")

    if ThisAxis==None:
        for F in PlotFormat:
            matplotlib.pyplot.savefig(PlotFile+F)      
    # matplotlib.pyplot.show()
    # matplotlib.pyplot.clf()
            
    return ax

def plot_data(
        ThisAxis = None, 
        PlotFile = "",        
        PlotFormat = ["png",],
        PlotTitle = "",
        PlotSize = [8.],
        System  = "aem05",
        Data = [],
        Errs = [],
        Labels=[],
        Linecolor=["k", "r", "g", "b", "c", "m"],
        Linetype=["-", ":", ";"],
        Linewidth=[1., 1.5, 2.],
        Markers = ["o", "s"],
        Markersize =[4],
        Fontsizes=[10,10,12],
        DataTrans=0,
        YLimits=[],
        YLabel = " ppm (-)",
        XLimits=[],        
        XLabel = "frequency (kHz)",
        Invalid=1.e30,
        Maxlines=30,
        Legend=True):

    cm = 1/2.54  # centimeters to inches
 
    if numpy.size(Data)==0:
        error("No data ensemble given!! Exit.")
        # nodat=True
   
    ax = ThisAxis

    if ThisAxis==None:
        nplots = 1
        fig, ax = matplotlib.pyplot.subplots(1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[0]*cm),
                                          gridspec_kw={"height_ratios": [1]})
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])
        Legend = True

    _, NN, _, _, Pars =  aesys.get_system_params(System)

    
    
    if "aem05" in System.lower(): 
        XLabel = "frequency (kHz)"
        YLabel = "Q/I (ppm)"
        XAxis = Pars[0]/1000.
        
                
        Qd =   Data[0:4]        
           
        if len(Errs) == 0:
          
            ax.plot(XAxis, Qd,
                         linestyle=Linetype[0], marker=Markers[0],
                         color=Linecolor[1], linewidth=Linewidth[0],
                         markersize=Markersize[0], label="Q, observed" )                
        else:
            Qe =   Errs[:,0:4]        
            ax.errorbar(XAxis, Qd, yerr=Qe,
                         linestyle=Linetype[0], marker=Markers[0],
                         color=Linecolor[1], linewidth=Linewidth[0],
                         markersize=Markersize[0], label="Q, observed" )

        Id =   Data[4:8]
        
        if len(Errs) == 0:
          
            ax.plot(XAxis, Id,
                         linestyle=Linetype[0], marker=Markers[1],
                         color=Linecolor[2],linewidth=Linewidth[0],
                         markersize=Markersize[0], label="I, observed" )                
        else:
            Ie =   Errs[:,4:8]        
            ax.errorbar(XAxis, Id, yerr=Ie,
                         linestyle=Linetype[0], marker=Markers[1],
                         color=Linecolor[2], linewidth=Linewidth[0],
                         markersize=Markersize[0], label="I, observed" )
 
    if "gene" in System.lower(): 
        
        XLabel = "time (ms)"
        YLabel = "H/Z (ppm)"
        XAxis = Pars[0]
        
        Hd =   Data[0:11]        
        Zd =   Data[11:22]
        
        if DataTrans==2:
             S = numpy.amin(Data)        
             Hd = Hd, _, _ = inverse.transform_data(d_vec=Hd, d_trn=DataTrans, scale=S)
             Zd = Zd, _, _ = inverse.transform_data(d_vec=Zd, d_trn=DataTrans, scale=S)
             
             YLabel = "atanh H/Z (-)"
             ax.set_yscale("linear")
             
        elif DataTrans==1:
           YLabel = "H/Z (ppm)"              
           ax.set_yscale("symlog")
        else:
           YLabel = "H/Z (ppm)"   
           ax.set_yscale("log", nonpositive="clip")
  


        if len(Errs) == 0:
          
            ax.plot(XAxis, Hd,
                         linestyle=Linetype[0], marker=Markers[0],
                         color=Linecolor[1], linewidth=Linewidth[0],
                         markersize=Markersize[0], label="H, observed" )                
        else:
            He =   Errs[:,0:4]        
            ax.errorbar(XAxis, Hd, yerr=He,
                         linestyle=Linetype[0], marker=Markers[0],
                         color=Linecolor[1], linewidth=Linewidth[0],
                         markersize=Markersize[0], label="H, observed" )

      
        if len(Errs) == 0:
          
            ax.plot(XAxis, Zd,
                         linestyle=Linetype[0], marker=Markers[1],
                         color=Linecolor[2],linewidth=Linewidth[0],
                         markersize=Markersize[0], label="Z, observed" )                
        else:
            Ze =   Errs[:,4:8]        
            ax.errorbar(XAxis, Zd, yerr=Ze,
                         linestyle=Linetype[0], marker=Markers[1],
                         color=Linecolor[2], linewidth=Linewidth[0],
                         markersize=Markersize[0], label="Z, observed" )
 
  
    
                
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

        ax.grid(True)
        ax.invert_yaxis()
        ax.grid("major", "both", linestyle=":", lw=0.3)
      
    if Legend:
        ax.legend(fontsize=Fontsizes[0], loc="best")

    if ThisAxis==None:
        for F in PlotFormat:
            matplotlib.pyplot.savefig(PlotFile+F)      
        # matplotlib.pyplot.show()
        # matplotlib.pyplot.clf()
            
    return ax


def plot_model(
        ThisAxis = None, 
        PlotFile = "",        
        PlotFormat = ["png",],
        PlotTitle = "",
        PlotSize = [8.],
        System  = "aem05",
        Model = [],
        Depth = [],
        Fillcolor=["0.8", "0.4"],
        Alphas = [0.3 , 0.6],
        Labels=[],
        Linecolor=["k", "r", "g", "b", "c", "m"],
        Linetype=["-", ":", ";"],
        Linewidth=[1., 1.5, 2.],
        Markers = ["v"],
        Markersize =[4],
        Fontsizes=[10,10,12],
        XLimits= [],
        YLimits=[],
        PlotStrng="",
        Invalid=1.e30,
        Legend=True):

    """
    
    """
    cm = 1/2.54  # centimeters to inches
    
    ax = ThisAxis
    
    if numpy.size(Model)==0:
       error("No parameter set given!! Exit.")
       # nopar=True

    if ThisAxis==None:
        nplots = 1
        fig, ax = matplotlib.pyplot.subplots(1,
                                          figsize=(PlotSize[0]*cm, nplots*PlotSize[0]*cm),
                                          gridspec_kw={"height_ratios": [1]})
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])
        Legend = True
        

    ax.step(Model, Depth,
                      where="pre", color=Linecolor[0],
                      linewidth=Linewidth[0],  linestyle=Linetype[0],           
                      label=Labels[0])

        # ax.step(medval, Depth,step="pre",
        #     linewidth=Linewidth[0]+2, color= Linecolor[1],
        #     label="medval")

    ax.set_xscale("log")
    ax.set_xlim(XLimits)
    ax.set_xlabel("resistivity ($\Omega$m)",fontsize=Fontsizes[0])
    ax.set_ylim(YLimits)
    ax.set_ylabel("depth (m)",fontsize=Fontsizes[0])
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position("both")
    ax.tick_params(labelsize=Fontsizes[1])

    ax.grid(True)
    ax.invert_yaxis()
    ax.grid("major", "both", linestyle=":", lw=0.3)
    
    if Legend: 
        ax.legend(fontsize=Fontsizes[0], loc="best")
        
    if ThisAxis==None:        
        for F in PlotFormat:
            matplotlib.pyplot.savefig(PlotFile+F)
    # matplotlib.pyplot.show()
    # matplotlib.pyplot.clf()
        
    return ax
