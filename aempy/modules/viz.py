# # -*- coding: utf-8 -*-
# """
# Created on Sun Dec 27 17:23:34 2020

# @author: vrath
# """

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings

import math
import numpy
import matplotlib
import matplotlib.pyplot
import matplotlib.colors
import matplotlib.cm
import mpl_toolkits.axes_grid1
import matplotlib.ticker
import cycler

import numpy
import numpy.ma
import simplekml

import aesys
import inverse
import util

def plot_depth_prof(
        ThisAxis = None, 
        PlotFile = None,
        PlotFormat = ["png",],                    
        PlotTitle = None,                    
        FigSize = [8.5*0.3937, 8.5*0.3937],
        Depth = [],        
        DLimits = [],        
        DLabel = " Depth (m)",
        Params = [],           
        Partyp = "",
        PLabel = "",  
        PLimits = [],
        Shade = [ 0.25 ],
        XScale = "log",
        PlotType = "steps",
        Legend = [],
        Linecolor = ["r", "g", "b", "m", "y"],
        Linetypes =  ["-","-","-","-","-"],
        Linewidth =  [1, 1, 1, 1,1,],
        Marker = ["v"],
        Markersize =[4],
        Fillcolor = [[0.7, 0.7, 0.7]],
        Logplot = True,
        Fontsizes =[10, 10, 12],
        PlotStrng="",
        StrngPos=[0.05,0.05],
        Save = True, 
        Invalid=1.e30):
    """
    General plot of (multiple) depth profiiles 
    
    Parameters
    ----------
    
    Depth :  np.array
        DESCRIPTION. The default is [].
    Params : np.array
        DESCRIPTION. The default is [].
    PLabels : list of strings, optional
        DESCRIPTION. The default is []. 

    XScale: string, optional
        "linear", "log", "symlog", "asinh"
        Last two need further parameters, e.g
        ax.set_yscale("asinh", linear_width=a0)
        x.set_yscale("symlog", linthresh=2,)
    Ptype : string, optional
        Proy type. The default is "steps".
    
    ALabels: string, optional
        Axis Labels for Params and Depth. The default is "", and "(m}.
    PLimits,  DLimits : lists, optional
        Limits for Params and Depth The default is [].    
    Errors : TYPE, optional
        DESCRIPTION. The default is [].
    PlotFile : string, optional
        Plot file name without extension
        The default is None.
    PlotTitle : string, optional
        Plot title
        he default is None.
    PlotFormat : string, optional
        List of output formats. The default is ["png",].    
     Linecolor : TYPE, optional
        DESCRIPTION. The default is ["y", "r", "g", "b", "m"].
    Linetypes : TYPE, optional
        DESCRIPTION. The default is "".
    Fontsizes : TYPE, optional
        DESCRIPTION. The default is [10, 10, 12].
    PlotStrng : string, optional
        Annotation. The default is "".
    StrngPos : TYPE, optional
        Annotation proition. The default is [0.05,0.05].

    Returns
    -------
    ax
       
    Created May 1, 2023
    @author: vrath

    """

    cm = 1/2.54  # centimeters in inches
    
    if ThisAxis==None:
        fig, ax =  matplotlib.pyplot.subplots(1, 1, figsize=(FigSize))
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])
    
    

    for iparset in range(len(Params)):
        
        P = Params[iparset]
        D = Depth[iparset]
        np = numpy.shape(P)[0]
        nd = numpy.shape(D)[0]

        df = D[-1] + 3*numpy.abs(D[-1]-D[-2])
 
 
        if Partyp=="":
    
            if "steps" in PlotType.lower():
                d = D
                for pp in numpy.arange(np):
                    print("PPPP ",P[pp])
                    p = numpy.append(P[pp],P[pp][-1])
                    ax.step(p , d , 
                         where='pre',
                         c=Linecolor[pp],
                         ls=Linetypes[pp], lw=Linewidth[pp])        
            else: 
                d = D
                for pp in numpy.arange(np):
                    p = P[pp]                       
                    p = numpy.append(P[pp],P[pp][-1])
                    ax.plot(p , d,
                            c=Linecolor[pp],                                  
                            ls=Linetypes[pp], lw=Linewidth[pp])
                    
                    
        if "sens" in Partyp.lower():
            if "steps" in PlotType.lower():
                d = D[:-1]
                for pp in numpy.arange(np):
                    p = P[pp]
                    print(numpy.shape(d),numpy.shape(p))
                    ax.step(p , d , 
                         where='pre',
                         c=Linecolor[pp],
                         ls=Linetypes[pp], lw=Linewidth[pp])        
            else: 
                d = D[:-1]
                for pp in numpy.arange(np):
                    p = P[pp]
                    ax.plot(p , d,
                            c=Linecolor[pp],                                  
                            ls=Linetypes[pp], lw=Linewidth[pp])
                    

        if "model" in Partyp.lower():
            
            d = numpy.append(D, df)
            
            if np==3: 

             
                p = numpy.append(P[0],P[0][-1])
                ep = numpy.append(P[1],P[1][-1])
                em = numpy.append(P[2],P[2][-1])
                
                if "fill" in PlotType.lower():    
                    ax.fill_betweenx(d, em, ep, 
                                step='post', 
                                color=Fillcolor[0],                                  
                                ls=Linetypes[0], lw=Linewidth[0],
                                alpha=Shade)
                                    
                ax.step(p , d , 
                         where='pre',
                         c=Linecolor[0],
                         ls=Linetypes[0], lw=Linewidth[0])        
                ax.plot(p[-1] , d[-1],
                        c=Linecolor[0],  ls=Linetypes[0], lw=0,
                        marker=Marker[0], markersize=Markersize[0])     
                ax.step(em , d , 
                         where='pre',
                         c=Linecolor[1],
                         ls=Linetypes[1], lw=Linewidth[1])                         
                ax.step(ep , d , 
                         where='pre',
                         c=Linecolor[1],
                         ls=Linetypes[1], lw=Linewidth[1]) 

            else:
   
                if "step" in PlotType.lower():
                    for pp in numpy.arange(np):
                        p = P[pp]
                        ax.step(p , d, 
                                where='pre',
                                color=Linecolor[pp],
                                    ls=Linetypes[pp],lw=Linewidth[0])
                else: 
                    for pp in numpy.arange(np):
                        p = P[pp]   
                        ax.plot(p , d,
                                c=Linecolor[pp],                                  
                                ls=Linetypes[pp], lw=Linewidth[0])
        
        ax.set_xlabel(PLabel, fontsize=Fontsizes[1])
        ax.set_ylabel(DLabel, fontsize=Fontsizes[1])
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("both")
        ax.tick_params(labelsize=Fontsizes[0])    
        
        if PLimits != []:
            ax.set_xlim(PLimits)
        if DLimits != []:
            ax.set_ylim(DLimits)
        
        if "lin" not in XScale:
            ax.set_xscale(XScale)
    
            
            
        ax.legend(Legend, fontsize=Fontsizes[1]-2, loc="best", ncol=1)
        
        if PLimits != []:
            ax.set_xlim(PLimits)
        if DLimits != []:
            ax.set_ylim(DLimits)
    
        ax.invert_yaxis()
    
        ax.grid("major", "both", linestyle=":", lw=0.3)
        ax.text(StrngPos[0], StrngPos[1],
                 PlotStrng, fontsize=Fontsizes[1]-1,transform=ax.transAxes,
                 bbox=dict(facecolor="white", alpha=0.5) )
        
        if ThisAxis==None:
            for F in PlotFormat:
                 matplotlib.pyplot.savefig(PlotFile+F)
        
            matplotlib.pyplot.show()
            matplotlib.pyplot.clf()
        
        return ax

def plot_matrix(
        ThisAxis = None, 
        PlotFile = None,
        PlotTitle = None,
        PlotFormat = ["png",],
        FigSize = [8.5*0.3937, 8.5*0.3937],
        Matrix = [],
        TickStr="",
        AxLabels = ["layer #", "layer #"],
        AxTicks = [[], []], 
        AxTickLabels = [[], []], 
        ColorMap="viridis",
        Fontsizes=[10,10,12],
        Unit = "",        
        PlotStrng="",
        StrngPos=[0.05,0.05],
        Aspect = "auto",
        Save = True, 
        Invalid=1.e30):
    """
    Plots jacobians, covariance and resolution matrices.
    
  
    Parameters
    ----------
    PlotFile : TYPE, optional
        DESCRIPTION. The default is None.
    PlotTitle : TYPE, optional
        DESCRIPTION. The default is None.
    PlotFormat : TYPE, optional
        DESCRIPTION. The default is ["png",].


    Returns
    -------
    ax
    
    
    Created April 30, 2023
    @author: vrath

    """
    cm = 1/2.54  # centimeters in inches
    nn = numpy.shape(Matrix)
    np = nn[0]
    if Matrix.ndim==1:
        np =math.isqrt(nn[0])
        Matrix = Matrix.reshape((np,np))
        
    if ThisAxis==None:
        fig, ax =  matplotlib.pyplot.subplots(1, 1, figsize=(FigSize))
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])
    

    im = ax.imshow(Matrix, cmap=ColorMap, origin="upper")
    
    xticks = AxTicks[0]
    xlabels = AxTickLabels[0]
    # print(xticks)
    # print(xlabels)
    ax.set_xticks(xticks, xlabels) #, minor=False)
    ax.set_xlabel(AxLabels[0], fontsize=Fontsizes[1])
    ax.xaxis.set_ticks_position("top") 
    ax.xaxis.set_label_position("top") 
    
    yticks = AxTicks[1]
    ylabels = AxTickLabels[1]    
    # print(yticks)
    # print(ylabels)
    ax.set_yticks(yticks, ylabels) #, minor=False)
    ax.set_ylabel(AxLabels[1], fontsize=Fontsizes[1])
                  
    if Aspect == "equal":
        ax.set_aspect("equal","box")
    else:
        ax.set_aspect(Aspect)
    
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cb = matplotlib.pyplot.colorbar(im, cax=cax)
    cb.ax.set_title(Unit)
    
    if PlotStrng != "":
        props = dict(facecolor="white", alpha=0.9) # boxstyle="round"
        ax.text(StrngPos[0], StrngPos[1], PlotStrng, 
                transform=ax.transAxes, 
                fontsize=Fontsizes[1],
                verticalalignment="center", bbox=props)
 

    if ThisAxis==None:
        for F in PlotFormat:
             matplotlib.pyplot.savefig(PlotFile+F)
    
        matplotlib.pyplot.show()
        matplotlib.pyplot.clf()
    
    return ax


def plot_data_genesis(
        PlotFile = None,
        PlotTitle = None,
        PlotFormat = ["png",],
        DPI = 400,
        Data=[],
        Errors=[],
        Position=None,
        Labels=[],
        Linecolor=[],
        Linetypes=["-", ":",],
        Linewidth=[1],
        Markers = ["v"],
        Markersize =[4],
        Fontsizes=[10,10,12],
        YLimits= [],
        PlotStrng="",
        StrngPos=[0.05,0.05],
        LogPlot = True,
        SymLog=False,
        TimeLog =True,
        DataTrans = 2, 
        Save = True,
        Invalid=1.e30):


    cm = 1/2.54  # centimeters in inches

    AEM_system = "genesis"
    _,_,_,_,pars = aesys.get_system_params(System=AEM_system)
    wincent = pars[0]


    nn = numpy.shape(Data)

    

    Prefix = ""
    if DataTrans == 1:
        Prefix = "log10 "
        LogPlot = False
    if DataTrans == 2:
        LogPlot = False
        Prefix = "arcsinh "

    if Save:
        fig, (ax1, ax2) =  matplotlib.pyplot.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])

    if len(Errors) == 0:
        for ll in numpy.arange(0, nn[0]):
            print(ll, Linetypes[ll], Markers[ll], Linecolor[ll])
            ax1.plot(wincent, Data[ll, 0:11].T,
                  linestyle=Linetypes[ll], marker=Markers[ll], color=Linecolor[ll],
                  linewidth=Linewidth,markersize=Markersize[ll],label = Labels[ll])
    else:
        for ll in numpy.arange(0, nn[0]):
            ax1.errorbar(wincent, Data[ll, 0:11].T,yerr=Errors[ll, 0:11].T,
                         linestyle=Linetypes[ll], marker=Markers[ll],
                         color=Linecolor[ll-1],linewidth=Linewidth,
                         markersize=Markersize[ll], label = Labels[ll])
    ax1.set_xlabel("Time (ms)", fontsize=Fontsizes[0])
    ax1.set_ylabel(Prefix+" B$_{inline}$ (ppm)", fontsize=Fontsizes[0])
    ax1.xaxis.set_tick_params(labelsize=Fontsizes[1])
    ax1.yaxis.set_tick_params(labelsize=Fontsizes[1])
    ax1.grid("major", "both", linestyle=":", lw=0.3)
    ax1.legend(Labels[:], fontsize=Fontsizes[1]-2, loc="best", ncol=1)
    
    if TimeLog:
        ax1.set_xscale("log")
    else:
        ax1.set_xscale("linear")


    if LogPlot:
        if SymLog:
            ax1.set_yscale("symlog")
        else:
            ax1.set_yscale("log", nonpositive="clip")
    else:
        ax1.set_yscale("linear")

    if len(YLimits) !=0:
        ax1.set_ylim(YLimits)


    if len(Errors) == 0:
        for ll in numpy.arange(0, nn[0]-1):            
            print(ll, Linetypes[ll], Markers[ll], Linecolor[ll])
            ax2.plot(wincent, Data[ll, 11:22].T,
                 linestyle=Linetypes[ll], marker=Markers[ll], color=Linecolor[ll],
                 linewidth=Linewidth ,markersize=Markersize[ll])
    else:
        for ll in numpy.arange(0, nn[0]-1):
            ax2.errorbar(wincent, Data[ll, 11:22].T,yerr=Errors[ll, 11:22].T,
                         linestyle=Linetypes[ll], marker=Markers[ll],
                         color=Linecolor[ll],linewidth=Linewidth,
                         markersize=Markersize[ll])


    ax2.set_xlabel("Time (ms)", fontsize=Fontsizes[0])
    ax2.set_ylabel(Prefix+"B$_{vert}$ (ppm)", fontsize=Fontsizes[0])
    ax2.xaxis.set_tick_params(labelsize=Fontsizes[1])
    ax2.yaxis.set_tick_params(labelsize=Fontsizes[1])
    ax2.grid("major", "both", linestyle=":", lw=0.3)

    if TimeLog:
        ax2.set_xscale("log")
    else:
        ax2.set_xscale("linear")

    if LogPlot:
        if SymLog:
            ax2.set_yscale("symlog")
        else:
            ax2.set_yscale("log", nonpositive="clip")
    else:
        ax2.set_yscale("linear")

    ax2.text(StrngPos[0], StrngPos[1],
             PlotStrng, fontsize=Fontsizes[1]-1,transform=ax2.transAxes,
             bbox=dict(facecolor="white", alpha=0.5))



    if Save:
        for F in PlotFormat:
             matplotlib.pyplot.savefig(PlotFile+F, dpi=DPI)
             

        matplotlib.pyplot.show()
        matplotlib.pyplot.clf()
    
    return fig,  ax1, ax2

def plot_data_aem05(
        PlotFile = None,
        PlotTitle = None,
        PlotFormat = ["png",],
        DPI = 400,
        DataTrans=None,
        Data=[],
        Errors=[],
        Position=None,
        Labels=[],
        Linecolor=None,
        Linetypes=["-", ":",],
        Linewidth=[1],
        Markers=['o'],
        Markersize = [4],
        Fontsizes=[10,10,12],
        YLimits= [],
        PlotStrng="",
        StrngPos=[0.05,0.05],
        LogPlot = False,
        SymLog=False, 
        Invalid=1.e30):

    cm = 1/2.54  # centimeters in inches

    AEM_system = "aem05"
    FwdCall,nD,_,_,freq = aesys.get_system_params(System=AEM_system)

    nn = numpy.shape(Data)

    Prefix=""
    if DataTrans == 1:
        Prefix = "log10 "
        LogPlot = False
        Data, Errors, _=inverse.transform_data(d_vec=Data, e_vec=Errors, d_trn=DataTrans)



    if Save:
        fig, (ax1, ax2) =  matplotlib.pyplot.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
        fig.suptitle(PlotTitle, fontsize=Fontsizes[2])

    if Errors.size == 0:
        for ll in numpy.arange(0, nn[0]-1):
            ax1.plot(freq, Data[ll, 0:4].T,
                  linestyle=Linetypes[ll], marker="o", color=Linecolor[ll],
                  linewidth=Linewidth, markersize=Markersize[ll],)
    else:
        for ll in numpy.arange(0, nn[0]-1):
            ax1.errorbar(freq, Data[ll, 0:4].T,yerr=Errors[ll, 0:4].T,
                         linestyle=Linetypes[ll], marker="o",
                         color=Linecolor[ll],linewidth=Linewidth,
                         markersize=Markersize[ll])

    ax1.set_xlabel(" Frequency (Hz)", fontsize=Fontsizes[0])
    ax1.set_ylabel(Prefix +" In-Phase (ppm)", fontsize=Fontsizes[0])
    ax1.xaxis.set_tick_params(labelsize=Fontsizes[1])
    ax1.yaxis.set_tick_params(labelsize=Fontsizes[1])
    ax1.grid("major", "both", linestyle=":", lw=0.3)
    ax1.legend(Labels, fontsize=Fontsizes[1]-2, loc="best", ncol=1)
    ax1.set_xscale("log", nonpositive="clip")
    ax1.set_yscale("linear")
    if len(YLimits) !=0:
        ax1.set_ylim(YLimits)


    if Errors.size == 0:
        for ll in numpy.arange(0, nn[0]-1):
            ax2.plot(freq, Data[ll, 4:8].T,
                  linestyle=Linetypes[ll], marker="o", color=Linecolor[ll],
                  linewidth=Linewidth, markersize=Markersize[ll])

    else:
        for ll in numpy.arange(0, nn[0]-1):
            ax2.errorbar(freq, Data[ll, 4:8].T,yerr=Errors[ll, 4:8].T,
                         linestyle=Linetypes[ll], marker="o",
                         color=Linecolor[ll], linewidth=Linewidth,
                         markersize=Markersize[ll])

    ax2.set_xlabel(" Frequency (Hz)", fontsize=Fontsizes[0])
    ax2.set_ylabel(Prefix +" Quadrature (ppm)", fontsize=Fontsizes[0])

    ax2.xaxis.set_tick_params(labelsize=Fontsizes[1])
    ax2.yaxis.set_tick_params(labelsize=Fontsizes[1])
    ax2.grid("major", "both", linestyle=":", lw=0.3)
    ax2.set_xscale("log", nonpositive="clip")
    ax2.set_yscale("linear")
    ax2.text(StrngPos[0], StrngPos[1],
             PlotStrng, fontsize=Fontsizes[1]-1,transform=ax2.transAxes,
             bbox=dict(facecolor="white", alpha=0.5))



    if Save:
        for F in PlotFormat:
             matplotlib.pyplot.savefig(PlotFile+F, dpi=DPI)

        matplotlib.pyplot.show()
        matplotlib.pyplot.clf()
    
    return fig, ax1, ax2




def plot_flightline_aem05(
        PlotName = None,
        PlotDir="./",
        PlotFmt = [".png"],
        PlotSize = [6., 18.],
        IncludePlots = ["qdata", "idata",],
        DataObs=None,
        DataCal=None,
        DataTrans=None,
        DLimits = [],
        QLimits = [],
        ILimits = [],
        HLimits = [],        
        ProfLabel = "profile distsnce (m)",
        PosDegrees=False,
        EPSG=32629,
        Linecolor=None,
        Linetypess=["-", ":",],
        Linewidth=[1,2,],
        Fontsizes=[12,10,16],
        Grey=[0.7],
        Logparams=[False,False,10.],
        PlotStrng="",
        PlotAltDiff=False,
        PlotPLM = False,
        PLimits = [0., 5.], 
        Save=True,
        Invalid=1.e30):

    cm = 1./2.54  # centimeters to inches
    ierr = 0

    Fsize = Fontsizes[0]
    Lsize = Fontsizes[1]
    Tsize = Fontsizes[2]

    Greyval = Grey[0]

    LogPlot= Logparams[0]
    LogSym = Logparams[1]
    LinThresh = Logparams[2]




    fline = DataObs[:, 0]
    site_x = DataObs[:, 1]
    site_y = DataObs[:, 2]
    site_gps = DataObs[:, 3]
    site_alt = DataObs[:, 4]
    site_dem = DataObs[:, 5]
    data_obs = DataObs[:, 6:14]
    site_plm = DataObs[:, 14]

    # print("alt")
    # print(site_alt)
    # print("plm")
    # print(site_plm)

    beg_pos = [site_x[0],  site_y[0]]
    end_pos = [site_x[-1], site_y[-1]]
    if PosDegrees:
        beg_pos = util.project_utm_to_latlon(beg_pos[0], beg_pos[1])
        end_pos = util.project_utm_to_latlon(end_pos[0], end_pos[1])

    # ####### Plot data along flight line #########
    # Name = str(fline[0])
    if "dist" in ProfLabel.lower():
        sx = site_x-site_x[0]
        sy = site_y-site_y[0]
        D = numpy.sqrt(sx * sx + sy * sy)
    else:
        D = numpy.arange(numpy.shape(site_x)[0])
        
    D_min, D_max = numpy.amin(D), numpy.amax(D)

    IData = data_obs[:, 0:4]
    QData = data_obs[:, 4:8]

    I_min, I_max = numpy.nanmin(IData), numpy.nanmax(IData)
    Q_min, Q_max = numpy.nanmin(QData), numpy.nanmax(QData)

    print("Dmin,max = "+str(D_min)+",  "+str(D_max))
    print("Imin,max = "+str(I_min)+",  "+str(I_max))
    print("Qmin,max = "+str(Q_min)+",  "+str(Q_max))

    nplots = numpy.size(IncludePlots)

    fig =  matplotlib.pyplot.figure(figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm))

    ii = 0

    for splot in IncludePlots:

        ii = ii+1

        print("Plotting "+splot+" at position "+str(nplots)+str(1)+str(ii))
        if ii == 1:
            ax =fig.add_subplot(nplots, 1, ii)
            first_ax= ax
        else:
            ax =fig.add_subplot(nplots, 1, ii, sharex = first_ax)


        if "alt"in splot.lower():


            if PlotAltDiff:
                ax.plot( D[:],
                           site_alt[:]-(site_gps[:] - site_dem[:]),
                    "r-", linewidth=Linewidth[0])
                ax.set_title("Altitude Differences", fontsize=Fsize, y=1.0, pad=-(Fsize+2))

            else:
                ax.plot(
                    D[:-1], site_alt[:-1],
                    "r", linewidth=Linewidth[0])
                ax.plot(
                    D[:-1], site_gps[:-1]-site_dem[:-1],
                    "g:", linewidth=Linewidth[0])
                ax.legend([" Radar alt", "GPS alt"], fontsize=Lsize, loc="upper left")
                if HLimits:
                    ax.set_ylim(HLimits)

            ax.set_ylabel(" altitude (m)", fontsize=Fsize)
            if ii == nplots:
                ax.set_xlabel(ProfLabel, fontsize=Lsize)
            else:
                ax.tick_params(labelbottom=False)

            ax.tick_params(labelsize=Lsize)
            ax.grid(True)


            if PlotPLM:
                ax2=ax.twinx()   # make a plot with different y-axis using second axis object
                ax2.plot(D[:],site_plm[:],color="blue",linewidth=Linewidth[0])
                ax2.set_ylabel("plm (nT)", fontsize=Lsize)
                ax2.legend([" powerline monitor"], fontsize=Lsize, loc="upper right")
                if PLimits:
                    ax2.set_ylim(PLimits)

        if "qdata"in splot.lower():

            ax.plot(
                D[:],        QData[:, 0],        "r",
                D[:],        QData[:, 1],        "g",
                D[:],        QData[:, 2],        "b",
                D[:],        QData[:, 3],        "y",
                linewidth=Linewidth[0])
            ax.set_title("Quadrature", fontsize=Fsize, y=1.0, pad=-(Fsize+2))
            ax.set_ylim([Q_min, Q_max])
            if numpy.size(QLimits)>0:
                ax.set_ylim(QLimits)

            ax.set_ylabel("(ppm)", fontsize=Lsize)
            if ii == nplots:
                ax.set_xlabel(ProfLabel, fontsize=Lsize)
            else:
                ax.tick_params(labelbottom=False)

            ax.tick_params(labelsize=Lsize)
            ax.grid(True)


            legend =ax.legend(
                [" 0.9 kHz", "3 kHz", "12 kHz", "24.5 kHz"],
                fontsize=Lsize-2, loc="best", ncol=2)
            legend.set_title("Frequency", prop={"size":Lsize})

            if LogPlot:
                if LogSym:
                    ax.set_yscale("symlog", linthresh=LinThresh)
                else:
                    ax.set_yscale("log")
            else:
                ax.set_yscale("linear")

        if "idata"in splot.lower():


            ax.plot(
                D[:],        IData[:, 0],        "r",
                D[:],        IData[:, 1],        "g",
                D[:],        IData[:, 2],        "b",
                D[:],        IData[:, 3],        "y",
                linewidth=Linewidth[0])

            ax.set_title("In-Phase", fontsize=Fsize, y=1.0, pad=-(Fsize+2))
            ax.set_ylim([I_min, I_max])
            if numpy.size(ILimits)>0:
                ax.set_ylim(ILimits)
            ax.set_ylabel("(ppm)", fontsize=Fsize)

            if ii == nplots:
                ax.set_xlabel(ProfLabel, fontsize=Lsize)
            else:
                ax.tick_params(labelbottom=False)



            ax.grid(True)

            legend =ax.legend(
                [" 0.9 kHz", "3 kHz", "12 kHz", "24.5 kHz"],
                fontsize=Lsize-2, loc="best", ncol=2)
            legend.set_title("Frequency", prop={"size":Lsize})

            if LogPlot:
                if LogSym:
                    ax.set_yscale("symlog", linthresh=LinThresh)
                else:
                    ax.set_yscale("log")
            else:
                ax.set_yscale("linear")

            if PosDegrees:
                beg_strng  = "Lat: "+str(beg_pos[0])+"\nLon: "+str(beg_pos[1])
                end_strng  = "Lat: "+str(end_pos[0])+"\nLon: "+str(end_pos[1])
            else:
                beg_strng  = "start:"+"\nX: "+str(round(beg_pos[0],2))+"\nY: "+str(round(beg_pos[1],2))+" (EPSG="+str(EPSG)+")"
                end_strng  = "end:  "+"\nX: "+str(round(end_pos[0],2))+"\nY: "+str(round(end_pos[1],2))

            fig.text(0.10, 0.004, beg_strng, fontsize = Fsize-2, ha="left")
            fig.text(0.85, 0.004, end_strng, fontsize = Fsize-2, ha="left")

    fig.suptitle("Flightline: "+PlotName+PlotStrng,
                 #x=0.5, y=0.9, fontsize=Tsize+1, ha="center") #, fontweight="bold")
                 x=0.5, y=0.95, fontsize=Tsize+0, ha="center", va ="top") #, fontweight="bold")

    if Save:
        for F in PlotFmt:
             matplotlib.pyplot.savefig(PlotDir+PlotName+F)

        if matplotlib.get_backend()!="cairo":
            matplotlib.pyplot.show()
            
        matplotlib.pyplot.clf()

    return fig, ax



def plot_flightline_genesis(
        PlotName = None,
        PlotDir="./",
        PlotFmt = [".png"],
        PlotSize = [18., 8.],
        IncludePlots = ["xdata", "zdata",],
        DataObs=None,
        DataCal=None,
        XLimits =[],
        ZLimits =[],
        DLimits = [],
        HLimits = [],
        ProfLabel = "profile distsnce (m)",
        Linecolor=None,
        Linetypess=["-", ":",],
        Linewidth=[1,2],
        Fontsizes=[12,10,12],
        Grey=[0.7],
        DataTrans=1,
        PosDegrees=False,
        EPSG=32629,
        Logparams=[True,True,0,1],
        PlotStrng="",
        PlotAltDiff=False, 
        Save=True, 
        Invalid=1.e30):

    cm = 1./2.54  # centimeters to inches
    ierr = 0

    Fsize = Fontsizes[0]
    Lsize = Fontsizes[1]
    Tsize = Fontsizes[2]
    
    Greyval = Grey[0]

    LogPlot= Logparams[0]
    LogSym = Logparams[1]
    LinThresh = Logparams[2]


    # fline = DataObs[:, 0]

    site_x = DataObs[:, 1]
    site_y = DataObs[:, 2]

    site_gps = DataObs[:, 3]
    site_alt = DataObs[:, 4]
    site_dem = DataObs[:, 5]
    XData = DataObs[:, 6:17]
    ZData = DataObs[:, 17:28]

    beg_pos = [site_x[0],  site_y[0]]
    end_pos = [site_x[-1], site_y[-1]]
    if PosDegrees:
        beg_pos = util.project_utm_to_latlon(beg_pos[0], beg_pos[1], utm_zone=EPSG)
        end_pos = util.project_utm_to_latlon(end_pos[0], end_pos[1], utm_zone=EPSG)

    dunit = "(ppm)"
    if DataTrans==2:
        
        XData = inverse.transform_data() numpy.arcsinh(XData)
        ZData = numpy.arcsinh(ZData)
        LogPlot = False
        dunit = "(-)"

    # Name = str(fline[0])
    if "dist" in ProfLabel.lower():
        sx = site_x-site_x[0]
        sy = site_y-site_y[0]
        D = numpy.sqrt(sx * sx + sy * sy)
    else:
        D = numpy.arange(numpy.shape(site_x)[0])

    D_min, D_max = numpy.amin(D), numpy.amax(D)

    X_min, X_max = numpy.amin(XData), numpy.amax(XData)
    Z_min, Z_max = numpy.amin(ZData), numpy.amax(ZData)

    print("Dmin,max = "+str(D_min)+",  "+str(D_max))
    print("Xmin,max = "+str(X_min)+",  "+str(X_max))
    print("Zmin,max = "+str(Z_min)+",  "+str(Z_max))




    if LogPlot:
        if not LogSym:
            X_minpos = numpy.amin(XData[XData>0.])
            Z_minpos = numpy.amin(ZData[ZData>0.])
            if LinThresh == []:
                LinThresh = numpy.amin([Z_minpos])
            XData[XData<=numpy.amin([X_minpos, Z_minpos])] = numpy.nan
            ZData[ZData<=numpy.amin([X_minpos, Z_minpos])] = numpy.nan
            print("Xminpos,Zminpos  = "+str(X_minpos)+",  "+str(Z_minpos))
            print("LinThresh = "+str( LinThresh))

    nplots = numpy.size(IncludePlots)

    fig =  matplotlib.pyplot.figure(figsize=(PlotSize[0]*cm, nplots*PlotSize[1]*cm))

    ii = 0

    for splot in IncludePlots:

        ii = ii+1

        if ii == 1:
            ax =fig.add_subplot(nplots, 1, ii)
            first_ax= ax
        else:
            ax =fig.add_subplot(nplots, 1, ii, sharex = first_ax)

        print("Plotting "+splot+" at position "+str(nplots)+str(1)+str(ii))


        if "alt"in splot.lower():
            if PlotAltDiff:
                ax.plot(
                    D[:],
                    site_alt[:]-(site_gps[:] - site_dem[:]),"r-")
                ax.set_title("Altitude Differences", fontsize=Fsize, y=1.0, pad=-(Fsize+2))

            else:
                ax.plot(
                    D[:], site_alt[:], "r",
                    D[:], site_gps[:] - site_dem[:], "g:")
                ax.set_title("Altitudes", fontsize=Fsize, y=1.0, pad=-(Fsize+2))
                ax.legend([" Radar alt", "GPS alt"], fontsize=Lsize, loc="best")

                if HLimits:
                    ax.set_ylim(HLimits)

            ax.set_title("Altitude", fontsize=Fsize, y=1.0, pad=-(Fsize+2))
            ax.set_ylabel("(m)", fontsize=Fsize)
            if ii == nplots:
                ax.set_xlabel(ProfLabel, fontsize=Lsize)
            else:
                ax.tick_params(labelbottom=False)

            ax.tick_params(labelsize=Lsize)
            ax.grid(True)


        if "xdata"in splot.lower():

            title="In-Line"
            if "asinh"in DataTrans.lower():
                title=title+" / arcsinh(X)"
                dunit = "(-)"

            ax.plot(D[:], XData[:,:], linewidth=Linewidth[0])
            ax.set_title(title, fontsize=Fsize, y=1.0, pad=-(Fsize+2))
            if numpy.size(XLimits)>0:
                ax.set_ylim(XLimits)
            ax.set_ylabel(dunit, fontsize=Fsize)

            if ii == nplots:
                ax.set_xlabel(ProfLabel, fontsize=Lsize)
            else:
                ax.tick_params(labelbottom=False)

            ax.grid(True)

            ax.tick_params(labelsize=Lsize)
            leg1 = ax.legend(
                [r"0.009 ms", r"0.026 ms", r"0.052 ms", r"0.095 ms",
                  r"0.156 ms", r"0.243 ms", r"0.365 ms", r"0.547 ms",
                  r"0.833 ms", r"1.259 ms", r"1.858 ms"],
                fontsize=Lsize-3, loc="best",ncol=3)
            leg1.set_title("Window (center)", prop={"size":Lsize-2})

            if LogPlot:
                if LogSym:
                    ax.set_yscale("symlog", linthresh=LinThresh)
                else:
                    ax.set_yscale("log")
            else:
                ax.set_yscale("linear")


        if "zdata"in splot.lower():

            title="Vertical"
            if "asinh"in DataTrans.lower():
                title=title+" / arcsinh(Z)"
                dunit = "(-)"

            ax.plot(D[:], ZData[:,:], linewidth=Linewidth[0])
            ax.set_title(title, fontsize=Fsize, y=1.0, pad=-(Fsize+2))
            if numpy.size(ZLimits)>0:
                ax.set_ylim(ZLimits)
            ax.set_ylabel(dunit, fontsize=Fsize)

            if ii == nplots:
                ax.set_xlabel(ProfLabel, fontsize=Lsize)
            else:
                ax.tick_params(labelbottom=False)

            ax.grid(True)
            ax.tick_params(labelsize=Lsize)
            leg2 = ax.legend(
                [r"0.009 ms", r"0.026 ms", r"0.052 ms", r"0.095 ms",
                  r"0.156 ms", r"0.243 ms", r"0.365 ms", r"0.547 ms",
                  r"0.833 ms", r"1.259 ms", r"1.858 ms"],
                fontsize=Lsize-3, loc="best",ncol=3,title="Window (center)")
            leg2.set_title("Window (center)", prop={"size":Lsize-2})
            if LogPlot:
                if LogSym:
                    ax.set_yscale("symlog", linthresh=LinThresh)
                else:
                    ax.set_yscale("log")
            else:
                ax.set_yscale("linear")

            if PosDegrees:
                beg_strng  = "Lat: "+str(beg_pos[0])+"\nLon: "+str(beg_pos[1])
                end_strng  = "Lat: "+str(end_pos[0])+"\nLon: "+str(end_pos[1])
            else:
                beg_strng  = "start:"+"\nX: "+str(round(beg_pos[0],2))+"\nY: "+str(round(beg_pos[1],2))+" (EPSG="+str(EPSG)+")"
                end_strng  = "end:  "+"\nX: "+str(round(end_pos[0],2))+"\nY: "+str(round(end_pos[1],2))

    fig.text(0.10, 0.004, beg_strng, fontsize = Fsize-2, ha="left")
    fig.text(0.85, 0.004, end_strng, fontsize = Fsize-2, ha="left")

    fig.set_tight_layout
    fig.suptitle("Flightline: "+PlotName+PlotStrng,
                 x=0.5, y=0.95, fontsize=Tsize+0, ha="center", va ="top") #, fontweight="bold")


    if Save:
        for F in PlotFmt:
             matplotlib.pyplot.savefig(PlotDir+PlotName+F)

        if matplotlib.get_backend()!="cairo":
            matplotlib.pyplot.show()
        matplotlib.pyplot.clf()
        
    return fig, ax 
    




def make_pdf_catalog(PDFList= None, FileName="Catalog.pdf"):
    """
    Make pdf catalog from list of pdf plots

    Parameters
    ----------
    PDFList : List of strings
        List of Filenames
    Filename : string
        Catalog file


    Returns
    -------
    None.

    """
    # error("not in 3.9! Exit")

    import fitz

    catalog = fitz.open()

    for pdf in PDFList:
        with fitz.open(pdf) as mfile:
            catalog.insert_pdf(mfile)

    catalog.save(FileName, garbage=4, clean = True, deflate=True)
    catalog.close()

    print("\n"+str(numpy.size(PDFList))+" files collected to "+FileName)

def save_geotiff(filename=None,
                 x_grid=None, y_grid=None, v_grid=None,
                 epsg="EPSG:32629",
                 outinfo=True):
    """
    Save interpolated data pr model parameter to geotiff

    Parameters
    ----------

    Filename : string
        Filename

    Returns
    -------
    None.

    """
    import rasterio
    import rasterio.transform
    import rasterio.crs

    error("Not yest worling! Exit.")
    #definition of the raster transform array
    nv = numpy.shape(v_grid)

    dx =numpy.abs(x_grid[1,0]-x_grid[0,0])
    dy =numpy.abs(y_grid[1,0]-y_grid[0,0])
    transform = rasterio.transform.Affine.translation(x_grid[0][0]-dx/2, y_grid[0][0]-dy/2)\
                *rasterio.transform.Affine.scale(dx,dy)

    #get crs as wkt
    crs = rasterio.crs.CRS.from_string(epsg)
    if outinfo: print(" CRS for input image: ", crs)
    #definition, register and close of interpolated raster
    rast_geo = rasterio.open(filename, "w", driver="GTiff",
                                    height=nv[0], width=nv[1],
                                    count=1,
                                    dtype=v_grid.dtype,
                                    crs=crs,
                                    transform=transform,
                                    )
    rast_geo.write(v_grid,1)
    rast_geo.close()
    if outinfo: print("Geotiff file written to ", filename)



def make_kml(kmzfile,
             llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""
    # from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
    #                    AltitudeMode, Camera)
    kml = simplekml.Kml()
    altitude = kw.pop("altitude", 2e7)
    roll = kw.pop("roll", 0)
    tilt = kw.pop("tilt", 0)
    altitudemode = kw.pop("altitudemode", simplekml.AltitudeMode.relativetoground)
    camera = simplekml.Camera(latitude=numpy.mean([urcrnrlat, llcrnrlat]),
                    longitude=numpy.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name="GroundOverlay")
        ground.draworder = draworder
        ground.visibility = kw.pop("visibility", 1)
        ground.name = kw.pop("name", "overlay")
        ground.color = kw.pop("color", "9effffff")
        ground.atomauthor = kw.pop("author", "ocefpaf")
        ground.latlonbox.rotation = kw.pop("rotation", 0)
        ground.description = kw.pop("description", "Matplotlib figure")
        ground.gxaltitudemode = kw.pop("gxaltitudemode",
                                       "clampToSeaFloor")
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name="ScreenOverlay")
        screen.icon.href = colorbar
        screen.overlayxy = simplekml.OverlayXY(x=0, y=0,
                                     xunits=simplekml.Units.fraction,
                                     yunits=simplekml.Units.fraction)
        screen.screenxy = simplekml.ScreenXY(x=0.015, y=0.075,
                                   xunits=simplekml.Units.fraction,
                                   yunits=simplekml.Units.fraction)
        screen.rotationXY = simplekml.RotationXY(x=0.5, y=0.5,
                                       xunits=simplekml.Units.fraction,
                                       yunits=simplekml.Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = simplekml.Units.fraction
        screen.size.yunits = simplekml.Units.fraction
        screen.visibility = 1

    # kmzfile = kw.pop("kmzfile", "overlay.kmz")
    kml.savekmz(kmzfile)

def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = numpy.cos(numpy.mean([llcrnrlat, urcrnrlat]) * numpy.pi/180.0)
    xsize = numpy.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = numpy.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        matplotlib.pyplot.ioff()  # Make `True` to prevent the KML components from poping-up.

    fig = matplotlib.pyplot.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax
