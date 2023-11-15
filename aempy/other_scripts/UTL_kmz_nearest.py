#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---

"""
@author: vr May 2022
"""

# Import required modules

import os
import sys
from sys import exit as error
import csv
import warnings
from time import process_time
from datetime import datetime
import simplekml
import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]

for pth in mypath:
    if pth not in sys.path:
        # sys.path.append(pth)
        sys.path.insert(0,pth)


import aesys
import util
import viz

from version import versionstrg


OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]


version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Generate KML for GENESIS-AEM05 pairs "+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


# Define the path to your data-files
data_dir = AEMPYX_DATA+"/Nearest/"
print(" data files read from: %s" % data_dir)
indatdir = data_dir+"/raw/"
inpltdir = data_dir+"/plot/"

if (not os.path.isdir(indatdir)) or (not os.path.isdir(inpltdir)):
    error(" File: %s or %s does not exist! Exit." % (indatdir, inpltdir))

# Fline_files = [data_dir+"TD_TB_Nearest_FDLines.txt",
#                data_dir+"TD_A1_Nearest_FDLines.txt",
#                data_dir+"TD_A2_Nearest_FDLines.txt",]

Fline_files = [data_dir+"TD_Nearest_FDLines.txt"]
Fline_files = [data_dir+"TD_Nearest_Somaye_2.txt",]

TellusAng = 345.
AngLimits = [TellusAng-5., TellusAng+5. ]
InvertDirection = True

for ii in range(len(Fline_files)):
    lines = numpy.loadtxt(Fline_files[ii], dtype=object)

    if ii==0:
        FLines =lines
    else:
        FLines = numpy.vstack((FLines, lines))


# Determine what is added to the KML-tags:

plots_td_include = False
plots_fd_include = False

plots_fmt = ".png"

kml_dir = data_dir+"/kmz/"
kml_file = kml_dir+"Nearest_Somaye_2_NoPlots"




kml = False
kmz = True

# Define the path for saving  kml files

icon_dir = AEMPYX_ROOT+"/aempy/share/icons/"

pltstep = 20
PlotMax = 999999999

line_icon =  icon_dir + "star.png"
data_icon_fd =  icon_dir + "triangle.png"
data_icon_td =  icon_dir + "square.png"

line_tcolor = simplekml.Color.white  # "#555500" #
line_tscale = 2.5  # scale the text

line_iscale = 2.
line_icolor = simplekml.Color.yellow
line_rcolor = simplekml.Color.yellow

data_iscale = 0.8
data_tscale = 1.2
data_icolor_fd = simplekml.Color.blue
data_rcolor_fd = simplekml.Color.blue

data_icolor_td = simplekml.Color.red
data_rcolor_td = simplekml.Color.red

# simplekml.Color.rgb(0, 0, 255)
# "ffff0000"


# Determine which geographical info is added to the KML-tags:
# define empty list
kml = simplekml.Kml(open=1)


line_iref = kml.addfile(line_icon)
fdata_iref = kml.addfile(data_icon_fd)
tdata_iref = kml.addfile(data_icon_td)

plot_files = os.listdir(inpltdir)

num_fl =-1
for f in FLines:
    num_fl =num_fl+1
    if num_fl>PlotMax:
        break
    print("Processing pair #"+str(num_fl))


    AEM_system = "aem05"
    fd_name,  ext = os.path.splitext(f[1])
    fd_file = indatdir+f[1]
    print("fd_file = "+fd_file)
    Data, _ = aesys.read_aempy(File=fd_file, System=AEM_system, OutInfo=False)
    if InvertDirection:
        nd =numpy.shape(Data)[0]
        spoint = [Data[round(nd*0.3),1], Data[round(nd*0.3),2]]
        epoint = [Data[round(nd*0.6),1], Data[round(nd*0.6),2]]
        ang, _ = util.get_direction_angle(spoint, epoint)
        if ang < AngLimits[0] or ang > AngLimits[1]:
            Data = numpy.flipud(Data)
            print(" Angle = "+str(round(ang,1))
                +" not in interval "
                +str(round(AngLimits[0],1))+" - "
                +str(round(AngLimits[1],1)))
            print("FD flightline direction has been reversed.")
        else:
            print("FD flightline direction is approx. 345 degrees")
    fdata = Data
    nfd =numpy.shape(fdata)[0]
    flat, flon = util.project_utm_to_latlon(fdata[:,1], fdata[:,2])
    fline = str(round(fdata[0,0],2)).replace(".","-")

    AEM_system = "genesis"
    td_file = indatdir+f[0]
    td_name,  ext = os.path.splitext(f[0])
    print("td_file = "+td_file)
    Data, _ = aesys.read_aempy(File=td_file, System=AEM_system, OutInfo=False)
    if InvertDirection:
        nd =numpy.shape(Data)[0]
        spoint = [Data[round(nd*0.3),1], Data[round(nd*0.3),2]]
        epoint = [Data[round(nd*0.6),1], Data[round(nd*0.6),2]]
        ang, _ = util.get_direction_angle(spoint, epoint)
        if ang < AngLimits[0] or ang > AngLimits[1]:
            Data = numpy.flipud(Data)
            print(" Angle = "+str(round(ang,1))
                +" not in interval "
                +str(round(AngLimits[0],1))+" - "
                +str(round(AngLimits[1],1)))
            print("TD flightline direction has been reversed.")
        else:
            print("TD flightline direction is approx. 345 degrees")
    tdata = Data
    ntd = numpy.shape(tdata)[0]
    tlat, tlon = util.project_utm_to_latlon(tdata[:,1], tdata[:,2])
    tline = str(round(tdata[0,0],2)).replace(".","-")

    folder_pair = kml.newfolder(name="FL"+tline+"_FL"+fline)

    if plots_fd_include:
        fd_plot = inpltdir+fd_name+ plots_fmt
        fd = folder_pair.newpoint(name=str(fdata[round(nfd/2),0]))
        fd.coords = [(flon[round(nfd/2)], flat[round(nfd/2)])]
        fd.style.labelstyle.color = data_icolor_fd
        fd.style.labelstyle.scale = data_tscale
        fd.style.iconstyle.icon.href = line_iref
        fd.style.iconstyle.scale = line_iscale
        fd.style.iconstyle.color = data_icolor_fd
        src= kml.addfile(fd_plot)
        fd.description = ('<img width="1200" align="left" src="' + src + '"/>')

    if plots_td_include:
        td_plot = inpltdir+td_name+ plots_fmt
        td = folder_pair.newpoint(name=str(tdata[round(ntd/2),0]))
        td.coords = [(tlon[round(ntd/2)], tlat[round(ntd/2)])]
        td.style.labelstyle.color = data_icolor_td
        td.style.labelstyle.scale = data_tscale
        td.style.iconstyle.icon.href = line_iref
        td.style.iconstyle.scale = line_iscale
        td.style.iconstyle.color = data_icolor_td
        src= kml.addfile(td_plot)
        td.description = ('<img width="1200" align="left" src="' + src + '"/>')


    folder_line = folder_pair.newfolder(name="FL"+fline)

    fd = folder_line.newpoint(name="S:"+str(fdata[0,0]))
    fd.style.labelstyle.color = data_icolor_fd
    fd.style.labelstyle.scale = data_tscale
    fd.style.iconstyle.icon.href = fdata_iref
    fd.style.iconstyle.scale = data_iscale*1.99
    fd.style.iconstyle.color = data_icolor_fd
    fd.coords = [(flon[0], flat[0])]

    fd = folder_line.newpoint(name="E:"+str(fdata[0,0]))
    fd.style.labelstyle.color = data_icolor_fd
    fd.style.labelstyle.scale = data_tscale
    fd.style.iconstyle.icon.href = fdata_iref
    fd.style.iconstyle.scale = data_iscale*1.99
    fd.style.iconstyle.color = data_icolor_fd
    fd.coords = [(flon[nfd-1], flat[nfd-1])]

    for ifd in numpy.arange(round(pltstep/2), nfd):
        if numpy.mod(ifd, pltstep) == 0:
            fd = folder_line.newpoint()
            fd.coords = [(flon[ifd], flat[ifd])]
            fd.style.iconstyle.icon.href = fdata_iref
            fd.style.iconstyle.scale = data_iscale
            fd.style.iconstyle.color = data_icolor_fd
            fd.description = AEM_system.upper()+"\nFlightline: "+str(fdata[ifd,0])


    folder_line = folder_pair.newfolder(name="FL"+tline)

    td = folder_line.newpoint(name="S:"+str(tdata[0,0]))
    td.style.labelstyle.color = data_icolor_td
    td.style.labelstyle.scale = data_tscale
    td.style.iconstyle.icon.href = tdata_iref
    td.style.iconstyle.scale = data_iscale*1.99
    td.style.iconstyle.color = data_icolor_td
    td.coords = [(tlon[0], tlat[0])]

    td = folder_line.newpoint(name="E:"+str(tdata[0,0]))
    td.style.labelstyle.color = data_icolor_td
    td.style.labelstyle.scale = data_tscale
    td.style.iconstyle.icon.href = tdata_iref
    td.style.iconstyle.scale = data_iscale*1.99
    td.style.iconstyle.color = data_icolor_td
    td.coords = [(tlon[ntd-1], tlat[ntd-1])]

    for itd in numpy.arange(ntd):

        if numpy.mod(itd, pltstep) == 0:
            td = folder_line.newpoint()
            td.coords = [(tlon[itd], tlat[itd])]
            td.style.iconstyle.icon.href = tdata_iref
            td.style.iconstyle.scale = data_iscale
            td.style.iconstyle.color = data_icolor_td
            td.description = AEM_system.upper()+"\nFlightline: "+str(tdata[itd,0])


    # Compressed kmz file:

kml.savekmz(kml_file + ".kmz")
