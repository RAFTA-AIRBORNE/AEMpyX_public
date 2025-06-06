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
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

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
print("Generate KML for boreholes & geophysics "+"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


# Define the path to your data-files
site_dir = AEMPYX_DATA+"/Intersection/boreholes/"
print(" site files read from: %s" % site_dir)
# site_files = ["Verified_Boreholes_utm_NM_A1.npz",
#                 "Verified_Boreholes_utm_NM_A2.npz",
#                 "Verified_Boreholes_utm_NM_TB.npz"]
MinDepth = 50.
site_files = ["Verified_Boreholes_utm_NM_A1.npz",
                "Verified_Boreholes_utm_NM_A2.npz",
                "Verified_Boreholes_utm_NM_TB.npz"]


SearchRadius=500.
data_dir = site_dir+"/data/"
infmt = ".asc"
print(" data files read from: %s" % data_dir)


# Determine what is added to the KML-tags:

plots_include = True
plots_search = "png"

kml = False
kmz = True

# Define the path for saving  kml files
kml_dir = data_dir

icon_dir = AEMPYX_ROOT+"/aempy/share/icons/"
site_icon =  icon_dir + "star.png"
data_icon_fd =  icon_dir + "circle.png"
data_icon_td =  icon_dir + "square.png"

site_tcolor = simplekml.Color.white  # "#555500" #
site_tscale = 2.0  # scale the text

site_iscale = 2.
site_icolor = simplekml.Color.yellow
site_rcolor = simplekml.Color.yellow

data_iscale = 0.75
data_tscale = 1.
data_icolor_fd = simplekml.Color.blue
data_rcolor_fd = simplekml.Color.blue

data_icolor_td = simplekml.Color.red
data_rcolor_td = simplekml.Color.red

# simplekml.Color.rgb(0, 0, 255)
# "ffff0000"



# Determine which geographical info is added to the KML-tags:
# define empty list
sites = []
for f in site_files:
    print("reading "+site_dir+f)
    sdata = numpy.load(site_dir+f, allow_pickle=True)["Data"]
    pname=sdata[:,0]

    slat, slon = util.project_utm_to_latlon(sdata[:,1], sdata[:,2])
    num_site =-1
    for p in pname:
        this_site = "/Borehole_"+p+"/"
        num_site =num_site+1
        if sdata[num_site,3]<=MinDepth:
            print("\n\n\n*****Borehole depth = "+str(sdata[num_site,3])+" for "+p+" is < "+str(MinDepth)+"m*****\n\n\n")
            continue

        kml = simplekml.Kml(open=1)
        site = kml.newpoint(name=p)
        site_iref = kml.addfile(site_icon)
        site.coords = [(slon[num_site], slat[num_site])]
        site.style.labelstyle.color = site_tcolor
        site.style.labelstyle.scale = site_tscale
        site.style.iconstyle.icon.href = site_iref
        site.style.iconstyle.scale = site_iscale
        site.style.iconstyle.color = site_icolor
        site.description = "Borehole Identification:  "+p+ "\nDepth ="+ str(sdata[num_site,3])
        kml_file = kml_dir+this_site + "Borehole_"+p+ "_S"+ str(SearchRadius)+"m_M"+ str(MinDepth)+"m"


        fdata_iref = kml.addfile(data_icon_fd)
        tdata_iref = kml.addfile(data_icon_td)

        indir = data_dir+"/Borehole_"+p+"/"
        if not os.path.isdir(indir):
            error(" File: %s does not exist! Exit." % indir)

        if plots_include:
            description = ""
            files = os.listdir(indir)
            for entry in files:
                if plots_search in entry:
                    png_name = indir+entry
                    srcfile = kml.addfile(png_name)
                    descr = (
                        '<img width="1200" align="left" src="' + srcfile + '"/>'
                        )
                    description = description + descr

        site.description = description


        AEM_system = "aem05"
        fd_file = AEM_system.upper()+"_Borehole"+p+"_Datafile"+"_SearchRadius"+str(round(SearchRadius))+"m"
        print(fd_file)
        fdata, _ = aesys.read_aempy(File=indir+fd_file+infmt, System=AEM_system, OutInfo=False)
        fs = numpy.shape(fdata)

        flines = numpy.unique(fdata[:,0])
        print("FD fligthlines:")
        print(flines)
        pos = []
        for i in range(len(flines)):
            # Iterate over list items by index pos
            for ipos in range(fs[0]):
                # Check if items matches the given element
                if fdata[ipos,0] == flines[i]:
                    pos.append(ipos)
                    print("Index of first occurence of "+str(flines[i])
                      +" in the list is: "+str(ipos))

                    break

        nfd = fs[0]
        numfline = fdata[:,0]
        flat, flon = util.project_utm_to_latlon(fdata[:,1], fdata[:,2])
        for ifd in numpy.arange(nfd):
            # if

            if ifd in pos:
                fd = kml.newpoint(name=str(fdata[ifd,0]))
                fd.style.labelstyle.color = data_icolor_fd
                fd.style.labelstyle.scale = data_tscale
            else:
                fd = kml.newpoint()

            fd.coords = [(flon[ifd], flat[ifd])]
            fd.style.iconstyle.icon.href = tdata_iref
            fd.style.iconstyle.scale = data_iscale
            fd.style.iconstyle.color = data_icolor_fd
            fd.description = AEM_system.upper()+"\nFlightline: "+str(fdata[ifd,0])


        AEM_system = "genesis"
        td_file = AEM_system.upper()+"_Borehole"+p+"_Datafile"+"_SearchRadius"+str(round(SearchRadius))+"m"
        print(td_file)
        tdata, _ = aesys.read_aempy(File=indir+td_file+infmt, System=AEM_system, OutInfo=False)
        ts = numpy.shape(tdata)

        tlines = numpy.unique(tdata[:,0])
        print("TD fligthlines:")
        print(tlines)
        pos = []
        for i in range(len(tlines)):
            # Iterate over list items by index pos
            for ipos in range(ts[0]):
                # Check if items matches the given element
                if tdata[ipos,0] == tlines[i]:
                    pos.append(ipos)
                    print("Index of first occurence of "+str(tlines[i])
                      +" in the list is: "+str(ipos))
                    break

        ntd = ts[0]
        tlat, tlon = util.project_utm_to_latlon(tdata[:,1], tdata[:,2])

        for itd in numpy.arange(ntd):

            if itd in pos:
                td = kml.newpoint(name=str(tdata[itd,0]))
                td.style.labelstyle.color = data_icolor_td
                td.style.labelstyle.scale = data_tscale
            else:
                td = kml.newpoint()


            td.coords = [(tlon[itd], tlat[itd])]
            td.style.iconstyle.icon.href = tdata_iref
            td.style.iconstyle.scale = data_iscale
            td.style.iconstyle.color = data_icolor_td
            td.description = AEM_system.upper()+"\nFlightline: "+str(tdata[itd,0])

        # Compressed kmz file:
        kml.savekmz(kml_file + ".kmz")

#     if plots_2:
#         nam_2 = name
#         srcfile_2 = kml.addfile(plots_dir + nam_2 + '.png')
#         description_2 = (
#             '<img width="900" align="left" src="' + srcfile_2 + '"/>'
#         )
#         description = description + description_2
