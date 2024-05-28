#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 20:18:33 2022

@author: vrath
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
import warnings
import csv
import simplekml

import numpy
import shapely.geometry as shg

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]

for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0, pth)

from version import versionstrg
import util
import aesys

warnings.simplefilter(action="ignore", category=FutureWarning)


AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus process Clara ERT "
      + "".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")

InputDir = "/home/vrath/work/Clara/work/"
OutputDir = InputDir

if not os.path.isdir(OutputDir):
    print(" File: %s does not exist, but will be created" % OutputDir)
    os.mkdir(OutputDir)

TellusAng = 345.
AngLimits = [TellusAng-15., TellusAng+15. ]


Profiles = "Clara_elec.csv"
Data = []
Name = []
X = []
Y = []
Z = []
with open(InputDir+Profiles) as fd:
      rows = csv.reader(fd, delimiter=',', )
      h = next(rows)
      if h != None:
          for row in rows:
             Data.append(row)
             Name.append( row[0])
             X.append(float(row[1]))
             Y.append(float(row[2]))
             Z.append(float(row[3]))

Data = numpy.array(Data, dtype = "object")
Xitm = numpy.array(X)
Yitm = numpy.array(Y)
Z = numpy.array(Z)

EPSG_in = 29903
Lat, Lon= util.project_utm_to_latlon(Xitm, Yitm,
                                    utm_zone=EPSG_in)
EPSG_in = 29903
EPSG_out = 32629
Xutm, Yutm = util.project_utm_to_utm(Xitm, Yitm,
                                    utm_zone_in=EPSG_in,
                                    utm_zone_out=EPSG_out)

Data[:,1] = Xutm
Data[:,2] = Yutm
Data[:,3] = Z
Data[:,4] = numpy.zeros_like(Z, dtype=int)


P1 = []
P2 = []
P3 = []
Profiles = "Clara_elec_utm29N.csv"
with open(InputDir+Profiles, 'w') as fd:
    print(*h, sep=',', file=fd)
    for row in Data:
        print(*row, sep=',', file=fd)
        if row[0].split()[0]=="1":
            P1.append(numpy.asarray(row[1:-1]))
        if row[0].split()[0]=="2":
            P2.append(numpy.asarray(row[1:-1]))
        if row[0].split()[0]=="3":
            P3.append(numpy.asarray(row[1:-1]))



P1 = numpy.asarray(P1, dtype=float)
ex, ey = P1[:, 0], P1[:, 1]
nd =numpy.shape(P1)[0]
spoint = [ex[round(nd*0.3)], ey[round(nd*0.3)]]
epoint = [ex[round(nd*0.6)], ey[round(nd*0.6)]]
ang, _ = util.get_direction_angle(spoint, epoint)
if ang < AngLimits[0] or ang > AngLimits[1]:
    P1 = numpy.flipud(P1)
    print(" Angle = "+str(round(ang,1))
        +" not in interval "
        +str(round(AngLimits[0],1))+" - "
        +str(round(AngLimits[1],1)))
    print("ERT line direction has been reversed.")
else:
    print("ERT line direction is approx. 345 degrees")

P2 = numpy.asarray(P2, dtype=float)
ex, ey = P2[:, 0], P2[:, 1]
nd =numpy.shape(P2)[0]
spoint = [ex[round(nd*0.3)], ey[round(nd*0.3)]]
epoint = [ex[round(nd*0.6)], ey[round(nd*0.6)]]
ang, _ = util.get_direction_angle(spoint, epoint)
if ang < AngLimits[0] or ang > AngLimits[1]:
    P2 = numpy.flipud(P2)
    print(" Angle = "+str(round(ang,1))
        +" not in interval "
        +str(round(AngLimits[0],1))+" - "
        +str(round(AngLimits[1],1)))
    print("ERT line direction has been reversed.")
else:
    print("ERT line direction is approx. 345 degrees")


P3 = numpy.asarray(P3, dtype=float)
ex, ey = P3[:, 0], P3[:, 1]
nd =numpy.shape(P3)[0]
spoint = [ex[round(nd*0.3)], ey[round(nd*0.3)]]
epoint = [ex[round(nd*0.6)], ey[round(nd*0.6)]]
ang, _ = util.get_direction_angle(spoint, epoint)
if ang < AngLimits[0] or ang > AngLimits[1]:
    P3 = numpy.flipud(P3)
    print(" Angle = "+str(round(ang,1))
        +" not in interval "
        +str(round(AngLimits[0],1))+" - "
        +str(round(AngLimits[1],1)))
    print("ERT line direction has been reversed.")
else:
    print("ERT line direction is approx. 345 degrees")



AEM_system = "genesis"
td_file = "NM_A1_intersection_FL13490-0.asc"
filein = InputDir+"/AEM_TD/"+td_file
print("\n Reading file " + filein)
tddata, tdheader = aesys.read_aempy(File=filein,
                               System=AEM_system, OutInfo=False)

nd = numpy.shape(tddata)[0]

tdx, tdy = tddata[:,1], tddata[:,2]

spoint = [tdx[round(nd*0.3)], tdy[round(nd*0.3)]]
epoint = [tdx[round(nd*0.6)], tdy[round(nd*0.6)]]
ang, _ = util.get_direction_angle(spoint, epoint)
if ang < AngLimits[0] or ang > AngLimits[1]:
    tddata = numpy.flipud(tddata)
    print(" Angle = "+str(round(ang,1))
        +" not in interval "
        +str(round(AngLimits[0],1))+" - "
        +str(round(AngLimits[1],1)))
    print("TD line direction has been reversed.")
else:
    print("TD line direction is approx. 345 degrees")

tdx, tdy = tddata[:,1], tddata[:,2]

ertdata = P1
erx, ery = ertdata[:,0], ertdata[:,1]
p0 = util.get_nearest_point(point=[erx[0], ery[0]], line=[tdx, tdy])[0]
p1 = util.get_nearest_point(point=[erx[-1], ery[-1]], line=[tdx, tdy])[0]
tddata = tddata[p0:p1,:]

name, ext = os.path.splitext(filein)
fileout = name+"_prof1"+ext
aesys.write_aempy(File=fileout, System=AEM_system, OutInfo=False,
                  Header=tdheader, Data=tddata)


AEM_system = "aem05"
fd_file = "A1_NM_intersection_FL11126-0_delete_PLM10_k5.asc"
filein = InputDir+"/AEM_FD/"+fd_file
print("\n Reading file " + filein)
fddata, fdheader = aesys.read_aempy(File=filein,
                               System=AEM_system, OutInfo=False)
nd = numpy.shape(fddata)[0]
fdx, fdy = fddata[:,1], fddata[:,2]
spoint = [fdx[round(nd*0.3)], fdy[round(nd*0.3)]]
epoint = [fdx[round(nd*0.6)], fdy[round(nd*0.6)]]
ang, _ = util.get_direction_angle(spoint, epoint)
if ang < AngLimits[0] or ang > AngLimits[1]:
    fddata = numpy.flipud(fddata)
    print(" Angle = "+str(round(ang,1))
        +" not in interval "
        +str(round(AngLimits[0],1))+" - "
        +str(round(AngLimits[1],1)))
    print("FD line direction has been reversed.")
else:
    print("FD line direction is approx. 345 degrees")

fdx, fdy = fddata[:,1], fddata[:,2]

ertdata = P1
erx, ery = ertdata[:,0], ertdata[:,1]
p0 = util.get_nearest_point(point=[erx[0], ery[0]], line=[fdx, fdy])[0]
p1 = util.get_nearest_point(point=[erx[-1], ery[-1]], line=[fdx,fdy])[0]
fddata = fddata[p0:p1,:]

name, ext = os.path.splitext(filein)
fileout = name+"_prof1"+ext
aesys.write_aempy(File=fileout, System=AEM_system, OutInfo=False,
                  Header=fdheader, Data=fddata)





#  Determine what is added to the KML-tags:
kml = simplekml.Kml(open=1)
# plots_include = True
# plots_search = "png"



# Define the path for saving  kml files
kml_dir = InputDir
kml_file = kml_dir+"ClaraData"
writekml = False
writekmz = True
icon_dir = AEMPYX_ROOT+"/aempy/share/icons/"

ert_icon =  icon_dir + "triangle.png"
ert_tcolor = simplekml.Color.yellow  # "#555500" #
ert_tscale = 1.5  # scale the text
ert_icolor = simplekml.Color.yellow
ert_iscale = 1.0

aem_icon_fd =  icon_dir + "circle.png"
aem_icon_td =  icon_dir + "square.png"
aem_icolor_fd = simplekml.Color.blue
aem_tcolor_fd = simplekml.Color.blue
aem_icolor_td = simplekml.Color.red
aem_tcolor_td = simplekml.Color.red
aem_iscale = 1.0
aem_tscale = 1.0

# simplekml.Color.rgb(0, 0, 255)  # "ffff0000"
ert_iref = kml.addfile(ert_icon)
aem_iref_fd = kml.addfile(aem_icon_fd)
aem_iref_td = kml.addfile(aem_icon_td)


EPSG_in = 32629
# Determine which geographical info is added to the KML-tags:
# define empty list
iprof= -1
profnames = ["P1", "P2", "P2"]
for ertdata in [P1, P2, P3]:
    iprof = iprof+1
    pname=profnames[iprof]
    print(pname)
    lat, lon= util.project_utm_to_latlon(ertdata[:,0], ertdata[:,1],
                                    utm_zone=32629)
    for ii in numpy.arange(len(lat)):
        # print(ii)
        if ii==0:
            clara_data = kml.newpoint(name=pname)
        else:
            clara_data = kml.newpoint()

        clara_data.coords = [(lon[ii], lat[ii])]
        clara_data.style.labelstyle.color = ert_tcolor
        clara_data.style.labelstyle.scale = ert_tscale
        clara_data.style.iconstyle.icon.href = ert_iref
        clara_data.style.iconstyle.scale = ert_iscale
        clara_data.style.iconstyle.color = ert_icolor
        clara_data.description = "" #pname+str(ii)


pname ="FL13490"
lat, lon= util.project_utm_to_latlon(tddata[:,1], tddata[:,2],
                                utm_zone=32629)

for ii in numpy.arange(len(lat)):
    if ii==len(lat)-1:
        print(ii)
        clara_data = kml.newpoint(name=pname)
    else:
        clara_data = kml.newpoint()

    clara_data.coords = [(lon[ii], lat[ii])]
    clara_data.style.labelstyle.color = aem_tcolor_td
    clara_data.style.labelstyle.scale = aem_tscale
    clara_data.style.iconstyle.icon.href = aem_iref_td
    clara_data.style.iconstyle.scale = aem_iscale
    clara_data.style.iconstyle.color = aem_icolor_td
    clara_data.description = "" #pname+str(ii)

pname ="FL11126"
lat, lon= util.project_utm_to_latlon(fddata[:,1], fddata[:,2],
                                utm_zone=32629)

for ii in numpy.arange(len(lat)):

    if ii==0:
        print(ii)

        clara_data = kml.newpoint(name=pname)
    else:
        clara_data = kml.newpoint()

    clara_data.coords = [(lon[ii], lat[ii])]
    clara_data.style.labelstyle.color = aem_tcolor_fd
    clara_data.style.labelstyle.scale = aem_tscale
    clara_data.style.iconstyle.icon.href = aem_iref_fd
    clara_data.style.iconstyle.scale = aem_iscale
    clara_data.style.iconstyle.color = aem_icolor_fd
    clara_data.description = "" #pname+str(ii)


# Compressed kmz file:
kml.savekmz(kml_file + ".kmz")
