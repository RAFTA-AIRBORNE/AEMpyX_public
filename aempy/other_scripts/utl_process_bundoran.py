#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---

"""
@author: vr April 2023
"""

# Import required modules
#!/usr/bin/env python3

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
import scipy
import scipy.spatial


AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg

import aesys
import util
import inverse
import alg

warnings.simplefilter(action="ignore", category=FutureWarning)

AEMPYX_DATA = os.environ["AEMPYX_DATA"]

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
AEM_system = "genesis"
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


# Bonduran A (east)
Point0= [543011.18, 6033721.42]
Point1= [542480.00, 6035601.00,]

# Bonduran B (center)
Point0= [542908.48, 6033692.60 ]
Point1= [542378.97, 6035576.55]

# Bonduran C
# Point0= [542908.48, 6033692.60 ]
# Point1= [542378.97, 6035576.55

FileList = "search"  # "search", "read"
InDatDir =  AEMPYX_DATA + "/Projects/Compare_systems/data/"
# SearchStrng = "CGG*FL90*.asc"
SearchStrng = "SGL*.npz"

if "set" in FileList.lower():
    print("Data files read from dir:  %s" % InDatDir)
    dat_files = []

else:
    # how = ["search", SearchStrng, InDatDir]
    # how = ["read", FileList, InDatDir]
    dat_files = util.get_data_list(how=["search", SearchStrng, InDatDir],
                              out= True, sort=True)
    ns = numpy.size(dat_files)

ns = numpy.size(dat_files)
if ns ==0:
    error("No files set!. Exit.")


"""
Output format is ".npz"
"""
OutFileFmt = ".npz"
OutDatDir =  AEMPYX_DATA + "/Projects/Compare_systems/data_reduced/"
print("Models written to dir: %s " % OutDatDir)
if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


TellusAng = 345.
AngLimits = [TellusAng-5., TellusAng+5. ]
CorrectDirection = True

setE = setW = setC = []
avgE0 = avgW0 = avgC0= numpy.array([0., 0.])
avgE1 = avgW1 = avgC1= numpy.array([0., 0.])
ndatE = ndatW = ndatC = 0

for file in dat_files:
    

    name,  ext = os.path.splitext(file)

    input_file = InDatDir+file
    print(" Data input from = "+input_file)
    
    output_file = OutDatDir+name+"_reduced"+OutFileFmt
    print(" Data output to = "+output_file)
    
    Data, Header, _ = aesys.read_aempy(File=input_file, 
                                  System=AEM_system, 
                                  OutInfo=False)
    
    if CorrectDirection:
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
    
    
    Prof = Data[:,1:3]
    # print(numpy.shape(Prof))
    
    D_tree=scipy.spatial.KDTree(Prof, leafsize=10,
                                compact_nodes=True,
                                copy_data=True,
                                balanced_tree=True,
                                boxsize=None)
    dist0, index0 = D_tree.query(Point0, k=1)
    dist1, index1 = D_tree.query(Point1, k=1)
    
    Data = Data[index0:index1,:]
    Prof = Prof[index0:index1,:]
 
    print(numpy.shape(Prof))
    
    Dist0 =  Point0 -Prof[0,:]
    # print ( Point0)
    # print ( Prof[0,:])
    # print ( dist0, index0)
    print ("Start point: ", Dist0)
        
    Dist1 = Point1 - Prof[-1,:]    
    # print (Point1)
    # print ( Prof[-1,:])
    # print ( dist1, index1)
    print ("End point: ", Dist1 )
    
    nd = numpy.shape(Prof)[0]
    print(numpy.shape(Prof))
    prof0 = numpy.array([Prof[round(nd*0.3),0], Prof[round(nd*0.3),1]]).reshape(-1,1)
    prof1 = numpy.array([Prof[round(nd*0.6),0], Prof[round(nd*0.6),1]]).reshape(-1,1)
    if   Dist1[0]<-10.:
        setW.append(output_file)
        avgW0 = numpy.column_stack((avgW0, prof0))
        avgW1 = numpy.column_stack((avgW1, prof1))
        ndatW = ndatW+1
    elif Dist1[0]>10.:
   
        setE.append(output_file)        
        avgE0 = numpy.column_stack((avgE0, prof0))
        avgE1 = numpy.column_stack((avgE1, prof1))
        ndatE = ndatE+1
    else:
        setC.append(output_file)        
        avgC0 = numpy.column_stack((avgC0, prof0))
        avgC1 = numpy.column_stack((avgC1, prof1))
        ndatC = ndatC+1

    
    NewHeader = aesys.grow_header(
                Header,"Startindex="+str(index0)+", Endindex="+str(index1))

    aesys.write_aempy(File=output_file, Data=Data,
                            System=AEM_system, Header=NewHeader, OutInfo=False)
    print("Data written to File: " + output_file)

    
    DelSamp = 10.
    Prof = Prof-Prof[0]
    # Dist = numpy.sqrt(numpy.power(Prof[:,0], 2.0) + numpy.power(Prof[:,1], 2.0))
    
    """
    Now project to new profile line for interpolaation
    """
    
    LinePoint0 = [Data[round(nd*0.3),0], Data[round(nd*0.3),1]]
    LinePoint1 = [Data[round(nd*0.6),0], Data[round(nd*0.6),1]]
    
    # util.project_to_line(x, y, line)
    
    # NewHeader = aesys.grow_header(
    #             NewHeader,"Interpolated, sample = "+str(DelSamp)+"m")
    
    

print(numpy.shape(avgE0))
print(numpy.shape(avgE1))

avgC0 = numpy.sum(avgC0, axis=1)/ndatC
avgC1 = numpy.sum(avgC1, axis=1)/ndatC
avgE0 = numpy.sum(avgE0, axis=1)/ndatE
avgE1 = numpy.sum(avgE1, axis=1)/ndatE
avgW0 = numpy.sum(avgW0, axis=1)/ndatW
avgW1 = numpy.sum(avgW1, axis=1)/ndatW

numpy.savez_compressed(
    file=AEMPYX_DATA + "/Projects/Compare_systems/BundoranSubsets.npz", 
    sets=[setE, setC, setW], avgs=[avgE0, avgE1, avgC0, avgC1, avgW0, avgW1] )
        