#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
# ---


"""
Created on 2020/11/11

This script present a work flow for processing data for inversions with
INV.py scripts

@author: vrath Oct 2020
"""

import os
import sys
from sys import exit as error
from time import process_time
from datetime import datetime
from random import randrange

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [AEMPYX_ROOT+"/aempy/modules/", AEMPYX_ROOT+"/aempy/scripts/"]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import prep
import aesys

OutInfo = True
AEMPYX_DATA = os.environ["AEMPYX_DATA"]

version, _ = versionstrg()
now = datetime.now()
Strng = "AEMpyX Version "+version
print("\n\n"+Strng)
print("Tellus NM CGG-GENESIS "
      +"\n"+"".join("Date " + now.strftime("%m/%d/%Y, %H:%M:%S")))
print("\n\n")


DatErr_add = 45.
DatErr_mult = 0.
kmax = 5
NanOut = False

AEM_system = "genesis"
FwdCall,NN,_,_,_, = aesys.get_system_params(System=AEM_system)
nD = NN[0]

impute = ["noise", 1000.]
# impute = ["delete", 0.]
"""
Input formats are 'npz','nc4','ascii'
"""
Filelist = "set" # "search", "read"
InFileFmt = ".npz"
OutFileFmt = ".npz"

# InDatDir = AEMPYX_DATA+"/Blocks/NM/raw/"
# OutDatDir = AEMPYX_DATA+"/Blocks/NM/proc_"+impute[0]+"/"
# SearchStrng = "NM_*FL*"


InDatDir = AEMPYX_DATA+"/RAFTA/overlaps_raw/data_asc/"
OutDatDir = AEMPYX_DATA+"/RAFTA/overlaps_"+impute[0]+"/"
SearchStrng = "NM_*FL*"

print("Data read from dir: %s " % InDatDir)
print("Data written to dir: %s " % OutDatDir)

dat_files = []

if "read" in FileList.lower()):
    dat_files = []
    with open("Filelist.txt", 'r') as file:
        for line in file:
            dat_files.append(line[:-1])

if "search" in FileList.lower()):
    dat_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InDatDir)


dat_files = util.get_filelist(searchstr=[SearchStrng], searchpath=InDatDir)
dat_files = sorted(dat_files)
ns = numpy.size(dat_files)
if ns ==0:
    error("No files corresponding to searchstring <"+SearchStrng+"> found!. Exit.")

run_number =str(randrange(10000))
files_to_do = dat_files.copy()
with open(InDatDir+"data_files_"+run_number+".txt", "w") as file:
    for item in files_to_do:
        file.write('%s\n' % item)


if not os.path.isdir(OutDatDir):
    print("File: %s does not exist, but will be created" % OutDatDir)
    os.mkdir(OutDatDir)


start = process_time()
num_sites = 0
bad_files = 0
for filename in dat_files:
    name, ext = os.path.splitext(filename)

    filein = InDatDir + filename
    print("\n Preprocessing file " + filein)
    Data, Header = aesys.read_aempy(File=filein, Format=InFileFmt,
                                   System=AEM_system, OutInfo=False)

    now = datetime.now()
    Header = aesys.grow_header(Header, ["Processed " + now.strftime("%m/%d/%Y, %H:%M:%S")])

    if numpy.size(Data)<=nD:
        print("Not enough data! Not processed")
        continue

    D = Data[:, :]

    D = numpy.where(D<1.e30,D, numpy.nan)
    nN = numpy.count_nonzero(numpy.isnan(D))
    print (str(nN)+" NaN in Data Block")
    if nN >0:
        bad_files = bad_files+1
        continue


    fline = Data[:, 0]
    sizedat = numpy.shape(D)
    nvars = sizedat[1]
    last = nvars - 1
    print("Flightline Data block on input has shape: ", numpy.shape(D))

    action = "lowpass filter"
    print("\n Proc action: " + action)
    columns = [4, 5]
    print(" dcolumns: ", columns)
    Header = aesys.grow_header(Header, ["LPF, IIR n=8"])
    D, _ = prep.filter_column(D, columns, method=["butter", 4, 1.0 / 20.0])
    print(" data block now has shape: ", numpy.shape(D))

    action = "greater than"
    threshval = 150.0
    columns = [4, 4]
    print("\n Proc action: " + action)
    print(" columns: ", columns)
    print(" thresh = ", threshval)
    Header = aesys.grow_header(
        Header, ["ALT, threshold = " + str(threshval)])
    D, nanindex = prep.insert_flag(D, action, threshval, columns,
                                   System = AEM_system)
    if NaNOut:
        head = Header
        filout = OutDatDir + name +"_nan"+OutFileFmt
        aesys.write_aempy(File=filout, Data=D, System=AEM_system,
                                Header=head, OutInfo=False)
        print("Data written to File: " + filout + "." + OutFileFmt.lower()[0:3])
        print("Info:")
        print(head)


    columns = [6, 28]
    print(" dcolumns: ", columns)
    D = prep.handle_gaps(D, columns, Impute=impute)
    print(" data block now has shape: ", numpy.shape(D))

    columns = []
    D = prep.handle_gaps(D, columns, Impute=impute, System=AEM_system)
    print(" data block now has shape: ", numpy.shape(D))

    if numpy.shape(D)[0] == 0:
        continue


    filout = OutDatDir + name+"_"+impute[0]+OutFileFmt
    aesys.write_aempy(File=filout, Data=D, System=AEM_system,
                        Header=head, OutInfo=False)

    print("Imputed data written to File: " + filout + "." + OutFileFmt.lower()[0:3])
    print("Info:")
    print(head)
    print("time taken = ", process_time() - start, "s \n")



    nDfinal = numpy.shape(D)
    num_sites = num_sites + nDfinal[0]
     if nDfinal[0]==0:
        break

    print("\nRunning pca filter ")
    columns = [6, 28]
    ncols = numpy.size(range(columns[0], columns[1] + 1))
    F = numpy.zeros(kmax)

    k = 0
    M = numpy.zeros(kmax)
    SVals = numpy.nan * numpy.ones((0, 22))
    MVals = numpy.nan * numpy.ones((0, 22))
    while k < kmax:


        k = k + 1
        print(" N pca: ", k)

        Data_k, U, S, V, MSE, FRO = prep.calc_svd_decomp(D, columns, k=k,
                                                  out_full=True)
        S = S / S[0]
        F[k-1] = FRO
        if OutInfo:
            print("TSVD: "+" k="+str(k)+" S(rel)="+str(S)+" FRO="+str(FRO))



        head = aesys.grow_header(
            Header,
            [   "\nTSVD: "+"SVals(norm) ="+str(S),
                "\nTSVD: "+" k="+str(k)+" S(rel)="+str(S)+" FRO="+str(FRO)])


        filout = OutDatDir + name+"_"+impute[0]+"_k" + str(k)+OutFileFmt
        aesys.write_aempy(File=filout, Data=Data_k, System=AEM_system,
                        Header=head, OutInfo=False)

        print("Data written to File: " + filout + "." + OutFileFmt.lower()[0:3])
        print("Info:")
        print(head)

    print("time taken = ", process_time() - start, "s \n")

    files_to_do.pop(0)
    with open(InDatDir+"data_files_to_do"+run_number+".txt", "w") as file:
        for item in files_to_do:
            file.write('%s\n' % item)


elapsed = process_time() - start
print(" Used %7.4f sec for %6i lines  - %8i sites\n" % (elapsed, ns, num_sites))
print(str(bad_files)+" bad files containing NaNs found.")
