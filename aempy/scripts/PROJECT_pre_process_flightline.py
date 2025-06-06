#!/usr/bin/env python3
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py,ipynb
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
# This script controls preprocessing of data required ore advantageous for
# subsequent inversions.
#

import os
import sys
from sys import exit as error
import copy
import getpass
import inspect

from time import process_time
from datetime import datetime
# from random import randrange
# import warnings

import numpy

AEMPYX_ROOT = os.environ["AEMPYX_ROOT"]
mypath = [os.path.join(AEMPYX_ROOT, "aempy/modules/")]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import prep
import aesys

# +

AEMPYX_DATA = os.environ["AEMPYX_DATA"]
# -

version, _ = versionstrg()
script = inspect.getfile(inspect.currentframe())  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+"\n\n")
Header = titstrng

# Define some parameters required for the different systems.

# +
# AEM_system = "genesis"
AEM_system = "aem05"

if "aem05" in AEM_system.lower():
    _, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]

if "genes" in AEM_system.lower():
    _, NN, _, _, _, = aesys.get_system_params(System=AEM_system)
    nD = NN[0]


# +
OutInfo = True
OutNaN = False
OutRes = False

SingValMax = 5


# +
InFileFmt = ".npz"
OutFileFmt = ".npz"

# -


##############################################################################
# StGormans
##############################################################################
AEMPYX_DATA = "/home/vrath/work/A1_StGormans/"

Filelist = "search" # "set", "read"


# rectangle case
# InputDataDir =  AEMPYX_DATA + "/raw/"
# OutputDataDir =  AEMPYX_DATA + "/proc/"
# SearchStrng = "*FL*data.npz"

# lines case
InputDataDir =  AEMPYX_DATA + "/lines/"
OutputDataDir =  AEMPYX_DATA + "/lines/"
SearchStrng = "*FL*data.npz"

print("\n\n")
print("Data read from dir:  %s" % InputDataDir)
print("Search flightline ID string: %s " % SearchStrng)
print("Processed data  written to dir: %s " % OutputDataDir)



if not os.path.isdir(OutputDataDir):
    print("File: %s does not exist, but will be created" % OutputDataDir)
    os.mkdir(OutputDataDir)


print("Data read from dir: %s " % InputDataDir)

dat_files = util.get_data_list(how=["search", SearchStrng, InputDataDir],
                              out= True, sort=True)
ns = numpy.size(dat_files)
if ns ==0:
    error("No files corresponding to searchstring <"+SearchStrng+"> found!. Exit.")

start = process_time()
num_sites = 0
num_files = 0
bad_files = 0

for filename in dat_files:
    num_files = num_files+1
    name, ext = os.path.splitext(filename)
    filein = InputDataDir + filename
    print("\n Preprocessing file " + filein)
    Data, Header, _ = aesys.read_aempy(File=filein,
                                   System=AEM_system, OutInfo=False)


    now = datetime.now()
    aesys.print_header(Header)
    Header = aesys.grow_header(Header, titstrng)


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

    action = "alt lowpass filter"
    print("\n Proc action: " + action)
    columns = [4, 5]
    print(" dcolumns: ", columns)
    Header = aesys.grow_header(Header, "LPF, IIR n=8")
    D, comment = prep.filter_column(D, columns, method=["butter", 4, 1.0 / 20.0])
    print(" data block now has shape: ", numpy.shape(D))

    action = " plm lowpass filter"
    print("\n Proc action: " + action)
    columns = [14, 14]
    print(" dcolumns: ", columns)
    Header = aesys.grow_header(Header, "LPF, IIR n=8")
    D, comment = prep.filter_column(D, columns, method=["butter", 4, 1.0 / 20.0])
    print(" data block now has shape: ", numpy.shape(D))

    action = "plm threshold "
    plmthresh = 0.25
    threshval = plmthresh
    columns = [14, 14]
    print("\n Proc action: " + action)
    print(" columns: ", columns)
    print(" thresh = ", threshval)
    Header = aesys.grow_header(
        Header, "PLM, threshold = " + str(threshval))
    D, nanindex = prep.insert_flag(D, action, threshval, columns,
                                    System=AEM_system)

    action = "less than"
    threshval = 0.0
    columns = [6, 14]
    print("\n Proc action: " + action)
    print(" columns: ", columns)
    print(" thresh = ", threshval)
    Header = aesys.grow_header(
        Header, "DAT, threshold = " + str(threshval))
    D, nanindex = prep.insert_flag(D, action, threshval, columns,
                                   System=AEM_system)
    action = "greater than"
    threshval = 120.0
    columns = [4, 4]
    print("\n Proc action: " + action)
    print(" columns: ", columns)
    print(" thresh = ", threshval)
    Header = aesys.grow_header(
        Header, "ALT, threshold = " + str(threshval))
    D, nanindex = prep.insert_flag(D, action, threshval, columns,
                                   System=AEM_system)

    print("Info:")
    print(Header)
    print("time taken = ", process_time() - start, "s \n")

    if OutNaN:
        OutNameStrng = name.replace("_data", "_proc_delete_nan_data")
        filout = OutputDataDir + OutNameStrng + OutFileFmt
        aesys.write_aempy(File=filout, Data=D, System=AEM_system,
                            Header=Header, OutInfo=False)
        print("Data with NaN (for plotting) written to File: " + filout)


    # impute = ["noise", 100.]
    action = "handle_gaps"
    impute = ["delete", 0.]
    print("Impute method:")
    print(impute)
    columns = [6, 14]
    Header = aesys.grow_header(
        Header, "GAP, method = " + impute[0])
    D = prep.handle_gaps(D, columns, Impute=impute, System=AEM_system)
    print(" data block now has shape: ", numpy.shape(D))

    if numpy.shape(D)[0] == 0:
        continue


    OutNameStrng = name.replace("_data", "_proc_"+impute[0]+"_data")
    filout = OutputDataDir + OutNameStrng + OutFileFmt
    aesys.write_aempy(File=filout, Data=D, System=AEM_system,
                    Header=Header, OutInfo=False)

    print("Imputed data written to File: " + filout)
    print("Info:")
    print(Header)
    print("time taken = ", process_time() - start, "s \n")

    nDfinal = numpy.shape(D)
    num_sites = num_sites + nDfinal[0]

    print("\nRunning pca ")
    columns = [6, 13]
    ncols = numpy.size(range(columns[0], columns[1] + 1))
    F = numpy.zeros(SingValMax)

    k = 0
    M = numpy.zeros(SingValMax)
    SVals = numpy.nan * numpy.ones((0, 8))
    MVals = numpy.nan * numpy.ones((0, 8))
    while k < SingValMax:

        k = k + 1
        print(" N pca: ", k)
        Data_k, U, S, V, MSE, FRO = prep.calc_svd_decomp(D, columns, k=k,
                                                  out_full=True)
        S = S / S[0]
        F[k-1] = FRO
        if OutInfo:
            print("TSVD: "+" k="+str(k)+" S(rel)="+str(S)+" FRO="+str(FRO))


        head = aesys.grow_header(
            Header,"TSVD: "+" k="+str(k)+" S(rel)="+str(S)+" FRO="+str(FRO))

        OutNameStrng = name.replace("_data", "_proc_"+impute[0]+"_k" + str(k) + "_data")
        filout = OutputDataDir + OutNameStrng+ OutFileFmt
        aesys.write_aempy(File=filout, Data=Data_k,
                        System=AEM_system, Header=head, OutInfo=False)
        print("Data written to File: " + filout)
        print("Info:")
        print(head)
        print("time taken = ", process_time() - start, "s \n")

        if OutRes:
            D_res = copy.deepcopy(D)
            nd1 = NN[1]
            nd2 = NN[1]+ NN[2]
            D_res[:,nd1:nd2] = D[:,nd1:nd2]-Data_k[:,nd1:nd2]
            print("min, max = ",numpy.amin(D_res[:,nd1:nd2]), ", ",numpy.amax(D_res[:,nd1:nd2]))
            print("std = ",numpy.std(D_res[:,nd1:nd2]))
            head = aesys.grow_header(
                Header,"TSVD: "+" k="+str(k)+" Res = min "+str(numpy.amin(D_res[:,nd1:nd2]))
                                                +" / max "+str(numpy.amax(D_res[:,nd1:nd2]))
                                                +" / std "+str(numpy.std(D_res[:,nd1:nd2])))

            OutNameStrng = name.replace("_data","_proc_"+impute[0]+"_k" + str(k)+"_res" + "_data")
            filout = OutputDataDir + OutNameStrng+ OutFileFmt
            aesys.write_aempy(File=filout, Data=D_res,
                            System=AEM_system, Header=head, OutInfo=False)
            print("Data written to File: " + filout)
            print("Info:")
            print(head)
            print("time taken = ", process_time() - start, "s \n")


print("\nAll done!")

elapsed = process_time() - start
print(" Used %7.4f sec for %6i lines  - %8i sites\n" % (elapsed, num_files + 1, num_sites))
