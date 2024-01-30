# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 16:05:35 2020

@author: vrath
"""
import sys
import os
from sys import exit as error
from datetime import datetime


import numpy
import scipy.signal
import scipy.linalg
import scipy.interpolate

#from scipy.signal import medfilt, decimate

import aesys
import util
import inverse



def filter_column(
        M=None,
        columns=None,
        method=["butter", 4, 1./10],
        OutInfo=True):
    """
     Filter columns in data.

    """
    cols = range(columns[0], columns[1] + 1)

    meth = method[0]
    ford = method[1]
    cutf = method[2]
    headstring = "filter columns " + str(columns) + ", method = " + str(method)
    print(meth.lower)
    if meth.lower() == "butter":
        bb, ab = scipy.signal.butter(ford, cutf)
    else:
        error("Method " + method[0] + " not yet implemented! Nothing done.")

    for icol in cols:
        tmp_in = M[:, icol]

        tmp_out = scipy.signal.filtfilt(
            bb,
            ab,
            tmp_in,
            method="gust",
            padtype="odd",
            padlen=24)
        comment = "filter column " + str(icol) + ", method = " + str(method)
        # print(numpy.isnan(tmp_out))

        M[:, icol] = tmp_out
        if OutInfo:
            comment = "filter column " + \
                str(icol) + ", method = " + str(method)
            print(comment)

    return M, headstring


def insert_flag(
        M,
        criterion="neg",
        threshval=0.,
        columns=[0, 0],
        flag=None,
        incol=None,
        System = "aem05"):
    """
    Replace bad values in data by flag (nan)

    Returns  array with replacement and a boolean array for bad values

    Last change vr  july 19, 2021

    """
    AEM_system = System
    _,nD,_,_,_ = aesys.get_system_params(System=AEM_system)

    if flag is None:
        flag = numpy.nan

    cols = numpy.arange(columns[0], columns[1] + 1)

    T = M[:, cols]

    oldflags = numpy.where(numpy.isnan(T))

    if "neg" in criterion.lower():
        T[oldflags] = -9999.
        T[T <= 0] = flag
    elif "less" in criterion.lower():
        T[oldflags] = threshval - 9999.
        T[T <= threshval] = flag
    elif "great" in criterion.lower():
        T[oldflags] = threshval + 9999.
        T[T >= threshval] = flag
    elif "plm" in criterion.lower():

        T[oldflags] = threshval + 9999.
        T[numpy.abs(T) > threshval] = flag

    elif "nan" in criterion.lower():
        if incol is None:
            error(" no column  given!")
        newindex = numpy.isnan(M[:, incol])
        T[newindex, :] = flag

    T[oldflags] = flag
    M[:, cols] = T
    NaNindex = numpy.where(numpy.isnan(M))

    return M, NaNindex


def handle_gaps(M, columns=[0, 0], Impute=["delete"], System = "aem05"):
    """
    Replace bad values in data (nans).

    Valid options are nan, delete, average, median, noise, spl, and rbf.

    if spl:
     interpolates bad values marked by nans in M.
     method: (str or int) specifies the kind of interpolation as a string
     (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’) where
     ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline interpolation of
     first, second or third order) or as an integer specifying the order
     of the spline interpolator to use. Default is ‘linear’

    Last change vr  Sep 1, 2021

    """
    AEM_system = System
    _,nD,_,_,_ = aesys.get_system_params(System=AEM_system)

    method = Impute.copy()

    if columns != []:
        cols = range(columns[0], columns[1] + 1)
    else:
        method[0] = "delete"
        # print(cols)

    nanindex = numpy.where(numpy.isnan(M))
    nanrows = numpy.ravel(nanindex[0])
    nancols = numpy.ravel(nanindex[1])

    if "nan" in method[0].lower():
        print("Method: " + method[0] + " - nothing to do.")
        pass


    elif "del" in method[0].lower():
        print("Method: " + method[0])
        M = numpy.delete(M, nanrows, 0)


    elif "ave" in method[0].lower():
        print("Method: " + method[0])
        T = M[:, cols]
        val = numpy.nanmean(T, 0)
        for irow in nanrows:
            for icol in nancols:
                T[irow, icol] = val[icol]
        M[:, cols] = T


    elif "med" in method[0].lower():
        print("Method: " + method[0])
        T = M[:, cols]
        val = numpy.nanmedian(T, 0)
        for irow in nanrows:
            for icol in nancols:
                T[irow, icol] = val[icol]
        M[:, cols] = T

    elif "noi" in method[0].lower():
        print("Method: "+method[0])
        print("Noise of "+str(method[1])+" ppm added.")
        T = M[:, cols]
        shapeT = numpy.shape(T)
        v_blk = numpy.nanmean(T, 0)
        v_std = Impute[1]
        # print(numpy.shape(v_blk))
        # print(numpy.shape(v_std))
        val = v_blk + v_std * numpy.random.default_rng().normal(shapeT[0], shapeT[1])
        for irow in nanrows:
            for icol in nancols:
                T[irow, icol] = val[irow, icol]
        M[:, cols] = T

    elif "spl" in method[0].lower():
        print("Method: " + method[0])
        xd = M[:, 1]
        yd = M[:, 2]
        rd = numpy.sqrt(xd * xd + yd * yd)
        for var in range[cols]:
            var = M[:, var]
            if not numpy.all(numpy.isnan(var)):
                spl = scipy.interpolate.UnivariateSpline(xd, yd,
                                                        w=numpy.ones_like(xd),
                                                        k=Impute[1])
                T[:, var] = spl(rd)
                if "spln" in method[0].lower():

                   if Impute[2] == 0.:
                       v0 = numpy.nanstd(var)
                   else:
                       v0 = Impute(2)

                   T[:, var] = T[:, var] + v0 * numpy.random.standard_normal(numpy.shape(var))
                   print("Noise of "+str(v0)+" ppm added.")

            else:
                print("Channel "+str(var)+" all nan! nothing done.")

        M[:, cols] = T

    else:
        error("Method " + method + " not yet implemented! Nothing done.")

    return M


def calc_svd_decomp(D=None, columns=None, ErrD=None, k=1, thresh_MSE=0.,
              out_full=False, OutInfo=False):
    """
     PCA filter similar to Minsley et al.(2012)

     IN:
        M block of data (usually flightline)
        k number of PCs < ndata

        Last change vr Nov 20, 2020

    """

    if (D is None) or (columns is None):
        error("pcafilter: no data definded! ")

    cols = range(columns[0], columns[1] + 1)

    if ErrD is None:
        SqCovDi =  SqCovD =numpy.diag(numpy.ones_like(cols))
    else:
        SqCovDi = numpy.diag(numpy.ones_like(cols)/ErrD)
        SqCovD  = numpy.diag(numpy.ones_like(cols)*ErrD)


    if k > numpy.size(cols):
        print("Number of SVs k ="+str(k)+" > "+str(numpy.size(cols))+". k assumed!")
        k = numpy.size(cols)


    P = D.copy()
    T = D[:, cols]
    print(numpy.shape(T))
    T_blk = numpy.nanmean(T, 0)
    X = numpy.dot(T - T_blk,SqCovDi)
    U, S, Vt = scipy.linalg.svd(X, full_matrices=False)
    SS = numpy.diag(S)
    V = Vt.T

    Tk = numpy.dot(U[:, :k], numpy.dot(SS[:k, :k], V[:, :k].T))

    MSE = numpy.sqrt(numpy.mean((X - Tk)**2))

    FRO = numpy.linalg.norm(X - Tk, "fro")/numpy.linalg.norm(X, "fro")

    if OutInfo:
        print(" using " + str(k) + " PCs, |X-Xk|/|X| = %.6G" % (FRO))

    P[:, cols] = numpy.dot(Tk,SqCovD) + T_blk

    if out_full:
        return P, U, S, V, MSE, FRO
    else:
        return P


# def reduce_data0(Data=None, System="aem05", Method=["median", 5], OutInfo = True):
#     """
#      Averages data into blocks of size blocksize.

#      Last change vr  Sep 11, 2021

#     """
#     if Data.size==0:
#         error("No data! Exit.")

#     _, nD, _, _, _, = aesys.get_system_params(System)
#     nVal =nD[1]+nD[2]+nD[3]


#     meth, blocksize = Method
#     comment = "Method = " + meth + " with blocksize = " + str(blocksize)


#     if blocksize != 1:
#         print(comment)
#         start=int(blocksize/2)
#         step = blocksize
#         for col in numpy.arange(0,nVal):

#            tmp, comment  = process_datacolumn(Data[:, col], Method=Method)
#            if col==0:
#                data_out = tmp
#            else:
#                data_out = numpy.column_stack((data_out,tmp))

#     else:
#         print("Blocksize is 1! Nothing done.")

#     Data = data_out[start:-1:step]
#     comment = "method = " + meth + ", blocksize = " + str(blocksize)
#     if OutInfo:
#         print(comment)
#         print("New data shape: ", numpy.shape(Data))


#     return Data, comment


# def process_datacolumn0(Data, Method=["mean", 3], OutInfo =False):
#     """
#      Averages data into blocks of size blocksize.

#      Last change vr Nov 20, 2020

#     """
#     meth = Method[0]

#     if "mean" in meth.lower():
#         blocksize = Method[1]
#         exestr = "numpy.nanmean(tmp[i - blocksize:i])"
#     elif "med" in meth.lower():
#         blocksize = Method[1]
#         exestr = "numpy.nanmedian(tmp[i - blocksize:i])"
#     elif "dec" in meth.lower():
#         qf = Method[1]
#         ftype = "iir"
#         exestr = "scipy.signal.decimate(tmp,"+qf+", ftype=" + ftype + ",zero_phase=True)"
#     else:
#         error("method " + meth.lower() + " not implemented! Exit.")

#     comment = "method = " + Method[0] + " parameters = " + str(Method[1:])
#     tmp = Data
#     data_out = []
#     for i, d in enumerate(tmp[:]):
#         if (i % blocksize) == 0 and i != 0:
#             val = eval(exestr)
#             data_out.append(val)

#     return data_out, comment


def reduce_data(datavec=numpy.array([]), System="aem05",
                      Method=["mean", 5],
                      ErrEst = False,
                      OutInfo = True):
    """
      Averages datavec into blocks of size blocksize.

      Last change: vr  Feb 21, 2023

    """
    if datavec.size==0:
        error("No datavec! Exit.")

    dat_out = numpy.array([])
    err_out = numpy.array([])

    _, nD, _, _, _, = aesys.get_system_params(System)

    meth, blocksize  = Method

    if blocksize != 1:

        if blocksize%2==0:
            blocksize = blocksize+1
            Method[1] = blocksize



        for col in numpy.arange(nD[0]):

            if col==0:
                dat, err, comment  = process_column(datavec[:, col],
                                                   Method=Method)
                dat_out = dat

            else:

                dat, err, comment  = process_column(datavec[:, col],
                                   Method=Method, ErrEst=ErrEst)

                dat_out = numpy.column_stack((dat_out,dat))
                if ErrEst:
                    if col==1:
                        err_out = err
                    else:
                        err_out = numpy.column_stack((err_out,err))

    else:
        print("Blocksize is 1! Nothing done.")


    # dat_out = dat_out[(blocksize-1)//2:-1:blocksize, :]

    if ErrEst:
        err_out = numpy.nanmean(err_out, axis=1)
        # err_out = err_out[(blocksize-1)//2:-1:blocksize, :]

    comment = "method = " + meth + ", blocksize = " + str(blocksize)


    if OutInfo:
        print(comment)
        print("New dat shape: ", numpy.shape(dat_out))
        if ErrEst:
            print("Err shape: ", numpy.shape(err_out))

    if ErrEst:
        return dat_out, err_out , comment
    else:
        return dat_out, comment


def process_column(datavec=numpy.array([]),
                       Method=["mean", 3],
                       ErrEst = False,
                       OutInfo =False):
    """
     Averages datavec into blocks of size blocksize.

     Last change vr Nov 20, 2020

    """
    # if not ErrEst:
    #     exestr_err="continue"

    comment = "method = " + Method[0] + " parameters = " + str(Method[1:])

    if "mean" in Method[0].lower():
        blocksize = Method[1]
        if blocksize%2 == 0: blocksize=blocksize+1
        exestr_avg = "numpy.nanmean(tmp[i - blocksize:i])"
        if ErrEst:
            exestr_err = "numpy.nanstd(tmp[i - blocksize:i])"

    elif "med" in Method[0].lower():
        blocksize = Method[1]
        if blocksize%2 == 0: blocksize=blocksize+1
        exestr_avg = "numpy.nanmedian(tmp[i - blocksize:i])"
        if ErrEst:
            exestr_err, _ = "inverse.calc_mad(tmp[i - blocksize:i], valavg)"

    elif "dec" in Method[0].lower():
        qf = Method[1]
        ftype = "iir"
        exestr_avg = "scipy.signal.decimate(tmp,"+qf+", ftype=" + ftype + ",zero_phase=True)"
        if ErrEst:
            error("For method " + Method[0].lower() + " no error esimation possible! Exit.")
    else:
        error("method " + Method[0].lower() + " not implemented! Exit.")


    data_head = numpy.flipud(datavec[1:blocksize//2+1])
    data_tail = numpy.flipud(datavec[-blocksize//2:len(datavec)-1])
    tmp = numpy.concatenate([data_head, datavec, data_tail], axis=0)
    dat_out = []
    err_out = []
    for i, d in enumerate(tmp[:]):
        if (i % blocksize) == 0 and i != 0:
            valavg = eval(exestr_avg)
            dat_out.append(valavg)
            if ErrEst:
                valerr = eval(exestr_err)
                err_out.append(valerr)

    dat_out = numpy.asarray(dat_out)
    dat_err = numpy.array([])
    if ErrEst:
        dat_err = numpy.asarray(err_out)
    else:
        dat_err = numpy.array([])

    comment = "method = " + Method[0] + " parameters = " + str(Method[1:])



    return dat_out, dat_err,  comment
