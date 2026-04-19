# -*- coding: utf-8 -*-
'''
prep.py - Data preparation utilities for AEM data processing.

Provenance
----------
AEMpyX project.

@authors: Duygu Kiyan (DIAS), Volker Rath (DIAS)
With support of Claude (Anthropic, 2026)

Functions
---------
filter_column       : Apply Butterworth (or other) filter to selected columns.
insert_flag         : Replace bad/out-of-range values with a flag (NaN).
handle_gaps         : Impute or remove NaN-flagged rows/values.
calc_svd_decomp     : PCA/SVD-based filter after Minsley et al. (2012).
reduce_data         : Block-average a data array by column.
process_column      : Block-average a single column vector.

Author
------
vrath

Created
-------
Sat Nov 21 16:05:35 2020

Last change
-----------
vr  Apr 2026

'''

import sys

import numpy
import scipy.signal
import scipy.linalg
import scipy.interpolate

import aesys


def filter_column(
        M=None,
        Columns=None,
        Method=['butter', 4, 1./10],
        OutInfo=True):
    '''
    Apply a digital filter to selected columns of a data matrix.

    Parameters
    ----------
    M : numpy.ndarray, required
        Input data matrix (nsamples x ncolumns). Modified in-place.
    Columns : list of int [col_start, col_end], required
        Inclusive column range to filter.
    Method : list, optional
        Filter specification: [name, order, cutoff_frequency].
        Currently supported: ['butter', order, cutoff].
        The default is ['butter', 4, 1./10].
    OutInfo : bool, optional
        Print per-column progress to stdout. The default is True.

    Returns
    -------
    M : numpy.ndarray
        Data matrix with filtered columns.
    headstring : str
        Summary string describing the operation.

    Provenance
    ----------
    AEMpyX project.

    Last change: vr  Apr 2026

    '''
    cols = range(Columns[0], Columns[1] + 1)

    meth = Method[0]
    ford = Method[1]
    cutf = Method[2]
    headstring = 'filter columns ' + str(Columns) + ', method = ' + str(Method)

    if meth.lower() == 'butter':
        bb, ab = scipy.signal.butter(ford, cutf)
    else:
        sys.exit('Method ' + Method[0] + ' not yet implemented! Nothing done.')

    for icol in cols:
        tmp_in = M[:, icol]
        tmp_out = scipy.signal.filtfilt(
            bb,
            ab,
            tmp_in,
            method='gust',
            padtype='odd',
            padlen=24)

        M[:, icol] = tmp_out
        if OutInfo:
            comment = 'filter column ' + str(icol) + ', method = ' + str(Method)
            print(comment)

    return M, headstring


def insert_flag(
        M,
        Criterion='neg',
        ThreshVal=0.,
        Columns=[0, 0],
        Flag=None,
        InCol=None,
        System='aem05'):
    '''
    Replace bad or out-of-range values in selected columns with a flag.

    Parameters
    ----------
    M : numpy.ndarray, required
        Input data matrix. Modified in-place.
    Criterion : str, optional
        Flagging criterion. Valid options:
          'neg'   - flag values <= 0.
          'less'  - flag values <= ThreshVal.
          'great' - flag values >= ThreshVal.
          'plm'   - flag values where |value| > ThreshVal.
          'nan'   - propagate NaNs from column InCol to all Columns.
        The default is 'neg'.
    ThreshVal : float, optional
        Threshold value used with 'less', 'great', and 'plm'.
        The default is 0.
    Columns : list of int [col_start, col_end], optional
        Inclusive column range to check. The default is [0, 0].
    Flag : float, optional
        Replacement value for bad entries. The default is numpy.nan.
    InCol : int, optional
        Reference column for the 'nan' criterion. The default is None.
    System : str, optional
        AEM system identifier. The default is 'aem05'.

    Returns
    -------
    M : numpy.ndarray
        Data matrix with bad values replaced by Flag.
    NaNindex : tuple
        Indices (row, col) of all NaN entries in M after flagging.

    Provenance
    ----------
    AEMpyX project.

    Last change: vr  Jul 2021

    '''
    AEM_system = System
    _, nD, _, _, _ = aesys.get_system_params(System=AEM_system)

    if Flag is None:
        Flag = numpy.nan

    cols = numpy.arange(Columns[0], Columns[1] + 1)
    T = M[:, cols]
    oldflags = numpy.where(numpy.isnan(T))

    if 'neg' in Criterion.lower():
        T[oldflags] = -9999.
        T[T <= 0] = Flag
    elif 'less' in Criterion.lower():
        T[oldflags] = ThreshVal - 9999.
        T[T <= ThreshVal] = Flag
    elif 'great' in Criterion.lower():
        T[oldflags] = ThreshVal + 9999.
        T[T >= ThreshVal] = Flag
    elif 'plm' in Criterion.lower():
        T[oldflags] = ThreshVal + 9999.
        T[numpy.abs(T) > ThreshVal] = Flag
    elif 'nan' in Criterion.lower():
        if InCol is None:
            sys.exit('insert_flag: no reference column (InCol) given! Exit.')
        newindex = numpy.isnan(M[:, InCol])
        T[newindex, :] = Flag

    T[oldflags] = Flag
    M[:, cols] = T
    NaNindex = numpy.where(numpy.isnan(M))

    return M, NaNindex


def handle_gaps(M, Columns=[0, 0], Impute=['delete'], System='aem05'):
    '''
    Impute or remove NaN-flagged rows/values in selected columns.

    Valid imputation methods (Impute[0]):
      'nan'    - do nothing, keep NaNs.
      'delete' - remove all rows containing NaNs.
      'ave'    - replace NaNs with column means.
      'med'    - replace NaNs with column medians.
      'noi'    - replace NaNs with mean + Impute[1]*noise.
      'spl'    - spline interpolation along profile distance;
                 Impute[1] sets spline order.
                 'spln' adds zero-mean Gaussian noise of std
                 Impute[2] (or nanstd if Impute[2]==0).

    Parameters
    ----------
    M : numpy.ndarray, required
        Input data matrix. Modified in-place or rows deleted.
    Columns : list of int [col_start, col_end], optional
        Inclusive column range to process. The default is [0, 0].
    Impute : list, optional
        Imputation method and optional parameters. The default is ['delete'].
    System : str, optional
        AEM system identifier. The default is 'aem05'.

    Returns
    -------
    M : numpy.ndarray
        Data matrix after gap handling.

    Provenance
    ----------
    AEMpyX project.

    Last change: vr  Sep 2021

    '''
    AEM_system = System
    _, nD, _, _, _ = aesys.get_system_params(System=AEM_system)

    method = Impute.copy()

    if Columns != []:
        cols = range(Columns[0], Columns[1] + 1)
    else:
        method[0] = 'delete'

    nanindex = numpy.where(numpy.isnan(M))
    nanrows = numpy.ravel(nanindex[0])
    nancols = numpy.ravel(nanindex[1])

    if 'nan' in method[0].lower():
        print('Method: ' + method[0] + ' - nothing to do.')
        pass

    elif 'del' in method[0].lower():
        print('Method: ' + method[0])
        M = numpy.delete(M, nanrows, 0)

    elif 'ave' in method[0].lower():
        print('Method: ' + method[0])
        T = M[:, cols]
        val = numpy.nanmean(T, 0)
        for irow in nanrows:
            for icol in nancols:
                T[irow, icol] = val[icol]
        M[:, cols] = T

    elif 'med' in method[0].lower():
        print('Method: ' + method[0])
        T = M[:, cols]
        val = numpy.nanmedian(T, 0)
        for irow in nanrows:
            for icol in nancols:
                T[irow, icol] = val[icol]
        M[:, cols] = T

    elif 'noi' in method[0].lower():
        print('Method: ' + method[0])
        print('Noise of ' + str(method[1]) + ' ppm added.')
        T = M[:, cols]
        shapeT = numpy.shape(T)
        v_blk = numpy.nanmean(T, 0)
        v_std = Impute[1]
        val = v_blk + v_std * numpy.random.default_rng().normal(shapeT[0], shapeT[1])
        for irow in nanrows:
            for icol in nancols:
                T[irow, icol] = val[irow, icol]
        M[:, cols] = T

    elif 'spl' in method[0].lower():
        print('Method: ' + method[0])
        xd = M[:, 1]
        yd = M[:, 2]
        rd = numpy.sqrt(xd * xd + yd * yd)
        for var in range(cols):
            var = M[:, var]
            if not numpy.all(numpy.isnan(var)):
                spl = scipy.interpolate.UnivariateSpline(
                    xd, yd,
                    w=numpy.ones_like(xd),
                    k=Impute[1])
                T[:, var] = spl(rd)
                if 'spln' in method[0].lower():
                    if Impute[2] == 0.:
                        v0 = numpy.nanstd(var)
                    else:
                        v0 = Impute[2]
                    T[:, var] = (T[:, var]
                                 + v0 * numpy.random.standard_normal(numpy.shape(var)))
                    print('Noise of ' + str(v0) + ' ppm added.')
            else:
                print('Channel ' + str(var) + ' all nan! Nothing done.')
        M[:, cols] = T

    else:
        sys.exit('Method ' + method[0] + ' not yet implemented! Nothing done.')

    return M


def calc_svd_decomp(D=None, Columns=None, ErrD=None, K=1, ThreshMSE=0.,
                    OutFull=False, OutInfo=False):
    '''
    PCA/SVD filter after Minsley et al. (2012).

    Parameters
    ----------
    D : numpy.ndarray, required
        Input data matrix (nsamples x nfeatures). The default is None.
    Columns : list of int [col_start, col_end], required
        Inclusive column range defining the data block to decompose.
        The default is None.
    ErrD : numpy.ndarray, optional
        Per-column error (standard deviation) used to weight the
        decomposition. If None, an identity weighting is applied.
        The default is None.
    K : int, optional
        Number of principal components to retain. Must be <= len(Columns).
        The default is 1.
    ThreshMSE : float, optional
        MSE threshold (currently unused; reserved for future filtering).
        The default is 0.
    OutFull : bool, optional
        If True, return the full SVD factors in addition to the
        reconstructed matrix. The default is False.
    OutInfo : bool, optional
        Print relative Frobenius-norm reconstruction error.
        The default is False.

    Returns
    -------
    P : numpy.ndarray
        Reconstructed (filtered) data matrix.
    U, S, V : numpy.ndarray
        Left singular vectors, singular values, right singular vectors.
        Only returned when OutFull is True.
    MSE : float
        Root-mean-square error of the K-rank approximation.
        Only returned when OutFull is True.
    FRO : float
        Relative Frobenius-norm error ||X - X_k|| / ||X||.
        Only returned when OutFull is True.

    References
    ----------
    Minsley, B.J. et al. (2012). Airborne electromagnetic filtering ...

    Provenance
    ----------
    AEMpyX project.

    Last change: vr  Nov 2020

    '''
    if (D is None) or (Columns is None):
        sys.exit('calc_svd_decomp: no data defined! Exit.')

    cols = range(Columns[0], Columns[1] + 1)

    if ErrD is None:
        SqCovDi = SqCovD = numpy.diag(numpy.ones_like(cols))
    else:
        SqCovDi = numpy.diag(numpy.ones_like(cols) / ErrD)
        SqCovD  = numpy.diag(numpy.ones_like(cols) * ErrD)

    if K > numpy.size(cols):
        print('Number of SVs K=' + str(K) + ' > ' + str(numpy.size(cols)) + '. K assumed.')
        K = numpy.size(cols)

    P = D.copy()
    T = D[:, cols]
    print(numpy.shape(T))
    T_blk = numpy.nanmean(T, 0)
    X = numpy.dot(T - T_blk, SqCovDi)
    U, S, Vt = scipy.linalg.svd(X, full_matrices=False)
    SS = numpy.diag(S)
    V = Vt.T

    Tk  = numpy.dot(U[:, :K], numpy.dot(SS[:K, :K], V[:, :K].T))
    MSE = numpy.sqrt(numpy.mean((X - Tk) ** 2))
    FRO = numpy.linalg.norm(X - Tk, 'fro') / numpy.linalg.norm(X, 'fro')

    if OutInfo:
        print(' using ' + str(K) + ' PCs, |X-Xk|/|X| = %.6G' % (FRO))

    P[:, cols] = numpy.dot(Tk, SqCovD) + T_blk

    if OutFull:
        return P, U, S, V, MSE, FRO
    else:
        return P


def reduce_data(DataVec=numpy.array([]), System='aem05',
                Method=['mean', 5],
                ErrEst=False,
                OutInfo=True):
    '''
    Block-average a multi-column data array.

    Parameters
    ----------
    DataVec : numpy.ndarray, required
        Input data array (nsamples x ncolumns). The default is an empty array.
    System : str, optional
        AEM system identifier used to determine the active column count.
        The default is 'aem05'.
    Method : list, optional
        ['method_name', blocksize]. Blocksize is forced odd if even.
        Supported methods: 'mean', 'med' (median), 'dec' (decimate).
        The default is ['mean', 5].
    ErrEst : bool, optional
        If True, also return per-block error estimates. The default is False.
    OutInfo : bool, optional
        Print summary of output shapes. The default is True.

    Returns
    -------
    dat_out : numpy.ndarray
        Block-averaged data array.
    err_out : numpy.ndarray
        Block error estimates. Only returned when ErrEst is True.
    comment : str
        Summary string describing method and blocksize.

    Provenance
    ----------
    AEMpyX project.

    Last change: vr  Feb 2023

    '''
    if DataVec.size == 0:
        sys.exit('reduce_data: no data given! Exit.')

    dat_out = numpy.array([])
    err_out = numpy.array([])

    _, nD, _, _, _ = aesys.get_system_params(System)
    print(nD)

    meth, blocksize = Method

    if blocksize != 1:

        if blocksize % 2 == 0:
            blocksize = blocksize + 1
            Method[1] = blocksize

        for col in numpy.arange(nD[0]):

            if col == 0:
                dat, err, comment = process_column(DataVec[:, col],
                                                   Method=Method)
                dat_out = dat
            else:
                dat, err, comment = process_column(DataVec[:, col],
                                                   Method=Method,
                                                   ErrEst=ErrEst)
                dat_out = numpy.column_stack((dat_out, dat))
                if ErrEst:
                    if col == 1:
                        err_out = err
                    else:
                        err_out = numpy.column_stack((err_out, err))
    else:
        print('Blocksize is 1! Nothing done.')

    if ErrEst:
        err_out = numpy.nanmean(err_out, axis=1)

    comment = 'method = ' + meth + ', blocksize = ' + str(blocksize)

    if OutInfo:
        print(comment)
        print('New data shape: ', numpy.shape(dat_out))
        if ErrEst:
            print('Err shape: ', numpy.shape(err_out))

    if ErrEst:
        return dat_out, err_out, comment
    else:
        return dat_out, comment


def process_column(DataVec=numpy.array([]),
                   Method=['mean', 3],
                   ErrEst=False,
                   OutInfo=False):
    '''
    Block-average a single column (1-D) data vector.

    Parameters
    ----------
    DataVec : numpy.ndarray, required
        1-D input data vector. The default is an empty array.
    Method : list, optional
        ['method_name', blocksize]. Blocksize is forced odd if even.
        Supported methods: 'mean', 'med' (median), 'dec' (decimate).
        The default is ['mean', 3].
    ErrEst : bool, optional
        If True, compute per-block error (std for 'mean', MAD for 'med').
        Not available for 'dec'. The default is False.
    OutInfo : bool, optional
        Reserved for future verbose output. The default is False.

    Returns
    -------
    dat_out : numpy.ndarray
        Block-averaged values.
    dat_err : numpy.ndarray
        Block error estimates, or empty array if ErrEst is False.
    comment : str
        Summary string describing method and parameters.

    Provenance
    ----------
    AEMpyX project.

    Last change: vr  Nov 2020

    '''
    comment = 'method = ' + Method[0] + ' parameters = ' + str(Method[1:])

    if 'mean' in Method[0].lower():
        blocksize = Method[1]
        if blocksize % 2 == 0:
            blocksize = blocksize + 1
        exestr_avg = 'numpy.nanmean(tmp[i - blocksize:i])'
        if ErrEst:
            exestr_err = 'numpy.nanstd(tmp[i - blocksize:i])'

    elif 'med' in Method[0].lower():
        blocksize = Method[1]
        if blocksize % 2 == 0:
            blocksize = blocksize + 1
        exestr_avg = 'numpy.nanmedian(tmp[i - blocksize:i])'
        if ErrEst:
            exestr_err = 'inverse.calc_mad(tmp[i - blocksize:i], valavg)'

    elif 'dec' in Method[0].lower():
        qf = Method[1]
        ftype = 'iir'
        exestr_avg = ('scipy.signal.decimate(tmp,' + str(qf)
                      + ', ftype=' + ftype + ',zero_phase=True)')
        if ErrEst:
            sys.exit('process_column: error estimation not available '
                     'for method ' + Method[0].lower() + '! Exit.')
    else:
        sys.exit('process_column: method ' + Method[0].lower()
                 + ' not implemented! Exit.')

    data_head = numpy.flipud(DataVec[1:blocksize // 2 + 1])
    data_tail = numpy.flipud(DataVec[-blocksize // 2:len(DataVec) - 1])
    tmp = numpy.concatenate([data_head, DataVec, data_tail], axis=0)

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
    if ErrEst:
        dat_err = numpy.asarray(err_out)
    else:
        dat_err = numpy.array([])

    comment = 'method = ' + Method[0] + ' parameters = ' + str(Method[1:])

    return dat_out, dat_err, comment
