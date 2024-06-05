# -*- coding: utf-8 -*-
import os
import sys
from sys import exit as error
import numpy
from numpy.linalg import norm
from scipy.io import FortranFile
from scipy.ndimage import laplace, convolve, gaussian_gradient_magnitude
from scipy.ndimage import uniform_filter, gaussian_filter, median_filter


# import scipy.sparse as scp

import netCDF4 as nc
# import h5netcdf as hc
def decode_h2(strng):
    """
    Decode header2 string from ModEM Jacobian (old style).

    ----------
    strng : string
       header string

    Returns
    -------
    i1, i2, i3 : integer
        frequency, dattype, site numbers

    """

    s = strng.replace(";","").split()
    i1 = int(s[3])
    i2 = int(s[5])
    i3 = int(s[7])

    ivals = [i1, i2, i3]
    return ivals

def read_jac(JacFile=None, OutInfo=False):
    """
    Read Jacobian from ModEM OutInfoput.

    author: vrath
    last changed: Feb 10, 2021
    """
    if OutInfo:
        print("Opening and reading " + JacFile)

    fjac = FortranFile(JacFile, "r")
    tmp1 = []
    tmp2 = []

    _ = fjac.read_record(numpy.byte)
    # h1 = ''.join([chr(item) for item in header1])
    # print(h1)
    _ = fjac.read_ints(numpy.int32)
    # nAll = fjac.read_ints(numpy.int32)
    # print("nAll"+str(nAll))
    nTx = fjac.read_ints(numpy.int32)
    # print("ntx"+str(nTx))
    for i1 in range(nTx[0]):
        nDt = fjac.read_ints(numpy.int32)
        # print("nDt"+str(nDt))
        for i2 in range(nDt[0]):
            nSite = fjac.read_ints(numpy.int32)
            # print("nSite"+str(nSite))
            for i3 in range(nSite[0]):
                # header2
                header2 = fjac.read_record(numpy.byte)
                h2 = ''.join([chr(item) for item in header2])
                tmp2.append(decode_h2(h2))
                # print(h2)
                # print(i1,i2,i3)
                nSigma = fjac.read_ints(numpy.int32)
                # print("nSigma"+str(nSigma))
                for i4 in range(nSigma[0]):
                    # paramType
                    _ = fjac.read_ints(numpy.byte)
                    # p = ''.join([chr(item) for item in paramType])
                    # print(p)
                    # dims
                    _ = fjac.read_ints(numpy.int32)
                    # print(dims)
                    # dx
                    _ = fjac.read_reals(numpy.float64)
                    # dy
                    _ = fjac.read_reals(numpy.float64)
                    # dz
                    _ = fjac.read_reals(numpy.float64)
                    # AirCond
                    _ = fjac.read_reals(numpy.float64)
                    # ColJac = fjac.read_reals(numpy.float64).flatten(order="F")
                    ColJac = fjac.read_reals(numpy.float64).flatten()
                    # print(numpy.shape(CellSens))
                    # ColJac =  CellSens.flatten(order='F')
                    tmp1.append(ColJac)
                    # print(numpy.shape(tmp1))
                    # tmp2.append()
    Jac = numpy.asarray(tmp1)
    Inf = numpy.asarray(tmp2)
#    Inf = numpy.asarray(tmp2,dtype=object)

    fjac.close()

    if OutInfo:
        print("...done reading " + JacFile)

    return Jac, Inf  #, Site, Freq, Comp


def read_data_jac(DatFile=None, OutInfo=True):
    """
    Read ModEM input data.

    author: vrath
    last changed: Feb 10, 2021
    """
    Data = []
    Site = []
    Comp = []
    Head = []

    with open(DatFile) as fd:
        for line in fd:
            if line.startswith("#") or line.startswith(">"):
                Head.append(line)
                continue

            t = line.split()

            if "PT" in t[7] or "RH" in t[7] or "PH" in t[7]:
                tmp1 = [
                    float(t[0]),
                    float(t[2]),
                    float(t[3]),
                    float(t[4]),
                    float(t[5]),
                    float(t[6]),
                    float(t[8]),
                    float(t[9]),
                ]
                Data.append(tmp1)
                Site.append([t[1]])
                Comp.append([t[7]])
            else:
                tmp1 = [
                    float(t[0]),
                    float(t[2]),
                    float(t[3]),
                    float(t[4]),
                    float(t[5]),
                    float(t[6]),
                    float(t[8]),
                    float(t[10]),
                ]
                Data.append(tmp1)
                tmp2 = [
                    float(t[0]),
                    float(t[2]),
                    float(t[3]),
                    float(t[4]),
                    float(t[5]),
                    float(t[6]),
                    float(t[9]),
                    float(t[10]),
                ]
                Data.append(tmp2)
                Comp.append([t[7] + "R", t[7] + "I"])
                Site.append([t[1], t[1]])

    Site = [item for sublist in Site for item in sublist]
    Site = numpy.asarray(Site, dtype=object)
    Comp = [item for sublist in Comp for item in sublist]
    Comp = numpy.asarray(Comp, dtype=object)

    Data = numpy.asarray(Data)
    Freq = Data[:,0]

    nD = numpy.shape(Data)
    if OutInfo:
        print("readDat: %i data read from %s" % (nD[0], DatFile))

    return Data, Site, Freq, Comp, Head


def write_jac_ncd(NCFile=None, Jac=None, Dat=None, Site=None, Comp=None,
               zlib_in=True, shuffle_in=True, OutInfo=True):
    """
    Write Jacobian from ModEM OutInfoput to NETCDF/HDF5 file.

    author: vrath
    last changed: July 25, 2020
    """
    JacDim = numpy.shape(Jac)
    DatDim = numpy.shape(Dat)

    if JacDim[0] != DatDim[0]:
        print(
            "Error:  Jac dim="
            + str(JacDim[0])
            + " does not match Dat dim="
            + str(DatDim[0])
        )
        sys.exit(1)

    ncOutInfo = nc.Dataset(NCFile, "w", format="NETCDF4")
    ncOutInfo.createDimension("data", JacDim[0])
    ncOutInfo.createDimension("param", JacDim[1])

    S = ncOutInfo.createVariable(
        "site", str, ("data"), zlib=zlib_in, shuffle=shuffle_in)
    C = ncOutInfo.createVariable(
        "comp", str, ("data"), zlib=zlib_in, shuffle=shuffle_in)

    Per = ncOutInfo.createVariable(
        "Per", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    Lat = ncOutInfo.createVariable(
        "Lat", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    Lon = ncOutInfo.createVariable(
        "Lon", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    X = ncOutInfo.createVariable(
        "X", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    Y = ncOutInfo.createVariable(
        "Y", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    Z = ncOutInfo.createVariable(
        "Z", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    Val = ncOutInfo.createVariable(
        "Val", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)
    Err = ncOutInfo.createVariable(
        "Err", "float64", ("data"), zlib=zlib_in, shuffle=shuffle_in)

    J = ncOutInfo.createVariable(
        "Jac", "float64", ("data", "param"), zlib=zlib_in, shuffle=shuffle_in)

    S[:] = Site[
        :,
    ]
    C[:] = Comp[
        :,
    ]
    Per[:] = Dat[:, 0]
    Lat[:] = Dat[:, 1]
    Lon[:] = Dat[:, 2]
    X[:] = Dat[:, 3]
    Y[:] = Dat[:, 4]
    Z[:] = Dat[:, 5]
    Val[:] = Dat[:, 6]
    Err[:] = Dat[:, 7]
    J[:] = Jac

    ncOutInfo.close()

    if OutInfo:
        print(
            "writeJacNC: data written to %s in %s format" %
            (NCFile, ncOutInfo.data_model)
        )



def read_data(DatFile=None, OutInfo=True):
    """
    Read ModEM input data.

    author: vrath
    last changed: Feb 10, 2021
    """
    Data = []
    Site = []
    Comp = []
    Head = []

    with open(DatFile) as fd:
        for line in fd:
            if line.startswith("#") or line.startswith(">"):
                Head.append(line)

                continue

            t = line.split()

            if "PT" in t[7] or "RH" in t[7] or "PH" in t[7]:
                tmp = [
                    float(t[0]),
                    float(t[2]),
                    float(t[3]),
                    float(t[4]),
                    float(t[5]),
                    float(t[6]),
                    float(t[8]),
                    float(t[9]),
                    0.,
                ]
                Data.append(tmp)
                Site.append([t[1]])
                Comp.append([t[7]])
            else:
                tmp = [
                    float(t[0]),
                    float(t[2]),
                    float(t[3]),
                    float(t[4]),
                    float(t[5]),
                    float(t[6]),
                    float(t[8]),
                    float(t[9]),
                    float(t[10]),
                ]
                Data.append(tmp)
                Comp.append([t[7]])
                Site.append([t[1]])


    Site = [item for sublist in Site for item in sublist]
    Site = numpy.asarray(Site, dtype=object)
    Comp = [item for sublist in Comp for item in sublist]
    Comp = numpy.asarray(Comp, dtype=object)
    Data = numpy.asarray(Data)


    nD = numpy.shape(Data)
    if OutInfo:
        print("readDat: %i data read from %s" % (nD[0], DatFile))

    return Site, Comp, Data, Head


def write_data(DatFile=None, Dat=None, Site=None, Comp=None, Head = None,
               OutInfo=True):
    """
    Write ModEM input data file.

    author: vrath
    last changed: Feb 10, 2021
    """
    datablock =numpy.column_stack((Dat[:,0], Site[:], Dat[:,1:6], Comp[:], Dat[:,6:10]))
    nD, _ = numpy.shape(datablock)

    hlin = 0
    nhead = len(Head)
    nblck = int(nhead/8)
    print(str(nblck)+" blocks will be written.")

    with open(DatFile,"w") as fd:

        for ib in numpy.arange(nblck):
            blockheader = Head[hlin:hlin+8]
            hlin = hlin + 8
            for ii in numpy.arange(8):
                fd.write(blockheader[ii])

            if "Impedance" in blockheader[2]:

                fmt = "%14e %14s"+"%15.6f"*2+" %15.1f"*3+" %14s"+" %14e"*3

                indices = []
                block = []
                for ii in numpy.arange(len(Comp)):
                    if ("ZX" in Comp[ii]) or ("ZY" in Comp[ii]):
                        indices.append(ii)
                        block.append(datablock[ii,:])

                if OutInfo:
                    print('Impedances')
                    print(numpy.shape(block))

            elif "Vertical" in blockheader[2]:

                fmt = "%14e %14s"+"%15.6f"*2+" %15.1f"*3+" %14s"+" %14e"*3

                indices = []
                block = []
                for ii in numpy.arange(len(Comp)):
                    if ("TX" == Comp[ii]) or ("TY" == Comp[ii]):
                        indices.append(ii)
                        block.append(datablock[ii,:])

                if OutInfo:
                    print('Tipper')
                    print(numpy.shape(block))

            elif "Tensor" in blockheader[2]:

                fmt = "%14e %14s"+"%15.6f"*2+" %15.1f"*3+" %14s"+" %14e"*3

                indices = []
                block = []
                for ii in numpy.arange(len(Comp)):
                    if ("PT" in Comp[ii]):
                        indices.append(ii)
                        block.append(datablock[ii,:])

                if OutInfo:
                    print('Phase Tensor')
                    print(numpy.shape(block))

            else:
                error("Data type "+blockheader[3]+'not implemented! Exit.')

            numpy.savetxt(fd,block, fmt = fmt)


def write_data_ncd(
        NCFile=None, Dat=None, Site=None, Comp=None,
        zlib_in=True, shuffle_in=True, OutInfo=True
        ):
    """
    Write Jacobian from ModEM OutInfoput to NETCDF file.

    author: vrath
    last changed: July 24, 2020
    """
    try:
        NCFile.close
    except BaseException:
        pass

    DatDim = numpy.shape(Dat)

    ncOutInfo = nc.Dataset(NCFile, "w", format="NETCDF4")
    ncOutInfo.createDimension("data", DatDim[0])

    S = ncOutInfo.createVariable(
        "site", str, ("data",), zlib=zlib_in, shuffle=shuffle_in)
    C = ncOutInfo.createVariable(
        "comp", str, ("data",), zlib=zlib_in, shuffle=shuffle_in)

    Per = ncOutInfo.createVariable(
        "Per", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    Lat = ncOutInfo.createVariable(
        "Lat", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    Lon = ncOutInfo.createVariable(
        "Lon", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    X = ncOutInfo.createVariable(
        "X", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    Y = ncOutInfo.createVariable(
        "Y", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    Z = ncOutInfo.createVariable(
        "Z", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    Val = ncOutInfo.createVariable(
        "Val", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )
    Err = ncOutInfo.createVariable(
        "Err", "float64", ("data",), zlib=zlib_in, shuffle=shuffle_in
    )

    S[:] = Site[
        :,
    ]
    C[:] = Comp[
        :,
    ]
    Per[:] = Dat[:, 0]
    Lat[:] = Dat[:, 1]
    Lon[:] = Dat[:, 2]
    X[:] = Dat[:, 3]
    Y[:] = Dat[:, 4]
    Z[:] = Dat[:, 5]
    Val[:] = Dat[:, 6]
    Err[:] = Dat[:, 7]

    ncOutInfo.close()

    if OutInfo:
        print(
            "writeDatNC: data written to %s in %s format"
            % (NCFile, ncOutInfo.data_model)
        )


def write_model_ncd(
    NCFile=None,
    x=None,
    y=None,
    z=None,
    Mod=None,
    Sens=None,
    Ref=None,
    trans="LINEAR",
    zlib_in=True,
    shuffle_in=True,
    OutInfo=True,
):
    """
    Write Model from ModEM OutInfoput to NETCDF/HDF5 file.

    author: vrath
    last changed: Jan 21, 2021
    """
    ModDim = numpy.shape(Mod)

    ncOutInfo = nc.Dataset(NCFile, "w", format="NETCDF4")

    ncOutInfo.createDimension("msiz", ModDim)
    ncOutInfo.createDimension("nx", ModDim[0])
    ncOutInfo.createDimension("ny", ModDim[1])
    ncOutInfo.createDimension("nz", ModDim[2])

    ncOutInfo.createDimension("ref", (3))

    X = ncOutInfo.createVariable(
        "x", "float64", ("nx"), zlib=zlib_in, shuffle=shuffle_in)
    Y = ncOutInfo.createVariable(
        "y", "float64", ("ny"), zlib=zlib_in, shuffle=shuffle_in)
    Z = ncOutInfo.createVariable(
        "z", "float64", ("nz"), zlib=zlib_in, shuffle=shuffle_in)
    X[:] = x[:]
    Y[:] = y[:]
    Z[:] = z[:]

    trans = trans.upper()

    if trans == "LOGE":
        Mod = numpy.log(Mod)
        if OutInfo:
            print("resistivities to " + NCFile + " transformed to: " + trans)
    elif trans == "LOG10":
        Mod = numpy.log10(Mod)
        if OutInfo:
            print("resistivities to " + NCFile + " transformed to: " + trans)
    elif trans == "LINEAR":
        pass
    else:
        print("Transformation: " + trans + " not defined!")
        sys.exit(1)

    M = ncOutInfo.createVariable(
        "model", "float64", ("msiz"), zlib=zlib_in, shuffle=shuffle_in
    )
    M[:, :, :] = Mod[:, :, :]

    if Sens is not None:
        S = ncOutInfo.createVariable(
            "sens", "float64", ("msiz"), zlib=zlib_in, shuffle=shuffle_in
        )
        S[:, :, :] = Sens[:, :, :]

    if Ref is not None:
        R = ncOutInfo.createVariable(
            "ref", "float64", ("ref"), zlib=zlib_in, shuffle=shuffle_in
        )
        R[:, :, :] = Ref[:, :, :]

    ncOutInfo.close()

    if OutInfo:
        print(
            "writeModNC: data written to %s in %s format"
            % (NCFile, ncOutInfo.data_model)
        )


# def write_model_vtk(ModFile=None, dx=None, dy=None, dz=None, rho=None, reference=None,
#                  OutInfo=True):
#     """
#     write ModEM model input in .

#     Expects rho in physical units

#     author: vrath
#     last changed: Mar 13, 2021

#     """
#     dims = numpy.shape(rho)
#     nx = dims[0]
#     ny = dims[1]
#     nz = dims[2]

#     with open(ModFile, "w") as f:
#         numpy.savetxt(
#             f, [" # 3D MT model written by ModEM in WS format"], fmt="%s")
#         # line = numpy.array([nx, ny,nz, dummy, trans],dtype=('i8,i8,i8,i8,U10'))
#         numpy.savetxt(f, dx.reshape(1, dx.shape[0]), fmt="%12.3f")
#         numpy.savetxt(f, dy.reshape(1, dy.shape[0]), fmt="%12.3f")
#         numpy.savetxt(f, dz.reshape(1, dz.shape[0]), fmt="%12.3f")
#         # write OutInfo the layers from resmodel
#         for zi in range(dz.size):
#             f.write("\n")
#             for yi in range(dy.size):
#                 # line = rho[::-1, yi, zi]
#                 # line = numpy.flipud(rho[:, yi, zi])
#                 line = rho[:, yi, zi]
#                 numpy.savetxt(f, line.reshape(1, nx), fmt="%12.5e")

#         f.write("\n")

#         cnt = numpy.asarray(reference)
#         numpy.savetxt(f, cnt.reshape(1, cnt.shape[0]), fmt="%10.1f")
#         f.write("%10.2f  \n" % (0.0))


def write_model(ModFile=None, dx=None, dy=None, dz=None, mval=None, reference=None,
                trans=None, aircells = None, mvalair = 1.e17, blank = 1.e17, OutInfo=True):
    """
    Write ModEM model input.

    Expects mval in physical units (linear).

    author: vrath
    last changed: Aug 18, 2020


    In Fortran:

    DO iz = 1,Nz
        DO iy = 1,Ny
            DO ix = Nx,1,-1
                READ(10,*) mval(ix,iy,iz)
            ENDDO
        ENDDO
    ENDDO

    """
    dims = numpy.shape(mval)

    nx = dims[0]
    ny = dims[1]
    nz = dims[2]
    dummy = 0

    if trans is not None:

        trans = trans.upper()

        if trans == "LOGE":
            mval = numpy.log(mval)
            mvalair = numpy.log(mvalair)
            if OutInfo:
                print("values to " + ModFile + " transformed to: " + trans)
        elif trans == "LOG10":
            mval = numpy.log10(mval)
            mvalair = numpy.log10(mvalair)
            if OutInfo:
                print("values to " + ModFile + " transformed to: " + trans)
        elif trans == "LINEAR":
            pass

        else:
            print("Transformation: " + trans + " not defined!")
            sys.exit(1)


    else:
        trans = "LINEAR"

    if not aircells == None:
        mval.reshape(dims)[aircells] = mvalair

    if not blank == None:
        blanks = numpy.where(~numpy.isfinite(mval))
        mval.reshape(dims)[blanks] = mvalair

    trans = numpy.array(trans)
    with open(ModFile, "w") as f:
        numpy.savetxt(
            f, ["# 3D MT model written by ModEM in WS format"], fmt="%s")
        line = numpy.array([nx, ny,nz, dummy, trans],dtype="object")
        # line = numpy.array([nx, ny, nz, dummy, trans])
        # numpy.savetxt(f, line.reshape(1, 5), fmt="   %s"*5)
        numpy.savetxt(f, line.reshape(1, 5), fmt =["  %i","  %i","  %i","  %i", "  %s"])

        numpy.savetxt(f, dx.reshape(1, dx.shape[0]), fmt="%12.3f")
        numpy.savetxt(f, dy.reshape(1, dy.shape[0]), fmt="%12.3f")
        numpy.savetxt(f, dz.reshape(1, dz.shape[0]), fmt="%12.3f")
        # write OutInfo the layers from resmodel
        for zi in range(dz.size):
            f.write("\n")
            for yi in range(dy.size):
                line = mval[::-1, yi, zi]
                # line = numpy.flipud(mval[:, yi, zi])
                # line = mval[:, yi, zi]
                numpy.savetxt(f, line.reshape(1, nx), fmt="%12.5e")

        f.write("\n")

        cnt = numpy.asarray(reference)
        numpy.savetxt(f, cnt.reshape(1, cnt.shape[0]), fmt="%10.1f")
        f.write("%10.2f  \n" % (0.0))


def read_model(ModFile=None, trans="LINEAR", volumes=False, OutInfo=True):
    """
    Read ModEM model input.

    Returns mval in physical units

    author: vrath
    last changed: Aug 18, 2020

    In Fortran:

    DO iz = 1,Nz
        DO iy = 1,Ny
            DO ix = Nx,1,-1
                READ(10,*) mval(ix,iy,iz)
            ENDDO
        ENDDO
    ENDDO

    """
    with open(ModFile, "r") as f:
        lines = f.readlines()

    lines = [line.split() for line in lines]
    dims = [int(sub) for sub in lines[1][0:3]]
    nx, ny, nz = dims
    trns = lines[1][4]
    dx = numpy.array([float(sub) for sub in lines[2]])
    dy = numpy.array([float(sub) for sub in lines[3]])
    dz = numpy.array([float(sub) for sub in lines[4]])

    mval = numpy.array([])
    for line in lines[5:-2]:
        line = numpy.flipud(line) #  numpy.fliplr(line)
        mval = numpy.append(mval, numpy.array([float(sub) for sub in line]))

    if OutInfo:
        print("values in " + ModFile + " are: " + trns)
    if trns == "LOGE":
        mval = numpy.exp(mval)
    elif trns == "LOG10":
        mval = numpy.power(10.0, mval)
    elif trns == "LINEAR":
        pass
    else:
        print("Transformation: " + trns + " not defined!")
        sys.exit(1)

    # here mval should be in physical units, not log...
    if "loge" in trans.lower() or "ln" in trans.lower():
        mval = numpy.log(mval)
        if OutInfo:
            print("values transformed to: " + trans)
    elif "log10" in trans.lower():
        mval = numpy.log10(mval)
        if OutInfo:
            print("values transformed to: " + trans)
    else:
        if OutInfo:
            print("values transformed to: " + trans)
        pass

    mval = mval.reshape(dims, order="F")

    reference = [float(sub) for sub in lines[-2][0:3]]

    if OutInfo:
        print(
            "readMod: %i x %i x %i model read from %s" % (nx, ny, nz, ModFile))

    if volumes:
        vcell = numpy.zeros_like(mval)
        for ii in numpy.arange(len(dx)):
            for jj in numpy.arange(len(dy)):
                for kk in numpy.arange(len(dz)):
                    vcell[ii,jj,kk] = dx[ii]*dy[jj]*dz[kk]

        if OutInfo:
            print(
                "readMod: %i x %i x %i cell volumes calculated" % (nx, ny, nz))

        return dx, dy, dz, mval, reference, trans, vcell



    else:
        return dx, dy, dz, mval, reference, trans


def linear_interpolation(p1, p2, x0):
    """
    Function that receives as arguments the coordinates of two points (x,y)
    and returns the linear interpolation of a y0 in a given x0 position. This is the
    equivalent to obtaining y0 = y1 + (y2 - y1)*((x0-x1)/(x2-x1)).
    Look into https://en.wikipedia.org/wiki/Linear_interpolation for more
    information.

    Parameters
    ----------
    p1     : tuple (floats)
        Tuple (x,y) of a first point in a line.
    p2     : tuple (floats)
        Tuple (x,y) of a second point in a line.
    x0     : float
        X coordinate on which you want to interpolate a y0.

    Return float (interpolated y0 value)
    """
    y0 = p1[1] + (p2[1] - p1[1]) * ((x0 - p1[0]) / (p2[0] - p1[0]))

    return y0


def clip_model(x, y, z, rho,
               pad=[0, 0, 0], centers=False, scale=[1., 1., 1.]):
    """
    Clip model to ROI.

    Parameters
    ----------
    x, y, z : float
        Node coordinates
    rho : float
        resistivity/sensitivity/diff values.
    pad : integer, optional
        padding in x/y/z. The default is [0, 0, 0].
    centers: bool, optional
        nodes or centers. The default is False (i.e. nodes).
    scale: float
        scling, e.g. to km (1E-3). The default is [1., 1.,1.].

    Returns
    -------
    xn, yn, zn, rhon

    """
    if numpy.size(scale) == 1:
        scale = [scale, scale, scale]

    p_x, p_y, p_z = pad
    s_x, s_y, s_z = scale

    xn = s_x * x[p_x:-p_x]
    yn = s_y * y[p_y:-p_y]
    zn = s_z * z[0:-p_z]
    rhon = rho[p_x:-p_x, p_y:-p_y, 0:-p_z]

    if centers:
        print("cells3d returning cell center coordinates.")
        xn = 0.5 * (xn[:-1] + xn[1:])
        yn = 0.5 * (yn[:-1] + yn[1:])
        zn = 0.5 * (zn[:-1] + zn[1:])

    return xn, yn, zn, rhon



def mt1dfwd(freq, sig, d, inmod="res", outdat="both"):
    """
    1D magnetotelluric forward modelling
    based on A. Pethik's script at www.digitalearthlab.com
    Last change vr Feb 28, 2023
    """

    mu0 = 4.E-7 * numpy.pi   		# Magnetic Permeability (H/m)

    sig = numpy.array(sig)
    freq = numpy.array(freq)
    d = numpy.array(d)

    if "cond" in inmod.lower():
        sig = numpy.array(sig)
    if "res" in inmod.lower():
        sig = 1. / numpy.array(sig)

    if sig.ndim > 1:
        error('IP not yet implemented')

    n = numpy.size(sig)

    Z = numpy.zeros_like(freq) + 1j * numpy.zeros_like(freq)
    w = numpy.zeros_like(freq)

    ifr = -1
    for f in freq:
        ifr = ifr + 1
        w[ifr] = 2. * numpy.pi * f
        imp = numpy.array(range(n)) + numpy.array(range(n)) * 1j

        # compute basement impedance
        imp[n - 1] = numpy.sqrt(1j * w[ifr] * mu0 / sig[n - 1])

        for layer in range(n - 2, -1, -1):
            sl = sig[layer]
            dl = d[layer]
            # 3. Compute apparent rho from top layer impedance
            # Step 2. Iterate from bottom layer to top(not the basement)
            #   Step 2.1 Calculate the intrinsic impedance of current layer
            dj = numpy.sqrt(1j * w[ifr] * mu0 * sl)
            wj = dj / sl
            #   Step 2.2 Calculate Exponential factor from intrinsic impedance
            ej = numpy.exp(-2 * dl * dj)

            #   Step 2.3 Calculate reflection coeficient using current layer
            #          intrinsic impedance and the below layer impedance
            impb = imp[layer + 1]
            rj = (wj - impb) / (wj + impb)
            re = rj * ej
            Zj = wj * ((1 - re) / (1 + re))
            imp[layer] = Zj

        Z[ifr] = imp[0]
        # print(Z[ifr])

    if "imp" in outdat.lower():
        return Z

    if "rho" in outdat.lower():
        absZ = numpy.abs(Z)
        rhoa = (absZ * absZ) / (mu0 * w)
        phase = numpy.rad2deg(numpy.arctan(Z.imag / Z.real))
        return rhoa, phase
    else:
        absZ = numpy.abs(Z)
        rhoa = (absZ * absZ) / (mu0 * w)
        phase = numpy.rad2deg(numpy.arctan(Z.imag / Z.real))
        return Z, rhoa, phase

def insert_body(
    dx=None,
    dy=None,
    dz=None,
    rho_in=None,
    body=None,
    pad=[0, 0, 0],
    smooth=None,
    scale=1.0,
    Out=True,
):
    """
    Insert 3d ellipsoid or box into given model.

    Created on Sun Jan 3 10:35:28 2021
    @author: vrath
    """
    xpad = pad[0]
    ypad = pad[1]
    zpad = pad[2]

    xc, yc, zc = cells3d(dx, dy, dz, otype='c')

    modcenter = [0.5 * numpy.sum(dx), 0.5 * numpy.sum(dy), 0.0]

    xc = xc - modcenter[0]
    yc = yc - modcenter[1]
    zc = zc - modcenter[2]

    nx = numpy.shape(xc)[0]
    ny = numpy.shape(yc)[0]
    nz = numpy.shape(zc)[0]

    rho_OutInfo = numpy.log10(rho_in)

    geom = body[0]
    action = body[1]
    rhoval = body[2]
    bcent = body[3:6]
    baxes = body[6:9]
    bangl = body[9:12]

    if action[0:3] == "rep":
        actstring = "rhoval"
    elif action[0:3] == "add":
        actstring = "rho_OutInfo[ii,jj,kk] + rhoval"
    else:
        error("Action" + action + " not implemented! Exit.")

    if Out:
        print(
            "Body type   : " + geom + ", " + action + " rho =",
            str(numpy.power(10.0, rhoval)) + " Ohm.m",
        )
        print("Body center : " + str(bcent))
        print("Body axes   : " + str(baxes))
        print("Body angles : " + str(bangl))
        print("Smoothed with " + smooth[0] + " filter")

    if geom[0:3] == "ell":

        for kk in numpy.arange(0, nz - zpad - 1):
            zpoint = zc[kk]
            for jj in numpy.arange(ypad + 1, ny - ypad - 1):
                ypoint = yc[jj]
                for ii in numpy.arange(xpad + 1, nx - xpad - 1):
                    xpoint = xc[ii]
                    position = [xpoint, ypoint, zpoint]
                    # if Out:
                    # print('position')
                    # print(position)
                    # print( bcent)
                    if in_ellipsoid(position, bcent, baxes, bangl):
                        rho_OutInfo[ii, jj, kk] = eval(actstring)
                        # if Out:
                        #     print("cell %i %i %i" % (ii, jj, kk))

    if geom[0:3] == "box":

        for kk in numpy.arange(0, nz - zpad - 1):
            zpoint = zc[kk]
            for jj in numpy.arange(ypad + 1, ny - ypad - 1):
                ypoint = yc[jj]
                for ii in numpy.arange(xpad + 1, nx - xpad - 1):
                    xpoint = xc[ii]
                    position = [xpoint, ypoint, zpoint]
                    # if Out:
                    # print('position')
                    # print(position)
                    # print( bcent)

                    if in_box(position, bcent, baxes, bangl):
                        rho_OutInfo[ii, jj, kk] = eval(actstring)
                        # if Out:
                        #     print("cell %i %i %i" % (ii, jj, kk))

    if smooth is not None:
        if smooth[0][0:3] == "uni":
            fsize = smooth[1]
            rho_OutInfo = uniform_filter(rho_OutInfo, fsize)

        elif smooth[0][0:3] == "gau":
            gstd = smooth[1]
            rho_OutInfo = gaussian_filter(rho_OutInfo, gstd)

        else:
            error("Smoothing filter  " + smooth[0] + " not implemented! Exit.")

    rho_OutInfo = numpy.power(10.0, rho_OutInfo)

    return rho_OutInfo


def cells3d(dx, dy, dz, center=False, reference=[0., 0., 0.]):
    """
    Define cell coordinates.

    dx, dy, dz in m,
    Created on Sat Jan 2 10:35:28 2021

    @author: vrath

    """
    x = numpy.append(0.0, numpy.cumsum(dx))
    y = numpy.append(0.0, numpy.cumsum(dy))
    z = numpy.append(0.0, numpy.cumsum(dz))

    x = x + reference[0]
    y = y + reference[1]
    z = z + reference[2]

    if center:
        print("cells3d returning cell center coordinates.")
        xc = 0.5 * (x[:-1] + x[1:])
        yc = 0.5 * (y[:-1] + y[1:])
        zc = 0.5 * (z[:-1] + z[1:])
        return xc, yc, zc

    else:
        print("cells3d returning node coordinates.")
        return x, y, z


def in_ellipsoid(
    point=None,
    cent=[0.0, 0.0, 0.0],
    axs=[1.0, 1.0, 1.0],
    ang=[0.0, 0.0, 0.0],
    find_inside=True,
):
    """
    Find points inside arbitrary box.

    Defined by the 3-vectors cent, axs, ang
    vr dec 2020

    """
    # subtract center
    p = numpy.array(point) - numpy.array(cent)
    # rotation matrices
    rz = rotz(ang[2])
    p = numpy.dot(rz, p)
    ry = roty(ang[1])
    p = numpy.dot(ry, p)
    rx = rotx(ang[0])
    p = numpy.dot(rx, p)
    # R = rz*ry*rx
    # p = R*p

    # position in ellipsoid coordinates

    p = p / axs

    t = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] < 1.0
    # print(p,t)
    if not find_inside:
        t = not t

    return t


def in_box(
    point=None,
    cent=[0.0, 0.0, 0.0],
    axs=[1.0, 1.0, 1.0],
    ang=[0.0, 0.0, 0.0],
    find_inside=True,
):
    """
    Find points inside arbitrary ellipsoid.

    Defined by the 3-vectors cent, axs, ang
    vr dec 2020

    """
    # subtract center
    p = numpy.array(point) - numpy.array(cent)
    # rotation matrices
    rz = rotz(ang[2])
    p = numpy.dot(rz, p)
    ry = roty(ang[1])
    p = numpy.dot(ry, p)
    rx = rotx(ang[0])
    p = numpy.dot(rx, p)
    # R = rz*ry*rx
    # p = R*p

    # position in ellipsoid coordinates

    p = p / axs

    t = (
        p[0] <= 1.0
        and p[0] >= -1.0
        and p[1] <= 1.0
        and p[1] >= -1.0
        and p[2] <= 1.0
        and p[2] >= -1.0
    )
    # print(p,t)

    if not find_inside:
        t = not t

    return t


def rotz(theta):
    """
    Calculate 3x3 rotation matriz for rotation around z axis.

    vr dec 2020
    """
    t = numpy.radians(theta)
    s = numpy.sin(t)
    c = numpy.cos(t)

    M = numpy.array([c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0]).reshape(3, 3)

    return M


def roty(theta):
    """
    Calculate 3x3 rotation matrix for rotationa around y axis.

    vr dec 2020
    """
    t = numpy.radians(theta)
    s = numpy.sin(t)
    c = numpy.cos(t)

    M = numpy.array([c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c]).reshape(3, 3)

    return M


def rotx(theta):
    """
    Calculate 3x3 rotation matriz for rotation around x axis.

    vr dec 2020
    """
    t = numpy.radians(theta)
    s = numpy.sin(t)
    c = numpy.cos(t)

    M = numpy.array([1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c]).reshape(3, 3)

    return M


def medfilt3D(
        M,
        kernel_size=[3, 3, 3], boundary_mode="nearest", maxiter=1, Out=True):
    """
    Run iterated median filter in nD.

    vr  Jan 2021
    """
    tmp = M.copy()
    for it in range(maxiter):
        if Out:
            print("iteration: " + str(it))
        tmp = median_filter(tmp, size=kernel_size, mode=boundary_mode)

    G = tmp.copy()

    return G


def anidiff3D(
        M,
        ckappa=50, dgamma=0.1, foption=1, maxiter=30, Out=True):
    """
    Apply anisotropic nonlinear diffusion in nD.

    vr  Jan 2021
    """
    tmp = M.copy()

    tmp = anisodiff3D(
        tmp,
        niter=maxiter,
        kappa=ckappa,
        gamma=dgamma,
        step=(1.0, 1.0, 1.0),
        option=foption)

    G = tmp.copy()

    return G


def anisodiff3D(
        stack,
        niter=1, kappa=50, gamma=0.1, step=(1.0, 1.0, 1.0), option=1,
        plotit=False):
    """
    Apply 3D Anisotropic diffusion.

    Usage:
    stackOutInfo = anisodiff(stack, niter, kappa, gamma, option)

    Arguments:
            stack  - input stack
            niter  - number of iterations
            kappa  - conduction coefficient 20-100 ?
            gamma  - max value of .25 for stability
            step   - tuple, the distance between adjacent pixels in (z,y,x)
            option - 1 Perona Malik diffusion equation No 1
                     2 Perona Malik diffusion equation No 2
            plot.n - if True, the middle z-plane will be plot.ed on every
                     iteration

    Returns:
            stackOutInfo   - diffused stack.

    kappa controls conduction as a function of gradient.  If kappa is low
    small intensity gradients are able to block conduction and hence diffusion
    across step edges.  A large value reduces the influence of intensity
    gradients on conduction.

    gamma controls speed of diffusion (you usually want it at a maximum of
    0.25)

    step is used to scale the gradients in case the spacing between adjacent
    pixels differs in the x,y and/or z axes

    Diffusion equation 1 favours high contrast edges over low contrast ones.
    Diffusion equation 2 favours wide regions over smaller ones.

    Reference:
    P. Perona and J. Malik.
    Scale-space and edge detection using ansotropic diffusion.
    IEEE Transactions on Pattern Analysis and Machine Intelligence,
    12(7):629-639, July 1990.

    Original MATLAB code by Peter Kovesi
    School of Computer Science & Software Engineering
    The University of Western Australia
    pk @ csse uwa edu au
    <http://www.csse.uwa.edu.au>

    Translated to Python and optimised by Alistair Muldal
    Department of Pharmacology
    University of Oxford
    <alistair.muldal@pharm.ox.ac.uk>

    June 2000  original version.
    March 2002 corrected diffusion eqn No 2.
    July 2012 translated to Python
    Jan 2021 slightly adapted python3 VR
    """
    # initialize OutInfoput array

    stackOutInfo = stack.copy()

    # initialize some internal variables
    deltaS = numpy.zeros_like(stackOutInfo)
    deltaE = deltaS.copy()
    deltaD = deltaS.copy()
    NS = deltaS.copy()
    EW = deltaS.copy()
    UD = deltaS.copy()
    gS = numpy.ones_like(stackOutInfo)
    gE = gS.copy()
    gD = gS.copy()

    # create the plot.figure, if requested
    if plotit:
        import pylab as pl
        from time import sleep

        showplane = stack.shape[0] // 2

        fig = pl.figure(figsize=(20, 5.5), num="Anisotropic diffusion")
        ax1, ax2 = fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)

        ax1.imshow(
            stack[showplane, ...].squeeze(),
            interpolation="nearest")
        ih = ax2.imshow(
            stackOutInfo[showplane, ...].squeeze(),
            interpolation="nearest", animated=True
        )
        ax1.set_title("Original stack (Z = %i)" % showplane)
        ax2.set_title("Iteration 0")

        fig.canvas.draw()

    for ii in range(niter):

        # calculate the diffs
        deltaD[:-1, :, :] = numpy.diff(stackOutInfo, axis=0)
        deltaS[:, :-1, :] = numpy.diff(stackOutInfo, axis=1)
        deltaE[:, :, :-1] = numpy.diff(stackOutInfo, axis=2)

        # conduction gradients (only need to compute one per dim!)
        if option == 1:
            gD = numpy.exp(-((deltaD / kappa) ** 2.0)) / step[0]
            gS = numpy.exp(-((deltaS / kappa) ** 2.0)) / step[1]
            gE = numpy.exp(-((deltaE / kappa) ** 2.0)) / step[2]
        elif option == 2:
            gD = 1.0 / (1.0 + (deltaD / kappa) ** 2.0) / step[0]
            gS = 1.0 / (1.0 + (deltaS / kappa) ** 2.0) / step[1]
            gE = 1.0 / (1.0 + (deltaE / kappa) ** 2.0) / step[2]

        # update matrices
        D = gD * deltaD
        E = gE * deltaE
        S = gS * deltaS

        # subtract a copy that has been shifted 'Up/North/West' by one
        # pixel. don't as questions. just do it. trust me.
        UD[:] = D
        NS[:] = S
        EW[:] = E
        UD[1:, :, :] -= D[:-1, :, :]
        NS[:, 1:, :] -= S[:, :-1, :]
        EW[:, :, 1:] -= E[:, :, :-1]

        # update the image
        stackOutInfo += gamma * (UD + NS + EW)

        if plotit:
            iterstring = "Iteration %i" % (ii + 1)
            ih.set_data(stackOutInfo[showplane, ...].squeeze())
            ax2.set_title(iterstring)
            fig.canvas.draw()
            # sleep(0.01)

    return stackOutInfo


def shock3d(
        M,
        dt=0.2, maxiter=30, filt=[3, 3, 3, 0.5],
        boundary_mode="nearest", signfunc=None):
    """
    Apply shock filter in nD.

    vr  Jan 2021
    """
    if signfunc is None or signfunc == "sign":
        signcall = "-numpy.sign(L)"

    elif signfunc[0] == "sigmoid":
        scale = 1.0
        signcall = "-1./(1. + numpy.exp(-scale *L))"

    else:
        error("sign func " + signfunc + " not defined! Exit.")

    kersiz = (filt[0], filt[1], filt[2])
    kerstd = filt[3]
    K = gauss3D(kersiz, kerstd)
    # print(numpy.sum(K.flat))
    G = M

    for it in range(maxiter):

        G = convolve(G, K, mode=boundary_mode)

        g = numpy.gradient(G)
    #         print(numpy.shape(g))
    #         normg=norm(g)
    #         normg=numpy.sqrt(g[0])
    #         print(numpy.shape(normg))
    #         L = laplace(G)

    #         S = eval(signcall)

    #         G=G+dt*normg*S

    return G


def gauss3D(Kshape=(3, 3, 3), Ksigma=0.5):
    """
    Define 2D gaussian mask.

    Should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])

    vr  Jan 2021
    """
    k, m, n = [(ss - 1) / 2 for ss in Kshape]
    x, y, z = numpy.ogrid[-n:n+1, -m:m+1, -k:k+1]
    h = numpy.exp(-(x * x + y * y + z * z) / (2.0 * Ksigma * Ksigma))
    h[h < numpy.finfo(h.dtype).eps * h.max()] = 0
    s = h.sum()
    if s != 0:
        h /= s

    K = h

    return K


def prepare_model(rho, rhoair=1.0e17):
    """
    Prepare model for filtering etc.

    Mainly redefining the boundaries (in the case of topograpy)
    Air domain is filed with vertical surface value
    Created on Tue Jan  5 11:59:42 2021

    @author: vrath
    """
    nn = numpy.shape(rho)

    rho_new = rho

    for ii in range(nn[0]):
        for jj in range(nn[1]):
            tmp = rho[ii, jj, :]
            na = numpy.argwhere(tmp < rhoair / 100.0)[0]
            # print(' orig')
            # print(tmp)
            tmp[: na[0]] = tmp[na[0]]
            # print(' prep')
            # print(tmp)
            rho_new[ii, jj, :] = tmp

    return rho_new


def calc_rhoa_phas(freq=None, Z=None):

    mu0 = 4.0e-7 * numpy.pi  # Magnetic Permeability (H/m)
    omega = 2.*numpy.pi*freq

    rhoa = numpy.power(numpy.abs(Z), 2) / (mu0 * omega)
    # phi = numpy.rad2deg(numpy.arctan(Z.imag / Z.real))
    phi = numpy.angle(Z, deg=True)

    return rhoa, phi

def mt1dfwd(freq, sig, d, inmod="r", out="imp", magfield="b"):
    """
    Calulate 1D magnetotelluric forward response.

    based on A. Pethik's script at www.digitalearthlab.com
    Last change vr Nov 20, 2020
    """
    mu0 = 4.0e-7 * numpy.pi  # Magnetic Permeability (H/m)

    sig = numpy.array(sig)
    freq = numpy.array(freq)
    d = numpy.array(d)

    if "c" in inmod[0].lower():
        sig = numpy.array(sig)
    else:
        sig = 1.0 / numpy.array(sig)

    if sig.ndim > 1:
        error("IP not yet implemented")

    n = numpy.size(sig)

    Z = numpy.zeros_like(freq) + 1j * numpy.zeros_like(freq)
    w = numpy.zeros_like(freq)

    ifr = -1
    for f in freq:
        ifr = ifr + 1
        w[ifr] = 2.0 * numpy.pi * f
        imp = numpy.array(range(n)) + numpy.array(range(n)) * 1j

        # compute basement impedance
        imp[n - 1] = numpy.sqrt(1j * w[ifr] * mu0 / sig[n - 1])

        for layer in range(n - 2, -1, -1):
            sl = sig[layer]
            dl = d[layer]
            # 3. Compute apparent rho from top layer impedance
            # Step 2. Iterate from bottom layer to top(not the basement)
            #   Step 2.1 Calculate the intrinsic impedance of current layer
            dj = numpy.sqrt(1j * w[ifr] * mu0 * sl)
            wj = dj / sl
            #   Step 2.2 Calculate Exponential factor from intrinsic impedance
            ej = numpy.exp(-2 * dl * dj)

            #   Step 2.3 Calculate reflection coeficient using current layer
            #          intrinsic impedance and the below layer impedance
            impb = imp[layer + 1]
            rj = (wj - impb) / (wj + impb)
            re = rj * ej
            Zj = wj * ((1 - re) / (1 + re))
            imp[layer] = Zj

        Z[ifr] = imp[0]
        # print(Z[ifr])

    if "imp" in out.lower():

        if "b" in magfield.lower():
            return Z/mu0
        else:
            return Z

    elif "rho" in out.lower():
        absZ = numpy.abs(Z)
        rhoa = (absZ * absZ) / (mu0 * w)
        phase = numpy.rad2deg(numpy.arctan(Z.imag / Z.real))

        return rhoa, phase
    else:
        absZ = numpy.abs(Z)
        rhoa = (absZ * absZ) / (mu0 * w)
        phase = numpy.rad2deg(numpy.arctan(Z.imag / Z.real))
        return Z, rhoa, phase

