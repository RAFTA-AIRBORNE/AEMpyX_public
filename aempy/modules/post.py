# -*- coding: utf-8 -*-
"""
Created on Apr 4, 2021
@author: vrath
"""
import sys
import os
from sys import exit as error
from datetime import datetime

import numpy
import scipy.linalg
import scipy.sparse 
import scipy.signal
import scipy.interpolate
from scipy.signal import medfilt, decimate
from scipy.ndimage import laplace, convolve
from scipy.ndimage import uniform_filter, gaussian_filter, median_filter
import differint 
import pylops

def fractrans(m=None, x=None , a=0.5):
    """
    Caklculate fractional derivative of m.

    VR Apr 2021
    """
    # import differint as df

    if m == None or x == None:
        error("No vector for diff given! Exit.")

    if numpy.size(m) != numpy.size(x):
        error("Vectors m and x have different length! Exit.")

    x0 = x[0]
    x1 = x[-1]
    npnts = numpy.size(x)
    mm = differint.differint(a, m, x0, x1, npnts)

    return mm

def crossgrad(m1=numpy.array([]),
                m2=numpy.array([]),
                mesh= [numpy.array([]), numpy.array([]), numpy.array([])],
                Out=True):
    """
    Crossgrad function
    
    
    See:
    Rosenkjaer GK, Gasperikova E, Newman, GA, Arnason K, Lindsey NJ (2015) 
        Comparison of 3D MT inversions for geothermal exploration: Case studies 
        for Krafla and Hengill geothermal systems in Iceland 
        Geothermics , Vol. 57, 258-274
        
    Schnaidt, S. (2015) 
        Improving Uncertainty Estimation in Geophysical Inversion Modelling 
        PhD thesis, University of Adelaide, AU
        
    vr  July 2023
    """
    sm = numpy.shape(m1)
    dm = m1.dim
    if dm==1:
        error("crossgrad: For dim="+str(dm)+" no crossgrad! Exit.")
    elif dm==2:
        cgdim = 1
    else:
        cgdim = 3

    gm1 = numpy.gradient(m1)
    gm2 = numpy.gradient(m2)
    sgm = numpy.shape(gm1)

    g1 = numpy.ravel(gm1)
    g2 = numpy.ravel(gm2)

    cgm = numpy.zeros_like(g1,cgdim)
    for k in numpy.arange(numpy.size(g1)):
        cgm[k,:] = numpy.cross (g1[k], g2[k])

    cgm =numpy.reshape(cgm,(sm+cgdim))

    
    cgnm = numpy.abs(cgm)/(numpy.abs(gm1)*numpy.abs(gm2))

    return cgm, cgnm

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
        ploton=False):
    """
    Apply 3D Anisotropic diffusion.

    Usage:
    stackout = anisodiff(stack, niter, kappa, gamma, option)

    Arguments:
            stack  - input stack
            niter  - number of iterations
            kappa  - conduction coefficient 20-100 ?
            gamma  - max value of .25 for stability
            step   - tuple, the distance between adjacent pixels in (z,y,x)
            option - 1 Perona Malik diffusion equation No 1
                     2 Perona Malik diffusion equation No 2
            ploton - if True, the middle z-plane will be plotted on every
                     iteration

    Returns:
            stackout   - diffused stack.

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

    Transcipy.linalgted to Python and optimised by Alistair Muldal
    Department of Pharmacology
    University of Oxford
    <alistair.muldal@pharm.ox.ac.uk>

    June 2000  original version.
    March 2002 corrected diffusion eqn No 2.
    July 2012 transcipy.linalgted to Python
    Jan 2021 slightly adapted for python3 VR
    """
    # initialize output array
    stackout = stack.copy()

    # initialize some internal variables
    deltaS = numpy.zeros_like(stackout)
    deltaE = numpy.zeros_like(stackout)
    deltaD = numpy.zeros_like(stackout)
    NS = numpy.zeros_like(stackout)
    EW = numpy.zeros_like(stackout)
    UD = numpy.zeros_like(stackout)
    gS = numpy.ones_like(stackout)
    gE = numpy.ones_like(stackout)
    gD = numpy.ones_like(stackout)

    # create the plot figure, if requested

    for ii in range(niter):

        # calculate the diffs
        deltaD[:-1, :, :] = numpy.diff(stackout, axis=0)
        deltaS[:, :-1, :] = numpy.diff(stackout, axis=1)
        deltaE[:, :, :-1] = numpy.diff(stackout, axis=2)

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
        stackout += gamma * (UD + NS + EW)

    return stackout


# def shock3d(
#         M,
#         dt=0.2, maxiter=30, filt=[3, 3, 3, 0.5],
#         boundary_mode="nearest", signfunc=None):
#     """
#     Apply shock filter in nD.

#     vr  Jan 2021
#     """
#     if signfunc is None or signfunc == "sign":
#         signcall = "-numpy.sign(L)"

#     elif signfunc[0] == "sigmoid":
#         scale = 1.0
#         signcall = "-1./(1. + numpy.exp(-scale *L))"

#     else:
#         error("sign func " + signfunc + " not defined! Exit.")

#     kersiz = (filt[0], filt[1], filt[2])
#     kerstd = filt[3]
#     K = gauss3D(kersiz, kerstd)
#     # print(numpy.sum(K.flat))
#     G = M

#     for it in range(maxiter):

#         G = convolve(G, K, mode=boundary_mode)

#         g = numpy.gradient(G)
#     #         print(numpy.shape(g))
#     #         normg=norm(g)
#     #         normg=numpy.sqrt(g[0])
#     #         print(numpy.shape(normg))
#     #         L = laplace(G)

#     #         S = eval(signcall)

#     #         G=G+dt*normg*S

#     return G


def gauss3D(Kshape=(3, 3, 3), Ksigma=0.5):
    """
    Define 3D gauscipy.signalan mask.

    Should give the same result as MATLAB's
    fspecial('gauscipy.signalan',[shape],[sigma])

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

def get_good_sites(q_val=None, q_thresh=None, out=True):
    """
    Clean models based on error or datafit

    """
    if (q_val is None) or (q_thresh is None):
        error("mod_qc: required input missing! Exit.")
        
        
    nsites= numpy.shape(q_val)[0]
    
    good = numpy.where(q_val<q_thresh)
    
    if out:
        print("number of good sites:", 
              numpy.count_nonzero(good),"from", nsites )
        
    return good
