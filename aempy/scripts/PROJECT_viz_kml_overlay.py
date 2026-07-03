# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:light,ipynb
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (Spyder)
#     language: python3
#     name: python3
# ---

# +
#!/usr/bin/env python3
# -

# This script interpolates flight-line data onto a regular grid and
# produces georeferenced raster overlays (borderless PNG) plus the
# accompanying KML/KMZ GroundOverlay wrapper, for display in Google Earth.

# +
import os
import sys
import shutil
import zipfile

from time import process_time
from datetime import datetime
import getpass
import inspect

import numpy
import matplotlib
import matplotlib.pyplot

import scipy.interpolate
import scipy.spatial

import pyproj

AEMPYX_ROOT = os.environ['AEMPYX_ROOT']
mypath = [os.path.join(AEMPYX_ROOT, 'aempy/modules/')]
for pth in mypath:
    if pth not in sys.path:
        sys.path.insert(0,pth)

from version import versionstrg
import util
import aesys
import inverse

# +
OutInfo = True
cm = 1/2.54
AEMPYX_DATA = os.environ['AEMPYX_DATA']

version, _ = versionstrg()
# script = 'Tutorial1_VIZ_kml_overlay.py'
script = inspect.getfile(inspect.currentframe())  # this only works in python, not jupyter notebook
titstrng = util.print_title(version=version, fname=script, out=False)
print(titstrng+'\n\n')
Header = titstrng
# -

# The following cell gives values to AEM-system related settings.
#
# Data transformation is activated by the variable _DataTrans_. Currently
# three possible options are allowed: _DataTrans = 0_: No transformation,
# i.e., the raw data are used. _DataTrans = 1_: The natural log of data
# is taken, only allowed for strictly positive values. _DataTrans = 2_:
# If data scale logarithmically, an _asinh_ transformation (introduced by
# Scholl, 2000) is applied. It allows negatives, which may occur in TDEM,
# when IP effects are present.

# +
# AEM_system = 'genesis'
AEM_system = 'aem05'
if 'aem05' in AEM_system.lower():
    _, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 0
    DatErr_add =  50.
    DatErr_mult = 0.05
    data_active = numpy.ones(NN[2], dtype='int8')

    CompDict = Misc[3]
    CompLabl = list(CompDict.keys())
    print(CompLabl)

if 'genes' in AEM_system.lower():
    _, NN, _, _, Misc, = aesys.get_system_params(System=AEM_system)
    nL = NN[0]
    ParaTrans = 1
    DataTrans = 2
    DatErr_add = 100.
    DatErr_mult = 0.01
    data_active = numpy.ones(NN[2], dtype='int8')
    data_active[0:11]=0  # only vertical component
    CompDict =Misc[2]
    CompLabl = list(CompDict.keys())

# +

InFileFmt = '.npz'


##############################################################################
# Somaye A5
##############################################################################
# AEMPYX_DATA  = AEMPYX_ROOT+"/aempy/examples/A1_StGormans/"
AEMPYX_DATA  = "/home/vrath/work/Data_Somaye/"
FileList = "search"
# un/comment according to which data  you want to plot

# # processed data
SearchStrng = "*FL*k2.npz"# "search", "read"
InDatDir = AEMPYX_DATA+"/proc_delete_PLM3s/"
PlotDir = AEMPYX_DATA +"/plots/"
KmlDir = AEMPYX_DATA
PlotStrng = "proc_k2"


PlotName = "Somaye_"+PlotStrng


print('Data read from dir: %s ' % InDatDir)
print('Overlays written to dir: %s ' % KmlDir)
print('Overlay basename: %s ' % PlotName)

ListName = ''
# -

if 'set' in FileList.lower():
    dat_files = []

if 'read' in FileList.lower():
    print('File names read from : '+ListName)
    how = ['read', ListName, InDatDir]
    dat_files = util.get_data_list(how=how,
                              out= True, sort=True)

if 'search' in FileList.lower():
    print('Searchstring is : '+SearchStrng)
    how = ['search', SearchStrng, InDatDir]
    dat_files = util.get_data_list(how=how, fullpath=True,
                              out= True, sort=True)

ns = numpy.size(dat_files)
if ns ==0:
    sys.exit('No files set!. Exit.')


MergeData = True
DataMergeFile = InDatDir+PlotName+'_merged.npz'

# The source coordinate system of Easting/Northing in the data files.
# _EPSGCode_ must match the projected CRS used when the flight-line data
# were processed (e.g. Irish Transverse Mercator = 'epsg:2157',
# UTM zone 29N = 'epsg:32629'). The overlay is always reprojected to
# WGS84 lon/lat ('epsg:4326'), as required by the KML LatLonBox.

EPSGCode = 'epsg:32629'
WGS84Code = 'epsg:4326'
Transformer = pyproj.Transformer.from_crs(EPSGCode, WGS84Code, always_xy=True)

XYUnits = '(m)'
XYFact = 1.0   # keep native metres; only used internally for the mesh

#
# Kernel functions for RBF:
#     The radial basis function, based on the radius, r,
#     given by the norm (default is Euclidean distance); the default is 'multiquadric':
#         'linear' : -r
#         'thin_plate_spline' : r**2 * log(r)
#         'cubic' : r**3
#         'quintic' : -r**5
#
# Methods for griddata:
#         'nearest'       data point closest to the point of interpolation
#         'linear'        tessellate the input point set to N-D simplices
#                         and interpolate linearly on each simplex
#         'cubic'         piecewise cubic, continuously differentiable (C1)
#

step = 1

InterpMethod = ['griddata','linear']
# InterpMethod = ['griddata', 'cubic']
# InterpMethod = ['rbf', 'linear', 0.0]
# InterpMethod = ['rbf', 'thin_plate_spline', 0.0]

numIndexes = [301, 351]
smooth = 0.
MaskNeg = True
MaskPoly = False
MaskDist = True

if MaskDist:
    DistMask = 100.*XYFact

if MaskPoly:
    PolyDir = AEMPYX_DATA+'/Blocks/polygons/'
    PolyFiles = [PolyDir+'A5_2019_utm.npz']
    Polygon= numpy.load(PolyFiles[0], allow_pickle=True)['Poly'][0]

# The following cell determines the settings for individual components.
# Each sublist associated to a component contains the name, followed by
# a list of parameters determining the data limits ([vmin, vmax, step]),
# and, for some components, an additional threshold value.

# +
CompList=[
    ['P1', []],
    ['Q1', []],
    ['P2', []],
    ['Q2', []],
    ['P3', []],
    ['Q3', []],
    ['P4', []],
    ['Q4', []],
    ['PLM', [], 0.5],      # PLMthresh = 0.5
    ['ALT', [40., 120., 20.], 300.]   # ALTthresh = 300
]
# -

# Below, some graphic parameters are set. The overlay image itself is
# always rendered without axes, frame, ticks, or margin, since it must
# align exactly with the geographic bounding box declared in the KML.

matplotlib.pyplot.style.use('seaborn-v0_8-paper')
matplotlib.rcParams['figure.dpi'] = 400
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['savefig.facecolor'] = 'none'

OverlayDpi = 300
OverlayWidthInches = 8.

#
# Determine colormap. A comprehensive list of colormaps can be found at
# https://matplotlib.org/stable/gallery/color/colormap_reference.html
#

# +
Cmap ='viridis'
# Cmap = 'jet_r'
cmp = matplotlib.colormaps[Cmap]

KmzBundle = True    # bundle all component overlays into a single .kmz
KmzName = PlotName+'_overlays.kmz'
# -

if not os.path.isdir(KmlDir):
    print('File: %s does not exist, but will be created' % KmlDir)
    os.mkdir(KmlDir)


def make_overlay_png(XI, YI, DI, cmap, vmin, vmax, outfile,
                      dpi=OverlayDpi, width_in=OverlayWidthInches):
    """
    Render an interpolated grid as a borderless, transparent PNG whose
    canvas exactly matches the [XI, YI] bounding box, for use as a KML
    GroundOverlay Icon.

    Parameters
    ----------
    XI, YI : 2D arrays
        Mesh coordinates (native projected CRS, e.g. metres).
    DI : 2D array
        Interpolated data values, same shape as XI, YI.
    cmap : matplotlib colormap
    vmin, vmax : float or None
        Colour limits.
    outfile : str
        Output PNG path.

    Returns
    -------
    x_min, x_max, y_min, y_max : float
        Bounding box in the native projected CRS (metres).

    Provenance
    ----------
    AEMpyX project. @authors: Duygu Kiyan (DIAS), Volker Rath (DIAS).
    With support of Claude (Anthropic, 2026)

    Last change
    -----------
    2026-07
    """
    x_min, x_max = numpy.amin(XI), numpy.amax(XI)
    y_min, y_max = numpy.amin(YI), numpy.amax(YI)

    aspect = (y_max-y_min)/(x_max-x_min)
    fig = matplotlib.pyplot.figure(figsize=(width_in, width_in*aspect), dpi=dpi)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.axis('off')

    ax.pcolormesh(XI, YI, DI, cmap=cmap, vmin=vmin, vmax=vmax,
                  shading='auto')

    fig.savefig(outfile, dpi=dpi, transparent=True, pad_inches=0)
    matplotlib.pyplot.close(fig)

    return x_min, x_max, y_min, y_max


def make_lookat(lon_min, lon_max, lat_min, lat_max):
    """
    Build a KML <LookAt> element centred on a lon/lat bounding box, so
    Google Earth flies the camera to the overlay on load instead of
    leaving it un-navigated in the Places panel.

    Provenance
    ----------
    AEMpyX project. @authors: Duygu Kiyan (DIAS), Volker Rath (DIAS).
    With support of Claude (Anthropic, 2026)

    Last change
    -----------
    2026-07
    """
    lon_c = 0.5*(lon_min+lon_max)
    lat_c = 0.5*(lat_min+lat_max)
    # crude range estimate: degrees -> metres -> viewing range
    span_deg = max(lon_max-lon_min, lat_max-lat_min, 1e-4)
    rng = max(span_deg*150000., 500.)
    return f"""  <LookAt>
    <longitude>{lon_c:.6f}</longitude>
    <latitude>{lat_c:.6f}</latitude>
    <altitude>0</altitude>
    <range>{rng:.0f}</range>
    <tilt>0</tilt>
    <heading>0</heading>
    <altitudeMode>clampToGround</altitudeMode>
  </LookAt>
"""


def write_kml_groundoverlay(png_name, lon_min, lon_max, lat_min, lat_max,
                             name, description=''):
    """
    Build a single KML <GroundOverlay> element (as a string).

    Parameters
    ----------
    png_name : str
        Icon href, relative to the KML/KMZ root.
    lon_min, lon_max, lat_min, lat_max : float
        Bounding box in WGS84 lon/lat.
    name : str
        Overlay display name (shown in Google Earth's layer list).
    description : str, optional

    Returns
    -------
    kml_fragment : str

    Provenance
    ----------
    AEMpyX project. @authors: Duygu Kiyan (DIAS), Volker Rath (DIAS).
    With support of Claude (Anthropic, 2026)

    Last change
    -----------
    2026-07
    """
    return f"""  <GroundOverlay>
    <name>{name}</name>
    <description>{description}</description>
    <Icon>
      <href>{png_name}</href>
    </Icon>
    <LatLonBox>
      <north>{lat_max:.6f}</north>
      <south>{lat_min:.6f}</south>
      <east>{lon_max:.6f}</east>
      <west>{lon_min:.6f}</west>
      <rotation>0</rotation>
    </LatLonBox>
  </GroundOverlay>
"""


if MergeData:
    Data = inverse.merge_data_sets(infile_list=dat_files,
                                outfile_name=DataMergeFile,
                                aem_system='aem05', out= True)
    dat_files = [DataMergeFile]

kml_folders = []
png_files = []
agg_lon_min, agg_lon_max = numpy.inf, -numpy.inf
agg_lat_min, agg_lat_max = numpy.inf, -numpy.inf

for filein in dat_files:
    start = process_time()
    print('\nData read from: %s' % filein)
    Data, header, _ = aesys.read_aempy(File=filein, System=AEM_system, OutInfo=False)

    E = Data[:,1][::step]*XYFact
    E_min = numpy.amin(E)
    E_max = numpy.amax(E)
    N = Data[:,2][::step]*XYFact
    N_min = numpy.amin(N)
    N_max = numpy.amax(N)

    xi = numpy.linspace(E_min, E_max, numIndexes[0])
    yi = numpy.linspace(N_min, N_max, numIndexes[1])
    dx = numpy.around(numpy.diff(xi)[0], decimals=0)
    dy = numpy.around(numpy.diff(yi)[0], decimals=0)
    print('Interpolation mesh, dx = '+ str(dx)+' m, dy ='+ str(dy)+' m')
    XI, YI = numpy.meshgrid(xi, yi, indexing='ij')
    Pnts = numpy.stack([ E.ravel(),  N.ravel()], -1)
    Mesh = numpy.stack([XI.ravel(), YI.ravel()], -1)

    if MaskDist:
        D_tree=scipy.spatial.KDTree(Pnts, leafsize=10,
                                    compact_nodes=True,
                                    copy_data=True,
                                    balanced_tree=True,
                                    boxsize=None)
        mindist, _ = D_tree.query(Mesh, k=1)
        blankdist = mindist>=DistMask

    if MaskPoly:
        blankpoly=[]
        for ipnt in numpy.arange(numpy.size(Mesh[:,0])):
            outside = not util.point_inside_polygon(Mesh[ipnt,0], Mesh[ipnt,1],
                                                    Polygon)
            blankpoly.append(outside)

    for nc in numpy.arange(len(CompList)):

        Comp = CompList[nc][0]
        comp = CompDict[Comp][0]
        indx = CompLabl.index(Comp)

        titl = CompLabl[indx]+CompDict[Comp][2]+': '\
            +PlotStrng+'/'+str(DataTrans)\
            +'/'+InterpMethod[0]+'/'+InterpMethod[1]

        print('Building overlay for component '+titl)
        D = Data[:,comp][::step]

        if ('Z' in Comp) or ('H' in Comp) or ('P' in Comp) or ('Q' in Comp):
            Unit = 'ppm'
            if DataTrans ==1:
                Unit = '-'
                D = numpy.log10(D)
            if DataTrans ==2:
                S = 100.
                Unit = '-'
                if not numpy.isfinite(S):
                   S = inverse.get_S(D)
                D= numpy.arcsinh(D/S)

        if ('PL' in Comp):
            Unit = '(-)'
            PLMthresh= CompList[nc][2]

        if ('A' in Comp):
            Unit = 'm'
            ALTthresh= CompList[nc][2]

        if ('R' in Comp):
            Unit = 'Ohm.m'

        Dats = D.flatten()
        if 'grid' in InterpMethod[0].lower():
            DI = scipy.interpolate.griddata(Pnts, Dats, Mesh,
                                            method=InterpMethod[1].lower())
            DI = numpy.reshape(DI,(len(xi), len(yi)))

        elif 'rbf' in InterpMethod[0].lower():
            RBF = scipy.interpolate.RBFInterpolator(
                        Pnts, Dats,
                        kernel=InterpMethod[1], smoothing=InterpMethod[2])
            DI = RBF(Mesh)
            DI = numpy.reshape(DI,(len(xi), len(yi)))

        if ('PL' in Comp):
            DI[DI<=PLMthresh]=numpy.nan
        if ('A' in Comp):
            DI[DI>=ALTthresh]=numpy.nan
        if MaskNeg:
            DI[DI<=0.]= numpy.nan
        if MaskPoly:
            DIF = DI.flatten().reshape(-1,1)
            DIF[blankpoly] = numpy.nan
            DI = numpy.reshape(DIF,(len(xi), len(yi)))
        if MaskDist:
            DIF = DI.flatten().reshape(-1,1)
            DIF[blankdist] = numpy.nan
            DI = numpy.reshape(DIF,(len(xi), len(yi)))

        D_min = numpy.nanmin(DI)
        D_max = numpy.nanmax(DI)
        print('Data, interpolated   min='+str( D_min)+'   max='+str( D_max))

        if len(CompList[nc][1])==0:
            vmin, vmax = None, None
        else:
            vmin, vmax, _ = CompList[nc][1]

        pngfile = KmlDir+PlotName+'_'+AEM_system+'_kmloverlay_'+Comp+'.png'
        x_min, x_max, y_min, y_max = make_overlay_png(
            XI, YI, DI, cmp, vmin, vmax, pngfile)

        lon_min, lat_min = Transformer.transform(x_min, y_min)
        lon_max, lat_max = Transformer.transform(x_max, y_max)
        print('  bbox lon/lat: '+str(round(lon_min,4))+', '+str(round(lat_min,4))
              +'  to  '+str(round(lon_max,4))+', '+str(round(lat_max,4)))

        agg_lon_min = min(agg_lon_min, lon_min)
        agg_lon_max = max(agg_lon_max, lon_max)
        agg_lat_min = min(agg_lat_min, lat_min)
        agg_lat_max = max(agg_lat_max, lat_max)

        kmlfile = KmlDir+PlotName+'_'+AEM_system+'_kmloverlay_'+Comp+'.kml'
        overlay_name = AEM_system.upper()+' '+Comp+' ('+Unit+')'
        fragment = write_kml_groundoverlay(
            os.path.basename(pngfile), lon_min, lon_max, lat_min, lat_max,
            name=overlay_name, description=titl)

        lookat = make_lookat(lon_min, lon_max, lat_min, lat_max)
        kml_doc = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
{lookat}{fragment}</Document>
</kml>
"""

        with open(kmlfile, 'w') as f:
            f.write(kml_doc)
        print('Overlay written to '+pngfile+' / '+kmlfile)

        kml_folders.append(fragment)
        png_files.append(pngfile)

if KmzBundle:
    folders_joined = ''.join(kml_folders)
    agg_lookat = make_lookat(agg_lon_min, agg_lon_max, agg_lat_min, agg_lat_max)
    print('\nAggregate overlay bbox (lon/lat): '
          +str(round(agg_lon_min,4))+', '+str(round(agg_lat_min,4))
          +'  to  '+str(round(agg_lon_max,4))+', '+str(round(agg_lat_max,4)))
    print('If this does not fall on your survey area, EPSGCode is wrong.')

    doc_kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>{PlotName}</name>
{agg_lookat}{folders_joined}</Document>
</kml>
"""

    doc_kml_path = KmlDir+'doc.kml'
    with open(doc_kml_path, 'w') as f:
        f.write(doc_kml)

    kmz_path = KmlDir+KmzName
    with zipfile.ZipFile(kmz_path, 'w', zipfile.ZIP_DEFLATED) as kmz:
        kmz.write(doc_kml_path, arcname='doc.kml')
        for png in png_files:
            kmz.write(png, arcname=os.path.basename(png))

    print('\nKMZ bundle written to '+kmz_path)
    d = {'Title': KmzName,
         'Author': getpass.getuser(),
         'CreationDate': datetime.now().strftime('%d/%m/%Y %H:%M:%S')}
    print(d)
