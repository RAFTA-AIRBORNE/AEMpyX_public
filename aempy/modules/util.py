# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 17:08:06 2020

@author: vrath
"""

import os
import sys
import ast
import warnings
from sys import exit as error
import fnmatch
from datetime import datetime

import numpy
import random
import pyproj
from pyproj import CRS, Transformer

import shapely
import scipy.ndimage.filters
from scipy.fftpack import dct, idct
import scipy.spatial

import matplotlib
import matplotlib.pyplot


# import skimage.filters

import aesys
import inverse


def check_env(envar="CONDA_PREFIX", action="error"):
    """
    Check if environment variable exists

    Parameters
    ----------
    envar : strng, optional
        The default is ["CONDA_PREFIX"].

    Returns
    -------
    None.

    """
    act_env = os.environ[envar]    
    if len(act_env)>0:
        print("\n\n")
        print("Active conda Environment  is:  " + act_env)
        print("\n\n")
    else:
        if "err" in action.lower():
            error("Environment "+ act_env+"is not activated! Exit.")
            
def sample_list(in_list= [], method = ["sample", 10], out= True): 
    """
   

    Parameters
    ----------
    inlist : list of items 
        The default is [].
    method : list
        Determines samples . The default is ["sample", Nsample].
    out : boolean
        Output to stdout. The default is True.

    Returns
    -------
    outlist : list of items
    
    Created  Dec 2023
    @author: vrath

    """     
    if len(in_list)==0:
        error("sample_list: list empty! Exit.")

    if "pass" in method[0].lower() or method[0].lower()=="":
        out_list = in_list
        if out:
           print("sample_list: return original list")       
    
    if "rand" in method[0].lower():
        nsamples = method[1]
        out_list = random.sample(range(len(in_list)), nsamples)
        out_list = sorted(out_list)
        if out:
            print("sample_list: random samples = ", nsamples)
    
    if "step" in method[0].lower():

        start, step, stop = method[1:]
        out_list = in_list[start:stop:step]
        if out:
            print("sample_list: reduced list with start/stop/step = ", start,stop,step)
    
    return out_list 


def get_data_list(how=["search", ".npz", "./"],
                  sort=True, fullpath=False, out= True):
    """
    constructs a list of data files

    Parameters
    ----------
    how : list, required
             ["search". ".npz", directory].
             ["read", filename, directory]
    out :   boolean, optional
            Default is True.

    Returns
    -------
    data_files : list of filenames to be processed

    Created on Sun Jan 22 12:18:58 2023

    @author: vrath

    """
    if how[0].lower() !="search" and how[0].lower() !="read":
        error("No method given! Exit.")


    if "read" in how[0].lower():
        list_file = how[1]
        indir = how[2]
        print("Data files read from dir:  %s" % indir)
        print("Data list file:  "+list_file)
        dat_files = []
        with open(list_file, "r") as file:
            for line in file:
                dat_files.append(line[:-1])
        ns = numpy.size(dat_files)
        if ns ==0:
            error("No files in list <"+list_file+"> found!. Exit.")

    if "search" in how[0].lower():
        indir = how[2]
        print("Data files read from dir:  %s" % indir)
        SearchStrng = how[1]
        print("Search string: %s " % SearchStrng)
        dat_files = get_filelist(searchstr=[SearchStrng], searchpath=indir, fullpath=fullpath)
        ns = numpy.size(dat_files)
        if ns ==0:
            error("No files corresponding to searchstring <"+SearchStrng+"> found!. Exit.")


    #if fullpath:
       #dat_files = [os.path.join(searchpath,dat_files[ii]) for ii in range(len(dat_files))]
    #else:
       #dat_files = [os.path.join(searchpath,dat_files[ii]) for ii in range(len(dat_files))]


    if sort:
        dat_files = sorted(dat_files)

    if out:
        print(str(ns)+" data files found:")
        print(dat_files)

    return dat_files

def get_filebase(file=""):
    if len(file)==0:
        error("get_filebase: No file!. Exit")
        
        name, ext =os.path.splitext(os.path.basename(file))
        
        return  name, ext



def get_filelist(searchstr=["*"], searchpath="./", sortedlist =True, fullpath=False):
    """
    Generate filelist from path and unix wildcard list.

    author: VR 3/20

    last change 4/23
    """

    filelist = fnmatch.filter(os.listdir(searchpath), "*")
    for sstr in searchstr:
        filelist = fnmatch.filter(filelist, sstr)

    filelist = [os.path.basename(f) for f in filelist]

    if sortedlist:
        filelist = sorted(filelist)

    if fullpath:
       filelist = [os.path.join(searchpath,filelist[ii]) for ii in range(len(filelist))]

    return filelist


def merge_data_sets(infile_list=None, outfile_name="./tmp.npz", aem_system="aem05", dictout=False, out=False):
    """
    Merge data sets from file list
    Created on Sun Jan 22 12:18:58 2023
    @author: vrath
    """
    dateform="%m/%d/%Y, %H:%M:%S"
    _,NN, _, _, _, = aesys.get_system_params(System=aem_system)

    k = 0
    for file in infile_list:
        k = k+1
        if out: print("\nData read from: %s" % file)
        data_k, header, _ = aesys.read_aempy(File=file, System=aem_system, OutInfo=False)
        if k == 1:
            merged_data = data_k
        else:
            merged_data = numpy.vstack((merged_data, data_k))
    if outfile_name!=None:
        header = "merged data set:"+"".join("Date " + datetime.now().strftime(dateform))
        aesys.write_aempy(File=outfile_name, Data=merged_data, System=aem_system,
                Header=header, OutInfo=out)

    if dictout:
        merged_data = {
                "fnam": merged_data[:,0],
                "x": merged_data[:,1],
                "y": merged_data[:,2],
                "gps": merged_data[:,3],
                "alt": merged_data[:,4],
                "dem": merged_data[:,5],
                "data": merged_data[:,6:6+NN[2]]
                        }

    return merged_data

def merge_model_sets(infile_list=None, outfile_name="./tmp.npz", dictout= True, out=False):
    """
    Merge models from file list
    Created on Sun Jan 22 12:18:58 2023
    @author: vrath
    """
    dateform="%m/%d/%Y, %H:%M:%S"
    # _,NN, _, _, _, = aesys.get_system_params(System=aem_system)


    k = 0
    for file in infile_list:
        k = k+1
        print("\nData read from: %s" % file)
        results = numpy.load(file, allow_pickle=True)
        ctrl = results["ctrl"]

        mod_ref = results["mod_ref"]
        mod_act = results["mod_act"]

        site_mod = results["site_modl"]
        site_sns = results["site_sens"]
        site_rms = results["site_nrms"]

        site_y = results["site_y"]
        site_x = results["site_x"]
        site_gps = results["site_gps"]
        site_alt = results["site_alt"]
        site_dem = results["site_dem"]


        site_d = numpy.zeros_like(site_mod)
        site_z = numpy.zeros_like(site_mod)

        nlyr = inverse.get_nlyr(mod_ref)
        tmp = numpy.cumsum(mod_ref[6*nlyr:7*nlyr-1])
        tmp = numpy.insert(tmp, 0, 0.)
        tmp = numpy.append(tmp, tmp[-1])

        for isite in numpy.arange(numpy.shape(site_mod)[0]):
            site_d[isite] = 0.5*(tmp[0:len(tmp)-1]+tmp[1:len(tmp)])
            site_z[isite] = site_dem[isite] - site_d[isite]

        # print(numpy.shape(site_x))
        # print(numpy.shape(site_z))

        if k == 1:
            merged_mod = site_mod
            merged_sns = site_sns
            merged_rms = site_rms
            merged_x = site_x.reshape(-1,1)
            merged_y = site_y.reshape(-1,1)
            merged_z = site_z
            merged_d = site_d
            merged_gps = site_gps.reshape(-1,1)
            merged_alt = site_alt.reshape(-1,1)
            merged_dem = site_dem.reshape(-1,1)

        else:
            merged_mod = numpy.vstack((merged_mod, site_mod))
            merged_sns = numpy.vstack((merged_sns, site_sns))
            merged_rms = numpy.vstack((merged_rms, site_rms))
            merged_x = numpy.vstack((merged_x, site_x.reshape(-1,1)))
            merged_y = numpy.vstack((merged_y, site_y.reshape(-1,1)))
            merged_z = numpy.vstack((merged_z, site_z))
            merged_d = numpy.vstack((merged_d, site_d))
            merged_gps = numpy.vstack((merged_gps, site_gps.reshape(-1,1)))
            merged_alt = numpy.vstack((merged_alt, site_alt.reshape(-1,1)))
            merged_dem = numpy.vstack((merged_dem, site_dem.reshape(-1,1)))


    if outfile_name!=None:
        if not ".npz" in os.path.splitext(outfile_name)[1]:
            error("merge_model_sets: Only npz format implemented.! Exit.")

        header = "merged model set:"+"".join("Date " + datetime.now().strftime(dateform))
        numpy.savez_compressed(
            file=outfile_name,
            header=header,
            ctrl = ctrl,
            mod_act=mod_act,
            mod_ref=mod_ref,
            mod=merged_mod,
            sns=merged_sns,
            rms=merged_rms,
            x=merged_x,
            y=merged_y,
            z=merged_z,
            d=merged_d,
            gps=merged_gps,
            alt=merged_alt,
            dem=merged_dem)

    if dictout:
        merged_models = {
                "ctrl": ctrl,
                "mod_ref": mod_ref,
                "mod_act": mod_act,

                "mod": merged_mod,
                "sns": merged_sns,
                "rms": merged_rms,
                "x": merged_x,
                "y": merged_y,
                "z": merged_z,
                "d": merged_d,
                "gps": merged_gps,
                "alt": merged_alt,
                "dem": merged_dem,
                        }
    else:
        k=0
        for nd in numpy.arange(numpy.size(merged_x)):
            k=k+1
            tmp = numpy.array([merged_x[nd],
                               merged_y[nd],
                               merged_z[nd,:],
                               merged_d[nd,:],
                               merged_mod[nd,:],
                               merged_sns[nd,:],
                               merged_rms[nd],
                               merged_gps[nd],
                               merged_alt[nd],
                               merged_dem[nd]
                               ])

            if k ==1:
                merged_models=tmp
            else:
                merged_models=numpy.vstack((merged_models,tmp))

    return merged_models


def get_nearest_point(point=None,line=None, tol = 1.e-3):
    """
    Find nearest point in profile

    Parameters
    ----------
    point : float
        array/tuple of point coordinates.
    linedata : float
        array of point coordinates.

    Returns
    -------
    nearest : float
        position of nearest ppoint in profile

    """

    x0 = point[0]
    y0 = point[1]
    # dist =[]

    x1 = line[0]
    y1 = line[1]
    # print(numpy.abs(x1-x0))
    # print(numpy.abs(y1-y0))
    r = numpy.sqrt((x1-x0)**2+(y1-y0)**2)

    r_min = numpy.amin(r)
    # print(r)
    # print(numpy.isclose(r_min, r, 1.e-3))
    nearest =  numpy.where(numpy.isclose(r_min, r, tol))[0]
    # next((x for x in test_list if x.value == value), None)
    # numpy.where(r==r.min())
    # print(nearest)
    return nearest


def get_direction_angle(p1=None, p2=None):
    """
    Determine direction of vector from p1 to p2

    Parameters
    ----------
    p1 p2: float
        point coordinates, [Northing, Easting]

    Returns
    -------
    ang : float
        direction angle .

    """
    if (p1 == None) or (p1 == None):
        error("No angle can be determied as points  are not given! Exit.")
    v1 = numpy.asarray(p1)
    v2 = numpy.asarray(p2)
    dist = v2-v1
    length = numpy.linalg.norm(dist)
    d = dist/length
    rad = numpy.arctan2(d[0],d[1])
    ang = numpy.rad2deg(rad)%360.
    return ang, length

def list_functions(filename):
    """
    Generate list of functions in module.

    author: VR 3/21
    """

    print(filename)
    tree = parse_ast(filename)
    for func in find_functions(tree.body):
        print("  %s" % func.name)


def parse_ast(filename):
    with open(filename, "rt") as file:

        return ast.parse(file.read(), filename=filename)


def find_functions(body):
    return (f for f in body if isinstance(f, ast.FunctionDef))


def project_wgs_to_geoid(lat, lon, alt, geoid=3855 ):
    """
    transform ellipsoid heigth to geoid, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 09/21

    """

    geoidtrans =pyproj.crs.CompoundCRS(name="WGS 84 + EGM2008 height", components=[4979, geoid])
    wgs = pyproj.Transformer.from_crs(
            pyproj.CRS(4979), geoidtrans, always_xy=True)
    lat, lon, elev = wgs.transform(lat, lon, alt)

    return lat, lon, elev

def project_utm_to_geoid(utm_x, utm_y, utm_z, utm_zone=32629, geoid=3855):
    """
    transform ellipsoid heigth to geoid, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 09/21

    """

    geoidtrans =pyproj.crs.CompoundCRS(name="UTM + EGM2008 height", components=[utm_zone, geoid])
    utm = pyproj.Transformer.from_crs(
            pyproj.CRS(utm_zone), geoidtrans, always_xy=True)
    utm_x, utm_y, elev = utm.transform(utm_x, utm_y, utm_z)

    return utm_x, utm_y, elev



def project_gk_to_latlon(gk_x, gk_y, gk_zone=5684):
    """
    transform utm to latlon, using pyproj
    Look for other EPSG at https://epsg.io/
    VR 04/21
    """
    prj_wgs = pyproj.CRS("epsg:4326")
    prj_gk = pyproj.CRS("epsg:" + str(gk_zone))
    latitude, longitude = pyproj.transform(prj_gk, prj_wgs, gk_x, gk_y)
    return latitude, longitude


# def get_utm_zone(latitude=None, longitude=None):
#     """
#     Find EPSG from position in lat/lon, using pyproj

#     VR 04/21
#     """
#     prj_wgs = pyproj.CRS("epsg:4326")
#     utm_list = pyproj.query_utm_crs_info(
#         datum_name="WGS 84",
#         area_of_interest=pyproj.AreaOfInterest(
#         west_lon_degree=longitude,
#         south_lat_degree=latitude,
#         east_lon_degree=longitude,
#         north_lat_degree=latitude, ), )

#     return utm_list
 
def get_utm_zone(latitude=None, longitude=None):
    """
    Find EPSG from position, using pyproj

    VR 08/23
    """
    from pyproj.aoi import AreaOfInterest
    from pyproj.database import query_utm_crs_info
    utm_list = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
        west_lon_degree=longitude,
        south_lat_degree=latitude,
        east_lon_degree=longitude,
        north_lat_degree=latitude, ), )
    utm_crs =CRS.from_epsg(utm_list[0].code)
    EPSG = CRS.to_epsg(utm_crs)

    return EPSG, utm_crs

def project_latlon_to_utm(latitude, longitude, utm_zone=32629):
    """
    transform latlon to utm , using pyproj
    Look for other EPSG at https://epsg.io/

    VR 08/23
    """
    prj_wgs = CRS("epsg:4326")
    prj_utm = CRS("epsg:" + str(utm_zone))
    transformer = Transformer.from_crs(prj_wgs, prj_utm)
    utm_x, utm_y = transformer.transform(latitude, longitude)
    # transfor
    # prj_wgs = CRS("epsg:4326")
    # prj_utm = CRS("epsg:" + str(utm_zone))
    # utm_x, utm_y = pyproj.transform(prj_wgs, prj_utm, latitude, longitude)

    return utm_x, utm_y

def project_utm_to_latlon(utm_e, utm_n, utm_zone=32629):
    """
    transform latlon to utm , using pyproj
    Look for other EPSG at https://epsg.io/
    
    VR 08/23
    """    
    prj_wgs = CRS("epsg:4326")
    prj_utm = CRS("epsg:" + str(utm_zone))
    transformer = Transformer.from_crs(prj_utm, prj_wgs)
    latitude, longitude = transformer.transform(utm_e, utm_n)
    
    return latitude, longitude


def project_latlon_to_itm(longitude, latitude):
    """
    transform latlon to itm , using pyproj
    Look for other EPSG at https://epsg.io/

    VR 08/23
    """
    prj_wgs = CRS("epsg:4326")
    prj_itm = CRS("epsg:2157")
    transformer = Transformer.from_crs(prj_wgs, prj_itm)
    itm_e, itm_n = transformer.transform(latitude, longitude)

    return itm_e, itm_n


def project_itm_to_latlon(itm_e, itm_n):
    """
    transform itm to latlon, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 08/23
    """
    prj_wgs = CRS("epsg:4326")
    prj_itm = CRS("epsg:2157")
    transformer = Transformer.from_crs(prj_itm, prj_wgs)
    latitude, longitude = transformer.transform(itm_e, itm_n)

    return latitude, longitude


def project_itm_to_utm(itm_x, itm_y, utm_zone=32629):
    """
    transform itm to utm, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 08/23
    """
    prj_utm = CRS("epsg:" + str(utm_zone))
    prj_itm = CRS("epsg:2157") 
    transformer = Transformer.from_crs(prj_itm, prj_utm)
    utm_e, utm_n = transformer.transform(itm_x, itm_y)
    
    return utm_e, utm_n


def project_utm_to_itm(utm_e, utm_n, utm_zone=32629):
    """
    transform utm to itm, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 08/23
    """
    prj_utm = CRS("epsg:" + str(utm_zone))
    prj_itm = CRS("epsg:2157")
    
    transformer = Transformer.from_crs(prj_utm, prj_itm)
    itm_e, itm_n = transformer.transform(utm_e, utm_n)
    
    return itm_e, itm_n


def project_utm_to_utm(utm_e_in, utm_n_in, utm_zone_in=32629, utm_zone_out=32629):
    """
    transform utm to utm, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 08/23

    """
    if utm_zone_in==utm_zone_out:
        return  utm_e_in, utm_n_in
        
    prj_utm_in = CRS("epsg:" + str(utm_zone_in))
    prj_utm_out = CRS("epsg:" + str(utm_zone_out))
    
    transformer = Transformer.from_crs(prj_utm_in, prj_utm_out)
    utm_e, utm_n = transformer.transform(utm_e_in, utm_n_in)
    
    return utm_e, utm_n


def project_wgs_to_geoid(lat, lon, alt, geoid=3855 ):
    """
    transform ellipsoid heigth to geoid, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 09/21

    """

    geoidtrans =pyproj.crs.CompoundCRS(name="WGS 84 + EGM2008 height", components=[4979, geoid])
    wgs = pyproj.Transformer.from_crs(
            pyproj.CRS(4979), geoidtrans, always_xy=True)
    lat, lon, elev = wgs.transform(lat, lon, alt)

    return lat, lon, elev

def project_utm_to_geoid(utm_x, utm_y, utm_z, utm_zone=32629, geoid=3855):
    """
    transform ellipsoid heigth to geoid, using pyproj
    Look for other EPSG at https://epsg.io/

    VR 09/21

    """

    geoidtrans =pyproj.crs.CompoundCRS(name="UTM + EGM2008 height", components=[utm_zone, geoid])
    utm = pyproj.Transformer.from_crs(
            pyproj.CRS(utm_zone), geoidtrans, always_xy=True)
    utm_x, utm_y, elev = utm.transform(utm_x, utm_y, utm_z)

    return utm_x, utm_y, elev

def modify_polygon(Polygons=None,
                    Operator="intersection",
                    Params=[], Out=True):
    """
    VR 8/21
    """
    if not Polygons :
        error("No Polygons given!")

    print("Operator is " + Operator)

    PDims= numpy.size(Polygons)
    # print(PDims)
    if PDims == 1:

        Poly0 = Polygons[0]
        if ("rot" in Operator.lower()):
            Angle = Params[0]
            if numpy.size(Params)==1:
                Poly = shapely.affinity.rotate(Poly0, Angle, "center")
            else:
                Center = Params[1]
                Poly = shapely.affinity.rotate(Poly0, Angle, "center")

    else:
        Poly0 = Polygons[0]
        Poly1 = Polygons[1]

        if ("int" in Operator.lower()):
            Poly = Poly0.intersection(Poly1)

        if ("uni" in Operator.lower()):
            Poly = Poly0.union(Poly1)


    return Poly


def extract_data_poly(Data=None, PolyPoints=None, method=None, Out=True):
    """
     Chooses polygon area from aempy data set, given
     PolyPoints = [[X1 Y1],...[XN YN]]. First and last points will
     be connected for closure.

     VR 7/21
    """
    if method == None:
        error("extract_data_poly: No method given")

    if Data.size == 0:
        error("No Data given!")

    Ddims = numpy.shape(Data)
    if Out:
        print("data matrix input: " + str(Ddims))


    if "env" in method.lower():
        poly = numpy.column_stack((PolyPoints[:,0], PolyPoints[:,1]))
        poly = shapely.geometry.Multipoint((PolyPoints[:,0], PolyPoints[:,1])).envelope
    if "con" in method.lower():
        poly = numpy.column_stack((PolyPoints[:,0], PolyPoints[:,1]))
        poly = shapely.geometry.Multipoint((PolyPoints[:,0], PolyPoints[:,1])).convex_hull
    if "shp" in method.lower():
        poly = PolyPoints

    DPoly = []
    for row in numpy.arange(Ddims[0] - 1):
        if point_inside_polygon(Data[row, 1], Data[row, 2], poly):
            DPoly.append(Data[row, :])

    DPoly = numpy.asarray(DPoly, dtype=float)

    if Out:
        Ddims = numpy.shape(DPoly)
        print("data matrix output: " + str(Ddims))

    return DPoly


def point_inside_polygon(x, y, poly, method = "shapely"):
    """
    Determine if a point is inside a given polygon or not, where
    the polygon is given as a list of (x,y) pairs.
    Returns True  when point (x,y) ins inside polygon poly, False otherwise
    Based on shapelyp.geometry module(method="shapely") or pure python
    VR 8/21
    """
    if method==None:
        error("point_inside_polygon: No method given")

    if "sha" in method.lower():

        polygon = poly #shapely.geometry.Polygon(poly) # create polygon
        point = shapely.geometry.Point(x,y) # create point
        # contains= polygon.contains(point) # check if polygon contains point
        inside = point.within(polygon)      # check if point is in polygon

    else:
        n = len(poly)
        inside = False
        pv = numpy.column_stack((poly[:,0], poly[:,1]))
        p1x, p1y = pv[0]
        for i in range(n + 1):
            p2x, p2y = pv[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

    return inside


def extract_data_rect(Data=None, Corners=None, Out=True):
    """
     Chooses rectangular area from aempy data set, given
     the left lower and right uper corners in m as [minX maxX minY maxY]
     This only produces axis-parallel selections. If a rotation is desired,
     One could take the envelope (see modify_polygons) of  Corners, i.e.,
     [[X1 Y1], [X2 Y2]] treating this as a polygon.

    VR 7/21

    """

    if Data.size == 0:
        error("No Data given!")
    if not Corners:
        error("No Rectangle given!")
    print(Corners)
    Ddims = numpy.shape(Data)
    if Out:
        print("data matrix input: " + str(Ddims))
    Rect = []

    

    if Out:
        Emin = numpy.amin(Data[:, 1])
        Emax = numpy.amax(Data[:, 1])
        Nmin = numpy.amin(Data[:, 2])
        Nmax = numpy.amax(Data[:, 2])
        print("Easting:  "+str(Emin)+"-"+str(Emax))
        print("Northing: "+str(Nmin)+"-"+str(Nmax))
        
    X = [Corners[0], Corners[2]]
    Y = [Corners[1], Corners[3]]
    Xll = numpy.amin(X)
    Xur = numpy.amax(X)
    Yll = numpy.amin(Y)
    Yur = numpy.amax(Y)
    
    if Out:
        print("Rect lower left : "+str(Xll)+", "+str(Yll))
        print("Rect upper right: "+str(Xur)+", "+str(Yur))
    

    for row in numpy.arange(Ddims[0] - 1):
        if (Data[row, 1] > Xll and Data[row, 1] < Xur and
            Data[row, 2] > Yll and Data[row, 2] < Yur):
            Rect.append(Data[row, :])
            
    Rect = numpy.asarray(Rect, dtype=float)
    if Out:
        Ddims = numpy.shape(Rect)
        print("data matrix output: " + str(Ddims))

    return Rect


def project_to_line(x, y, line):
    """
    Projects a point onto a line, where line is represented by two arbitrary
    points. as an array

    VR 02/21
    """
    x1 = line[0, 0]
    x2 = line[1, 0]
    y1 = line[0, 1]
    y2 = line[1, 1]
    m = (y2 - y1) / (x2 - x1)
    b = y1 - (m * x1)

    xn = (m * y + x - m * b) / (m * m + 1.)
    yn = (m * m * y + m * x + b) / (m * m + 1.)

    return xn, yn


def gen_searchgrid(Points=None,
                   XLimits=None, dX=None, YLimits=None, dY=None, Out=False):
    """
    Generate equidistant grid for searching (in m).

    VR 02/21
    """
    small = 0.1

    datax = Points[:, 0]
    datay = Points[:, 1]
    nD = numpy.shape(Points)[0]

    X = numpy.arange(numpy.min(XLimits), numpy.max(XLimits) + small, dX)
    nXc = numpy.shape(X)[0]-1
    Y = numpy.arange(numpy.min(YLimits), numpy.max(YLimits)+ small, dY)
    nYc = numpy.shape(Y)[0]-1
    if Out:
        print("Mesh size: "+str(nXc)+"X"+str(nYc)
              +"\nCell sizes: "+str(dX)+"X"+str(dY)
              +"\nNuber of data = "+str(nD))


    p = numpy.zeros((nXc, nYc), dtype=object)
    # print(numpy.shape(p))


    xp= numpy.digitize(datax, X, right=False)
    yp= numpy.digitize(datay, Y, right=False)

    for ix in numpy.arange (nXc):
        incol = numpy.where(xp == ix)[0]
        for iy in numpy.arange (nYc):
            rlist = incol[numpy.where(yp[incol] == iy)[0]]
            p[ix,iy]=rlist

    # pout = numpy.array(p,dtype=object)

            if Out:
               print("mesh cell: "+str(ix)+" "+str(iy))

    return p



def splitall(path):
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path:  # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def get_files(SearchString=None, SearchDirectory="."):
    """
    FileList = get_files(Filterstring) produces a list
    of files from a searchstring (allows wildcards)

    VR 11/20
    """
    FileList = fnmatch.filter(os.listdir(SearchDirectory), SearchString)

    return FileList



def unique(list, out=False):
    """
    find unique elements in list/array

    VR 9/20
    """

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    if out:
        for x in unique_list:
            print(x)

    return unique_list



def strcount(keyword=None, fname=None):
    """
    count occurences of keyword in file
     Parameters
    ----------
    keywords : TYPE, optional
        DESCRIPTION. The default is None.
    fname : TYPE, optional
        DESCRIPTION. The default is None.

    VR 9/20
    """
    with open(fname, "r") as fin:
        return sum([1 for line in fin if keyword in line])
    # sum([1 for line in fin if keyword not in line])


def strdelete(keyword=None, fname_in=None, fname_out=None, out=True):
    """
    delete lines containing on of the keywords in list

    Parameters
    ----------
    keywords : TYPE, optional
        DESCRIPTION. The default is None.
    fname_in : TYPE, optional
        DESCRIPTION. The default is None.
    fname_out : TYPE, optional
        DESCRIPTION. The default is None.
    out : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    VR 9/20
    """
    nn = strcount(keyword, fname_in)

    if out:
        print(str(nn) + " occurances of <" + keyword + "> in " + fname_in)

    # if fname_out == None: fname_out= fname_in
    with open(fname_in, "r") as fin, open(fname_out, "w") as fou:
        for line in fin:
            if keyword not in line:
                fou.write(line)


def strreplace(key_in=None, key_out=None, fname_in=None, fname_out=None):
    """
    replaces key_in in keywords by key_out

    Parameters
    ----------
    key_in : TYPE, optional
        DESCRIPTION. The default is None.
    key_out : TYPE, optional
        DESCRIPTION. The default is None.
    fname_in : TYPE, optional
        DESCRIPTION. The default is None.
    fname_out : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    VR 9/20

    """

    with open(fname_in, "r") as fin, open(fname_out, "w") as fou:
        for line in fin:
            fou.write(line.replace(key_in, key_out))


def change_filename(old_filename="", how = ["append",""]):
    """
        Changes filename from template name.

        Parameters
        ----------
        old_filename : string
            input filename as string. The default is "".
        how : list
            One of the following:
                ["append", "append_string"]
                ["prepend"," prepend_string"]
                ["replace", "string_in","string_out"]"
            The default is ["append",""].

        Returns
        -------
       new_filename: string

    """

    if len(old_filename)==0:
        error("change_filename: No input file given! Exit.")

    new_filename = old_filename


    how[0] = how[0].lower()

    if "app" in how[0]:
        if len(how[1]) == 0:
            print("change_filename: No append string given!, nothing done")
        else:
            head, tail = os.path.split(old_filename)
            name, ext = os.path.splitext(tail)
            name =name+how[1]+ext
            new_filename = os.path.join(head,name)

    if "prep" in how[0]:
        if len(how[1])==0:
            print("change_filename: No prepend string given!, nothing done")
        else:
            head, tail = os.path.split(old_filename)
            name, ext = os.path.splitext(tail)
            name =how[1]+name+ext
            new_filename = os.path.join(head,name)

    if "repl" in how[0]:
        if  len(how[1])==0 or len(how[2])==0:
            print("change_filename: No replace strings given!, nothing done")
        else:
            new_filename = old_filename.replace(how[1], how[2])

    return new_filename

def gen_grid_latlon(
        LatLimits=None,
        nLat=None,
        LonLimits=None,
        nLon=None,
        out=True):
    """
     Generates equidistant 1-d grids in latLong.

     VR 11/20
    """
    small = 0.000001
# LonLimits = ( 6.275, 6.39)
# nLon = 31
    LonStep = (LonLimits[1] - LonLimits[0]) / nLon
    Lon = numpy.arange(LonLimits[0], LonLimits[1] + small, LonStep)

# LatLimits = (45.37,45.46)
# nLat = 31
    LatStep = (LatLimits[1] - LatLimits[0]) / nLat
    Lat = numpy.arange(LatLimits[0], LatLimits[1] + small, LatStep)

    return Lat, Lon


def gen_grid_utm(XLimits=None, nX=None, YLimits=None, nY=None, out=True):
    """
     Generates equidistant 1-d grids in m.

     VR 11/20
    """

    small = 0.000001
# LonLimits = ( 6.275, 6.39)
# nLon = 31
    XStep = (XLimits[1] - XLimits[0]) / nX
    X = numpy.arange(XLimits[0], XLimits[1] + small, XStep)

# LatLimits = (45.37,45.46)
# nLat = 31
    YStep = (YLimits[1] - YLimits[0]) / nY
    Y = numpy.arange(YLimits[0], YLimits[1] + small, YStep)

    return X, Y

def fractrans(m=None, x=None , a=0.5):
    """
    Caklculate fractional derivative of m.

    VR Apr 2021
    """
    import differint as df

    if m == None or x == None:
        error("No vector for diff given! Exit.")

    if numpy.size(m) != numpy.size(x):
        error("Vectors m and x have different length! Exit.")

    x0 = x[0]
    x1 = x[-1]
    npnts = numpy.size(x)
    mm = df.differint(a, m, x0, x1, npnts)

    return mm


def nearly_equal(a,b,sig_fig=6):
    return (a==b or int(a*10**sig_fig) == int(b*10**sig_fig))


def make_pdf_catalog(WorkDir="./", PdfList= None, FileName=None):
    """
    Make pdf catalog from site-plot.

    Parameters
    ----------
    Workdir : string
        Working directory.
    Filename : string
        Filename. Files to be appended must begin with this string.

    Returns
    -------
    None.

    """
    error("not in 3.9! Exit")
    import fitz

    catalog = fitz.open()

    for pdf in PdfList:
        with fitz.open(pdf) as mfile:
            catalog.insert_pdf(mfile)

    catalog.save(FileName, garbage=4, clean = True, deflate=True)
    catalog.close()

    print("\n"+str(numpy.size(PdfList))+" files collected to "+FileName)


def list_functions_in(this_module):

    from inspect import getmembers, isfunction
    import this_module

    funclist = getmembers(this_module, isfunction)
    for f in funclist:
        print(f[0])


def segment_distance(p, p1, p2, axis=None, return_t=False, segment=True):
    r"""
    Find the distance between an N-dimensional point and a line or line
    segment.
    The distance from a point to a line in N dimensions is the length
    of a normal dropped to the line. Using the fact that the dot
    product of orthogonal vectors is we can find the point
    :math:`\vec{p}_0` on the line that corresponds to this location.
    First, parametrize the points on the line through parameter
    :math:`t` as
    .. math::
       \vec{\ell} = \vec{p}_1 + t (\vec{p}_2 - \vec{p}_1)
    Then set up the equation with dot-products and solve for :math:`t`:
    .. math::
       (\vec{p} - \vec{p}_0) \cdot (\vec{p}_2 - \vec{p}_1) = 0
    .. math::
       (\vec{p} - \vec{p}_1 - t (\vec{p}_2 - \vec{p}_1)) \cdot (\vec{p}_2 - \vec{p}_1) = 0
    .. math::
       t (\vec{p}_2 - \vec{p}_1) \cdot (\vec{p}_2 - \vec{p}_1) = (\vec{p} - \vec{p}_1) \cdot (\vec{p}_2 - \vec{p}_1)
    .. math::
       t = \frac{(\vec{p} - \vec{p}_1) \cdot (\vec{p}_2 - \vec{p}_1)}{\left\lVert\vec{p}_2 - \vec{p}_1\right\rVert ^ 2}
    .. math::
       \vec{p}_0 = \vec{p}_1 + \frac{(\vec{p} - \vec{p}_1) \cdot (\vec{p}_2 - \vec{p}_1)}{\left\lVert\vec{p}_2 - \vec{p}_1\right\rVert ^ 2} (\vec{p}_2 - \vec{p}_1)
    The value of :math:`t` represents the location of :math:`\vec{p}_0`
    in relationship to :math:`\vec{p}_1` and :math:`\vec{p}_2`: values
    in the range :math:`[0, 1]` are on the line segment, negative
    values are on the side closer to :math:`\vec{p}_1`, and values
    greater than one are on the side of the line closer to
    :math:`\vec{p}_2`.
    The value of :math:`t` at the closest approach can be returned by
    setting ``return_t=True``. The value returned in this case applies
    to the entire line, even if ``segment == True`` and the closest
    point is one of the endpoints of the line segment.
    Parameters
    ----------
    p : array-like
        The target point. Must broadscast to `p1` and `p2`.
    p1 : array-like
        The start of the line segment. Must broadcast to the same shape
        as `p` and `p2`.
    p2 : array-like
        The end of the line segment. Must broadcast to the same shape
        as `p` and `p1`.
    axis : int or None
        The axis corresponding to the point vectors in the broadcasted
        arrays. If `None`, all point arrays are raveled.
    return_t: bool
        If `True`, return an additional value indicating the parameter
        :math:`t` at the distance of closest approach along the line.
        This will be the same regargless of `segment`.
    segment: bool
        If `True`, find the nearest point on the line segment bounded
        by `p1` and `p2` rather than the line passing between them.
    Returns
    -------
    dist : float or ~numpy.ndarray
        Distance from `p` to the line or line segment passing through
        `p1` and `p2`. The shape of the result is the broadcasted shape
        of the inputs, collapsed along `axis`.
        Scalar if `axis is None` or the inputs are all one-dimensional.
    t : float or ~numpy.ndarray
        An array of the same shape as `dist` containing the value of
        parameter :math:`t` for each line. The parameter is the
        location of the normal from `p` to the line passing through
        `p1` and `p2`, regardless if the distance is to the line
        segment or the line.
        Returned only if `return_t` is set.


        Created on Sun Jul 17 07:33:05 2022

    Based on:
    haggis: a library of general purpose utilities
    Copyright (C) 2019  Joseph R. Fox-Rabinovitz <jfoxrabinovitz at gmail dot com>

    VR July 2022

    """
    p, p1, p2 = numpy.broadcast_arrays(p, p1, p2)
    seg = p2 - p1
    norm2_seg = (seg * seg).sum(axis=axis, keepdims=True)
    t = ((p - p1) * seg).sum(axis=axis, keepdims=True) / norm2_seg
    p0 = p1 + t * seg

    if segment:
        dist = numpy.empty_like(p0)
        mask1 = t < 0
        mask2 = t > 1
        numpy.subtract(p0, p, where=~(mask1 | mask2), out=dist)
        numpy.subtract(p1, p, where=mask1, out=dist)
        numpy.subtract(p2, p, where=mask2, out=dist)
    else:
        dist = p0 - p
    dist = numpy.square(dist, out=dist).sum(axis=axis, keepdims=True)
    numpy.sqrt(dist, out=dist)

    if axis is None or seg.ndim == 1:
        dist = dist.item()
        t = t.item()
    else:
        dist = dist.squeeze(axis)
        t = t.squeeze(axis)

    if return_t:
        return dist, t

def find_nearest(site0=(0.,0.), sitevec=numpy.array([])):
    """
    Find smallest distance between points - brute force
    """
    if numpy.size(sitevec)==0:
        error("find_nearest: No vector uiven! exit.")

    x0 = site0[0]
    y0 = site0[1]
    dist = numpy.full_like(sitevec, numpy.nan)
    for site1 in numpy.arange(numpy.shape(sitevec)[0]):
        x1 = sitevec[site1,0]
        y1 = sitevec[site1,1]
        dist[site1] = numpy.sqrt((x1-x0)**2+(y1-y0)**2)

    minpos = numpy.argmin(dist)


    return minpos

def add_object_npz(filein=None,
                   xkey = [], xobject=numpy.array([]),
                   fileout = None):
    """
    Add object to .npz file.

    filein, fileout: str
        Valid .npz file.
    fileout: str
        Filename for output, default is filein
    xobject:
        Object to add to file (numpy.array)
    xkey: str
        Name of object

    vr Oct 31, 2022
    """

    if filein==None:
        error("File not given! Exit.")
    if fileout==None:
        fileout=filein

    if len(xobject)==0 or len(xkey)==0:
        error("Objects/keys not not given! Exit.")
    if len(xobject) != len(xkey):
        print("Object/key sizes do mot match! Set t0: "+xkey)

    tmp = numpy.load(filein, allow_pickle=True)
    tmp = dict(tmp)
    for iobj in numpy.arange(len(xkey)):
        tmp[xkey[iobj]] = xobject[iobj]
        print("Item "+xkey[iobj]+" added to "+fileout)

    numpy.savez_compressed(fileout,**tmp)



def del_object_npz(filein=None,
                   xkey = None, xobject=numpy.array([]),
                   fileout = None):
    """
    delete object from .npz file.

    filein:  str
        Valid .npz file.
    fileout: str
        Filename for output, default is filein

    xkey: list of str
        Names of objects to be deleted.

    vr Nov 12, 2022
    """

    if filein==None:
        error("File not given! Exit.")

    if xkey==None:
        error("No key  not given! Exit.")

    if fileout==None:
        fileout=filein

    tmp = numpy.load(filein, allow_pickle=True)
    tmp = dict(tmp)
    for iobj in numpy.arange(len(xkey)):
        if xkey[iobj] in tmp:
            tmp.pop(xkey[iobj])

    else:
        print("Key {"+xkey+"} is not in the dictionary")
        return

    numpy.savez_compressed(fileout,**tmp)


def stack_ragged(array_list, axis=0):
    """
    Stack ragged arrays.

    based on:
        https://tonysyu.github.io/ragged-arrays.html#.Y2D0UjPMIVM

    Parameters
    ----------
    array_list : list
        list of np.arrays.
    axis : TYPE, optional
        The default is 0.

    Returns
    -------
    stacked : np.array
        stacked array.
    idx : np.array
        points to end of stacked arrays.

    """
    lengths = [numpy.shape(a)[axis] for a in array_list]
    idx = numpy.cumsum(lengths[:-1])
    stacked = numpy.concatenate(array_list, axis=axis)
    return stacked, idx


def save_stacked_array(fname, array_list, axis=0, compressed =True):
    """
    Stack ragged arrays.

    based on:
        https://tonysyu.github.io/ragged-arrays.html#.Y2D0UjPMIVM

    Parameters
    ----------
    fnane: str
        filename
    array_list : list
        list of np.arrays.
    axis : TYPE, optional
        The default is 0.

    Returns
    -------
    stacked : np.array
        stacked array.
    idx : np.array
        points to end of stacked arrays.
    """
    stacked, idx = stack_ragged(array_list, axis=axis)
    if compressed:
        numpy.savez_compressed(fname, stacked_array=stacked, stacked_index=idx)
    else:
        numpy.savez(fname, stacked_array=stacked, stacked_index=idx)


def load_stacked_arrays(fname, axis=0):
    """
    Stack ragged arrays.

    based on:
        https://tonysyu.github.io/ragged-arrays.html#.Y2D0UjPMIVM

    Parameters
    ----------
    fname: str
        filename
    axis : TYPE, optional
        The default is 0.

    Returns
    -------
    stacked : np.array
        stacked array.
    idx : np.array
        points to end of stacked arrays.
    """
    npzfile = numpy.load(fname)
    idx = npzfile['stacked_index']
    stacked = npzfile['stacked_array']
    return numpy.split(stacked, idx, axis=axis)


def print_title(version="", fname="", form="%m/%d/%Y, %H:%M:%S", out=True):
    """
    Print version, calling filename, and modification date.
    """

    import os.path
    from datetime import datetime

    title = ""

    if len(version)==0:
        print("No version string given! Not printed to title.")
        tstr = ""
    else:
       ndat = "\n"+"".join("Date " + datetime.now().strftime(form))
       tstr =  "AEMpyX Version "+version+ndat+ "\n"

    if len(fname)==0:
        print("No calling filenane given! Not printed to title.")
        fstr = ""
    else:
        fnam = os.path.basename(fname)
        mdat = datetime.fromtimestamp((os.path.getmtime(fname))).strftime(form)
        fstr = fnam+", modified "+mdat+"\n"
        fstr = fstr + fname

    title = tstr+ fstr

    if out:
        print(title)

    return title


def anisodiff(img,niter=1,kappa=50,gamma=0.1,step=(1.,1.),option=1):
    """
    Anisotropic diffusion.

    Usage:
    imgout = anisodiff(im, niter, kappa, gamma, option)

    Arguments:
            img    - input image
            niter  - number of iterations
            kappa  - conduction coefficient 20-100 ?
            gamma  - max value of .25 for stability
            step   - tuple, the distance between adjacent pixels in (y,x)
            option - 1 Perona Malik diffusion equation No 1
                     2 Perona Malik diffusion equation No 2
            ploton - if True, the image will be plotted on every iteration

    Returns:
            imgout   - diffused image.

    kappa controls conduction as a function of gradient.  If kappa is low
    small intensity gradients are able to block conduction and hence diffusion
    across step edges.  A large value reduces the influence of intensity
    gradients on conduction.

    gamma controls speed of diffusion (you usually want it at a maximum of
    0.25)

    step is used to scale the gradients in case the spacing between adjacent
    pixels differs in the x and y axes

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

    April 2019 - Corrected for Python 3.7    -  AvW
    January 2022 simplified for application  -  VR

    """

    # ...you could always diffuse each color channel independently if you
    # really want
    if img.ndim == 3:
        warnings.warn("Only grayscale images allowed, converting to 2D matrix")
        img = img.mean(2)

    # initialize output array
    imgout = img.copy()

    # initialize some internal variables
    deltaS = numpy.zeros_like(imgout)
    deltaE = deltaS.copy()
    NS = deltaS.copy()
    EW = deltaS.copy()
    gS = numpy.ones_like(imgout)
    gE = gS.copy()

    for ii in range(niter):

        # calculate the diffs
        deltaS[:-1,: ] = numpy.diff(imgout,axis=0)
        deltaE[: ,:-1] = numpy.diff(imgout,axis=1)

        # conduction gradients (only need to compute one per dim!)
        if option == 1:
            gS = numpy.exp(-(deltaS/kappa)**2.)/step[0]
            gE = numpy.exp(-(deltaE/kappa)**2.)/step[1]
        elif option == 2:
            gS = 1./(1.+(deltaS/kappa)**2.)/step[0]
            gE = 1./(1.+(deltaE/kappa)**2.)/step[1]

        # update matrices
        E = gE*deltaE
        S = gS*deltaS

        # subtract a copy that has been shifted 'North/West' by one
        # pixel. don't as questions. just do it. trust me.
        NS[:] = S
        EW[:] = E
        NS[1:,:] -= S[:-1,:]
        EW[:,1:] -= E[:,:-1]

        # update the image
        imgout += gamma*(NS+EW)

    return imgout

def anisodiff3(stack,niter=1,kappa=50,gamma=0.1,step=(1.,1.,1.),option=1):
    """
    3D Anisotropic diffusion.

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

    Translated to Python and optimised by Alistair Muldal
    Department of Pharmacology
    University of Oxford
    <alistair.muldal@pharm.ox.ac.uk>

    June 2000  original version.
    March 2002 corrected diffusion eqn No 2.
    July 2012 translated to Python
    """

    # ...you could always diffuse each color channel independently if you
    # really want
    if stack.ndim == 4:
        warnings.warn("Only grayscale stacks allowed, converting to 3D matrix")
        stack = stack.mean(3)

    # initialize output array
    stackout = stack.copy()

    # initialize some internal variables
    deltaS = numpy.zeros_like(stackout)
    deltaE = deltaS.copy()
    deltaD = deltaS.copy()
    NS = deltaS.copy()
    EW = deltaS.copy()
    UD = deltaS.copy()
    gS = numpy.ones_like(stackout)
    gE = gS.copy()
    gD = gS.copy()

    for ii in range(niter):

        # calculate the diffs
        deltaD[:-1,: ,:  ] = numpy.diff(stackout,axis=0)
        deltaS[:  ,:-1,: ] = numpy.diff(stackout,axis=1)
        deltaE[:  ,: ,:-1] = numpy.diff(stackout,axis=2)

        # conduction gradients (only need to compute one per dim!)
        if option == 1:
            gD = numpy.exp(-(deltaD/kappa)**2.)/step[0]
            gS = numpy.exp(-(deltaS/kappa)**2.)/step[1]
            gE = numpy.exp(-(deltaE/kappa)**2.)/step[2]
        elif option == 2:
            gD = 1./(1.+(deltaD/kappa)**2.)/step[0]
            gS = 1./(1.+(deltaS/kappa)**2.)/step[1]
            gE = 1./(1.+(deltaE/kappa)**2.)/step[2]

        # update matrices
        D = gD*deltaD
        E = gE*deltaE
        S = gS*deltaS

        # subtract a copy that has been shifted 'Up/North/West' by one
        # pixel. don't as questions. just do it. trust me.
        UD[:] = D
        NS[:] = S
        EW[:] = E
        UD[1:,: ,: ] -= D[:-1,:  ,:  ]
        NS[: ,1:,: ] -= S[:  ,:-1,:  ]
        EW[: ,: ,1:] -= E[:  ,:  ,:-1]

        # update the image
        stackout += gamma*(UD+NS+EW)

    return stackout

def ms2ohmm(ms):
    s = 1.e-3*ms
    return 1./s
