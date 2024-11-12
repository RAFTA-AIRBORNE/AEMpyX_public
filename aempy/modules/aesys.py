"""
Module controlling data IO.

Created on Sun Nov  1 17:08:06 2020

@author: vrath
@author: duygu - edited on 12 June 2021
@author: duygu - edited on 6 April 2023
"""
import os
from sys import exit as error
import numpy
from datetime import datetime
import util
import netCDF4 as nc

# from numba import njit

def get_system_params(System="aem05", OutInfo = True):
    """
    Get system parameters, format and columnheaders, for AEM systems.

    Parameters
    ----------
    System : string
        The string should be a unique indicator for the system.
        Should be consistent with naming in the Fortran source code
        in core library.

    Returns
    -------
    fwdcall :  string
        Forward modeling call
    NN  : list
        This is a list of lengths of different parameter sets in the
        internal aempy format.
        [line, meta information, data,other]
        meta information includes flight line, position (UTM), height asl,
        radar altitude, and DTM.
        other could be a power line monitor where available, and
        other available quality related information.
    fmt :  string
        Output format for ASCII data in aempy
    col :  string
        Fcolumn names for data in aempy
    miscpars: List of objects
        Other pars

    Author: vrath, 2021/02/15

    """
    if System.lower() == "aem05":

        n_meta = 6
        n_data = 8
        n_optn = 3
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  # 6+8+1+2

        fwdcall = "core1d.aemfwd1d_aem05(nlyr, m, alt)"

        fmt = \
            "%10.2f "+"%14.2f "*2+"   "+" %10.2f "*3+"   "+"%10.2f "*11
        cols = [
            "line,",
            "x,", "y,", "gps,", "alt,", "dem,",
            "p1,", "p2,", "p3,", "p4,", "q1,", "q2,", "q3,", "q4,",
            "pli,", "qflag,","pflag"]
        col = ""
        for field in cols:
            col = col + "{:>6}".format(field[:6])

        freqs = numpy.asarray([912., 3005., 11962.,24510.])
        xunits = "(Hz)"
        dunits= "(ppm)"
        compdict=\
        dict([
            ("P1",  [6, freqs[0], " - 912 Hz"]), ("Q1", [10, freqs[1], " - 912 Hz"]),
            ("P2",  [7, freqs[2],  " - 3005 Hz"]), ("Q2", [11, freqs[3]," - 3005 Hz"]),
            ("P3",  [8, freqs[0],  " - 11962 Hz"]), ("Q3", [12, freqs[1],  " - 11962 Hz"]),
            ("P4",  [9, freqs[2],  " - 24510 Hz"]), ("Q4", [13, freqs[3],  " - 24510 Hz"]),
            ("PLM", [14, 0., ""]),("ALT", [4, 0., ""])
            ])
        # CompLabl = list(compdict.keys())
        miscpars = [freqs, xunits, dunits, compdict]

    elif System.lower() == "genesis":
        n_meta = 6
        n_data = 22
        n_optn = 0
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  #6+22+2

        fwdcall = "core1d.aemfwd1d_genesis(nlyr, m, alt)"

        fmt = \
            "%10.1f " + "%14.2f "*2 + "%10.2f "*3 + "%10.2f "*24
        cols = [
            "line,",
            "x,", "y,", "gps,", "alt,", "dem,",
            "x1,", "x2,", "x3,", "x4,", "x5,",
            "x6,", "x7,", "x8,", "x9,", "x10,", "x11,",
            "z1,", "z2,", "z3,", "z4,", "z5,",
            "z6,", "z7,", "z8,", "z9,", "z10,", "z11,", "qflag,","pflag"]
        col = ""
        for field in cols:
            col = col + "{:>6}".format(field[:6])


        w0 = [0.01,  0.017, 0.035, 0.069, 0.122, 0.191, 0.295, 0.434, 0.660, 1.007, 1.51]
        w1 = [0.017, 0.035, 0.069, 0.122, 0.191, 0.295, 0.434, 0.660, 1.007, 1.51, 2.205]
        # wct = 0.5 *(numpy.log10(w0) + numpy.log10(w1))
        wct = 0.5 *(numpy.array(w0) + numpy.array(w1))
        # wct =numpy.power(10.,wct)
        xunits = "(ms)"
        dunits= "(ppm)"
        compdict =\
            dict([
                ("H1", [ 6, wct[0], " - 0.0135 ms"]), ("H2", [ 7, wct[1], " - 0.026 ms"]),
                ("H3", [ 8, wct[2], " - 0.052 ms"]),  ("H4", [ 9, wct[3], " - 0.0955 ms"]),
                ("H5", [10, wct[4], " - 0.1565 ms"]), ("H6", [11, wct[5], " - 0.243 ms"]),
                ("H7", [12, wct[6], " - 0.3645 ms"]), ("H8", [13, wct[7], " - 0.547 ms"]),
                ("H9", [14, wct[8], " - 0.8335 ms"]), ("H10", [15, wct[9], " - 1.2585 ms"]),
                ("H11", [16, wct[10], " - 1.8575 ms"]),
                ("Z1", [17, wct[0], " - 0.0135 ms"]), ("Z2", [18, wct[1]," - 0.026 ms"]),
                ("Z3", [19, wct[2]," - 0.052 ms"]),   ("Z4", [20, wct[3], " - 0.0955 ms"]),
                ("Z5", [21, wct[4], " - 0.1565 ms"]), ("Z6", [22, wct[5], " - 0.243 ms"]),
                ("Z7", [23, wct[6], " - 0.3645 ms"]), ("Z8", [24, wct[7], " - 0.547 ms"]),
                ("Z9", [25, wct[8], " - 0.8335 ms"]), ("Z10", [26, wct[9], " - 1.2585 ms"]),
                ("Z11", [27, wct[10], " - 1.8575 ms"]),
                ("ALT",  [4, 0., ""])
               ])
        # CompLabl = list(compdict.keys())

        miscpars = [wct, xunits, dunits, compdict]

    elif System.lower() == "sglo":
        """
        special version for overlap a1/nm
        """


        n_meta = 6
        n_data = 8
        n_optn = 3
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  # 6+8+1+2

        fwdcall = "core1d.aemfwd1d_aem05(nlyr, m, alt)"

        fmt = \
            "%10.2f "+"%14.2f "*2+"   "+" %10.2f "*3+"   "+"%10.2f "*11
        cols = [
            "line,",
            "x,", "y,", "gps,", "alt,", "dem,",
            "p1,", "p2,", "p3,", "p4,", "q1,", "q2,", "q3,", "q4,",
            "pli,", "qflag,","pflag"]
        col = ""
        for field in cols:
            col = col + "{:>6}".format(field[:6])

        freqs = numpy.asarray([912., 3005., 11962.,24510.])
        xunits = "(Hz)"
        dunits= "(ppm)"
        compdict=\
        dict([
            ("P1",  [6, freqs[0], " - 912 Hz"]), ("Q1", [10, freqs[1], " - 912 Hz"]),
            ("P2",  [7, freqs[2],  " - 3005 Hz"]), ("Q2", [11, freqs[3]," - 3005 Hz"]),
            ("P3",  [8, freqs[0],  " - 11962 Hz"]), ("Q3", [12, freqs[1],  " - 11962 Hz"]),
            ("P4",  [9, freqs[2],  " - 24510 Hz"]), ("Q4", [13, freqs[3],  " - 24510 Hz"]),
            ("PLM", [14, 0., ""]),("ALT", 0., [4, ""])
            ])
        # CompLabl = list(compdict.keys())
        miscpars = [freqs, xunits, dunits, compdict]

    elif System.lower() == "cggo":
        """
        special version for overlap a1/nm
        """

        n_meta = 6
        n_data = 22
        n_optn = 2
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  #6+22+2

        fwdcall = "core1d.aemfwd1d_genesis(nlyr, m, alt)"

        fmt = \
            "%10.1f " + "%16.8g "*2 + "%10.2f "*3 + "%10.2f "*24
        cols = [
            "line,",
            "x,", "y,", "gps,", "alt,", "dem,",
            "x1,", "x2,", "x3,", "x4,", "x5,",
            "x6,", "x7,", "x8,", "x9,", "x10,", "x11,",
            "z1,", "z2,", "z3,", "z4,", "z5,",
            "z6,", "z7,", "z8,", "z9,", "z10,", "z11,", "qflag,","pflag"]
        col = ""
        for field in cols:
            col = col + "{:>6}".format(field[:6])

        n_meta = 6
        n_data = 22
        n_optn = 2
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  #6+22+2

        w0 = [0.01,  0.017, 0.035, 0.069, 0.122, 0.191, 0.295, 0.434, 0.660, 1.007, 1.51]
        w1 = [0.017, 0.035, 0.069, 0.122, 0.191, 0.295, 0.434, 0.660, 1.007, 1.51, 2.205]
        # wct = 0.5 *(numpy.log10(w0) + numpy.log10(w1))
        wct = 0.5 *(numpy.array(w0) + numpy.array(w1))
        # wct =numpy.power(10.,wct)
        xunits = "(ms)"
        dunits= "(fT)"
        compdict =\
            dict([
                ("H1", [ 6, wct[0], " - 0.0135"]), ("H2", [ 7, wct[1], " - 0.026"]),
                ("H3", [ 8, wct[2], " - 0.052"]),  ("H4", [ 9, wct[3], " - 0.0955"]),
                ("H5", [10, wct[4], " - 0.1565"]), ("H6", [11, wct[5], " - 0.243"]),
                ("H7", [12, wct[6], " - 0.3645"]), ("H8", [13, wct[7], " - 0.547"]),
                ("H9", [14, wct[8], " - 0.8335"]), ("H10", [15, wct[9], " - 1.2585"]),
                ("H11", [16, wct[10], " - 1.8575"]),
                ("Z1", [17, wct[0], " - 0.0135"]), ("Z2", [18, wct[1]," - 0.026"]),
                ("Z3", [19, wct[2]," - 0.052"]),   ("Z4", [20, wct[3], " - 0.0955"]),
                ("Z5", [21, wct[4], " - 0.1565"]), ("Z6", [22, wct[5], " - 0.243"]),
                ("Z7", [23, wct[6], " - 0.3645"]), ("Z8", [24, wct[7], " - 0.547"]),
                ("Z9", [25, wct[8], " - 0.8335"]), ("Z10", [26, wct[9], " - 1.2585"]),
                ("Z11", [27, wct[10], " - 1.8575"]),
                ("ALT",  [4, 0., ""])
               ])
        # CompLabl = list(compdict.keys())


        miscpars = [wct, xunits, dunits, compdict]

    elif System.lower() == "sglt":
        """
        special version for testlines DKV 2062
        """

        n_meta = 6
        n_data = 8
        n_optn = 3
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  # 6+8+1+1

        fwdcall = "core1d.aemfwd1d_aem05(nlyr, mm, alt)"

        fmt = \
            "%10.2f "+"%14.2f "*2+"   "+" %10.2f "*3+"   "+"%10.2f "*11
        cols = [
            "line,",
            "x,", "y,", "gps,", "alt,", "dem,",
            "p1,", "p2,", "p3,", "p4,", "q1,", "q2,", "q3,", "q4,",
            "pli,", "qflag,","pflag"]
        col = ""
        for field in cols:
            col = col + "{:>6}".format(field[:6])
        freqs = numpy.asarray([912., 3005., 11962.,24510.])
        xunits = "(Hz)"
        dunits= "(ppm)"
        compdict=\
        dict([
            ("P1",  [6, freqs[0], " - 912 Hz"]), ("Q1", [10, freqs[1], " - 912 Hz"]),
            ("P2",  [7, freqs[2],  " - 3005 Hz"]), ("Q2", [11, freqs[3]," - 3005 Hz"]),
            ("P3",  [8, freqs[0],  " - 11962 Hz"]), ("Q3", [12, freqs[1],  " - 11962 Hz"]),
            ("P4",  [9, freqs[2],  " - 24510 Hz"]), ("Q4", [13, freqs[3],  " - 24510 Hz"]),
            ("PLM", [14, 0., ""]),("ALT", 0., [4, ""])
            ])
        # CompLabl = list(compdict.keys())

        miscpars = [freqs, xunits, dunits, compdict]

    elif System.lower() == "cggt":
        """
        special version for testlines
        """

        n_meta = 6
        n_data = 22
        n_optn = 2
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  #6+22+1

        fwdcall = "core1d.aemfwd1d_genesis(nlyr, mm, alt)"

        fmt = \
            "%10.1f " + "%16.8g "*2 + "%10.2f "*3 + "%10.2f "*24
        cols = [
            "line,",
            "x,", "y,", "gps,", "alt,", "dem,",
            "x1,", "x2,", "x3,", "x4,", "x5,",
            "x6,", "x7,", "x8,", "x9,", "x10,", "x11,",
            "z1,", "z2,", "z3,", "z4,", "z5,",
            "z6,", "z7,", "z8,", "z9,", "z10,", "z11,", "qflag, ","pflag"]
        col = ""
        for field in cols:
            col = col + "{:>6}".format(field[:6])

        n_meta = 6
        n_data = 22
        n_optn = 2
        nn =[n_meta+n_data+n_optn, n_meta, n_data, n_optn]  #6+22+2

        w0 = [0.01,  0.017, 0.035, 0.069, 0.122, 0.191, 0.295, 0.434, 0.660, 1.007, 1.51]
        w1 = [0.017, 0.035, 0.069, 0.122, 0.191, 0.295, 0.434, 0.660, 1.007, 1.51, 2.205]
        # wct = 0.5 *(numpy.log10(w0) + numpy.log10(w1))
        wct = 0.5 *(numpy.array(w0) + numpy.array(w1))
        # wct =numpy.power(10.,wct)
        xunits = "(ms)"
        dunits= "(fT)"

        compdict =\
            dict([
                ("H1", [ 6, wct[0], " - 0.0135"]), ("H2", [ 7, wct[1], " - 0.026"]),
                ("H3", [ 8, wct[2], " - 0.052"]),  ("H4", [ 9, wct[3], " - 0.0955"]),
                ("H5", [10, wct[4], " - 0.1565"]), ("H6", [11, wct[5], " - 0.243"]),
                ("H7", [12, wct[6], " - 0.3645"]), ("H8", [13, wct[7], " - 0.547"]),
                ("H9", [14, wct[8], " - 0.8335"]), ("H10", [15, wct[9], " - 1.2585"]),
                ("H11", [16, wct[10], " - 1.8575"]),
                ("Z1", [17, wct[0], " - 0.0135"]), ("Z2", [18, wct[1]," - 0.026"]),
                ("Z3", [19, wct[2]," - 0.052"]),   ("Z4", [20, wct[3], " - 0.0955"]),
                ("Z5", [21, wct[4], " - 0.1565"]), ("Z6", [22, wct[5], " - 0.243"]),
                ("Z7", [23, wct[6], " - 0.3645"]), ("Z8", [24, wct[7], " - 0.547"]),
                ("Z9", [25, wct[8], " - 0.8335"]), ("Z10", [26, wct[9], " - 1.2585"]),
                ("Z11", [27, wct[10], " - 1.8575"]),
                ("ALT",  [4, 0., ""])
               ])
        # CompLabl = list(compdict.keys())

        miscpars = [wct, xunits, dunits, compdict]

    elif System.lower() == "geotem" or System.lower() == "tempest" :
        fwdcall = ""
        error("System " + System + " not yet implemented! Exit.")

    else:
        error("System " + System + " not implemented! Exit.")

    if OutInfo:
        print("\nAEM system is "+System)
        print("Forward model call: "+fwdcall)
        print("Data:"+str(nn))

    return fwdcall, nn, fmt, col, miscpars

def read_survey_data(DatFile=None, Survey="A5", OutInfo=False, Invalid=numpy.nan,
                     EPSG_in=2157,
                     EPSG_out=32629):
    """
    Read raw data.

    """
    print("Survey is "+Survey.lower())

    # Invalid could also be fac*numpy.finfo(float).max

    invalid_strng = str(Invalid)

    if any(s in Survey.lower() for s in ["a1", "a2", "a3", "a4", "wf", "tb", ]):

        """
        ===================================================================================================================
        Surveys A1 - A4, TB, WF (NEW, 2022):

        ====================================================================================================================
        Project Code: GSI-D_17.IRL Eire
        Project:      Airborne Magnetic, Radiometric and Frequency Domain Electromagnetic Survey in the Republic of Ireland
        Client:       Geological Survey of Ireland (GSI)
        Delivery:     DLV2088 FEM + DLV2424 PLM_nT data
        Date:         11th March 2022
        Comment:      GSI merge of DLV2088 FEM data with DLV2424 PLM_nT data
        ====================================================================================================================


        File Name:  GSI-D_17.IRL_DLV2088_FEM_PLM_A1.xyz
        File Name:  GSI-D_17.IRL_DLV2088_FEM_PLM_A2.xyz
        File Name:  GSI-D_17.IRL_DLV2088_FEM_PLM_A3.xyz
        File Name:  GSI-D_17.IRL_DLV2088_FEM_PLM_A4.xyz
        File Name:  GSI-D_17.IRL_DLV2088_FEM_PLM_Tellus_Border.xyz
        File Name:  GSI-D_17.IRL_DLV2088_FEM_PLM_Waterford.xyz

        Translation table from .xyz files to internal data format:

        Name              Storage        Pos        units              Description
        Rerelease2022:
        ITM_X            1               0         m       	X coordinate, IRENET95 ITM
        ITM_Y            2               1	       m       	Y coordinate, IRENET95 ITM
        DATE           	                 2	                    Date YYYYMMDD
        LINE             0               3                     Line Number  BBLLLL.SR
        LONG         		             4      degree    	    Longitude, WGS-84
        LAT          		             5      degree    	    Latitude, WGS-84
        MSLHGT           3               6         m		    GPS Elevation above Mean Sea Level
        RADAR		      4	              7         m		    Clearance above Terrain from Radar
        DEM              5	              8 		m           Digital Elevation Model with respect to Mean Sea Level from Radar Clearance
        P09lev        	  6	              9         ppm		    Levelled and filtered in-phase 912 Hz
        Q09lev        	  7	             10         ppm		    Levelled and filtered quadrature 912 Hz
        P3lev         	  8	             11         ppm		    Levelled and filtered in-phase 3005 Hz
        Q3lev         	  9	             12         ppm		    Levelled and filtered quadrature 3005 Hz
        P12lev        	 10	             13         ppm		    Levelled and filtered in-phase 11962 Hz
        Q12lev        	 11	             14         ppm	    	Levelled and filtered quadrature 11962 Hz
        P25lev        	 12	             15         ppm	    	Levelled and filtered in-phase 24510 Hz
        Q25lev        	 13	             16         ppm		    Levelled and filtered quadrature 24510 Hz
        PLM_nT           14 	         17         nT		    Power line monitor                                                                (B=Block, L=line, S=segment, R=reflight)

        DLV2088:
        LINE             0               0                              Line Number  BBLLLL.SR
        DATE                             1           -                  Date YYYYMMDD
        DAY                              2           -                  Day of year
        TIME                             3          sec                 Fiducial Seconds
        LAT                              4         degree               Latitude, WGS-84
        LONG                             5         degree               Longitude, WGS-84
        ITM-X             1              6           m                  X coordinate, IRENET95 ITM
        ITM-Y             2              7           m                  Y coordinate, IRENET95 ITM
        UTM-X                            8           m                  X coordinate, WGS-84 UTM 29N
        UTM-Y                            9           m                  Y coordinate, WGS-84 UTM 29N
        UTM-Z                            10           m                 GPS Elevation above WGS-84 Ellipsoid
        MSLHGT            3              11           m                 GPS Elevation above Mean Sea Level
        RADAR             4              12           m                 Radar Altimeter
        P09lev            5              13          ppm                Levelled and filtered in-phase 912 Hz
        Q09lev            9              14          ppm                Levelled and filtered quadrature 912 Hz
        P3lev             6              15          ppm                Levelled and filtered in-phase 3005 Hz
        Q3lev            10              16          ppm                Levelled and filtered quadrature 3005 Hz
        P12lev            7              17          ppm                Levelled and filtered in-phase 11962 Hz
        Q12lev           11              18          ppm                Levelled and filtered quadrature 11962 Hz
        P25lev            8              19          ppm               Levelled and filtered in-phase 24510 Hz
        Q25lev           12              20          ppm                Levelled and filtered quadrature 24510 Hz
        PLM_nT           13              21           nT                Power line monitor


        """
        # ncol =  [0, 6, 7,    11, 12,     13, 15, 17, 19,  14, 16, 18, 20, 21]       # DLV2088
        ncol =  [3,  0, 1,    6, 7, 8,     9, 11, 13, 15,  10, 12, 14, 16, 17]        #2022

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[0].lower().startswith("#")
                    or line[0].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie" in line[:24].lower()):
                    continue

                t = line.split()

                t = [w.replace("*", invalid_strng) for w in t]
                # print(t[:22])
                tmp = [t[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        EPSG_in = 2157
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                          utm_zone_in=EPSG_in,
                                                          utm_zone_out=EPSG_out)

        # Data = numpy.insert(Data, 5, Data[:, 3]-Data[:, 4], axis=1)
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))


    if any(s in Survey.lower() for s in [ "cav","cv", ]):
        """
        =======================================================================
        Surveys cavan2006 :
        =======================================================================
        0       1    X              Easting coordinate (Irish Grid)
        1       2    Y              Northing coordinate (Irish Grid)
        2   LAT    WGS84 Latitude
        3   LON    WGS84 Longitude
        4   FID    Fiducial
        5       0    FLIGHT	        Flight number
        6   DATE    Survey date
        7   DAY	    Julian day
        8   TIME    Recording time in HHMMSS.SSS format
        9   DIR	    Flight direction (from North)
        10      4    RALT       Flight altitude
        11      3    GPS_H      GPS altitude
        12      5    DTM        Terrain model (elevation from WGS-84 ellipsoid surface)
        13      14   PLM        Powerline monitor
        14      6    RE09       Real component, 912 Hz
        15      10   M09        Imaginary component, 912 Hz
        16      7    RE3        Real component, 3005 Hz
        17      11   IM3        Imaginary component, 3005 Hz
        18      8    RE12       Real component, 11962 Hz
        19      12   IM12       Imaginary component, 11962 Hz
        20      9    RE25       Real component, 24510 Hz
        21      13   IM25       Imaginary component, 24510 Hz
   """
        ncol =  [5, 0, 1,    11, 10, 12,    14, 16, 18, 20, 15, 17, 19, 21, 13]

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[0].lower().startswith("#")
                    or line[0].lower().startswith("/")):
                    continue

                if ("line" in line.lower() or "tie" in line.lower()):
                    t = line.split()
                    fl = numpy.float(t[1])
                    print(t)
                    continue

                t = line.split()

                t = [w.replace("*", invalid_strng) for w in t]
                # print(t[:22])
                tmp = [t[ii] for ii in ncol]
                tmp[0] = fl
                tmp[1] = numpy.float(tmp[1])
                tmp[2] = numpy.float(tmp[2])
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        EPSG_in = 29903
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                          utm_zone_in=EPSG_in,
                                                          utm_zone_out=EPSG_out)

        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))


    if any(s in Survey.lower() for s in [ "a5","a6", "a7", ]):

        """
        ===================================================================================================================
        Surveys A5 - A7 :
        ===================================================================================================================

        Translation table from .xyz files to internal data format:

        AEMPY       N           Name            Units        Description
        0           0           LINE                         Line Number - LLLL.SR
                                                            (L=line, S=segment, R=reflight)
                    1           FLT                          Flight Number
                    2           DATE  	      -        -     Date YYYYMMDD
                    3           DAY		      -        -     Day of year
                    4           TIME  	        sec          Fiducial Seconds
                    5           LAT		        degree       Latitude, WGS-84
                    6           LONG  	        degree       Longitude, WGS-84
        1           7           ITM-X 	        m            X coordinate, IRENET95 ITM (EPSG 2157)
        2           8           ITM-Y 	        m            Y coordinate, IRENET95 ITM (EPSG 2157)
                    9           UTM-X 	        m            X coordinate, WGS-84 UTM 29N
                   10           UTM-Y 	        m            Y coordinate, WGS-84 UTM 29N
                   11           UTM-Z 	        m            GPS Elevation above WGS-84 Ellipsoid
        3          12           MSLHGT	        m            GPS Elevation above Mean Sea Level
        4          13           CLEARANCE	    m            Clearance above Terrain from Laser
        5          14           DEM		        m            DEM from Laser with respect
                                                             to Mean Sea Level
                   15           TEMP  	        degree C     Temperature
                   16-31                                     Raw/filtered data
        6          32           P09lev          ppm          Levelled and filtered in-phase 912 Hz
       10          33           Q09lev          ppm          Levelled and filtered quadrature 912 Hz
        7          34           P3lev           ppm          Levelled and filtered in-phase 3005 Hz
       11          35           Q3lev           ppm          Levelled and filtered quadrature 3005 Hz
        8          36           P12lev          ppm          Levelled and filtered in-phase 11962 Hz
       12          37           Q12lev          ppm          Levelled and filtered quadrature 11962 Hz
        9          38           P25lev          ppm          Levelled and filtered in-phase 24510 Hz
       13          39           Q25lev          ppm          Levelled and filtered quadrature 24510 Hz
                   40           Radio_Flag      -            Radio call flag
       14          41           PLM_nT          nT           Power line monitor

       """
        ncol = [0, 7, 8,  12, 13, 14,   32, 34, 36, 38, 33, 35, 37, 39,   41]

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[0].lower().startswith("#")
                    or line[0].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie" in line[:24].lower()):
                    continue

                t = line.split()

                t = [w.replace("*", invalid_strng) for w in t]
                # print(t[:22])
                tmp = [t[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        EPSG_in = 2157
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                          utm_zone_in=EPSG_in,
                                                          utm_zone_out=EPSG_out)

        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))



    if any(s in Survey.lower() for s in [ "a8","a9", ]):

        """
       ===================================================================================================================
       Surveys A8 - A9 (NEW):
       ===================================================================================================================

         Translation table from .xyz files to internal data format:

         AEMPY  N        Name          Units     	Description
         1      0       ITM_X          m       	X coordinate, IRENET95 ITM
         2      1       ITM_Y          m       	Y coordinate, IRENET95 ITM
                2       DATE           -        Date YYYYMMDD
         0      3       LINE		   -		Line number
                4       LONG           degree   Longitude, WGS-84
                5       LAT            degree   Latitude, WGS-84
         3      6       MSLHGT         m		GPS Elevation above Mean Sea Level
         4      7       CLEARANCE      m		Clearance above Terrain from Laser
         5      8       DEM            m		Digital Elevation Model with respect
                                                to Mean Sea Level from Laser Clearance
                9       TEMP           degree C	Temperature
         6      10      P09lev         ppm		Levelled and filtered in-phase 912 Hz
        10      11      Q09lev         ppm		Levelled and filtered quadrature 912 Hz
         7      12      P3lev          ppm		Levelled and filtered in-phase 3005 Hz
        11      13      Q3lev          ppm		Levelled and filtered quadrature 3005 Hz
         8      14      P12lev         ppm		Levelled and filtered in-phase 11962 Hz
        12      15      Q12lev         ppm		Levelled and filtered quadrature 11962 Hz
         9      16      P25lev         ppm		Levelled and filtered in-phase 24510 Hz
        13      17      Q25lev         ppm		Levelled and filtered quadrature 24510 Hz
        14      18      PLM_nT         nT		Power line monitor

        """

        ncol = [3, 0, 1, 6, 7, 8,   10, 12, 14, 16, 11, 13, 15, 17,   18]

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[0].lower().startswith("#")
                    or line[0].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie" in line[:24].lower()):
                    continue

                t = line.split()

                t = [w.replace("*", invalid_strng) for w in t]
                tmp = [t[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        EPSG_in = 2157
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                           utm_zone_in=EPSG_in,
                                                           utm_zone_out=EPSG_out)

        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))

    if Survey.lower() in ["sglt2062", ]:
        """
        ===================================================================================================================
        Testlines SGL (Bundoran)    DLV 2062 (2019)  File: GSI_TESTLINE_FEM.xyz
        ===================================================================================================================
        GSNI__11.IRL, GSI___15.IRL, GSI___16.IRL, GSI___17.IRL AND GSI___18.IRL  (GSI_TESTLINE_FEM.xyz)

            AEM  column  Size        Name  Units    Description

            0    01     11          LINE    -       Line number -
                                                    LHHHDAV.YY
                                                    L = 1, 2 or 3 (from west to east),
                                                    H = nominal height,
                                                    D = direction, either 0 (for northbound or "forward")
                                                    or 1 (for southbound or "backward"),
                                                    A = attempt, starts at 1, increments for each new test flight,
                                                    V = version, starts at 0, increments for each pass,
                                                    Y = year}
                 02     6         FLI`GHT    -       Flight Number
                 03     10          DATE    -       Date YYYYMMDD
                 04     5            DAY    -       Day of year
                 05     11          TIME    s       Fiducial Seconds
            1    06     13         ITM-X    m       X coordinate, IRENET95 ITM
            2    07     13         ITM-Y    m       Y coordinate, IRENET95 ITM
                 08     13         UTM-X    m       X coordinate, WGS-84 UTM29N
                 09     13         UTM-Y    m       Y coordinate, WGS-84 UTM29N
                 10     13         UTM-Z    m       GPS Elevation (above WGS-84 Ellipsoid)
            3    11     13        MSLHGT    m       Mean Sea Level Altitude
                 12     16           LAT  degree    Latitude, WGS-84
                 13     16          LONG  degree    Longitude, WGS-84
            4    14     11         RADAR    m       Radar
                 15     11          TEMP   °C       Temperature
                 16     I11       P09ppm   ppm      In-phase 912 Hz
                 17     I11       Q09ppm   ppm      Quadrature 912 Hz
                 18     I11        P3ppm   ppm      In-phase 3005 Hz
                 19     I11        Q3ppm   ppm      Quadrature 3005 Hz
                 20     I11       P12ppm   ppm      In-phase 11962 Hz
                 21     I11       Q12ppm   ppm      Quadrature 11962 Hz
                 22     I11       P25ppm   ppm      In-phase 24510 Hz
                 23     I11       Q25ppm   ppm      Quadrature 24510 Hz
            6    24     I11      P09filt   ppm      Filtered in-phase 912 Hz
            10   25     I11      Q09filt   ppm      Filtered quadrature 912 Hz
            7    26     I11       P3filt   ppm      Filtered in-phase 3005 Hz
            11   27     I11       Q3filt   ppm      Filtered quadrature 3005 Hz
            8    28     I11      P12filt   ppm      Filtered in-phase 11962 Hz
            12   29     I11      Q12filt   ppm      Filtered quadrature 11962 Hz
            9    30     I11      P25filt   ppm      Filtered in-phase 24510 Hz
            13   31     I11      Q25filt   ppm      Filtered quadrature 24510 Hz
                 32     F11.3    912_res ohm-m      912Hz Apparent Resistivity
                 33     F11.3   3005_res ohm-m      3005Hz Apparent Resistivity
                 34     F11.3  11962_res ohm-m      11962Hz Apparent Resistivity
                 35     F11.3  24510_res ohm-m      24510Hz Apparent Resistivity
                 36     F11.2     PLM_mV    mV      Power line monitor
            14   37     F11.2     PLM_nT    nT      Magnetometer Derived Power line monitor



"""
        print(Survey.lower(), EPSG_in,EPSG_out)

        # filtered:
        ncol = [ 0,  5, 6,  10, 13,  23, 25, 27, 29, 24, 26, 28, 30, 36, ]

        # raw
        # ncol = [ 0,  5,  6, 10, 13, 15, 17, 19, 21, 16, 18, 20, 22, 36,]

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[:24].lower().startswith("#")
                    or line[:24].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie"  in line[:24].lower()):
                    continue


                t = line.split()



                tmp = [w.replace("*", invalid_strng) for w in t]
                tmp = [tmp[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        # EPSG_in = 2157
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                           utm_zone_in=EPSG_in,
                                                           utm_zone_out=EPSG_out)

        Data = numpy.insert(Data, 5, Data[:, 3]-Data[:, 4], axis=1)

        Data = numpy.where(Data<Invalid/1.e2,Data, Invalid)
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))


    if Survey.lower() in ["sglt2402", ]:
        """
        ===================================================================================================================
        Testlines SGL (Bundoran)    DLV 2062 (2019)  File: GSI_TESTLINE_FEM.xyz
        ===================================================================================================================
        GSNI__11.IRL, GSI___15.IRL, GSI___16.IRL, GSI___17.IRL AND GSI___18.IRL  (GSI_TESTLINE_FEM.xyz)

            AEM  column  Size        Name  Units    Description

            0    01     11          LINE    -       Line number -
                                                    LHHHDAV.YY
                                                    L = 1, 2 or 3 (from west to east),
                                                    H = nominal height,
                                                    D = direction, either 0 (for northbound or "forward")
                                                    or 1 (for southbound or "backward"),
                                                    A = attempt, starts at 1, increments for each new test flight,
                                                    V = version, starts at 0, increments for each pass,
                                                    Y = year}
                 02     6         FLI`GHT    -       Flight Number
                 03     10          DATE    -       Date YYYYMMDD
                 04     5            DAY    -       Day of year
                 05     11          TIME    s       Fiducial Seconds
            1    06     13         ITM-X    m       X coordinate, IRENET95 ITM
            2    07     13         ITM-Y    m       Y coordinate, IRENET95 ITM
                 08     13         UTM-X    m       X coordinate, WGS-84 UTM29N
                 09     13         UTM-Y    m       Y coordinate, WGS-84 UTM29N
                 10     13         UTM-Z    m       GPS Elevation (above WGS-84 Ellipsoid)
            3    11     13        MSLHGT    m       Mean Sea Level Altitude
                 12     16           LAT  degree    Latitude, WGS-84
                 13     16          LONG  degree    Longitude, WGS-84
            4    14     11         RADAR    m       Radar
                 15     11          TEMP   °C       Temperature
                 16     I11       P09ppm   ppm      In-phase 912 Hz
                 17     I11       Q09ppm   ppm      Quadrature 912 Hz
                 18     I11        P3ppm   ppm      In-phase 3005 Hz
                 19     I11        Q3ppm   ppm      Quadrature 3005 Hz
                 20     I11       P12ppm   ppm      In-phase 11962 Hz
                 21     I11       Q12ppm   ppm      Quadrature 11962 Hz
                 22     I11       P25ppm   ppm      In-phase 24510 Hz
                 23     I11       Q25ppm   ppm      Quadrature 24510 Hz
            6    24     I11      P09filt   ppm      Filtered in-phase 912 Hz
            10   25     I11      Q09filt   ppm      Filtered quadrature 912 Hz
            7    26     I11       P3filt   ppm      Filtered in-phase 3005 Hz
            11   27     I11       Q3filt   ppm      Filtered quadrature 3005 Hz
            8    28     I11      P12filt   ppm      Filtered in-phase 11962 Hz
            12   29     I11      Q12filt   ppm      Filtered quadrature 11962 Hz
            9    30     I11      P25filt   ppm      Filtered in-phase 24510 Hz
            13   31     I11      Q25filt   ppm      Filtered quadrature 24510 Hz
                 32     F11.3    912_res ohm-m      912Hz Apparent Resistivity
                 33     F11.3   3005_res ohm-m      3005Hz Apparent Resistivity
                 34     F11.3  11962_res ohm-m      11962Hz Apparent Resistivity
                 35     F11.3  24510_res ohm-m      24510Hz Apparent Resistivity
            14   36     F11.2     PLM_nT    nT      Magnetometer Derived Power line monitor



"""

        # filtered:
        ncol = [ 0,  5, 6,  10, 13,  23, 25, 27, 29, 24, 26, 28, 30, 35, ]

        # raw
        # ncol = [ 0,  5,  6, 10, 13, 15, 17, 19, 21, 16, 18, 20, 22, 35,]

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[:24].lower().startswith("#")
                    or line[:24].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie"  in line[:24].lower()):
                    continue


                t = line.split()



                tmp = [w.replace("*", invalid_strng) for w in t]
                tmp = [tmp[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        EPSG_in = EPSG_out
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                           utm_zone_in=EPSG_in,
                                                           utm_zone_out=EPSG_out)
        Data = numpy.insert(Data, 5, Data[:, 3]-Data[:, 4], axis=1)
        Data = numpy.where(Data<Invalid/1.e2,Data, Invalid)
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))


    if any(s in Survey.lower() for s in ["nm", ]):
        """
        ===================================================================================================================
        Survey NM (GENESIS TDEM system):
        ===================================================================================================================

        Translation table from .xyz files to internal data format:

        0   0       -           Line Number
            1       -           Flight Number
            2       yyyymmdd    Date of Survey Flight
            3       sec         Universal Time (Seconds Since Midnight)
            4       degrees     Point by Point Bearing
            5       degrees     Latitude in WGS84
            6       degrees     Longitude in WGS84
        1   7       m           Easting (X) in WGS84 UTM Zone 29N
        2   8       m           Northing (Y) in WGS84 UTM Zone 29N
            9       m           Easting (X) in IRENET95 Irish Transverse Mercato
            10      m           Northing (Y) in IRENET95 Irish Transverse Mercator
        3   11      m           GPS Elevation (Referenced to Mean Sea Level)
        4   12      m           Radar Altimeter
        5   13      m           Terrain (Referenced to Mean Sea Level)
            14      A           Transmitter Current
            15      ppm         X-Coil Late Time N294.51     194.9ormalization Channel
            16      ppm         Z-Coil Late Time Normalization Channel
            17-27   ppm         Raw X-Coil Channels 01 to 11
            28-38   ppm         Raw Z-Coil Channels 01 to 11
            39-49   ppm         Levelled X-Coil Channels 01 to 11
            50-60   ppm         Levelled Z-Coil Channels 01 to 11
            61-71   ppm         Levelled and Height Corrected X-Coil Channels 01 to 11
            18-28 72-82   ppm         Levelled and Height Corrected Z-Coil Channels 01 to 11

        ncol = [0, 7, 8, 11, 12, 13,
                    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
                    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82]

        """
        ncol = [0, 7, 8, 11, 12, 13,
                61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
                72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82]

        # print(numpy.shape(ncol))
        # Data = []
        # iline = 0
        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[:24].lower().startswith("#")
                    or line[:24].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie"  in line[:24].lower()
                    or "*" in line):
                    continue

                t = line.split()
                tmp = [w.replace("*", invalid_strng) for w in t]
                tmp = [tmp[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        EPSG_in = EPSG_out
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                       utm_zone_in=EPSG_in,
                                                       utm_zone_out=EPSG_out)



    if any(s in Survey.lower() for s in [ "nmx", ]):
        """
        ===================================================================================================================
        Survey NM (GENESIS TDEM system):
        ===================================================================================================================

        Translation table from .xyz files to internal data format:

        0   0       -           Line Number
            1       -           Flight Number
            2       yyyymmdd    Date of Survey Flight
            3       sec         Universal Time (Seconds Since Midnight)
            4       degrees     Point by Point Bearing
            5       degrees     Latitude in WGS84
            6       degrees     Longitude in WGS84
        1   7       m           Easting (X) in WGS84 UTM Zone 29N
        2   8       m           Northing (Y) in WGS84 UTM Zone 29N
            9       m           Easting (X) in IRENET95 Irish Transverse Mercato
            10      m           Northing (Y) in IRENET95 Irish Transverse Mercator
        3   11      m           GPS Elevation (Referenced to Mean Sea Level)
        4   12      m           Radar Altimeter
        5   13      m           Terrain (Referenced to Mean Sea Level)
            14      A           Transmitter Current
            15      ppm         X-Coil Late Time N294.51     194.9ormalization Channel
            16      ppm         Z-Coil Late Time Normalization Channel
            17-27   ppm         Raw X-Coil Channels 01 to 11
            28-38   ppm         Raw Z-Coil Channels 01 to 11
            39-49   ppm         Levelled X-Coil Channels 01 to 11
            50-60   ppm         Levelled Z-Coil Channels 01 to 11
            61-71   ppm         Levelled and Height Corrected X-Coil Channels 01 to 11
            18-28 72-82   ppm         Levelled and Height Corrected Z-Coil Channels 01 to 11

        ncol = [0, 7, 8, 11, 12, 13,
                    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
                    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82]

        """
        ncol = [0, 7, 8, 11, 12, 13,
                61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
                72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82]

        # print(numpy.shape(ncol))
        # Data = []
        # iline = 0
        Data = []
        Norm = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[:24].lower().startswith("#")
                    or line[:24].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie"  in line[:24].lower()
                    or "*" in line):
                    continue

                t = line.split()
                tmp = [w.replace("*", invalid_strng) for w in t]
                tmp = [tmp[ii] for ii in ncol]

                Data.append(tmp)
                Norm.append([line[59], line[60]])

        Data = numpy.asarray(Data, dtype=float)
        Norm = numpy.asarray(Norm, dtype=float)

        Data[:, 6:6+11] = Norm[:,0]*Data[:, 6:6+11]
        Data[:, 18:18+11] = Norm[:,1]*Data[:, 18:18+11]

        EPSG_in = EPSG_out
        if EPSG_in != EPSG_out:
            Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                       utm_zone_in=EPSG_in,
                                                       utm_zone_out=EPSG_out)


    if Survey.lower() in ["cggt", ]:
        """
        ===================================================================================================================
        Testlines GENESIS (Bundoran)
        ===================================================================================================================

        Translation table from .xyz files to internal data format:

        0   0       -           Line Number
            1       -           Flight Number
            2      m            Easting (X) in IRENET95 Irish Transverse Mercator
            3      m            Northing (Y) in IRENET95 Irish Transverse Mercator
        3   4      m            GPS Elevation (Referenced to Mean Sea Level)
        4   5      m            Radar Altimeter
        5   6      m            Terrain (Referenced to Mean Sea Level)
            10-21               Levelled and Height Corrected X-Coil Channels 01 to 11
            22-33               Levelled and Height Corrected Z-Coil Channels 01 to 11



        """

        ncol = [0, 1, 2, 3, 4, 5,
                7,  8,  9,  10, 11, 12, 13, 14, 15, 16,
                17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]

        Data = []
        iline = 0
        with open(DatFile) as fd:
            for line in fd:
                iline = iline + 1
                if (line[:24].lower().startswith("#")
                    or line[:24].lower().startswith("/")
                    or "line" in line[:24].lower()
                    or "tie"  in line[:24].lower()
                    or "*" in line):
                    continue


                t = line.split()



                tmp = [w.replace("*", invalid_strng) for w in t]
                tmp = [tmp[ii] for ii in ncol]
                Data.append(tmp)

        Data = numpy.asarray(Data, dtype=float)

        Data[:,1], Data[:,2] = util.project_utm_to_utm(Data[:,1], Data[:,2],
                                                       utm_zone_in=2157,
                                                       utm_zone_out=32629)

        Data1 = Data[:,:16]
        Data2 = Data[:,16:]
        Data = numpy.column_stack((Data1,numpy.nan*numpy.ones_like(Data[:,0])))
        Data = numpy.column_stack((Data, Data2))
        Data = numpy.where(Data<Invalid/1.e2, Data, Invalid)
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))
        Data = numpy.column_stack((Data,numpy.zeros_like(Data[:,0])))

    nn = numpy.shape(Data)

    if numpy.size(Data)==0:
        error("Survey format for " + Survey + " not found!.Exit")

    if OutInfo:
        print("readDat: %i x %i data read from %s" % (nn[0], nn[1], DatFile))

    return Data

def get_header(file, headstr="# ", OutInfo=False):
    """
    Get header lines from file.

    VR 11/20

    """
    if OutInfo:
        print("\nHeader:")
    header = []
    with open(file) as fd:
        for line in fd:
            if line.startswith(headstr) and not "line" in line:
                # line=line[2:]
                header.append(line.replace(headstr, ""))

                if OutInfo:
                    print(line[1:-1])
    if OutInfo:
        print("\n")
        print(header)
    Header = header
    return Header


def grow_header(Header=None, Addstr=None, OutInfo=False):
    """
    Grow header string.

    author: vrath
    last changed: 2022/11/4
    """
    if OutInfo:
        print("\nOld Header:")
        print(Header)

        print("New lines:")
        print(Addstr)

    if type(Header) == str:
            Header=[Header]
    if Addstr is None:
        pass
    else:
        Header.append(Addstr)
            # Header = Header+line

    # workaround for old files only...
    try:
        Header = " | ".join(Header)
    except:
        Header =str(Header).replace(", ", " | ")
        Header = Header.replace("'", " ")
        Header = Header.replace("[", "")
        Header = Header.replace("]", "")

    if OutInfo:
        print("\nNew Header:")
        # print(" \n".join(Header))
        print(Header)

    return Header


def print_header(Header):
    """
    Grow header string.

    author: vrath
    last changed: 2020/12/4
    """
    print("Current Header:")
    # print("\n".join(Header))
    print(Header)


def read_aempy(File=None, Format=None, System="aem05", OutInfo=False):
    """
    Read Tellus Data.

    author: vrath
    last changed:  nov 2021

    """
    if File == None:
        error("No filename given! Exit.")

    if Format == None:
        name, ext =  os.path.splitext(File)
        ff = ext
    else:
        ff = Format.lower()


    if   "npz" in ff:
        Data, Header, System = read_aempy_npz(File, OutInfo=OutInfo)
    elif "asc" in ff:
        Data, Header, System = read_aempy_asc(File, OutInfo=OutInfo)
    elif "nc4" in ff:
        Data, Header, _System= read_aempy_ncd(File, OutInfo=OutInfo)
    else:
        error("Input format " + ff + " not implemented! Exit.")

    return Data, Header, System

def write_aempy(File=None, Data=None, Format=None, System="aem05",
                Header="", OutInfo=False):
    """
    Write Tellus Data.

    author: vrath
    last changed:  sep 2021
    """
    if File == None:
        error("No filename given! Exit.")

    if Format == None:
        name, ext =  os.path.splitext(File)
        ff = ext
    else:
        ff = "."+Format.lower()
        File= File+ff

    print("Output file format is "+ff)

    if  "npz" in ff:
        write_aempy_npz(File , Data=Data, Header=Header, OutInfo=OutInfo)
        if OutInfo:
            print("Data written to File: "+File)

    elif "asc" in ff:
        write_aempy_asc(File , Data=Data, System=System,
                           Header=Header, OutInfo=OutInfo)
        if OutInfo:
            print("Data written to File: "+File)

    elif "nc4" in ff:
        write_aempy_ncd(File, Data=Data,
                            Header=Header, OutInfo=OutInfo)
        if OutInfo:
            print("Data written to File: "+File)

    else:
        error("Output format " + ff + " not implemented! Exit.")


def write_aempy_asc(File=None, Data=None,
                    System="aem05", Header="", OutInfo=False):
    """
    Write Tellus AEM05 Data to ASCII file.

    author: vrath
    last changed: 2021/02/22
    """
    if os.path.exists(File):
        os.remove(File)

    _, nn, outfmt, columns, _ = get_system_params(System)
    head = grow_header(Header, " |"+System+"| "+columns)
    # head = "".join(head).replace("|", "\n")
    print(outfmt)
    print(numpy.shape(Data))
    numpy.savetxt(File, Data, header=head, fmt=outfmt)


def read_aempy_asc(File=None, System="unknown", OutInfo=False):
    """
    Read Tellus  Data to ASCII file.

    author: vrath
    last changed: 2020/11/15
    """
    # print(File)

    if not os.path.isfile(File):
        error("File %s does not exist! Exit" % File)

    Header = get_header(File, headstr="# ", OutInfo=OutInfo)

    if "aem05" in Header:
        System = "aem05"
    if "gen" in Header:
        System = "genesis"

    Data = numpy.loadtxt(File, comments="#")
    # print(numpy.shape(Data))

    return Data, Header, System


def write_aempy_npz(File=None, Data=None, Header="", System="aem05", OutInfo=False):
    """
    Write Tellus AEM05 Data to compressed pickle file.

    author: vrath
    last changed: 2020/11/15
    """
    if os.path.exists(File):
        os.remove(File)

    Head = numpy.array(Header)
    numpy.savez_compressed(File, Data=Data, Header=Head, System=System)
    # print(numpy.shape(Data))
    if OutInfo:
        print("Data written to " + File)
        print(numpy.shape(Data))


def read_aempy_npz(File=None, OutInfo=False):
    """
    Read Data from compressed pickle file.

    author: vrath
    last changed: 2020/11/15
    """
    if not os.path.isfile(File):
        error("File %s does not exist! Exit" % File)

    with numpy.load(File, allow_pickle=True) as data:
        print(data.files)
        Data = data["Data"]
        Header = data["Header"]
        System = data["System"]




    if isinstance(Header,numpy.ndarray):
        Header = Header.tolist()

    if OutInfo:
        print("Data read from " + File)
        print("Header:")
        print(Header)
    # print(numpy.shape(Header))
    # print(type(Header))

    return Data, Header, System


def write_aempy_ncd(File=None, Data=None, Header="",System="aem05",
                    zlib_in=True, shuffle_in=True, OutInfo=False):
    """
    Write data to NETCDF4 file.

    author: vrath
    last changed: 2020/11/15
    """
    # try:
    #     File.close
    # except:
    #     pass

    if os.path.exists(File):
        os.remove(File)

    # if os.path.isfile(File):

    #     try:
    #         os.remove(File)
    #         print("File %s has been removed successfully" % File)
    #     except BaseException:
    #         error("File %s can not be removed" % File)

    # _, columns = getTellusfmts()
    Head = Header
    DataDim = numpy.shape(Data)
    Head = numpy.array(Header, dtype=object)
    HeadDim = numpy.shape(Head)

    ncout = nc.Dataset(File, "w", format="NETCDF4")

    now = datetime.now()
    ncout.description = "\n".join(Header)
    ncout.history = "Created " + now.strftime("%m/%d/%Y, %H:%M:%S")

    ncout.createDimension("data", DataDim[0])
    ncout.createDimension("obs", 8)
    ncout.createDimension("head", HeadDim[0])

    h = ncout.createVariable("header", "str", "head")
    s = ncout.createVariable("system", "str", 1)

    f = ncout.createVariable(
        "line",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    f.units = "none"
    x = ncout.createVariable(
        "x",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    x.units = "m"
    y = ncout.createVariable(
        "y",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    y.units = "m"
    g = ncout.createVariable(
        "gps",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    g.units = "m"
    c = ncout.createVariable(
        "clr",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    c.units = "m"
    e = ncout.createVariable(
        "dem",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    e.units = "m"

    p = ncout.createVariable(
        "pli",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    p.units = "none"

    q = ncout.createVariable(
        "qfl",
        "float64",
        "data",
        zlib=zlib_in,
        shuffle=shuffle_in)
    q.units = "none"

    d = ncout.createVariable(
        "data",
        "float64",
        ("data",
         "obs"),
        zlib=zlib_in,
        shuffle=shuffle_in)
    d.units = "ppm"

    f[:] = Data[:, 0]
    x[:] = Data[:, 1]
    y[:] = Data[:, 2]
    g[:] = Data[:, 3]
    c[:] = Data[:, 4]
    e[:] = Data[:, 5]
    p[:] = Data[:, 14]
    q[:] = Data[:, 15]

    d[:, :] = Data[:, 6:14]

    h[:] = Head[:]
    s[:] = System

    ncout.close()

    if OutInfo:
        print("Data written to %s in %s format" % (File, ncout.data_model))


def read_aempy_ncd(File=None, Data=None, Split=False, OutInfo=False):
    """
    Read data from NETCDF4 file.

    author: vrath
    last changed: 2021/02/15
    """
    if not os.path.isfile(File):
        error("File %s does not exist! Exit" % File)

    ncin = nc.Dataset(File, "r", format="NETCDF4")
    print(ncin.file_format)

    f = ncin.variables("line")[:]
    x = ncin.variables("x")
    y = ncin.variables("y")
    g = ncin.variables("gps")
    c = ncin.variables("clr")
    e = ncin.variables("dem")
    p = ncin.variables("pli")
    q = ncin.variables("qfl")
    d = ncin.variables("data")

    h = ncin.variables("header")
    s = ncin.variables("system")

    System = s

    Header = h

    Data = numpy.column_stack((f, x, y, g, c, e, d, p))

    if OutInfo:
        print(Header)
        print(numpy.shape(Data))

    ncin.close()

    if OutInfo:
        print("Data read from %s in %s format" % (File, ncin.data_model))

    if Split:
        return f, x, y, g, c, e, p, q, d, h, s
    else:
        return Data, Header, System


def merge_data_files(File_list=[],
                     Merged="merged_data", MergedHeader=None,
                     AEM_system="aem05", OutInfo=False):
    """


    Parameters
    ----------
    f_list : list of strings
        List of data files (npz), with full paths.
    merged : string, optional
        Merged data files. The default is merged.npz.

    Returns
    -------
    None.

    """
    nf = len (File_list)
    if nf==0:
        error("merge_data_files: file list is empty! Exit.")



    for ifile in numpy.arange(nf):
        filein = File_list[ifile]
        print("\n Reading file " + filein)
        data_in, header, system = read_aempy_npz(File=filein, OutInfo=False)
        # print(numpy.shape(data_in))
        if ifile==0:
            system_out = system
            data_out = data_in.copy()
            # print(numpy.shape(data_out))
        else:
            data_out = numpy.vstack((data_out, data_in))
            # print(numpy.shape(data_out))

    print("\n\n", numpy.shape(data_out), " data written to ", Merged )

    write_aempy_npz(File=Merged, Data=data_out, Header=MergedHeader,
                        System=system_out , OutInfo=False)

    return data_out
