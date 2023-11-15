!*==main.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      PROGRAM airBeo
!------------
 
!*** Calls DCPRM_FD, DCPRM_TD, HSBOSS_TD, HSBOSS_FD, NLSQ2,
!          SET_NORM_TD, SET_SOURCE, WRITE_FD, WRITE_TD, WRITE_MODEL,
!          WRITE_LOG_FILE
 
!*** Calls from INPUT_DATA_FOR_AIRBEO:
!          READ_SYSTEM_AND_SURVEY, READ_MODEL, READ_INVRT_CNTRL_AND_DATA
!          SET_SURVEY, SET_TRP, WRITE_NW1_INITIAL
 
 !   PROGRAM AirBeo       Version 4.7.0    16 March 2007       Language: ANSI Standard Fortran 95
!   ======= ======       =============    =============       ----------------------------------
!
!
!            Developed by:  Art Raiche
!                     For:  AMIRA project P223F
!
!
!================
!   DESCRIPTION |
!================
!
!  AirBeo is an AEM inversion program which fits a layered-earth model to
!  frequency or time-domain data.  The flight line is descretised into
!  NSTAT stations.  A separate layered earth model computed for each station.
!
!  In frequency domain, only the maximally coupled component is fit.  This
!  can be mixed horizontal coplanar, vertical co-axial or vertical co-planar
!  or a new option, vertical coplanar broadside.  The data are assumed to be
!  normalised; ie, ppm, ppb ppt or percent.
!
!  In time domain, the user can invert for the vertical (Z), component only,
!  the horizontal in-line (X) componentant only, both components or the total
!  component.  These can be expressed as dB/dt, magnetic field, or normalised
!  data.
!
!  The NLYR model consists of NLYR-1 layers over basement.  The user can
!  express the model in terms of layer thicknesses or depth to the base of
!  each layer.  The latter is best used when depth to basement is important.
!
!  The layer lithologies are composed of resistivities, relative dielectric
!  constants, relative magnetic permeabilities and the Cole-Cole parameters.
!
!  This version of AirBeo inverts for resistivities & either layer thicknesses
!  or depth to basement of each layer.  Lithology components other than
!  resistivity do not vary during inversion.
!
!*********************
!      Systems
!*********************
!
!  AirBeo can be used to model any existing AEM system in frequency or time-
!  domain mode.  Currently the transmitter is modelled as a magnetic dipole
!  whose axis is in the vertical plane along the flight path except for the
!  new vertical coplanar broadside option.  The other exception to this is
!  for co-axial time-domain HEM where the system normalisation assumes a
!  horizontal loop of finite radius.
!
!   Time-domain waveform sign convention
!   ------------------------------------
!
!  If the user specifies waveform excitation as either the transmitter
!  current or primary B at the receiver, the program assumes that current or
!  magnetic field will start at 0 or some low value, and rise to a positive
!  maximum before going to zero or oscillating about zero with small
!  magnitudes.  In this case, dI/dt will be computed using the negative I or B
!  so that the early off-time response is positive for vertical fields.
!
!  If the user specifies the excitation waveform as dB/dt at the receiver,
!  the program assumes that the response will rise to a positive maximum
!  followed by a negative maximum before oscillating about zero with smaller
!  amplitudes.  In this case dI/dt is derived by reversing the sign of the
!  input primary dB/dt and dividing by the geometric coupling factor.  Again,
!  this procedure is designed so that the early off-time response is positive
!  for vertical fields.
!
!*********************
!   FILE CONVENTIONS
!*********************
!
!   INPUT FILES:
!   -----------
!
!   The input control file, named AirBeo.cfl is read
!   from logical unit number NR = 3.
!
!   For inversion the data to be inverted must be contained in
!   AirBeo.inv, read from logical unit NRI = 13
!
!
!   VERBOSE-OUTPUT FILES:
!   --------------------
!
!   The AirBeo.out is written to logical unit NW = 4
!
!   Messages about data or runtime errors are written in a file
!   called AirBeo.log on logical unit NLG = 9
!
!
!   OUTPUT FILES FOR PLOTTING:
!   --------------------------
!
!   Terse inversion output for plotting is written to
!   logical unit NI = 14
!
!   NMD = 15 contains a very terse file containing resistivities & thicknesses
!            as a function of position.  The resistivities for all layers,
!            starting from the top to bottom are written out followed by
!            the thicknesses starting from the top to bottom.
!
!
!    UNIT #   UNIT ID      FILE ID      Function
!    ------   -------      -------      --------
!       3       NR       AirBeo.cfl   Input control file
!       4       NW       AirBeo.out   Verbose output data
!       9       NLG      AirBeo.log   Data error messages
!      13       NRI      AirBeo.inv   Data to be inverted
!      14       NW1      AirBeo.mv1   Tabular inversion output
!      15       MD1      AirBeo.mdl   Resistivities + layer depths & thicknesses
!
!-------------------------------------------------------------------------------
!
!****************************************************
!     DESCRIPTION OF DATA RECORDS for AIrBeo.CFL
!****************************************************
!
!  All records are in list directed (free) format except for
!  the TITLE in RECORD 1.
!
!    NOTE:  All distances and locations are to be specified in metres.
!    ----   All angles are to be specified in degrees.
!
!      In plan view, the coordinate system used for input and output
!      is (East, North, Depth).   Depth is positive downwards.
!
!
!**  RECORD 1:    TITLE - up to 120 characters
!
!
!             MODEL CONTROL, UNITS & PRINT PARAMETERS
!             ---------------------------------------
!
!**  RECORD 2:  TDFD, DO1D, PRFL, ISTOP
!
!
!      TDFD = 1 => time-domain modelling
!           = 2 => frequency-domain modelling
!
!      DO1D = 1 or -1  Inversion using the same starting model for each station
!           = 2 or -2  Inversion using the final model each station as the
!                      model for the following station. (Care required.)
!           = 0 =>  Forward modelling only
!
!
!      PRFL =  1 prints response in profile mode using the default units and
!                line tag.  Each column contains the response at all stations
!                for a channel or frequency.
!
!           =  0 prints responses in temporal or frequency mode, using the
!                default units and line tag.  Each column contains the
!                responses for a single position for all frequencies (TDFD=2)
!                or delay times (TDFD=1).
!
!           = -1 same as PRFL = 1 except for new feature activation
!
!     ISTOP = 0  read the input data and run the specified models.
!           = 1  read the input data, print the model description
!                and STOP so that the model description can be verified.
!
!            REMEMBER to change ISTOP to 0 once the models have been verified.
!
!
!  RECORDS 3 & 4 for FREQUENCY-DOMAIN AEM SYSTEM INFORMATION
!  -------------------------------------------------------------
!
!**  RECORD 3:  NFREQ, CMP, NPPF
!
!      NFREQ - number of frequencies
!
!    INVERSION:   data must be in normalised form either in PPM, PPB, PPT or percent
!    ---------    Invert on single component response in direction of
!                 transmitter orientation; eg co-planar or co-axial.
!
!      CMP = 1 => Transmitter dipole axis azimuth is oriented along flight path
!                 Transmitter dipole axis inclinations are specified in RECORD 5
!
!      CMP = -1 for VCPB (vertical co-planar broadside array)
!               Transmitter dipole axis azimuth is 90 degrees to flight line
!               Transmitter dipole axis inclination is fixed at 90 degrees
!               regardless of RECORD 5
!
!
!    FORWARD MODELLING
!    -----------------
!      CMP = 1 or -1 => compute single component response in direction of
!                       transmitter orientation; eg co-planar or co-axial.
!                       Use CMP = -1 only for VCBB array
!          = 2 => compute vertical and horizontal in-line response.
!          = 3 => compute 3 component response.
!
!           If CMP = 1, output is in normalised units specified by NPPF
!           If CMP > 1, output is in picoTeslas
!
!
!     NPPF = 1 => normalisation is in percent (parts per hundred)
!          = 2 => normalisation in parts per thousand
!          = 3 => normalisation in parts per million
!          = 4 => normalisation in parts per billion
!
!
!**  RECORDS 4: (FREQ(J), ZRX(J), XRX(J), YRX(J), TXCLN(J), J = 1, NFRQ)
!            (enter a separate record for each of the NFRQ offsets in metres
!
!      FREQ(J) - frequency in Hz
!
!      ZRX(J) - vertical Rx offset : positive if Rx is below Tx
!      XRX(J) - in-line Rx offset :  positive if Rx is behind Tx
!      YRX(J) - transverse offset :  positive if Rx is left of Tx
!
!      TXCLN(J) - inclination angle in degrees that transmitter dipole axis makes
!                 with the vertical.
!                 TXCLN = 0 for a horizontal loop in level flight.
!                 TXCLN is positive if the front of the loop is
!                       above the rear of the loop
!
!      If CMP = -1, the default TXCLN = 90.  In this case, entries for
!      TXCLN in RECORD 4 are ignored and needn't be specified.
!
!     DATA ENTRY CONTINUES WITH FLIGHT PATH SPECIFICATION IN RECORD 9
!
!  _____________________________________________________________________
!
!  RECORDS 3, 4, & 5 for TIME-DOMAIN AEM SYSTEM INFORMATION  (ISW > 0)
!  ---------------------------------------------------------------------
!
!   The user will only specify the waveform for 1/2 cycle.  The program
!   will then fill in the other half as the negative of the first half.
!   Any non-negative off-time value can be specified.
!
!**  RECORD 3:  ISW, NSX, STEP, UNITS, NCHNL, KRXW, OFFTIME
!
!      ISW =  1 =>  RECORD 4 will contain a transmitter current waveform in amps
!      ISW = -1 =>  Aerotem triangular current input with on-time channels
!
!      ISW =  4 =>  Step B system where the Tx current waveform is a
!                   bipolar square wave (no off-time) whose frequency
!                   will be read in RECORD 4.  Output is B in femtoteslas.
!
!                   For ISW = 4, set OFFTIME = 0  and NSX = 1
!
!                   Towed Bird ONLY options
!                   ------------------------
!       ISW = 10  => RECORD 4 will contain the horizontal IN-LINE dB/dt receiver
!                    calibration waveform in dB/dt UNITS
!
!       ISW = 11 =>  RECORD 4 will contain the horizontal IN-LINE B receiver
!                    calibration waveform in B UNITS.
!
!       ISW = 30 =>  RECORD 4 will contain the VERTICAL dB/dt receiver
!                    calibration waveform in dB/dt UNITS
!
!       ISW = 31 =>  RECORD 4 will contain the VERTICAL B receiver
!                    calibration waveform in B UNITS.
!
!                   Central loop ONLY options
!                   -------------------------
!       ISW = 130 =>  RECORD 4 will contain the VERTICAL dB/dt receiver
!                     calibration waveform in dB/dt UNITS
!       ISW = 131 =>  RECORD 4 will contain the VERTICAL B receiver
!                     calibration waveform in B UNITS.
!
!            Geotem / Questem Stripping Option
!            ---------------------------------
!
!            Geotem & Questem data are processed by using an algorithm to strip
!            primary field and bird motion from the output waveform.  Users
!            have the option to apply this algorithm to AirBeo model output.
!            This can be done by specifying ISW as a negative integer for
!            receiver waveform options.
!
!            In other words, setting ISW = -10, -11, -30 or -31 instead of
!            10, 11, 30 or 31 respectively will produce stripped output.
!
!
!      STEP = O =>  Compute dB/dt in nanovolts per unit area  (nV/m^2)
!                    (same as nanoteslas per second (nT/s) )
!
!           = 1 =>  Compute B field in picoteslas (pT)
!
!      UNITS apply to  - waveform input except for ISW = 1
!                      - un-normalised forward time-domain model output
!                      - un-normalised time-domain data to be inverted
!
!                    STEP=0   STEP=1
!                    dB/dt      B
!                    -----    ------
!      UNITS = 1     nT/s       nT
!              2     pT/s       pT
!              3     fT/s       fT
!
!      NCHNL - number of receiver channels
!
!      NSX - number of digitised points in waveform including endpoints
!
!      KRXW = 1 => receiver channels will be read in terms of start and end
!                  times in ms.
!      KRXW = 2 => receiver channels will be read in terms of midpoints and
!                  channel widths.
!
!      OFFTIME - time (milliseconds) between end of one pulse and the start of
!                the next pulse (of opposite sign) since a bipolar waveform is
!                assumed.  For systems which have a signal which is always on,
!                OFFTIME = 0.
!
!
!**  RECORD 4  for ISW /= 4:  (TXON(J), WAVEFORM(J), J = 1,NSX)
!              -------------
!
!      TXON(J) = digitised time (in milliseconds)
!                In most cases, TXON(1) = 0, TXON(NSX) = pulse on-time
!
!      WAVEFORM(J) = transmitter current (in amps) at time TXON(J) if ISW = 1
!                  = vertical dB/dt receiver waveform if ISW = 30 or 130
!                  = vertical B receiver waveform if ISW = 31 or 131
!                  = horizontal in-line dB/dt receiver waveform if ISW = 10
!                  = horizontal in-line B receiver waveform if ISW = 11
!
!**  RECORD 4  for ISW = 4 ONLY:  FREQ, TXAMPL
!              ----------------
!      FREQ = frequency of bipolar square wave
!      TXAMPL = peak to peak current amplitude
!
!
!         In RECORDS 5 (& 6), the receiver windows can be specified in terms of
!         start & end times if KRXW = 1; or centres & widths if KRXW = 2.
!
!   If KRXW = 1, Specify start and end times (in ms) of receiver windows.
!                   (measured from start of signal turn on.)
!   ---------------------------------------------------------------------
!**  RECORD 5:  (TOPN(J), TCLS(J), J=1, NCHNL)
!
!         No RECORD 6 if KRXW = 2!  Go to RECORD 7
!   ---------------------------------------------------------------------
!
!   If KRXW = 2,  Specify centres and widths of receiver windows.
!   ------------------------------------------------------------
!**  RECORD 5:  (TMS(J), J=1,NCHNL) - centre of receiver window gates in ms.
!                                     measured from start of signal turn on.
!
!**  RECORD 6:  (WIDTH(J), J=1,NCHNL)  -  width of Rx window gates in ms.
!
!**  RECORD 7:  TXCLN, CMP, KPPM
!
!         TXCLN - angle in degrees that transmitter dipole axis makes
!                 with the vertical for level flight.
!                 TXDIP = 0 for a horizontal loop in level flight.
!                 TXCLN > 0 if the front of the loop is above the rear of the loop
!
!
!    INVERSION:
!    ---------
!         CMP = 11 => invert on horizontal in-line component only
!             = 13 => invert on vertical component only
!             = 2  => joint inversion on vertical & horizontal in-line components
!             = 3  => invert on all 3 components
!             = 4  => invert on total field
!
!    FORWARD MODELLING
!    -----------------
!       CMP = 11 => print horizontal in-line component only
!           = 13 => print vertical component only
!           = 2  => print vertical & horizontal in-line components
!           = 3  => print all three components
!
!      KPPM = 0  => No normalisation (automatic if ISW = 4)
!
!      KPPM > 0 => normalisation based on maximum dI/dt response for dB/dt output.
!                  In this case, an entry for KPPF in RECORD 7.1 is necessary.
!                  ---------------------------------------------------------------
!
!      KPPM = 1   => All components are normalised to in-line primary field
!      KPPM = 3   => All components are normalised to vertical primary field
!      KPPM = 123 => Vertical component normalised to the vertical primary &
!                    In-line component to the in-line primary
!                    Transverse component to the transverse primary
!                     (For null coupling, total field is used.)
!
!      KPPM = 4 =>   all components are normalised to total primary field
!
!      -------------------------
!      only if KPPM > 0
!**    RECORD 7.1:  NPPF
!      -------------------------
!
!           NPPF = 1 => normalisation is in percent (parts per hundred)
!                = 2 => normalisation in parts per thousand
!                = 3 => normalisation in parts per million
!                = 4 => normalisation in parts per billion
!
!      -------------------------
!      only if ISW = 1 or ISW > 100
!**    RECORD 7.2:  TXAREA, NTRN
!      -------------------------
!
!      TXAREA - transmitter loop area (sq. metres)
!        NTRN - number of transmitter loop turns
!
!      TXAREA & NTRN need be specified ONLY for ISW = 1 or ISW > 100,
!                    because the program ignores them for all other ISW options.
!
!**  RECORDS 8: ZRX0, XRX0, YRX0
!
!      ZRX0 - initial vertical Rx offset : positive if Rx is below Tx
!      XRX0 - initial in-line Rx offset :  positive if Rx is behind Tx
!      YRX0 - initial transverse offset :  positive if Rx is left of Tx
!
!             These are the high altitude calibration offsets.
!             If SURVEY < 3 in RECORD 9, these offsets are used for each station
!
!             Flight Path Information
!             -----------------------
!
!  ======================================================================
!  RECORD 9.1 and 9.2 are used only for forward modelling, DO1D = 0
!
!  For inversion, make a single entry for RECORD 9.0 and go on to RECORD 10
!
!========================================================================
!----------------------------------------------------------------------
!           only for inversion - NPASS should be set to 0 or 1
!           NPASS > 1 will the cause program to ignore the next NPASS records
!
!**  RECORD 9.0: NPASS
!                      GO TO RECORD 10
!----------------------------------------------------------------------
!======================================================================
!
!    RECORDS 9 FOR FORWARD MODELLING (DO1D = 0)
!    ------------------------------------------
!
!**  RECORD 9.0:  NSTAT, SURVEY, BAROMTRC, LINE_TAG
!
!        NSTAT - total number of transmitter data positions for all lines
!
!       SURVEY = 1 constant altitude, course & Tx-Rx geometry
!              = 2 variable altitude & course but constant Tx-Rx geometry
!              = 3 variable altitude, course, transmitter pitch
!                          and receiver offset (time-domain only)
!
!     BAROMTRC = 0 => altitudes are ground clearance in metres.
!              = 1 => altitudes are barometric; ie, metres above sea level.
!
!     LINE_TAG = 0 => default line number (1000) will be used for all stations
!              = 1 => a line number must be specified for each station
!
!-------------------------------------------------------------------------------
!   ONLY IF SURVEY = 1    Constant course and altitude
!   ------------------
!
!    IF LINE_TAG = 0
!**  RECORD 9.1:  EAST(1), NORTH(1), ALT(1), BEARING, DSTAT
!
!    IF LINE_TAG = 1
!**  RECORD 9.1:  LINE(1), EAST(1), NORTH(1), ALT(1), BEARING, DSTAT
!
!      LINE(1) - integer line number
!      EAST(1) - Initial East coordinate of transmitter
!     NORTH(1) - Initial North coordinate of transmitter
! ---------------------------------------------------------------------------------
!         In frequency domain, these are interpreted as the Tx-Rx midpoint.
!         The Tx position for all frequencies are based on the offset for
!         the first frequency.  The Rx position for each frequency is computed
!         using the offset for that frequency relative to the Tx position for
!         frequency 1.
! ---------------------------------------------------------------------------------
!
!       ALT(1) - Altitude of transmitter
!      BEARING - Flight path angle with respect to North (in degrees)
!                Due North = 0;  due East = 90.
!        DSTAT - Distance (metres) between successive transmitter positions.
!
!-------------------------------------------------------------------------------
!
!   IF SURVEY > 1    Variable course and altitude  (NSTAT records to follow)
!   -------------
!
!        Enter  RECORD 9.J,  J = 1,NSTAT  - ie one for each station
!
!   If SURVEY = 2  and  LINE_TAG = 0
!   ---------------------------------
!**  RECORD 9.J:  EAST(J), NORTH(J), ALT(J)
!
!
!   If SURVEY = 3  and  LINE_TAG = 0    (Time-domain only)
!   ---------------------------------
!**  RECORD 9.J:  EAST(J), NORTH(J), ALT(J), TXCLN(J), ZRX(J), XRX(J), YRX(J)
!
!
!   If SURVEY = 2  and  LINE_TAG = 1
!   ---------------------------------
!**  RECORD 9.J:  LINE(J), EAST(J), NORTH(J), ALT(J)
!
!
!   If SURVEY = 3  and  LINE_TAG = 1    (Time-domain only)
!   ---------------------------------
!**  RECORD 9.J:  LINE(J), EAST(J), NORTH(J), ALT(J), TXCLN(J), ZRX(J), XRX(J), YRX(J)
!
!
!      LINE(1) - integer line number
!      EAST(J) - East coordinate of transmitter at station J
!     NORTH(J) - North coordinate of transmitter at station J
! ---------------------------------------------------------------------------------
!         In frequency domain, these are interpreted as the Tx-Rx midpoint.
!         The Tx position for all frequencies are based on the offset for
!         the first frequency.  The Rx position for each frequency is computed
!         using the offset for that frequency relative to the Tx position for
!         frequency 1.
! ---------------------------------------------------------------------------------
!
!       ALT(J) - Altitude of transmitter at station J
!
!       ZRX(J) - Vertical offset of receiver at station J
!                positive if Rx is below Tx
!       XRX(J) - In-line offset of receiver at station J
!                positive if Rx is behing Tx
!       YRX(J) - Transverse offset of receiver at station J
!                positive if Rx is left of Tx

!
!          LITHOLOGY & STRUCTURE FOR AIrBeo
!          ================================
!
!** RECORD 10:  NLYR, QLYR, NLITH, GND_LVL
!
!          NLYR - number of layers including basement.
!
!          QLYR = 1 => structure specified using thickness of each layer.
!
!         NLITH - number of layer lithologies.
!                 Any number of lithologies may be defined.
!
!       GND_LVL - Relative level of flat surface (m)
!
!
!          DEFINE LITHOLOGIES
!          ------------------
!
!** RECORD 11.1: RES(1), SIG_T(1), RMU(1), REPS(1), CHRG(1), CTAU(1), CFREQ(1)
!** RECORD 11.2: RES(2), SIG_T(2), RMU(2), REPS(2), CHRG(2), CTAU(2), CFREQ(2)
!     .
!
!** RECORD 11.N: RES(N), SIG_T(N), RMU(N), REPS(N), CHRG(N), CTAU(N), CFREQ(N)
!
!           N = NLITH
!      RES(I) - layer resistivity
!    SIG_T(I) - Conductance (conductivity-thickness product)
!               (not used but must be entered)
!      RMU(I) - relative layer magnetic permeability for LITH_INDEX(I)
!     REPS(I) - relative layer dielectric constant (permittivity for LITH_INDEX(I)
!     CHRG(I) - Cole-Cole layer chargeability for LITH_INDEX(I)
!     CTAU(I) - Cole-Cole layer time constant for LITH_INDEX(I)
!    CFREQ(I) - Cole-Cole layer frequency constant for LITH_INDEX(I)
!
!    Default values:  RMU = 1   REPS = 1   CHRG = 0   CTAU = 0   CFREQ = 1
!
!    The default means no magnetic permeability contrast (MU = 4 PI * 10^(-7))
!                      no dielectric constant contrast  (EPSILON = 8.854215E-12)
!                      and no IP effects (no Cole-Cole)
!
!          LAYERED EARTH STRUCTURE
!          -----------------------
!
!** RECORD 12.1: LITH(1), THK(1)
!
!** RECORD 12.J: LITH(J), THK(J)
!
!** RECORD 12.NLYR: LITH (NLYR) - basement lithology
!
!      LITH(J) = integer which assigns the resistivity and other
!                physical properties from the list of RECORD 11
!                to layer J.
!
!      THK(J) = thickness of layer J
!
!___________________________________________
!
!     END OF DATA ENTRY IF NO INVERSION
!___________________________________________
!========================================================
!
!             INVERSION CONTROL
!             =================
!
!    The inversion will run until one of the following occurs:
!
!    1. The maximum number of iterations, specified by MAXITS in
!       RECORD 16 has been completed.  Hopefully this won't
!       happen since the sugested value is MAXITS = 90
!
!    2. The RMS error is below that specified by PCTCNV.
!       Default = 1 percent
!
!    3. The inversion is unable to reduce the error any further.
!
!    4. The predicted error decrease goes to zero
!
!    5. After 10 iterations, the combined error reduction of the previous
!       two iterations is less than 1 percent.
!
!    In view of these criteria, it is best not to set MAXITS too low.  During
!    developmental testing, the number of iterations required for convergence
!    rarely exceeded 25.  Often 10 iterations were sufficient.
!
!    MAXITS = 90 will allow the inversion to work to its full capability
!    and is suggested unless the user is in a hurry.
!
!** RECORD 16: MAXITS, CNVRG, NFIX, MV1PRT, OUTPRT
!
!      MAXITS - Upper limit on the maximum number of iterations
!               It is unusual not to achieve convergence of some sort before
!               20 iterations.  The suggested value for MAXITS is 90.
!               Use a low value of MAXITS to limit of change the inversion.
!
!      CNVRG = 1 => iterations will proceed unless the error can no
!                   longer be reduced.
!
!            = 2 => iterations will not proceed any further if the
!                   RMS (root mean square) error is less than a user
!                   specified percent (PCTCNV in RECORD 16.1).
!
!      NFIX - the number of parameters that will be constrained.  If all
!             parameters are free to vary without restriction, NFIX = 0
!
!             If NFIX > 0, NFIX records describing the constraint must be
!             entered as RECORDS 17
!
!      MV1PRT refers to the output print level in AirBeo.mv1
!      OUTPRT refers to the output print level in AirBeo.OUT
!
!            =  0  No output DURING inversion.  The final model set AFTER inversion,
!                  but NOT the final model data, is written to output files.
!
!            =  1  as above plus final model data
!
!            =  2  as above plus intermediate model sets after each iteration
!
!            =  3 as above plus intermediate model data after each iteration
!
!     ------------------
!     only if CNVRG = 2
!**   RECORD 16.1: PCTCNV - terminate inversion if SQRT {SUMSQ / WSUM } <n PCTCNV
!     -------------------
!
!        SUMSQ - sum of the squares of the weighted symmetric error
!                 SERR(J) at each data point.
!
!        WSUM - sum of weights for all data points
!             = NDATA if all data points are equally weighted.
!
!                    W(J) * ABS ( M(J) - D(J) )
!        SERR(J) =   --------------------------
!                    (ABS (M(J)) + ABS (D(J)) ) / 2.
!
!        M(J), D(J) & W(J) are the model value, data value & weight
!        respectively at CHANNEL J.
!
!     ______________________________________________________________
!
!     only if NFIX > 0  (NFIX records)
!
!**   RECORD 17: CTYPE, LYR_INDX, KPAR, ELAS, LBND, UBND
!     ______________________________________________________________
!
!     CTYPE = 1 : parameter is fixed to initial model value
!                 only LYR_INDX & KPAR need be specified.
!                 ELAS(KPAR), LBND(KPAR), UBND(KPAR) are not read in
!
!           = 2 : frictional restraint
!                 only LYR_INDX, KPAR & ELAS(KPAR)need be specified.
!                 LBND(KPAR), UBND(KPAR) are not read in
!
!           = 3 : Buffered boundaries
!                 LYR_INDX, KPAR, ELAS(KPAR), LBND(KPAR), UBND(KPAR) must be specified.
!
!     LYR_INDX = constained layer number - increasing with depth
!     LYR_INDX = 1 or -1 refers to the top layer
!     LYR_INDX = NLYR or -NLYR refers to the basement.
!
!     The fact that LYR_INDX can be positive or negative is to allow
!     consistency with LeroiAir where plate indices are posiive and
!     layer indices are negative.  Positive indices are mor eintuitive,
!     especially for a layered earth program.
!
!     KPAR = 1 => layer resistivity
!            2 => layer thickness
!
!
!     ELAS(KPAR) - elasticity  of parameter KPAR  (0 < ELAS < 1)
!           If ELAS > 0.95, it is set to ELAS = 1
!           If ELAS < 0.05, it is set to ELAS = 0
!
!     LBND(KPAR) - lower bound of parameter KPAR
!     UBND(KPAR) - upper bound of parameter KPAR
!
!        -------------------------------
!    After each iteration, the inversion will compute a proposed step change
!    for each parameter.  Call this STEP(KPAR)
!
!     For CTYPE = 1, the allowed STEP change = 0
!
!     For CTYPE = 2 - restrained step
!        The maximum allowed parameter change will be ELAS * STEP(KPAR)
!        ELAS serves like a frictional restraint or rubber band preventing a
!        parameter from making the full change suggested by the inversion
!        process.  It is used when a parameter value is thought to be known
!        but allows more latitude than a hard fix.
!
!     For CTYPE = 3 - buffered bounds
!        Suppose the B1(KPAR) is the bound in the direction proposed by STEP(KPAR)
!        Define D1 as the distance between the current parameter value and B1
!        The maximum allowed step for the current iteration S1 = ELAS * D1
!        The actual parameter change will be the minimum of S1 & STEP(KPAR)
!
!
!=========================================
!                                        =
!     END OF ENTRIES FOR AirBeo.cfl      =
!         (Logical unit NR = 3)          =
!                                        =
!     BEGIN DESCRIPTON FOR AirBeo.inv    =
!         (Logical unit NRI = 13)        =
!                                        =
!=========================================
!
!      Any number of comment lines can be inserted at the beginning of the
!      .inv file by using either / or \ as the first character of the line.
!
!      AirBeo will start reading at the first line not containing / or \
!      as the first character.  Thereafter, any line beginning with  / or \
!      will cause an eror and an end to execution.
!
!
!             DATA DESCRIPTION & WEIGHTING  - read from AirBeo.inv
!             ============================
!
!    In what follows, the channel reference term PCHNL refers to each
!    component for each frequency or time-domain channel.
!
!    NPCHNL = the number of PCHNLs.
!    For example if there are 2 components of 10 CHANNEL time-domain data, NPCHNL = 20
!    For the 5 frequency DIGHEM system NPCHNL = 10 corresponding
!    to 5 in-phase and 5 quadrature data
!
!         Frequency-domain:
!         ----------------
!           Channels are odered from the lowest to the highest frequency.
!
!              PCHNL = 1 to NFRQ for the in-phase and
!                    = NFRQ+1 to 2*NFRQ for the quadrature.
!
!         Time-domain:
!         -----------
!            Regardless of the order in which the data is read in, for weighting
!            PCHNL = 1 to NCHNL refer to
!                               Vertical response if CMP = 13, 2 or 3
!                               In-Line response of CMP = 11
!                               Transverse response of CMP = 12
!
!            PCHNL = NCHNL+1 to 2*NCHNL refer to in-line response If CMP = 2
!            PCHNL = 2*NCHNL+1 to 3*NCHNL refer to transverse response of CMP = 3
!
!
!**  RECORD 1: NSTAT, SURVEY, BAROMTRC, KCMP, ORDER
!
!      NSTAT - number of stations to be read in from LeroiAir.inv
!
!        SURVEY = 2 variable altitude & course but constant Tx-Rx geometry
!        SURVEY = 3 variable altitude, course, transmitter pitch
!                          and receiver offset (time-domain only)
!
!      BAROMTRC = 0 => altitudes are ground clearance in metres.
!               = 1 => altitudes are barometric; ie, metres above sea level.
!
!      (time-domain)
!      KCMP = 1   : only X (in-line component) data will be read in
!           = 3   : only Z (vertical component) data will be read in
!           = 13  : X data followed by Z data
!           = 31  : Z data followed by X data
!           = 123 : X data followed by Y data followed by Z data
!           = 312 : Z data followed by X data followed by Y data
!
!      (time-domain)
!      ORDER = 1 the first data component for all times is followed
!                by the next data component for all times in the order
!                specified by KCMP.  (for example X1 X2 Y1 Y2 Z1 Z2)
!
!            = 2 all of the data components for each channel are followed by
!                all of the data components for the next channel in the order
!                specified by KCMP   (for example X1 Y1 Z1 X2 Y2 Z2)
!
!
!      (frequency-domain)
!      KCMP(1:NFRQ) : specify components for each frequency.  These must be in the
!                     order as the frequencies in LeroiAir.cfl
!
!                  1 : HCP - horizontal coplanar
!                  2 : VCA - vertical coaxial
!                  3 : VCP - vertical coplanar
!                  4 : VCB - vertical coplanar broadside
!
!     Dighem example:       KCMP(1:5) = 1 2 2 1 1  for 5 frequencies starting from the lowest.
!     GTK wingtip example:  KCMP(1:2) = 3 3        for 2 frequencies starting from the lowest.
!
!      (frequency-domain)
!      ORDER = 1122 : all inphase data will be followed by all quadrature data
!            = 1212 : data consists of paired inphase and quadrature for each frequency
!            = 2211 : all quadrature data will be followed by all inphase data
!            = 2121 : data consists of paired quadrature and inphase for each frequency
!
!**  RECORD 2: DATA_FLOOR
!
!      Any data value whose absolute magnitude is less than DATA_FLOOR will be
!      weighted to zero.
!
!      For TIME-DOMAIN only one value is required.
!
!      For FREQUENCY-DOMAIN, 2 * NFRQ values must be entered.
!
!        The first NFRQ values refer to inphase measurements for all frequencies
!        followed by NFRQ values for all quadrature measurements.
!        These must be read in order from the lowest to the highest frequency
!
!
!**  RECORD 3: N0STAT, N0CHNL, N0PTS
!
!      N0STAT - number of stations for which all the data will be weighted to zero
!      N0CHNL - number of PCHNLs for which all the data will be weighted to zero
!      N0PTS  - number of data points not covered by K0STAT & K0CHNL which will be
!               weighted to zero
!
!     ------------------
!     only if N0STAT > 0
!**   RECORD 4:  K0STAT(1:K0STAT) - indices of stations for which all data
!     ------------------               will be weighted to 0
!
!     ------------------
!     only if N0STAT < 0
!**   RECORD 4:  K0STAT(1:K0STAT) - indices of stations for which all data
!     ------------------               will NOT be weighted to 0
!
!     ------------------
!     only if N0CHNL > 0
!**   RECORD 5:  K0CHNL(1:N0CHNL) - indices of PCHNLs for which all data
!     ------------------               will be weighted to 0
!
!     ------------------
!     only if N0PTS > 0
!**   RECORD 6:  (J0CH(I),J0ST(I)), I = 1,N0PTS)
!     ------------------
!         PCHNL and station indices of individual points to be weighted to 0.
!         using the above PCHNL ordering convention
!
!
!             DATA ENTRY   (read from AirBeo.inv)
!             ==========
!
!      NSTAT records -
!
!    for SURVEY = 2 or -2
!**  RECORD 7.J: LINE(J), EAST(J), NORTH(J), ALT(J), DATA(1:NPCHNL, J)
!
!    for SURVEY = 3 or -3 (time domain)
!**  RECORD 7.J: LINE(J), EAST(J), NORTH(J), ALT(J), TXCLN(J), ZRX(J), XRX(J), YRX(J), DATA(1:NPCHNL, J)
!                                       ZRX(J), XRX(J),
!        DATA(I,J) = datum of PCHNL I at station J
!          Note that unnormalised B data is in pT and unnormalised dB/dt data is in nT/s
!          unless otherwise specified in RECORD 2.2
!
!        For frequency-domain inversion or if KPPM > 0, normalised data are expressed
!        in ppm unless otherwise specified using NPPF /= 3.
!
!
!        ZRX(J), XRX(J) = vertical & inline offsets at station J
!
!        TXCLN(J) = transmitter inclination
!
!==============================
!  END DATA ENTRY DESCRIPTIOM |
!==============================
!
!============================================================================



      USE input_data_for_airBeo
 
      IMPLICIT NONE
      
      INTEGER ider , knrm , mprnt
      REAL rmserr

      REAL , ALLOCATABLE , DIMENSION(:) :: norm
      REAL , ALLOCATABLE , DIMENSION(:,:,:) :: btd
      COMPLEX , ALLOCATABLE , DIMENSION(:,:,:) :: bfd
      REAL cmp_start , cmp_end , elapsed
      LOGICAL wrt_nw1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmlfd/ nfrq,freq,txcln,txa90,nstat,sz,zrx,xrx,         &
     &                 yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,       &
     &                 bfd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmltd/ js,step,ider,nsx,swx,swy,npuls,pulse,           &
     &                        ntypls,ntyrp,trp,nchnl,topn,tcls,txcln,   &
     &                        nstat,sz,zrx,xrx,yrx,nlyr,res,reps,rmu,   &
     &                        thk,calf,ctau,cfreq,gstrp,astrp,btd

      NAMELIST /mynmlinv/  js,nw,nw1,nlg,mv1prt,outprt,maxits,cnvrg,pctcnv, &
     &                 ndata,xdata,xmodl,xwts,npar,cxpar,elas,lbnd,ubnd,&
     &                 tdfd,cmp,knrm,norm,step,ider,nsx,swx,swy,ntyrp,  &
     &                 trp,npuls,pulse,ntypls,nchnl,topn,tcls,gstrp,    &
     &                 astrp,nfrq,freq,txcln,txa90,nrx,nrxst,xrx,yrx,   &
     &                 zrx,nstat,sx,sy,sz,rx,ry,rz,gnd_lvl,title,line,  &
     &                 nlyr,thk,res,rmu,reps,calf,ctau,cfreq,mpar
 
      
      CALL cpu_time(cmp_start)

      CALL read_system_and_survey
      CALL read_model
 
      IF ( invert ) THEN
         OPEN (nri,FILE='AirBeo.inv',STATUS='OLD')
         OPEN (nw1,FILE='AirBeo.mv1',STATUS='REPLACE')
         OPEN (md1,FILE='AirBeo.mdl',STATUS='REPLACE')
         CALL read_invrt_cntrl_and_data
      ELSE
         mprnt = 100
         js = 0
         CALL write_model(nw,mprnt,js,nlyr,thk,res,chrg,ctau,cfreq,rmu, &
     &                    reps)
         OPEN (nw1,FILE='AirBeo.mf1',STATUS='REPLACE')
      ENDIF
 
      CALL set_survey
      wrt_nw1 = .TRUE.
      IF ( tdfd==2 .AND. cmp>1 ) wrt_nw1 = .FALSE.
      IF ( wrt_nw1 ) CALL write_nw1_initial
 
      IF ( mxerr==0 ) THEN
         WRITE (*,'(/T3,A//T3,A//)')                                    &
     &                              'Control file passed initial tests.'&
     &                              , 'Computation begins.'
      ELSEIF ( mxerr==1 ) THEN
         WRITE (*,'(/T3,A//T3,A//)') 'Computation begins.' ,            &
     &                          'Look at warning messages in AirBeo.log'
      ELSEIF ( mxerr==2 ) THEN
         WRITE (*,'(/T3,A//T3,A//T3,A)')                                &
     &           'FATAL INPUT DATA ERRORS IN AirBeo.cfl' ,              &
     &          'Refer to messages in AirBeo.log' ,                     &
     &          'Execution will not occur until these are corrected.'
         STOP
      ENDIF
!============================================================================
      IF ( istop==1 ) STOP
!============================================================================
 
      IF ( tdfd==1 ) THEN
                       ! Time-Domain
 
! For time-domain, set up frequencies, interpolation times
! For time-domain, call SET_SOURCE to compute dI/dt at the transmitter using
! the DC coupling if waveform at the receiver has been specified.  Then
! (for time-domain) convert PRM_TD to the peak primary dB/dt in NT if
! impulse response is output or B in pT for step response.
! SWY will be in amps/s  * Tx area * NTRN
 
! IDER = 0 => that current derivative must be computed: ISW = 1, 11 or 31
!      = 1 => that current derivative has specified through voltage
!             calibration: ISW = 10 or 30
!      = 4 => ISW = 4  (pure rectangular pulse)
 
         ider = 0
         IF ( isw==10 .OR. isw==30 .OR. isw==130 ) ider = 1
         IF ( isw==4 ) ider = 4
         CALL set_trp
         knrm = 3
         ALLOCATE (norm(knrm))
         IF ( isw==4 ) THEN
            norm = 1.E6
         ELSE
            
            CALL dcprm_td(xrx0,yrx0,zrx0,txcln0,txarea,prm_td)
            CALL set_source(step,isw,bffac,waveform,nsx,swx,swy,prm_td)
            write(*,*) nw,bunit,bffac,kppm,punit, ppfac,prm_td
            IF ( invert ) CALL set_norm_td(nw,bunit,bffac,kppm,punit,   &
     &           ppfac,prm_td,norm)
             write(*,*) nw,bunit,bffac,kppm,punit,ppfac,prm_td,norm
         ENDIF
      ELSE
         knrm = nfrq
         ALLOCATE (norm(knrm))
         CALL dcprm_fd(nfrq,xrx,yrx,zrx,txcln,txa90,prm_fd,ppfac,norm)
      ENDIF
 
!============================================================================
      IF ( invert ) THEN
 
         DO js = 1 , nstat
            IF ( do1d==1 ) THEN
               res = res0
               thk = thk0
            ENDIF
            DO jt = 1 , ndata
               xdata(jt) = rdata(jt,js)
               xwts(jt) = rwts(jt,js)
            ENDDO
            IF ( maxval(xwts)<1 ) THEN
               WRITE (nw,99004) js
               WRITE (*,99004) js
               CYCLE
            ENDIF

!            write(*,*) 'THIS IS mynmlinv'
!            write(*,nml=mynmlinv)
            !write(*,*) size(title)
            !write(*,*) sizeof(title)

            CALL nlsq2(js,nw,nw1,nlg,mv1prt,outprt,maxits,cnvrg,pctcnv, &
     &                 ndata,xdata,xmodl,xwts,npar,cxpar,elas,lbnd,ubnd,&
     &                 tdfd,cmp,knrm,norm,step,ider,nsx,swx,swy,ntyrp,  &
     &                 trp,npuls,pulse,ntypls,nchnl,topn,tcls,gstrp,    &
     &                 astrp,nfrq,freq,txcln,txa90,nrx,nrxst,xrx,yrx,   &
     &                 zrx,nstat,sx,sy,sz,rx,ry,rz,gnd_lvl,title,line,  &
     &                 nlyr,thk,res,rmu,reps,calf,ctau,cfreq,mpar,      &
     &                 rmserr)
 
            res(1:nlyr) = mpar(1:nlyr)
            thk(1:nlyr-1) = mpar(nlyr+1:npar)
            IF(nlyr>1) THEN
               CALL cnvrt2_depth(nlyr,thk,depth)
               WRITE (md1,99001) js , syd(js) , sxd(js) , rmserr, res(1:nlyr) ,    &
     &                        depth(1:nlyr-1) , thk(1:nlyr-1)
            ELSE
               WRITE (md1,99001) js , syd(js) , sxd(js) , rmserr, res(1:nlyr)
            ENDIF
 
99001       FORMAT (I5,2F12.1,100G13.4)
 
         ENDDO
!===================================================================================
 
      ELSE    ! FORWARD MODEL OPTION
 
         CLOSE (nr)
 
         IF ( tdfd==1 ) THEN
                         ! Time-Domain.
            ALLOCATE (btd(nchnl,nstat,3))
            btd = 0.
 
!  Compute BTD, the layered earth response convolved with the excitation waveform
!  as dB/dt in nT/s if STEP = 0;  or as B in pT if STEP = 1
 
            DO js = 1 , nstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
            WRITE(*,*) 'THIS IS mynmltd before hsboss_td'
            WRITE(*,nml=mynmltd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
               CALL hsboss_td(js,step,ider,nsx,swx,swy,npuls,pulse,     &
     &                        ntypls,ntyrp,trp,nchnl,topn,tcls,txcln,   &
     &                        nstat,sz,zrx,xrx,yrx,nlyr,res,reps,rmu,   &
     &                        thk,calf,ctau,cfreq,gstrp,astrp,btd)
            ENDDO
 
!  Write out the results.
 
            CALL write_td(nw,nw1,title,nstat,line,sxd,syd,sz,txdeg,rxd, &
     &                    ryd,rz,xrx,yrx,zrx,nchnl,tms,prfl,qunit,bunit,&
     &                    ppfac,bffac,prm_td,cmp,kppm,btd)
 
         ELSE
         !  Construct the frequency-domain response.
 
            ALLOCATE (bfd(nfrq,nstat,3))
            bfd = (0.,0.)
 
            DO js = 1 , nstat
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!           WRITE(*,*) 'THIS IS mynmlfd before hsboss_fd '
!           WRITE(*,nml=mynmlfd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
               CALL hsboss_fd(js,nfrq,freq,txcln,txa90,nstat,sz,zrx,xrx,&
     &                        yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,&
     &                        bfd)
            ENDDO
            CALL write_fd(nw,nw1,title,nstat,line,txcln,txa90,sxd,syd,  &
     &                    sz,rxd,ryd,rz,config,nfrq,freq,prfl,qunit,    &
     &                    ppfac,bffac,prm_fd,cmp,bfd)
         ENDIF
      ENDIF
 
      CALL date_and_time(date,time,zone,qqdt)
      qqhms(1:2) = qqdt(5:6)
 
      CALL cpu_time(cmp_end)
      elapsed = cmp_end - cmp_start
 
      IF ( invert ) THEN
         WRITE (nw,99006) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,       &
     &                    qqdt(1) , elapsed
         WRITE (*,99006) qqhms(1:2) , qqdt(3) , month(qqdt(2)) , qqdt(1)&
     &                   , elapsed
         WRITE (nw1,99002) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,      &
     &                     qqdt(1) , elapsed
99002    FORMAT (T1,'/'/'/ AirBeo inversion completed at ',I2.2,':',    &
     &           I2.2,' on',I3.2,1X,A,I5/T1,'/ Computation time = ',    &
     &           F10.2,' seconds.')
      ELSE
         WRITE (nw,99005) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,       &
     &                    qqdt(1) , elapsed
         WRITE (*,99005) qqhms(1:2) , qqdt(3) , month(qqdt(2)) , qqdt(1)&
     &                   , elapsed
         WRITE (nw1,99003) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,      &
     &                     qqdt(1) , elapsed
99003    FORMAT (T1,'/'/'/ AirBeo forward model completed at ',I2.2,':',&
     &           I2.2,' on',I3.2,1X,A,I5/T1,'/ Computation time = ',    &
     &           F10.2,' seconds.')
      ENDIF
 
      CLOSE (nw)
      STOP
99004 FORMAT (//T14,'=============================='/T15,               &
     &        'No inversion for station',I4/T14,                        &
     &        '==============================')
99005 FORMAT (//T3,'AirBeo forward model completed at ',I2.2,':',I2.2,  &
     &        ' on',I3.2,1X,A,I5//T3,'Computation time = ',F10.2,       &
     &        ' seconds.'//)
99006 FORMAT (//T3,'AirBeo inversion completed at ',I2.2,':',I2.2,' on',&
     &        I3.2,1X,A,I5//T3,'Computation time = ',F10.2,             &
     &        ' seconds.'//)
      
      
      END PROGRAM airBeo
