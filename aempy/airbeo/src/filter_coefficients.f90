!*==filter_coefficients.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
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
 
 
      MODULE filter_coefficients
!  --------------------------
 
      IMPLICIT NONE
 
      INTEGER , PARAMETER :: Jnlo = -250 , Jnhi = 150 , Ndec_jn = 15
      INTEGER j9
      REAL shftjn , wj0(Jnlo:Jnhi) , wj1(Jnlo:Jnhi)
      SAVE 
 
!  Filter restored to original LeroiAir 7 February, 2000 (artificial shift removed)
 
!  J0 filter coefficients computed from the Niels Christensen program, FILCOA
!  for the following parameters:
!
!   ANY =  0      AMY =  0      NDEC = 15       NLO = -250        NHI =  150
!   IOPT = 1   ISHIFT = 0      OMEGA = .3 PI    EPS = 1.0E-12      SC = 3.257209
!      A = 0.162875             DEL0 = 0.14314998               ERROR =  1.4032E-08
 
      DATA shftjn/0.14314998/
      DATA (wj0(j9),j9=-250,-161)/2.86608135867E-18 ,                   &
     &      3.34160553102E-18 , 3.89602601168E-18 , 4.54243283439E-18 , &
     &      5.29608785801E-18 , 6.17478510356E-18 , 7.19927087644E-18 , &
     &      8.39373359283E-18 , 9.78637487555E-18 , 1.14100754027E-17 , &
     &      1.33031712306E-17 , 1.55103589191E-17 , 1.80837508313E-17 , &
     &      2.10841055215E-17 , 2.45822622636E-17 , 2.86608135867E-17 , &
     &      3.34160553102E-17 , 3.89602601168E-17 , 4.54243283439E-17 , &
     &      5.29608785801E-17 , 6.17478510356E-17 , 7.19927087644E-17 , &
     &      8.39373359283E-17 , 9.78637487555E-17 , 1.14100754027E-16 , &
     &      1.33031712306E-16 , 1.55103589191E-16 , 1.80837508313E-16 , &
     &      2.10841055215E-16 , 2.45822622636E-16 , 2.86608135867E-16 , &
     &      3.34160553102E-16 , 3.89602601168E-16 , 4.54243283439E-16 , &
     &      5.29608785801E-16 , 6.17478510356E-16 , 7.19927087644E-16 , &
     &      8.39373359283E-16 , 9.78637487555E-16 , 1.14100754027E-15 , &
     &      1.33031712306E-15 , 1.55103589191E-15 , 1.80837508313E-15 , &
     &      2.10841055215E-15 , 2.45822622636E-15 , 2.86608135867E-15 , &
     &      3.34160553102E-15 , 3.89602601168E-15 , 4.54243283439E-15 , &
     &      5.29608785801E-15 , 6.17478510356E-15 , 7.19927087644E-15 , &
     &      8.39373359283E-15 , 9.78637487555E-15 , 1.14100754027E-14 , &
     &      1.33031712306E-14 , 1.55103589191E-14 , 1.80837508313E-14 , &
     &      2.10841055215E-14 , 2.45822622636E-14 , 2.86608135867E-14 , &
     &      3.34160553102E-14 , 3.89602601168E-14 , 4.54243283439E-14 , &
     &      5.29608785801E-14 , 6.17478510356E-14 , 7.19927087644E-14 , &
     &      8.39373359283E-14 , 9.78637487555E-14 , 1.14100754027E-13 , &
     &      1.33031712306E-13 , 1.55103589191E-13 , 1.80837508313E-13 , &
     &      2.10841055215E-13 , 2.45822622636E-13 , 2.86608135867E-13 , &
     &      3.34160553102E-13 , 3.89602601168E-13 , 4.54243283439E-13 , &
     &      5.29608785801E-13 , 6.17478510356E-13 , 7.19927087644E-13 , &
     &      8.39373359283E-13 , 9.78637487555E-13 , 1.14100754027E-12 , &
     &      1.33031712306E-12 , 1.55103589191E-12 , 1.80837508313E-12 , &
     &      2.10841055215E-12 , 2.45822622636E-12/
      DATA (wj0(j9),j9=-160,-71)/2.86608135867E-12 , 3.34160553102E-12 ,&
     &      3.89602601168E-12 , 4.54243283439E-12 , 5.29608785801E-12 , &
     &      6.17478510356E-12 , 7.19927087644E-12 , 8.39373359283E-12 , &
     &      9.78637487555E-12 , 1.14100754027E-11 , 1.33031712306E-11 , &
     &      1.55103589191E-11 , 1.80837508313E-11 , 2.10841055215E-11 , &
     &      2.45822622636E-11 , 2.86608135867E-11 , 3.34160553102E-11 , &
     &      3.89602601168E-11 , 4.54243283439E-11 , 5.29608785801E-11 , &
     &      6.17478510356E-11 , 7.19927087644E-11 , 8.39373359283E-11 , &
     &      9.78637487555E-11 , 1.14100754027E-10 , 1.33031712306E-10 , &
     &      1.55103589191E-10 , 1.80837508313E-10 , 2.10841055215E-10 , &
     &      2.45822622636E-10 , 2.86608135867E-10 , 3.34160553102E-10 , &
     &      3.89602601168E-10 , 4.54243283439E-10 , 5.29608785801E-10 , &
     &      6.17478510356E-10 , 7.19927087644E-10 , 8.39373359283E-10 , &
     &      9.78637487555E-10 , 1.14100754027E-09 , 1.33031712306E-09 , &
     &      1.55103589191E-09 , 1.80837508313E-09 , 2.10841055215E-09 , &
     &      2.45822622636E-09 , 2.86608135867E-09 , 3.34160553102E-09 , &
     &      3.89602601168E-09 , 4.54243283439E-09 , 5.29608785801E-09 , &
     &      6.17478510356E-09 , 7.19927087644E-09 , 8.39373359283E-09 , &
     &      9.78637487555E-09 , 1.14100754027E-08 , 1.33031712306E-08 , &
     &      1.55103589191E-08 , 1.80837508313E-08 , 2.10841055215E-08 , &
     &      2.45822622636E-08 , 2.86608135867E-08 , 3.34160553102E-08 , &
     &      3.89602601168E-08 , 4.54243283439E-08 , 5.29608785801E-08 , &
     &      6.17478510356E-08 , 7.19927087644E-08 , 8.39373359283E-08 , &
     &      9.78637487555E-08 , 1.14100754027E-07 , 1.33031712306E-07 , &
     &      1.55103589191E-07 , 1.80837508313E-07 , 2.10841055215E-07 , &
     &      2.45822622635E-07 , 2.86608135866E-07 , 3.34160553102E-07 , &
     &      3.89602601167E-07 , 4.54243283438E-07 , 5.29608785799E-07 , &
     &      6.17478510354E-07 , 7.19927087640E-07 , 8.39373359277E-07 , &
     &      9.78637487545E-07 , 1.14100754026E-06 , 1.33031712304E-06 , &
     &      1.55103589187E-06 , 1.80837508307E-06 , 2.10841055205E-06 , &
     &      2.45822622620E-06/
      DATA (wj0(j9),j9=-70,19)/2.86608135842E-06 , 3.34160553063E-06 ,  &
     &      3.89602601105E-06 , 4.54243283340E-06 , 5.29608785643E-06 , &
     &      6.17478510107E-06 , 7.19927087248E-06 , 8.39373358656E-06 , &
     &      9.78637486561E-06 , 1.14100753870E-05 , 1.33031712056E-05 , &
     &      1.55103588795E-05 , 1.80837507685E-05 , 2.10841054221E-05 , &
     &      2.45822621060E-05 , 2.86608133369E-05 , 3.34160549143E-05 , &
     &      3.89602594894E-05 , 4.54243273495E-05 , 5.29608770041E-05 , &
     &      6.17478485378E-05 , 7.19927048056E-05 , 8.39373296541E-05 , &
     &      9.78637388116E-05 , 1.14100738267E-04 , 1.33031687328E-04 , &
     &      1.55103549604E-04 , 1.80837445571E-04 , 2.10840955776E-04 , &
     &      2.45822465035E-04 , 2.86607886087E-04 , 3.34160157229E-04 , &
     &      3.89601973751E-04 , 4.54242289050E-04 , 5.29607209800E-04 , &
     &      6.17476012564E-04 , 7.19923128912E-04 , 8.39367085119E-04 , &
     &      9.78627543681E-04 , 1.14099178031E-03 , 1.33029214523E-03 , &
     &      1.55099630479E-03 , 1.80831234191E-03 , 2.10831111434E-03 , &
     &      2.45806862870E-03 , 2.86583158466E-03 , 3.34120966900E-03 , &
     &      3.89539861933E-03 , 4.54143849891E-03 , 5.29451197347E-03 , &
     &      6.17228756167E-03 , 7.19531268313E-03 , 8.38746058912E-03 , &
     &      9.77643350230E-03 , 1.13943208262E-02 , 1.32782050079E-02 , &
     &      1.54707967971E-02 , 1.80210634703E-02 , 2.09847837166E-02 , &
     &      2.44249145050E-02 , 2.84115778193E-02 , 3.30213524808E-02 , &
     &      3.83353639832E-02 , 4.44353673090E-02 , 5.13965627145E-02 , &
     &      5.92752031985E-02 , 6.80880607240E-02 , 7.77794366644E-02 , &
     &      8.81696149649E-02 , 9.88766639298E-02 , 1.09202052802E-01 , &
     &      1.17971700371E-01 , 1.23332521049E-01 , 1.22530035854E-01 , &
     &      1.11753240889E-01 , 8.62569960973E-02 , 4.11899187108E-02 , &
     &      -2.61456504772E-02 , -1.11691705121E-01 ,                   &
     &      -1.97411432453E-01 , -2.44254055664E-01 ,                   &
     &      -1.95918893763E-01 , -1.49300191739E-02 ,                   &
     &      2.33634698676E-01 , 3.13582629541E-01 , -4.47760615930E-03 ,&
     &      -3.86535797015E-01 , -3.87589109967E-03 ,                   &
     &      4.18653972543E-01 , -4.16298788795E-01/
      DATA (wj0(j9),j9=20,109)/2.34448877498E-01 , -9.52158343728E-02 , &
     &      3.09020778713E-02 , -8.49535839509E-03 , 2.06835506815E-03 ,&
     &      -4.67185821059E-04 , 1.02086153218E-04 ,                    &
     &      -2.20830053233E-05 , 4.76413760468E-06 ,                    &
     &      -1.02705545675E-06 , 2.21421979164E-07 ,                    &
     &      -4.77750910705E-08 , 1.03340738634E-08 ,                    &
     &      -2.25102276694E-09 , 4.99715357680E-10 ,                    &
     &      -1.16500471179E-10 , 3.03986897639E-11 ,                    &
     &      -9.72611811870E-12 , 3.99994042396E-12 ,                    &
     &      -2.00348565820E-12 , 1.11608417099E-12 ,                    &
     &      -6.50767639555E-13 , 3.86180817012E-13 ,                    &
     &      -2.30659587418E-13 , 1.38093695980E-13 ,                    &
     &      -8.27455585993E-14 , 4.95961642994E-14 ,                    &
     &      -2.97302965597E-14 , 1.78224472343E-14 ,                    &
     &      -1.06841897105E-14 , 6.40498685290E-15 ,                    &
     &      -3.83968417568E-15 , 2.30182896520E-15 ,                    &
     &      -1.37991039489E-15 , 8.27234374391E-16 ,                    &
     &      -4.95913890248E-16 , 2.97292643817E-16 ,                    &
     &      -1.78222228351E-16 , 1.06841401468E-16 ,                    &
     &      -6.40497544674E-17 , 3.83968128138E-17 ,                    &
     &      -2.30182807939E-17 , 1.37991004842E-17 ,                    &
     &      -8.27234560136E-18 , 4.95913797287E-18 ,                    &
     &      -2.97292590016E-18 , 1.78222272891E-18 ,                    &
     &      -1.06841382487E-18 , 6.40497431324E-19 ,                    &
     &      -3.83968224515E-19 , 2.30182767120E-19 ,                    &
     &      -1.37990980321E-19 , 8.27234414081E-20 ,                    &
     &      -4.95914134387E-20 , 2.97292537295E-20 ,                    &
     &      -1.78222241286E-20 , 1.06841455108E-20 ,                    &
     &      -6.40497317742E-21 , 3.83968156424E-21 ,                    &
     &      -2.30182923671E-21 , 1.37990955793E-21 ,                    &
     &      -8.27234267383E-22 , 4.95914046240E-22 ,                    &
     &      -2.97292739490E-22 , 1.78222209690E-22 ,                    &
     &      -1.06841436161E-22 , 6.40497753124E-23 ,                    &
     &      -3.83968088314E-23 , 2.30182784256E-23 ,                    &
     &      -1.37991049701E-23 , 8.27234475022E-24 ,                    &
     &      -4.95913958682E-24 , 2.97292559305E-24 ,                    &
     &      -1.78222330828E-24 , 1.06841371450E-24 ,                    &
     &      -6.40497639510E-25 , 3.83968184851E-25 ,                    &
     &      -2.30182842033E-25 , 1.37990966066E-25 ,                    &
     &      -8.27234682962E-26 , 4.95914083158E-26 ,                    &
     &      -2.97292634049E-26 , 1.78222222810E-26 ,                    &
     &      -1.06841489841E-26 , 6.40497251344E-27 ,                    &
     &      -3.83968281228E-27 , 2.30182702533E-27 ,                    &
     &      -1.37991000702E-27 , 8.27234181627E-28 , -4.95914207635E-28/
      DATA wj0(110:150)/2.97292963477E-28 , -1.78222420371E-28 ,        &
     &     1.06841425086E-28 , -6.40497412376E-29 , 3.83968377606E-29 , &
     &     -2.30182957681E-29 , 1.37991153609E-29 , -8.27235098582E-30 ,&
     &     4.95914332316E-30 , -2.97292528486E-30 , 1.78222312353E-30 , &
     &     -1.06841451903E-30 , 6.40498122076E-31 , -3.83968474142E-31 ,&
     &     2.30183015458E-31 , -1.37991188353E-31 , 8.27234597206E-32 , &
     &     -4.95914031749E-32 , 2.97292858145E-32 , -1.78222357152E-32 ,&
     &     1.06841478804E-32 , -6.40498282844E-33 , 3.83968570659E-33 , &
     &     -2.30182876031E-33 , 1.37991104718E-33 , -8.27234805187E-34 ,&
     &     4.95914156225E-34 , -2.97292932767E-34 , 1.78222401887E-34 , &
     &     -1.06841414093E-34 , 6.40497895409E-35 , -3.83968338099E-35 ,&
     &     2.30182933903E-35 , -1.37991139355E-35 , 8.27235013127E-36 , &
     &     -4.95914281087E-36 , 2.97292752582E-36 , -1.78222294016E-36 ,&
     &     1.06841440910E-36 , -6.40498056176E-37 , 3.83968434477E-37/
 
!  J1 filter coefficients computed from the Niels Christensen program, FILCOA
!  for the following parameters:
!
!   ANY =  1      AMY =  0      NDEC = 15       NLO = -250        NHI =  150
!   IOPT = 1   ISHIFT = 0      OMEGA = .3 PI    EPS = 1.0E-12      SC = 3.257209
!      A = 0.162875             DEL0 = 0.14314998               ERROR =  1.4032E-08
 
      DATA (wj1(j9),j9=-250,-161)/2.67560875879E-35 ,                   &
     &      3.63710586576E-35 , 4.94412310292E-35 , 6.72082533724E-35 , &
     &      9.13599687416E-35 , 1.24190757379E-34 , 1.68819499732E-34 , &
     &      2.29485865865E-34 , 3.11953078380E-34 , 4.24055410750E-34 , &
     &      5.76442432690E-34 , 7.83590704850E-34 , 1.06517903247E-33 , &
     &      1.44795792522E-33 , 1.96829085937E-33 , 2.67560875879E-33 , &
     &      3.63710586576E-33 , 4.94412310292E-33 , 6.72082533724E-33 , &
     &      9.13599687416E-33 , 1.24190757379E-32 , 1.68819499732E-32 , &
     &      2.29485865865E-32 , 3.11953078380E-32 , 4.24055410750E-32 , &
     &      5.76442432690E-32 , 7.83590704850E-32 , 1.06517903247E-31 , &
     &      1.44795792522E-31 , 1.96829085937E-31 , 2.67560875879E-31 , &
     &      3.63710586576E-31 , 4.94412310292E-31 , 6.72082533724E-31 , &
     &      9.13599687416E-31 , 1.24190757379E-30 , 1.68819499732E-30 , &
     &      2.29485865865E-30 , 3.11953078380E-30 , 4.24055410750E-30 , &
     &      5.76442432690E-30 , 7.83590704850E-30 , 1.06517903247E-29 , &
     &      1.44795792522E-29 , 1.96829085937E-29 , 2.67560875879E-29 , &
     &      3.63710586576E-29 , 4.94412310292E-29 , 6.72082533724E-29 , &
     &      9.13599687416E-29 , 1.24190757379E-28 , 1.68819499732E-28 , &
     &      2.29485865865E-28 , 3.11953078380E-28 , 4.24055410750E-28 , &
     &      5.76442432690E-28 , 7.83590704850E-28 , 1.06517903247E-27 , &
     &      1.44795792522E-27 , 1.96829085937E-27 , 2.67560875879E-27 , &
     &      3.63710586576E-27 , 4.94412310292E-27 , 6.72082533724E-27 , &
     &      9.13599687416E-27 , 1.24190757379E-26 , 1.68819499732E-26 , &
     &      2.29485865865E-26 , 3.11953078380E-26 , 4.24055410750E-26 , &
     &      5.76442432690E-26 , 7.83590704850E-26 , 1.06517903247E-25 , &
     &      1.44795792522E-25 , 1.96829085937E-25 , 2.67560875879E-25 , &
     &      3.63710586576E-25 , 4.94412310292E-25 , 6.72082533724E-25 , &
     &      9.13599687416E-25 , 1.24190757379E-24 , 1.68819499732E-24 , &
     &      2.29485865865E-24 , 3.11953078380E-24 , 4.24055410750E-24 , &
     &      5.76442432690E-24 , 7.83590704850E-24 , 1.06517903247E-23 , &
     &      1.44795792522E-23 , 1.96829085937E-23/
      DATA (wj1(j9),j9=-160,-71)/2.67560875879E-23 , 3.63710586576E-23 ,&
     &      4.94412310292E-23 , 6.72082533724E-23 , 9.13599687416E-23 , &
     &      1.24190757379E-22 , 1.68819499732E-22 , 2.29485865865E-22 , &
     &      3.11953078380E-22 , 4.24055410750E-22 , 5.76442432690E-22 , &
     &      7.83590704850E-22 , 1.06517903247E-21 , 1.44795792522E-21 , &
     &      1.96829085937E-21 , 2.67560875879E-21 , 3.63710586576E-21 , &
     &      4.94412310292E-21 , 6.72082533724E-21 , 9.13599687416E-21 , &
     &      1.24190757379E-20 , 1.68819499732E-20 , 2.29485865865E-20 , &
     &      3.11953078380E-20 , 4.24055410750E-20 , 5.76442432690E-20 , &
     &      7.83590704850E-20 , 1.06517903247E-19 , 1.44795792522E-19 , &
     &      1.96829085937E-19 , 2.67560875879E-19 , 3.63710586576E-19 , &
     &      4.94412310292E-19 , 6.72082533724E-19 , 9.13599687416E-19 , &
     &      1.24190757379E-18 , 1.68819499732E-18 , 2.29485865865E-18 , &
     &      3.11953078380E-18 , 4.24055410750E-18 , 5.76442432690E-18 , &
     &      7.83590704850E-18 , 1.06517903247E-17 , 1.44795792522E-17 , &
     &      1.96829085937E-17 , 2.67560875879E-17 , 3.63710586576E-17 , &
     &      4.94412310292E-17 , 6.72082533724E-17 , 9.13599687416E-17 , &
     &      1.24190757379E-16 , 1.68819499732E-16 , 2.29485865865E-16 , &
     &      3.11953078380E-16 , 4.24055410750E-16 , 5.76442432690E-16 , &
     &      7.83590704850E-16 , 1.06517903247E-15 , 1.44795792522E-15 , &
     &      1.96829085937E-15 , 2.67560875879E-15 , 3.63710586576E-15 , &
     &      4.94412310292E-15 , 6.72082533724E-15 , 9.13599687416E-15 , &
     &      1.24190757379E-14 , 1.68819499732E-14 , 2.29485865865E-14 , &
     &      3.11953078380E-14 , 4.24055410750E-14 , 5.76442432690E-14 , &
     &      7.83590704849E-14 , 1.06517903247E-13 , 1.44795792522E-13 , &
     &      1.96829085938E-13 , 2.67560875878E-13 , 3.63710586577E-13 , &
     &      4.94412310288E-13 , 6.72082533728E-13 , 9.13599687406E-13 , &
     &      1.24190757380E-12 , 1.68819499729E-12 , 2.29485865868E-12 , &
     &      3.11953078372E-12 , 4.24055410758E-12 , 5.76442432666E-12 , &
     &      7.83590704871E-12 , 1.06517903240E-11 , 1.44795792527E-11 , &
     &      1.96829085917E-11/
      DATA (wj1(j9),j9=-70,19)/2.67560875891E-11 , 3.63710586515E-11 ,  &
     &      4.94412310317E-11 , 6.72082533541E-11 , 9.13599687462E-11 , &
     &      1.24190757324E-10 , 1.68819499736E-10 , 2.29485865695E-10 , &
     &      3.11953078363E-10 , 4.24055410221E-10 , 5.76442432542E-10 , &
     &      7.83590703194E-10 , 1.06517903172E-09 , 1.44795791998E-09 , &
     &      1.96829085611E-09 , 2.67560874206E-09 , 3.63710585268E-09 , &
     &      4.94412304898E-09 , 6.72082528725E-09 , 9.13599669890E-09 , &
     &      1.24190755523E-08 , 1.68819493996E-08 , 2.29485859113E-08 , &
     &      3.11953059487E-08 , 4.24055386543E-08 , 5.76442370102E-08 , &
     &      7.83590618983E-08 , 1.06517882412E-07 , 1.44795762309E-07 , &
     &      1.96829016283E-07 , 2.67560770231E-07 , 3.63710352883E-07 , &
     &      4.94411942636E-07 , 6.72081747305E-07 , 9.13598412795E-07 , &
     &      1.24190492063E-06 , 1.68819059152E-06 , 2.29484968860E-06 , &
     &      3.11951559104E-06 , 4.24052372735E-06 , 5.76437203602E-06 , &
     &      7.83580400571E-06 , 1.06516106220E-05 , 1.44792293329E-05 , &
     &      1.96822917833E-05 , 2.67548981332E-05 , 3.63689436167E-05 , &
     &      4.94371845248E-05 , 6.72010067340E-05 , 9.13461935181E-05 , &
     &      1.24165945005E-04 , 1.68772580859E-04 , 2.29400955289E-04 , &
     &      3.11793204874E-04 , 4.23764974965E-04 , 5.75897507579E-04 , &
     &      7.82597702990E-04 , 1.06332133421E-03 , 1.44456435715E-03 , &
     &      1.96195766368E-03 , 2.66401748131E-03 , 3.61551958902E-03 , &
     &      4.90456094796E-03 , 6.64729428357E-03 , 9.00112880743E-03 , &
     &      1.21689223295E-02 , 1.64231258930E-02 , 2.20996958736E-02 , &
     &      2.96400942278E-02 , 3.95385050500E-02 , 5.24078149405E-02 , &
     &      6.87615215337E-02 , 8.91013723344E-02 , 1.13192375541E-01 , &
     &      1.40192739735E-01 , 1.66618485339E-01 , 1.87030308669E-01 , &
     &      1.89612379729E-01 , 1.61380285157E-01 , 8.29859362099E-02 , &
     &      -4.46335736689E-02 , -2.01737898138E-01 ,                   &
     &      -2.84006740802E-01 , -1.90854624427E-01 ,                   &
     &      1.45861570853E-01 , 3.42338340245E-01 , 5.72930699760E-02 , &
     &      -4.71068534718E-01 , 2.63969067746E-01 , 8.25956507901E-02/
      DATA (wj1(j9),j9=20,109)/ - 2.22236420794E-01 ,                   &
     &      2.04428998525E-01 , -1.44401888321E-01 , 9.24618900674E-02 ,&
     &      -5.69896615248E-02 , 3.45697730305E-02 ,                    &
     &      -2.08227940873E-02 , 1.25054653306E-02 ,                    &
     &      -7.50178808640E-03 , 4.49828025678E-03 ,                    &
     &      -2.69688071237E-03 , 1.61678766116E-03 ,                    &
     &      -9.69249547051E-04 , 5.81052166908E-04 ,                    &
     &      -3.48332124427E-04 , 2.08819730575E-04 ,                    &
     &      -1.25184162926E-04 , 7.50459390809E-05 ,                    &
     &      -4.49888596104E-05 , 2.69701130091E-05 ,                    &
     &      -1.61681580285E-05 , 9.69255610555E-06 ,                    &
     &      -5.81053473294E-06 , 3.48332405883E-06 ,                    &
     &      -2.08819791213E-06 , 1.25184175990E-06 ,                    &
     &      -7.50459418954E-07 , 4.49888602168E-07 ,                    &
     &      -2.69701131398E-07 , 1.61681580566E-07 ,                    &
     &      -9.69255611161E-08 , 5.81053473425E-08 ,                    &
     &      -3.48332405911E-08 , 2.08819791219E-08 ,                    &
     &      -1.25184175991E-08 , 7.50459418957E-09 ,                    &
     &      -4.49888602168E-09 , 2.69701131398E-09 ,                    &
     &      -1.61681580566E-09 , 9.69255611161E-10 ,                    &
     &      -5.81053473425E-10 , 3.48332405911E-10 ,                    &
     &      -2.08819791219E-10 , 1.25184175991E-10 ,                    &
     &      -7.50459418957E-11 , 4.49888602168E-11 ,                    &
     &      -2.69701131398E-11 , 1.61681580566E-11 ,                    &
     &      -9.69255611161E-12 , 5.81053473425E-12 ,                    &
     &      -3.48332405911E-12 , 2.08819791219E-12 ,                    &
     &      -1.25184175991E-12 , 7.50459418957E-13 ,                    &
     &      -4.49888602168E-13 , 2.69701131398E-13 ,                    &
     &      -1.61681580566E-13 , 9.69255611161E-14 ,                    &
     &      -5.81053473425E-14 , 3.48332405911E-14 ,                    &
     &      -2.08819791219E-14 , 1.25184175991E-14 ,                    &
     &      -7.50459418957E-15 , 4.49888602168E-15 ,                    &
     &      -2.69701131398E-15 , 1.61681580566E-15 ,                    &
     &      -9.69255611161E-16 , 5.81053473425E-16 ,                    &
     &      -3.48332405911E-16 , 2.08819791219E-16 ,                    &
     &      -1.25184175991E-16 , 7.50459418957E-17 ,                    &
     &      -4.49888602168E-17 , 2.69701131398E-17 ,                    &
     &      -1.61681580566E-17 , 9.69255611161E-18 ,                    &
     &      -5.81053473425E-18 , 3.48332405911E-18 ,                    &
     &      -2.08819791219E-18 , 1.25184175991E-18 ,                    &
     &      -7.50459418957E-19 , 4.49888602168E-19 ,                    &
     &      -2.69701131398E-19 , 1.61681580566E-19 ,                    &
     &      -9.69255611161E-20 , 5.81053473425E-20 ,                    &
     &      -3.48332405911E-20 , 2.08819791219E-20 ,                    &
     &      -1.25184175991E-20 , 7.50459418957E-21/
      DATA wj1(110:150)/ - 4.49888602168E-21 , 2.69701131398E-21 ,      &
     &     -1.61681580566E-21 , 9.69255611161E-22 , -5.81053473425E-22 ,&
     &     3.48332405911E-22 , -2.08819791219E-22 , 1.25184175991E-22 , &
     &     -7.50459418957E-23 , 4.49888602168E-23 , -2.69701131398E-23 ,&
     &     1.61681580566E-23 , -9.69255611161E-24 , 5.81053473425E-24 , &
     &     -3.48332405911E-24 , 2.08819791219E-24 , -1.25184175991E-24 ,&
     &     7.50459418957E-25 , -4.49888602168E-25 , 2.69701131398E-25 , &
     &     -1.61681580566E-25 , 9.69255611161E-26 , -5.81053473425E-26 ,&
     &     3.48332405911E-26 , -2.08819791219E-26 , 1.25184175991E-26 , &
     &     -7.50459418957E-27 , 4.49888602168E-27 , -2.69701131398E-27 ,&
     &     1.61681580566E-27 , -9.69255611161E-28 , 5.81053473425E-28 , &
     &     -3.48332405911E-28 , 2.08819791219E-28 , -1.25184175991E-28 ,&
     &     7.50459418957E-29 , -4.49888602168E-29 , 2.69701131398E-29 , &
     &     -1.61681580566E-29 , 9.69255611161E-30 , -5.81053473425E-30/
 
      END MODULE filter_coefficients
