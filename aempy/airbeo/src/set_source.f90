!*==set_source.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE set_source(step,isw,bffac,waveform,nsx,swx,swy,prm_td)
!------------------------------------------------------------------
 
! For time-domain, SET_SOURCE computes dI/dt at the transmitter using
! the DC coupling if waveform at the receiver has been specified.  Then
! (for time-domain) it converts PRM_TD to the peak primary dB/dt in nT/s
! if impulse response is required or B in nT for step response.
!
! SWY will be in amps / sec * Tx area * NTRN
!  Computes SWY to be TXMNT * dI(t)/dt & PKSX, the peak response.
!  The units of SWY are amps * m^2 / s
 
!*** Called by: MAIN
 
!             INPUT
!             -----
!
!       STEP = 0 for impulse response (dB/dt)
!            = 1 for step response (B)
!        ISW - waveform indicator
!      BUNIT & BFAC are set in SUBROUTINE READ_SYSTEM_SETUP
!      BUNIT can have values nT/s, nT, pT, pT/s, fT or fT/s
!      BFFAC is the conversion factor needed to achieve this from nT or nT/s
!            = 1, 1000, or 1E6
!   WAVEFORM - amps * TXMNT at the transmitter if ISW = 1
!            - vertical dB/dt if ISW = 30
!            - vertical B if ISW = 31
!            - horizontal dB/dt if ISW = 10
!            - horizontal B if ISW = 11
!        NSX - number of source points in waveform
!        SWX - time abscissae for input waveform in seconds
!
!   PRM_TD(I) = Ith component of B in nT per unit dipole moment
!               at the receiver
!
!        I = 1, 2 & 3 are the in-line, transverse & vertical components.
!
!             OUTPUT
!             ------
!
!        SWY - TRXMNT * dI(t)/dt  amps/s * m^2
!     PRM_TD - peak primary dB/dt in nT/s if STEP = 0 or
!              peak primary B (in nT if STEP = 1)
 
      IMPLICIT NONE
      INTEGER isw , step , nsx , jt
      REAL bffac , swx(nsx) , waveform(nsx) , swy(nsx,3) , pksx , delt ,&
     &     coupling , prm_td(3)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmlso/ step,isw,bffac,waveform,nsx,swx,swy,prm_td
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) ' start set_norm_source'  
        write(*,nml=mynmlso)

 
      IF ( swx(2)-swx(1)<0.5E-7 ) swx(2) = swx(1) + 0.5E-7
      IF ( swx(nsx)-swx(nsx-1)<0.5E-7 ) swx(nsx-1) = swx(nsx) - 0.5E-7
 
!  Remove the receiver coupling if the receiver voltage or magnetic field is
!  to be used to get dI/dt.  Ensure that the input waveform is converted to
!  nT or nT/s
 
      swy = 0.
      IF ( isw/=1 ) waveform = waveform/bffac
      swy(1:nsx,3) = waveform(1:nsx)
                                 !  Store original waveform for Geotem stripping.
 
! Remove receiver coupling factor if ISW isn't given in amps.
! Compensate 1.E9 for the fact that PRM_TD is in nT
 
      coupling = 1.
      IF ( isw==30 .OR. isw==31 .OR. isw==130 .OR. isw==131 ) THEN
         coupling = prm_td(3)
                        !  Waveform derived from vertical measurement
      ELSEIF ( isw==10 .OR. isw==11 ) THEN
         coupling = prm_td(1)
                        ! Waveform derived from in-line measurement
      ENDIF
      coupling = abs(coupling)
 
! The current is obtained by dividing the waveform by the coupling.  At this
! stage, the coupling, PRM_TD, is expressed in nT for unit dipole moment.
! The waveform is in nT/s or nT so coupling meeds to be multiplied
! by BFFAC to also be in UNITS thus get the current in amps and SWY in amps / s.
 
      IF ( isw/=1 ) waveform = waveform/coupling
                                               ! TXMNT * dI/dt or I
 
!  Compute the source waveform for the program from the input data.
!  This is dI/dt * the Tx-Rx moment.
 
      IF ( isw/=30 .AND. isw/=10 .AND. isw/=130 ) THEN
 
! Compute delta I in SWY(*,2) and dI/dt if it exists in SW(*,1).
! Compute the negative derivative so that dI/dt will be positive
! just before signal turn-off.  Store on right node (27.01.00)
 
         DO jt = 2 , nsx
            swy(jt,2) = waveform(jt-1) - waveform(jt)
            delt = swx(jt) - swx(jt-1)
            IF ( delt>1.0E-7 ) swy(jt,1) = swy(jt,2)/delt
         ENDDO
      ELSEIF ( step==0 ) THEN
         swy(1:nsx,1) = -waveform(1:nsx)
                                      ! Reverse negatibe derivative at signal end
      ELSE
         swy(1:nsx,1) = waveform(1:nsx)
      ENDIF
 
      pksx = 0
           ! Compute peak I or dI/dt for step normalisation.
 
      IF ( step==0 ) THEN
         pksx = maxval(abs(swy(1:nsx-1,1)))
                                           ! Compute peak dI/dt.
      ELSEIF ( step==1 ) THEN
         pksx = maxval(abs(waveform))
      ENDIF
 
      prm_td = pksx*prm_td
 
      END SUBROUTINE set_source
