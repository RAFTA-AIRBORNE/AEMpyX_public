!*==dcprm_td.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE dcprm_td(xrx0,yrx0,zrx0,txcln0,txarea,prm_td)
!---------------------------------------------------------
 
!***  Called by: READ_INVERSION_CONTROL
 
! For time-domain, PRM_TD is the 3 component Tx-Rx dc coupling factor
! per unit dipole moment, expressed in NANOTESLAS per unit amp
! Multiply it by dI/dt and get nanovolts per m^2 which is the
! same as nT/s.  Multiply it by current and get nT.
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  SIGN CONVENTION:
!  ----------------
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!                               INPUT
!                               -----
!      TXCLN0 - angle in degrees that the transmitting dipole makes with vertical
!              (climb = positive for VMD transmitter)
!      TXAREA - transmitter area in sq. metres
!      ZRX0, XRX0 & YRX0 are the initial vertical, in-line and transverse offsets
!                        of the receiver relative to transmitter
!                        below = + ;  behind = + ;  port = +
!
!                                 OUTPUT (time domain)
!                                 --------------------
!     PRM_TD = primary field coupling factor for B in nT (unit TX moment)
!     PRM_TD(1) = in-line B
!     PRM_TD(2) = transverse B
!     PRM_TD(3) = vertical B for station
 
      IMPLICIT NONE
      REAL , PARAMETER :: Pi = 3.141592654
      REAL sntx , cstx , xbd , ybd , zbd , rbrd , rsq , rsq1 , bfac ,   &
     &     fac , faczx , txarea , txcln0 , theta , xrx0 , yrx0 , zrx0 , &
     &     inline , vert , trans , prm_td(3)
 
!       write(*,*) ' Now in dcprm_td'
!  BFAC = 1.0E9 * MU / (4 * PI) for time domain  ! NANOTESLAS
 
      prm_td = 0.
 
!----------------------------------------------------------------------------
! In-loop time-domain HEM
 
      IF ( txarea>1. ) THEN
         rbrd = abs(xrx0) + abs(yrx0)
         IF ( rbrd<1. ) THEN
            zbd = sqrt(zrx0**2+txarea/Pi)
            prm_td(3) = 200./zbd**3  ! 1.0E9 * MU / (2 * PI) = 200. (nT)
            RETURN
         ENDIF
      ENDIF
!----------------------------------------------------------------------------
 
      bfac = 100.
      theta = txcln0*Pi/180.
      sntx = sin(theta)
      cstx = cos(theta)
 
      xbd = -xrx0
              !  XRX is defined as positive behind the TX.
      ybd = yrx0
      rbrd = sqrt(xbd**2+ybd**2)
      zbd = zrx0
      rsq = zbd**2 + rbrd**2
      rsq1 = sqrt(rsq)
      fac = bfac/rsq1**5
      faczx = 3.*fac*xbd*zbd
      vert = fac*cstx*(3.*zbd**2-rsq) + sntx*faczx
      inline = fac*sntx*(3.*xbd**2-rsq) + cstx*faczx
      trans = 3.*fac*ybd*((cstx*zbd)+(sntx*xbd))
      prm_td(1) = inline
      prm_td(2) = trans
      prm_td(3) = vert
 
      END SUBROUTINE dcprm_td
