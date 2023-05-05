!*==dcprm_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE dcprm_fd(nfrq,xrx,yrx,zrx,txcln,txa90,prm_fd,ppfac,    &
     &                    norm)
!--------------------------------------------------------------------
 
!***  Called by: MAIN
 
!  In frequency-domain, it computes the maximally coupled component of B at each
!  receiver location for each frequency assuming unit dipoles and current
!  transmitters are co-oriented.
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
!       NFRQ - number of frequencies
!        ZRX - vertical offset of RX relative to transmitter (below = + ).
!        XRX - in-line offset of RX relative to transmitter (behind = + ).
!        YRX - transverse offset of RX relative to transmitter (port = + ).
!
!       PPFAC = 100  => parts per hundred (percent) pct
!            or 1000 => parts per thousand (ppt)
!            or 1.e6 => parts per million (ppm)
!            or 1.e9 => parts per billion (ppb)
!
!
!                                 OUTPUT (frequency domain)
!                                 -------------------------
!     PRM_FD(1:NFRQ)  = primary B (nT) per unit dipole moment at each frequency.
!       NORM(1:NFRQ)  = PPM normalisation factor for fields expresed in nT
 
      IMPLICIT NONE
      INTEGER nfrq , jf
      REAL sntx , cstx , xbd , ybd , zbd , rbrd , rsq , rsq1 , bfac ,   &
     &     fac , faczx , inline , vert , ppfac
      REAL , DIMENSION(nfrq) :: txcln , xrx , yrx , zrx , prm_fd , norm
      LOGICAL coplanar , txa90
 
!  BFAC = 1.0E9 * MU / (4 * PI)  ! NANOTESLAS
 
      prm_fd = 0.
      bfac = 100.
 
      DO jf = 1 , nfrq
         sntx = sin(txcln(jf))
         cstx = cos(txcln(jf))
 
 
         xbd = -xrx(jf)
                   !  XRX is defined as positive behind the TX.
         ybd = yrx(jf)
         rbrd = sqrt(xbd**2+ybd**2)
         zbd = zrx(jf)
         rsq = zbd**2 + rbrd**2
         rsq1 = sqrt(rsq)
         coplanar = .FALSE.
         IF ( abs(sntx)<.01 ) coplanar = .TRUE.
         IF ( txa90 ) coplanar = .TRUE.
         IF ( coplanar ) THEN
            prm_fd(jf) = -bfac/rsq1**3
         ELSE
            fac = bfac/rsq1**5
            faczx = 3.*fac*xbd*zbd
            vert = fac*cstx*(3.*zbd**2-rsq) + sntx*faczx
            inline = fac*sntx*(3.*xbd**2-rsq) + cstx*faczx
            prm_fd(jf) = cstx*vert + sntx*inline
         ENDIF
         norm(jf) = ppfac/abs(prm_fd(jf))
      ENDDO
 
      END SUBROUTINE dcprm_fd
