SUBROUTINE aemfwd1d_aem05( nlyr, mvec, alt,  calc_data)
!f2py intent(in)      nlyr,alt,mvec
!f2py depend(nlyr)    mvec
!f2py intent (out)    calc_data
!f2py threadsafe
!----------------------------------------------------------------------
!  Uses AEM1D_FD to compute  the frequency-domain layered earth
!  H field for a dipole of unit moment and current.
!.
!                             INPUT
!                             -----
!        JS = station reference
!      FREQ - array of NFRQ frequencies
!     TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!     TXA90 - true for vertical co-planar briadside array
!     NSTAT - number of stations in survey line.
!       ALT - array of transmitter altitudes
!       ZRX - vertical offset of each receiver from transmitter  (below = +)
!       XRX - in-line horizontal offset of RX J;  (behind = +)
!       YRX - transverse horizontal offset of RX J (port = +)
!      NLYR - number of layers
!      RES - layer resistivities (nlyr)
!      REPS - array of relative dislectric constants (nlyr)
!      RMUX - mu(i) / mu(0)   (nlyr)
!      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters. (nlyr)
!      THK - array of layer thicknesses (nlyr-1)
!
!--------------       ----------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   1 Dec 2020   VR
!------------------------------------------------------------------------
!

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nlyr
      REAL(KIND=8), INTENT(in)  :: mvec(7*nlyr)
      REAL(KIND=8), INTENT(out) :: calc_data(8)

!       REAL(KIND=8), ALLOCATABLE :: mcurrent(:)
!       REAL(KIND=8):: mcurrent(1000)
      INTEGER jf
      REAL(KIND=8), DIMENSION(nlyr) ::  thk,res,reps,rmu,calf,ctau,cfreq
      COMPLEX(KIND=8) bfd(4,3),xbfd

      INTEGER, PARAMETER :: nfrq = 4
      REAL(KIND=8)  :: ppfac=1.d6, alt
      REAL(KIND=8), DIMENSION(nfrq) :: zrx(1:nfrq) = 0.D0, xrx(1:nfrq) = 0.D0, yrx(1:nfrq) = 21.36d0
      REAL(KIND=8), DIMENSION(nfrq) :: freq (1:nfrq) = (/ 912.d0, 3005.d0, 11962.d0,24510.d0/)
      REAL(KIND=8), DIMENSION(nfrq) :: txcln(1:nfrq) = 1.57079637    ! 90. degrees
      REAL(KIND=8), DIMENSION(nfrq) :: cstx, sntx, norm, prm_fd
      LOGICAL :: txa90=.false.
      LOGICAL :: debug=.false.



! Set system-specific parameters for 4 freq GTK system

      cstx(1:nfrq) = dcos(txcln(1:nfrq))
      sntx(1:nfrq) = dsin(txcln(1:nfrq))

      !write(*,'(/A, i7)') ' transform= ',transform
      CALL setup_aem05(nfrq,xrx,yrx,zrx,txcln,txa90,prm_fd,ppfac,norm)
!       WRITE(*,'(  A/,7(/12G12.4))') 'trvec  ', mvec(0*nlyr+1:7*nlyr)
      CALL unpack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mvec(0*nlyr+1:7*nlyr))
      !WRITE(*,'( A/,7(/3G12.4))') 'unvec  ', mcurrent
      IF (debug) THEN
         WRITE(*,*) 'FWD-AEM05'
         WRITE(*,'(  A,3I7)')      ' nlyr, nfrq', nlyr, nfrq
         WRITE(*,'(  A,24G14.6)') ' res  ', res
         WRITE(*,'(  A,24G14.6)') ' reps ', reps
         WRITE(*,'(  A,24G14.6)') ' rmu  ', rmu
         WRITE(*,'(  A,24G14.6)') ' calf ', calf
         WRITE(*,'(  A,24G14.6)') ' ctau ', ctau
         WRITE(*,'(  A,24G14.6)') ' cfrq ', cfreq
         WRITE(*,'(  A,24G14.6)') ' thk  ', thk
         WRITE(*,'(  A,24G14.6)') ' alt  ', alt
     ENDIF

!      DEALLOCATE (mcurrent)

     CALL aem1d_fd(nfrq,freq,txcln,txa90,alt,zrx,xrx,yrx,  &
     &                  nlyr,res,reps,rmu,thk,calf,ctau,cfreq,bfd)


!  Store maximally coupled components in XMODL
      DO jf = 1 , nfrq
!            IF ( txa90 ) THEN
!               xbfd = norm(jf)*bfd(jf,2)
!            ELSE
               xbfd = norm(jf)                                          &
     &                *(bfd(jf,1)*sntx(jf)+bfd(jf,3)*cstx(jf))
!            ENDIF
            calc_data(jf)      = realpart(xbfd)
            calc_data(jf+nfrq) = imagpart(xbfd)
      ENDDO


       IF (debug) THEN
         WRITE(*,*) 'FWD-AEM05'
         WRITE(*,'(/A/4G14.6/4G14.6)') ' result is: ',calc_data
!         WRITE(*,*) ' result is: ',bfd
      ENDIF

END SUBROUTINE aemfwd1d_aem05

SUBROUTINE setup_aem05(nfrq,xrx,yrx,zrx,txcln,txa90,prm_fd,ppfac,norm)
!--------------------------------------------------------------------
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
!------------------------------------------------------------------------
! LAST CHANGE:   11 May 2016   VR
!------------------------------------------------------------------------
!

      IMPLICIT NONE
      INTEGER nfrq , jf
      REAL(KIND=8) sntx , cstx , xbd , ybd , zbd , rbrd , rsq , rsq1 , bfac ,   &
     &     fac , faczx , inline , vert , ppfac
      REAL(KIND=8) , DIMENSION(nfrq) :: txcln , xrx , yrx , zrx , prm_fd , norm
      LOGICAL coplanar , txa90

!  BFAC = 1.0E9 * MU / (4 * PI)  ! NANOTESLAS
      bfac = 100.

      prm_fd = 0.


      DO jf = 1 , nfrq
         sntx = dsin(txcln(jf))
         cstx = dcos(txcln(jf))


         xbd = -xrx(jf)
                   !  XRX is defined as positive behind the TX.
         ybd = yrx(jf)
         rbrd = dsqrt(xbd**2+ybd**2)
         zbd = zrx(jf)
         rsq = zbd**2 + rbrd**2
         rsq1 = dsqrt(rsq)
         coplanar = .FALSE.
         IF ( dabs(sntx)<.01 ) coplanar = .TRUE.
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
         norm(jf) = ppfac/dabs(prm_fd(jf))
      ENDDO

END SUBROUTINE setup_aem05
