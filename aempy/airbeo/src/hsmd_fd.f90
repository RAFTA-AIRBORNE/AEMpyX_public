!*==hsmd_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE hsmd_fd(nfrq,freq,alt,nrxf,txclnd,txa90,zrxd,xrxd,yrxd,&
     &                   nlyr,thk,res,reps,rmu,calf,ctau,cfreq,tdfd,    &
     &                   bfdd)
!------------------------------------------------------------------------
 
!***  Called by: HSBOSS_TD, HSBOSS_FD
!***      Calls: HSMD_HNK
 
!  Computes the frequency-domain layered earth magnetic field for a dipole of
!  unit moment and current.
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  SIGN CONVENTION:
!  ----------------
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 
!                             INPUT
!                             -----
!      FREQ - array of NFRQ frequencies
!       ALT - transmitter altitude above earth
!      NRXF - receiver offset dimension = NFRQ in frequency-domain and 1 in time-domain
!    TXCLND - angle in radians (QL) that TX dipole makes with vertical (climb = +)
!     TXA90 - true for vertical co-planar briadside array
!      ZRXD - vertical receiver offset (QL) from transmitter, (below = +)
!      XRXD - in-line horizontal offset (QL)  of RX J;  (behind = +)
!      YRXD - transverse horizontal offset (QL) of RX J (port = +)
!      NLYR - number of layers
!       THK - array of layer thicknesses
!       RES - array of layer resistivities
!      REPS - array of relative dielectric constants for each layer
!       RMU - array of relative magnetic permeabilities for each layer
!     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
!     TDFD  = 1 for time domain;  = 2 for frequency-domain
!
!                             OUTPUT
!                             ------
!   BFDD(1:NFRQ,3) - vertical magnetic field in nT
!   BFDD(1:NFRQ,1) - in-line magnetic field in nT
!   BFDD(1:NFRQ,2) - transverse magnetic field in nT
 
 
      IMPLICIT NONE
 
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80)
      REAL , PARAMETER :: Twopi = 6.2831853 , C_light = 2.99793E8
      COMPLEX(KIND=Ql) , PARAMETER :: Zero = (0._QL,0._QL)
      INTEGER nrxf , nfrq , nlyr , icole(nlyr) , jf , jl , tdfd , jq
      REAL w , res(nlyr) , freq(nfrq) , alt
      REAL(KIND=Ql) sig0(nlyr) , thkd(nlyr-1) , rmux(nlyr) , sntx ,     &
     &              cstx , zrfd , rhod , xbrq , xbrq2 , ybrq
      REAL(KIND=Ql) , DIMENSION(nrxf) :: txclnd , zrxd , xrxd , yrxd
      REAL , DIMENSION(nlyr) :: reps , rmu , ctau , cfreq , calf
      REAL , DIMENSION(nlyr-1) :: thk
      COMPLEX(KIND=Ql) iw , dispd , hlyr(3) , bfdd(nfrq,3) , vert ,     &
     &                 inline , trans
      LOGICAL txa90
 
      icole = 0
      DO jl = 1 , nlyr
         IF ( cfreq(jl)>1.E-3 .AND. ctau(jl)>1.E-12 ) icole(jl) = 1
      ENDDO
 
! Set extended precision variables
 
      sig0(1:nlyr) = real(1./res(1:nlyr),kind=Ql)
      thkd = 0._QL
      thkd(1:nlyr-1) = real(thk(1:nlyr-1),kind=Ql)
      rmux = real(rmu,kind=Ql)
 
!  Compute layered earth fields BLE_LYR at first station for each different altitude.
 
      bfdd = Zero
 
      jq = 1
      DO jf = 1 , nfrq
         IF ( tdfd==2 ) jq = jf
         zrfd = real(2.*alt,kind=Ql) - zrxd(jq)
                                              ! Reflected distance from TX to ground to RX
         rhod = sqrt(xrxd(jq)**2+yrxd(jq)**2)
         IF ( rhod>.01_QL ) THEN
            xbrq = -real(xrxd(jq))/rhod
                                     !  XRXD is defined + in negative direction
            ybrq = real(yrxd(jq))/rhod
         ELSE
            rhod = .01_QL
            xbrq = 0._QL
            ybrq = 0._QL
         ENDIF
         xbrq2 = xbrq**2
 
         w = Twopi*freq(jf)
         iw = cmplx(0.D0,w,kind=Ql)
         dispd = (iw/C_light)**2
         CALL hsmd_hnk(iw,nlyr,thkd,dispd,sig0,reps,rmux,icole,calf,    &
     &                 ctau,cfreq,zrfd,rhod,hlyr)
 
         IF ( txa90 ) THEN
            bfdd(jf,2) = hlyr(3)
 
         ELSE
 
            sntx = sin(txclnd(jq))
            cstx = cos(txclnd(jq))
 
            vert = (cstx*hlyr(1)) + (xbrq*sntx*hlyr(2))
 
            inline = sntx*((1._QL-2._QL*xbrq2)*hlyr(3)+xbrq2*hlyr(1))   &
     &               - cstx*xbrq*hlyr(2)
 
            trans = sntx*xbrq*ybrq*(hlyr(1)-2._QL*hlyr(3))              &
     &              - cstx*ybrq*hlyr(2)
            bfdd(jf,3) = vert
            bfdd(jf,1) = inline
            bfdd(jf,2) = trans
 
         ENDIF
      ENDDO
 
      END SUBROUTINE hsmd_fd
