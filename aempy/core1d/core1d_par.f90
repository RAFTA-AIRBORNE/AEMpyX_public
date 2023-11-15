
SUBROUTINE aem1d_fd(nfrq,freq,txcln,txa90,alt,zrx,xrx,   &
     &                     yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,   &
     &                     bfd)
!----------------------------------------------------------------------
!  Computes the frequency-domain layered earth H field for a dipole of
!  unit moment and current.
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
!       RES - layer resistivities
!      REPS - array of relative dislectric constants
!      RMUX - mu(i) / mu(0)
!       THK - array of layer thicknesses
!     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
!
!                             OUTPUT
!                             ------
!     BFD(JF,JS,1) - the in-line component of the layered earth response at
!                    time JT, station JS. (nT)
!     BFD(JF,JS,2) - the transverse component
!     BFD(JF,JS,3) - the vertical component
!
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   23 May 2016   VR
!------------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER nfrq , nlyr , tdfd, nrxf
      REAL(KIND=8)   alt
      REAL(KIND=8) , DIMENSION(nfrq) :: freq , txcln , zrx , xrx , yrx
      REAL(KIND=8) , DIMENSION(nlyr) :: res , rmu , reps , ctau , cfreq , calf
      REAL(KIND=8) , DIMENSION(nlyr) :: thk
      COMPLEX(KIND=8) bfd(nfrq,3)
      LOGICAL txa90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      NAMELIST /mynmlfd/ nfrq,freq,txcln,txa90,alt,zrx,xrx,   &
!     &                     yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,   &
!     &                     bfd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Compute layered earth fields

      tdfd = 2
      nrxf = nfrq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           WRITE(*,*) 'THIS IS mynmlfd before hsmd_fd '
!           WRITE(*,nml=mynmlfd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL hsmd_fd(nfrq,freq,alt,nrxf,txcln,txa90,zrx,xrx,yrx,nlyr, &
     &             thk,res,reps,rmu,calf,ctau,cfreq,tdfd,bfd)


      END SUBROUTINE aem1d_fd

SUBROUTINE aem1d_td(step,ider,nsx,swx,swy,npuls,pulse,ntypls,    &
     &                     ntyrp,trp,nchnl,topn,tcls,txcln,            &
     &                     alt,zrx, xrx,yrx,                           &
     &                     nlyr,res,reps,rmu,thk,calf,ctau, cfreq,     &
     &                     gstrp,astrp,td_out)
!------------------------------------------------------------------------------------
!  Computes BTD, the time-domain layered earth response convolved with the
!  excitation waveform and the receiver channels per unit receiver area.
!  For impulse response, it computes dB/dt in nT / s which is the same as
!  nanovolts per unit area.
!
!                             INPUT
!                             -----
!         JS = station reference
!       STEP = 1 iff step response is to be computed
!       IDER = 1 if source waveform was dB/dt; = 0 if amps pr B
!        NSX - number of points used to discretise transmitter signal
!        SWX - abscissae (seconds) of current waveform
!        SWY - dI/dt * Tx moment & nanotesla conversion at times SWX
!      NPULS - number of bipolar pulses of length PULSE
!      PULSE - length of half-cycle on pulse plus off-time
!      NTYRP - number of values in TRP for total signal length: 2 * NPULS *PULSE
!        TRP - array of time values for FD -> TD transformations
!     NTYPLS - number of TRP values in 1 PULSE
!      NCHNL - number of channels
!       TOPN - time at which receiver channel I opens.
!       TCLS - time at which receiver channel I closes.
!      TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!      NSTAT - number of stations in survey line.
!         SZ - array of transmitter altitudes
!        ZRX - vertical offset of RX at each station from transmitter  (below = +)
!        XRX - in-line horizontal offset of RX at each stationJ;  (behind = +)
!        YRX - transverse offset of RX at each stationJ (port = +)
!       NLYR - number of layers
!        RES - array of layer resistivities
!       REPS - relative dielectric constant
!        RMU - mu(i) / mu(0)
!        THK - array of layer thicknesses
!     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
!
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight patREAL*8, INTENT(in)h, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   25 Sept 2016   VR
!------------------------------------------------------------------------
!
      USE frequency_select

      IMPLICIT NONE
      INTEGER, PARAMETER :: nfrq = nf_6pde
      INTEGER, PARAMETER :: nrxf =1
      REAL(KIND=8) , PARAMETER :: twopi = 6.283185307D0
      INTEGER step , ider , nsx , npuls , ntypls , ntyrp , nchnl ,            &
     &        nlyr , tdfd , gstrp , astrp ,  jf , jt , jc
      REAL(KIND=8) swx(nsx) , swy(nsx,3) , pulse , trp(ntyrp) , alt,           &
     &     t , yprm(4,ntyrp) , yfrq(4,nfrq), costrn , td_out(nchnl,3)

      REAL(KIND=8) , DIMENSION(nfrq) :: freq , wf
      REAL(KIND=8) , DIMENSION(nrxf) :: txcln , zrx , xrx , yrx
      REAL(KIND=8), DIMENSION(nchnl) :: topn , tcls, ycum
      REAL(KIND=8) , DIMENSION(nlyr) :: res , reps , ctau , cfreq , calf , rmu
      REAL(KIND=8) , DIMENSION(nlyr) :: thk
      !REAL(KIND=8)  :: zrx , xrx , yrx
      COMPLEX(KIND=8) bfdd(nfrq,3)
      LOGICAL txa90

      tdfd = 1
      txa90 = .FALSE.


      freq(1:nfrq) = frq_6pde(1:nfrq)
      wf(1:nfrq) = dlog(twopi*freq(1:nfrq))


      CALL hsmd_fd(nfrq,freq,alt,nrxf,txcln,txa90,zrx,xrx,yrx,nlyr, &
     &             thk,res,reps,rmu,calf,ctau,cfreq,tdfd,bfdd)

!  Compute BTD, the 'observed' layered earth response by folding the BLEXT,
!  the extended response over NPULS bipolar cycles into 1 PULSE and then
!  convolving this with the TX waveform.  It is during the convolution that we
!  shift from teslas to nanoteslas or nT/s.

      yfrq = 0.d0
      DO jc = 1 , 3
         DO jf = 1 , nfrq
            yfrq(1,jf) = dimag(bfdd(jf,jc))/(twopi*freq(jf))
         ENDDO
         CALL cubspl(wf,yfrq,nfrq,0,0)

         yprm = 0.
         DO jt = 1 , ntyrp
                      !  Convert to step-function time-domain.
            t = trp(jt)
            yprm(1,jt) = costrn(wf,yfrq,nfrq,nfrq,t)

         ENDDO

         CALL fold_and_convolve(step,ider,nsx,swx,swy,npuls,pulse,trp,  &
     &                          ntypls,ntyrp,nchnl,topn,tcls,yprm,gstrp,&
     &                          astrp,ycum)

         td_out(1:nchnl,jc) = ycum(1:nchnl)
        ! WRITE(*,'(/i8/11G12.4/)') jc,ycum(1:nchnl)
      ENDDO
!       WRITE(*,*) 'SUBROUTINE aem1d_td'
!       WRITE(*,*) 'in-line: '
!       WRITE(*,*)  bfdd(:,1)
!       WRITE(*,*)'cross-line: '
!       WRITE(*,*)  bfdd(:,2)
!       WRITE(*,*) 'vertical: '
!       WRITE(*,*)  bfdd(:,3)
!       WRITE(*,*) 'bfdd',
!       WRITE(*,*)  bfdd


      END SUBROUTINE aem1d_td

SUBROUTINE hsmd_fd(nfrq,freq,alt,nrxf,txcln,txa90,zrx,xrx,yrx,&
     &                   nlyr,thk,res,reps,rmu,calf,ctau,cfreq,tdfd,    &
     &                   bfdd)
!------------------------------------------------------------------------

!***  Called by: HSBOSS_TD, HSBOSS_FD
!***      Calls: HSMD_HNK

!  Computes the frequency-domain layered earth magnetic field for a dipole of
!  unit moment and current.
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
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   14 May 2016   VR
!------------------------------------------------------------------------
!
      USE OMP_LIB


      IMPLICIT NONE

      INTEGER :: NT = 4

      REAL(KIND=8) , PARAMETER :: Twopi = 6.2831853D0 , c_light = 2.99793D8
      INTEGER  nfrq , nlyr , icole(nlyr) , jf , jl , tdfd , jq, nrxf, thrd
      REAL(KIND=8)  w , alt, sntx ,cstx , zrfd , rhod , xbrq , xbrq2 , ybrq
      REAL(KIND=8) , DIMENSION(nrxf) :: txcln , zrx , xrx , yrx
      REAL(KIND=8) , DIMENSION(nfrq) :: freq
      REAL(KIND=8) , DIMENSION(nlyr) :: res, sig, reps , rmu,     &
     &                  ctau , cfreq , calf
      REAL(KIND=8) , DIMENSION(nlyr) :: thk
      COMPLEX(KIND=8) :: iw , dispd , hlyr(3) , bfdd(nfrq,3) , vert ,     &
     &                 inline , trans, Zero
      LOGICAL txa90

      Zero = complex(0.D0,0.D0)

      icole = 0
      DO jl = 1 , nlyr
         IF ( cfreq(jl)>1.D-3 .AND. ctau(jl)>1.D-12 ) icole(jl) = 1
      ENDDO

! Set sigma
      sig(1:nlyr) = 1./res(1:nlyr)


!  Compute layered earth fields

      bfdd = Zero

      WRITE(*,*) 'freq # ',nfrq
      jq = 1
      !$omp parallel do num_threads(NT)
      !!!!private(jf, jq, xrx, yrx, zrx,iw, dispd, zrfd, rhod, hlyr, bfdd, sntx, cstx, vert, inline, trans)
      DO jf = 1 , nfrq
         thrd = omp_get_thread_num()
!          WRITE(*,*) 'thread # ',thrd

         IF ( tdfd==2 ) jq = jf
         ! write(*,*) jq
         zrfd = 2.d0*alt - zrx(jq)     ! Reflected distance from TX to ground to RX
         rhod = dsqrt(xrx(jq)**2+yrx(jq)**2)
         IF ( rhod>1.D-2 ) THEN
            xbrq = -xrx(jq)/rhod       !  XRXD is defined + in negative direction
            ybrq =  yrx(jq)/rhod
         ELSE
            rhod = 1.d-2
            xbrq = 0.d0
            ybrq = 0.d0
         ENDIF
         xbrq2 = xbrq**2

         w = Twopi*freq(jf)
         iw = dcmplx(0.D0,w)
         dispd = (iw/c_light)**2
         CALL hsmd_hnk(iw,nlyr,thk,dispd,sig,reps,rmu,icole,calf,    &
     &                 ctau,cfreq,zrfd,rhod,hlyr)

         IF ( txa90 ) THEN
            bfdd(jf,2) = hlyr(3)
         ELSE

            sntx = dsin(txcln(jq))
            cstx = dcos(txcln(jq))

            vert = (cstx*hlyr(1)) + (xbrq*sntx*hlyr(2))

            inline = sntx*((1.D0-2.D0*xbrq2)*hlyr(3)+xbrq2*hlyr(1))   &
     &               - cstx*xbrq*hlyr(2)

            trans = sntx*xbrq*ybrq*(hlyr(1)-2.D0*hlyr(3))              &
     &              - cstx*ybrq*hlyr(2)
            bfdd(jf,3) = vert
            bfdd(jf,1) = inline
            bfdd(jf,2) = trans


         ENDIF
      ENDDO

      END SUBROUTINE hsmd_fd

SUBROUTINE hsmd_hnk(iw,nlyr,thkd,dispd,sig0,reps,rmux,icole,calf, &
     &                    ctau,cfreq,zrfd,rhod,hlyr)

!  Magnetic dipole transmitter & receiver above or on earth surface
!  VALID FOR ANY NUMBER OF LAYERS
!  Computes transform integrals HLYR(3) which are used to compute vertical
!  and horizontal frequency-domain magnetic field components at the RX from
!  VMD and HMD sources.  It evaluates the Hankel transform integral using a
!  15 points per decade filter coefficient set derived from Christensen's
!  FLTGEN program.

!      IW - iw  angular frequency *(0.,1.)
!    RHOD - horizontal TX -> RX distance.
!     KER - stores kernel values from HSMD_KER
!
!
!    OUTPUT is HLYR(1:3)  forward model components
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!
!------------------------------------------------------------------------
! LAST CHANGE:   14 May 2016   VR
!------------------------------------------------------------------------
!

      USE filt_coef_q

      IMPLICIT NONE
      REAL(KIND=8) , PARAMETER :: Vfac0 = 100.D0
      INTEGER nlyr , i , icole(nlyr)
      REAL(KIND=8) , DIMENSION(nlyr) :: reps , ctau , cfreq , calf
      REAL(KIND=8) del_jn , rho_jn , y , lmbda , rhod , zrfd ,         &
     &              sig0(nlyr) , thkd(nlyr) , rmux(nlyr)
      COMPLEX(KIND=8) iw , hlyr(3) , dispd , qfd
      LOGICAL jump

      del_jn = dlog(10.D0)/dble(Ndec_jn)
      rho_jn = -dlog(rhod) - shftjn

      hlyr = (0.D0,0.D0)

      DO i = -50 , Jnhi       ! Start at I = -50 to pick up low values.
         y = rho_jn + dble(i)*del_jn
         lmbda = dexp(y)
         CALL hsmd_ker(iw,lmbda,nlyr,thkd,dispd,sig0,reps,rmux,icole,   &
     &                 calf,ctau,cfreq,zrfd,qfd)
         CALL hs_frq_jmp
         IF ( jump .AND. i>-40 ) EXIT
      ENDDO

      jump = .FALSE.      ! Finish off the low end for RHOTRP(1)
      DO i = -51 , Jnlo , -1
         y = rho_jn + dble(i)*del_jn
         lmbda = dexp(y)
         CALL hsmd_ker(iw,lmbda,nlyr,thkd,dispd,sig0,reps,rmux,icole,   &
     &                 calf,ctau,cfreq,zrfd,qfd)

         CALL hs_frq_jmp
         IF ( jump .AND. i<-60 ) EXIT
      ENDDO

      hlyr = Vfac0*hlyr/rhod

      CONTAINS

      SUBROUTINE hs_frq_jmp
!  ---------------------

!***  Called by: HSMD_HNK

!  Accumulates function calls for the Hankel transformation &
!  checks convergence.

      REAL(KIND=8) , PARAMETER :: Tol = 1.D-6 , Tol2 = 1.D-35
      INTEGER jint
      REAL(KIND=8) qr , qi
      COMPLEX(KIND=8) fw(3)

      fw(1) = wj0(i)*qfd*lmbda
      fw(2) = wj1(i)*qfd*lmbda
      fw(3) = wj1(i)*qfd/rhod

      hlyr = hlyr + fw

      jump = .TRUE.
      DO jint = 1 , 3
         qr = dabs(realpart(hlyr(jint)))
         qi = dabs(imagpart(hlyr(jint)))
         IF ( qr>Tol2 .AND. dabs(realpart(fw(jint)))>Tol*qr ) jump = .FALSE.
         IF ( qi>Tol2 .AND. dabs(imagpart(fw(jint)))>Tol*qi ) jump = .FALSE.
      ENDDO

      END SUBROUTINE hs_frq_jmp

      END SUBROUTINE hsmd_hnk

SUBROUTINE hsmd_ker(iw,lmbda,nlyr,thkd,dispd,sig0,reps,rmux,icole,&
     &                    calf,ctau,cfreq,zrfd,qfd)
!-----------------------------------------------------------------------------

!***  Called by: HSMD_SNTR

!  Kernel for dipole transmitter and receiver above earth.
!
!          Input
!          -----
!      IW - iw = complex angular frequency
!   LMBDA = Hankel transform variable
!    NLYR - number of layers
!    SIG0 = real conductivities ( 1. / RES)
!    REPS - array of relative dielectric constants
!    RMUX - mu(i) / mu(0)
!    THKD - layer thicknesses
!   ICOLE - C-C layer indicator array. (0 for pure real, 1 for C-C layer)
!    CALF - complementary chargeability array; ie., CALF(I) = 1.0 - CHRG(I)
!    CTAU - array of layer relaxation times (sec).
!   CFREQ - array of layer frequency parameters.
!    ZRFD - reflected distance from transmitter to earth to receiver
!
!          Output
!          ------
!  QFD  forward model kernel
!
!------------------------------------------------------------------------
! LAST CHANGE:   14 May 2016   VR
!------------------------------------------------------------------------


      IMPLICIT NONE

      REAL(KIND=8)  , PARAMETER :: Eps0 = 8.854156D-12 ,                &
     &                       Mu0 = 12.56637D-7 , Exp_tol = 80.D0
      INTEGER nlyr , icole(nlyr) , j
      REAL(KIND=8) , DIMENSION(nlyr) :: calf , ctau , cfreq , reps
      REAL(KIND=8)   xp0 , sig0(nlyr) , thkd(nlyr) , lmbda , zrfd ,    &
     &              rmux(nlyr) , rmusq(nlyr)
      COMPLEX(KIND=8) dispd , s0 , t0 , iw , p , lmbsq , ep ,          &
     &                 sigl(nlyr) , t(nlyr) , qfd, One, Zero
      COMPLEX(KIND=8) , DIMENSION(nlyr) :: f , xp , ksq , s , e


      One  = complex(1.D0,0.D0)
      Zero = complex(0.D0,0.D0)

      ksq = Zero
      e = Zero
      xp = Zero
      t = Zero

      lmbsq = dcmplx(lmbda*lmbda,0.D0)
      s0 = cdsqrt(lmbsq-dispd)
      xp0 = -lmbda*dexp(-lmbda*zrfd)

      sigl(1:nlyr) = dcmplx(sig0(1:nlyr),0.D0)
      DO j = nlyr , 1 , -1
         rmusq(j) = rmux(j)*rmux(j)
         p = (iw*ctau(j))**cfreq(j)
         p = icole(j)*p
         sigl(j) = sigl(j)*(One+p)/(One+calf(j)*p)
         sigl(j) = sigl(j) + iw*Eps0*reps(j)
                                            !  Add in displacement term
         ksq(j) = iw*Mu0*rmux(j)*sigl(j)
         s(j) = cdsqrt(ksq(j)+lmbsq)

         IF ( j==nlyr ) CYCLE

         ep = 2.D0*s(j)*thkd(j)
         IF ( realpart(ep)<Exp_tol ) e(j) = cdexp(-ep)
         t(j) = ((rmusq(j+1)-rmusq(j))*lmbsq+rmusq(j+1)*ksq(j)-rmusq(j) &
     &          *ksq(j+1))/(rmux(j+1)*s(j)+rmux(j)*s(j+1))**2
      ENDDO

      t0 = ((rmusq(1)-1.D0)*lmbsq-ksq(1))/(rmux(1)*s0+s(1))**2
      xp(1:nlyr-1) = e(1:nlyr-1)
      f = Zero
      DO j = nlyr - 1 , 1 , -1
         f(j) = xp(j)*(t(j)+f(j+1))/(One+t(j)*f(j+1))
      ENDDO

      qfd = xp0*(t0+f(1))/(One+t0*f(1))

      END SUBROUTINE hsmd_ker

SUBROUTINE fold_and_convolve(step,ider,nsx,swx,swy,npuls,pulse,   &
     &                             trp,ntypls,ntyrp,nchnl,topn,tcls,    &
     &                             ypls,gstrp,astrp,ycum)
!-------------------------------------------------------------------------------

!  Computes the "observed" response YCUM by convolving the splined earth
!  response function, YPLS, with the TX waveform.

!***  Called by: HSBOSS_TD, TDEM3D
!***      Calls: CUBDER, CUBVAL, CUBSPL, TXCNVD, TXCNVL, TQSTRIP

!     IDER - derivative indicator
!      NSX - number of points used to discretise transmitter signal
!      SWX - abscissae (seconds) of current waveform
!      SWY - dI/dt at times SWX
!    NPULS - number of bipolar pulses of length PULSE
!    PULSE - length single on pulse plus off-time
!    NTYRP - number of values in TRP for total signal length: 2 * NPULS *PULSE
!      TRP - array of time values for FD -> TD transformations
!   NTYPLS - number of TRP values in 1 PULSE
!    NCHNL - number of channels
!     TOPN - time at which receiver channel I opens.
!     TCLS - time at which receiver channel I closes.
!    GSTRP = 1 if Geotem / Questem stripping is to be applied

      IMPLICIT NONE
      INTEGER jt , ntyrp , ntypls , npuls , ider , step , nsx , nchnl , &
     &        jgl , jp , gstrp , astrp , mxcnv , ipl
      REAL(KIND=8) pulse , trp(ntyrp) , swx(nsx) , swy(nsx,3) , topn(nchnl) , &
     &     tcls(nchnl) , t1 , t2 , width , tf , tfh , hwidth , yc1 ,    &
     &     tc(3) ,ypls(4,ntyrp) , x , xp ,                              &
     &     ycum(nchnl) ,  wt , fold(ntypls) , ycnv(4,nsx)
      REAL(KIND=8) :: cubval , cubder , txcnvl , txcnvd
      REAL(KIND=8), DIMENSION(3) ::  glx = (/ -.7745967d0 , 0.d0 , .7745967d0/)
      REAL(KIND=8), DIMENSION(3) ::  glw = (/.5555556d0 , .8888889d0 , .5555556d0/)

      INTENT (in)ider , nsx , swx , swy , npuls , pulse , trp , ntypls ,&
     &        ntyrp , nchnl , topn , tcls , gstrp
      INTENT (inout) ypls
      INTENT (out)   ycum

!  Accumulate the results of NPULS bipolar cycles by splining the instantaneous
!  response and folding the positive and negative parts of each cycle back
!  into a single pulse.


      CALL cubspl(trp,ypls,ntyrp,0,0)

      fold = 0.d0
      ipl = 1
      IF ( npuls==1 ) ipl = 0

      IF ( step==1 ) THEN
         DO jt = 1 , ntypls
            x = trp(jt)
            xp = x + pulse
            fold(jt) = cubval(trp,ypls,ntyrp,x)                         &
     &                 - ipl*cubval(trp,ypls,ntyrp,xp)
            DO jp = 2 , npuls
               x = xp + pulse
               xp = x + pulse
               fold(jt) = cubval(trp,ypls,ntyrp,x)                      &
     &                    - cubval(trp,ypls,ntyrp,xp) + fold(jt)
            ENDDO
         ENDDO
      ELSE
         DO jt = 1 , ntypls
            x = trp(jt)
            xp = x + pulse
            fold(jt) = ipl*cubder(trp,ypls,ntyrp,xp)                    &
     &                 - cubder(trp,ypls,ntyrp,x)
            DO jp = 2 , npuls
               x = xp + pulse
               xp = x + pulse
               fold(jt) = cubder(trp,ypls,ntyrp,xp)                     &
     &                    - cubder(trp,ypls,ntyrp,x) + fold(jt)
            ENDDO
         ENDDO
      ENDIF

      ypls = 0.d0
      ypls(1,1:ntypls) = fold(1:ntypls)
      CALL cubspl(trp,ypls,ntypls,0,0)
      ycum = 0.d0

!  Begin convolution.  If Geotem / Questem primary field stripping is required
!  the convolution must be done for all points in the waveform.
!  Otherwise, convolve only for those points needed in the windows.

!  The layered earth field is in IMPULSE form if dB/dt is desired
!  or in STEP form if B is to be computed.

      mxcnv = ntypls + nsx
      tf = swx(nsx)
      tfh = 0.5d0*tf

      IF ( gstrp==1 ) CALL tqstrip(ider,ntypls,trp,ypls,nsx,swx,swy,    &
     &                             ycnv)
      DO jt = 1 , nchnl
         t1 = topn(jt)
         t2 = tcls(jt)
         width = t2 - t1
         hwidth = width/2.d0

! Step response for step input or impulse response response for impulse input
! Average the response over receiver windows using 3 point Gaussian integration.

         tc(2) = (tcls(jt)+topn(jt))/2.
         tc(1) = tc(2) + hwidth*glx(1)
         tc(3) = tc(2) + hwidth*glx(3)

         DO jgl = 1 , 3
            t1 = tc(jgl)
            wt = glw(jgl)/2.d0
            yc1 = 0.d0

            IF ( gstrp==1 .AND. t1<tf ) THEN
               yc1 = cubval(swx,ycnv,nsx,t1)
            ELSEIF ( ider==0 ) THEN
                                ! Waveform input as I or B field (derived dI/dt)
               yc1 = txcnvl(t1,ntypls,trp,ypls,nsx,swx,swy)
               IF ( astrp==1 .AND. t1<tf ) THEN
                                            ! Asymmetric stripping
                  t2 = t1 - tfh
                  yc1 = yc1 + txcnvl(t2,ntypls,trp,ypls,nsx,swx,swy)
               ENDIF

            ELSEIF ( ider==1 ) THEN  ! Waveform input as voltage (known dI/dt)
               yc1 = txcnvd(mxcnv,t1,ntypls,trp,ypls,nsx,swx,swy)
            ELSEIF ( ider==4 ) THEN                ! pure on-time step
               yc1 = swy(1,1)*cubval(trp,ypls,ntypls,t1)
            ENDIF
            ycum(jt) = ycum(jt) + (wt*yc1)
         ENDDO
      ENDDO

      END SUBROUTINE fold_and_convolve

REAL(KIND=8) FUNCTION costrn(wf,yfrq,nfrq,kfrq,t)
!------------------------------------------

!***  Called by: HSBOSS_TD, TDEM_3D
!***      Calls: CUBVAL

! LAST MODIFICATION DATE: October, 2001

! Produces time-domain value at time T by cosine transformation of NFRQ
! frequency-domain values contained in cubic spline array YFRQ.
! KFRQ is the high frequency cutoff, less than or equal to NFRQ.
! Array WF contains the LOG (base e) of the angular frequency values.

! The routine uses filter coefficients derived from the Niels Christensen
! fast Hankel transform routine FILCOA at a spacing of 12 points per decade
! and omega = 0.3.  Various filters were tested using a vertical magnetic
! dipole receiver in a very large circular for which accurate frequency
! and time-domain solutions were programmed.  This particular filter gave
! the overall best accuracy for 1/2 spaces ranging in resistivity from
! .1 to 10,000 ohm-m for times ranging from .01 to 50 msec.


!  K(W,T) = (2/PI) * F(W) * COS(WT) dW

! Letting X = WT, the above becomes
!
!  K(W,T) = (2/PI*T) * F(X/T) * COS(X) dX
!
! From Abramowitz and Stegun, COS(X) = SQRT(X*PI/2) * J(-1/2:X).
! Filter Coefficients are used to represent X**(1/2) * J(-1/2:X)
!
!  COSTRN = SQRT (2/PI) * SUM(i) { WCOS(i) * F [X(i) /T] }

! The accumulation is done using 12 digit precision


      USE filt_coef_q

      IMPLICIT NONE
      INTEGER , PARAMETER :: Ndec_cos = 12 , Kflow = -200 , Kfhigh = 99
      REAL(KIND=8) , PARAMETER :: Fac = .7978846d0 , Tol = 1.0d-6
      INTEGER j1 , nfrq , kfrq
      REAL(KIND=8) wf(nfrq) , yfrq(4,nfrq)
      REAL(KIND=8) delta , y1 , y , ytym , val, t , cubval , v1

      INTENT (in)  wf , yfrq , nfrq , t


      delta = dlog(10.d0)/dble(Ndec_cos)
      ytym = 0.
      y1 = -dlog(t) - delcos

! Begin right side convolution at weight 0.
! Stop when frequency domain array is exhausted.

      MOVE_HIGH:DO j1 = 0 , Kfhigh

         y = y1 + j1*delta
         IF ( y>wf(kfrq) ) EXIT MOVE_HIGH
         IF ( y<wf(1) ) y = wf(1)
         v1 = cubval(wf,yfrq,nfrq,y)
         val = wcos(j1)*v1
         ytym = ytym + val
      ENDDO MOVE_HIGH

      y = y1

! Begin left side convolution at weight -1.
! When log angular frequency is less than WF(3), check convergence.
! Continue left using the fact that impulse B is inversely proportional to
! frequency as freq -> 0; i.e., step response B is constant.

      MOVE_LOW:DO j1 = -1 , Kflow , -1

         y = y1 + j1*delta
         IF ( y>wf(kfrq) ) CYCLE MOVE_LOW
         IF ( y<wf(1) ) y = wf(1)
         v1 = cubval(wf,yfrq,nfrq,y)
         val = wcos(j1)*v1
         ytym = ytym + val
         IF ( (y<wf(3)) ) THEN
            IF ( dabs(val)<Tol*dabs(ytym) ) EXIT MOVE_LOW
         ENDIF
      ENDDO MOVE_LOW

      costrn = Fac*ytym/t

      END FUNCTION costrn

SUBROUTINE cubspl(xnot,c,n,ibcbeg,ibcend)
! ----------------------------------------------

!***  Called by: EGT_CSPL, FOLD_AND_CONVOLVE, HSBOSS_TD, INTER_EGT_CSPL, MGTBS,
!                PRM_BOSS, TDEM_3D, TQSTRIP, TXCNVD,
!

!  Calculates coefficients for cubic spline interpolation.
!  Call function CUBVAL to evaluate function values after interpolation.
!  From  * A PRACTICAL GUIDE TO SPLINES *  by Carl de Boor.

!             INPUT
!             -----
!
!     N = number of data points. assumed to be > 1.
!
!  (XNOT(I), C(1,I), I=1,...,N) = abscissae and ordinates of the data points.
!                                 XNOT is assumed to be strictly increasing.
!
!     IBCBEG, IBCEND = boundary condition indicators, and
!     C(2,1), C(2,N) = boundary condition information. Specifically,
!
!     IBCBEG = 0  No boundary condition at XNOT(1) is given.  In this case,
!                 the not-a-knot condition is used, i.e. the jump in the
!                 third derivative across XNOT(2) is forced to zero.  Thus
!                 first and the second cubic polynomial pieces are made to
!                 coincide.
!     IBCBEG = 1  the slope at XNOT(1) is made to equal C(2,1),
!                 supplied by input.
!     IBCBEG = 2  the second derivative at XNOT(1) is made to equal C(2,1),
!                 supplied by input.
!
!     IBCEND = 0, 1, or 2 has analogous meaning concerning the boundary
!                 condition at XNOT(n), with the additional information
!                 taken from C(2,n).
!
!          OUTPUT
!          ------
!
!     C(J,I), J=1,...,4; I=1,...,L (= N-1) = the polynomial coefficients
!         of the cubic interpolating spline with interior knots (or joints)
!         XNOT(2), ..., XNOT(N-1).
!
!        In the interval: (XNOT(I) - XNOT(I+1)), the spline F is given by:
!
!        F(X) = C(1,I) + H* (C(2,I) + H* (C(3,I) + H* C(4,I)/3.) /2.)
!
!     where H = X - XNOT(I).  FUNCTION  *CUBVAL* may be
!     used to evaluate F or its derivatives from XNOT,C, L = N-1,
!     AND K=4.
!******************************************************************************
!******************************************************************************
      IMPLICIT NONE
      INTEGER ibcbeg , ibcend , n , i , j , l , m
      REAL(KIND=8) c(4,n) , xnot(n) , divdf1 , divdf3 , dxnot , g

      INTENT (in)    xnot , n , ibcbeg , ibcend
      INTENT (inout) c
      SAVE

!  A tridiagonal linear system for the unknown slopes S(I) of F at
!  XNOT(I), I=1,...,N, is generated and then solved by Gauss elimination,
!  with S(I) ending up in C(2,I), ALL I.
!  C(3,.) AND C(4,.) are used initially for temporary storage.

!  Compute first differences of XNOT sequence and store in C(3,.).
!  Also, compute first divided difference of data and store in C(4,.).

      l = n - 1
      DO m = 2 , n
         c(3,m) = xnot(m) - xnot(m-1)
         c(4,m) = (c(1,m)-c(1,m-1))/c(3,m)
      ENDDO

!  Construct first equation from the boundary condition, of the form
!  C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)

      IF ( ibcbeg<1 ) THEN
         IF ( n>2 ) THEN

!  Not-a-knot condition at left end and N > 2.

            c(4,1) = c(3,3)
            c(3,1) = c(3,2) + c(3,3)
            c(2,1) = ((c(3,2)+2.*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))&
     &               /c(3,1)
            GOTO 100
         ELSE

!  No condition at left end and N = 2.

            c(4,1) = 1.
            c(3,1) = 1.
            c(2,1) = 2.*c(4,2)
            GOTO 400
         ENDIF
      ELSEIF ( ibcbeg==1 ) THEN

!  Slope prescribed at left end.

         c(4,1) = 1.d0
         c(3,1) = 0.d0
      ELSE

!  Second derivative prescribed at left end.

         c(4,1) = 2.d0
         c(3,1) = 1.d0
         c(2,1) = 3.d0*c(4,2) - c(3,2)*c(2,1)/2.
      ENDIF
      IF ( n==2 ) GOTO 400

!  if there are interior knots, generate the corresponding equations and
!  perform the forward pass of Gauss elimination, after which the M-TH
!  equation reads    C(4,M)*S(M) + C(3,M)*S(M+1) = C(2,M).

 100  DO m = 2 , l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
         c(4,m) = g*c(3,m-1) + 2.*(c(3,m)+c(3,m+1))
      ENDDO

!  Construct last equation from the second boundary condition, of the form
!  (-G*C(4,N-1))*S(N-1) + C(4,N)*S(N) = C(2,N)
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since C array happens to be set up just right for it
!  at this point.

      IF ( ibcend<1 ) THEN
         IF ( n/=3 .OR. ibcbeg/=0 ) THEN

!  Not-a-knot and N > 2, and either N > 3 or also not-a-knot at
!  left end point.

            g = c(3,n-1) + c(3,n)
            c(2,n) = ((c(3,n)+2.*g)*c(4,n)*c(3,n-1)+c(3,n)**2*(c(1,n-1)-&
     &               c(1,n-2))/c(3,n-1))/g
            g = -g/c(4,n-1)
            c(4,n) = c(3,n-1)
            GOTO 500
         ENDIF
      ELSEIF ( ibcend==1 ) THEN
         GOTO 600
      ELSE
         GOTO 300
      ENDIF

!  Either (N=3 and not-a-knot also at left) or (N=2 and not not-a-
!  knot at left end point).

 200  c(2,n) = 2.d0*c(4,n)
      c(4,n) = 1.d0
      g = -1./c(4,n-1)
      GOTO 500

!  Second derivative prescribed at right endpoint.

 300  c(2,n) = 3.d0*c(4,n) + c(3,n)*c(2,n)/2.
      c(4,n) = 2.d0
      g = -1./c(4,n-1)
      GOTO 500
 400  IF ( ibcend<1 ) THEN
         IF ( ibcbeg>0 ) GOTO 200

!  Not-a-knot at right endpoint and at left endpoint and N = 2.

         c(2,n) = c(4,n)
         GOTO 600
      ELSEIF ( ibcend==1 ) THEN
         GOTO 600
      ELSE
         GOTO 300
      ENDIF

!  Complete forward pass of Gauss elimination.

 500  c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)

!  Perform back substitution.

 600  j = l
      DO
         c(2,j) = (c(2,j)-c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         IF ( j<=0 ) THEN

!  Generate cubic coefficients in each interval, i.e., the derivatives at its
!  left endpoint, from value and slope at its endpoints.

            DO i = 2 , n
               dxnot = c(3,i)
               divdf1 = (c(1,i)-c(1,i-1))/dxnot
               divdf3 = c(2,i-1) + c(2,i) - 2.*divdf1
               c(3,i-1) = 2.*(divdf1-c(2,i-1)-divdf3)/dxnot
               c(4,i-1) = (divdf3/dxnot)*(6./dxnot)
            ENDDO
            EXIT
         ENDIF
      ENDDO
      END SUBROUTINE cubspl

REAL(KIND=8) FUNCTION cubint(xknot,coef,knot,x1,x2)
! ------------------------------------------------
!
!***  Called by:  EGT_BOSS TXCNVD, TXCNVL
!***      Calls: INTERV.  On exit from INTERV

!  Integrates a function from X1 to X2 using its cubic spline representation.

!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range

!      KNOT - total number of knots including endpoints.
!
!     XKNOT(I), I = 1,KNOT - Location of the knots.  The rightmost data
!                            point used to calculate coefficients is not
!                            included.
!
!     COEF(J,I), J = 1,4; I = 1,KNOT
!
!              The coefficients of the cubic spline represent the
!              indefinite integral of F, on the I'th interval, as:
!
!       INTGR [ F(X) ] = COEF(4,I)/24 * H**4  +  COEF(3,I)/6 * H**3  +
!                        COEF(2,I)/2 * H**2  +  COEF(1,I) * H
!
!                          WITH  H = X - XKNOT(K)
!
!  This is a modification of the FUNCTION PPVALU in the book
!  "A PRACTICAL GUIDE TO SPLINES"  by C. DE BOOR

!*********************************************************************

      IMPLICIT NONE
      INTEGER i , i1 , i2 , mflag , knot
      REAL(KIND=8) h , h1 , h2 , x1 , x2 , xknot(knot) , coef(4,knot)

!  Find the indices I1 and I2 of largest breakpoints to the left of X1
!  and X2 respectively.
!
      CALL interv(xknot,knot-1,x1,i1,mflag)
      CALL interv(xknot,knot-1,x2,i2,mflag)
      h1 = x1 - xknot(i1)
      IF ( mflag==-1 ) h1 = 0.d0

      h2 = x2 - xknot(i2)
      cubint = (((coef(4,i2)*h2/4.d0+coef(3,i2))*h2/3.d0+coef(2,i2))      &
     &         *h2/2.d0+coef(1,i2))                                      &
     &         *h2 - (((coef(4,i1)*h1/4.d0+coef(3,i1))*h1/3.d0+coef(2,i1))&
     &         *h1/2.d0+coef(1,i1))*h1

!  Include integrals over intervening intervals.

      IF ( i2>i1 ) THEN
         DO i = i1 , i2 - 1
            h = xknot(i+1) - xknot(i)
            cubint = cubint +                                           &
     &               (((coef(4,i)*h/4.d0+coef(3,i))*h/3.d0+coef(2,i))     &
     &               *h/2.d0+coef(1,i))*h
         ENDDO
      ENDIF

      END FUNCTION cubint

REAL(KIND=8) FUNCTION cubval(xknot,coef,knot,x1)
! --------------------------------------------
!***  Called by: COSTRN, EGT_BOSS, FOLD_AND_CONVOLVE, INTER_EGT_BOSS,
!                MGTV1, PRM_BOSS, TXCNVD, TXCNVL
!
!***      Calls: INTERV.

!  On exit from INTERV,
!  Evaluates a function at X1 from from its cubic spline representation.

!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range

!      KNOT - total number of knots including endpoints.
!
!     XKNOT(I), I = 1,KNOT - location of the knots.  The rightmost data
!                            point used to calculate coefficients is not
!                            included.
!
!     COEF(J,I), J = 1,4; I = 1,KNOT
!
! The coefficients of the cubic spline on the I'th interval represent F as:
!
!                F(X) = COEF(4,I)/6 * H**3  +  COEF(3,I)/2 * H**2  +
!                       COEF(2,I) * H  +  COEF(1,I)
!
!                          with  H = X - XKNOT(I)
!
!  This is a modification of the FUNCTION PPVALU in the book
!  "A PRACTICAL GUIDE TO SPLINES"  by C. DE Boor
!
!             METHOD
!             ------
!
!  The interval index I, appropriate for X, is found through a call to INTERV.
!  The formula for F is evaluated using nested multiplication.
!******************************************************************************

      IMPLICIT NONE
      INTEGER i , mflag , knot
      REAL(KIND=8) xknot(knot) , coef(4,knot) , x1 , h

      INTENT (in)  xknot , coef , knot , x1
!
!  Find index I of largest breakpoint to the left of X1.
!
      CALL interv(xknot,knot-1,x1,i,mflag)
      h = x1 - xknot(i)
      IF ( mflag==-1 ) h = 0.d0

      cubval = ((coef(4,i)*h/3.d0+coef(3,i))*0.5d0*h+coef(2,i))            &
     &         *h + coef(1,i)

      END FUNCTION cubval

REAL(KIND=8) FUNCTION cubder(xknot,coef,knot,x1)
! --------------------------------------------

!***  Called by: FOLD_AND_CONVOLVE
!***      Calls: INTERV.  On exit from INTERV

!  Evaluates the first derivative of a function from its cubic spline
!  interpolation.

!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range

!         KNOT - total number of knots including endpoints.
!     XKNOT(I), I = 1,KNOT - location of the knots.  The rightmost data
!                            point used to calculate coefficients is not
!                            used.
!     COEF(J,I), J = 1,4; I = 1,KNOT = Jth derivative at H = 0
!                                      where  H = X - XKNOT(I)
!******************************************************************************
!******************************************************************************

      IMPLICIT NONE
      INTEGER i , mflag , knot
      REAL(KIND=8) xknot(knot) , coef(4,knot) , x1 , h

!  Find index i of largest breakpoint to the left of X1.

      CALL interv(xknot,knot-1,x1,i,mflag)
      h = x1 - xknot(i)
      IF ( mflag==-1 ) h = 0.d0

      cubder = (coef(4,i)*h/2.+coef(3,i))*h + coef(2,i)

      END FUNCTION cubder

REAL(KIND=8) FUNCTION linval(nx,xval,yval,x1,ider)
!-------------------------------------------

!***  Called by: TXCNVD
!***      Calls: INTERV

!  Evaluates a function at X1 from from its linear representation.
!
!           On exit from INTERV
!
!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range

!
!     XVAL(1:NX) - location of the abscissa knots.  The rightmost data point
!                  used to calculate coefficients is not included.
!
!     YVAL(1:NX,1) = function values.
!     YVAL(1:NX,2)   may be populated but aren't used.
!
!     If IDER = 0, the value in the interval is that of the leftmost knot.
!                  because the derivative has been computed using two knot
!                  values and stored at the left node.
!
!     If IDER = 1, the value is a linear interpolation between the knots.
!
!             METHOD
!             ------
!
!  The interval index I, appropriate for X, is found through a call to INTERV.
!  The formula for F is evaluated using nested multiplication.
!******************************************************************************

      IMPLICIT NONE
      INTEGER i , mflag , nx , ider
      REAL(KIND=8) xval(nx) , yval(nx,3) , x1 , h

      INTENT (in)nx , xval , yval , x1 , ider
!
!  Find index I of largest breakpoint to the left of X1.
!
      CALL interv(xval,nx-1,x1,i,mflag)

      IF ( ider==0 ) THEN !  Computed derivative values stored at right node (26.01.00)
         linval = yval(i+1,1)
      ELSE
         h = x1 - xval(i)
         IF ( mflag==-1 ) h = 0.d0

         linval = yval(i,1) + h*(yval(i+1,1)-yval(i,1))                 &
     &            /(xval(i+1)-xval(i))
      ENDIF

      END FUNCTION linval

SUBROUTINE interv(xt,lxt,x,left,mflag)
!-------------------------------------------

!***   Called by: CUBVAL, CUBINT, CUBDER, LINVAL

!---  Restructured April, 1997

!  from  * A PRACTICAL GUIDE TO SPLINES *  by C. DE BOOR
!  computes  LEFT = MAX( I , 1 <= I <= LXT  .AND.  XT(I) <= X )  .
!
!             INPUT
!             -----
!       XT - a real sequence, of length  LXT, assumed to be non-decreasing.
!      LXT - number of terms in the sequence  XT .
!        X - the point whose location with respect to the sequence XT is
!            to be determined.
!
!             OUTPUT
!             ------
!      LEFT, MFLAG.....are both integers, whose value is:
!
!        1     -1      IF               X <  XT(1)
!        I      0      IF   XT(I)  <= X < XT(I+1)
!       LXT     1      IF  XT(LXT) <= X
!
!        In particular, MFLAG = 0 is the 'usual' case.  MFLAG /= 0
!        indicates that X  lies outside the halfopen interval
!        XT(1) <= Y < XT(LXT) . The asymmetric treatment of the
!        interval is due to the decision to make all pp functions
!        continuous from the right.
!
!             METHOD
!             ------
!
!  The program is designed to be efficient in the common situation that
!  it is called repeatedly, with  X  taken from an increasing or decreasing
!  sequence. This will happen, e.g., when a pp function is to be grapged.
!  The first guess for  LEFT  is therefore taken to be the value returned at
!  the previous call and stored in the  L O C A L  variable ILO. A first
!  check ascertains that  ILO < LXT (This is necessary since the present
!  call may have nothing to do with the previous call).
!  Then, if XT(ILO) <= XT(ILO+1),
!  we set  LEFT = ILO  and are done after just three comparisons.
!  Otherwise, we repeatedly double the difference  ISTEP = IHI - ILO
!  while also moving  ILO  AND  IHI  in the direction of  X , until
!                      XT(ILO) <= X < XT(IHI) ,
!  after which we use bisection to get, in addition, ILO+1 = IHI .
!  LEFT = ILO  is then returned.
!******************************************************************************
!******************************************************************************

      IMPLICIT NONE
      INTEGER left , lxt , mflag , ihi , ilo , istep , middle , j1
      REAL(KIND=8) x , xt(lxt)
      SAVE ilo

      DATA ilo/1/

!***********************************************************
!  Trivial returns when X is not in the range.

      IF ( (x<=xt(1)) .OR. (lxt<=1) ) THEN
         mflag = -1
         left = 1
         RETURN
      ENDIF

      IF ( x>=xt(lxt) ) THEN
         mflag = 1
         left = lxt
         RETURN
      ENDIF

      mflag = 0
      IF ( ilo>=lxt ) ilo = lxt - 1
      ihi = ilo + 1

!  Trivial return when X is already in the interval.

      IF ( (x<=xt(ihi)) .AND. (x>=xt(ilo)) ) THEN
         left = ilo
         RETURN
      ENDIF
!***********************************************************

      IF ( x<=xt(ilo) ) THEN
                          ! decrease ILO  to capture X.
         istep = 1
         DO j1 = 1 , lxt
            ihi = ilo
            ilo = ihi - istep
            ilo = max(1,ilo)
            IF ( (x>=xt(ilo)) .OR. (ilo==1) ) EXIT
            istep = istep*2
         ENDDO

      ELSEIF ( x>=xt(ihi) ) THEN
                                ! increase IHI to capture X

         istep = 1
         DO j1 = 1 , lxt
            ilo = ihi
            ihi = ilo + istep
            ihi = min(ihi,lxt)
            IF ( (x<=xt(ihi)) .OR. (ihi==lxt) ) EXIT
            istep = istep*2
         ENDDO

      ENDIF

!  Now XT(ILO) <= X < XT(IHI) . Narrow the interval.

      DO j1 = 1 , lxt
         middle = (ilo+ihi)/2
         IF ( middle==ilo ) EXIT
         IF ( x<xt(middle) ) THEN
            ihi = middle
         ELSE
            ilo = middle
         ENDIF
      ENDDO

! Task complete

      left = ilo

      END SUBROUTINE interv

SUBROUTINE txcmrg(mxcnv,x1,y1,n1,x2,y2,n2,xcnv,ycnv,ncnv)
!----------------------------------------------------------

!***  Called by: TXCNVD

!  Merges two previously sorted list pairs X1, Y1 of length N1 and X2, Y2 of
!  length N2 into list pair XCNV, YCNV of length NCNV into ascending values of
!  XCNV.

      IMPLICIT NONE
      REAL(KIND=8) , PARAMETER :: Tol = 1.d-3
      INTEGER mxcnv , n1 , n2 , ncnv , k1 , k2 , n , j1
      REAL(KIND=8) delt , tl1 , xcnv(mxcnv) , x1(mxcnv) , y1(mxcnv) , x2(mxcnv) &
     &     , y2(mxcnv) , ycnv(4,mxcnv)
      LOGICAL list1 , list2

      INTENT (in)  mxcnv , x1 , y1 , n1 , x2 , y2 , n2
      INTENT (out) xcnv , ycnv , ncnv

      list1 = .TRUE.
      list2 = .TRUE.
      k1 = 1
      k2 = 1
      n = n1 + n2

      DO j1 = 1 , n
         IF ( list1 .AND. list2 ) THEN
            IF ( x1(k1)<x2(k2) ) THEN
               xcnv(j1) = x1(k1)
               ycnv(1,j1) = y1(k1)
               k1 = k1 + 1
               IF ( k1>n1 ) list1 = .FALSE.
            ELSE
               xcnv(j1) = x2(k2)
               ycnv(1,j1) = y2(k2)
               k2 = k2 + 1
               IF ( k2>n2 ) list2 = .FALSE.
            ENDIF
         ELSEIF ( list1 ) THEN
            xcnv(j1) = x1(k1)
            ycnv(1,j1) = y1(k1)
            k1 = k1 + 1
            IF ( k1>n1 ) list1 = .FALSE.
         ELSEIF ( list2 ) THEN
            xcnv(j1) = x2(k2)
            ycnv(1,j1) = y2(k2)
            k2 = k2 + 1
            IF ( k2>n2 ) list2 = .FALSE.
         ENDIF
      ENDDO

      ncnv = 1 !  Clean up list
      DO j1 = 2 , n
         delt = xcnv(j1) - xcnv(ncnv)
         tl1 = Tol*xcnv(j1)
         IF ( delt>tl1 ) THEN
            ncnv = ncnv + 1
            xcnv(ncnv) = xcnv(j1)
            ycnv(1,ncnv) = ycnv(1,j1)
         ENDIF
      ENDDO

      END SUBROUTINE txcmrg

REAL(KIND=8) FUNCTION txcnvd(mxcnv,t,ntypls,trp,ypls,nsx,swx,swy)
!----------------------------------------------------------
!
!***  Called by: FOLD_AND_CONVOLVE, TQSTRIP
!***      Calls: CUBINT, CUBSPL, CUBVAL, LINVAL, TXCMRG

!  Convolves impulse B (step dB/dt) earth response function (ERF) with the
!  specified derivative of the source waveform at NSX points to produce
!  the system dB/dt response of the earth.
!
!       MXCNV = NTYPLS + NSX
!           T - convolution time in sec measured from the beginning
!               of the source waveform.
!   TRP, YPLS - abscissa & ordinate values of earth response function to
!               be convolved.
!      NTYPLS - number of values in TRP and YPLS
!         SWX - abscissa of time values of source waveform in sec.
!         SWY - dI/dt values derived from receiver dB/dt.
!         NSX - number of points in SWX & in each waveform stored in SWY
!
!  Defining  T1 = MIN {T, signal length,}, the convolution is formally
!  computed as
!
!   TXCNVD (T) = INT (T0 -> T) { YPLS (tau) * SWY (T-tau)  d tau }

!  where T0 = MAX { TRP(1), T - SWX (NSX)}
!
!       ONTIME RESPONSE
!       ---------------
!  For response in the on-time period, ( T < signal length) a correction to
!  account for the response from 0 -> T0 is needed.  Analysis and subsequent
!  numerical experiments confirm that as T -> 0, step dB/dt -> A * T**(-1/2).
!  Thus ERFINT, the integral of YPLS from 0 to TRP(1), is simply
!  2 * TRP(1) * YPLS (TRP(1)) if TRP(1) is chosen sufficiently early.
!  The convolution correction factor is SWY(T) * ERFINT.

!  Alternatively, we can difference the step B field from 0 to TRP(1) which
!  is a lot easier since the step B field at T = 0 is simply the DC field due
!  to a transmitter image buried at z = ALT; i.e., the z+z' term.  In this case,
!  the bigger TRP(1) is, the more accurate the difference in B but this must be
!  sufficiently small so that the change in dI/dt is negligable.  Thus, TRP(1)
!  is chosen to be .1 microsecond.

      IMPLICIT NONE
      INTEGER mxcnv , ntypls , nsx , n1 , j1 , n2 , j2 , ncnv
      REAL(KIND=8) t , tc , t0 , trp(ntypls) , ypls(4,ntypls) , swx(nsx) ,      &
     &     swy(nsx,3) , ycnv(4,mxcnv) , xcnv(mxcnv) , x1(mxcnv) ,       &
     &     y1(mxcnv) , x2(mxcnv) , y2(mxcnv) , cubval , cubint , linval

      INTENT (in) mxcnv , t , ntypls , trp , ypls , nsx , swx , swy

!  Set up X1,Y1, the N1 values of SWX, SWY * YPLS for signal ontime < T.
!  where X1, the conjugate signal time, contains T-SWX values.
!  Set up X2,Y2, the N2 values of TRP, YPLS * SWY for ERF points  <= T.

!  Return TXCNVD = 0 if N1 + N2 < 4 or if NCNV < 4

      txcnvd = 0.0d0
      n1 = 0
      DO j1 = nsx , 1 , -1
         tc = t - swx(j1)
         IF ( tc<0. ) CYCLE
         n1 = n1 + 1
         x1(n1) = tc
         y1(n1) = swy(j1,1)*cubval(trp,ypls,ntypls,tc)
      ENDDO

      t0 = t - swx(nsx)
      t0 = max(t0,trp(1))/1.0001

      n2 = 0
      DO j2 = 1 , ntypls
         IF ( (trp(j2)>t0) .AND. (trp(j2)<t) ) THEN
            n2 = n2 + 1
            x2(n2) = trp(j2)
            tc = t - trp(j2)
            y2(n2) = ypls(1,j2)*linval(nsx,swx,swy,tc,1)
         ENDIF
      ENDDO

!  Merge the two lists into XCNV, YCNV of length NCNV.
!  Then spline and integrate

!+++++++++++++++++++++++++++++++++
      IF ( n1+n2<4 ) RETURN
!+++++++++++++++++++++++++++++++++

      CALL txcmrg(mxcnv,x1,y1,n1,x2,y2,n2,xcnv,ycnv,ncnv)

!+++++++++++++++++++++++++++++++++
      IF ( ncnv<4 ) RETURN
!+++++++++++++++++++++++++++++++++

      CALL cubspl(xcnv,ycnv,ncnv,0,0)
      txcnvd = cubint(xcnv,ycnv,ncnv,t0,t)

      END FUNCTION txcnvd

REAL(KIND=8) FUNCTION txcnvl(t,ntypls,trp,ypls,nsx,swx,swy)
!----------------------------------------------------

!***  Called by: FOLD_AND_CONVOLVE, TQSTRIEP
!***      Calls: CUBINT, CUBVAL

!  Computes the system dB/dt response by convolving the computed dI/dt with
!  the impulse B response of the earth.  For step current drops, system dB/dt
!  is computed asthe product of instantaneous current drop times the
!  earth step dB/dt.

!  This routine assumes that the source waveform is composed of NSX linear
!  segments.  Thus NSX-1 constant dI/dt values are contained in SWY(*,1).

!  The input earth response function (step dB/dt or equivalently, impulse B)
!  must be contained in a splined array of NTYPLS values of time (abscissa) TRP
!  and ordinate YPLS.  System dB/dt is computed by integrating YPLS between
!  the SWX points of constant dI/dt segments.

!              T - convolution time in sec measured from the beginning
!                  of the source waveform.
!      TRP, YPLS - abscissa & ordinate values of earth response function to
!                  be convolved.
!         NTYPLS - number of values in TRP and YPLS
!            SWX - abscissa of time values of source waveform in sec.
!       SWY(*,1) - dI/dt if it exists (0 otherwise)
!       SWY(*,2) - first difference values of source waveform
!                  (-delta I) in amps.
!            NSX - number of points in SWX & WAVEFORM

      IMPLICIT NONE
      REAL(KIND=8) , PARAMETER :: T0_min = 1.d-7
      INTEGER ntypls , nsx , jt
      REAL(KIND=8) t , tf , cnv , tb , delt , seg , trp(ntypls) , ypls(4,ntypls), &
    &                        swx(nsx) , swy(nsx,3) , tend , cubint , cubval
      LOGICAL der

      tf = t - trp(1)
      cnv = 0.d0
      DO jt = 2 , nsx
         IF ( swx(jt-1)>tf ) EXIT
         tb = t - min(tf,swx(jt))
         delt = swx(jt) - swx(jt-1)
         der = .FALSE.
         IF ( delt>T0_min ) THEN
            tend = t - swx(jt-1)
            der = .TRUE.
         ENDIF

!  For an instantaneous step drop in current, SEG is YPLS times SWY(*,2),
!  since YPLS is already the dB/dt step response.  Otherwise SEG is the
!  integral of YPLS * constant dI/dt SWY(*,1) since YPLS is also impulse B.

         IF ( der ) THEN
            seg = swy(jt,1)*cubint(trp,ypls,ntypls,tb,tend)
         ELSE
            seg = swy(jt,2)*cubval(trp,ypls,ntypls,tb)
         ENDIF
         cnv = cnv + seg
      ENDDO
      txcnvl = cnv

      END FUNCTION txcnvl

SUBROUTINE tqstrip(ider,ntypls,trp,ypls,nsx,swx,swy,ycnv)
!----------------------------------------------------------
!
!***  Called by: FOLD_AND_CONVOLVE
!***      Calls: CUBSPL, TXCNVD, TXCNVL

!  A stripped down version of TXCNVD is used to convolve the earth response
!  with the receiver waveform for for every point on that waveform.
!  The Geotem / Questem correlation is used to strip remnant primary field.
!  The result is sent back for binning into NCHNL receiver windows.
!
!   TRP, YPLS - abscissa & ordinate values of earth response function to
!               be convolved.
!        IDER - derivative indicator
!      NTYPLS - number of values in TRP and YPLS
!         SWX - abscissa of time values of source waveform in sec.
!         SWY - dI/dt values derived from receiver dB/dt. + raw waveform.
!         NSX - number of points in SWX & in each waveform stored in SWY
!        YCNV - the stripped convolved waveform
!
!  Defining  T1 = MIN {T, signal length,}, the convolution is formally
!  computed as
!
!   TXCNVD (T) = INT (T0 -> T) { YPLS (tau) * SWY (T-tau)  d tau }

!  where T0 = MAX { TRP(1), T - SWX (NSX)}
!

      IMPLICIT NONE
      REAL(KIND=8) , PARAMETER :: T0_min = 1.d-7
      INTEGER ider , ntypls , nsx , jt , mxcnv
      REAL(KIND=8) t , trp(ntypls) , ypls(4,ntypls) , swx(nsx) , swy(nsx,3) ,   &
     &     ycnv(4,nsx) , txcnvl , txcnvd , a1 , b1 , alpha

      INTENT (in)  ider , ntypls , trp , ypls , nsx , swx , swy
      INTENT (out) ycnv

      a1 = 0.d0
      b1 = 0.d0
      ycnv = 0.d0
      mxcnv = ntypls + nsx

      DO jt = 2 , nsx       !  Convolve NSW points using the derived waveform
         t = swx(jt)
         IF ( t<T0_min ) CYCLE
         IF ( ider==0 ) THEN
                           ! Waveform input as I or B field (derived dI/dt)
            ycnv(1,jt) = txcnvl(t,ntypls,trp,ypls,nsx,swx,swy)
         ELSE              ! Waveform input as voltage (known dI/dt)
            ycnv(1,jt) = txcnvd(mxcnv,t,ntypls,trp,ypls,nsx,swx,swy)
         ENDIF

         a1 = a1 + ycnv(1,jt)*swy(jt,3)
                                      !  Compute correlation
         b1 = b1 + swy(jt,3)*swy(jt,3)
      ENDDO

      alpha = a1/b1
      DO jt = 1 , nsx
         ycnv(1,jt) = ycnv(1,jt) - alpha*swy(jt,3)
      ENDDO

      CALL cubspl(swx,ycnv,nsx,0,0)

      END SUBROUTINE tqstrip
