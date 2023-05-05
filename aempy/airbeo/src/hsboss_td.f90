!*==hsboss_td.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE hsboss_td(js,step,ider,nsx,swx,swy,npuls,pulse,ntypls, &
     &                     ntyrp,trp,nchnl,topn,tcls,txcln,nstat,sz,zrx,&
     &                     xrx,yrx,nlyr,res,reps,rmu,thk,calf,ctau,     &
     &                     cfreq,gstrp,astrp,btd)
!------------------------------------------------------------------------------------
 
!***  Called by: MAIN, GET_FWD_MODL
!***      Calls: COSTRN, CUBSPL HSMD_FD, FOLD_AND_CONVOLVE
 
!  Computes BTD, the time-domain layered earth response convolved with the
!  excitation waveform and the receiver channels per unit receiver area.
!  For impulse response, it computes dB/dt in nT / s which is the same as
!  nanovolts per unit area.
 
!  For step response, it computes B in nanoteslas.
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  SIGN CONVENTION:
!  ----------------
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 
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
!                             OUTPUT
!                             ------
!     BTD(JT,JS,1) - the in-line component of the layered earth response at
!                    time JT, station JS.
!     BTD(JT,JS,2) - the horizontal transverse component
!     BTD(JT,JS,3) - the vertical component
 
      USE frequency_select
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Nfrq = Nf_6pde , Nrxf = 1 ,                &
     &                       Ql = selected_real_kind(12,80)
      REAL , PARAMETER :: Twopi = 6.283185307
      INTEGER step , ider , nsx , npuls , ntypls , ntyrp , nchnl ,      &
     &        nstat , nlyr , tdfd , gstrp , astrp , js , jf , jt , jc
      REAL swx(nsx) , swy(nsx,3) , pulse , trp(ntyrp) , sz(nstat) ,     &
     &     alt , t , yprm(4,ntyrp) , ycum(nchnl) , yfrq(4,Nfrq) ,       &
     &     freq(Nfrq) , wf(Nfrq) , costrn , btd(nchnl,nstat,3)
      REAL , DIMENSION(nchnl) :: topn , tcls
      REAL , DIMENSION(nlyr) :: res , reps , ctau , cfreq , calf , rmu
      REAL , DIMENSION(nlyr-1) :: thk
      REAL , DIMENSION(nstat) :: txcln , zrx , xrx , yrx
      REAL(KIND=Ql) , DIMENSION(Nrxf) :: xrxd , yrxd , zrxd , txclnd
      COMPLEX(KIND=Ql) bfdd(Nfrq,3)
      COMPLEX bfd
      LOGICAL txa90
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmltd/ js,step,ider,nsx,swx,swy,npuls,pulse,ntypls, &
     &                     ntyrp,trp,nchnl,topn,tcls,txcln,nstat,sz,zrx,&
     &                     xrx,yrx,nlyr,res,reps,rmu,thk,calf,ctau,     &
     &                     cfreq,gstrp,astrp,btd
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmlfd/freq,alt,txclnd,txa90,zrxd,xrxd,yrxd,nlyr, &
     &             thk,res,reps,rmu,calf,ctau,cfreq,tdfd,bfdd
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      tdfd = 1
      txa90 = .FALSE.
 
      freq(1:Nfrq) = frq_6pde(1:Nfrq)
      wf(1:Nfrq) = log(Twopi*freq(1:Nfrq))
 
      txclnd(1) = real(txcln(js),kind=Ql)
      alt = sz(js)
      zrxd(1) = real(zrx(js),kind=Ql)
      xrxd(1) = real(xrx(js),kind=Ql)
      yrxd(1) = real(yrx(js),kind=Ql)
 
!      WRITE(*,*) 'THIS IS mynmlfd in boss_td'
!      WRITE(*,nml=mynmlfd)
    
      CALL hsmd_fd(Nfrq,freq,alt,Nrxf,txclnd,txa90,zrxd,xrxd,yrxd,nlyr, &
     &             thk,res,reps,rmu,calf,ctau,cfreq,tdfd,bfdd)
 
!    WRITE(*,*) 'THIS IS mynmltg in boss_td'
!     WRITE(*,nml=mynmltd)

!  Compute BTD, the 'observed' layered earth response by folding the BLEXT,
!  the extended response over NPULS bipolar cycles into 1 PULSE and then
!  convolving this with the TX waveform.  It is during the convolution that we
!  shift from teslas to nanoteslas or nT/s.
 
      yfrq = 0.
      DO jc = 1 , 3
         DO jf = 1 , Nfrq
            bfd = cmplx(bfdd(jf,jc))
            yfrq(1,jf) = aimag(bfd)/(Twopi*freq(jf))
         ENDDO
         CALL cubspl(wf,yfrq,Nfrq,0,0)
!         WRITE(*,*) ' wf ',jc
!         WRITE(*,*)  wf
 
         yprm = 0.
         DO jt = 1 , ntyrp
                      !  Convert to step-function time-domain.
            t = trp(jt)
            yprm(1,jt) = costrn(wf,yfrq,Nfrq,Nfrq,t)
         ENDDO
!         WRITE(*,*) ' yprm: '
!         WRITE(*,*)  yprm

         CALL fold_and_convolve(step,ider,nsx,swx,swy,npuls,pulse,trp,  &
     &                          ntypls,ntyrp,nchnl,topn,tcls,yprm,gstrp,&
     &                          astrp,ycum)
 
         btd(1:nchnl,js,jc) = ycum(1:nchnl)
         !WRITE(*,'(/i8/11G12.4/)') jc,ycum(1:nchnl)
      ENDDO
 
      END SUBROUTINE hsboss_td
