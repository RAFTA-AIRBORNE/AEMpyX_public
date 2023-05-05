!*==forjac.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE forjac(ndata,xdata,xmodl,xwts,npar,xpar,cxpar,sumsq,   &
     &                  jcbn,a,verr,tdfd,cmp,step,ider,nsx,swx,swy,     &
     &                  ntyrp,trp,npuls,pulse,ntypls,nchnl,topn,tcls,   &
     &                  gstrp,astrp,nfrq,freq,knrm,norm,js,txcln,txa90, &
     &                  xrx,yrx,zrx,nrxst,nstat,sz,nlyr,rmu,reps,calf,  &
     &                  ctau,cfreq)
!--------------------------------------------------------------------------------------
 
!  Sets up and calls for forward model computation.
!  It also calculates the Jacobian and error vector if required
!  New convention: Dec, 2003: VERR is now VD - VM
!  Thus DELPAR is now added rather than subtracted during updates.
 
!*** Called by: NLSQ2
!***     Calls: SET_CELLS, SHIFT_CELLS, GET_FWD_MODL, CNVRT2_MPAR
 
!             General Inversion Input Variables
!             ---------------------------------
!
!          JS - station index
!       NDATA - dimension of vector to be inverted: = NCMP * NCHNL or 2*NFRQ.
!       XDATA - data to be inverted in user-specified units
!       XMODL - model data in user-specified units
!        XWTS - weights for XDATA  (0 or 1)
!        NPAR - number of parameters to be inverted (nominally 2*NLYR - 1)
!        XPAR - array of transformed model parameters
!       CXPAR = 0 => parameter is completely free to vary as dictated by inversion step
!             = 1 => parameter is fixed
!             = 2 => parameter is constrained by elasticity.
!             = 3 => parameter bounds are buffered.
!       SUMSQ - sum squared scaled symmetric error
!        JCBN - true if Jacobian required
!           A - a large array which carries the Jacobian out of this routine
!     VERR(J) - scaled symmetric error in channel J
!             AEM System Input Variables
!             --------------------------
!
!        TDFD = 1 for TD; 2 for FD
!         CMP - (time-domain only) component to be inverted
!               11: in-line; 13: vertical; 2: joint vertical & in-line
!               3 component; 4 total field
!
!        KNRM - dimension of NORM: = 3 for time-domain,  = NFRQ for frequency-domain
!        NORM - PPM conversion
!        FREQ - array of NFRQ frequencies
!        STEP = 1 iff step response is to be computed
!        IDER = 1 if source waveform was dB/dt; = 0 if amps pr B
!         NSX - number of points used to discretise transmitter signal
!         SWX - abscissae (seconds) of current waveform
!         SWY - dI/dt * Tx moment & nanotesla conversion at times SWX
!       NTYRP - number of values in TRP for total signal length: 2 * NPULS *PULSE
!         TRP - array of time values for FD -> TD transformations
!      NTYPLS - number of TRP values in 1 PULSE
!       NCHNL = number of time domain channels
!        TOPN - time at which receiver channel I opens.
!        TCLS - time at which receiver channel I closes.
!       GSTRP = 1 => apply Questem-Geotem stripping algorithm
!       ASTRP = 1 => apply Aerotem stripping algorithm
!        FREQ = array of NFRQ frequencies
!       TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!       TXA90 - true for vertical co-planar briadside array
!       NRXST - dimension for receiver offset & transmitter tilt
!             = NFRQ in FD;  = NSTAT in TD
!         ZRX - vertical receiver offset for each frequency; below = positive
!         XRX - in-line receiver offset for each frequency;  behind = positive
!         YRX - transverse receiver offset for each frequency; Port = positive.
!
!             AEM Survey Input Variables
!             --------------------------
!
!       NSTAT - number of stations in survey line.
!          SZ - altitude (re gnd level) of transmitter
!
!             Model Description Input Variables
!             ---------------------------------
!
!        NLYR - number of layers (1 or 2)
!         THK - layer thicknesses
!         RMU - mu(i) / mu(0)
!        REPS - relative dielectric constant
!        CALF, CTAU & CFREQ are the layered earth Cole-Cole parameters.
!
!
!          PHYSICAL PARAMETER 0RDERING
!          ---------------------------
!
!          1 : NLYR      - resistivities
!          NLYR+1 : NPAR - thicknesses
!
!
!          TRANSFORMED PARAMETERS FOR INVERSION  - Logarithms
!          ------------------------------------
 
 
      IMPLICIT NONE
      INTEGER ndata , npar , nlyr , cxpar(npar) , xwts(ndata) , nchnl , &
     &        tdfd , cmp , knrm , nfrq , nstat , nrxst , step , ider ,  &
     &        nsx , npuls , ntypls , ntyrp , gstrp , astrp , js , jd ,  &
     &        jp , jp1 , jf , jt , nc1 , nc2 , nc3
      REAL sumsq , a(ndata,npar+1) , xpar(npar) , sz(nstat) , norm(knrm)&
     &     , freq(nfrq) , swx(nsx) , swy(nsx,3) , pulse , topn(nchnl) , &
     &     tcls(nchnl) , trp(ntyrp) , x2 , x3 , vm , vj , vd , denom ,  &
     &     delta , xp0 , parfac
      REAL , DIMENSION(ndata) :: xmodl , xmodl0 , xdata , verr
      REAL , DIMENSION(nrxst) :: txcln , xrx , yrx , zrx
      REAL , DIMENSION(nlyr) :: res , reps , rmu , calf , ctau , cfreq
      REAL , DIMENSION(nlyr-1) :: thk
      REAL , DIMENSION(nchnl,nstat,3) :: btd
      COMPLEX xbfd
      COMPLEX , DIMENSION(nfrq,nstat,3) :: bfd
      LOGICAL jcbn , txa90
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynml1/ ndata,xdata,xmodl,xwts,npar,xpar,cxpar,sumsq,    &
     &                  jcbn,a,verr,tdfd,cmp,step,ider,nsx,swx,swy,     &
     &                  ntyrp,trp,npuls,pulse,ntypls,nchnl,topn,tcls,   &
     &                  gstrp,astrp,nfrq,freq,knrm,norm,js,txcln,txa90, &
     &                  xrx,yrx,zrx,nrxst,nstat,sz,nlyr,rmu,reps,calf,  &
     &                  ctau,cfreq
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynml2/ js,step,ider,nsx,swx,swy,npuls,pulse,ntypls,    &
     &                  ntyrp,trp,nchnl,topn,tcls,txcln,nstat,sz,zrx,   &
     &                  xrx,yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,  &
     &                  gstrp,astrp,btd

! Compute initial forward model & compute error
 
      CALL cnvrt2_mpar(npar,nlyr,xpar,res,thk)
 
      nc1 = nchnl
      nc2 = 2*nchnl
      nc3 = 3*nchnl
      jp1 = 0
      parfac = 1.0
      CALL get_fwd_modl(jp1,parfac)
      xmodl0 = xmodl 
  
      sumsq = 0.
      DO jd = 1 , ndata
         verr(jd) = 0.
         IF ( xwts(jd)>0 ) THEN
            vd = xdata(jd)
            vm = xmodl0(jd)
            denom = sqrt((vm*vm+vd*vd)/2.0)
            verr(jd) = xwts(jd)*(vd-vm)/denom
            sumsq = sumsq + verr(jd)**2
         ENDIF
      ENDDO
 
!  Initialise and then compute the Jacobian as the derivative of log(volts) wrt
!  log(parameter) for a three percent step.  Skip over held parameters
 
      IF ( jcbn ) THEN
         a = 0.
         parfac = 0.03
 
         DO jp1 = 1 , nlyr          ! Vary resistivities.
            IF ( cxpar(jp1)/=1 ) THEN
               xp0 = res(jp1)
               res(jp1) = (1.+parfac)*res(jp1)
               CALL get_fwd_modl(jp1,parfac)
               res(jp1) = xp0
            ENDIF
         ENDDO
         DO jp = 1 , nlyr - 1
            jp1 = nlyr + jp           ! Vary layer dimesnsion.
            IF ( cxpar(jp1)/=1 ) THEN
               xp0 = thk(jp)
               thk(jp) = (1.+parfac)*thk(jp)
               CALL get_fwd_modl(jp1,parfac)
               thk(jp) = xp0
            ENDIF
         ENDDO
      ENDIF
 
      CONTAINS
 
      SUBROUTINE get_fwd_modl(jp1,parfac)
!  ------------------------------------
 
!***  Called by: FORJAC
!***      Calls: TEM_3D, HSBOSS_TD, HSBOSS_FD
 
!  If JP1 = 0, performs forward model computation using existing parameters
!  If JP1 > 0, performs forward model computation where parameter JP1 is multiplied by
!             PARFAC and Jacobian column JP1 isconstructed.
 
      INTEGER jp1
      REAL parfac , cstx(nrxst) , sntx(nrxst)
 
      IF ( tdfd==1 ) THEN
     WRITE(*,*) 'THIS IS mynml in forjac td'
     WRITE(*,nml=mynml2)

    
         CALL hsboss_td(js,step,ider,nsx,swx,swy,npuls,pulse,ntypls,    &
     &                  ntyrp,trp,nchnl,topn,tcls,txcln,nstat,sz,zrx,   &
     &                  xrx,yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,  &
     &                  gstrp,astrp,btd)
 
         write (*,*) cmp,  norm 
         IF ( cmp==13 .OR. cmp==2 .OR. cmp==3 ) xmodl(1:nc1) = norm(3)  &
     &        *btd(1:nchnl,js,3)
 
         IF ( cmp==2 .OR. cmp==3 ) xmodl(nc1+1:nc2) = norm(1)           &
     &        *btd(1:nchnl,js,1)
 
         IF ( cmp==3 ) xmodl(nc2+1:nc3) = norm(2)*btd(1:nchnl,js,2)
 
         IF ( cmp==11 ) xmodl(1:nc1) = norm(1)*btd(1:nchnl,js,1)
 
         IF ( cmp==4 .OR. cmp==42 ) THEN
            DO jt = 1 , nchnl
               x2 = btd(jt,js,1)**2 + btd(jt,js,3)**2
               x3 = x2 + btd(jt,js,2)**2
               IF ( cmp==42 ) THEN
                  xmodl(jt) = sqrt(x2)
               ELSE
                  xmodl(jt) = sqrt(x3)
               ENDIF
            ENDDO
         ENDIF
 
      ELSE
     WRITE(*,*) 'THIS IS mynml in forjac fd'
     WRITE(*,nml=mynml1)

         CALL hsboss_fd(js,nfrq,freq,txcln,txa90,nstat,sz,zrx,xrx,yrx,  &
     &                  nlyr,res,reps,rmu,thk,calf,ctau,cfreq,bfd)
 
!  Store maximally coupled components in XMODL
 
         cstx(1:nfrq) = cos(txcln(1:nfrq))
         sntx(1:nfrq) = sin(txcln(1:nfrq))
 
         DO jf = 1 , nfrq
            IF ( txa90 ) THEN
               xbfd = norm(jf)*bfd(jf,js,2)
            ELSE
               xbfd = norm(jf)                                          &
     &                *(bfd(jf,js,1)*sntx(jf)+bfd(jf,js,3)*cstx(jf))
            ENDIF
            xmodl(jf) = real(xbfd)
            xmodl(jf+nfrq) = aimag(xbfd)
         ENDDO
      ENDIF
 
      IF ( jp1>0 ) THEN
         DO jd = 1 , ndata
            vm = xmodl0(jd)
            vj = xmodl(jd)
            delta = xwts(jd)*(vj-vm)
            denom = sqrt((vm**2+vj**2)/2.)
            IF ( abs(delta)>1.E-8*denom ) a(jd,jp1)                     &
     &           = delta/(parfac*denom)
         ENDDO
      ENDIF
 
      END SUBROUTINE get_fwd_modl
      END SUBROUTINE forjac
