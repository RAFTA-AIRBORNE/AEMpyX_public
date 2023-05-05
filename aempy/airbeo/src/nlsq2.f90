      SUBROUTINE nlsq2(js,nw,nw1,nlg,mv1prt,outprt,maxits,cnvrg,pctcnv, &
     &                 ndata,xdata,xmodl,xwts,npar,cxpar,elas,lbnd,ubnd,&
     &                 tdfd,cmp,knrm,norm,step,ider,nsx,swx,swy,ntyrp,  &
     &                 trp,npuls,pulse,ntypls,nchnl,topn,tcls,gstrp,    &
     &                 astrp,nfrq,freq,txcln,txa90,nrx,nrxst,xrx,yrx,   &
     &                 zrx,nstat,sx,sy,sz,rx,ry,rz,gnd_lvl,title,line,  &
     &                 nlyr,thk,res,rmu,reps,calf,ctau,cfreq,mpar,      &
     &                 rmserr)
!---------------------------------------------------------------------------------------
 
!***  Called by: MAIN
!***      Calls: CNVRT_BOUNDS, CNVRT2_XPAR, ESVD, FORJAC, INDEX_MPAR,
!                NOISE_2_SIGR, PARAMETER_SENSITIVITY, SOLVE2, WRITE_MODEL,
!                WRITE_MDATA, WRITE_MISFIT
 
!  SVD non-linear least square inversion using a modified
!  Gauss-Newton-Marquardt method as described in Jupp & Vozoff, 1975,
!  Stable iterative methods for the inversion of geophysical data,
!  Geophys. J. R. Astr. Soc.,42
 
!  New convergence criteria based on pre-specified RMS threshold
!  added December, 1998
!  New convention: Dec, 2003: VERR is now VD - VM
!  Thus DELPAR is now added rather than subtracted during updates.
!=========================================================================
 
!***  Called by: MAIN
!***      Calls: ESVD, FORJAC, WRITE_MODEL, WRITE_MISFIT, SOLVE2
!
!
!   Output: MPAR:  NLYR resistivities followed by NLYR-1 layer thicknesses
!
!
!
!             General Inversion Input Variables
!             ---------------------------------
!
!              JS - station index
!              NW - verbose output unit number
!             NW1 = plot file unit number
!             NLG = output unit for Airbeo.log
!          MV1PRT & OUTPRT are print options for the MV1 & OUT files respectively.
!                 =  0  No output DURING inversion.  The final model set AFTER inversion,
!                       but NOT the final model data, is written to output files.
!
!                 =  1  as above plue plus final model data
!
!                 =  2  as above plus intermediate model sets after each iteration
!
!                 =  3 as above plus intermediate model data after each iteration
!
!          MAXITS - maximum permitted iterations.
!           CNVRG = 1 => converge on predicted decrease etc
!                 = 2=> stop when RMS error < PCTCNV
!           NDATA = NCMP NCHNL or 2 * NFRQ
!           NDATA = total number of readings per station to be inverted(TD or FD)
!                 = the number of rows of the Jacobian matrix
!                 = NCHNL for time-domain when CMP = 11, 13, 4, 42, 43
!                 = 2* NCHNL for time-domain when CMP = 2
!                 = 3* NCHNL for time-domain when CMP = 3
!
!                 = 2 * NFRQ for frequency-domain
!           XDATA - Array of survey data points to be inverted.
!           XMODL - Array of model data points.
!            XWTS - Array of weights for XDATA (integer 0 or 1)
!           NPAR = 2*NLYR -1 = number of parameters
!           XPAR - On call, XPAR contains the log (initial parameter estimates).
!                  On exit it contains the final values.
!          CXPAR = 0 => parameter is completely free to vary as dictated by inversion step
!                = 1 => parameter is fixed
!                = 2 => parameter is constrained by elasticity.
!                = 3 => parameter bounds are buffered.
!           ELAS - [0 to 1] fraction of proposed step actually taken.
!     LBND, UBND - lower & upper bound of parameter
!
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
!         NRX - number of receivers = NFRQ in FD; = 1 in TD
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
!       SX, SY, SZ: north, east and altitude (re gnd level) of transmitter
!       RX, RY, RZ: north, east and altitude (re gnd level) of receiver(s)
!
!     GND_LVL - ground level
!        LINE - line number
!
!             Model Description Input Variables
!             ---------------------------------
!
!        NLYR - number of layers (1 or 2)
!         THK - layer thicknesses
!         RES - array of layer resistivities
!         RMU - mu(i) / mu(0)
!        REPS - relative dielectric constant
!        CALF, CTAU & CFREQ are the layered earth Cole-Cole parameters.
!
!              Other Stuff
!              -----------
!
!              A -  Jacobian matrix
!            WSP - working space.
!           VERR - error vector
!         RMSERR - RMS error (in percent)
!            BND - estimate of noise to signal ratio, used for damping limit.
!            WSP - working space.
!     SV, UMAT, and VMAT are the decomposition of A.  SV contains eigenvalues.
 
!      The model set consists of the model parameters, parameter importance,
!      standard error and RSVT, the relative singular value threshold.
!
!      The model data (distinct from survey data) is defined as the fields or
!      ppm values for each channel or frequency for each station; ie the
!      forward model response for a given model.
!
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ip = 1 , Ql = selected_real_kind(12,80)
      REAL , PARAMETER :: Bnd = 0.001 , Expnd = 2. , Rsvt0 = 0.1 ,       &
     &                    Eta = 1.E-7 , Tol = .5E-31 ,                  &
     &                    Rad2deg = 180./3.141592654
      INTEGER js , nw , nw1 , nlg , mv1prt , outprt , its , cnvrg ,     &
     &        ndata , npar , cxpar(npar) , xwts(ndata) , tdfd , cmp ,   &
     &        step , ider , nsx , npuls , ntypls , ntyrp , nchnl ,      &
     &        gstrp , astrp , nfrq , nstat , knrm , nrxst , nrx , icnt ,&
     &        maxits , nsv , fits , jp , line(nstat) , mprnt , j1 ,     &
     &        mxerr
      REAL pctcnv , a(ndata,npar+1) , umat(ndata,npar) , vmat(npar,npar)&
     &     , wsp(3*npar) , nsr , drms , rmserr , norm(knrm) , swx(nsx) ,&
     &     swy(nsx,3) , trp(ntyrp) , pulse , freq(nfrq) , gnd_lvl ,     &
     &     zero , delt , rsvt , fnm , gcrit , gtest , pdre , ssqnew ,   &
     &     sumsq , wgtsum , b1 , b2 , bmid
      REAL , ALLOCATABLE :: rms1(:)
      REAL , DIMENSION(nrxst) :: txcln , txdeg , xrx , yrx , zrx
      REAL , DIMENSION(nstat,nrx) :: rx , ry , rz
      REAL , DIMENSION(npar) :: sv , xpar , mpar , delpar , gxpar ,     &
     &                          import , elas , lbnd , ubnd , xlbnd ,   &
     &                          xubnd
      REAL , DIMENSION(nchnl) :: topn , tcls
      REAL , DIMENSION(ndata) :: verr , xmodl , xdata , verr_fill
      REAL , DIMENSION(nstat) :: sx , sy , sz , sz0 , rz0
      REAL(KIND=Ql) , DIMENSION(nstat) :: sxd , syd , rxd , ryd
      LOGICAL , PARAMETER :: Withu = .TRUE. , Withv = .TRUE.
      LOGICAL jcbn , txa90
      CHARACTER(LEN=120) cvar , title
      CHARACTER(LEN=80) ql0 , ql1
 
! Model Specific
      INTEGER nlyr
      REAL thk(nlyr-1)
      REAL , DIMENSION(nlyr) :: res , reps , rmu , chrg , calf , ctau , &
     &                          cfreq
 
!  Preset threshold parameters and index workspace, but return if problem
!  size in error.  First write initial model
 
      

      txdeg = Rad2deg*txcln
      sxd = real(sx,Ql)
      syd = real(sy,Ql)
      sz0 = sz + gnd_lvl
 
      rxd(js) = real(rx(js,1),Ql)
      ryd(js) = real(ry(js,1),Ql)
      rz0(js) = rz(js,1) + gnd_lvl
 
      verr_fill = 0.
      IF ( mv1prt==3 ) CALL write_xdata
 
      mprnt = 0
      chrg = 1. - calf
      CALL write_model(nw,mprnt,js,nlyr,thk,res,chrg,ctau,cfreq,rmu,    &
     &                 reps)
 
      zero = Bnd*Bnd
      IF ( zero<Eta ) zero = Eta
      gcrit = sqrt(Eta)
 
      a = 0. !  Initialise arrays
      verr = 0.
      import = 0.
 
      wgtsum = real(sum(xwts))
      its = 0
      rsvt = max(Bnd,Rsvt0)
                           ! Initialise eigenvalue damping at 10 percent
      WRITE (nw,'(/T3,A)') trim(title)
      WRITE (*,'(/T3,A)') trim(title)
      IF ( tdfd==2 ) THEN
         WRITE (nw,99011) js
         WRITE (*,99011) js
      ELSE
         IF ( cmp==11 ) cvar = 'in-line component'
         IF ( cmp==13 ) cvar = 'vertical component'
         IF ( cmp==2 ) cvar = 'joint vertical and in-line components'
         IF ( cmp==3 ) cvar =                                           &
     &               'joint vertical, in-line and transverse components'
         IF ( cmp==4 ) cvar = 'total field using all 3 components'
 
         WRITE (nw,99012) trim(cvar) , js
         WRITE (*,99012) trim(cvar) , js
      ENDIF
      WRITE (nw,99013) maxits
      WRITE (*,99013) maxits
 
!----------------------------------------------------------------------
!  Start of main loop.  Call FORJAC to get error and Jacobian.
!  First, tranform physical model into transformed parameters.
!  Call ESVD to find the S.V.D. of the Jacobian.
!----------------------------------------------------------------------
 
      CALL cnvrt2_xpar(npar,nlyr,res,thk,xpar)
      CALL cnvrt_bounds(npar,lbnd,ubnd,cxpar,xlbnd,xubnd)
 
      ALLOCATE (rms1(maxits))
      drms = 0.
      ITER_LOOP:DO its = 1 , maxits
 
         jcbn = .TRUE.
         CALL forjac(ndata,xdata,xmodl,xwts,npar,xpar,cxpar,sumsq,jcbn, &
     &               a,verr,tdfd,cmp,step,ider,nsx,swx,swy,ntyrp,trp,   &
     &               npuls,pulse,ntypls,nchnl,topn,tcls,gstrp,astrp,    &
     &               nfrq,freq,knrm,norm,js,txcln,txa90,xrx,yrx,zrx,    &
     &               nrxst,nstat,sz,nlyr,rmu,reps,calf,ctau,cfreq)
 
         rmserr = 100.*sqrt(sumsq/wgtsum)
         rms1(its) = rmserr
         IF ( its==1 ) THEN
            fits = 0
            mpar(1:npar) = exp(xpar(1:npar))
            IF ( mv1prt>1 ) THEN
               WRITE (nw1,99020) fits , js , rmserr , rsvt
               WRITE (nw1,99021) fits , mpar(1:npar)
            ENDIF
            IF ( outprt>1 ) WRITE (nw,99014) rmserr , rsvt
            IF ( mv1prt==3 ) CALL write_mdata(fits)
            IF ( outprt==3 ) CALL write_misfit(nw,fits,ndata,tdfd,nfrq, &
     &           freq,cmp,nchnl,verr,xmodl,xdata)
            WRITE (*,99014) rmserr , rsvt
         ENDIF
         fnm = 0.01*sqrt(sumsq)
         fnm = max(fnm,Eta)
 
!  Load the error vector into the NPAR+1 column of A.  On return from ESVD,
!  this column will contain the transformed error vector; i.e.,
!  VERR = U * VERR
 
         a(1:ndata,npar+1) = verr(1:ndata)

         CALL esvd(a,ndata,npar,Ip,Withu,Withv,sv,umat,vmat,wsp,Eta,Tol,&
     &             nw)
 
         IF ( abs(sv(1))<1.E-30 ) THEN
            CALL write_log_file(nlg,100,mxerr,2)
            WRITE (*,99022) js
            WRITE (nlg,99022) js
            RETURN
         ENDIF
 
         verr(1:ndata) = a(1:ndata,npar+1)
 
 
!  Solve for the correction vector, and test for convergence on
!  predicted decrease.  Loop over internal iterations.
 
         ICNT_LOOP:DO icnt = 0 , maxits
 
            CALL solve2(npar,rsvt,zero,vmat,verr,sv,nsv,delpar,wsp,pdre)
            delt = sqrt(pdre)
 
! If the predicted residual decrease < 1 percent of RMS error,
! terminate iterations.  Inversion won't improve.
 
            IF ( delt<fnm ) THEN
               WRITE (nw,99015) its
               WRITE (*,99015) its
               WRITE (nw,99016) rmserr , rsvt
               WRITE (*,99016) rmserr , rsvt
               GOTO 100
            ENDIF
 
            DO jp = 1 , npar
               SELECT CASE (cxpar(jp))
               CASE (0)                        ! No constraints
                  gxpar(jp) = xpar(jp) + delpar(jp)
               CASE (1)                        ! Fixed parameter
                  gxpar(jp) = xpar(jp)
               CASE (2)
                  gxpar(jp) = xpar(jp) + elas(jp)*delpar(jp)
               CASE (3)
                  b1 = xlbnd(jp)
                  b2 = xubnd(jp)
                  bmid = 0.5*(b1+b2)
                  IF ( xpar(jp)<b1 .OR. xpar(jp)>b2 ) xpar(jp) = bmid
                  b2 = xpar(jp) + elas(jp)*(b2-xpar(jp))
                  b1 = xpar(jp) + elas(jp)*(b1-xpar(jp))
                  gxpar(jp) = xpar(jp) + delpar(jp)
                  gxpar(jp) = min(gxpar(jp),b2)
                  gxpar(jp) = max(gxpar(jp),b1)
               END SELECT
            ENDDO
 
!  Get the error for model with corrected parameters.
!  Test for improvement (decrease) in residual with Goldstein condition
!  on ratio of actual to computed decrease in error.  If it fails, reduce
!  step and try again.  Give up and return after MAXITS. "internal" iterations.
 
            jcbn = .FALSE.
 
            CALL forjac(ndata,xdata,xmodl,xwts,npar,gxpar,cxpar,ssqnew, &
     &                  jcbn,a,verr,tdfd,cmp,step,ider,nsx,swx,swy,     &
     &                  ntyrp,trp,npuls,pulse,ntypls,nchnl,topn,tcls,   &
     &                  gstrp,astrp,nfrq,freq,knrm,norm,js,txcln,txa90, &
     &                  xrx,yrx,zrx,nrxst,nstat,sz,nlyr,rmu,reps,calf,  &
     &                  ctau,cfreq)
 
            gtest = sumsq - ssqnew
 
            IF ( gtest>gcrit*pdre ) THEN
                                    !  Error reduced using smaller step
               IF ( icnt==0 ) THEN
                  rsvt = rsvt/Expnd !  Decrease eigenvalue threshold damping
                  rsvt = max(Bnd,rsvt)
               ENDIF
               EXIT ICNT_LOOP       !  Start next iteration
            ENDIF
            rsvt = rsvt*Expnd       !  No error decrease. Raise threshold
            IF ( icnt==maxits ) THEN
               WRITE (nw,99001) icnt
                                    !  No improvement possible.  Maybe another starting guess?
99001          FORMAT (/T3,'The solution is trapped.  ICNT = ',I2/T3,   &
     &               'Another starting guess may yield a better result.'&
     &               )
               GOTO 100
            ENDIF
         ENDDO ICNT_LOOP
 
 
!  Error reduced.  Accept step, write out summary, and test convergence.
 
         xpar(1:npar) = gxpar(1:npar)
         delt = sqrt(gtest)
         sumsq = ssqnew
         rmserr = 100.*sqrt(sumsq/wgtsum)
         fnm = 0.01*sqrt(sumsq)
         fnm = max(fnm,Eta)
 
!  If the predicted residual decrease < 1 percent of RMS error,
!  If the error decrease from the new iteration < 1 percent of the previous error,
!  claim convergence and terminate iterations because inversion won't improve.
 
!  Else, write out the current model and continue iterating up until ITS = MAXITS.
 
         IF ( delt<fnm ) THEN
            WRITE (nw,99015) its
            WRITE (*,99015) its
            WRITE (nw,99016) rmserr , rsvt
            WRITE (*,99016) rmserr , rsvt
            EXIT ITER_LOOP
         ENDIF
         IF ( cnvrg==2 .AND. rmserr<pctcnv ) THEN
            WRITE (nw,99017) pctcnv , its
            WRITE (*,99017) pctcnv , its
            WRITE (nw,99016) rmserr , rsvt
            WRITE (*,99016) rmserr , rsvt
            EXIT ITER_LOOP
         ENDIF
         IF ( its==maxits ) THEN
            WRITE (nw,99018) its
            WRITE (*,99018) its
            WRITE (nw,99016) rmserr , rsvt
            WRITE (*,99016) rmserr , rsvt
            EXIT ITER_LOOP
         ENDIF
 
         mpar(1:npar) = exp(xpar(1:npar))
         IF ( mv1prt>1 ) THEN
            WRITE (nw1,99020) its , js , rmserr , rsvt
            WRITE (nw1,99021) its , mpar(1:npar)
         ENDIF
         IF ( mv1prt>2 ) CALL write_mdata(its)
 
         WRITE (*,99002) its , rmserr , rsvt
99002    FORMAT (/I4,' Iterations completed.  Symmetric RMS error =',   &
     &           F8.2,' percent.'/T44,'RSVT =',F8.3)
         IF ( outprt>1 ) THEN
            res(1:nlyr) = mpar(1:nlyr)
            thk(1:nlyr-1) = mpar(nlyr+1:npar)
            CALL write_model(nw,its,js,nlyr,thk,res,chrg,ctau,cfreq,rmu,&
     &                       reps)
            WRITE (nw,'(2X)')
            WRITE (*,'(2X)')
            WRITE (nw,99016) rmserr , rsvt
            WRITE (*,99016) rmserr , rsvt
 
            IF ( outprt==3 .AND. its<maxits )                           &
     &           CALL write_misfit(nw,its,ndata,tdfd,nfrq,freq,cmp,     &
     &           nchnl,verr,xmodl,xdata)
         ENDIF
         rms1(its) = rmserr
         IF ( its>3 ) drms = rms1(its-2) - rms1(its)
         IF ( its>10 .AND. drms<1. ) THEN
            WRITE (nw,99019) its
            WRITE (*,99019) its
            EXIT ITER_LOOP
         ENDIF
      ENDDO ITER_LOOP
                    !  END OF MAIN LOOP.  Write final model and exit.
 100  DEALLOCATE (rms1)
 
      CALL parameter_sensitivity(npar,vmat,sv,Bnd,import)
      CALL noise_2_sigr(npar,ndata,xmodl,xdata,xwts,nsr)
 
      WRITE (nw,99003)
99003 FORMAT (/T3,50('='))
      WRITE (nw1,99004) its , js , rmserr , rsvt
99004 FORMAT (T1,'/'/T1,'/ FINAL_ITERATION  ',I2.2,4X,'Station',I4/T1,  &
     &        '/ PERCENT_RMS_ERROR:',F9.2,3X,'RSVT:',F8.3)
      mpar(1:npar) = exp(xpar(1:npar))
      WRITE (nw1,99005) mpar(1:npar)
99005 FORMAT ('/ FINAL_MODEL',84G13.4)
      WRITE (nw1,99006) import(1:npar)
99006 FORMAT (T1,'/ IMPORTANCE:  ',2X,84F6.2)
      IF ( mv1prt==1 .OR. mv1prt==2 ) CALL write_xdata
      IF ( mv1prt>0 ) CALL write_mdata(-1)
 
      mprnt = -its
      res(1:nlyr) = mpar(1:nlyr)
      thk(1:nlyr-1) = mpar(nlyr+1:npar)
      CALL write_model(nw,mprnt,js,nlyr,thk,res,chrg,ctau,cfreq,rmu,    &
     &                 reps)
      WRITE (nw,99007)
99007 FORMAT (/T3,'Parameter Importance'/T3,'--------------------')
      WRITE (nw,99008,ADVANCE='NO') 'RES_01' , ('    RES_',j1,j1=2,nlyr)
99008 FORMAT (T3,A,30(A,I2.2))
      WRITE (nw,99009) ('  THICK_',j1,j1=1,nlyr-1)
99009 FORMAT (1X,30(A,I2.2))
      WRITE (nw,'(F6.2,30F10.2)') import(1) , (import(j1),j1=2,2*nlyr-1)
 
      WRITE (nw,99010) rmserr , nsr
99010 FORMAT (/T23,'Symmetric RMS error =',F9.2,' percent.'/T21,        &
     &        'Noise to signal ratio =',F9.3)
 
      fits = -1
      sxd = real(sx,Ql)
      syd = real(sy,Ql)
      sz0 = sz + gnd_lvl
      IF ( outprt>0 ) CALL write_misfit(nw,fits,ndata,tdfd,nfrq,freq,   &
     &                                  cmp,nchnl,verr,xmodl,xdata)
 
99011 FORMAT (/T3,'Begin frequency-domain inversion for station',I4)
99012 FORMAT (/T3,'Begin time-domain inversion on ',A,' for station',I4)
99013 FORMAT (T3,'Maximum iterations =',I3)
99014 FORMAT (/T3,'Initial symmetric root mean square error =',F8.2,    &
     &        ' percent.'/T3,                                           &
     &        'Initial RSVT (Relative Singular Value Threshold) =',F8.3)
99015 FORMAT (//T3,'Convergence on predicted decrease after',I3,        &
     &        ' iterations.')
99016 FORMAT (T3,'Symmetric RMS error =',F8.2,' percent.',3X,'RSVT =',  &
     &        F7.3)
99017 FORMAT (/T3,'Convergence within RMS error threshold of',F7.2,     &
     &        ' after',I3,' iterations.')
99018 FORMAT (/T3,'Inversion finished after the maximum',I3,            &
     &        ' Iterations')
99019 FORMAT (/T3,                                                      &
     &'-----------------------------------------------------------------&
     &---'/T3,'Inversion terminated after',I3,                          &
     &' iterations because, the reduction in'/T3,                       &
     &'the RMS error from the last two iterations was less than 1 percen&
     &t.'/T3,                                                           &
     &'_________________________________________________________________&
     &___')
99020 FORMAT (T1,'/'/T1,'/ ITERATION  ',I2.2,4X,'Station',I4/T1,        &
     &        '/ PERCENT_RMS_ERROR:',F9.2,3X,'RSVT:',F8.3)
99021 FORMAT (T1,'/ MODEL_',I2.2,84G13.4)
99022 FORMAT (/T3,'INVERSION HALTED for station',I4,                    &
     &        ' due to singular value decomposition failure.')
 
      CONTAINS
 
      SUBROUTINE write_mdata(kts)
!  --------------------------
 
!***  Called by NLSQ2
 
      INTEGER kts
 
      IF ( kts==-1 ) THEN ! Write Final Model Data
         WRITE (ql0,99001) line(js)
 
99001    FORMAT (T2,I10,'_ZFNL')
         READ (ql0,'(A)') ql1
         WRITE (nw1,99007)
         WRITE (nw1,99006) trim(adjustl(ql1)) , js
      ELSEIF ( kts==-2 ) THEN  ! Write Failure
         WRITE (ql0,99002) line(js)
99002    FORMAT (T2,I10,'_Failed_Inversion')
         READ (ql0,'(A)') ql1
         WRITE (nw1,99007)
         WRITE (nw1,99006) trim(adjustl(ql1)) , js
      ELSE                     ! Write Model Data after iteration KTS
         WRITE (ql0,99003) line(js) , kts
99003    FORMAT (T2,I10,'_I',I2.2)
         READ (ql0,'(A)') ql1
         WRITE (nw1,99006) trim(adjustl(ql1)) , js
      ENDIF
 
      IF ( tdfd==1 ) THEN
         WRITE (nw1,99004) js , syd(js) , sxd(js) , sz0(js) , txdeg(js) &
     &                     , ryd(js) , rxd(js) , rz0(js) ,              &
     &                     xmodl(1:ndata) , 100.*verr(1:ndata)
99004    FORMAT (T1,I5,2F12.1,F8.1,F6.1,2F12.1,F8.1,350G13.4:)
      ELSE
         WRITE (nw1,99005) js , syd(js) , sxd(js) , sz0(js) , ryd(js) , &
     &                     rxd(js) , rz0(js) , xmodl(1:ndata) ,         &
     &                     100.*verr(1:ndata)
99005    FORMAT (T1,I5,2F12.1,F8.1,2F12.1,F8.1,300G13.4)
      ENDIF
99006 FORMAT (T2,'Line ',A,4X,'Station',I4)
99007 FORMAT (T1,'/')
 
      END SUBROUTINE write_mdata
 
      SUBROUTINE write_xdata
!  ----------------------
 
!***  Called by NLSQ2
 
! Writes Survey Data
 
      WRITE (ql0,99001) line(js)
 
99001 FORMAT (T2,I10,'_Survey_Data')
      READ (ql0,'(A)') ql1
      WRITE (nw1,99002)
99002 FORMAT (T1,'/')
      WRITE (nw1,99003) trim(adjustl(ql1)) , js
99003 FORMAT (T2,'Line ',A,4X,'Station',I4)
 
      IF ( tdfd==1 ) THEN
         WRITE (nw1,99004) js , syd(js) , sxd(js) , sz0(js) , txdeg(js) &
     &                     , ryd(js) , rxd(js) , rz0(js) ,              &
     &                     xdata(1:ndata) , verr_fill(1:ndata)
99004    FORMAT (T1,I5,2F12.1,F8.1,F6.1,2F12.1,F8.1,350G13.4:)
      ELSE
         WRITE (nw1,99005) js , syd(js) , sxd(js) , sz0(js) , ryd(js) , &
     &                     rxd(js) , rz0(js) , xdata(1:ndata) ,         &
     &                     verr_fill(1:ndata)
99005    FORMAT (T1,I5,2F12.1,F8.1,2F12.1,F8.1,300G13.4)
      ENDIF
 
      END SUBROUTINE write_xdata
 
      END SUBROUTINE nlsq2
