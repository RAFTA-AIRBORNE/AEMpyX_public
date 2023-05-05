!*==main.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      PROGRAM airBeo
!------------
 
!*** Calls DCPRM_FD, DCPRM_TD, HSBOSS_TD, HSBOSS_FD, NLSQ2,
!          SET_NORM_TD, SET_SOURCE, WRITE_FD, WRITE_TD, WRITE_MODEL,
!          WRITE_LOG_FILE
 
!*** Calls from INPUT_DATA_FOR_AIRBEO:
!          READ_SYSTEM_AND_SURVEY, READ_MODEL, READ_INVRT_CNTRL_AND_DATA
!          SET_SURVEY, SET_TRP, WRITE_NW1_INITIAL
 
      USE input_data_for_airBeo
 
      IMPLICIT NONE
      
      INTEGER ider , knrm , mprnt
      REAL rmserr

      REAL , ALLOCATABLE , DIMENSION(:) :: norm
      REAL , ALLOCATABLE , DIMENSION(:,:,:) :: btd
      COMPLEX , ALLOCATABLE , DIMENSION(:,:,:) :: bfd
      REAL cmp_start , cmp_end , elapsed
      LOGICAL wrt_nw1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmlfd/ nfrq,freq,txcln,txa90,nstat,sz,zrx,xrx,         &
     &                 yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,       &
     &                 bfd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmltd/ js,step,ider,nsx,swx,swy,npuls,pulse,           &
     &                        ntypls,ntyrp,trp,nchnl,topn,tcls,txcln,   &
     &                        nstat,sz,zrx,xrx,yrx,nlyr,res,reps,rmu,   &
     &                        thk,calf,ctau,cfreq,gstrp,astrp,btd

      NAMELIST /mynmlinv/  js,nw,nw1,nlg,mv1prt,outprt,maxits,cnvrg,pctcnv, &
     &                 ndata,xdata,xmodl,xwts,npar,cxpar,elas,lbnd,ubnd,&
     &                 tdfd,cmp,knrm,norm,step,ider,nsx,swx,swy,ntyrp,  &
     &                 trp,npuls,pulse,ntypls,nchnl,topn,tcls,gstrp,    &
     &                 astrp,nfrq,freq,txcln,txa90,nrx,nrxst,xrx,yrx,   &
     &                 zrx,nstat,sx,sy,sz,rx,ry,rz,gnd_lvl,title,line,  &
     &                 nlyr,thk,res,rmu,reps,calf,ctau,cfreq,mpar
 
      
      CALL cpu_time(cmp_start)

      CALL read_system_and_survey
      CALL read_model
 
      IF ( invert ) THEN
         OPEN (nri,FILE='AirBeo.inv',STATUS='OLD')
         OPEN (nw1,FILE='AirBeo.mv1',STATUS='REPLACE')
         OPEN (md1,FILE='AirBeo.mdl',STATUS='REPLACE')
         CALL read_invrt_cntrl_and_data
      ELSE
         mprnt = 100
         js = 0
         CALL write_model(nw,mprnt,js,nlyr,thk,res,chrg,ctau,cfreq,rmu, &
     &                    reps)
         OPEN (nw1,FILE='AirBeo.mf1',STATUS='REPLACE')
      ENDIF
 
      CALL set_survey
      wrt_nw1 = .TRUE.
      IF ( tdfd==2 .AND. cmp>1 ) wrt_nw1 = .FALSE.
      IF ( wrt_nw1 ) CALL write_nw1_initial
 
      IF ( mxerr==0 ) THEN
         WRITE (*,'(/T3,A//T3,A//)')                                    &
     &                              'Control file passed initial tests.'&
     &                              , 'Computation begins.'
      ELSEIF ( mxerr==1 ) THEN
         WRITE (*,'(/T3,A//T3,A//)') 'Computation begins.' ,            &
     &                          'Look at warning messages in AirBeo.log'
      ELSEIF ( mxerr==2 ) THEN
         WRITE (*,'(/T3,A//T3,A//T3,A)')                                &
     &           'FATAL INPUT DATA ERRORS IN AirBeo.cfl' ,              &
     &          'Refer to messages in AirBeo.log' ,                     &
     &          'Execution will not occur until these are corrected.'
         STOP
      ENDIF
!============================================================================
      IF ( istop==1 ) STOP
!============================================================================
 
      IF ( tdfd==1 ) THEN
                       ! Time-Domain
 
! For time-domain, set up frequencies, interpolation times
! For time-domain, call SET_SOURCE to compute dI/dt at the transmitter using
! the DC coupling if waveform at the receiver has been specified.  Then
! (for time-domain) convert PRM_TD to the peak primary dB/dt in NT if
! impulse response is output or B in pT for step response.
! SWY will be in amps/s  * Tx area * NTRN
 
! IDER = 0 => that current derivative must be computed: ISW = 1, 11 or 31
!      = 1 => that current derivative has specified through voltage
!             calibration: ISW = 10 or 30
!      = 4 => ISW = 4  (pure rectangular pulse)
 
         ider = 0
         IF ( isw==10 .OR. isw==30 .OR. isw==130 ) ider = 1
         IF ( isw==4 ) ider = 4
         CALL set_trp
         knrm = 3
         ALLOCATE (norm(knrm))
         IF ( isw==4 ) THEN
            norm = 1.E6
         ELSE
            
            CALL dcprm_td(xrx0,yrx0,zrx0,txcln0,txarea,prm_td)
            CALL set_source(step,isw,bffac,waveform,nsx,swx,swy,prm_td)
            write(*,*) nw,bunit,bffac,kppm,punit, ppfac,prm_td
            IF ( invert ) CALL set_norm_td(nw,bunit,bffac,kppm,punit,   &
     &           ppfac,prm_td,norm)
             write(*,*) nw,bunit,bffac,kppm,punit,ppfac,prm_td,norm
         ENDIF
      ELSE
         knrm = nfrq
         ALLOCATE (norm(knrm))
         CALL dcprm_fd(nfrq,xrx,yrx,zrx,txcln,txa90,prm_fd,ppfac,norm)
      ENDIF
 
!============================================================================
      IF ( invert ) THEN
 
         DO js = 1 , nstat
            IF ( do1d==1 ) THEN
               res = res0
               thk = thk0
            ENDIF
            DO jt = 1 , ndata
               xdata(jt) = rdata(jt,js)
               xwts(jt) = rwts(jt,js)
            ENDDO
            IF ( maxval(xwts)<1 ) THEN
               WRITE (nw,99004) js
               WRITE (*,99004) js
               CYCLE
            ENDIF

!            write(*,*) 'THIS IS mynmlinv'
!            write(*,nml=mynmlinv)
            !write(*,*) size(title)
            !write(*,*) sizeof(title)

            CALL nlsq2(js,nw,nw1,nlg,mv1prt,outprt,maxits,cnvrg,pctcnv, &
     &                 ndata,xdata,xmodl,xwts,npar,cxpar,elas,lbnd,ubnd,&
     &                 tdfd,cmp,knrm,norm,step,ider,nsx,swx,swy,ntyrp,  &
     &                 trp,npuls,pulse,ntypls,nchnl,topn,tcls,gstrp,    &
     &                 astrp,nfrq,freq,txcln,txa90,nrx,nrxst,xrx,yrx,   &
     &                 zrx,nstat,sx,sy,sz,rx,ry,rz,gnd_lvl,title,line,  &
     &                 nlyr,thk,res,rmu,reps,calf,ctau,cfreq,mpar,      &
     &                 rmserr)
 
            res(1:nlyr) = mpar(1:nlyr)
            thk(1:nlyr-1) = mpar(nlyr+1:npar)
            IF(nlyr>1) THEN
               CALL cnvrt2_depth(nlyr,thk,depth)
               WRITE (md1,99001) js , syd(js) , sxd(js) , rmserr, res(1:nlyr) ,    &
     &                        depth(1:nlyr-1) , thk(1:nlyr-1)
            ELSE
               WRITE (md1,99001) js , syd(js) , sxd(js) , rmserr, res(1:nlyr)
            ENDIF
 
99001       FORMAT (I5,2F12.1,100G13.4)
 
         ENDDO
!===================================================================================
 
      ELSE    ! FORWARD MODEL OPTION
 
         CLOSE (nr)
 
         IF ( tdfd==1 ) THEN
                         ! Time-Domain.
            ALLOCATE (btd(nchnl,nstat,3))
            btd = 0.
 
!  Compute BTD, the layered earth response convolved with the excitation waveform
!  as dB/dt in nT/s if STEP = 0;  or as B in pT if STEP = 1
 
            DO js = 1 , nstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
            WRITE(*,*) 'THIS IS mynmltd before hsboss_td'
            WRITE(*,nml=mynmltd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
               CALL hsboss_td(js,step,ider,nsx,swx,swy,npuls,pulse,     &
     &                        ntypls,ntyrp,trp,nchnl,topn,tcls,txcln,   &
     &                        nstat,sz,zrx,xrx,yrx,nlyr,res,reps,rmu,   &
     &                        thk,calf,ctau,cfreq,gstrp,astrp,btd)
            ENDDO
 
!  Write out the results.
 
            CALL write_td(nw,nw1,title,nstat,line,sxd,syd,sz,txdeg,rxd, &
     &                    ryd,rz,xrx,yrx,zrx,nchnl,tms,prfl,qunit,bunit,&
     &                    ppfac,bffac,prm_td,cmp,kppm,btd)
 
         ELSE
         !  Construct the frequency-domain response.
 
            ALLOCATE (bfd(nfrq,nstat,3))
            bfd = (0.,0.)
 
            DO js = 1 , nstat
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!           WRITE(*,*) 'THIS IS mynmlfd before hsboss_fd '
!           WRITE(*,nml=mynmlfd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
               CALL hsboss_fd(js,nfrq,freq,txcln,txa90,nstat,sz,zrx,xrx,&
     &                        yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,&
     &                        bfd)
            ENDDO
            CALL write_fd(nw,nw1,title,nstat,line,txcln,txa90,sxd,syd,  &
     &                    sz,rxd,ryd,rz,config,nfrq,freq,prfl,qunit,    &
     &                    ppfac,bffac,prm_fd,cmp,bfd)
         ENDIF
      ENDIF
 
      CALL date_and_time(date,time,zone,qqdt)
      qqhms(1:2) = qqdt(5:6)
 
      CALL cpu_time(cmp_end)
      elapsed = cmp_end - cmp_start
 
      IF ( invert ) THEN
         WRITE (nw,99006) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,       &
     &                    qqdt(1) , elapsed
         WRITE (*,99006) qqhms(1:2) , qqdt(3) , month(qqdt(2)) , qqdt(1)&
     &                   , elapsed
         WRITE (nw1,99002) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,      &
     &                     qqdt(1) , elapsed
99002    FORMAT (T1,'/'/'/ AirBeo inversion completed at ',I2.2,':',    &
     &           I2.2,' on',I3.2,1X,A,I5/T1,'/ Computation time = ',    &
     &           F10.2,' seconds.')
      ELSE
         WRITE (nw,99005) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,       &
     &                    qqdt(1) , elapsed
         WRITE (*,99005) qqhms(1:2) , qqdt(3) , month(qqdt(2)) , qqdt(1)&
     &                   , elapsed
         WRITE (nw1,99003) qqhms(1:2) , qqdt(3) , month(qqdt(2)) ,      &
     &                     qqdt(1) , elapsed
99003    FORMAT (T1,'/'/'/ AirBeo forward model completed at ',I2.2,':',&
     &           I2.2,' on',I3.2,1X,A,I5/T1,'/ Computation time = ',    &
     &           F10.2,' seconds.')
      ENDIF
 
      CLOSE (nw)
      STOP
99004 FORMAT (//T14,'=============================='/T15,               &
     &        'No inversion for station',I4/T14,                        &
     &        '==============================')
99005 FORMAT (//T3,'AirBeo forward model completed at ',I2.2,':',I2.2,  &
     &        ' on',I3.2,1X,A,I5//T3,'Computation time = ',F10.2,       &
     &        ' seconds.'//)
99006 FORMAT (//T3,'AirBeo inversion completed at ',I2.2,':',I2.2,' on',&
     &        I3.2,1X,A,I5//T3,'Computation time = ',F10.2,             &
     &        ' seconds.'//)
      
      
      END PROGRAM airBeo
