        MODULE input_data_for_airBeo
 
! CONTAINS: READ_SYSTEM_AND_SURVEY, READ_MODEL, READ_INVRSION_CNTRL, , READ_INVRSION_DATA
 
      IMPLICIT NONE
 
! General Airborne & Layered Earth Dimensions
 
      INTEGER , PARAMETER :: Nprop = 7 , Ql = selected_real_kind(12,80)
      REAL , PARAMETER :: T0_min = 1.E-7 , Pi = 3.141592654 ,           &
     &                    Data_tol = 1.E-24 , Turn = Pi/10.
      INTEGER nr , nw , nd , nlg , nri , nw1 , md1 , ks , tdfd , step , &
     &        do1d , isw , prfl , nchnl , kppm , istop , nlith , nsx ,  &
     &        nsx1 , jf , jt , js , jl , jr , nfrq , ntyrp , nstat ,    &
     &        nlyr , nrx , nrxst , ntrn , mtxrx , npuls , ntypls ,      &
     &        survey , cmp , krxw , msg , mxerr , j , gstrp , astrp ,   &
     &        iunits , nppf , qqdt(8) , qqhms(2) , ndata , npar ,       &
     &        outprt , mv1prt , cnvrg , maxits , a2c , fvn , baromtrc , &
     &        line_tag
      INTEGER , ALLOCATABLE , DIMENSION(:) :: lith , kfix , cxpar ,     &
     &                                  xwts , line , cfg1
      INTEGER , ALLOCATABLE , DIMENSION(:,:) :: rwts
      REAL txampl , txfreq , pulse , pkcur , offtym , alf , delt ,      &
     &     alt1 , gnd_lvl , txcln0 , xrx0 , yrx0 , zrx0 , csf , snf ,   &
     &     dstat , prm_td(3) , bffac , ppfac , txarea , cmin , cmax ,   &
     &     t0 , qcnd , pctcnv , aline
      REAL , DIMENSION(:) , ALLOCATABLE :: freq , waveform , txon ,     &
     &                 swx , tms , wtms , topn , tcls , sz , sx , sy ,  &
     &                 zrx , xrx , yrx , bearing , fangle , depth ,     &
     &                 thk , res , rmu , reps , chrg , trp , calf ,     &
     &                 ctau , cfreq , txcln , txdeg , prm_fd , xdata ,  &
     &                 xmodl , ubnd , lbnd , elas , mpar , res0 , thk0
      REAL , DIMENSION(:,:) , ALLOCATABLE :: swy , rx , ry , rz , lyth ,&
     &                 rdata
      REAL(KIND=Ql) qfrq , qfrq6 , qfrq12 , fqq , east1 , north1
      REAL(KIND=Ql) , ALLOCATABLE :: sxd(:) , syd(:) , rxd(:,:) ,       &
     &                               ryd(:,:)
      LOGICAL invert , txa90
      CHARACTER(LEN=3) , ALLOCATABLE :: config(:)
      CHARACTER(LEN=1) tchr
      CHARACTER(LEN=120) inp , title
      CHARACTER(LEN=60) pvc
      CHARACTER(LEN=10) time , date , zone
      CHARACTER(LEN=3) month(12)
      CHARACTER(LEN=4) qunit , bunit , punit
      DATA month/'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , 'JUL' ,&
     &     'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'/
      DATA pvc , fvn/'AirBeoVR - Version 4.7.X    22 August 2016' , 450/
! Parameters for the 3D target
 
      CONTAINS
 
      SUBROUTINE read_system_and_survey
!  -----------------------=---------
 
!***  Called by: MAIN
!***      Calls: CALL CONFIG_ID, WRITE_LOG_FILE
 
!  If DO1D > 0, data to be inverted is read in this routine
 
      REAL bffacs(6) , ppfacs(4) , rho , delx , dely , a1
      CHARACTER(LEN=4) punits(4) , bunits(6)
      CHARACTER(LEN=19) wvfrm(3)
 
      DATA wvfrm/'Transmitter current' , 'Vertical receiver  ' ,        &
     &     'Horizontal receiver'/
      DATA bffacs/1. , 1000. , 1.E6 , 1. , 1000. , 1.E6/
      DATA ppfacs/1.E2 , 1.E3 , 1.E6 , 1.E9/
 
      nr = 3 !  Input unit number for AirBeo.cfl
      nw = 4 !  Output unit number for AirBeo.out
      nlg = 9
             !  Log file unit number for AirBeo.log
      nri = 13
             !  Inversion data for AirBeo.inv
      nw1 = 14
             !  Output unit number for AirBeo.mv1 (inversion) or AirBeo.mf1 (forward model)
      md1 = 15
             !  Output unit number for AirBeo.mdl
 
!      OPEN (nr,FILE=AirBeo_input//'.cfl',STATUS='OLD')
!      OPEN (nw,FILE=AirBeo_input//'.out',STATUS='REPLACE')
      
      OPEN (nr,FILE='AirBeo.cfl',STATUS='OLD')
      OPEN (nw,FILE='AirBeo.out',STATUS='REPLACE')

      CALL date_and_time(date,time,zone,qqdt)

      qqhms(1:2) = qqdt(5:6)
      WRITE (*,99031) qqhms(1:2) , qqdt(3) , month(qqdt(2)) , qqdt(1)
      WRITE (nw,99031) qqhms(1:2) , qqdt(3) , month(qqdt(2)) , qqdt(1)
 
!      Initialise some variables.
 
      mxerr = 0        !  Initialise input error flag
      nchnl = 1        !  Initialise dimensions
      nppf = 3         !  ppm default
      iunits = 5       !  pT default
      gstrp = 0
      astrp = 0
      txarea = 1.
      prm_td = 0.
      invert = .FALSE.
      txa90 = .FALSE.  !  Don't print scattered fields
 
!  Reproduce input data with no assignments and rewind file.
 
      WRITE (nw,99028) pvc
      WRITE (*,99028) pvc
      WRITE (nw,'(T11,A/T11,A/)') 'INPUT DATA' , '----------'
      REFLECT_DATA:DO jf = 1 , 10000
         READ (nr,'(A)',END=100) inp
         WRITE (nw,'(1X,A)') inp
      ENDDO REFLECT_DATA
 
 100  REWIND nr
      WRITE (nw,99001)
99001 FORMAT (T1,79('-'))
      
      title(1:120) = ' '
      READ (nr,'(A)') title
      WRITE (nw,'(/1X,A)') trim(title)
 
! Read model control & print parameters
 
      READ (nr,*) tdfd , do1d , prfl , istop
      do1d = abs(do1d)
      WRITE (nw,99002) tdfd , do1d , prfl , istop
99002 FORMAT (/T3,'TDFD =',I3,3X,'DO1D =',I3,3X,'PRFL =',I3,4X,         &
     &        'ISTOP =',I2)
      IF ( do1d>0 ) invert = .TRUE.
 
      IF ( prfl==0 .OR. prfl==2 .OR. prfl==10 ) THEN
         prfl = 2
      ELSE
         prfl = 1
      ENDIF
 
!   TDFD = 1 or 2 for TD or FD respectively.
!   DO1D = 1 => use seed model for all inversins
!        = 2 => use result of station J for seed for station J+1
!        = 0 => compute layered 1/2 space model only.
!   PRFL - indicates profile or decay curve output
!  ISTOP - read data and stop if ISTOP = 1
!        - used as a frequency setter later on in this routine.
 
      IF ( tdfd/=1 .AND. tdfd/=2 ) CALL write_log_file(nlg,1,mxerr,2)
      IF ( do1d>2 ) do1d = 1
 
      IF ( tdfd==1 ) THEN
 
! Transmitter system information
! ------------------------------
 
         nrx = 1
         nfrq = 1
         READ (nr,*) isw , nsx1 , step , iunits , nchnl , krxw , offtym
         WRITE (nw,99003) isw , nsx1 , step , iunits , nchnl , krxw ,   &
     &                    offtym
99003    FORMAT (T3,'ISW =',I4,T15,'NSX =',I4,T27,'STEP =',I2,T39,      &
     &           'UNITS =',I4,T52,'NCHNL =',I4,/T3,'KRXW =',I3,T15,     &
     &           'OFFTYM =',G12.4)
         IF ( iunits<1 .OR. iunits>3 ) THEN
            IF ( invert ) THEN
               CALL write_log_file(nlg,17,mxerr,2)
            ELSE
               CALL write_log_file(nlg,19,mxerr,1)
            ENDIF
         ENDIF
         IF ( step==0 .AND. iunits>3 ) iunits = 1       ! Default
         IF ( step==1 .AND. iunits<4 ) iunits = iunits + 3
         IF ( isw==4 ) THEN
            iunits = 6                  ! DefaulFt
            IF ( step==0 ) iunits = 3   ! Change from fT to fT/s
         ENDIF
         bunit = bunits(iunits)
         qunit = bunit
         npuls = 5
         IF ( isw==-1 ) THEN
            astrp = 1
            isw = 1
         ELSEIF ( isw==4 .OR. isw==40 ) THEN
            nsx = 1
            offtym = 0.
            IF ( isw==4 ) npuls = 1
            isw = 4
         ELSEIF ( isw==-10 .OR. isw==-30 ) THEN
            IF ( step==0 ) THEN
               gstrp = 1
               WRITE (nw,99004)
99004          FORMAT (/T3,                                             &
     &'Geotem / Questem stripping algorithm will be applied to computati&
     &ons.')
            ELSE
               WRITE (nw,99005)
99005          FORMAT (/T3,                                             &
     &'Geotem / Questem stripping algorithm is not applied for B field o&
     &utput.')
            ENDIF
         ELSEIF ( isw==-11 .OR. isw==-31 ) THEN
            WRITE (nw,99006)
99006       FORMAT (/T3,                                                &
     &'Geotem / Questem stripping algorithm is not applied for B field i&
     &nput.')
         ELSEIF ( isw==-1 .OR. isw==-4 ) THEN
            WRITE (nw,99007)
99007       FORMAT (/T3,                                                &
     &'Geotem / Questem stripping algorithm is not applied for ISW = 1 o&
     &r 4.')
         ENDIF
         isw = abs(isw)
 
         IF ( isw/=1 .AND. isw/=10 .AND. isw/=11 .AND. isw/=30 .AND.    &
     &        isw/=31 .AND. isw/=4 .AND. isw/=130 .AND. isw/=131 ) THEN
            CALL write_log_file(nlg,4,mxerr,2)
         ELSEIF ( isw==4 ) THEN
            ALLOCATE (swx(nsx),swy(nsx,3))
            IF ( step==1 ) WRITE (nw,99008) trim(bunit)
99008       FORMAT (//T10,'+-------------------------------------------'&
     &              /T10,'+   Airborne System Information'/T10,         &
     &              '+   100 Percent Duty Cycle STEP Response'/T10,     &
     &              '+   for Rectangular Waveform'/T10,                 &
     &              '+   B output will be in ',A/T10,                   &
     &              '+-------------------------------------------'/)
            IF ( step==0 ) WRITE (nw,99009) trim(bunit)
99009       FORMAT (//T10,'+-------------------------------------------'&
     &              /T10,'+   Airborne System Information'/T10,         &
     &              '+   100 Percent Duty Cycle Response'/T10,          &
     &              '+   for Rectangular Waveform'/T10,                 &
     &              '+   dB/dt output will be in ',A/T10,               &
     &              '+-------------------------------------------')
         ELSE
            IF ( step==0 ) WRITE (nw,99010) trim(bunit)
99010       FORMAT (//T10,                                              &
     &              '+----------------------------------------------'/  &
     &              T10,                                                &
     &              '+    Time-Domain AEM Impulse System Input Data '/  &
     &              T10,'+          dB/dt output will be in ',A/T10,    &
     &              '+----------------------------------------------')
            IF ( step==1 ) WRITE (nw,99011) trim(bunit)
99011       FORMAT (//T10,                                              &
     &              '+----------------------------------------------'/  &
     &              T10,                                                &
     &              '+    Time-Domain AEM Step System Input Data    '/  &
     &              T10,'+        B output will be in ',A/T10,          &
     &              '+----------------------------------------------')
            IF ( step/=1 .AND. step/=0 )                                &
     &           CALL write_log_file(nlg,5,mxerr,2)
         ENDIF
 
         IF ( krxw/=1 .AND. krxw/=2 ) CALL write_log_file(nlg,6,mxerr,2)
 
         ALLOCATE (txon(nsx1+1),waveform(nsx1+1),tms(nchnl),wtms(nchnl),&
     &             topn(nchnl),tcls(nchnl),freq(1))
 
         txon = 0.
         waveform = 0.
         tms = 0.
         wtms = 0.
         topn = 0.
         tcls = 0.
 
         IF ( isw==4 ) THEN
                         ! Step B response for full duty cycle rectangular pulse
            READ (nr,*) txfreq , txampl
            swx(1) = 0.
            pulse = .5/txfreq
            swy(1,1) = txampl
            IF ( npuls==1 ) THEN
               WRITE (nw,99012) txampl
99012          FORMAT (/T6,'Peak Current =',F6.1,' amps.'/T6,           &
     &                 'Single pulse response to step current turn-off.'&
     &                 )
            ELSE
               WRITE (nw,99013) txampl , txfreq , 1000.*pulse
99013          FORMAT (/T6,'Peak Current =',F6.1,' amps.'/T6,           &
     &                 'Step B System Frequency =',F6.1,' Hz.'/T6,      &
     &                 'Pulse On-Time =',F6.1,' ms.'/)
            ENDIF
         ELSE
            IF ( isw==1 ) THEN
               WRITE (nw,99029) wvfrm(1) , 'amps'
            ELSEIF ( isw==10 .OR. isw==11 ) THEN
                                             ! In-line component
               WRITE (nw,99029) wvfrm(3) , bunit
            ELSEIF ( isw>11 ) THEN           ! Vertical comonent
               WRITE (nw,99029) wvfrm(2) , bunit
            ENDIF
 
            READ (nr,*) (txon(j),waveform(j),j=1,nsx1)
                                                   ! Read in source waveform.
            nsx = nsx1
            IF ( txon(1)>1000.*T0_min ) THEN
                                        !  Fill in 0 point if not in original data
               nsx = nsx1 + 1
               DO jt = nsx , 2 , -1
                  txon(jt) = txon(jt-1)
                  waveform(jt) = waveform(jt-1)
               ENDDO
               txon(1) = 0.
               waveform(1) = 0.
            ENDIF
 
            DO j = 1 , nsx
               WRITE (nw,'(3X,I4,F13.3,5X,G13.4)') j , txon(j) ,        &
     &                waveform(j)
            ENDDO
 
            ALLOCATE (swx(nsx),swy(nsx,3))
            swx(1:nsx) = 1.E-3*txon(1:nsx)
            pulse = 1.E-3*(offtym+txon(nsx))
         ENDIF
         IF ( krxw==1 ) THEN
!            write(*,*) nchnl
            READ (nr,*) (topn(j),tcls(j),j=1,nchnl)
            tms = (topn+tcls)/2.
            wtms = tcls - topn
!            write (*,*) (topn(j),tcls(j),j=1,nchnl)
         ELSE
            READ (nr,*) tms(1:nchnl)
            READ (nr,*) wtms(1:nchnl)
            tcls = tms + wtms/2.
            topn = tms - wtms/2.
         ENDIF
         WRITE (nw,99014)
99014    FORMAT (/T10,'Receiver Window Specifications (ms)'/T10,        &
     &           '-----------------------------------'//T8,'Window',T19,&
     &           'Open',T31,'Close',T42,'Width',T53,'Centre'/T8,        &
     &           '------',T19,'----',T31,'-----',T42,'-----',T53,       &
     &           '------')
 
         DO jt = 1 , nchnl
            WRITE (nw,'(8X,I3,2F12.3,F11.3,F12.3)') jt , topn(jt) ,     &
     &             tcls(jt) , wtms(jt) , tms(jt)
            IF ( topn(jt)<=0 ) CALL write_log_file(nlg,7,mxerr,2)
         ENDDO
         topn = 1.E-3*topn
         tcls = 1.E-3*tcls
 
! Read in Tx area, turns and the number of receivers
         write(*,*) iunits,  bffacs(iunits)
         READ (nr,*) txcln0 , cmp , kppm
         IF ( kppm>0 ) THEN
            READ (nr,*) nppf
            IF ( nppf<1 .OR. nppf>4 ) nppf = 3
            kppm = abs(kppm)
         ENDIF
         punit = punits(nppf)
         bffac = bffacs(iunits)
                             !  Field unit conversion
         ppfac = ppfacs(nppf)
                             !  Normalisation conversion
         IF ( kppm>0 ) qunit = punit
 
         IF ( invert ) THEN
            IF ( cmp/=11 .AND. cmp/=13 .AND. cmp/=2 .AND. cmp/=3 .AND.  &
     &           cmp/=4 ) CALL write_log_file(nlg,8,mxerr,2)
 
         ELSEIF ( cmp/=11 .AND. cmp/=13 .AND. cmp/=2 .AND. cmp/=3 ) THEN
            cmp = 3
            CALL write_log_file(nlg,9,mxerr,1)
         ENDIF
         IF ( kppm/=0 .AND. kppm/=1 .AND. kppm/=3 .AND. kppm/=123 .AND. &
     &        kppm/=4 ) THEN
            kppm = 123
            CALL write_log_file(nlg,10,mxerr,2)
         ENDIF
 
! Normalisation isn't defined for step output and dB/dt waveform calibration
! or impulse output and B calibration
! It isn't used for the pure rectangular step waveform.
 
         IF ( (step==1) .AND. (isw==10 .OR. isw==30 .OR. isw==130) )    &
     &        THEN
            step = 0
            CALL write_log_file(nlg,11,mxerr,1)
         ELSEIF ( (step==0) .AND. (isw==11 .OR. isw==31 .OR. isw==131) )&
     &            THEN
            step = 1
            CALL write_log_file(nlg,12,mxerr,1)
         ENDIF
 
         !!!!!!!! IF ( isw==4 ) kppm = 0
         IF ( isw==4 ) kppm = 0
         IF ( kppm==123 ) THEN
            IF ( cmp==11 ) kppm = 1
            IF ( cmp==13 ) kppm = 3
         ENDIF
         WRITE (nw,99015) cmp , kppm , txcln0
99015    FORMAT (/T3,'CMP =',I3,4X,'KPPM =',I4/T3,                      &
     &           'Inclination angle of transmitter in level flight =',  &
     &           F5.1,' degrees (front up)')
 
         IF ( isw==1 ) THEN
            READ (nr,*) txarea , ntrn
            WRITE (nw,99016) nint(txarea) , ntrn
99016       FORMAT (/T3,'Tx area =',I8,' m^2;    NTRN =',I2)
         ELSEIF ( isw>100 ) THEN   !  ISW > 100 => central loop system
            READ (nr,*) txarea
            WRITE (nw,99017) nint(txarea)
99017       FORMAT (/T3,'Tx area =',I8)
         ENDIF
 
         IF ( isw==1 ) waveform = waveform*ntrn*txarea
 
         READ (nr,*) zrx0 , xrx0 , yrx0
         WRITE (nw,99018) zrx0 , xrx0 , yrx0
99018    FORMAT (/T3,'Initial Rx offset relative to Tx:',F7.1,' Below,',&
     &           F7.1,' Behind,',F6.1,' Left')
         rho = abs(xrx0) + abs(yrx0)
         IF ( rho<1. .AND. kppm>0 ) kppm = 3
 
      ELSEIF ( tdfd==2 ) THEN            ! Frequency-domain systems
 
         nchnl = 1
         nsx = 1
         ntyrp = 1
         ALLOCATE (swx(1),swy(1,3),trp(1),txon(1),waveform(1),tms(1),   &
     &             wtms(1),topn(1),tcls(1))
         WRITE (nw,99019)
99019    FORMAT (/10X,                                                  &
     &           '+------------------------------------------------+'/10&
     &           X,'+  Frequency-Domain Airborne System Information  +'/&
     &           10X,                                                   &
     &           '+------------------------------------------------+')
         IF ( iunits<4 ) iunits = iunits + 3
                                         ! Convert to B.  Default was 5
         bunit = bunits(iunits)
         READ (nr,*) nfrq , cmp , nppf
         IF ( nppf<1 .OR. nppf>4 ) nppf = 3
 
         IF ( cmp==-1 ) THEN
            txa90 = .TRUE.
            cmp = 1
            WRITE (nw,99020)
99020       FORMAT (/T3,                                                &
     &              'System orientation = vertical coplanar broadside')
            IF ( nppf<1 .OR. nppf>4 ) nppf = 3
         ENDIF
         IF ( cmp<1 .OR. cmp>3 ) CALL write_log_file(nlg,13,mxerr,1)
 
         IF ( invert .AND. cmp/=1 ) CALL write_log_file(nlg,14,mxerr,2)
         punit = punits(nppf)
         qunit = punit
         IF ( cmp>1 ) qunit = bunit
         ppfac = ppfacs(nppf)
         bffac = bffacs(iunits)
         write(*,*) iunits,bffac
         WRITE (nw,99021) nfrq , cmp , nppf , qunit
99021    FORMAT (/T3,'NFRQ =',I3,';  CMP =',I2,';  NPPF =',I2/T3,       &
     &           'Data will be expressed as ',A//T3,                    &
     &           'Frequencies, Tx Angles and Receiver Offset(s)'//T6,   &
     &           'Frequency  TXCLN  TXAZM   ZRX   XRX   YRX   CONFIG'/T6&
     &           ,'---------  -----  -----   ---   ---   ---   ------')
 
         nrx = nfrq
         nrxst = nfrq
 
         ALLOCATE (prm_fd(nfrq),freq(nfrq),zrx(nfrq),xrx(nfrq),yrx(nfrq)&
     &             ,txcln(nfrq),txdeg(nfrq),config(nfrq),cfg1(nfrq))
         zrx = 0.
         xrx = 0.
         yrx = 0.
         freq = 0
         txcln = 0.
         prm_fd = 0.
 
         DO jf = 1 , nfrq
            IF ( txa90 ) THEN
               txcln(jf) = 90.
               READ (nr,*) freq(jf) , zrx(jf) , xrx(jf) , yrx(jf)
            ELSE
               READ (nr,*) freq(jf) , zrx(jf) , xrx(jf) , yrx(jf) ,     &
     &                     txcln(jf)
            ENDIF
         ENDDO
 
         CALL config_id(nfrq,txcln,txa90,xrx,yrx,zrx,config,cfg1)
 
         a1 = 0.
         IF ( txa90 ) a1 = 90.
         DO jf = 1 , nfrq
            IF ( cmp<2 ) THEN
               WRITE (nw,99030) jf , freq(jf) , txcln(jf) , a1 , zrx(jf)&
     &                          , xrx(jf) , yrx(jf) , config(jf)
            ELSE
               WRITE (nw,99030) jf , freq(jf) , txcln(jf) , a1 , zrx(jf)&
     &                          , xrx(jf) , yrx(jf)
            ENDIF
         ENDDO
      ENDIF                 !  End frequency-domain specifics
 
! Flight path details for forward modelling only.  Convert FANGLE & TXCLN to radians
 
      IF ( invert ) THEN
         READ (nr,*) nstat
         IF ( nstat>1 ) THEN
            DO jf = 1 , nstat
               READ (nr,'(A)') inp
            ENDDO
         ENDIF
      ELSE
         READ (nr,*) nstat , survey , baromtrc , line_tag
         WRITE (nw,99022) nstat , survey , baromtrc , line_tag
99022    FORMAT (//T3,'NSTAT =',I4,3X,'SURVEY =',I2,3X,'BAROMTRC =',I2, &
     &           3X,'LINE_TAG =',I2)
         survey = abs(survey)
         IF ( survey<1 .OR. survey>3 ) THEN
            CALL write_log_file(nlg,16,mxerr,2)
         ELSEIF ( tdfd==2 .AND. abs(survey)>2 ) THEN
            CALL write_log_file(nlg,16,mxerr,2)
         ENDIF
 
!  NRXST is used to dimension the transmitter-receiver offsets and transmitter
!  orientation.  For time-domain systems, the offset is constant with frequency
!  but can vary with station => NRXST = NSTAT
 
!  With frequency-domain systems, the offset is constant along the survey but
!  can vary with frequency => NRXST = NFRQ
 
!  NRX is used to dimension absolute receiver locations as a function of
!  frequency, so NRX = NFRQ for frequency-domain systems but
!  NRX = 1 for time-domain systems
 
         IF ( tdfd==1 ) THEN
            nrxst = nstat
            ALLOCATE (zrx(nstat),xrx(nstat),yrx(nstat),txcln(nstat),    &
     &                txdeg(nstat))
            txcln = txcln0
            zrx = zrx0
            xrx = xrx0
            yrx = yrx0
         ENDIF
 
         ALLOCATE (line(nstat),sx(nstat),sy(nstat),sz(nstat),           &
     &             fangle(nstat),bearing(nstat),sxd(nstat),syd(nstat),  &
     &             rx(nstat,nrx),ry(nstat,nrx),rz(nstat,nrx),           &
     &             rxd(nstat,nrx),ryd(nstat,nrx))
 
         line = 1000
         sx = 0.
         sy = 0.
         sz = 0.
         fangle = 0.
         bearing = 0.
         sxd = 0.
         syd = 0.
         rx = 0.
         ry = 0.
         rz = 0.
         rxd = 0.
         ryd = 0.
 
! Read in course for forward modelling
! Read in course + data to be inverted for inversion
 
             ! Forward modelling only
 
         IF ( survey==1 ) THEN
            IF ( line_tag==1 ) THEN
               READ (nr,*) aline , syd(1) , sxd(1) , sz(1) , bearing(1) &
     &                     , dstat
               line(1) = floor(aline)
            ELSE
               READ (nr,*) syd(1) , sxd(1) , sz(1) , bearing(1) , dstat
            ENDIF
 
            WRITE (nw,99023) bearing(1)
99023       FORMAT (T3,'The flight path follows an angle of',F5.0,      &
     &              ' degrees East of North.')
            fangle(1:nstat) = bearing(1)*Pi/180.
            line(2:nstat) = line(1)
            sz(2:nstat) = sz(1)
 
            DO js = 2 , nstat
               sxd(js) = sxd(js-1) + real(cos(fangle(1))*dstat,8)
               syd(js) = syd(js-1) + real(sin(fangle(1))*dstat,8)
            ENDDO
         ELSE
            DO js = 1 , nstat
               IF ( abs(survey)==2 ) THEN
                  IF ( line_tag==1 ) THEN
                     READ (nr,*) aline , syd(js) , sxd(js) , sz(js)
                     line(js) = floor(aline)
                  ELSE
                     READ (nr,*) syd(js) , sxd(js) , sz(js)
                  ENDIF
               ELSEIF ( abs(survey)==3 ) THEN
                  IF ( line_tag==1 ) THEN
                     READ (nr,*) aline , syd(js) , sxd(js) , sz(js) ,   &
     &                           txcln(js) , zrx(js) , xrx(js) , yrx(js)
                     line(js) = floor(aline)
                  ELSE
                     READ (nr,*) syd(js) , sxd(js) , sz(js) , txcln(js) &
     &                           , zrx(js) , xrx(js) , yrx(js)
                  ENDIF
               ENDIF
               IF ( js>1 ) THEN
                  delx = real(sxd(js)-sxd(js-1))
                  dely = real(syd(js)-syd(js-1))
                  rho = sqrt(delx**2+dely**2)
                  IF ( rho>0.01 ) fangle(js) = atan2(dely,delx)
               ENDIF
            ENDDO
         ENDIF
 
         IF ( abs(survey)>1 ) THEN
                                !  Bearing correction for new line
            fangle(1) = fangle(2)
            DO js = 2 , nstat - 2
               IF ( abs(fangle(js+1)-fangle(js))>Turn ) fangle(js+1)    &
     &              = fangle(js+2)
            ENDDO
         ENDIF
         bearing = fangle*180./Pi
 
         txdeg = txcln
         txcln = txcln*Pi/180.
 
         IF ( tdfd==1 ) THEN
            WRITE (nw,99024) nstat
99024       FORMAT (/T7,I3,                                             &
     &              ' transmitter positions along the flight path'//T6, &
     &'Line   Stat     East       North       Alt      Bearing    Pitch &
     &  ZRX    XRX      YRX'/T6,                                        &
     &'----   ----     ----       -----       ---      -------    ----- &
     &  ---    ---      ---'/)
            DO js = 1 , nstat
               WRITE (nw,99025) line(js) , js , syd(js) , sxd(js) ,     &
     &                          sz(js) , bearing(js) , txcln(js) ,      &
     &                          zrx(js) , xrx(js) , yrx(js)
99025          FORMAT (T1,I9,I6,2F12.1,2F10.1,F9.1,3F8.1)
            ENDDO
         ELSE
            WRITE (nw,99026) nstat
99026       FORMAT (/T7,I3,                                             &
     &              ' transmitter positions along the flight path'//T6, &
     &       'Line   Stat       East        North       Alt     Bearing'&
     &       /T6,                                                       &
     &       '----   ----       ----        -----       ---     -------'&
     &       )
            DO js = 1 , nstat
               WRITE (nw,99027) line(js) , js , syd(js) , sxd(js) ,     &
     &                          sz(js) , bearing(js)
99027          FORMAT (T1,I9,I6,2X,2F12.1,2F10.1)
            ENDDO
         ENDIF
      ENDIF
 
99028 FORMAT (T25,A/T25,'Develped by: Art Raiche'/T33,                  &
     &        'for: AMIRA project P223F'///)
 
99029 FORMAT (//T27,A/T12,'TXON (ms)      waveform in ',A/T12,          &
     &        '---------      -----------------'/)
99030 FORMAT (I3,F9.0,F8.0,F7.0,F7.1,2F6.1,T51,A)
99031 FORMAT (//T3,'AirBeo task started at ',I2.2,':',I2.2,' on',I3.2,  &
     &        1X,A,I5//)
 
      END SUBROUTINE read_system_and_survey
 
      SUBROUTINE read_model
!  ---------------------
 
!***  Called by: MAIN
!***      Calls: WRITE_LOG_FILE
 
      IMPLICIT NONE
      INTEGER qlyr
 
!  Layered Model Specification
!  ---------------------------
 
      READ (nr,*) nlyr , qlyr , nlith , gnd_lvl
      WRITE (nw,99001) nlyr , nlith , gnd_lvl
 
99001 FORMAT (//T3,'NLAYER =',I3,';   NLITH =',I3,';   GND_LVL =',F8.2)
      IF ( qlyr/=1 ) THEN
         qlyr = 1
         CALL write_log_file(nlg,50,mxerr,1)
      ENDIF
      npar = 2*nlyr - 1
 
      ALLOCATE (lyth(nlith,Nprop),lith(nlyr),depth(nlyr-1),thk(nlyr-1), &
     &          res(nlyr),thk0(nlyr-1),res0(nlyr),chrg(nlyr),calf(nlyr),&
     &          ctau(nlyr),cfreq(nlyr),rmu(nlyr),reps(nlyr),mpar(npar))
 
      depth = 0.
      thk = 0
      res = 0
      chrg = 0
      calf = 1
      ctau = 0
      cfreq = 1
      rmu = 1
      reps = 1
      lith = 0
 
!  Initialise lithology list.
 
      lyth(1:nlith,1) = -1.
                          !  blank resistivity indicator
      lyth(1:nlith,2) = -1.
                          !  blank conductance (SIG_T) indicator
      lyth(1:nlith,3) = 1.
                          !  Relative magnetic permeabilities
      lyth(1:nlith,4) = 1.
                          !  Relative dielectric constants
      lyth(1:nlith,5) = 0.
                          !  Chargeabilities
      lyth(1:nlith,6) = 0.
                          !  CTAUs
      lyth(1:nlith,7) = 1.
                          !  CFREQs
 
      WRITE (nw,99002)
99002 FORMAT (//T27,'LITHOLOGY PROPERTIES'/T27,'--------------------'// &
     &        T35,'Relative   Relative     Cole-Cole Parameters'/T9,    &
     &'Resistivity  Conductance     MU     Dielectric   CHRG    CTAU    &
     &   CFREQ'/)
      DO j = 1 , nlith
         READ (nr,*) lyth(j,1:Nprop)
         WRITE (nw,'(I4,T8,G12.4,T22,F7.1,F12.3,F11.3,F10.2,G12.3,F8.2)'&
     &          ) j , lyth(j,1:Nprop)
         IF ( lyth(j,1)<0 .AND. lyth(j,2)<0 )                           &
     &        CALL write_log_file(nlg,53,mxerr,2)
 
         IF ( lyth(j,3)<0.01 ) lyth(j,3) = 1.
                                          ! Default RMU
         IF ( lyth(j,4)<0.01 ) lyth(j,4) = 1.
                                          ! Default REPS
 
 
         IF ( lyth(j,5)<1.E-3 .OR. lyth(j,6)<1.E-15 .OR. lyth(j,7)      &
     &        <1.E-6 ) THEN
            lyth(j,5) = 0
                     ! default CHRG
            lyth(j,6) = 0
                     ! default CTAU
            lyth(j,7) = 1
                     ! default CFRQ
         ENDIF
 
      ENDDO
 
      WRITE (nw,99003)
99003 FORMAT (//T3,'LAYERED EARTH INPUT DATA'/T3,                       &
     &        '------------------------'/)
      IF ( nlyr>1 ) THEN
         DO j = 1 , nlyr - 1
            READ (nr,*) lith(j) , thk(j)
            WRITE (nw,'(2I4,F7.1,T19,A)') j , lith(j) , thk(j) ,        &
     &             'J, LITH(J), THK(J)'
         ENDDO
      ENDIF
      READ (nr,*) lith(nlyr)
      WRITE (nw,'(2I4,T22,A)') nlyr , lith(nlyr) , 'Basement Lithology'
 
      DO jl = 1 , nlyr
         j = lith(jl)
 
         IF ( j<1 .OR. j>nlith ) THEN
            WRITE (nw,'(T3,A,I2,A,I4)') 'LITH(' , jl , ') =' , j
            CALL write_log_file(nlg,54,mxerr,2)
         ENDIF
 
         res(jl) = lyth(j,1)
         IF ( res(jl)<0 ) CALL write_log_file(nlg,55,mxerr,2)
 
         rmu(jl) = lyth(j,3)
         reps(jl) = lyth(j,4)
         chrg(jl) = lyth(j,5)
         ctau(jl) = lyth(j,6)
         cfreq(jl) = lyth(j,7)
 
         calf(jl) = 1. - chrg(jl)
      ENDDO
      res0 = res
      thk0 = thk
      mpar(1:nlyr) = res(1:nlyr)
      mpar(nlyr+1:npar) = thk(1:nlyr-1)
 
      END SUBROUTINE read_model
 
 
      SUBROUTINE read_invrt_cntrl_and_data
!  ------------------------------------
 
!***  Called by: MAIN
!***      Calls: WRITE_LOG_FILE
 
      INTEGER nfix , order , mdchnl , ndcmp , n0stat , n0chnl , n2b ,   &
     &        n2e , n3b , n3e , ctype , lyr_indx , kpar , nsta , n0pts ,&
     &        jp , jp1 , j1
      INTEGER , ALLOCATABLE , DIMENSION(:) :: kcmp , k0stat , k0chnl
      REAL tdata , e1 , e2 , e3 , a1 , a2 , delx , dely , rho
      REAL , ALLOCATABLE , DIMENSION(:) :: qdata , q2data , data_floor
      CHARACTER(LEN=1) tchr
      CHARACTER(LEN=3) kcmpc(0:5)
      CHARACTER(LEN=11) lyr_prm(2)
      DATA kcmpc/'   ' , 'HCP' , 'VCP' , 'VCA' , 'VCB' , 'HCA'/

! ----------------------------------------------------------------
! Set inversion dimensions:
!
!  NCHNL = number of time domain channels
!  NDATA = total number of readings per station to be inverted(TD or FD)
!        = NCHNL for time-domain when CMP = 11, 13, 4, 42, 43
!        = 2* NCHNL for time-domain when CMP = 2
!        = 3* NCHNL for time-domain when CMP = 3
!
!        = 2 * NFRQ for frequency-domain
!
!  RDATA & RWTS are data and weights in array form (NDATA, NSTAT)
!  XDATA & XWTS are data and weights in individual columns of RDATA & RWTS
!  RWTS  & XWTS are now restricted to integer values of 0 or 1 (reject or accep)
!  Note that data are normalised so that balance weighting is unnecessary.
!
!  XMODL contains model results in column form (NDATA)
!  NRX = 1 for TD & NFRQ for FD
!
!  For TD, RDATA is ordered as vertical components for all delay times
!  followed by all in-line components, followed by all transverse
!  components for each station.
!
!  For FD data, in-phase data is followed by quadrature data
!
!  XDATA are the RDATA stacked into a column, station by station
!  The same convention applies to XMODL, XWTS & RWTS
! ----------------------------------------------------------------
 
!  Set degree of constraint on each parameter
!  CXPAR = 0 => parameter is completely free to vary as dictated by inversion step
!        = 1 => parameter is fixed
!        = 2 => parameter is constrained by elasticity.
!        = 3 => parameter bounds are buffered.
 
 
      lyr_prm(1) = 'Resistivity'
      lyr_prm(2) = ' Thickness '
 
      IF ( tdfd==1 ) THEN
         ndata = nchnl
         IF ( cmp==2 ) ndata = 2*nchnl
         IF ( cmp==3 ) ndata = 3*nchnl
      ELSE
         ndata = 2*nfrq
      ENDIF
      ALLOCATE (elas(npar),lbnd(npar),ubnd(npar),cxpar(npar))
      cxpar = 0
      IF ( tdfd==1 ) ALLOCATE (data_floor(1),kcmp(1))
      IF ( tdfd==2 ) ALLOCATE (data_floor(2*nfrq),kcmp(nfrq))
 
      data_floor = 0.
      e1 = 1.
      elas = 1.
      e2 = 1.
      lbnd = 1.
      e3 = 1.
      ubnd = 1.
 
      WRITE (nw,99001)
 
99001 FORMAT (//T3,'-------------------------'/T3,                      &
     &        'Inversion Controls & Data'/T3,                           &
     &        '-------------------------')
      IF ( tdfd==1 ) THEN
         IF ( cmp==13 ) WRITE (nw,99002)
99002    FORMAT (/T3,'Inversion of Time-Domain Vertical Component Data')
         IF ( cmp==11 ) WRITE (nw,99003)
99003    FORMAT (/T3,'Inversion of Time-Domain In-Line Component Data')
         IF ( cmp==2 ) WRITE (nw,99004)
99004    FORMAT (/T3,                                                   &
     &'Joint Inversion of Time-Domain Vertical & In-Line Component Data'&
     &)
         IF ( cmp==3 ) WRITE (nw,99005)
99005    FORMAT (/T3,                                                   &
     &           'Joint Inversion of Time-Domain Three Component Data')
         IF ( kppm==0 ) THEN
            WRITE (nw,99006) trim(qunit)
99006       FORMAT (T3,'The data to be inverted is expressed as ',A)
         ELSE
            WRITE (nw,99043) trim(qunit)
         ENDIF
      ELSEIF ( tdfd==2 ) THEN
         WRITE (nw,99007)
99007    FORMAT (/T3,'Inversion of Frequency-Domain Data')
         WRITE (nw,99043) trim(qunit)
      ENDIF
      IF ( do1d==1 ) WRITE (nw,99008)
99008 FORMAT (/T3,                                                      &
     &        'The same starting model will be used for all inversions')
      IF ( do1d==2 ) WRITE (nw,99009)
99009 FORMAT (/T3,                                                      &
     &'The starting model for all inversions after the first will be the&
     &',/T3,'final model from the previous inversion.')
      WRITE (nw,99010) npar , ndata
99010 FORMAT (/T3,'NPAR =',I3,3X,'MCHNL =',I3)
      READ (nr,*) maxits , cnvrg , nfix , mv1prt , outprt
      WRITE (nw,99011) maxits , cnvrg , nfix , mv1prt , outprt
99011 FORMAT (T3,'MAXITS =',I3,3X,'CNVRG =',I2,3X,'NFIX =',I3,3X,       &
     &        'MV1PRT =',I2,3X,'OUTPRT =',I2)
      IF ( mv1prt<0 .OR. mv1prt>3 ) mv1prt = 1
      IF ( outprt<0 .OR. outprt>3 ) outprt = 1
      IF ( cnvrg/=1 .AND. cnvrg/=2 ) THEN
         cnvrg = 1
         CALL write_log_file(nlg,205,mxerr,1)
      ENDIF
 
      IF ( cnvrg==2 ) THEN
         READ (nr,*) pctcnv
         WRITE (nw,99012) pctcnv
99012    FORMAT (T3,                                                    &
     &       'AirBeo will finish when either the RMS error is less than'&
     &       ,F6.1,' percent'/T3,'or after',I3,                         &
     &       ' iterations, whichever comes first.')
      ELSE
         WRITE (nw,99013) maxits
99013    FORMAT (T3,                                                    &
     &'AirBeo will run until the error can no longer be reduced or after&
     &',I3,' iterations,'/T3,'whichever comes first.')
      ENDIF
 
      IF ( nfix>0 ) THEN
         WRITE (nw,99014)
99014    FORMAT (//T12,'Constrained Parameters'/T12,                    &
     &           '----------------------'//T5,                          &
     &           'Global  Layer  Parameter'/T5,                         &
     &'Index   Index   Index      Parameter      Elasticity  Lower Bound&
     &   Upper Bound   CXPAR'/T5,                                       &
     &'------  -----  ---------   ---------      ----------  -----------&
     &   -----------   -----')
         DO jp = 1 , nfix
            READ (nr,*) ctype 
            
            BACKSPACE nr
            SELECT CASE (ctype)     ! J1 is a dummy variable
            CASE (1)
               READ (nr,*) j1 , lyr_indx , kpar
               e1 = 0.
            CASE (2)
               READ (nr,*) j1 , lyr_indx , kpar , e1
            CASE (3)
               READ (nr,*) j1 , lyr_indx , kpar , e1 , e2 , e3
            END SELECT
            e1 = abs(e1)
            IF ( e1<0.05 ) e1 = 0.
                                 ! Hold for elasticities < 0.05
            IF ( e1>0.95 ) e1 = 1.
                                 ! Allow full freedom for elasticities > 0.95
 
            lyr_indx = abs(lyr_indx)
            IF ( kpar/=1 .AND. kpar/=2 ) CYCLE
            jp1 = lyr_indx
            IF ( kpar==2 ) jp1 = nlyr + lyr_indx
            IF ( jp1>npar ) CYCLE
 
            IF ( abs(e1)>0.95 ) e1 = 1.
                                       ! Allow fuLl freedom for elasticities > 0.95
            IF ( e2>e3 ) THEN          ! Switch if LB > UB
               a1 = e3
               e3 = e2
               e2 = a1
            ENDIF
            a1 = e3 - e2
            a2 = 0.005*(abs(e2)+abs(e3))
            IF ( a1<a2 ) e1 = 0.
 
            cxpar(jp1) = ctype
            IF ( abs(e1)<0.05 ) THEN
                                    ! Hold for elasticities < 0.05
               e1 = 0.
               cxpar(jp1) = 1
            ENDIF
            WRITE (nw,99015) jp1 , lyr_indx , kpar , lyr_prm(kpar) ,    &
     &                       e1 , e2 , e3 , cxpar(jp1)
99015       FORMAT (T7,I2,T14,I2,T23,I1,T31,A,T50,F4.2,T59,G12.4,T73,   &
     &              G12.4,T87,I3)
 
            elas(jp1) = e1
            lbnd(jp1) = e2
            ubnd(jp1) = e3
         ENDDO
         WRITE (nw,99016)
99016    FORMAT (/T3,90('-'))
      ELSE
         WRITE (nw,99017)
99017    FORMAT (/T3,                                                   &
     &   'All model parameters will be allowed to vary during inversion'&
     &   )
      ENDIF
 
!  Start reading from AirBeo.inv on UNIT NRI = 13
!  First skip over al comment lines
 
      DO
         READ (nri,'(A)') tchr
         IF ( tchr/='\' .AND. tchr/='/' ) EXIT
      ENDDO
      BACKSPACE (nri)
 
      IF ( tdfd==1 ) THEN
         READ (nri,*) nstat , survey , baromtrc , kcmp(1) , order
         WRITE (nw,99018) nstat , survey , baromtrc , kcmp(1) , order
99018    FORMAT (/T3,'NSTAT =',I4,3X,'SURVEY =',I2,3X,'BAROMTRC =',I2,  &
     &           3X,'KCMP =',I4,3X,'ORDER =',I2)
         SELECT CASE (cmp)
         CASE (11)
            IF ( kcmp(1)==3 ) CALL write_log_file(nlg,202,mxerr,2)
         CASE (13)
            IF ( kcmp(1)==1 ) CALL write_log_file(nlg,203,mxerr,2)
         CASE (2)
            IF ( kcmp(1)==1 ) CALL write_log_file(nlg,210,mxerr,2)
            IF ( kcmp(1)==3 ) CALL write_log_file(nlg,211,mxerr,2)
         CASE (3)
            IF ( kcmp(1)==1 ) CALL write_log_file(nlg,212,mxerr,2)
            IF ( kcmp(1)==3 ) CALL write_log_file(nlg,213,mxerr,2)
            IF ( kcmp(1)>4 .AND. kcmp(1)<100 )                          &
     &           CALL write_log_file(nlg,214,mxerr,2)
         END SELECT
 
      ELSEIF ( tdfd==2 ) THEN
         READ (nri,*) nstat , survey , baromtrc , kcmp(1:nfrq) , order
         WRITE (nw,99019) nstat , survey , baromtrc , order
99019    FORMAT (//T3,'NSTAT =',I4,3X,'SURVEY =',I2,3X,'BAROMTRC =',I2, &
     &           3X,'ORDER = ',I4)
         IF ( maxval(freq)<1.E5 ) THEN
            WRITE (nw,99020) (j,kcmpc(kcmp(j)),freq(j),j=1,nfrq)
99020       FORMAT (/T3,'Inversion components:',20(I5,': ',A,F8.1))
         ELSE
            WRITE (nw,99021) (j,kcmpc(kcmp(j)),freq(j),j=1,nfrq)
99021       FORMAT (/T3,'Inversion components:',20(I5,': ',A,F10.1))
         ENDIF
         DO jf = 1 , nfrq
            IF ( cfg1(jf)/=kcmp(jf) ) THEN
               CALL write_log_file(nlg,215,mxerr,2)
               WRITE (nlg,99022) kcmp(kcmp(1:nfrq))
99022          FORMAT (/T3,'Components from AirBeo.inv:',I4,20I5)
               WRITE (nlg,99046) kcmpc(kcmp(1:nfrq))
               WRITE (nlg,99023) cfg1(1:nfrq)
99023          FORMAT (/T3,'Components from AirBeo.cfl:',I4,20I5)
               WRITE (nlg,99046) kcmpc(cfg1(1:nfrq))
               EXIT
            ENDIF
         ENDDO
      ENDIF
 
      IF ( abs(survey)==1 ) CALL write_log_file(nlg,206,mxerr,2)
 
!        SET SYSTEM DIMENSIONS
!        ---------------------
!  NRXST is used to dimension the transmitter-receiver offsets and transmitter
!  orientation.  For time-domain systems, the offset is constant with frequency
!  but can vary with station => NRXST = NSTAT
 
!  With frequency-domain systems, the offset is constant along the survey but
!  can vary with frequency => NRXST = NFRQ
 
!  NRX is used to dimension absolute receiver locations as a function of
!  frequency, so NRX = NFRQ for frequency-domain systems but
!  NRX = 1 for time-domain systems
 
      IF ( tdfd==1 ) THEN
         nrxst = nstat
         ALLOCATE (zrx(nstat),xrx(nstat),yrx(nstat),txcln(nstat),       &
     &             txdeg(nstat))
         txcln = txcln0
         zrx = zrx0
         xrx = xrx0
         yrx = yrx0
      ENDIF
 
      ALLOCATE (line(nstat),sx(nstat),sy(nstat),sz(nstat),fangle(nstat),&
     &          bearing(nstat),sxd(nstat),syd(nstat),rx(nstat,nrx),     &
     &          ry(nstat,nrx),rz(nstat,nrx),rxd(nstat,nrx),             &
     &          ryd(nstat,nrx))
 
      sx = 0.
      sy = 0.
      rx = 0.
      ry = 0.
      rz = 0.
      rxd = 0.
      ryd = 0.
      fangle = 0.
 
      IF ( tdfd==1 ) THEN
         IF ( kcmp(1)/=1 .AND. kcmp(1)/=3 .AND. kcmp(1)/=13 .AND.       &
     &        kcmp(1)/=31 .AND. kcmp(1)/=123 .AND. kcmp(1)/=321 )       &
     &        CALL write_log_file(nlg,209,mxerr,2)
         IF ( kcmp(1)<100 .AND. cmp==3 ) THEN
            mxerr = 2
            WRITE (nlg,99024) kcmp(1)
99024       FORMAT (T3,                                                 &
     &         'Three cpomponent inversion has been specified (CMP = 3)'&
     &         /T3,'but three components are not read in: KCMP = ',I4)
         ELSEIF ( kcmp(1)<13 ) THEN
            IF ( cmp==2 ) THEN
               mxerr = 2
               WRITE (nlg,99025) kcmp(1)
99025          FORMAT (T3,                                              &
     &            'Two component inversion has been specified (CMP = 2)'&
     &            /T3,'but two components are not read in: KCMP = ',I4)
            ELSEIF ( cmp==13 .AND. kcmp(1)/=3 ) THEN
               mxerr = 2
               WRITE (nlg,99026) kcmp(1)
99026          FORMAT (T3,                                              &
     &      'Vertical component inversion has been specified (CMP = 13)'&
     &      /T3,'but this component is not read in: KCMP = ',I4)
            ELSEIF ( cmp==11 .AND. kcmp(1)/=1 ) THEN
               mxerr = 2
               WRITE (nlg,99027) kcmp(1)
99027          FORMAT (T3,                                              &
     &       'In-line component inversion has been specified (CMP = 11)'&
     &       /T3,'but this component is not read in: KCMP = ',I4)
            ENDIF
         ENDIF
         READ (nri,*) data_floor(1)
         WRITE (nw,99044) abs(data_floor(1)) , trim(qunit)
         ndcmp = 3
         IF ( kcmp(1)<100 ) ndcmp = 2
         IF ( kcmp(1)<13 ) ndcmp = 1
         mdchnl = ndcmp*nchnl
!               write(*,*) mdchnl, ndcmp,nchnl
      ELSEIF ( tdfd==2 ) THEN
         mdchnl = ndata
         READ (nri,*) data_floor(1:2*nfrq)
 
         WRITE (nw,99028) trim(qunit)
99028    FORMAT (/T8,'Frequency Data Floors (',A,')'/T8,                &
     &           '----------------------------'//T8,                    &
     &           'Freq     In-phase   Quadrature'/)
         DO jp = 1 , nfrq
            WRITE (nw,99029) jp , freq(jp) , abs(data_floor(jp)) ,      &
     &                       abs(data_floor(jp+nfrq))
99029       FORMAT (I4,F9.0,2G12.4)
         ENDDO
      ENDIF
 
! Invert one station at a time
      ALLOCATE (rdata(ndata,nstat),rwts(ndata,nstat),xdata(ndata),      &
     &          xwts(ndata),xmodl(ndata),qdata(mdchnl),q2data(mdchnl))
      xdata = 0.
      rdata = 0.
      xwts = 1
      rwts = 1
 
      READ (nri,*) n0stat , n0chnl , n0pts
      WRITE (nw,99030) n0stat , n0chnl , n0pts
99030 FORMAT (/T3,'N0STAT =',I4,3X,'N0CHNL =',I3,3X,'N0PTS =',I4)
 
      IF ( n0stat/=0 ) THEN
         nsta = abs(n0stat)
         ALLOCATE (k0stat(nsta))
         READ (nri,*) k0stat(1:nsta)
         IF ( n0stat>0 ) THEN
            WRITE (nw,99031) k0stat(1:nsta)
99031       FORMAT (/T3,                                                &
     &      'The data from the following stations will not be inverted:'&
     &      /T3,60I4)
            DO j1 = 1 , nsta
               rwts(1:ndata,k0stat(j1)) = 0
            ENDDO
         ELSE
            rwts = 0
            WRITE (nw,99032) k0stat(1:nsta)
99032       FORMAT (/T3,                                                &
     &         'Only data from the following stations will be inverted:'&
     &         /T3,60I4)
            DO j1 = 1 , nsta
               rwts(1:ndata,k0stat(j1)) = 1
            ENDDO
         ENDIF
         DEALLOCATE (k0stat)
      ENDIF
      IF ( n0chnl>0 ) THEN
         ALLOCATE (k0chnl(n0chnl))
         READ (nri,*) k0chnl(1:n0chnl)
         WRITE (nw,99033) k0chnl(1:n0chnl)
99033    FORMAT (/T3,                                                   &
     &        'Data from the following PCHNLs will be weighted to zero:'&
     &        /T3,60I4)
         DO j1 = 1 , n0chnl
            rwts(k0chnl(j1),1:nstat) = 0
         ENDDO
         DEALLOCATE (k0chnl)
      ENDIF
 
      IF ( n0pts>0 ) THEN
         ALLOCATE (k0stat(n0pts),k0chnl(n0pts))
         READ (nri,*) (k0chnl(j1),k0stat(j1),j1=1,n0pts)
         WRITE (nw,99034)
99034    FORMAT (/T3,                                                   &
     &'Data from the following (PCHNL, STAT) pairs will be weighted to z&
     &ero:')
         DO j1 = 1 , n0pts
            WRITE (nw,'(T3,2I4)') k0chnl(j1) , k0stat(j1)
            rwts(k0chnl(j1),k0stat(j1)) = 0
         ENDDO
         DEALLOCATE (k0stat,k0chnl)
      ENDIF
 
!======================
!      DATA ENTRY
!======================
 
      DO js = 1 , nstat
         IF ( abs(survey)==2 ) THEN
            READ (nri,*) aline , syd(js) , sxd(js) , sz(js) ,           &
     &                   qdata(1:mdchnl)
            line(js) = floor(aline)
         ELSEIF ( abs(survey)==3 ) THEN
            READ (nri,*) aline , syd(js) , sxd(js) , sz(js) , txcln(js) &
     &                   , zrx(js) , xrx(js) , yrx(js) , qdata(1:mdchnl)
            line(js) = floor(aline)
         ENDIF
 
         IF ( js>1 ) THEN
            delx = real(sxd(js)-sxd(js-1))
            dely = real(syd(js)-syd(js-1))
            rho = sqrt(delx**2+dely**2)
            IF ( rho>0.01 ) fangle(js) = atan2(dely,delx)
         ENDIF
 
 
!  Put data in the order of all Z followed by all X followed by all Y
!  depending upon which components are present.
 
         IF ( tdfd==1 ) THEN
            n2b = nchnl + 1
            n2e = 2*nchnl
            n3b = n2e + 1
            n3e = 3*nchnl
 
            IF ( kcmp(1)==1 .OR. kcmp(1)==3 ) q2data(1:nchnl)           &
     &           = qdata(1:nchnl)
            IF ( order==1 ) THEN
               IF ( kcmp(1)==31 .OR. kcmp(1)==312 ) q2data(1:n2e)       &
     &              = qdata(1:n2e)
               IF ( kcmp(1)==312 ) q2data(n3b:n3e) = qdata(n3b:n3e)
               IF ( kcmp(1)==13 ) THEN
                  q2data(1:nchnl) = qdata(n2b:n2e)
                  q2data(n2b:n2e) = qdata(1:nchnl)
               ELSEIF ( kcmp(1)==123 ) THEN
                  q2data(1:nchnl) = qdata(n3b:n3e)
                  q2data(n2b:n2e) = qdata(1:nchnl)
                  q2data(n3b:n3e) = qdata(n2b:n2e)
               ENDIF
            ELSEIF ( order==2 ) THEN
               DO jt = 1 , nchnl
                  IF ( kcmp(1)==31 ) THEN
                     q2data(jt) = qdata(2*jt-1)
                     q2data(jt+nchnl) = qdata(2*jt)
                  ELSEIF ( kcmp(1)==13 ) THEN
                     q2data(jt) = qdata(2*jt)
                     q2data(jt+nchnl) = qdata(2*jt-1)
                  ELSEIF ( kcmp(1)==123 ) THEN
                     q2data(jt) = qdata(3*jt)
                     q2data(jt+nchnl) = qdata(3*jt-2)
                     q2data(jt+2*nchnl) = qdata(3*jt-1)
                  ELSEIF ( kcmp(1)==312 ) THEN
                     q2data(jt) = qdata(3*jt-2)
                     q2data(jt+nchnl) = qdata(3*jt-1)
                     q2data(jt+2*nchnl) = qdata(3*jt)
                  ENDIF
               ENDDO
            ELSE
               CALL write_log_file(nlg,208,mxerr,2)
            ENDIF
            IF ( cmp==13 ) THEN
               rdata(1:nchnl,js) = q2data(1:nchnl)
            ELSEIF ( cmp==11 ) THEN
               IF ( ndcmp==1 ) rdata(1:nchnl,js) = q2data(1:nchnl)
               IF ( ndcmp>1 ) rdata(1:nchnl,js) = q2data(n2b:n2e)
            ELSEIF ( cmp==2 .OR. cmp==3 ) THEN
               rdata(1:ndata,js) = q2data(1:ndata)
            ELSEIF ( cmp==4 ) THEN
               DO jt = 1 , nchnl
                  tdata = q2data(jt)**2
                  IF ( ndcmp>1 ) tdata = tdata + q2data(jt+nchnl)**2
                  IF ( ndcmp==3 ) tdata = tdata + q2data(jt+2*nchnl)**2
                  rdata(jt,js) = sqrt(tdata)
               ENDDO
            ENDIF
         ELSEIF ( tdfd==2 ) THEN ! Store all in-phase data first and then all quadature data
            SELECT CASE (order)
            CASE (1122)
               rdata(1:2*nfrq,js) = qdata(1:2*nfrq)
            CASE (1212)
               DO jf = 1 , nfrq
                  rdata(jf,js) = qdata(2*jf-1)
                  rdata(jf+nfrq,js) = qdata(2*jf)
               ENDDO
            CASE (2211)
               rdata(1:nfrq,js) = qdata(nfrq+1:2*nfrq)
               rdata(nfrq+1:2*nfrq,js) = qdata(1:nfrq)
            CASE (2121)
               DO jf = 1 , nfrq
                  rdata(jf,js) = qdata(2*jf)
                  rdata(jf+nfrq,js) = qdata(2*jf-1)
               ENDDO
            CASE DEFAULT
               CALL write_log_file(nlg,220,mxerr,2)
            END SELECT
         ENDIF
      ENDDO
 
! Frequency-domain TXDEG & TXCLN already established in READ_SYSTEM_AND_SURVEY
! for both forward modelling and inversion.
 
      txdeg = txcln
      txcln = txcln*Pi/180.
 
      fangle(1) = fangle(2)
                         ! Bearing correction for start of new lines
      DO js = 2 , nstat - 2
         IF ( abs(fangle(js+1)-fangle(js))>Turn ) fangle(js+1)          &
     &        = fangle(js+2)
      ENDDO
      bearing = fangle*180./Pi
 
      DEALLOCATE (qdata,q2data)
 
! Write the data and weights in blocked format to make checking easier.
 
      IF ( tdfd==1 ) THEN
         WRITE (nw,99035)
99035    FORMAT (//T6,                                                  &
     &           'Line   Station   Bearing     East       North     Alt'&
     &           ,T68,'Vertical Component Data'/T6,                     &
     &           '----   -------   -------     ----       -----     ---'&
     &           ,T68,'-----------------------')
         DO js = 1 , nstat
            WRITE (nw,99045) line(js) , js , bearing(js) , syd(js) ,    &
     &                       sxd(js) , sz(js) , rdata(1:nchnl,js)                  ! First component
         ENDDO
         IF ( cmp==2 .OR. cmp==3 ) THEN
            n2b = nchnl + 1
            n2e = 2*nchnl
            WRITE (nw,99036)
99036       FORMAT (//T6,                                               &
     &           'Line   Station   Bearing     East       North     Alt'&
     &           ,T68,'Horizontal In-line Component Data'/T6,           &
     &           '----   -------   -------     ----       -----     ---'&
     &           ,T68,'----------------------------------')
            DO js = 1 , nstat
               WRITE (nw,99045) line(js) , js , bearing(js) , syd(js) , &
     &                          sxd(js) , sz(js) , rdata(n2b:n2e,js)                 ! TD in-line data
            ENDDO
         ENDIF
 
         IF ( cmp==3 ) THEN
            n3b = 2*nchnl + 1
            n3e = 3*nchnl
            WRITE (nw,99037)
99037       FORMAT (//T6,                                               &
     &           'Line   Station   Bearing     East       North     Alt'&
     &           ,T68,'Horizontal Transverse Component Data'/T6,        &
     &           '----   -------   -------     ----       -----     ---'&
     &           ,T68,'------------------------------------')
            DO js = 1 , nstat
               WRITE (nw,99045) line(js) , js , bearing(js) , syd(js) , &
     &                          sxd(js) , sz(js) , rdata(n3b:n3e,js)                 ! TD transverse data
            ENDDO
         ENDIF
 
!  Weights
 
         WRITE (nw,99044) abs(data_floor(1)) , trim(qunit)
         DO js = 1 , nstat
            DO jt = 1 , ndata
               IF ( abs(rdata(jt,js))<data_floor(1) ) rwts(jt,js) = 0
               IF ( isnan(rdata(jt,js))) rwts(jt,js) = 0
            ENDDO
         ENDDO
 
 
         WRITE (nw,99038)
99038    FORMAT (/T3,'Station   Vertical Weights'/T3,                   &
     &           '-------   ----------------')
         DO js = 1 , nstat
            WRITE (nw,99047) js , rwts(1:nchnl,js)  ! Solo or vertical weights
         ENDDO
 
         IF ( cmp==2 .OR. cmp==3 ) THEN
            n2b = nchnl + 1
            n2e = 2*nchnl
            WRITE (nw,99039)
99039       FORMAT (/T3,'Station   Horizontal In-line Weights'/T3,      &
     &              '-------   --------------------------')
            DO js = 1 , nstat
               WRITE (nw,99047) js , rwts(n2b:n2e,js) ! TD in-line weights
            ENDDO
         ENDIF
 
         IF ( cmp==3 ) THEN
            n3b = 2*nchnl + 1
            n3e = 3*nchnl
            WRITE (nw,99040)
99040       FORMAT (/T3,'Station   Horizontal Transverse Weights'/T3,   &
     &              '-------   -----------------------------')
            DO js = 1 , nstat
               WRITE (nw,99047) js , rwts(n3b:n3e,js) ! TD transverse weights
            ENDDO
         ENDIF
      ELSEIF ( tdfd==2 ) THEN
         WRITE (nw,99041)
99041    FORMAT (/T6,                                                   &
     &          'Line   Station   Bearing    East       North       Alt'&
     &          ,T68,'In-phase Data followed by Quadrature Data'/T6,    &
     &          '----   -------   -------    ----       -----       ---'&
     &          ,T68,'----------------------------------------')
         DO js = 1 , nstat
            WRITE (nw,99045) line(js) , js , bearing(js) , syd(js) ,    &
     &                       sxd(js) , sz(js) , rdata(1:2*nfrq,js)
         ENDDO
         WRITE (nw,99042)
99042    FORMAT (/T3,                                                   &
     &       'Station   In-phase Weights followed by Quadrature Weights'&
     &       /T3,                                                       &
     &       '-------   ----------------------------------------------')
         DO js = 1 , nstat
            DO jf = 1 , nfrq
               j1 = jf + nfrq
               IF ( abs(rdata(jf,js))<data_floor(jf) ) rwts(jf,js) = 0
               IF ( abs(rdata(j1,js))<data_floor(j1) ) rwts(j1,js) = 0
               IF ( isnan(rdata(jf,js))) rwts(jf,js) = 0
               IF ( isnan(rdata(j1,js))) rwts(j1,js) = 0
           ENDDO
            WRITE (nw,99047) js , rwts(1:2*nfrq,js)
         ENDDO
      ENDIF
99043 FORMAT (T3,'The data to be inverted has been normalised to ',A)
99044 FORMAT (/T3,'Time-Domain Data Floor =',G12.4,1X,A)
99045 FORMAT (I9,I7,F11.0,2F12.1,F9.1,300G13.4)
99046 FORMAT (/T32,20(2X,A))
99047 FORMAT (T1,I10,T13,300I2)
 
      END SUBROUTINE read_invrt_cntrl_and_data
 
      SUBROUTINE set_survey
!  ---------------------
 
!***  Called by: MAIN
 
 
      IMPLICIT NONE
 
      IF ( baromtrc==1 ) THEN
         sz = sz - gnd_lvl    ! Change barometric altitude to ground clearance
         IF ( abs(gnd_lvl)>0.01 ) WRITE (nw,99001) gnd_lvl
 
 
99001    FORMAT (/T3,                                                   &
     &        'Barometric altitude will be changed to ground clearance.'&
     &        /T3,'The vertical origin is shifted down by',F6.1)
      ENDIF
 
!  Compute receiver coordinates in both body centred and real world systems.
 
      DO js = 1 , nstat
         csf = cos(fangle(js))
         snf = sin(fangle(js))
 
         IF ( tdfd==1 ) THEN
            sx(js) = real(sxd(js))
            sy(js) = real(syd(js))
            rx(js,1) = sx(js) - xrx(js)*csf + yrx(js)*snf
            ry(js,1) = sy(js) - yrx(js)*csf - xrx(js)*snf
            rz(js,1) = sz(js) - zrx(js)
         ELSE
            sx(js) = real(sxd(js)) + 0.5*(xrx(1)*csf-yrx(1)*snf)
            sy(js) = real(syd(js)) + 0.5*(yrx(1)*csf+xrx(1)*snf)
            DO jf = 1 , nfrq
               rx(js,jf) = sx(js) - xrx(jf)*csf + yrx(jf)*snf
               ry(js,jf) = sy(js) - yrx(jf)*csf - xrx(jf)*snf
               rz(js,jf) = sz(js) - zrx(jf)
            ENDDO
         ENDIF
      ENDDO
      rxd = real(rx,kind=Ql)
      ryd = real(ry,kind=Ql)
      sxd = real(sx,kind=Ql)
                            !  SXD & SYD were midpoints in Freq domain
      syd = real(sy,kind=Ql)
                            !  Now they are Tx positions based on XRX(1), YRX(1)
 
      END SUBROUTINE set_survey
 
      SUBROUTINE set_trp
!  ------------------
 
!***  Called by: MAIN
 
!  Sets up interpolation times for FD -> TD transform which use the
!  exact 6 points per decade frequency-domain data plus 6 per decade
!  interpolated values.  These are based on a 12 point per decade
!  cosine filter derived from the Niels Christensen routine FILCOA
!  with OMEGA = .3 PI and shift 0.
 
!             OUTPUT
!             ------
 
!        TRP - array of time values for FD -> TD transformations
!      NTYRP - number of values in TRP
!     EXTENT - the latest time for which time-domain output is required.
!      PULSE - time length of one signal pulse
!     NTYPLS - number of TRP values in 1 PULSE
 
 
      REAL , PARAMETER :: Twopi = 6.2831853
      INTEGER mxtym , j1
      REAL t0 , extent
      REAL , ALLOCATABLE :: qqq(:)
      REAL(KIND=Ql) tbase , qtym , tq
 
      mxtym = 200
      ALLOCATE (qqq(mxtym))
      qqq = 0.
 
      qtym = log(10.D0)/12.D0
      qtym = exp(qtym)
      extent = 2.0*npuls*pulse
  
      t0 = minval(topn) - swx(nsx)
      t0 = max(t0,T0_min)
       
      tbase = 1.D0/dble(Twopi)
      DO j1 = 1 , mxtym
         IF ( tbase<t0 ) EXIT
         tbase = tbase/qtym
      ENDDO
      
      tq = tbase
      qqq(1) = real(tq)
      DO j1 = 2 , mxtym
         ntyrp = j1
         tq = tq*qtym
         qqq(j1) = real(tq)
         IF ( qqq(j1)<pulse ) ntypls = j1 + 2
         IF ( qqq(j1)>extent ) EXIT
      ENDDO
 
      ALLOCATE (trp(ntyrp))
      trp(1:ntyrp) = qqq(1:ntyrp)
      DEALLOCATE (qqq)
 
      END SUBROUTINE set_trp
 
      SUBROUTINE write_nw1_initial
!  ---------------------------
 
!***  Called by: MAIN
 
! Sets up the initial part of the output plotting file for inversion.
 
      INTEGER ncmp , j1
      CHARACTER(LEN=5) chz , chx , chy , cht , wez , wex , wey , wet ,  &
     &                 weq , wei
      CHARACTER(LEN=6) , DIMENSION(nfrq) :: qfrq , ifrq
 
      WRITE (nw1,99012) fvn , pvc , trim(title)
      WRITE (md1,99012) fvn , pvc , trim(title)
      WRITE (md1,99001)
99001 FORMAT (T1,                                                       &
     &'/ Station  EastTx  NorthTx  RES(1:NLYR)  DEPTH(1:NLYR-1)  THK(1:N&
     &LYR-1)')
      IF ( tdfd==1 ) THEN
         ncmp = cmp
         IF ( cmp>10 ) ncmp = 1
         chz = '  CHZ'
         wez = '  WEZ'
         chx = '  CHX'
         wex = '  WEX'
         chy = '  CHY'
         wey = '  WEY'
         cht = '  CHT'
         wet = '  WET'
 
         WRITE (nw1,99002) trim(qunit) , nstat , nchnl , ncmp
99002    FORMAT (T1,'/ UNITS=',A,3X,'NSTAT=',I3.3,3X,'NCH=',I3.3,3X,    &
     &           'NCMP=',I1)
         WRITE (nw1,99003) tms(1:nchnl)
99003    FORMAT (T1,'/ TIMES(ms)=',T17,100G13.4)
         WRITE (nw1,99004) wtms(1:nchnl)
99004    FORMAT (T1,'/ CHNL_WDTH(ms)=',T17,100G13.4)
         WRITE (nw1,99005)
99005    FORMAT (T1,'/ SURVEY=TD_AEM  PP=RX_POS')
         WRITE (nw1,99014) nlyr
         IF ( invert ) THEN
            IF ( cmp==11 ) WRITE (nw1,99013) (chx,jt,jt=1,nchnl) ,      &
     &                            (wex,jt,jt=1,nchnl)
            IF ( cmp==13 ) WRITE (nw1,99013) (chz,jt,jt=1,nchnl) ,      &
     &                            (wez,jt,jt=1,nchnl)
            IF ( cmp==2 ) WRITE (nw1,99013) (chz,jt,jt=1,nchnl) ,       &
     &                           (chx,jt,jt=1,nchnl) ,                  &
     &                           (wez,jt,jt=1,nchnl) ,                  &
     &                           (wex,jt,jt=1,nchnl)
            IF ( cmp==3 ) WRITE (nw1,99013) (chz,jt,jt=1,nchnl) ,       &
     &                           (chx,jt,jt=1,nchnl) ,                  &
     &                           (chy,jt,jt=1,nchnl) ,                  &
     &                           (wez,jt,jt=1,nchnl) ,                  &
     &                           (wex,jt,jt=1,nchnl) ,                  &
     &                           (wey,jt,jt=1,nchnl)
            IF ( cmp==4 ) WRITE (nw1,99013) (cht,jt,jt=1,nchnl) ,       &
     &                           (wet,jt,jt=1,nchnl)
         ELSE
            IF ( cmp==11 ) WRITE (nw1,99013) (chx,jt,jt=1,nchnl)
            IF ( cmp==13 ) WRITE (nw1,99013) (chz,jt,jt=1,nchnl)
            IF ( cmp==2 ) WRITE (nw1,99013) (chz,jt,jt=1,nchnl) ,       &
     &                           (chx,jt,jt=1,nchnl)
            IF ( cmp==3 ) WRITE (nw1,99013) (chz,jt,jt=1,nchnl) ,       &
     &                           (chx,jt,jt=1,nchnl) ,                  &
     &                           (chy,jt,jt=1,nchnl)
            IF ( cmp==4 ) WRITE (nw1,99013) (cht,jt,jt=1,nchnl)
         ENDIF
      ELSE
         weq = '  WEQ'
         wei = '  WEI'
         DO jf = 1 , nfrq
            qfrq(jf) = '  Q'//config(jf)
            ifrq(jf) = '  I'//config(jf)
         ENDDO
 
         WRITE (nw1,99006) trim(qunit) , nstat , nfrq
99006    FORMAT (T1,'/ UNITS=',A,3X,'NSTAT=',I3.3,3X,'NFRQ=',I2.2,3X,   &
     &           'NCMP=1')
         WRITE (nw1,99007) freq(1:nfrq)
99007    FORMAT (T1,'/ FREQS(Hz) =',60F13.2)
         WRITE (nw1,99008)
99008    FORMAT (T1,'/ SURVEY=FD_AEM  PP=RX_POS')
         IF ( invert ) THEN
            WRITE (nw1,99015) (ifrq(jf),jf,jf=1,nfrq) ,                 &
     &                        (qfrq(jf),jf,jf=1,nfrq) ,                 &
     &                        (wei,jf,jf=1,nfrq) , (weq,jf,jf=1,nfrq)
         ELSE
            WRITE (nw1,99015) (ifrq(jf),jf,jf=1,nfrq) ,                 &
     &                        (qfrq(jf),jf,jf=1,nfrq)
         ENDIF
         WRITE (nw1,99014) nlyr
      ENDIF
 
      WRITE (nw1,99009)
99009 FORMAT (T1,'/ MODEL_HEADER')
      WRITE (nw1,99010,ADVANCE='NO')
99010 FORMAT (T1,'/')
      WRITE (nw1,'(30(A,I2.2))') ('  RES_',j1,j1=1,nlyr) ,              &
     &                           (' THICK_',j1,j1=1,nlyr-1)
 
      IF ( .NOT.invert ) WRITE (nw1,99011) mpar(1:npar)
99011 FORMAT (T1,'/'/T1,'/ INITIAL_MODEL',84G13.4)
 
99012 FORMAT (T1,'/ ',I4.4,T15,'File version number'/T1,                &
     &        '/ PROGRAM_NAME=',A/T1,'/ TITLE: ',A)
99013 FORMAT (T1,'/ LINE_HEADER'/T1,                                    &
     &'/ Station  EastTx  NorthTx  AltTx  Txcln  EastRx  NorthRx  AltRx'&
     &,350(A,I3.3))
99014 FORMAT (T1,'/ LAYERS=',I2.2)
99015 FORMAT (T1,'/ LINE_HEADER'/T1,                                    &
     &     '/ Station  EastTx  NorthTx  AltTx  EastRx  NorthRx  AltRx  '&
     &     ,100(A,I2.2))
 
      END SUBROUTINE write_nw1_initial
 
      END MODULE input_data_for_airBeo
