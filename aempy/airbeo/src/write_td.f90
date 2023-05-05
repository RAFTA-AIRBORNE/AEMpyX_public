!*==write_td.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE write_td(nw,nw1,title,nstat,line,sxd,syd,sz,txdeg,rxd, &
     &                    ryd,rz,xrx,yrx,zrx,nchnl,tms,prfl,qunit,bunit,&
     &                    ppfac,bffac,prm_td,cmp,kppm,btd)
!----------------------------------------------------------------------------------------------
 
!***  Called by: MAIN
!***      Calls: WRSLV
!
!  Prints the results of TEM computations.
!
!                NW - output unit number
!          LINE_TAG - Line ID up to 20 alphanumeric characters
!             NSTAT - Number of stations
!            SZ(JS) - transmitter altitude at station JS
!           SXD(JS) - North coordinate of transmitter at station JS
!           SYD(JS) - East)coordinate of transmitter at station JS
!          RXD(I,1) - North coordinates for receiver at station I
!          RYD(I,1) - East coordinates for receiver at station I
!           RZ(I,1) - altitude of receiver J at station I
!          XRX (JS) - distance receiver is behind the transmitter at station JS
!          ZRX (JS) - distance receiver is below the transmitter at station JS
!          YRX (JS) - port offset of distance receiver at station JS
!             NCHNL - number of channels.
!               TMS - array of times at channel midpoints
!              PRFL = 1 for orofile output;  = 2 for decay output.
!             QUNIT = text units for B, dB/dt or normalisation
!             BFFAC = numeric factor for nT, nT/s, pT, fT, pT/s or fT/s output
!             PPFAC = numeric factor for output in pct, ppt, ppm or ppb
!       PRM_TD(1:3) - (1) in-line; (2) transverse; (3) vertical primary field
!                     nT/s for impulse, nT for step.
!            CMP = 11 => print horizontal in-line component only
!                = 13 => print vertical component only
!                =  2 => print vertical & horizontal in-line components
!                =  3 => print all three components
!
!           KPPM = 0   => no PPM normalisation
!           KPPM = 1   => all components are normalised to in-line primary field
!           KPPM = 3   => all components are normalised to vertical primary field
!           KPPM = 4   => all components are normalised to total primary field
!           KPPM = 123 => vertical component normalised to the vertical primary &
!                         in-line component to the in-line primary
!
!       BTD(JT,JS,JC) - total field for channel JT, source position JS,component JC
!          JC = 1 => in-line component; 2 => transverse component;  3 => vertical component
 
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80)
      REAL , PARAMETER :: Tol = 1.E-3
      INTEGER nw , nw1 , prfl , nstat , line(nstat) , nchnl , mchnl ,   &
     &        cmp , kppm , jc , js
      REAL tms(nchnl) , rz(nstat,1) , prm_td(3) , xbd , zbd , xrn ,     &
     &     norm(4) , ppfac , bffac , ytr(nchnl,nstat) , xq(4)
      REAL , DIMENSION(nstat) :: sz , xrx , yrx , zrx , txdeg
      REAL , DIMENSION(nchnl,nstat,3) :: btd
      REAL , ALLOCATABLE :: qdata(:)
      REAL(KIND=Ql) sxd(nstat) , syd(nstat) , rxd(nstat,1) ,            &
     &              ryd(nstat,1)
      LOGICAL wl , prtybd , prtx , prty , prtz
      CHARACTER(LEN=120) title
      CHARACTER(LEN=25) xloc(2) , xlc
      CHARACTER(LEN=10) cmp2 , ql0 , ql1
      CHARACTER(LEN=8) cmp3
      CHARACTER(LEN=7) cmp1
      CHARACTER(LEN=5) zloc(2) , zlc , ylc
      CHARACTER(LEN=4) qunit , bunit
 
      DATA zloc(1:2)/'below' , 'above'/
      DATA xloc(1:2)/'behind the transmitter.  ' ,                      &
     &     'ahead of the transmitter.'/
      DATA cmp1 , cmp2 , cmp3/'IN-LINE' , 'TRANSVERSE' , 'VERTICAL'/
 
!  Set up receiver locations and write TITLE.
 
! Normalisation isn't defined for step output and dB/dt waveform calibration
! It isn't used for the pure rectangular step waveform.
 
      prtx = .TRUE.
      prty = .TRUE.
      prtz = .TRUE.
      IF ( cmp/=3 ) prty = .FALSE.
      IF ( cmp==11 ) prtz = .FALSE.
      IF ( cmp==13 ) prtx = .FALSE.
 
      norm(1) = abs(prm_td(1))
      norm(2) = abs(prm_td(2))
      norm(3) = abs(prm_td(3))
      norm(4) = sqrt(norm(1)**2+norm(2)**2+norm(3)**2)
      IF ( kppm==0 ) THEN            !  Compute fields in requied units
         btd = bffac*btd
      ELSEIF ( kppm>0 ) THEN        !  Compute normalised response
         DO jc = 1 , 3
            xrn = norm(jc)/norm(4)
            IF ( xrn<Tol ) norm(jc) = norm(4)
         ENDDO
 
         IF ( kppm==4 ) norm(1:3) = norm(4)             ! TOTAL primary normalisation
         IF ( kppm==1 .OR. kppm==3 ) norm(1:3) = norm(kppm)
                                                        ! vertical or in-line normalisation
 
         DO jc = 1 , 3
            btd(1:nchnl,1:nstat,jc) = ppfac*btd(1:nchnl,1:nstat,jc)     &
     &                                /norm(jc)
         ENDDO
      ENDIF
 
      xbd = abs(xrx(1))
      zbd = abs(zrx(1))
      zlc = zloc(1)
      xlc = xloc(1)
      IF ( zrx(1)<0 ) zlc = zloc(2)
      IF ( xrx(1)<0 ) xlc = xloc(2)
 
      WRITE (nw,99001)
 
99001 FORMAT (//T3,'TIME-DOMAIN Airbeo OUTPUT'/T3,28('-')//T3,          &
     &'The IN-LINE component is defined as the horizontal component alon&
     &g'/T3,                                                            &
     &'the flight path.  It is positive in the forward flight direction.&
     &')
      IF ( cmp==3 ) WRITE (nw,99002)
99002 FORMAT (/T3,'The TRANSVERSE component is the horizontal component'&
     &        /T3,'perpendicular to the flight path.')
      prtybd = .FALSE.
      IF ( abs(yrx(1))>.5 ) THEN
         prtybd = .TRUE.
         ylc = 'left'
         IF ( yrx(1)<0. ) ylc = 'right'
      ENDIF
 
      xq(1:3) = bffac*prm_td(1:3)
      xq(4) = bffac*norm(4)
      IF ( kppm>0 ) THEN
                      !  Print out primary fields
         WRITE (nw,99003) xq(3) , bunit , xq(1) , bunit , xq(2) ,       &
     &                    bunit , xq(4) , bunit
99003    FORMAT (/T20,'  Vertical primary =',G12.4,1X,A/T20,            &
     &           '   In-line primary =',G12.4,1X,A/T20,                 &
     &           'Transverse primary =',G12.4,1X,A/T20,                 &
     &           '     Total primary =',G12.4,1X,A)
 
         IF ( kppm==1 ) WRITE (nw,99004)
99004    FORMAT (//T3,                                                  &
     &       'Each component is normalised to the in-line primary field'&
     &       )
         IF ( kppm==3 ) WRITE (nw,99005)
99005    FORMAT (/T3,                                                   &
     &      'Each component is normalised to the vertical primary field'&
     &      )
         IF ( kppm==4 ) WRITE (nw,99006)
99006    FORMAT (/T3,                                                   &
     &         'Each component is normalised to the total primary field'&
     &         )
         IF ( kppm==123 ) WRITE (nw,99007)
99007    FORMAT (/T3,                                                   &
     & 'Each component is normalised to its corresponding primary field'&
     & )
      ENDIF
 
      WRITE (nw,99008) zbd , zlc , xbd , xlc
99008 FORMAT (/T3,'The receiver is',F6.1,' metres ',A,' and',F6.1,      &
     &        ' metres ',A)
      IF ( prtybd ) WRITE (nw,99009) abs(yrx(1)) , ylc
99009 FORMAT (T3,'It is ',F6.1,' metres to ',A)
 
      IF ( prfl==2 ) THEN
         WRITE (nw,99010)
99010    FORMAT (/T57,'Altitude'/T3,                                    &
     &           'The first 3 rows of columns 3-8 are the transmitter', &
     &           T57,'East'/T57,'North   coordinates.'//T57,            &
     &           'Altitude'/T12,                                        &
     &           'below which are the corresponding receiver',T57,      &
     &           'East'/T57,'North   coordinates.'//T3,                 &
     &       'Underneath each receiver (Altitude, East, North) triplet,'&
     &       )
         WRITE (nw,99011) trim(qunit)
99011    FORMAT (T3,'the output is printed in ',A)
         WRITE (nw,99012)
99012    FORMAT (T3,                                                    &
     &'Channel times are in milliseconds from the start of the transmitt&
     &er pulse',/T3,'to the centre of the receiver window.')
      ENDIF
 
      IF ( prtz ) THEN
                 !  Total vertical component
         WRITE (nw,99013) trim(title)
99013    FORMAT (//T2,'TITLE:  ',A/T2,'-----')
         WRITE (nw,99015) cmp3 , trim(qunit)
         ytr(1:nchnl,1:nstat) = btd(1:nchnl,1:nstat,3)
         CALL wrslv(nw,prfl,rz,rxd,ryd,nstat,sz,sxd,syd,nchnl,tms,ytr)
      ENDIF
 
      IF ( prtx ) THEN
                  !  Total in-line component
         WRITE (nw,99015) cmp1 , trim(qunit)
         ytr(1:nchnl,1:nstat) = btd(1:nchnl,1:nstat,1)
         CALL wrslv(nw,prfl,rz,rxd,ryd,nstat,sz,sxd,syd,nchnl,tms,ytr)
      ENDIF
 
      IF ( prty ) THEN     !  Total transverse component
         WRITE (nw,99015) cmp2 , trim(qunit)
         ytr(1:nchnl,1:nstat) = btd(1:nchnl,1:nstat,2)
         CALL wrslv(nw,prfl,rz,rxd,ryd,nstat,sz,sxd,syd,nchnl,tms,ytr)
      ENDIF
 
!  Write out scattered voltages stripped of layered 1/2 space responses
!  if requested
!  ====================================================================
 
!  Finish writing Airbeo.mf1
 
      mchnl = nchnl
      IF ( cmp<4 ) mchnl = cmp*nchnl
      ALLOCATE (qdata(mchnl))
 
      DO js = 1 , nstat
         IF ( cmp/=11 ) THEN
            qdata(1:nchnl) = btd(1:nchnl,js,3)
            IF ( cmp<4 ) qdata(nchnl+1:2*nchnl) = btd(1:nchnl,js,1)
            IF ( cmp==3 ) qdata(2*nchnl+1:3*nchnl) = btd(1:nchnl,js,2)
         ELSE
            qdata(1:nchnl) = btd(1:nchnl,js,1)
         ENDIF
 
         wl = .TRUE.
         IF ( js>1 .AND. line(js)==line(js-1) ) wl = .FALSE.
         IF ( wl ) THEN
            WRITE (ql0,'(I10)') line(js)
            READ (ql0,'(A)') ql1
            WRITE (nw1,'(2A)') 'Line ' , trim(adjustl(ql1))
         ENDIF
         WRITE (nw1,99014) js , syd(js) , sxd(js) , sz(js) , txdeg(js) ,&
     &                     ryd(js,1) , rxd(js,1) , rz(js,1) ,           &
     &                     qdata(1:mchnl)
99014    FORMAT (I5,2F12.1,F8.1,F6.1,2F12.1,F8.1,200G13.4)
 
      ENDDO
99015 FORMAT (/T11,A,' COMPONENT - ',A/)
 
      END SUBROUTINE write_td
