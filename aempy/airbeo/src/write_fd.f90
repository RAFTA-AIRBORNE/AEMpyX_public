!*==write_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE write_fd(nw,nw1,title,nstat,line,txcln,txa90,sxd,syd,  &
     &                    sz,rxd,ryd,rz,config,nfrq,freq,prfl,qunit,    &
     &                    ppfac,bffac,prm_fd,cmp,bfd)
!------------------------------------------------------------------------------------------
 
!***  Called by: MAIN
!***      Calls: WRSLV_FD
 
!  Prints the results of TEM computations.
!
!                NW - output unit number
!               NW1 - unit number for mf1 file
!             NSTAT - Number of stations
!              LINE - Line number
!         TXCLN(JF) - transmitter orientation at frequency JF (in radians)
!             TXA90 - true for vertical co-planar briadside array
!              PRFL = 1 => profile output
!                   = 2 => spectral output
!             QUNIT = text units for B, dB/dt or normalisation
!             BFFAC = numeric factor for nT, nT/s, pT, fT, pT/s or fT/s output
!             PPFAC = numeric factor for output in pct, ppt, ppm or ppb
!              NFRQ - number of frequencies.
!              FREQ - array of frequencies.
!               CMP = 1: write single component response only in direction of transmitter
!                        orientation in normalised units
!                   = 2: write un-normalised vertical and horizontal in-line response. (pT)
!                   = 3: write un-normalised 3 component response. (pT)
!
!             JF, below refers to receiver location varying with frequency
!             since the offset can vary with frequency
!
!         RXD(I,JF) - North coordinates for receiver JF at station I
!         RYD(I,JF) - East coordinates for receiver JF at station I
!          RZ(I,JF) - altitude of receiver JF at station I
!
!        PRM_FD(JF) - DC (nT) response along transmitting dipole direction
!                     for Tx-Rx configuration for frequency JF.
!
!       BFD(JF,JS,JC) - layered + scattered earth magnetic field for frequency JF,
!                       source position JS, component JC (nT) in aircraft components
!
 
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80)
      INTEGER nw , nw1 , nstat , line(nstat) , nfrq , cmp , prfl ,      &
     &        txd(nfrq) , jf , js
      REAL freq(nfrq) , rz(nstat,nfrq) , sz(nstat) , mpz1(nstat) ,      &
     &     prm_fd(nfrq) , prm4 , norm(nfrq) , ytr(nfrq,nstat) ,         &
     &     txcln(nfrq) , ppfac , bffac
      REAL(KIND=Ql) rxd(nstat,nfrq) , ryd(nstat,nfrq) , mxd(nstat) ,    &
     &              myd(nstat) , sxd(nstat) , syd(nstat)
      COMPLEX bfd(nfrq,nstat,3) , bfd1(nfrq,nstat,4)
      LOGICAL txa90 , wl
      CHARACTER(LEN=120) title
      CHARACTER(LEN=11) qi , ql0 , ql1
      CHARACTER(LEN=9) qr
      CHARACTER(LEN=10) cmp2
      CHARACTER(LEN=3) config(nfrq)
      CHARACTER(LEN=4) qunit
      CHARACTER(LEN=8) cmp3
      CHARACTER(LEN=7) cmp1
      DATA cmp1 , cmp2 , cmp3/'IN-LINE' , 'TRANSVERSE' , 'VERTICAL'/
      DATA qr , qi/'IN-PHASE ' , 'QUADRATURE '/
 
! This routine assumes a single Tx-Rx separation for each frequency.
! Put all the results into  BFD1
 
!   BFD1(JF,JS,1:3) will contain the in-line(1), transverse(2) and
!                   vertical fields(3) respectively (in pT)
!                   for frequency JF and station JS
!   BFD1(JF,JS,4) will contain the maximally coupled field (ppm)
!                 normalised to the parallel primary component.
!   BFD1SC(JF,JS,1:4) contains the scattered fields
 
! Print results at Tx-Rx midpoint
 
      mxd(1:nstat) = (rxd(1:nstat,1)+sxd(1:nstat))/2.
      myd(1:nstat) = (ryd(1:nstat,1)+syd(1:nstat))/2.
      mpz1(1:nstat) = (rz(1:nstat,1)+sz(1:nstat))/2.
      bfd1 = (0.,0.)
      IF ( cmp==1 ) WRITE (nw,99001)
99001 FORMAT (//T4,'Frequency    Coupled Primary   Normalisation'/T4,   &
     &        'Index        Field (pT)        Factor'/T4,               &
     &        '---------    ---------------   ------------')
      DO jf = 1 , nfrq
         prm4 = bffac*abs(prm_fd(jf))                ! coupled primary in nT, pT or fT
         norm(jf) = ppfac/prm4                     ! for pct, ppt, ppm or ppb output
 
         IF ( cmp==1 ) WRITE (nw,'(I6,T17,G12.4,T34,g12.4)') jf , prm4 ,&
     &                        norm(jf)
 
         bfd1(jf,1:nstat,1:3) = bffac*bfd(jf,1:nstat,1:3)     ! total field in nT, pT or fT
      ENDDO
 
! Normalise them as indicated in input data file.
! For CMP = 1, compute component along Tx direction
 
      txd(1:nfrq) = nint(180.*txcln(1:nfrq)/3.1416)
      WRITE (nw,99002)
 
99002 FORMAT (//T3,'FREQUENCY-DOMAIN Airbeo OUTPUT'/T3,33('-'))
      IF ( cmp>1 ) WRITE (nw,99003)
99003 FORMAT (/T3,                                                      &
     &'The IN-LINE component is defined as the horizontal component alon&
     &g'/T3,                                                            &
     &'the flight path.  It is positive in the forward flight direction.&
     &')
      IF ( cmp>2 ) WRITE (nw,99004)
99004 FORMAT (/T3,                                                      &
     & 'The TRANSVERSE component is defined as the horizontal component'&
     & ,/T3,'perpendicular to the flight path.'/)
 
      DO jf = 1 , nfrq
 
!  maximally coupled response
         IF ( txa90 ) THEN
            bfd1(jf,1:nstat,4) = bfd1(jf,1:nstat,2)
         ELSE
            bfd1(jf,1:nstat,4) = bfd1(jf,1:nstat,1)*sin(txcln(jf))      &
     &                           + bfd1(jf,1:nstat,3)*cos(txcln(jf))
         ENDIF
 
         bfd1(jf,1:nstat,4) = norm(jf)*bfd1(jf,1:nstat,4)
      ENDDO
 
      IF ( cmp==1 ) THEN
         WRITE (nw,99005) trim(title)
99005    FORMAT (/T3,'TITLE:  ',A/T3,'-----')
         WRITE (nw,99006)
99006    FORMAT (/T3,                                                   &
     &    'SINGLE COMPONENT RESPONSE ALONG TRANSMITTER DIPOLE DIRECTION'&
     &    )
 
         WRITE (nw,99008) qr , trim(qunit)
         ytr(1:nfrq,1:nstat) = real(bfd1(1:nfrq,1:nstat,4))
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
 
         WRITE (nw,99008) qi , trim(qunit)
         ytr(1:nfrq,1:nstat) = aimag(bfd1(1:nfrq,1:nstat,4))
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
      ENDIF
 
      IF ( cmp>1 ) THEN
         WRITE (nw,'(/3X,A)') trim(title)
                                        !  Vertical component
         WRITE (nw,99009) qr , cmp3 , qunit
         ytr(1:nfrq,1:nstat) = real(bfd1(1:nfrq,1:nstat,3))
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
 
         WRITE (nw,99009) qi , cmp3 , qunit
         ytr(1:nfrq,1:nstat) = aimag(bfd1(1:nfrq,1:nstat,3))
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
 
         WRITE (nw,'(//3X,A)') trim(title)
                                         !  In-line component
         WRITE (nw,99009) qr , cmp1 , qunit
         ytr(1:nfrq,1:nstat) = real(bfd1(1:nfrq,1:nstat,1))
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
 
         WRITE (nw,99009) qi , cmp1 , qunit
         ytr(1:nfrq,1:nstat) = aimag(bfd1(1:nfrq,1:nstat,1))
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
      ENDIF
 
      IF ( cmp==3 ) THEN
         WRITE (nw,'(//3X,A)') trim(title)
                                         ! Transverse component
         ytr(1:nfrq,1:nstat) = real(bfd1(1:nfrq,1:nstat,2))
         WRITE (nw,99009) qr , cmp2 , qunit
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
 
         ytr(1:nfrq,1:nstat) = aimag(bfd1(1:nfrq,1:nstat,2))
         WRITE (nw,99009) qi , cmp2 , qunit
         CALL wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,    &
     &                 config,ytr)
      ENDIF
 
!  Finish writing Airbeo.mf1
 
      IF ( cmp>1 ) RETURN
      DO js = 1 , nstat
         wl = .TRUE.
         IF ( js>1 ) THEN
           IF (line(js)==line(js-1) ) wl = .FALSE.
         ENDIF
         IF ( wl ) THEN
            WRITE (ql0,'(I10)') line(js)
            READ (ql0,'(A)') ql1
            WRITE (nw1,'(2A)') 'Line ' , trim(adjustl(ql1))
         ENDIF
 
         WRITE (nw1,99007) js , syd(js) , sxd(js) , sz(js) , ryd(js,1) ,&
     &                     rxd(js,1) , rz(js,1) ,                       &
     &                     real(bfd1(1:nfrq,js,4)) ,                    &
     &                     aimag(bfd1(1:nfrq,js,4))
99007    FORMAT (I5,2F12.1,F8.1,2F12.1,F8.1,50G13.4)
      ENDDO
99008 FORMAT (//T10,A,'COMPONENT - ',A)
99009 FORMAT (/T10,2A,' COMPONENT - ',A)
 
      END SUBROUTINE write_fd
