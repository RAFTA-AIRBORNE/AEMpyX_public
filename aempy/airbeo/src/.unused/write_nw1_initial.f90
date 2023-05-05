!*==write_nw1_initial.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
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
