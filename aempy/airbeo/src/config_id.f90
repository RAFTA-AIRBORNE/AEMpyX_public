!*==config_id.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE config_id(nfrq,txcln,txa90,xrx,yrx,zrx,config,cfg1)
!---------------------------------------------------------------
 
!***  Called by: READ_SYSTEM_AND_SURVEY
 
!  Returns CONFIG = HCP, VCA, VCP, VCB, HCA or '   '
!            CFG1 =  1    2    3    4    5       0
!
!        NFRQ - number of frequencies
!       TXCLN - transmitter inclination in degrees
!       TXA90 - true for vertical co-planar briadside array
!         ZRX - vertical receiver offset for each frequency;   below = positive
!         XRX - in-line receiver offset for each frequency;    behind = positive
!         YRX - transverse receiver offset for each frequency; left = positive.
 
      INTEGER nfrq , jf , txo , cfg1(nfrq)
      REAL xabs , yabs , zabs , rabs
      REAL , DIMENSION(nfrq) :: txcln , xrx , yrx , zrx
      LOGICAL txa90
      CHARACTER(LEN=3) config(nfrq) , kcmpc(0:5)
      DATA kcmpc/'   ' , 'HCP' , 'VCA' , 'VCP' , 'VCB' , 'HCA'/
 
      IF ( txa90 ) THEN
         config = 'VCB'
         cfg1 = 4
         RETURN
      ENDIF
 
      cfg1 = 0
      DO jf = 1 , nfrq
         xabs = abs(xrx(jf))
         yabs = abs(yrx(jf))
         zabs = abs(zrx(jf))
         rabs = sqrt(xabs**2+yabs**2)
         txo = -1
         IF ( abs(abs(txcln(jf))-90.)<0.1 ) txo = 90
         IF ( abs(abs(txcln(jf))-0.)<0.1 ) txo = 0
         IF ( zabs<0.1 .AND. rabs>1. ) THEN
            IF ( txo==0 ) cfg1(jf) = 1
            IF ( txo==90 ) THEN
               IF ( xabs<0.01 ) cfg1(jf) = 3
               IF ( yabs<0.01 ) cfg1(jf) = 2
            ENDIF
         ELSEIF ( zabs>1. .AND. rabs>0.1 ) THEN
            cfg1(jf) = 5
         ENDIF
         config(jf) = kcmpc(cfg1(jf))
      ENDDO
 
      END SUBROUTINE config_id
