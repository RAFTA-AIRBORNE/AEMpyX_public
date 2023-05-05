!*==write_model.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
!==================================================================================
 
      SUBROUTINE write_model(nw,mprnt,js,nlyr,thk,res,chrg,ctau,cfreq,  &
     &                       rmu,reps)
!-------------------------------------------------------------------------------
 
!***  Called by: MAIN, NLSQ2
 
!      NW - output unit
!   MPRNT = 100 - forward model
!         =  0 - initial model before inversion
!         =  ITS - model after ITS iteations
!         = -ITS - final model after ITS iterations
!    NLYR - number of layers
!     THK - layer thicknesses
!     RES - array of layer resistivities
!    REPS - relative dielectric constant
!     RMU - mu(i) / mu(0)
!    CHRG -  chargeability
!    CTAU - layer time constants
!   CFREQ - layer frequency constants
 
      IMPLICIT NONE
      INTEGER nw , mprnt , nlyr , j , js
      REAL , DIMENSION(nlyr) :: res , ctau , cfreq , rmu , reps , chrg
      REAL , DIMENSION(nlyr-1) :: thk , depth , cnd
      LOGICAL full
 
      IF ( mprnt==100 ) THEN
         WRITE (nw,99001)
 
99001    FORMAT (//T9,'Model Description'/T9,'=================')
      ELSE
         IF ( mprnt==0 ) WRITE (nw,99007) js
         IF ( mprnt==0 ) WRITE (*,99007) js
         IF ( mprnt>0 ) WRITE (nw,99008) mprnt , js
         IF ( mprnt>0 ) WRITE (*,99008) mprnt , js
         IF ( mprnt<0 ) WRITE (nw,99009) abs(mprnt) , js
         IF ( mprnt<0 ) WRITE (*,99009) abs(mprnt) , js
      ENDIF
 
      IF ( nlyr>1 ) CALL cnvrt2_depth(nlyr,thk,depth)
      full = .FALSE.
      DO j = 1 , nlyr - 1
         cnd(j) = thk(j)/res(j)
      ENDDO
 
      IF ( maxval(chrg)>1.E-4 ) full = .TRUE.
      IF ( maxval(rmu)>1.0001 ) full = .TRUE.
      IF ( maxval(reps)>1.0001 ) full = .TRUE.
      IF ( full ) THEN
         IF ( nlyr>1 ) THEN
            WRITE (nw,99002)
99002       FORMAT (/T2,                                                &
     &'Layer  Resistivity  Depth  Thickness  Conductance   MU-R   EPS-R &
     &  CHRG    CTAU      CFREQ'/T2,                                    &
     &'-----  -----------  -----  ---------  -----------   ----   ----- &
     &  ----    ----      -----')
            DO j = 1 , nlyr - 1
               WRITE (nw,99010) j , res(j) , depth(j) , thk(j) , cnd(j) &
     &                          , rmu(j) , reps(j) , chrg(j) , ctau(j) ,&
     &                          cfreq(j)
               IF ( mprnt<100 ) WRITE (*,99010) j , res(j) , depth(j) , &
     &                                 thk(j) , cnd(j)
            ENDDO
            j = nlyr
            WRITE (nw,99011) j , res(j) , rmu(j) , reps(j) , chrg(j) ,  &
     &                       ctau(j) , cfreq(j)
            IF ( mprnt<100 ) WRITE (*,99011) j , res(j)
         ELSE
            WRITE (nw,99003)
99003       FORMAT (/T2,                                                &
     &            'Resistivity   MU-R   EPS-R   CHRG    CTAU      CFREQ'&
     &            /T2,                                                  &
     &            '-----------   ----   -----   ----    ----      -----'&
     &            )
            j = nlyr
            WRITE (nw,99004) res(j) , rmu(j) , reps(j) , chrg(j) ,      &
     &                       ctau(j) , cfreq(j)
99004       FORMAT (G12.4,2F7.2,F8.2,G11.2,F7.2)
         ENDIF
      ELSEIF ( nlyr>1 ) THEN
         WRITE (nw,99005)
99005    FORMAT (/T2,'Layer  Resistivity  Depth  Thickness  Conductance'&
     &           /T2,'-----  -----------  -----  ---------  -----------'&
     &           )
         DO j = 1 , nlyr - 1
            WRITE (nw,99010) j , res(j) , depth(j) , thk(j) , cnd(j)
            IF ( mprnt<100 ) WRITE (*,99010) j , res(j) , depth(j) ,    &
     &                              thk(j) , cnd(j)
         ENDDO
         j = nlyr
         WRITE (nw,99011) j , res(j)
         IF ( mprnt<100 ) WRITE (*,99011) j , res(j)
      ELSE
         WRITE (nw,99006) res(1)
99006    FORMAT (/T3,'Host resistivity =',G15.4)
      ENDIF
99007 FORMAT (//T9,'Initial Model Before Inversion for Station',I4,/T9, &
     &        '----------------------------------------------')
99008 FORMAT (//T9,'Model after',I3,' Iterations for Station',I4,/T9,   &
     &        '-----------------------------------------')
99009 FORMAT (//T9,'Final Model after',I3,' Iterations for Station',I4, &
     &        /T9,'===============================================')
99010 FORMAT (I4,G15.4,F7.1,F9.1,T39,G12.4,2F7.2,F8.2,G11.2,F7.2)
99011 FORMAT (I4,G15.4,28X,2F7.2,F8.2,G11.2,F7.2)
 
      END SUBROUTINE write_model
