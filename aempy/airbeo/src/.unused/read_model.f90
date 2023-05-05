!*==read_model.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
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
 
      ALLOCATE (lyth(nlith,nprop),lith(nlyr),depth(nlyr-1),thk(nlyr-1), &
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
         READ (nr,*) lyth(j,1:nprop)
         WRITE (nw,'(I4,T8,G12.4,T22,F7.1,F12.3,F11.3,F10.2,G12.3,F8.2)'&
     &          ) j , lyth(j,1:nprop)
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
