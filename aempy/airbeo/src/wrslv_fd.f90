!*==wrslv_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE wrslv_fd(nw,prfl,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp, &
     &                    config,ytr)
!------------------------------------------------------------------------------
 
!***  Called by: WRITE_FD
!***      Calls: (conditional) WRSLVS_FD
 
!  Writes frequency-domain output in profile form
!  if PRFL = 1  or in spectral output if PRFL = 2
 
!         TXD(JF) - transmitter orientation at frequency JF (in degrees)
!      YTR(JF,JS) - field at station JS for frequency JF.
 
!    All other variables defined in SUBROUTINE WRITE_FD
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80)
      INTEGER nw , nfrq , prfl , cmp , nstat , txd(nfrq) , krx , kry ,  &
     &        krz , js
      REAL mpz1(nstat) , freq(nfrq) , ytr(nfrq,nstat)
      REAL(KIND=Ql) mxd(nstat) , myd(nstat)
      CHARACTER(LEN=3) config(nfrq)
 
      IF ( prfl==2 ) THEN
         CALL wrslvs_fd(nw,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,ytr)
      ELSE
         IF ( cmp/=1 ) THEN
            WRITE (nw,99001) (freq(1:nfrq))
99001       FORMAT (/T10,'EAST     NORTH    ALT',30G14.4)
         ELSEIF ( maxval(freq)<100000. ) THEN
            WRITE (nw,99002) (freq(1:nfrq))
 
99002       FORMAT (/T10,'EAST     NORTH    ALT',F10.0,40F12.0)
         ELSE
            WRITE (nw,99003) (freq(1:nfrq))
99003       FORMAT (/T10,'EAST     NORTH    ALT',40G11.3)
         ENDIF
 
         IF ( cmp==1 ) WRITE (nw,99004) config(1:nfrq)
99004    FORMAT (T37,30(A,9X))
         WRITE (nw,'(3X)')
         DO js = 1 , nstat
            krx = nint(mxd(js))
            kry = nint(myd(js))
            krz = nint(mpz1(js))
            IF ( cmp==1 ) THEN
               WRITE (nw,99005) js , kry , krx , krz , ytr(1:nfrq,js)
99005          FORMAT (I3,2I10,I7,40F12.2)
            ELSE
               WRITE (nw,99006) js , kry , krx , krz , ytr(1:nfrq,js)
99006          FORMAT (I3,2I10,I7,40G14.4)
            ENDIF
         ENDDO
 
      ENDIF
 
      END SUBROUTINE wrslv_fd
