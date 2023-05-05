!*==wrslvs_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE wrslvs_fd(nw,mxd,myd,mpz1,nstat,nfrq,freq,txd,cmp,ytr)
!-------------------------------------------------------------------
 
!***  Called by: WRSLV_FD
 
!  Writes frequency-domain output in spectral form
 
!         TXD(JF) - transmitter orientation at frequency JF (in degrees)
!      YTR(JF,JS) - field at station JS for frequency JF.
 
!    All other variables defined in SUBROUTINE WRITE_FD
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80) , Ncol = 40
      INTEGER nw , nfrq , cmp , nstat , txd(nfrq) , nblks , j1 , j2 ,   &
     &        jb , jf
      REAL mpz1(nstat) , freq(nfrq) , ytr(nfrq,nstat)
      REAL(KIND=Ql) mxd(nstat) , myd(nstat)
 
      nblks = nstat/Ncol
      IF ( mod(nstat,Ncol)>0 ) nblks = nblks + 1
 
      DO j1 = 1 , nblks
         jf = j1*Ncol
         jb = jf - Ncol + 1
         jf = min(jf,nstat)
 
         WRITE (nw,'(/T14,A,F9.0,39F13.0)') 'Receiver  Z' , mpz1(jb:jf)
         WRITE (nw,'(T13,A,F9.0,39F13.0)') 'Positions  E' , myd(jb:jf)
         WRITE (nw,'(T24,A,F9.0,39F13.0)') 'N' , mxd(jb:jf)
         WRITE (nw,'(T7,A)') 'Freq      TXCLN'
         WRITE (nw,'(3X)')
 
         DO j2 = 1 , nfrq
            IF ( cmp==1 ) THEN
               WRITE (nw,'(I3,G13.5,I4,3X,40F13.2)') j2 , freq(j2) ,    &
     &                txd(j2) , ytr(j2,jb:jf)                                    ! PPM
            ELSE
               WRITE (nw,'(I3,G13.5,I4,3X,40G13.4)') j2 , freq(j2) ,    &
     &                txd(j2) , ytr(j2,jb:jf)                                    ! pT
            ENDIF
         ENDDO
         WRITE (nw,99001)
 
99001    FORMAT (85('-')/)
      ENDDO
 
      END SUBROUTINE wrslvs_fd
