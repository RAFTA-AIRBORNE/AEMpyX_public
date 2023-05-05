!*==wrslv.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE wrslv(nw,prfl,rz,rxd,ryd,nstat,sz,sxd,syd,nchnl,tms,   &
     &                 ytr)
!--------------------------------------------------------------------
 
!***  Called by: WRITE_TD
!***      Calls: WRSLVP
 
!  Writes time-domain output in temporal form for receiver
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80) , Ncol = 30
      INTEGER nw , prfl , nstat , nchnl , nblks , j1 , jb , jf , jt
      REAL rz(nstat,1) , sz(nstat) , tms(nchnl) , ytr(nchnl,nstat)
      REAL(KIND=Ql) sxd(nstat) , syd(nstat) , rxd(nstat,1) ,            &
     &              ryd(nstat,1)
 
      IF ( prfl==1 ) THEN
         CALL wrslvp(nw,sxd,syd,sz,nstat,nchnl,tms,ytr)
      ELSE
 
         nblks = nstat/Ncol
         IF ( mod(nstat,Ncol)>0 ) nblks = nblks + 1
 
         DO j1 = 1 , nblks
            jf = j1*Ncol
            jb = jf - Ncol + 1
            jf = min(jf,nstat)
 
            WRITE (nw,'(/T15,A,F10.0,29F13.0)') 'Z' , sz(jb:jf)
            WRITE (nw,'(T2,A,F10.0,29F13.0)') 'Transmitter  E' ,        &
     &             syd(jb:jf)
            WRITE (nw,'(T2,A,F10.0,29F13.0)') 'Positions    N' ,        &
     &             sxd(jb:jf)
            WRITE (nw,'(/T2,A,F10.0,29F13.0)') 'Rx Positions Z' ,       &
     &             rz(jb:jf,1)
            WRITE (nw,'(T15,A,F10.0,29F13.0)') 'E' , ryd(jb:jf,1)
            WRITE (nw,'(T2,A,F10.0,29F13.0)') 'CHNL   TIME  N' ,        &
     &             rxd(jb:jf,1)
            WRITE (nw,'(3X)')
 
            DO jt = 1 , nchnl
               WRITE (nw,'(I3,F9.3,3X,30G13.4)') jt , tms(jt) ,         &
     &                ytr(jt,jb:jf)
            ENDDO
            WRITE (nw,99001)
 
99001       FORMAT (85('-')/)
         ENDDO
      ENDIF
 
      END SUBROUTINE wrslv
