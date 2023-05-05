!*==write_xdata.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
      SUBROUTINE write_xdata
!  ----------------------
 
!***  Called by NLSQ2
 
! Writes Survey Data
 
      WRITE (ql0,99001) line(js)
 
99001 FORMAT (T2,I10,'_Survey_Data')
      READ (ql0,'(A)') ql1
      WRITE (nw1,99002)
99002 FORMAT (T1,'/')
      WRITE (nw1,99003) trim(adjustl(ql1)) , js
99003 FORMAT (T2,'Line ',A,4X,'Station',I4)
 
      IF ( tdfd==1 ) THEN
         WRITE (nw1,99004) js , syd(js) , sxd(js) , sz0(js) , txdeg(js) &
     &                     , ryd(js) , rxd(js) , rz0(js) ,              &
     &                     xdata(1:ndata) , verr_fill(1:ndata)
99004    FORMAT (T1,I5,2F12.1,F8.1,F6.1,2F12.1,F8.1,350G13.4:)
      ELSE
         WRITE (nw1,99005) js , syd(js) , sxd(js) , sz0(js) , ryd(js) , &
     &                     rxd(js) , rz0(js) , xdata(1:ndata) ,         &
     &                     verr_fill(1:ndata)
99005    FORMAT (T1,I5,2F12.1,F8.1,2F12.1,F8.1,300G13.4)
      ENDIF
 
      END SUBROUTINE write_xdata
