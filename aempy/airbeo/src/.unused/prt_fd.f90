!*==prt_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
      SUBROUTINE prt_fd
!  -----------------
 
!***  Called by: WRITE_MISFIT
 
      WRITE (nw,'(/40F11.0)') freq(1:nfrq)
      WRITE (nw,'(/40F11.2)') 100.*verr(n1:n2)
      WRITE (nw,'( 40F11.2)') xmodl(n1:n2)
      WRITE (nw,'( 40F11.2)') xdata(n1:n2)
      END SUBROUTINE prt_fd
