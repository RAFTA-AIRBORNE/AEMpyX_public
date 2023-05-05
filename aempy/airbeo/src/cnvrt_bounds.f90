!*==cnvrt_bounds.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
!==================================================================================
!  ROUTINES SPECIFIC FOR INVERSION
!  -------------------------------
 
      SUBROUTINE cnvrt_bounds(npar,lbnd,ubnd,cxpar,xlbnd,xubnd)
!----------------------------------------------------------
 
!*** Called by NLSQ2
 
! Converts the LBND & UBND specified by the user into the form used
! by the inversion, XLBND & XUBND.
 
! Bounds are only applied to layer parameters for CXPAR = 3
 
!  INPUT:  NPAR,LBND,UBND,CXPAR
! OUTPUT:  XLBND,XUBND
 
      INTEGER npar , cxpar(npar) , j1
      REAL , DIMENSION(npar) :: lbnd , ubnd , xlbnd , xubnd
 
      DO j1 = 1 , npar
         IF ( cxpar(j1)==3 ) THEN
            xlbnd(j1) = log(lbnd(j1))
            xubnd(j1) = log(ubnd(j1))
         ENDIF
      ENDDO
 
      END SUBROUTINE cnvrt_bounds
