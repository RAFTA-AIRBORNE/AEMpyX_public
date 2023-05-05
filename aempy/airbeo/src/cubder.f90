!*==cubder.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION cubder(xknot,coef,knot,x1)
! --------------------------------------------
 
!***  Called by: FOLD_AND_CONVOLVE
!***      Calls: INTERV.  On exit from INTERV
 
!  Evaluates the first derivative of a function from its cubic spline
!  interpolation.
 
!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range
 
!         KNOT - total number of knots including endpoints.
!     XKNOT(I), I = 1,KNOT - location of the knots.  The rightmost data
!                            point used to calculate coefficients is not
!                            used.
!     COEF(J,I), J = 1,4; I = 1,KNOT = Jth derivative at H = 0
!                                      where  H = X - XKNOT(I)
!******************************************************************************
!******************************************************************************
 
      IMPLICIT NONE
      INTEGER i , mflag , knot
      REAL xknot(knot) , coef(4,knot) , x1 , h
 
!  Find index i of largest breakpoint to the left of X1.
 
      CALL interv(xknot,knot-1,x1,i,mflag)
      h = x1 - xknot(i)
      IF ( mflag==-1 ) h = 0.
 
      cubder = (coef(4,i)*h/2.+coef(3,i))*h + coef(2,i)
 
      END FUNCTION cubder
