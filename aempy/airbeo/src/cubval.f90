!*==cubval.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION cubval(xknot,coef,knot,x1)
! --------------------------------------------
 
!***  Called by: COSTRN, EGT_BOSS, FOLD_AND_CONVOLVE, INTER_EGT_BOSS,
!                MGTV1, PRM_BOSS, TXCNVD, TXCNVL
!
!***      Calls: INTERV.
 
!  On exit from INTERV,
!  Evaluates a function at X1 from from its cubic spline representation.
 
!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range
 
!      KNOT - total number of knots including endpoints.
!
!     XKNOT(I), I = 1,KNOT - location of the knots.  The rightmost data
!                            point used to calculate coefficients is not
!                            included.
!
!     COEF(J,I), J = 1,4; I = 1,KNOT
!
! The coefficients of the cubic spline on the I'th interval represent F as:
!
!                F(X) = COEF(4,I)/6 * H**3  +  COEF(3,I)/2 * H**2  +
!                       COEF(2,I) * H  +  COEF(1,I)
!
!                          with  H = X - XKNOT(I)
!
!  This is a modification of the FUNCTION PPVALU in the book
!  "A PRACTICAL GUIDE TO SPLINES"  by C. DE Boor
!
!             METHOD
!             ------
!
!  The interval index I, appropriate for X, is found through a call to INTERV.
!  The formula for F is evaluated using nested multiplication.
!******************************************************************************
 
      IMPLICIT NONE
      INTEGER i , mflag , knot
      REAL xknot(knot) , coef(4,knot) , x1 , h
 
      INTENT (in)xknot , coef , knot , x1
!
!  Find index I of largest breakpoint to the left of X1.
!
      CALL interv(xknot,knot-1,x1,i,mflag)
      h = x1 - xknot(i)
      IF ( mflag==-1 ) h = 0.
 
      cubval = ((coef(4,i)*h/3.0+coef(3,i))*0.5*h+coef(2,i))            &
     &         *h + coef(1,i)
 
      END FUNCTION cubval
