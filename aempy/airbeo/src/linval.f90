!*==linval.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION linval(nx,xval,yval,x1,ider)
!-------------------------------------------
 
!***  Called by: TXCNVD
!***      Calls: INTERV
 
!  Evaluates a function at X1 from from its linear representation.
!
!           On exit from INTERV
!
!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range
 
!
!     XVAL(1:NX) - location of the abscissa knots.  The rightmost data point
!                  used to calculate coefficients is not included.
!
!     YVAL(1:NX,1) = function values.
!     YVAL(1:NX,2)   may be populated but aren't used.
!
!     If IDER = 0, the value in the interval is that of the leftmost knot.
!                  because the derivative has been computed using two knot
!                  values and stored at the left node.
!
!     If IDER = 1, the value is a linear interpolation between the knots.
!
!             METHOD
!             ------
!
!  The interval index I, appropriate for X, is found through a call to INTERV.
!  The formula for F is evaluated using nested multiplication.
!******************************************************************************
 
      IMPLICIT NONE
      INTEGER i , mflag , nx , ider
      REAL xval(nx) , yval(nx,3) , x1 , h
 
      INTENT (in)nx , xval , yval , x1 , ider
!
!  Find index I of largest breakpoint to the left of X1.
!
      CALL interv(xval,nx-1,x1,i,mflag)
 
      IF ( ider==0 ) THEN !  Computed derivative values stored at right node (26.01.00)
         linval = yval(i+1,1)
      ELSE
         h = x1 - xval(i)
         IF ( mflag==-1 ) h = 0.
 
         linval = yval(i,1) + h*(yval(i+1,1)-yval(i,1))                 &
     &            /(xval(i+1)-xval(i))
      ENDIF
 
      END FUNCTION linval
