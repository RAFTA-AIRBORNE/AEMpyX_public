!*==cubint.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION cubint(xknot,coef,knot,x1,x2)
! ------------------------------------------------
!
!***  Called by:  EGT_BOSS TXCNVD, TXCNVL
!***      Calls: INTERV.  On exit from INTERV
 
!  Integrates a function from X1 to X2 using its cubic spline representation.
 
!       MFLAG = -1  => X is to the left of interpolated range
!             =  1  => X is to the right of interpolated range
!             =  0  => X is in the interpolated range
 
!      KNOT - total number of knots including endpoints.
!
!     XKNOT(I), I = 1,KNOT - Location of the knots.  The rightmost data
!                            point used to calculate coefficients is not
!                            included.
!
!     COEF(J,I), J = 1,4; I = 1,KNOT
!
!              The coefficients of the cubic spline represent the
!              indefinite integral of F, on the I'th interval, as:
!
!       INTGR [ F(X) ] = COEF(4,I)/24 * H**4  +  COEF(3,I)/6 * H**3  +
!                        COEF(2,I)/2 * H**2  +  COEF(1,I) * H
!
!                          WITH  H = X - XKNOT(K)
!
!  This is a modification of the FUNCTION PPVALU in the book
!  "A PRACTICAL GUIDE TO SPLINES"  by C. DE BOOR
 
!*********************************************************************
 
      IMPLICIT NONE
      INTEGER i , i1 , i2 , mflag , knot
      REAL h , h1 , h2 , x1 , x2 , xknot(knot) , coef(4,knot)
 
!  Find the indices I1 and I2 of largest breakpoints to the left of X1
!  and X2 respectively.
!
      CALL interv(xknot,knot-1,x1,i1,mflag)
      CALL interv(xknot,knot-1,x2,i2,mflag)
      h1 = x1 - xknot(i1)
      IF ( mflag==-1 ) h1 = 0.
 
      h2 = x2 - xknot(i2)
      cubint = (((coef(4,i2)*h2/4.0+coef(3,i2))*h2/3.0+coef(2,i2))      &
     &         *h2/2.0+coef(1,i2))                                      &
     &         *h2 - (((coef(4,i1)*h1/4.0+coef(3,i1))*h1/3.0+coef(2,i1))&
     &         *h1/2.0+coef(1,i1))*h1
 
!  Include integrals over intervening intervals.
 
      IF ( i2>i1 ) THEN
         DO i = i1 , i2 - 1
            h = xknot(i+1) - xknot(i)
            cubint = cubint +                                           &
     &               (((coef(4,i)*h/4.0+coef(3,i))*h/3.0+coef(2,i))     &
     &               *h/2.0+coef(1,i))*h
         ENDDO
      ENDIF
 
      END FUNCTION cubint
