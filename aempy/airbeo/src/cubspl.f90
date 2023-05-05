!*==cubspl.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE cubspl(xnot,c,n,ibcbeg,ibcend)
! ----------------------------------------------
 
!***  Called by: EGT_CSPL, FOLD_AND_CONVOLVE, HSBOSS_TD, INTER_EGT_CSPL, MGTBS,
!                PRM_BOSS, TDEM_3D, TQSTRIP, TXCNVD,
!
 
!  Calculates coefficients for cubic spline interpolation.
!  Call function CUBVAL to evaluate function values after interpolation.
!  From  * A PRACTICAL GUIDE TO SPLINES *  by Carl de Boor.
 
!             INPUT
!             -----
!
!     N = number of data points. assumed to be > 1.
!
!  (XNOT(I), C(1,I), I=1,...,N) = abscissae and ordinates of the data points.
!                                 XNOT is assumed to be strictly increasing.
!
!     IBCBEG, IBCEND = boundary condition indicators, and
!     C(2,1), C(2,N) = boundary condition information. Specifically,
!
!     IBCBEG = 0  No boundary condition at XNOT(1) is given.  In this case,
!                 the not-a-knot condition is used, i.e. the jump in the
!                 third derivative across XNOT(2) is forced to zero.  Thus
!                 first and the second cubic polynomial pieces are made to
!                 coincide.
!     IBCBEG = 1  the slope at XNOT(1) is made to equal C(2,1),
!                 supplied by input.
!     IBCBEG = 2  the second derivative at XNOT(1) is made to equal C(2,1),
!                 supplied by input.
!
!     IBCEND = 0, 1, or 2 has analogous meaning concerning the boundary
!                 condition at XNOT(n), with the additional information
!                 taken from C(2,n).
!
!          OUTPUT
!          ------
!
!     C(J,I), J=1,...,4; I=1,...,L (= N-1) = the polynomial coefficients
!         of the cubic interpolating spline with interior knots (or joints)
!         XNOT(2), ..., XNOT(N-1).
!
!        In the interval: (XNOT(I) - XNOT(I+1)), the spline F is given by:
!
!        F(X) = C(1,I) + H* (C(2,I) + H* (C(3,I) + H* C(4,I)/3.) /2.)
!
!     where H = X - XNOT(I).  FUNCTION  *CUBVAL* may be
!     used to evaluate F or its derivatives from XNOT,C, L = N-1,
!     AND K=4.
!******************************************************************************
!******************************************************************************
      IMPLICIT NONE
      INTEGER ibcbeg , ibcend , n , i , j , l , m
      REAL c(4,n) , xnot(n) , divdf1 , divdf3 , dxnot , g
 
      INTENT (in)xnot , n , ibcbeg , ibcend
      INTENT (inout)c
      SAVE 
 
!  A tridiagonal linear system for the unknown slopes S(I) of F at
!  XNOT(I), I=1,...,N, is generated and then solved by Gauss elimination,
!  with S(I) ending up in C(2,I), ALL I.
!  C(3,.) AND C(4,.) are used initially for temporary storage.
 
!  Compute first differences of XNOT sequence and store in C(3,.).
!  Also, compute first divided difference of data and store in C(4,.).
 
      l = n - 1
      DO m = 2 , n
         c(3,m) = xnot(m) - xnot(m-1)
         c(4,m) = (c(1,m)-c(1,m-1))/c(3,m)
      ENDDO
 
!  Construct first equation from the boundary condition, of the form
!  C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)
 
      IF ( ibcbeg<1 ) THEN
         IF ( n>2 ) THEN
 
!  Not-a-knot condition at left end and N > 2.
 
            c(4,1) = c(3,3)
            c(3,1) = c(3,2) + c(3,3)
            c(2,1) = ((c(3,2)+2.*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))&
     &               /c(3,1)
            GOTO 100
         ELSE
 
!  No condition at left end and N = 2.
 
            c(4,1) = 1.
            c(3,1) = 1.
            c(2,1) = 2.*c(4,2)
            GOTO 400
         ENDIF
      ELSEIF ( ibcbeg==1 ) THEN
 
!  Slope prescribed at left end.
 
         c(4,1) = 1.
         c(3,1) = 0.
      ELSE
 
!  Second derivative prescribed at left end.
 
         c(4,1) = 2.
         c(3,1) = 1.
         c(2,1) = 3.*c(4,2) - c(3,2)*c(2,1)/2.
      ENDIF
      IF ( n==2 ) GOTO 400
 
!  if there are interior knots, generate the corresponding equations and
!  perform the forward pass of Gauss elimination, after which the M-TH
!  equation reads    C(4,M)*S(M) + C(3,M)*S(M+1) = C(2,M).
 
 100  DO m = 2 , l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
         c(4,m) = g*c(3,m-1) + 2.*(c(3,m)+c(3,m+1))
      ENDDO
 
!  Construct last equation from the second boundary condition, of the form
!  (-G*C(4,N-1))*S(N-1) + C(4,N)*S(N) = C(2,N)
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since C array happens to be set up just right for it
!  at this point.
 
      IF ( ibcend<1 ) THEN
         IF ( n/=3 .OR. ibcbeg/=0 ) THEN
 
!  Not-a-knot and N > 2, and either N > 3 or also not-a-knot at
!  left end point.
 
            g = c(3,n-1) + c(3,n)
            c(2,n) = ((c(3,n)+2.*g)*c(4,n)*c(3,n-1)+c(3,n)**2*(c(1,n-1)-&
     &               c(1,n-2))/c(3,n-1))/g
            g = -g/c(4,n-1)
            c(4,n) = c(3,n-1)
            GOTO 500
         ENDIF
      ELSEIF ( ibcend==1 ) THEN
         GOTO 600
      ELSE
         GOTO 300
      ENDIF
 
!  Either (N=3 and not-a-knot also at left) or (N=2 and not not-a-
!  knot at left end point).
 
 200  c(2,n) = 2.*c(4,n)
      c(4,n) = 1.
      g = -1./c(4,n-1)
      GOTO 500
 
!  Second derivative prescribed at right endpoint.
 
 300  c(2,n) = 3.*c(4,n) + c(3,n)*c(2,n)/2.
      c(4,n) = 2.
      g = -1./c(4,n-1)
      GOTO 500
 400  IF ( ibcend<1 ) THEN
         IF ( ibcbeg>0 ) GOTO 200
 
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
 
         c(2,n) = c(4,n)
         GOTO 600
      ELSEIF ( ibcend==1 ) THEN
         GOTO 600
      ELSE
         GOTO 300
      ENDIF
 
!  Complete forward pass of Gauss elimination.
 
 500  c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
 
!  Perform back substitution.
 
 600  j = l
      DO
         c(2,j) = (c(2,j)-c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         IF ( j<=0 ) THEN
 
!  Generate cubic coefficients in each interval, i.e., the derivatives at its
!  left endpoint, from value and slope at its endpoints.
 
            DO i = 2 , n
               dxnot = c(3,i)
               divdf1 = (c(1,i)-c(1,i-1))/dxnot
               divdf3 = c(2,i-1) + c(2,i) - 2.*divdf1
               c(3,i-1) = 2.*(divdf1-c(2,i-1)-divdf3)/dxnot
               c(4,i-1) = (divdf3/dxnot)*(6./dxnot)
            ENDDO
            EXIT
         ENDIF
      ENDDO
      END SUBROUTINE cubspl
