!*==dprod1.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION dprod1(n,n1,n2,a,b,c)
!------------------------------------
 
!***  Called by: ESVD, SOLVE2
 
!     Double precision inner product routine
!     DPROD = A + B * C
 
!         A = scalar
!       B,C = vectors (can be rows or columns of arrays)
!         N = length of vectors
!     N1,N2 = increment for b,c
!           = 1 if col of array
!           = col length (i.e. no. of rows) if row of array
 
!  DPROD must be declared external by any routine using it because there is
!  a standard intrinsic function with the same name.
!  If omitted compilation warnings result.
 
      IMPLICIT NONE
      INTEGER n , na , nb , n1 , n2 , i
      DOUBLE PRECISION z1 , z2 , z3
      REAL a , b(*) , c(*)
 
      INTENT (in)n , n1 , n2 , a , b , c
 
      z1 = a
      IF ( n>=1 ) THEN
         na = 1
         nb = 1
         DO i = 1 , n
            z2 = b(na)
            z3 = c(nb)
            z1 = z1 + z2*z3
            na = na + n1
            nb = nb + n2
         ENDDO
      ENDIF
      dprod1 = real(z1)
 
      END FUNCTION dprod1
