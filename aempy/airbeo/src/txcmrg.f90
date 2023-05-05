!*==txcmrg.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE txcmrg(mxcnv,x1,y1,n1,x2,y2,n2,xcnv,ycnv,ncnv)
!----------------------------------------------------------
 
!***  Called by: TXCNVD
 
!  Merges two previously sorted list pairs X1, Y1 of length N1 and X2, Y2 of
!  length N2 into list pair XCNV, YCNV of length NCNV into ascending values of
!  XCNV.
 
      IMPLICIT NONE
      REAL , PARAMETER :: Tol = 1.E-3
      INTEGER mxcnv , n1 , n2 , ncnv , k1 , k2 , n , j1
      REAL delt , tl1 , xcnv(mxcnv) , x1(mxcnv) , y1(mxcnv) , x2(mxcnv) &
     &     , y2(mxcnv) , ycnv(4,mxcnv)
      LOGICAL list1 , list2
 
      INTENT (in)mxcnv , x1 , y1 , n1 , x2 , y2 , n2
      INTENT (out)xcnv , ycnv , ncnv
 
      list1 = .TRUE.
      list2 = .TRUE.
      k1 = 1
      k2 = 1
      n = n1 + n2
 
      DO j1 = 1 , n
         IF ( list1 .AND. list2 ) THEN
            IF ( x1(k1)<x2(k2) ) THEN
               xcnv(j1) = x1(k1)
               ycnv(1,j1) = y1(k1)
               k1 = k1 + 1
               IF ( k1>n1 ) list1 = .FALSE.
            ELSE
               xcnv(j1) = x2(k2)
               ycnv(1,j1) = y2(k2)
               k2 = k2 + 1
               IF ( k2>n2 ) list2 = .FALSE.
            ENDIF
         ELSEIF ( list1 ) THEN
            xcnv(j1) = x1(k1)
            ycnv(1,j1) = y1(k1)
            k1 = k1 + 1
            IF ( k1>n1 ) list1 = .FALSE.
         ELSEIF ( list2 ) THEN
            xcnv(j1) = x2(k2)
            ycnv(1,j1) = y2(k2)
            k2 = k2 + 1
            IF ( k2>n2 ) list2 = .FALSE.
         ENDIF
      ENDDO
 
      ncnv = 1 !  Clean up list
      DO j1 = 2 , n
         delt = xcnv(j1) - xcnv(ncnv)
         tl1 = Tol*xcnv(j1)
         IF ( delt>tl1 ) THEN
            ncnv = ncnv + 1
            xcnv(ncnv) = xcnv(j1)
            ycnv(1,ncnv) = ycnv(1,j1)
         ENDIF
      ENDDO
 
      END SUBROUTINE txcmrg
