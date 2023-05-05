!*==cnvrt2_depth.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE cnvrt2_depth(nlyr,thk,depth)
!----------------------------------------
 
!*** Called by: WRITE_MODEL, MAIN
 
!   NLYR - number of layers
!    THK - array of layer thicknessess (input)
!  DEPTH  THK - array of layer depths  (output)
 
      INTEGER j1 , nlyr
      REAL , DIMENSION(nlyr-1) :: thk , depth
 
      depth(1) = thk(1)
      DO j1 = 2 , nlyr - 1
         depth(j1) = depth(j1-1) + thk(j1)
      ENDDO
 
      END SUBROUTINE cnvrt2_depth
