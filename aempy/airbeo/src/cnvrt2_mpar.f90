!*==cnvrt2_mpar.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE cnvrt2_mpar(npar,nlyr,xpar,res,thk)
!-----------------------------------------------
 
!*** Called by: NLSQ2
 
!  Converts inversion parameters to model parameters
!
!      NLYR = number of layers = 1 or 2
!      NPAR =  2*NLYR-1
!       RES = layer resistivities
!       THK - layer thicknesses
!      XPAR = transformed parameters
 
      INTEGER npar , nlyr
      REAL xpar(npar) , res(nlyr) , thk(nlyr-1)
 
      res(1:nlyr) = exp(xpar(1:nlyr))
      thk(1:nlyr-1) = exp(xpar(nlyr+1:npar))
 
      END SUBROUTINE cnvrt2_mpar
