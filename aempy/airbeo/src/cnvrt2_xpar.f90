!*==cnvrt2_xpar.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE cnvrt2_xpar(npar,nlyr,res,thk,xpar)
!-----------------------------------------------
 
!*** Called by: NLSQ2
 
!  Converts from model parameters to inversion parameters
 
!      NLYR = number of layers = 1 or 2
!      NPAR =  2*NLYR-1
!       RES = layer resistivities
!       THK - layer thicknesses
!      XPAR = transformed parameters
 
      INTEGER nlyr , npar
      REAL xpar(npar) , res(nlyr) , thk(nlyr-1)
 
      xpar(1:nlyr) = log(res(1:nlyr))
      xpar(nlyr+1:npar) = log(thk(1:nlyr-1))
 
      END SUBROUTINE cnvrt2_xpar
