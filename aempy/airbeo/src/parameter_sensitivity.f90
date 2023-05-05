!*==parameter_sensitivity.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE parameter_sensitivity(npar,vmat,sv,bnd,import)
!----------------------------------------------------------
 
!***  Called by NLSQ2
 
!  Compute IMPORTANCE = parameter sensitivity for each parameter.  This is simply
!  an RMS rotation of the damping factors rotated into physical parameter space
 
!***  Called by: NLSQ2
 
!     NPAR - number of parameters
!     VMAT - the NPAR * NPAR matrix from SVD of Jacobian
!       SV - ordered singular values
!      BND - set minimum singular value threshold
!   IMPORT - importance  (damping factors rotated into physical parameter space.
 
      IMPLICIT NONE
      REAL , PARAMETER :: Eta = 1.E-7
      INTEGER npar , j1 , j2
      REAL vmat(npar,npar) , bnd , svnr , svnr4 , cum , ep4 , epsq
      REAL , DIMENSION(npar) :: sv , import , dmpfac
 
!  Normalise the singular values and set the damping factors for the
!  second order Marquardt method.
!  The threshold level is set at that used in NLSQ2
 
      dmpfac = 0.0
      epsq = bnd*bnd
      epsq = max(epsq,Eta)
      ep4 = epsq*epsq
 
      DO j1 = 1 , npar
         svnr = sv(j1)/sv(1)
         svnr = max(svnr,Eta)
         svnr4 = svnr**4
         IF ( svnr>epsq ) dmpfac(j1) = svnr4/(svnr4+ep4)
      ENDDO
 
      DO j1 = 1 , npar
         cum = 0.0
         DO j2 = 1 , npar
            cum = cum + (vmat(j1,j2)*dmpfac(j2))**2
         ENDDO
         import(j1) = sqrt(cum)
      ENDDO
 
      END SUBROUTINE parameter_sensitivity
