!*==solve2.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE solve2(npar,rsvt,zero,vmat,r,sv,nsv,delpar,wsp,pdre)
!----------------------------------------------------------------
 
!  Calculates the new parameter changes by taking the product of
!  the transformed error vector times the damping factors (second order
!  Marquardt damping) times the inverse singular value matrix.
!  It also calculates the predicted error decrease.
 
!*** Called by: NLSQ2
!***     Calls: DPROD1
 
!       Input
!       -----
!     NPAR -  the number of unknown parameters
!     RSVT -  the relative singular value threshold for the current iteration.
!     ZERO -  rejection level for relative singular values; i.e., truncation
!             instead of damping for RSV less than zero.  in NLSQ2, ZERO is
!             set to the square of the minimum allowed value of RSVT.
!     VMAT -  is the (NPAR * NPAR) V matrix of parameter space eigenvectors
!             output by ESVD.
!        R -  the transformed error vector output by ESVD.  It is the product
!             of the transposed U matrix times the observed error vector.
!       SV -  the (ordered) array of NPAR singular values returned by ESVD.
!             SV(1) > SV(2) > SV(3) etc.
!
!       Output
!       ------
!
!     NSV -  the returned number of singular values used in the solution.
!  DELPAR -  the returned solution ( parameter changes)
!     WSP -  WSP (1:NPAR) contains the squared r vector
!             WSP (NPAR: 2*NPAR) contains the damping factors.
!    PDRE -  the returned predicted decrease in residual squared error
 
      INTEGER nsv , i , npar
      REAL dprod1 , eigpar(npar) , r(npar) , sv(npar) , delpar(npar) ,  &
     &     wsp(3*npar) , vmat(npar,npar) , q , dmpfac , pdre , rsvt ,   &
     &     zero
      EXTERNAL dprod1
 
      INTENT (in)npar , rsvt , zero , vmat , r , sv
      INTENT (out)nsv , delpar , wsp , pdre
 
      nsv = 0
      DO i = 1 , npar
         wsp(i) = 0.
         wsp(npar+i) = 0.
         eigpar(i) = 0.
         q = sv(i)/sv(1)
         IF ( q<=zero ) CYCLE
         nsv = nsv + 1
         IF ( q<rsvt ) THEN
            q = (q/rsvt)**4
            dmpfac = q/(1.0+q)
         ELSE
            dmpfac = 1.0/(1.0+(rsvt/q)**4)
         ENDIF
 
!  Eigenparameter calculation.  store the damping factors in WSP(NPAR+I)
 
         eigpar(i) = r(i)*(dmpfac/sv(i))
         wsp(npar+i) = dmpfac
      ENDDO
 
!  Calculate change in physical parameters from eigenparameters.
 
      DO i = 1 , npar
         delpar(i) = dprod1(npar,npar,1,0.0,vmat(i,1),eigpar)
         wsp(i) = r(i)*r(i)
      ENDDO
      pdre = dprod1(npar,1,1,0.0,wsp(1),wsp(npar+1))
 
      END SUBROUTINE solve2
