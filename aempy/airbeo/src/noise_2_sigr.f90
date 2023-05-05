!*==noise_2_sigr.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE noise_2_sigr(npar,ndata,xmodl,xdata,xwts,nsr)
!--------------------------------------------------------
 
!***  Called by NLSQ2
 
!  Calculates and prints the noise to signal ratio.
 
!***  Called by: MAIN
 
!      NPAR - the number of free variables
!     NDATA - the number of data values
!     XMODL - model data in microvolts
!     XDATA - data to be inverted in microvolts
!      XWTS - weights for data points (intger 0 or 1)
!       NSR - noise to signal ratio
 
      IMPLICIT NONE
      INTEGER npar , ndata , xwts(ndata) , j1
      REAL , PARAMETER :: Tol = .1E-5
      REAL xmodl(ndata) , xdata(ndata) , ym(ndata) , yd(ndata) , pk ,   &
     &     pkm , base , zm , zd , ybar , nsr , cumm , cumd , cwmd
 
      INTENT (in)npar , ndata , xmodl , xdata , xwts
      INTENT (out)nsr
 
 
!  Set up stitched log representation for model and data voltages.
!  Use a base of 6 orders of magnitude or the data dynamic range,
!  whichever is the lesser.
!  Accumulate means and variances in double precision.
 
      base = abs(xdata(1))
      pk = base
      DO j1 = 2 , ndata
         base = min(abs(xdata(j1)),base)
         pk = max(abs(xdata(j1)),pk)
      ENDDO
      pkm = Tol*pk
      base = max(pkm,base)
 
      cumm = 0.
      ym = 0.
      yd = 0.
      DO j1 = 1 , ndata
         IF ( xwts(j1)>0 ) THEN
            zm = abs(xmodl(j1))
            zd = abs(xdata(j1))
            IF ( zm>base ) THEN
               ym(j1) = log(zm/base)
               IF ( xmodl(j1)<0 ) ym(j1) = -ym(j1)
            ENDIF
            IF ( zd>base ) THEN
               yd(j1) = log(zd/base)
               IF ( xdata(j1)<0 ) yd(j1) = -yd(j1)
            ENDIF
            cumm = cumm + ym(j1)
         ENDIF
      ENDDO
 
      ybar = cumm/ndata
 
      cumm = 0.
      cumd = 0.
      DO j1 = 1 , ndata
         cumm = cumm + (ym(j1)-ybar)**2
         cwmd = xwts(j1)*(ym(j1)-yd(j1))
         cumd = cumd + cwmd**2
      ENDDO
 
      cumm = cumm/(ndata-1)
      cumd = cumd/(ndata-npar)
 
      nsr = sqrt(cumd/cumm)
 
      END SUBROUTINE noise_2_sigr
