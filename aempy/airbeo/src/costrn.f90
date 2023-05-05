!*==costrn.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION costrn(wf,yfrq,nfrq,kfrq,t)
!------------------------------------------
 
!***  Called by: HSBOSS_TD, TDEM_3D
!***      Calls: CUBVAL
 
! LAST MODIFICATION DATE: October, 2001
 
! Produces time-domain value at time T by cosine transformation of NFRQ
! frequency-domain values contained in cubic spline array YFRQ.
! KFRQ is the high frequency cutoff, less than or equal to NFRQ.
! Array WF contains the LOG (base e) of the angular frequency values.
 
! The routine uses filter coefficients derived from the Niels Christensen
! fast Hankel transform routine FILCOA at a spacing of 12 points per decade
! and omega = 0.3.  Various filters were tested using a vertical magnetic
! dipole receiver in a very large circular for which accurate frequency
! and time-domain solutions were programmed.  This particular filter gave
! the overall best accuracy for 1/2 spaces ranging in resistivity from
! .1 to 10,000 ohm-m for times ranging from .01 to 50 msec.
 
 
!  K(W,T) = (2/PI) * F(W) * COS(WT) dW
 
! Letting X = WT, the above becomes
!
!  K(W,T) = (2/PI*T) * F(X/T) * COS(X) dX
!
! From Abramowitz and Stegun, COS(X) = SQRT(X*PI/2) * J(-1/2:X).
! Filter Coefficients are used to represent X**(1/2) * J(-1/2:X)
!
!  COSTRN = SQRT (2/PI) * SUM(i) { WCOS(i) * F [X(i) /T] }
 
! The accumulation is done using 12 digit precision
 
 
      USE filt_coef_q
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ndec_cos = 12 , Kflow = -200 , Kfhigh = 99
      REAL , PARAMETER :: Fac = .7978846 , Tol = 1.0E-6
      INTEGER j1 , nfrq , kfrq
      REAL wf(nfrq) , yfrq(4,nfrq) , t , ys , cubval , v1
      REAL(KIND=Ql) delta , y1 , y , td , ytym , val
 
      INTENT (in)wf , yfrq , nfrq , t
 
 
      delta = log(10._QL)/real(Ndec_cos,kind=Ql)
      td = real(t,kind=Ql)
      ytym = 0.
      y1 = -log(td) - delcos
 
! Begin right side convolution at weight 0.
! Stop when frequency domain array is exhausted.
 
      MOVE_HIGH:DO j1 = 0 , Kfhigh
 
         y = y1 + j1*delta
         ys = real(y)
         IF ( ys>wf(kfrq) ) EXIT MOVE_HIGH
         IF ( ys<wf(1) ) ys = wf(1)
         v1 = cubval(wf,yfrq,nfrq,ys)
         val = wcos(j1)*real(v1,kind=Ql)
         ytym = ytym + val
      ENDDO MOVE_HIGH
 
      y = y1
 
! Begin left side convolution at weight -1.
! When log angular frequency is less than WF(3), check convergence.
! Continue left using the fact that impulse B is inversely proportional to
! frequency as freq -> 0; i.e., step response B is constant.
 
      MOVE_LOW:DO j1 = -1 , Kflow , -1
 
         y = y1 + j1*delta
         ys = real(y)
         IF ( ys>wf(kfrq) ) CYCLE MOVE_LOW
         IF ( ys<wf(1) ) ys = wf(1)
         v1 = cubval(wf,yfrq,nfrq,ys)
         val = wcos(j1)*real(v1,kind=Ql)
         ytym = ytym + val
         IF ( (y<wf(3)) ) THEN
            IF ( abs(val)<Tol*abs(ytym) ) EXIT MOVE_LOW
         ENDIF
      ENDDO MOVE_LOW
 
      costrn = Fac*real(ytym)/t
 
      END FUNCTION costrn
