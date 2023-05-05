!*==txcnvd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION txcnvd(mxcnv,t,ntypls,trp,ypls,nsx,swx,swy)
!----------------------------------------------------------
!
!***  Called by: FOLD_AND_CONVOLVE, TQSTRIP
!***      Calls: CUBINT, CUBSPL, CUBVAL, LINVAL, TXCMRG
 
!  Convolves impulse B (step dB/dt) earth response function (ERF) with the
!  specified derivative of the source waveform at NSX points to produce
!  the system dB/dt response of the earth.
!
!       MXCNV = NTYPLS + NSX
!           T - convolution time in sec measured from the beginning
!               of the source waveform.
!   TRP, YPLS - abscissa & ordinate values of earth response function to
!               be convolved.
!      NTYPLS - number of values in TRP and YPLS
!         SWX - abscissa of time values of source waveform in sec.
!         SWY - dI/dt values derived from receiver dB/dt.
!         NSX - number of points in SWX & in each waveform stored in SWY
!
!  Defining  T1 = MIN {T, signal length,}, the convolution is formally
!  computed as
!
!   TXCNVD (T) = INT (T0 -> T) { YPLS (tau) * SWY (T-tau)  d tau }
 
!  where T0 = MAX { TRP(1), T - SWX (NSX)}
!
!       ONTIME RESPONSE
!       ---------------
!  For response in the on-time period, ( T < signal length) a correction to
!  account for the response from 0 -> T0 is needed.  Analysis and subsequent
!  numerical experiments confirm that as T -> 0, step dB/dt -> A * T**(-1/2).
!  Thus ERFINT, the integral of YPLS from 0 to TRP(1), is simply
!  2 * TRP(1) * YPLS (TRP(1)) if TRP(1) is chosen sufficiently early.
!  The convolution correction factor is SWY(T) * ERFINT.
 
!  Alternatively, we can difference the step B field from 0 to TRP(1) which
!  is a lot easier since the step B field at T = 0 is simply the DC field due
!  to a transmitter image buried at z = ALT; i.e., the z+z' term.  In this case,
!  the bigger TRP(1) is, the more accurate the difference in B but this must be
!  sufficiently small so that the change in dI/dt is negligable.  Thus, TRP(1)
!  is chosen to be .1 microsecond.
 
      IMPLICIT NONE
      INTEGER mxcnv , ntypls , nsx , n1 , j1 , n2 , j2 , ncnv
      REAL t , tc , t0 , trp(ntypls) , ypls(4,ntypls) , swx(nsx) ,      &
     &     swy(nsx,3) , ycnv(4,mxcnv) , xcnv(mxcnv) , x1(mxcnv) ,       &
     &     y1(mxcnv) , x2(mxcnv) , y2(mxcnv) , cubval , cubint , linval
 
      INTENT (in)mxcnv , t , ntypls , trp , ypls , nsx , swx , swy
 
!  Set up X1,Y1, the N1 values of SWX, SWY * YPLS for signal ontime < T.
!  where X1, the conjugate signal time, contains T-SWX values.
!  Set up X2,Y2, the N2 values of TRP, YPLS * SWY for ERF points  <= T.
 
!  Return TXCNVD = 0 if N1 + N2 < 4 or if NCNV < 4
 
      txcnvd = 0.0
      n1 = 0
      DO j1 = nsx , 1 , -1
         tc = t - swx(j1)
         IF ( tc<0. ) CYCLE
         n1 = n1 + 1
         x1(n1) = tc
         y1(n1) = swy(j1,1)*cubval(trp,ypls,ntypls,tc)
      ENDDO
 
      t0 = t - swx(nsx)
      t0 = max(t0,trp(1))/1.0001
 
      n2 = 0
      DO j2 = 1 , ntypls
         IF ( (trp(j2)>t0) .AND. (trp(j2)<t) ) THEN
            n2 = n2 + 1
            x2(n2) = trp(j2)
            tc = t - trp(j2)
            y2(n2) = ypls(1,j2)*linval(nsx,swx,swy,tc,1)
         ENDIF
      ENDDO
 
!  Merge the two lists into XCNV, YCNV of length NCNV.
!  Then spline and integrate
 
!+++++++++++++++++++++++++++++++++
      IF ( n1+n2<4 ) RETURN
!+++++++++++++++++++++++++++++++++
 
      CALL txcmrg(mxcnv,x1,y1,n1,x2,y2,n2,xcnv,ycnv,ncnv)
 
!+++++++++++++++++++++++++++++++++
      IF ( ncnv<4 ) RETURN
!+++++++++++++++++++++++++++++++++
 
      CALL cubspl(xcnv,ycnv,ncnv,0,0)
      txcnvd = cubint(xcnv,ycnv,ncnv,t0,t)
 
      END FUNCTION txcnvd
