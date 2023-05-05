!*==tqstrip.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE tqstrip(ider,ntypls,trp,ypls,nsx,swx,swy,ycnv)
!----------------------------------------------------------
!
!***  Called by: FOLD_AND_CONVOLVE
!***      Calls: CUBSPL, TXCNVD, TXCNVL
 
!  A stripped down version of TXCNVD is used to convolve the earth response
!  with the receiver waveform for for every point on that waveform.
!  The Geotem / Questem correlation is used to strip remnant primary field.
!  The result is sent back for binning into NCHNL receiver windows.
!
!   TRP, YPLS - abscissa & ordinate values of earth response function to
!               be convolved.
!        IDER - derivative indicator
!      NTYPLS - number of values in TRP and YPLS
!         SWX - abscissa of time values of source waveform in sec.
!         SWY - dI/dt values derived from receiver dB/dt. + raw waveform.
!         NSX - number of points in SWX & in each waveform stored in SWY
!        YCNV - the stripped convolved waveform
!
!  Defining  T1 = MIN {T, signal length,}, the convolution is formally
!  computed as
!
!   TXCNVD (T) = INT (T0 -> T) { YPLS (tau) * SWY (T-tau)  d tau }
 
!  where T0 = MAX { TRP(1), T - SWX (NSX)}
!
 
      IMPLICIT NONE
      REAL , PARAMETER :: T0_min = 1.E-7
      INTEGER ider , ntypls , nsx , jt , mxcnv
      REAL t , trp(ntypls) , ypls(4,ntypls) , swx(nsx) , swy(nsx,3) ,   &
     &     ycnv(4,nsx) , txcnvl , txcnvd , a1 , b1 , alpha
 
      INTENT (in)ider , ntypls , trp , ypls , nsx , swx , swy
      INTENT (out)ycnv
 
      a1 = 0.
      b1 = 0.
      ycnv = 0.
      mxcnv = ntypls + nsx
 
      DO jt = 2 , nsx       !  Convolve NSW points using the derived waveform
         t = swx(jt)
         IF ( t<T0_min ) CYCLE
         IF ( ider==0 ) THEN
                           ! Waveform input as I or B field (derived dI/dt)
            ycnv(1,jt) = txcnvl(t,ntypls,trp,ypls,nsx,swx,swy)
         ELSE              ! Waveform input as voltage (known dI/dt)
            ycnv(1,jt) = txcnvd(mxcnv,t,ntypls,trp,ypls,nsx,swx,swy)
         ENDIF
 
         a1 = a1 + ycnv(1,jt)*swy(jt,3)
                                      !  Compute correlation
         b1 = b1 + swy(jt,3)*swy(jt,3)
      ENDDO
 
      alpha = a1/b1
      DO jt = 1 , nsx
         ycnv(1,jt) = ycnv(1,jt) - alpha*swy(jt,3)
      ENDDO
 
      CALL cubspl(swx,ycnv,nsx,0,0)
 
      END SUBROUTINE tqstrip
