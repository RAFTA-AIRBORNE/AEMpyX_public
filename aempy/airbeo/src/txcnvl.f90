!*==txcnvl.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      REAL FUNCTION txcnvl(t,ntypls,trp,ypls,nsx,swx,swy)
!----------------------------------------------------
 
!***  Called by: FOLD_AND_CONVOLVE, TQSTRIP
!***      Calls: CUBINT, CUBVAL
 
!  Computes the system dB/dt response by convolving the computed dI/dt with
!  the impulse B response of the earth.  For step current drops, system dB/dt
!  is computed asthe product of instantaneous current drop times the
!  earth step dB/dt.
 
!  This routine assumes that the source waveform is composed of NSX linear
!  segments.  Thus NSX-1 constant dI/dt values are contained in SWY(*,1).
 
!  The input earth response function (step dB/dt or equivalently, impulse B)
!  must be contained in a splined array of NTYPLS values of time (abscissa) TRP
!  and ordinate YPLS.  System dB/dt is computed by integrating YPLS between
!  the SWX points of constant dI/dt segments.
 
!              T - convolution time in sec measured from the beginning
!                  of the source waveform.
!      TRP, YPLS - abscissa & ordinate values of earth response function to
!                  be convolved.
!         NTYPLS - number of values in TRP and YPLS
!            SWX - abscissa of time values of source waveform in sec.
!       SWY(*,1) - dI/dt if it exists (0 otherwise)
!       SWY(*,2) - first difference values of source waveform
!                  (-delta I) in amps.
!            NSX - number of points in SWX & WAVEFORM
 
      IMPLICIT NONE
      REAL , PARAMETER :: T0_min = 1.E-7
      INTEGER ntypls , nsx , jt
      REAL t , tf , cnv , tb , delt , seg , trp(ntypls) , ypls(4,ntypls)&
     &     , swx(nsx) , swy(nsx,3) , tend , cubint , cubval
      LOGICAL der
 
      tf = t - trp(1)
      cnv = 0.
      DO jt = 2 , nsx
         IF ( swx(jt-1)>tf ) EXIT
         tb = t - min(tf,swx(jt))
         delt = swx(jt) - swx(jt-1)
         der = .FALSE.
         IF ( delt>T0_min ) THEN
            tend = t - swx(jt-1)
            der = .TRUE.
         ENDIF
 
!  For an instantaneous step drop in current, SEG is YPLS times SWY(*,2),
!  since YPLS is already the dB/dt step response.  Otherwise SEG is the
!  integral of YPLS * constant dI/dt SWY(*,1) since YPLS is also impulse B.
 
         IF ( der ) THEN
            seg = swy(jt,1)*cubint(trp,ypls,ntypls,tb,tend)
         ELSE
            seg = swy(jt,2)*cubval(trp,ypls,ntypls,tb)
         ENDIF
         cnv = cnv + seg
      ENDDO
      txcnvl = cnv
 
      END FUNCTION txcnvl
