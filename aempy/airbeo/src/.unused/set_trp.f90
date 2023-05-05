!*==set_trp.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
      SUBROUTINE set_trp
!  ------------------
 
!***  Called by: MAIN
 
!  Sets up interpolation times for FD -> TD transform which use the
!  exact 6 points per decade frequency-domain data plus 6 per decade
!  interpolated values.  These are based on a 12 point per decade
!  cosine filter derived from the Niels Christensen routine FILCOA
!  with OMEGA = .3 PI and shift 0.
 
!             OUTPUT
!             ------
 
!        TRP - array of time values for FD -> TD transformations
!      NTYRP - number of values in TRP
!     EXTENT - the latest time for which time-domain output is required.
!      PULSE - time length of one signal pulse
!     NTYPLS - number of TRP values in 1 PULSE
 
 
      REAL , PARAMETER :: Twopi = 6.2831853
      INTEGER mxtym , j1
      REAL t0 , extent
      REAL , ALLOCATABLE :: qqq(:)
      REAL(KIND=ql) tbase , qtym , tq
 
      mxtym = 200
      ALLOCATE (qqq(mxtym))
      qqq = 0.
 
      qtym = log(10.D0)/12.D0
      qtym = exp(qtym)
      extent = 2.0*npuls*pulse
 
 
      t0 = minval(topn) - swx(nsx)
      t0 = max(t0,t0_min)
      tbase = 1.D0/dble(Twopi)
      DO j1 = 1 , mxtym
         IF ( tbase<t0 ) EXIT
         tbase = tbase/qtym
      ENDDO
 
      tq = tbase
      qqq(1) = real(tq)
      DO j1 = 2 , mxtym
         ntyrp = j1
         tq = tq*qtym
         qqq(j1) = real(tq)
         IF ( qqq(j1)<pulse ) ntypls = j1 + 2
         IF ( qqq(j1)>extent ) EXIT
      ENDDO
 
      ALLOCATE (trp(ntyrp))
      trp(1:ntyrp) = qqq(1:ntyrp)
      DEALLOCATE (qqq)
 
      END SUBROUTINE set_trp
