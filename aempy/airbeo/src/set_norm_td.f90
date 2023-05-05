!*==set_norm_td.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE set_norm_td(nw,bunit,bffac,kppm,punit,ppfac,prm_td,    &
     &                       norm)
!--------------------------------------------------------------------
 
!*** Called by: MAIN
!
!   Airbeo computes all fields in nT or nT/s.  In order to match field data for
!   inversion, it computes factors (NORM) to convert computed data into pT, pT/s,
!   fT, fT/s, pct, ppt, ppm, ppb as required.
!
!             Input
!             -----
!      BUNIT, BFAC, PUNIT & PPFAC are set in SUBROUTINE READ_SYSTEM_SETUP
!      BUNIT can have values nT/s, nT, pT, pT/s, fT or fT/s
!      BFFAC is the conversion factor needed to achieve this from nT or nT/s
!            = 1, 1000, or 1E6
!
!
!
!      KPPM = 0   => No PPM normalisation (automatic if ISW = 4)
!      KPPM = 1   => All components are normalised to in-line primary field
!      KPPM = 3   => All components are normalised to vertical primary field
!      KPPM = 123 => Vertical component normalised to the vertical primary &
!      KPPM = 4   => All components are normalised to total primary field
!
!      For KPPM > 0:
!        PUNIT can have values pct, ppt, ppm or ppb; eg,
!           parts per hundred (percent)
!           parts per thousand
!           parts per million
!           parts per billion
!
!      PPFAC = 1e2, 1e3, 1e6 or 1e9 is the conversion factor to achieve this
!
!   PRM_TD(I) - peak primary dB/dt in nT/s if STEP = 0 or
!               peak primary B (in nT if STEP = 1)
!               I = 1, 2 & 3 are the in-line, transverse & vertical components.
!
!             Output
!             ------
!
!       NORM - Factor to convert time-domain response in nT or nT/s into relevant units.
 
      IMPLICIT NONE
      INTEGER nw , kppm , j
      REAL prm_td(3) , norm(3) , ptd(4) , ppfac , bffac , xq(4)
      LOGICAL no(3)
      CHARACTER(LEN=4) bunit , punit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAMELIST /mynmlno/ nw,bunit,bffac,kppm,punit,ppfac,prm_td,norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        write(*,*) ' start set_norm_td'  
!        write(*,nml=mynmlno)
      
        
      no = .FALSE.
      ptd(1:3) = abs(prm_td(1:3))
      ptd(4) = sqrt(ptd(1)**2+ptd(2)**2+ptd(3)**2)
      xq(1:3) = bffac*prm_td(1:3)
      xq(4) = bffac*ptd(4)
      DO j = 1 , 3
         IF ( ptd(j)<1.E-3*ptd(4) ) THEN
            ptd(j) = ptd(4)
            no(j) = .TRUE.
         ENDIF
      ENDDO
      write(*,*) ' ptd',ptd
      IF ( kppm>0 ) WRITE (nw,99001) punit , bunit , xq(3) , xq(1) ,    &
     &                               xq(2) , xq(4)
99001 FORMAT (//T31,'Normalisation in ',A/T31,                          &
     &        'Primary field units are ',A//T20,'  Vertical primary =', &
     &        G12.4/T20,'   In-line primary =',G12.4/T20,               &
     &        'Transverse primary =',G12.4/T20,'     Total primary =',  &
     &        G12.4)
 
!  Primary fields, PRM_TD, are computed in nT, NORM, the normalisation.  If the
!  field data are not in nT, NORM, must be converted to pT or fT as required.
 
      IF ( kppm==0 ) THEN
         norm = bffac
      ELSEIF ( kppm==1 ) THEN
         norm = ppfac/ptd(1)
         IF ( no(1) ) THEN
            WRITE (nw,99005) punit , norm(1)
         ELSE
            WRITE (nw,99002) punit , norm(1)
 
99002       FORMAT (/T11,                                               &
     &       'Each component is normalised to the in-line primary field'&
     &       /T26,'In-line ',A,' norm =',G12.4)
         ENDIF
      ELSEIF ( kppm==3 ) THEN
         norm = ppfac/ptd(3)
         IF ( no(3) ) THEN
            WRITE (nw,99005) punit , norm(3)
         ELSE
            WRITE (nw,99003) punit , norm(3)
99003       FORMAT (/T11,                                               &
     &      'Each component is normalised to the vertical primary field'&
     &      /T25,'Vertical ',A,' norm =',G12.4)
         ENDIF
      ELSEIF ( kppm==4 ) THEN
         norm = ppfac/ptd(4)
         WRITE (nw,99005) bunit , norm(1)
      ELSEIF ( kppm==123 ) THEN
         norm(1:3) = ppfac/ptd(1:3)
         WRITE (nw,99004) norm(1:3)
99004    FORMAT (/T11,                                                  &
     & 'Each component is normalised to its corresponding primary field'&
     & /T23,'   In-line norm =',G12.4/T23,'Transverse norm =',G12.4/T23,&
     & '  Vertical norm =',G12.4)
      ENDIF
99005 FORMAT (/T3,T11,                                                  &
     &        'Each component is normalised to the total primary field'/&
     &        T28,'Total ',A,' norm =',G12.4)
        
 !       write(*,*) ' end of set_norm_td'  
 !       write(*,nml=mynmlno)

      END SUBROUTINE set_norm_td
