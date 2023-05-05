!*==set_survey.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
      SUBROUTINE set_survey
!  ---------------------
 
!***  Called by: MAIN
 
 
      IMPLICIT NONE
 
      IF ( baromtrc==1 ) THEN
         sz = sz - gnd_lvl    ! Change barometric altitude to ground clearance
         IF ( abs(gnd_lvl)>0.01 ) WRITE (nw,99001) gnd_lvl
 
 
99001    FORMAT (/T3,                                                   &
     &        'Barometric altitude will be changed to ground clearance.'&
     &        /T3,'The vertical origin is shifted down by',F6.1)
      ENDIF
 
!  Compute receiver coordinates in both body centred and real world systems.
 
      DO js = 1 , nstat
         csf = cos(fangle(js))
         snf = sin(fangle(js))
 
         IF ( tdfd==1 ) THEN
            sx(js) = real(sxd(js))
            sy(js) = real(syd(js))
            rx(js,1) = sx(js) - xrx(js)*csf + yrx(js)*snf
            ry(js,1) = sy(js) - yrx(js)*csf - xrx(js)*snf
            rz(js,1) = sz(js) - zrx(js)
         ELSE
            sx(js) = real(sxd(js)) + 0.5*(xrx(1)*csf-yrx(1)*snf)
            sy(js) = real(syd(js)) + 0.5*(yrx(1)*csf+xrx(1)*snf)
            DO jf = 1 , nfrq
               rx(js,jf) = sx(js) - xrx(jf)*csf + yrx(jf)*snf
               ry(js,jf) = sy(js) - yrx(jf)*csf - xrx(jf)*snf
               rz(js,jf) = sz(js) - zrx(jf)
            ENDDO
         ENDIF
      ENDDO
      rxd = real(rx,kind=ql)
      ryd = real(ry,kind=ql)
      sxd = real(sx,kind=ql)
                            !  SXD & SYD were midpoints in Freq domain
      syd = real(sy,kind=ql)
                            !  Now they are Tx positions based on XRX(1), YRX(1)
 
      END SUBROUTINE set_survey
