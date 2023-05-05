!*==fold_and_convolve.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE fold_and_convolve(step,ider,nsx,swx,swy,npuls,pulse,   &
     &                             trp,ntypls,ntyrp,nchnl,topn,tcls,    &
     &                             ypls,gstrp,astrp,ycum)
!-------------------------------------------------------------------------------
 
!  Computes the "observed" response YCUM by convolving the splined earth
!  response function, YPLS, with the TX waveform.
 
!***  Called by: HSBOSS_TD, TDEM3D
!***      Calls: CUBDER, CUBVAL, CUBSPL, TXCNVD, TXCNVL, TQSTRIP
 
!     IDER - derivative indicator
!      NSX - number of points used to discretise transmitter signal
!      SWX - abscissae (seconds) of current waveform
!      SWY - dI/dt at times SWX
!    NPULS - number of bipolar pulses of length PULSE
!    PULSE - length single on pulse plus off-time
!    NTYRP - number of values in TRP for total signal length: 2 * NPULS *PULSE
!      TRP - array of time values for FD -> TD transformations
!   NTYPLS - number of TRP values in 1 PULSE
!    NCHNL - number of channels
!     TOPN - time at which receiver channel I opens.
!     TCLS - time at which receiver channel I closes.
!    GSTRP = 1 if Geotem / Questem stripping is to be applied
 
      IMPLICIT NONE
      INTEGER jt , ntyrp , ntypls , npuls , ider , step , nsx , nchnl , &
     &        jgl , jp , gstrp , astrp , mxcnv , ipl
      REAL pulse , trp(ntyrp) , swx(nsx) , swy(nsx,3) , topn(nchnl) ,   &
     &     tcls(nchnl) , t1 , t2 , width , tf , tfh , hwidth , yc1 ,    &
     &     tc(3) , glx(3) , glw(3) , ypls(4,ntyrp) , x , xp ,           &
     &     ycum(nchnl) , cubval , cubder , txcnvl , txcnvd , wt ,       &
     &     fold(ntypls) , ycnv(4,nsx)
      DATA glw(1:3)/.5555556 , .8888889 , .5555556/ , glx(1:3)          &
     &     / - .7745967 , 0. , .7745967/
 
      INTENT (in)ider , nsx , swx , swy , npuls , pulse , trp , ntypls ,&
     &        ntyrp , nchnl , topn , tcls , gstrp
      INTENT (inout)ypls
      INTENT (out)ycum
      
      
       NAMELIST /mynmlfac/ step,ider,nsx,swx,swy,npuls,pulse,   &
     &                             trp,ntypls,ntyrp,nchnl,topn,tcls,    &
     &                             ypls,gstrp,astrp,ycum
       
!  Accumulate the results of NPULS bipolar cycles by splining the instantaneous
!  response and folding the positive and negative parts of each cycle back
!  into a single pulse.
 
!     WRITE(*,nml=mynmlfac)
      
      
      CALL cubspl(trp,ypls,ntyrp,0,0)
      fold = 0.
      ipl = 1
      IF ( npuls==1 ) ipl = 0
 
      IF ( step==1 ) THEN
         DO jt = 1 , ntypls
            x = trp(jt)
            xp = x + pulse
            fold(jt) = cubval(trp,ypls,ntyrp,x)                         &
     &                 - ipl*cubval(trp,ypls,ntyrp,xp)
            DO jp = 2 , npuls
               x = xp + pulse
               xp = x + pulse
               fold(jt) = cubval(trp,ypls,ntyrp,x)                      &
     &                    - cubval(trp,ypls,ntyrp,xp) + fold(jt)
            ENDDO
         ENDDO
      ELSE
         DO jt = 1 , ntypls
            x = trp(jt)
            xp = x + pulse
            fold(jt) = ipl*cubder(trp,ypls,ntyrp,xp)                    &
     &                 - cubder(trp,ypls,ntyrp,x)
            DO jp = 2 , npuls
               x = xp + pulse
               xp = x + pulse
               fold(jt) = cubder(trp,ypls,ntyrp,xp)                     &
     &                    - cubder(trp,ypls,ntyrp,x) + fold(jt)
            ENDDO
         ENDDO
      ENDIF
 
      ypls = 0.
      ypls(1,1:ntypls) = fold(1:ntypls)
      CALL cubspl(trp,ypls,ntypls,0,0)
      ycum = 0.
 
!  Begin convolution.  If Geotem / Questem primary field stripping is required
!  the convolution must be done for all points in the waveform.
!  Otherwise, convolve only for those points needed in the windows.
 
!  The layered earth field is in IMPULSE form if dB/dt is desired
!  or in STEP form if B is to be computed.
 
      mxcnv = ntypls + nsx
      tf = swx(nsx)
      tfh = 0.5*tf
 
      IF ( gstrp==1 ) CALL tqstrip(ider,ntypls,trp,ypls,nsx,swx,swy,    &
     &                             ycnv)
      DO jt = 1 , nchnl
         t1 = topn(jt)
         t2 = tcls(jt)
         width = t2 - t1
         hwidth = width/2.
 
! Step response for step input or impulse response response for impulse input
! Average the response over receiver windows using 3 point Gaussian integration.
 
         tc(2) = (tcls(jt)+topn(jt))/2.
         tc(1) = tc(2) + hwidth*glx(1)
         tc(3) = tc(2) + hwidth*glx(3)
 
         DO jgl = 1 , 3
            t1 = tc(jgl)
            wt = glw(jgl)/2.
            yc1 = 0.
 
            IF ( gstrp==1 .AND. t1<tf ) THEN
               yc1 = cubval(swx,ycnv,nsx,t1)
            ELSEIF ( ider==0 ) THEN
                                ! Waveform input as I or B field (derived dI/dt)
               yc1 = txcnvl(t1,ntypls,trp,ypls,nsx,swx,swy)
               IF ( astrp==1 .AND. t1<tf ) THEN
                                            ! Asymmetric stripping
                  t2 = t1 - tfh
                  yc1 = yc1 + txcnvl(t2,ntypls,trp,ypls,nsx,swx,swy)
               ENDIF
 
            ELSEIF ( ider==1 ) THEN  ! Waveform input as voltage (known dI/dt)
               yc1 = txcnvd(mxcnv,t1,ntypls,trp,ypls,nsx,swx,swy)
            ELSEIF ( ider==4 ) THEN                ! pure on-time step
               yc1 = swy(1,1)*cubval(trp,ypls,ntypls,t1)
            ENDIF
            ycum(jt) = ycum(jt) + (wt*yc1)
         ENDDO
      ENDDO
 
      END SUBROUTINE fold_and_convolve
