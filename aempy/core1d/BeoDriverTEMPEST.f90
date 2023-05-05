      PROGRAM beodriver_tempest
!----------------------------------------------------------------------
!  Uses AEM1D_FD to compute  the frequency-domain layered earth
!  H field for a dipole of unit moment and current.
!.
!                             INPUT
!                             -----
!        JS = station reference
!      FREQ - array of NFRQ frequencies
!     TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!     TXA90 - true for vertical co-planar briadside array
!     NSTAT - number of stations in survey line.
!       ALT - array of transmitter altitudes
!       ZRX - vertical offset of each receiver from transmitter  (below = +)
!       XRX - in-line horizontal offset of RX J;  (behind = +)
!       YRX - transverse horizontal offset of RX J (port = +)
!      NLYR - number of layers
!       RES - layer resistivities
!      REPS - array of relative dislectric constants
!      RMUX - mu(i) / mu(0)
!       THK - array of layer thicknesses
!     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
!
!                             OUTPUT
!                             ------
!     BFD(JF,JS,1) - the in-line component of the layered earth respon6se at
!                    time JT, station JS. (nT)
!     BFD(JF,JS,2) - the transverse component
!     BFD(JF,JS,3) - the vertical component
!
!------------------------------------------------------------------------
!  SIGN CONVENTION:
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!------------------------------------------------------------------------
! LAST CHANGE:   26 Sep 2016   VR
!------------------------------------------------------------------------
!
      !USE core1d_module
      IMPLICIT NONE
      INTEGER , PARAMETER :: nmods=1, nlyr=3, ndata = 33
      INTEGER js, mode
      REAL*8 alt, &
     &    thk(nlyr),res(nlyr),reps(nlyr),rmu(nlyr),calf(nlyr),ctau(nlyr),cfreq(nlyr), &
     &    mvec(7*nlyr),cmp1,cmp2
      REAL*8, DIMENSION(ndata) :: calc_data


      alt            =  120.d0

      thk(1:nlyr-1)    = (/30.D0, 30.D0 /)
      res(1:nlyr)      = (/100.d0, 300.d0, 10.d0 /)
      reps(1:nlyr)     = 1.D0
      rmu(1:nlyr)      = 1.D0
      calf(1:nlyr)     = 0.D0
      ctau(1:nlyr)     = 0.D0
      cfreq(1:nlyr)    = 1.D0
       write(*,*)thk(1:nlyr-1)
       mvec(0*nlyr+1:1*nlyr)   = dlog10(res(1:nlyr))
       mvec(1*nlyr+1:2*nlyr)   = reps(1:nlyr)
       mvec(2*nlyr+1:3*nlyr)   = rmu(1:nlyr)
       mvec(3*nlyr+1:4*nlyr)   = calf(1:nlyr)
       mvec(4*nlyr+1:5*nlyr)   = ctau(1:nlyr)
       mvec(5*nlyr+1:6*nlyr)   = cfreq(1:nlyr)
       mvec(6*nlyr+1:7*nlyr-1)   = thk(1:nlyr-1)


      calc_data     = 0.D0


      CALL cpu_time(cmp1)
      mode =0
      DO js=0, nmods

             CALL aemfwd1d_tempest(mode,alt, nlyr,mvec,calc_data)

      ENDDO

      CALL cpu_time(cmp2)
      Write(*,*) 'cpu time for ',nmods+1,' 1-D models is ',cmp2-cmp1,' s'


      WRITE(*,'(/A/11G12.4/11G12.4/11G12.4)') ' result is: ',calc_data

      END
