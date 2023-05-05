!*==hsboss_fd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE hsboss_fd(js,nfrq,freq,txcln,txa90,nstat,sz,zrx,xrx,   &
     &                     yrx,nlyr,res,reps,rmu,thk,calf,ctau,cfreq,   &
     &                     bfd)
!----------------------------------------------------------------------
 
!  Computes the frequency-domain layered earth H field for a dipole of
!  unit moment and current.
 
!***  Called by: MAIN, GET_FWD_MODL
!***      Calls: HSMD_FD
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  SIGN CONVENTION:
!  ----------------
!  The normal layered earth field coordinate system used in this
!  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
!  positive to starboard, and Z (JC=3) positive down.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 
!                             INPUT
!                             -----
!        JS = station reference
!      FREQ - array of NFRQ frequencies
!     TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
!     TXA90 - true for vertical co-planar briadside array
!     NSTAT - number of stations in survey line.
!        SZ - array of transmitter altitudes
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
!     BFD(JF,JS,1) - the in-line component of the layered earth response at
!                    time JT, station JS. (nT)
!     BFD(JF,JS,2) - the transverse component
!     BFD(JF,JS,3) - the vertical component
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80)
      INTEGER nfrq , nstat , nlyr , js , tdfd , nrxf
      REAL sz(nstat) , alt
      REAL , DIMENSION(nfrq) :: freq , txcln , zrx , xrx , yrx
      REAL , DIMENSION(nlyr) :: res , rmu , reps , ctau , cfreq , calf
      REAL , DIMENSION(nlyr-1) :: thk
      REAL(KIND=Ql) , DIMENSION(nfrq) :: xrxd , yrxd , zrxd , txclnd
      COMPLEX(KIND=Ql) bfdd(nfrq,3)
      COMPLEX bfd(nfrq,nstat,3)
      LOGICAL txa90
 
!  Compute layered earth fields BLE_LYR at first station for each different altitude.
 
      tdfd = 2
      nrxf = nfrq
      zrxd(1:nfrq) = real(zrx(1:nfrq),kind=Ql)
      xrxd(1:nfrq) = real(xrx(1:nfrq),kind=Ql)
      yrxd(1:nfrq) = real(yrx(1:nfrq),kind=Ql)
      txclnd(1:nfrq) = real(txcln(1:nfrq),kind=Ql)
 
      alt = sz(js)
      CALL hsmd_fd(nfrq,freq,alt,nrxf,txclnd,txa90,zrxd,xrxd,yrxd,nlyr, &
     &             thk,res,reps,rmu,calf,ctau,cfreq,tdfd,bfdd)
 
      bfd(1:nfrq,js,1:3) = cmplx(bfdd(1:nfrq,1:3))
 
      END SUBROUTINE hsboss_fd
