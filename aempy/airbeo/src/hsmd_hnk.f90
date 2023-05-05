!*==hsmd_hnk.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE hsmd_hnk(iw,nlyr,thkd,dispd,sig0,reps,rmux,icole,calf, &
     &                    ctau,cfreq,zrfd,rhod,hlyr)
!-------------------------------------------------------------------------------------------
 
!***  Called by: HSMD_FD
!***      Calls: HS_JMP, HSMD_KER
 
!  VALID FOR ANY NUMBER OF LAYERS
!  Magnetic dipole transmitter & receiver above or on earth surface
 
!  Computes transform integrals HLYR(3) which are used to compute vertical
!  and horizontal frequency-domain magnetic field components at the RX from
!  VMD and HMD sources.  It evaluates the Hankel transform integral using a
!  15 points per decade filter coefficient set derived from Christensen's
!  FLTGEN program.
 
!      IW - iw  angular frequency *(0.,1.)
!    RHOD - horizontal TX -> RX distance.
!     KER - stores kernel values from HSMD_KER
 
!  NLYR,SIG0,REPS,RMUX,THK,ICOLE,CALF,CTAU,CFREQ,ZRFD
!  are described in HSMD_FD
!
!    OUTPUT is HLYR(1:3)  forward model components
 
      USE filt_coef_q
 
      IMPLICIT NONE
      REAL(KIND=Ql) , PARAMETER :: Vfac0 = 100._QL
      INTEGER nlyr , i , icole(nlyr)
      REAL , DIMENSION(nlyr) :: reps , ctau , cfreq , calf
      REAL(KIND=Ql) del_jn , rho_jn , y , lmbda , rhod , zrfd ,         &
     &              sig0(nlyr) , thkd(nlyr-1) , rmux(nlyr)
      COMPLEX(KIND=Ql) iw , hlyr(3) , dispd , qfd
      LOGICAL jump
 
      del_jn = log(10.D0)/dble(Ndec_jn)
      rho_jn = -log(rhod) - shftjn
 
      hlyr = (0._QL,0._QL)
 
      DO i = -50 , Jnhi       ! Start at I = -50 to pick up low values.
         y = rho_jn + dble(i)*del_jn
         lmbda = exp(y)
         CALL hsmd_ker(iw,lmbda,nlyr,thkd,dispd,sig0,reps,rmux,icole,   &
     &                 calf,ctau,cfreq,zrfd,qfd)
         CALL hs_frq_jmp
         IF ( jump .AND. i>-40 ) EXIT
      ENDDO
 
      jump = .FALSE.      ! Finish off the low end for RHOTRP(1)
      DO i = -51 , Jnlo , -1
         y = rho_jn + dble(i)*del_jn
         lmbda = exp(y)
         CALL hsmd_ker(iw,lmbda,nlyr,thkd,dispd,sig0,reps,rmux,icole,   &
     &                 calf,ctau,cfreq,zrfd,qfd)
 
         CALL hs_frq_jmp
         IF ( jump .AND. i<-60 ) EXIT
      ENDDO
 
      hlyr = Vfac0*hlyr/rhod
 
      CONTAINS
 
      SUBROUTINE hs_frq_jmp
!  ---------------------
 
!***  Called by: HSMD_HNK
 
!  Accumulates function calls for the Hankel transformation &
!  checks convergence.
 
      REAL(KIND=Ql) , PARAMETER :: Tol = 1.D-6 , Tol2 = 1.D-35
      INTEGER jint
      REAL(KIND=Ql) qr , qi
      COMPLEX(KIND=Ql) fw(3)
 
      fw(1) = wj0(i)*qfd*lmbda
      fw(2) = wj1(i)*qfd*lmbda
      fw(3) = wj1(i)*qfd/rhod
 
      hlyr = hlyr + fw
 
      jump = .TRUE.
      DO jint = 1 , 3
         qr = abs(real(hlyr(jint),kind=Ql))
         qi = abs(aimag(hlyr(jint)))
         IF ( qr>Tol2 .AND. abs(real(fw(jint)))>Tol*qr ) jump = .FALSE.
         IF ( qi>Tol2 .AND. abs(aimag(fw(jint)))>Tol*qi ) jump = .FALSE.
      ENDDO
 
      END SUBROUTINE hs_frq_jmp
 
      END SUBROUTINE hsmd_hnk
