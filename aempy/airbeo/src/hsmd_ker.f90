!*==hsmd_ker.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE hsmd_ker(iw,lmbda,nlyr,thkd,dispd,sig0,reps,rmux,icole,&
     &                    calf,ctau,cfreq,zrfd,qfd)
!-----------------------------------------------------------------------------
 
!***  Called by: HSMD_SNTR
 
!  Kernel for dipole transmitter and receiver above earth.
!
!          Input
!          -----
!      IW - iw = complex angular frequency
!   LMBDA = Hankel transform variable
!    NLYR - number of layers
!    SIG0 = real conductivities ( 1. / RES)
!    REPS - array of relative dielectric constants
!    RMUX - mu(i) / mu(0)
!    THKD - layer thicknesses
!   ICOLE - C-C layer indicator array. (0 for pure real, 1 for C-C layer)
!    CALF - complementary chargeability array; ie., CALF(I) = 1.0 - CHRG(I)
!    CTAU - array of layer relaxation times (sec).
!   CFREQ - array of layer frequency parameters.
!    ZRFD - reflected distance from transmitter to earth to receiver
!
!          Output
!          ------
!  QFD  forward model kernel
!
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80)
      REAL(KIND=Ql) , PARAMETER :: Eps0 = 8.854156D-12 ,                &
     &                             Mu0 = 12.56637D-7 , Exp_tol = 80.D0
      COMPLEX(KIND=Ql) , PARAMETER :: One = (1._QL,0._QL) ,             &
     &                                Zero = (0._QL,0._QL)
      INTEGER nlyr , icole(nlyr) , j
      REAL , DIMENSION(nlyr) :: calf , ctau , cfreq , reps
      REAL(KIND=Ql) xp0 , sig0(nlyr) , thkd(nlyr-1) , lmbda , zrfd ,    &
     &              rmux(nlyr) , rmusq(nlyr)
      COMPLEX(KIND=Ql) dispd , s0 , t0 , iw , p , lmbsq , ep ,          &
     &                 sigl(nlyr) , t(nlyr) , qfd
      COMPLEX(KIND=Ql) , DIMENSION(nlyr) :: f , xp , ksq , s , e
 
      ksq = Zero
      e = Zero
      xp = Zero
      t = Zero
 
      lmbsq = cmplx(lmbda*lmbda,0._QL,kind=Ql)
      s0 = sqrt(lmbsq-dispd)
      xp0 = -lmbda*exp(-lmbda*zrfd)
 
      sigl(1:nlyr) = cmplx(sig0(1:nlyr),0._QL,kind=Ql)
      DO j = nlyr , 1 , -1
         rmusq(j) = rmux(j)*rmux(j)
         p = (iw*ctau(j))**cfreq(j)
         p = icole(j)*p
         sigl(j) = sigl(j)*(One+p)/(One+calf(j)*p)
         sigl(j) = sigl(j) + iw*Eps0*reps(j)
                                            !  Add in displacement term
         ksq(j) = iw*Mu0*rmux(j)*sigl(j)
         s(j) = sqrt(ksq(j)+lmbsq)
 
         IF ( j==nlyr ) CYCLE
 
         ep = 2.D0*s(j)*thkd(j)
         IF ( real(ep)<Exp_tol ) e(j) = exp(-ep)
         t(j) = ((rmusq(j+1)-rmusq(j))*lmbsq+rmusq(j+1)*ksq(j)-rmusq(j) &
     &          *ksq(j+1))/(rmux(j+1)*s(j)+rmux(j)*s(j+1))**2
      ENDDO
 
      t0 = ((rmusq(1)-1._QL)*lmbsq-ksq(1))/(rmux(1)*s0+s(1))**2
      xp(1:nlyr-1) = e(1:nlyr-1)
      f = Zero
      DO j = nlyr - 1 , 1 , -1
         f(j) = xp(j)*(t(j)+f(j+1))/(One+t(j)*f(j+1))
      ENDDO
 
      qfd = xp0*(t0+f(1))/(One+t0*f(1))
 
      END SUBROUTINE hsmd_ker
