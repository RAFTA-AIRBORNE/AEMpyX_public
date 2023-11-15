SUBROUTINE unpack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mvec)
!f2py depend(nlyr)  mvec,res,reps,rmu,calf,ctau,cfreq,thk
!f2py intent(in) nlyr, mvec
!f2py intent (out) res,reps,rmu,calf,ctau,cfreq,thk
!f2py threadsafe
!------------------------------------------------------------------------
! Unpacks the parameter vector to physical properties
!
!      NLYR - number of layers
!      RES - layer resistivities (nlyr)
!      REPS - array of relative dielectric constants (nlyr)
!      RMUX - mu(i) / mu(0)   (nlyr)
!      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters. (nlyr)
!      THK - array of layer thicknesses (nlyr-1)
!
!------------------------------------------------------------------------
! LAST CHANGE:   1 Sep 2016   VR
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in ) :: nlyr
      REAL(KIND=8), INTENT(in) ::  mvec(7*nlyr)
      REAL(KIND=8), INTENT(out) :: thk(nlyr),res(nlyr),reps(nlyr),rmu(nlyr),calf(nlyr),ctau(nlyr),&
     & cfreq(nlyr)
     
       res(1:nlyr)      = mvec(0*nlyr+1:1*nlyr)
       rmu(1:nlyr)      = mvec(1*nlyr+1:2*nlyr)
       reps(1:nlyr)     = mvec(2*nlyr+1:3*nlyr)
       calf(1:nlyr)     = mvec(3*nlyr+1:4*nlyr)
       ctau(1:nlyr)     = mvec(4*nlyr+1:5*nlyr)
       cfreq(1:nlyr)    = mvec(5*nlyr+1:6*nlyr)
       thk(1:nlyr-1)    = mvec(6*nlyr+1:7*nlyr-1)

END SUBROUTINE unpack_mvec

SUBROUTINE pack_mvec(nlyr,res,reps,rmu,calf,ctau,cfreq,thk,mvec)
!f2py intent(out) mvec
!f2py depend(nlyr)  mvec,res,reps,rmu,calf,ctau,cfreq,thk
!f2py intent (in) nlyr,res,reps,rmu,calf,ctau,cfreq,thk
!f2py threadsafe

!------------------------------------------------------------------------
! packs the parameter vector to physical properties
!
!      NLYR - number of layers
!      RES - layer resistivities (nlyr)
!      REPS - array of relative dislectric constants (nlyr)
!      RMUX - mu(i) / mu(0)   (nlyr)
!      CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters. (nlyr)
!      THK - array of layer thicknesses (nlyr-1)
!
!------------------------------------------------------------------------
! LAST CHANGE:   1 Sep 2016   VR
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in ) :: nlyr
      REAL(KIND=8), INTENT(out) :: mvec(7*nlyr)
      REAL(KIND=8), INTENT(in ) :: thk(nlyr),res(nlyr),reps(nlyr),rmu(nlyr),calf(nlyr),ctau(nlyr),&
     & cfreq(nlyr)

       mvec(0*nlyr+1:1*nlyr) = res(1:nlyr)
       mvec(1*nlyr+1:2*nlyr) = rmu(1:nlyr)
       mvec(2*nlyr+1:3*nlyr) = reps(1:nlyr)
       mvec(3*nlyr+1:4*nlyr) = calf(1:nlyr)
       mvec(4*nlyr+1:5*nlyr) = ctau(1:nlyr)
       mvec(5*nlyr+1:6*nlyr) = cfreq(1:nlyr)
       mvec(6*nlyr+1:7*nlyr-1) = thk(1:nlyr-1)

END SUBROUTINE pack_mvec

SUBROUTINE trans_mvec(transform,nlyr,mvec,nvec)
!f2py intent(in)    mvec
!f2py intent(out)   nvec, delt
!f2py intent (in)   transform, nlyr
!f2py depend(nlyr)  mvec, nvec, delt
!f2py threadsafe
!------------------------------------------------------------------------,
! transforms the parameter vector to pre-defined functions for inversion
! generates appropriate DD pertutbations
!
!      NLYR - number of layers
!
!------------------------------------------------------------------------
! LAST CHANGE:   5 DEC 2020   VR
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)           :: nlyr, transform
      REAL(KIND=8), INTENT(in)      :: mvec(7*nlyr)      
      REAL(KIND=8), INTENT(out)     :: nvec(7*nlyr)


      nvec = mvec

!       
!       write(*,*) 'fbeg nlyr :',nlyr
!       write(*,*) 'fbeg mvec :', mvec
!       write(*,*) 'fbeg delt:', delt
      
      SELECT CASE (abs(transform))
          CASE (0)
            nvec(0*nlyr+1:1*nlyr)   = dlog(mvec(0*nlyr+1:1*nlyr))
            nvec(6*nlyr+1:7*nlyr-1) = dlog(mvec(6*nlyr+1:7*nlyr-1))
          CASE (1)
! all log Hoenig 2002
            nvec(0*nlyr+1:1*nlyr) = dlog(mvec(0*nlyr+1:1*nlyr))
            nvec(3*nlyr+1:4*nlyr) = dlog(mvec(3*nlyr+1:4*nlyr))
            nvec(4*nlyr+1:5*nlyr) = dlog(mvec(4*nlyr+1:5*nlyr))
            nvec(5*nlyr+1:6*nlyr) = dlog(mvec(5*nlyr+1:6*nlyr))
            nvec(6*nlyr+1:7*nlyr-1) = dlog(mvec(6*nlyr+1:7*nlyr-1))
        CASE (2)
! 10^mu/(10^mu + 1), INVERSE normalized chargeability Ghorbani 2007
            nvec(0*nlyr+1:1*nlyr) = dlog(mvec(0*nlyr+1:1*nlyr))
            nvec(3*nlyr+1:4*nlyr) = dlog(mvec(3*nlyr+1:4*nlyr)/(1D0-mvec(3*nlyr+1:4*nlyr)))
            nvec(4*nlyr+1:5*nlyr) = dlog(mvec(4*nlyr+1:5*nlyr))
            nvec(5*nlyr+1:6*nlyr) = dlog(mvec(5*nlyr+1:6*nlyr))
            nvec(6*nlyr+1:7*nlyr-1) = dlog(mvec(6*nlyr+1:7*nlyr-1))
    END SELECT
!     write(*,*) 'fend nlyr :',nlyr
!     write(*,*) 'fend mvec :',nvec
!     write(*,*) 'fend delt:', delt
            
END SUBROUTINE trans_mvec
        
SUBROUTINE untrans_mvec(transform,nlyr,mvec,nvec)
!f2py intent(in)    mvec,transform, nlyr
!f2py intent(out)   nvec
!f2py depend(nlyr)  mvec, nvec
!f2py threadsafe
!------------------------------------------------------------------------,
! transforms the parameter vector fropm pre-defined functions for inversion
! generates appropriate DD pertutbations
!
!      NLYR - number of layers
!
!------------------------------------------------------------------------
! LAST CHANGE:   5 DEC 2020   VR
!------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)              :: nlyr, transform
      REAL(KIND=8), INTENT(in)         :: mvec(7*nlyr)      
      REAL(KIND=8), INTENT(out)        :: nvec(7*nlyr)
      

      nvec = mvec  
   
      
      SELECT CASE (abs(transform))
          CASE (0)
      
            write(*,*) 'before:'
            write(*,'(10g12.4)') mvec(0*nlyr+1:1*nlyr)            
            write(*,'(10g12.4)') mvec(6*nlyr+1:7*nlyr)
            nvec(0*nlyr+1:1*nlyr)   = dexp(mvec(0*nlyr+1:1*nlyr))
            write(*,*) ' after:'
            write(*,'(10g12.4)') nvec(0*nlyr+1:1*nlyr)   
            write(*,'(10g12.4)') nvec(6*nlyr+1:7*nlyr)

            nvec(6*nlyr+1:7*nlyr-1) = dexp(mvec(6*nlyr+1:7*nlyr-1))
            write(*,'(10g12.4)') nvec(6*nlyr+1:7*nlyr)
          CASE (1)
! all log Hoenig 2002
            nvec(0*nlyr+1:1*nlyr) = dexp(mvec(0*nlyr+1:1*nlyr))
            nvec(3*nlyr+1:4*nlyr) = dexp(mvec(3*nlyr+1:4*nlyr))
            nvec(4*nlyr+1:5*nlyr) = dexp(mvec(4*nlyr+1:5*nlyr))
            nvec(5*nlyr+1:6*nlyr) = dexp(mvec(5*nlyr+1:6*nlyr))
            nvec(6*nlyr+1:7*nlyr-1) = dexp(mvec(6*nlyr+1:7*nlyr-1))

        CASE (2)
! 10^mu/(10^mu + 1), INVERSE normalized chargeability Ghorbani 2007
            nvec(0*nlyr+1:1*nlyr) = dexp(mvec(0*nlyr+1:1*nlyr))
            nvec(3*nlyr+1:4*nlyr) = dexp(mvec(3*nlyr+1:4*nlyr))/ (dexp(mvec(3*nlyr+1:4*nlyr))+1.d0)
            nvec(4*nlyr+1:4*nlyr) = dexp(mvec(4*nlyr+1:5*nlyr))
            nvec(5*nlyr+1:6*nlyr) = dexp(mvec(5*nlyr+1:6*nlyr))
            nvec(6*nlyr+1:7*nlyr-1) = dexp(mvec(6*nlyr+1:7*nlyr-1))
      END SELECT

      write(*,*) 'nvec :', nvec

END SUBROUTINE untrans_mvec
