!*==interv.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE interv(xt,lxt,x,left,mflag)
!-------------------------------------------
 
!***   Called by: CUBVAL, CUBINT, CUBDER, LINVAL
 
!---  Restructured April, 1997
 
!  from  * A PRACTICAL GUIDE TO SPLINES *  by C. DE BOOR
!  computes  LEFT = MAX( I , 1 <= I <= LXT  .AND.  XT(I) <= X )  .
!
!             INPUT
!             -----
!       XT - a real sequence, of length  LXT, assumed to be non-decreasing.
!      LXT - number of terms in the sequence  XT .
!        X - the point whose location with respect to the sequence XT is
!            to be determined.
!
!             OUTPUT
!             ------
!      LEFT, MFLAG.....are both integers, whose value is:
!
!        1     -1      IF               X <  XT(1)
!        I      0      IF   XT(I)  <= X < XT(I+1)
!       LXT     1      IF  XT(LXT) <= X
!
!        In particular, MFLAG = 0 is the 'usual' case.  MFLAG /= 0
!        indicates that X  lies outside the halfopen interval
!        XT(1) <= Y < XT(LXT) . The asymmetric treatment of the
!        interval is due to the decision to make all pp functions
!        continuous from the right.
!
!             METHOD
!             ------
!
!  The program is designed to be efficient in the common situation that
!  it is called repeatedly, with  X  taken from an increasing or decreasing
!  sequence. This will happen, e.g., when a pp function is to be grapged.
!  The first guess for  LEFT  is therefore taken to be the value returned at
!  the previous call and stored in the  L O C A L  variable ILO. A first
!  check ascertains that  ILO < LXT (This is necessary since the present
!  call may have nothing to do with the previous call).
!  Then, if XT(ILO) <= XT(ILO+1),
!  we set  LEFT = ILO  and are done after just three comparisons.
!  Otherwise, we repeatedly double the difference  ISTEP = IHI - ILO
!  while also moving  ILO  AND  IHI  in the direction of  X , until
!                      XT(ILO) <= X < XT(IHI) ,
!  after which we use bisection to get, in addition, ILO+1 = IHI .
!  LEFT = ILO  is then returned.
!******************************************************************************
!******************************************************************************
 
      IMPLICIT NONE
      INTEGER left , lxt , mflag , ihi , ilo , istep , middle , j1
      REAL x , xt(lxt)
      SAVE ilo
 
      DATA ilo/1/
 
!***********************************************************
!  Trivial returns when X is not in the range.
 
      IF ( (x<=xt(1)) .OR. (lxt<=1) ) THEN
         mflag = -1
         left = 1
         RETURN
      ENDIF
 
      IF ( x>=xt(lxt) ) THEN
         mflag = 1
         left = lxt
         RETURN
      ENDIF
 
      mflag = 0
      IF ( ilo>=lxt ) ilo = lxt - 1
      ihi = ilo + 1
 
!  Trivial return when X is already in the interval.
 
      IF ( (x<=xt(ihi)) .AND. (x>=xt(ilo)) ) THEN
         left = ilo
         RETURN
      ENDIF
!***********************************************************
 
      IF ( x<=xt(ilo) ) THEN
                          ! decrease ILO  to capture X.
         istep = 1
         DO j1 = 1 , lxt
            ihi = ilo
            ilo = ihi - istep
            ilo = max(1,ilo)
            IF ( (x>=xt(ilo)) .OR. (ilo==1) ) EXIT
            istep = istep*2
         ENDDO
 
      ELSEIF ( x>=xt(ihi) ) THEN
                                ! increase IHI to capture X
 
         istep = 1
         DO j1 = 1 , lxt
            ilo = ihi
            ihi = ilo + istep
            ihi = min(ihi,lxt)
            IF ( (x<=xt(ihi)) .OR. (ihi==lxt) ) EXIT
            istep = istep*2
         ENDDO
 
      ENDIF
 
!  Now XT(ILO) <= X < XT(IHI) . Narrow the interval.
 
      DO j1 = 1 , lxt
         middle = (ilo+ihi)/2
         IF ( middle==ilo ) EXIT
         IF ( x<xt(middle) ) THEN
            ihi = middle
         ELSE
            ilo = middle
         ENDIF
      ENDDO
 
! Task complete
 
      left = ilo
 
      END SUBROUTINE interv
