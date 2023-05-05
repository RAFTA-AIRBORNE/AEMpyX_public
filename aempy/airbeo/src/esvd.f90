!*==esvd.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE esvd(ajac,ndata,npar,kp,withu,withv,sv,umat,vmat,wsp,  &
     &                eta,tol,nw)
!----------------------------------------------------------------------------
 
!***  Called by: NLSQ2
!***      Calls: DPROD1, WRITE_MESSAGE
!
!==ESVD.spg  processed by SPAG 4.51A  at 17:56 on 17 Mar 1995
!
!  Singular value decomposition based on Golub's method with an option of
!  least squares reduction of RHS vectors for under or over determined
!  problems
!
!  REFERENCE:  G.H. Golub & C. Reinsch,
!             'Singular value decomposition and least squares solutions'
!              Num. Math.14,403-420(1970) ... (algol language) features both
!              SVD and minfit with accumulation in double precision
!
!     This routine is based on a program written by P. Businger
!     of Bell Telephone Lab., but incorporates a few modifications
!     due to R. Underwood of Stanford University, and D.L.B. Jupp
!==============================================================================
 
!          AJAC -  NDATA*(N+KP) array containing matrix to be decomposed
!         NDATA -  Number of data channels to be inverted
!          NPAR -  Number of parameters to be inverted.
!            KP -  If KP > 0 columns NPAR+1, ... ,NPAR+KP of AJAC contain KP
!                  'right hand sides'. these are multiplied by the transpose
!                  of UMAT for use in SOLVE2 (accompanying routine)
!   WITHU,WITHV -  logical control variables governing whether or not UMAT
!                  and VMAT respectively are constructed
!            SV -  on return SV(1) ... SV(N) contain the ordered singular
!                  values of AJAC. (SV(1) > SV(2) >  ...  > SV(NPAR))
!     UMAT,VMAT -  on return contain data and parameter space eigenvectors
!                  of AJAC. (AJAC =UMAT*SV*VMATtr)  depending upon the truth
!                  of WITHU and WITHV
!           WSP -  a workspace array of length 3*NPAR
!     ETA,TOL are machine dependent constants
!     ETA is the relative precision of real numbers
!     TOL is the smallest representable number divided by ETA
!
!     EXTERNAL REFERENCES ... DPROD1
!==================================================================
 
      IMPLICIT NONE
      INTEGER , PARAMETER :: Itmax = 30
      INTEGER i , im1 , iter , j , k , kk , kp , k1 , l , ll , l1 ,     &
     &        mkp , mp1 , ndata , nmk , npi , npk , npl , npar , np ,   &
     &        nw , n1 , n2
      REAL cs , eps , eta , f , fac , ftemp , g , h , q , r , sn , tol ,&
     &     w , x , y , z
      LOGICAL withu , withv
      REAL ajac(ndata,npar+1) , umat(ndata,npar) , vmat(npar,npar) ,    &
     &     sv(npar) , wsp(3*npar) , dprod1
      EXTERNAL dprod1
 
      np = npar + kp
      n1 = npar + 1
      n2 = npar + npar
      wsp(n1) = 0.e0
      k = 1
      DO
 
!  Householder reduction to upper bidiagonal form
 
         k1 = k + 1
         npk = npar + k1
 
!  Row transformation
 
         z = 0.0
         mkp = ndata - k + 1
         nmk = npar - k

         z = dprod1(mkp,1,1,z,ajac(k,k),ajac(k,k))
         wsp(k) = 0.0
         IF ( z>tol ) THEN
            z = sqrt(z)
            wsp(k) = z
            w = abs(ajac(k,k))
            fac = z + w
            q = 1.0
            IF ( ajac(k,k)<0.0 ) q = -1.0
            ajac(k,k) = q*fac
            IF ( k/=np ) THEN
               fac = z*fac
               DO j = k1 , np
                  q = 0.0
                  q = dprod1(mkp,1,1,q,ajac(k,k),ajac(k,j))
                  q = q/fac
                  ajac(k:ndata,j) = ajac(k:ndata,j) - q*ajac(k:ndata,k)
               ENDDO
 
!  Phase transformation
               IF ( ajac(k,k)>0.0 ) ajac(k,k1:np) = -ajac(k,k1:np)
            ENDIF
         ENDIF
 
         IF ( k==npar ) THEN
 
!  End of householder reduction.  Set tolerance for iteration.
 
            eps = 0.0
            DO k = 1 , npar
               npk = n2 + k
               npl = npar + k
               sv(k) = wsp(k)
               wsp(npk) = wsp(npl)
               eps = max(eps,sv(k)+wsp(npk))
            ENDDO
            eps = eps*eta
 
!  Set UMAT, and VMAT, to identity and preset pad of zero's
!  for case NDATA < NPAR.
 
            IF ( withu ) THEN
               umat = 0.
               DO j = 1 , npar
                  umat(j,j) = 1.0
               ENDDO
            ENDIF
 
            IF ( withv ) THEN
               vmat = 0.
               DO j = 1 , npar
                  vmat(j,j) = 1.0
               ENDDO
            ENDIF
 
            IF ( ndata<npar .AND. np>npar ) THEN
               mp1 = ndata + 1
               DO i = mp1 , npar
                  ajac(i,n1:np) = 0.
               ENDDO
            ENDIF
 
!  Main iteration loop on K ... Q-R algorithm due to Francis
 
            DO kk = 1 , npar
               k = n1 - kk
               npk = n2 + k
               iter = 0
 10            LOOP1:DO ll = 1 , k
                  l = k + 1 - ll
                  npl = n2 + l
                  IF ( abs(wsp(npl))<=eps ) GOTO 20
                  IF ( abs(sv(l-1))<=eps ) EXIT LOOP1
               ENDDO LOOP1
 
!  Cancellation
 
               cs = 0.0
               sn = 1.0
               l1 = l - 1
               LOOP2:DO i = l , k
                  npi = n2 + i
                  f = sn*wsp(npi)
                  wsp(npi) = cs*wsp(npi)
                  IF ( abs(f)<=eps ) EXIT LOOP2
                  h = sv(i)
                  IF ( abs(f)<=abs(h) ) THEN
                     sn = f/h
                     w = sqrt(sn*sn+1.0)
                     sv(i) = h*w
                     cs = 1.0/w
                     sn = -sn*cs
                  ELSE
                     cs = h/f
                     w = sqrt(cs*cs+1.0)
                     sv(i) = f*w
                     sn = -1.0/w
                     cs = -cs*sn
                  ENDIF
                  IF ( withu ) THEN
                     DO j = 1 , npar
                        x = umat(j,l1)
                        y = umat(j,i)
                        umat(j,l1) = x*cs + y*sn
                        umat(j,i) = y*cs - x*sn
                     ENDDO
                  ENDIF
                  IF ( np/=npar ) THEN
                     DO j = n1 , np
                        q = ajac(l1,j)
                        r = ajac(i,j)
                        ajac(l1,j) = q*cs + r*sn
                        ajac(i,j) = r*cs - q*sn
                     ENDDO
                  ENDIF
               ENDDO LOOP2
 
!  Test for convergence
 
 20            w = sv(k)
               IF ( l/=k ) THEN
 
!     TEST FOR MAX ITERATIONS
 
                  iter = iter + 1
                  IF ( iter<=Itmax ) THEN
 
!  Compute the implicit shift of origin from bottom 2x2 minor
 
                     x = sv(l)
                     y = sv(k-1)
                     g = wsp(npk-1)
                     h = wsp(npk)
                     f = (y-w)*(y+w) + (g-h)*(g+h)
                     ftemp = 2.0*h*y
                     IF ( abs(f)>abs(ftemp) ) THEN
                        f = ftemp/f
                        g = sqrt(f*f+1.0)
                        f = 1.0/f
                        g = g*f
                     ELSE
                        f = f/ftemp
                        g = sqrt(f*f+1.0)
                        IF ( f<0.0 ) g = -g
                     ENDIF
                     f = ((x-w)*(x+w)+(y/(f+g)-h)*h)/x
                     cs = 1.0
                     sn = 1.0
                     l1 = l + 1
 
!  Main loop Q-R transformation for SV(K)
 
                     DO i = l1 , k
                        im1 = i - 1
                        npi = n2 + i
                        g = wsp(npi)
                        y = sv(i)
                        h = sn*g
                        g = cs*g
 
!  Givens rotation from the right
 
                        IF ( abs(f)<=abs(h) ) THEN
                           cs = f/h
                           w = sqrt(cs*cs+1.0)
                           wsp(npi-1) = h*w
                           sn = 1.0/w
                           cs = cs*sn
                        ELSE
                           sn = h/f
                           w = sqrt(sn*sn+1.0)
                           wsp(npi-1) = f*w
                           cs = 1.0/w
                           sn = sn*cs
                        ENDIF
                        f = x*cs + g*sn
                        g = g*cs - x*sn
                        h = y*sn
                        y = y*cs
                        IF ( withv ) THEN
                           DO j = 1 , npar
                              x = vmat(j,im1)
                              w = vmat(j,i)
                              vmat(j,im1) = x*cs + w*sn
                              vmat(j,i) = w*cs - x*sn
                           ENDDO
                        ENDIF
 
!  Givens rotation from the left
 
                        IF ( abs(f)<=abs(h) ) THEN
                           cs = f/h
                           w = sqrt(cs*cs+1.0)
                           sv(im1) = h*w
                           sn = 1.0/w
                           cs = cs*sn
                        ELSE
                           sn = h/f
                           w = sqrt(sn*sn+1.0)
                           sv(im1) = f*w
                           cs = 1.0/w
                           sn = sn*cs
                        ENDIF
                        f = cs*g + sn*y
                        x = cs*y - sn*g
                        IF ( withu ) THEN
                           DO j = 1 , npar
                              y = umat(j,im1)
                              w = umat(j,i)
                              umat(j,im1) = y*cs + w*sn
                              umat(j,i) = w*cs - y*sn
                           ENDDO
                        ENDIF
                        IF ( npar/=np ) THEN
                           DO j = n1 , np
                              q = ajac(im1,j)
                              r = ajac(i,j)
                              ajac(im1,j) = q*cs + r*sn
                              ajac(i,j) = r*cs - q*sn
                           ENDDO
                        ENDIF
                     ENDDO
                     wsp(npl) = 0.0
                     wsp(npk) = f
                     sv(k) = x
                     GOTO 10
                  ELSE
                     WRITE (nw,'(//T3,A,I3)')                           &
     &' MESSAGE FROM ESVD: Maximum iterations exceeded for singular valu&
     &e.' , k
                  ENDIF
               ENDIF
 
!  Convergence ... if singular value negative, make it positive
 
               IF ( w<0.0e0 ) THEN
                  sv(k) = -w
                  IF ( withv ) vmat(1:npar,k) = -vmat(1:npar,k)
               ENDIF
            ENDDO
 
!  End of main loop.  Order singular values, UMAT, VMAT, and RHS'S.
 
            DO k = 1 , npar
               g = -1.0
               j = k
               DO i = k , npar
                  IF ( sv(i)>g ) THEN
                     g = sv(i)
                     j = i
                  ENDIF
               ENDDO
               IF ( j/=k ) THEN
                  sv(j) = sv(k)
                  sv(k) = g
                  IF ( withv ) THEN
                     DO i = 1 , npar
                        q = vmat(i,j)
                        vmat(i,j) = vmat(i,k)
                        vmat(i,k) = q
                     ENDDO
                  ENDIF
                  IF ( withu ) THEN
                     DO i = 1 , npar
                        q = umat(i,j)
                        umat(i,j) = umat(i,k)
                        umat(i,k) = q
                     ENDDO
                  ENDIF
                  IF ( npar/=np ) THEN
                     DO i = n1 , np
                        q = ajac(j,i)
                        ajac(j,i) = ajac(k,i)
                        ajac(k,i) = q
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
 
!  Update umat with stored Householder transformations
 
            IF ( withu ) THEN
               DO kk = 1 , npar
                  k = n1 - kk
                  IF ( abs(wsp(k))>tol ) THEN
                     IF ( k<=ndata ) THEN
                        mkp = ndata - k + 1
                        fac = abs(ajac(k,k))*wsp(k)
 
!  Undo the phase
                        IF ( ajac(k,k)>0.0 ) umat(k,1:npar)             &
     &                       = -umat(k,1:npar)
                        DO j = 1 , npar
                           q = 0.0
                           q = dprod1(mkp,1,1,q,ajac(k,k),umat(k,j))
                           q = q/fac
                           umat(k:ndata,j) = umat(k:ndata,j)            &
     &                        - q*ajac(k:ndata,k)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
 
!  Update VMAT with stored householder transformations
 
            IF ( withv ) THEN
               IF ( npar>=2 ) THEN
                  DO kk = 2 , npar
                     k = n1 - kk
                     k1 = k + 1
                     npk = npar + k1
                     IF ( abs(wsp(npk))>tol ) THEN
                        IF ( k<=ndata ) THEN
                           nmk = npar - k
                           fac = abs(ajac(k,k1))*wsp(npk)
 
!  Undo the phase
 
                           IF ( ajac(k,k1)>0.0 ) vmat(k1,1:npar)        &
     &                          = -vmat(k1,1:npar)
                           DO j = 1 , npar
                              q = 0.0
                              q = dprod1(nmk,ndata,1,q,ajac(k,k1),      &
     &                            vmat(k1,j))
                              q = q/fac
                              vmat(k1:npar,j) = vmat(k1:npar,j)         &
     &                           - q*ajac(k,k1:npar)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
            RETURN
         ELSE
 
!  Column transformation
 
            z = 0.0
            IF ( k<=ndata ) z = dprod1(nmk,ndata,ndata,z,ajac(k,k1),    &
     &                          ajac(k,k1))
            wsp(npk) = 0.0
            IF ( z>tol ) THEN
               z = sqrt(z)
               wsp(npk) = z
               w = abs(ajac(k,k1))
               fac = z + w
               q = 1.0
               IF ( ajac(k,k1)<0.0 ) q = -1.0
               ajac(k,k1) = q*fac
               IF ( ndata>k ) THEN
                  fac = z*fac
                  DO i = k1 , ndata
                     q = 0.0
                     q = dprod1(nmk,ndata,ndata,q,ajac(k,k1),ajac(i,k1))
                     q = q/fac
                     ajac(i,k1:npar) = ajac(i,k1:npar)                  &
     &                                 - q*ajac(k,k1:npar)
                  ENDDO
 
!  Phase transformation
                  IF ( ajac(k,k1)>0.0 ) ajac(k1:ndata,k1)               &
     &                 = -ajac(k1:ndata,k1)
               ENDIF
            ENDIF
            k = k1
         ENDIF
      ENDDO
 
      END SUBROUTINE esvd
