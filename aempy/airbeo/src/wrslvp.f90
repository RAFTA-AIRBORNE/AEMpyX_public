!*==wrslvp.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE wrslvp(nw,sxd,syd,sz,nstat,nchnl,tms,ytr)
!-----------------------------------------------------
 
!***  Called by: WRSLV
 
!  Writes time-domain output in profile form for receiver
      IMPLICIT NONE
      INTEGER , PARAMETER :: Ql = selected_real_kind(12,80) , Ncol = 30
      INTEGER nw , nstat , nchnl , nblks , j1 , jb , jf , krx , kry ,   &
     &        krz , js
      REAL sz(nstat,1) , tms(nchnl) , ytr(nchnl,nstat)
      REAL(KIND=Ql) sxd(nstat,1) , syd(nstat,1)
      CHARACTER(LEN=8) chn(150)
      DATA chn(1:150)/' CHNL 1 ' , ' CHNL 2 ' , ' CHNL 3 ' ,            &
     &     ' CHNL 4 ' , ' CHNL 5 ' , ' CHNL 6 ' , ' CHNL 7 ' ,          &
     &     ' CHNL 8 ' , ' CHNL 9 ' , 'CHNL 10 ' , 'CHNL 11 ' ,          &
     &     'CHNL 12 ' , 'CHNL 13 ' , 'CHNL 14 ' , 'CHNL 15 ' ,          &
     &     'CHNL 16 ' , 'CHNL 17 ' , 'CHNL 18 ' , 'CHNL 19 ' ,          &
     &     'CHNL 20 ' , 'CHNL 21 ' , 'CHNL 22 ' , 'CHNL 23 ' ,          &
     &     'CHNL 24 ' , 'CHNL 25 ' , 'CHNL 26 ' , 'CHNL 27 ' ,          &
     &     'CHNL 28 ' , 'CHNL 29 ' , 'CHNL 30 ' , 'CHNL 31 ' ,          &
     &     'CHNL 32 ' , 'CHNL 33 ' , 'CHNL 34 ' , 'CHNL 35 ' ,          &
     &     'CHNL 36 ' , 'CHNL 37 ' , 'CHNL 38 ' , 'CHNL 39 ' ,          &
     &     'CHNL 40 ' , 'CHNL 41 ' , 'CHNL 42 ' , 'CHNL 43 ' ,          &
     &     'CHNL 44 ' , 'CHNL 45 ' , 'CHNL 46 ' , 'CHNL 47 ' ,          &
     &     'CHNL 48 ' , 'CHNL 49 ' , 'CHNL 50 ' , 'CHNL 51 ' ,          &
     &     'CHNL 52 ' , 'CHNL 53 ' , 'CHNL 54 ' , 'CHNL 55 ' ,          &
     &     'CHNL 56 ' , 'CHNL 57 ' , 'CHNL 58 ' , 'CHNL 59 ' ,          &
     &     'CHNL 60 ' , 'CHNL 61 ' , 'CHNL 62 ' , 'CHNL 63 ' ,          &
     &     'CHNL 64 ' , 'CHNL 65 ' , 'CHNL 66 ' , 'CHNL 67 ' ,          &
     &     'CHNL 68 ' , 'CHNL 69 ' , 'CHNL 70 ' , 'CHNL 71 ' ,          &
     &     'CHNL 72 ' , 'CHNL 73 ' , 'CHNL 74 ' , 'CHNL 75 ' ,          &
     &     'CHNL 76 ' , 'CHNL 77 ' , 'CHNL 78 ' , 'CHNL 79 ' ,          &
     &     'CHNL 80 ' , 'CHNL 81 ' , 'CHNL 82 ' , 'CHNL 83 ' ,          &
     &     'CHNL 84 ' , 'CHNL 85 ' , 'CHNL 86 ' , 'CHNL 87 ' ,          &
     &     'CHNL 88 ' , 'CHNL 89 ' , 'CHNL 90 ' , 'CHNL 91 ' ,          &
     &     'CHNL 92 ' , 'CHNL 93 ' , 'CHNL 94 ' , 'CHNL 95 ' ,          &
     &     'CHNL 96 ' , 'CHNL 97 ' , 'CHNL 98 ' , 'CHNL 99 ' ,          &
     &     'CHNL 100' , 'CHNL 101' , 'CHNL 102' , 'CHNL 103' ,          &
     &     'CHNL 104' , 'CHNL 105' , 'CHNL 106' , 'CHNL 107' ,          &
     &     'CHNL 108' , 'CHNL 109' , 'CHNL 110' , 'CHNL 111' ,          &
     &     'CHNL 112' , 'CHNL 113' , 'CHNL 114' , 'CHNL 115' ,          &
     &     'CHNL 116' , 'CHNL 117' , 'CHNL 118' , 'CHNL 119' ,          &
     &     'CHNL 120' , 'CHNL 121' , 'CHNL 122' , 'CHNL 123' ,          &
     &     'CHNL 124' , 'CHNL 125' , 'CHNL 126' , 'CHNL 127' ,          &
     &     'CHNL 128' , 'CHNL 129' , 'CHNL 130' , 'CHNL 131' ,          &
     &     'CHNL 132' , 'CHNL 133' , 'CHNL 134' , 'CHNL 135' ,          &
     &     'CHNL 136' , 'CHNL 137' , 'CHNL 138' , 'CHNL 139' ,          &
     &     'CHNL 140' , 'CHNL 141' , 'CHNL 142' , 'CHNL 143' ,          &
     &     'CHNL 144' , 'CHNL 145' , 'CHNL 146' , 'CHNL 147' ,          &
     &     'CHNL 148' , 'CHNL 149' , 'CHNL 150'/
 
      nblks = nchnl/Ncol
      IF ( mod(nchnl,Ncol)>0 ) nblks = nblks + 1
 
      DO j1 = 1 , nblks
         jf = j1*Ncol
         jb = jf - Ncol + 1
         jf = min(jf,nchnl)
 
         WRITE (nw,99001) chn(jb:jf)
 
99001    FORMAT (T10,'TRANSMITTER POSITION',T35,30(A:,5X))
         WRITE (nw,99002) tms(jb:jf)
99002    FORMAT (T10,'EAST     NORTH    ALT',F10.3,29F13.3)
         WRITE (nw,'(3X)')
         DO js = 1 , nstat
            krx = 0
            kry = 0
            krz = nint(sz(js,1))
            krx = nint(sxd(js,1))
            kry = nint(syd(js,1))
            WRITE (nw,99003) js , kry , krx , krz , ytr(jb:jf,js)
99003       FORMAT (I3,2I10,I7,1X,30G13.4)
         ENDDO
         WRITE (nw,99004)
99004    FORMAT (85('-')/)
      ENDDO
 
      END SUBROUTINE wrslvp
