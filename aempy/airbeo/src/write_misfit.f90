!*==write_misfit.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE write_misfit(nw,itspr,ndata,tdfd,nfrq,freq,cmp,nchnl,  &
     &                        verr,xmodl,xdata)
!---------------------------------------------------------------------------------
 
!***  Called by: NLSQ2
 
!  Writes out the error misfit at any stage of the inversion
!
!     ITSPR > 0 => printout after ITSPR iterations on unit NW
!           < 0 => printout after final iterations on unit NW
!           = 0 => printout initial error structure unit NW
!
!     NDATA - total number of data points
!      TDFD = 1: time domain  or 2: frequency-domain
!      FREQ - array of NFRQ frequencies for FD inversion
!     XDATA - 1D array of NDATA measured data points
!     XMODL - 1D array of NDATA computed data points from most recent model
!      VERR - symmetric error at each data point
!       CMP = 11 => inversion on horizontal in-line component only
!           = 13 => inversion on vertical component only
!           = 2  => joint inversion on vertical & horizontal in-line components
!           = 3  => inversion on all 3 components
!           = 4 or 42 => inversion on total field
!
 
      INTEGER nw , itspr , ndata , tdfd , nfrq , cmp , nchnl , nchn ,   &
     &        n1 , n2
      REAL freq(nfrq)
      REAL , DIMENSION(ndata) :: xmodl , xdata , verr
      CHARACTER(LEN=8) chn(50)
      DATA chn(1:50)/' CHNL_1 ' , ' CHNL_2 ' , ' CHNL_3 ' , ' CHNL_4 ' ,&
     &     ' CHNL_5 ' , ' CHNL_6 ' , ' CHNL_7 ' , ' CHNL_8 ' ,          &
     &     ' CHNL_9 ' , 'CHNL_10 ' , 'CHNL_11 ' , 'CHNL_12 ' ,          &
     &     'CHNL_13 ' , 'CHNL_14 ' , 'CHNL_15 ' , 'CHNL_16 ' ,          &
     &     'CHNL_17 ' , 'CHNL_18 ' , 'CHNL_19 ' , 'CHNL_20 ' ,          &
     &     'CHNL_21 ' , 'CHNL_22 ' , 'CHNL_23 ' , 'CHNL_24 ' ,          &
     &     'CHNL_25 ' , 'CHNL_26 ' , 'CHNL_27 ' , 'CHNL_28 ' ,          &
     &     'CHNL_29 ' , 'CHNL_30 ' , 'CHNL_31 ' , 'CHNL_32 ' ,          &
     &     'CHNL_33 ' , 'CHNL_34 ' , 'CHNL_35 ' , 'CHNL_36 ' ,          &
     &     'CHNL_37 ' , 'CHNL_38 ' , 'CHNL_39 ' , 'CHNL_40 ' ,          &
     &     'CHNL_41 ' , 'CHNL_42 ' , 'CHNL_43 ' , 'CHNL_44 ' ,          &
     &     'CHNL_45 ' , 'CHNL_46 ' , 'CHNL_47 ' , 'CHNL_48 ' ,          &
     &     'CHNL_49 ' , 'CHNL_50 '/
 
!  Put data into matrix form
 
      IF ( itspr==0 ) THEN
         WRITE (nw,99001)
 
99001    FORMAT (//T7,'ERROR STRUCTURE OF INITIAL MODEL'/T7,            &
     &           '---------------------------------')
      ELSEIF ( itspr>0 ) THEN
         WRITE (nw,99002) itspr
99002    FORMAT (//T7,'ERROR STRUCTURE OF MODEL AFTER',I3,              &
     &           ' ITERATIONS'/T7,                                      &
     &           '--------------------------------------------')
      ELSEIF ( itspr<0 ) THEN
         WRITE (nw,99003)
99003    FORMAT (//T7,'ERROR STRUCTURE OF FINAL MODEL'/T7,              &
     &           '------------------------------')
      ENDIF
      WRITE (nw,99004)
99004 FORMAT (/T3,'For each station:'/T3,'----------------'/T5,         &
     &        'The first line is the percent symmetric error.'/T5,      &
     &        'The second line is the model response.'/T5,              &
     &        'The third line is the data.')
 
      IF ( tdfd==1 ) THEN
         nchn = min(nchnl,50)
         IF ( cmp==13 .OR. cmp==2 .OR. cmp==3 ) THEN
            WRITE (nw,99005)
99005       FORMAT (//T7,'VERTICAL COMPONENT STRUCTURE'/T7,             &
     &              '----------------------------')
            n1 = 1
            n2 = nchn
            CALL prt_td
         ENDIF
 
         IF ( cmp==11 .OR. cmp==2 .OR. cmp==3 ) THEN
            WRITE (nw,99006)
99006       FORMAT (//T7,'IN-LINE COMPONENT STRUCTURE'/T7,              &
     &              '---------------------------')
            n1 = 1
            n2 = nchn
            IF ( cmp==2 .OR. cmp==3 ) THEN
               n1 = nchnl + 1
               n2 = n1 - 1 + nchn
            ENDIF
            CALL prt_td
         ENDIF
 
         IF ( cmp==3 ) THEN
            WRITE (nw,99007)
99007       FORMAT (//T7,'TRANSVERSE COMPONENT STRUCTURE'/T7,           &
     &              '------------------------------')
            n1 = 2*nchnl + 1
            n2 = n1 - 1 + nchn
            CALL prt_td
         ENDIF
 
         IF ( cmp==4 .OR. cmp==42 ) THEN
            WRITE (nw,99008)
99008       FORMAT (//T7,'TOTAL COMPONENT STRUCTURE'/T7,                &
     &              '-------------------------')
            n1 = 1
            n2 = nchn
            CALL prt_td
         ENDIF
 
      ELSE
         WRITE (nw,99009)
99009    FORMAT (//T7,'IN-PHASE COMPONENT STRUCTURE'/T7,                &
     &           '----------------------------')
         n1 = 1
         n2 = nfrq
         CALL prt_fd
         WRITE (nw,99010)
99010    FORMAT (//T7,'QUADRATURE COMPONENT STRUCTURE'/T7,              &
     &           '------------------------------')
         n1 = nfrq + 1
         n2 = 2*nfrq
         CALL prt_fd
      ENDIF
 
      CONTAINS
 
      SUBROUTINE prt_td
!  -----------------
 
!***  Called by: WRITE_MISFIT
 
      WRITE (nw,'(/6X,30(A:,5X))') chn(1:nchn)
      WRITE (nw,'(/50G13.4)') 100.*verr(n1:n2)
      WRITE (nw,'( 50G13.4)') xmodl(n1:n2)
      WRITE (nw,'( 50G13.4)') xdata(n1:n2)
      END SUBROUTINE prt_td
 
 
      SUBROUTINE prt_fd
!  -----------------
 
!***  Called by: WRITE_MISFIT
 
      WRITE (nw,'(/40F11.0)') freq(1:nfrq)
      WRITE (nw,'(/40F11.2)') 100.*verr(n1:n2)
      WRITE (nw,'( 40F11.2)') xmodl(n1:n2)
      WRITE (nw,'( 40F11.2)') xdata(n1:n2)
      END SUBROUTINE prt_fd
 
      END SUBROUTINE write_misfit
