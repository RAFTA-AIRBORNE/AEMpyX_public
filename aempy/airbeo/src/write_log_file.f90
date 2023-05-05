!*==write_log_file.f90  processed by SPAG 6.55Rc at 14:57 on 13 Oct 2014
 
      SUBROUTINE write_log_file(nlg,msg,mxerr,err_lvl)
! -------------------------------------------------
 
!***  Called by: READ_SYSTEM_AND_SURVEY, READ_MODEL, READ_INVRT_CNTRL_AND_DATA
 
! This subroutine prints out warning and fatal error messages on the LOG file.
!
! NLG = Airbeo.log unit index
! MSG refers to error message index
! ERR_LVL = 1 for warnings;  = 2 for fatal errors
! MXERR = MAX (ERR_LVL)
 
      INTEGER err_lvl , msg , nlg , mxerr
 
      IF ( mxerr==0 ) OPEN (nlg,FILE='Airbeo.log',STATUS='REPLACE')
 
      mxerr = max(err_lvl,mxerr)
      IF ( err_lvl==1 ) WRITE (nlg,99001)
 
99001 FORMAT (/T2,'WARNING'/T2,'-------')
      IF ( err_lvl==2 ) WRITE (nlg,99002)
99002 FORMAT (/T2,'FATAL ERROR'/T2,'----- -----')
      IF ( msg==100 ) WRITE (nlg,'(2X)')
 
      IF ( msg==1 ) WRITE (nlg,99003)
 
99003 FORMAT (/T3,'The value for TDFD is outside the permitted range.'/ &
     &        T3,                                                       &
     &'The allowed values are: 1 for time-domain or 2 for frequency doma&
     &in.'/)
      IF ( msg==4 ) WRITE (nlg,99004)
99004 FORMAT (/T3,'The value for ISW is outside the permitted range.')
      IF ( msg==5 ) WRITE (nlg,99005)
99005 FORMAT (/T3,'The value for STEP is outside the permitted range.'/ &
     &        T3,'The allowed values are: 0 or 1.'/)
      IF ( msg==6 ) WRITE (nlg,99006)
99006 FORMAT (/T3,'The value for KRXW is outside the permitted range.'/ &
     &        T3,'The allowed values are: 1 or 2.'/)
      IF ( msg==7 ) WRITE (nlg,99007)
99007 FORMAT (/T3,'This value for TOPN is outside the permitted range.'/&
     &        T3,'It must be > 0.'/)
      IF ( msg==8 ) WRITE (nlg,99008)
99008 FORMAT (/T3,'For inversion CMP must be 11, 13, 2, 3 or 4.')
      IF ( msg==9 ) WRITE (nlg,99009)
99009 FORMAT (/T3,'CMP must be 11, 13, 2 or 3.  It has been reset to 3')
      IF ( msg==10 ) WRITE (nlg,99010)
99010 FORMAT (/T3,                                                      &
     &  'KPPM is outside allowed range.  It must be 0, 1, 3, 123, or 4.'&
     &  )
      IF ( msg==11 ) WRITE (nlg,99011)
99011 FORMAT (/T3,                                                      &
     &        'Input calibration and output must have the same units'/T3&
     &        ,'Output will be in dB/dt units')
      IF ( msg==12 ) WRITE (nlg,99012)
99012 FORMAT (/T3,                                                      &
     &        'Input calibration and output must have the same units'/T3&
     &        ,'Output will be in B units')
      IF ( msg==13 ) WRITE (nlg,99013)
99013 FORMAT (//T3,                                                     &
     &        'CMP must = 1, 2, or 3 for frequency domain modelling'/T3,&
     &        'CMP has been set to 1')
      IF ( msg==14 ) WRITE (nlg,99014)
99014 FORMAT (//T3,                                                     &
     &'Frequency-domain inversion is based on the maximally coupled comp&
     &onent only.'/T3,'CMP must = 1 or -1')
      IF ( msg==16 ) WRITE (nlg,99015)
99015 FORMAT (/T3,'SURVEY must be either 1, 2, or 3 for time-domain.'/T3&
     &        ,'SURVEY = 3 cannot be used for frequency-domain because'/&
     &        T3,                                                       &
     &        'Tx-Rx offset must be constant as a function of position')
      IF ( msg==17 ) WRITE (nlg,99016)
99016 FORMAT (/T3,'IUNITS must be 1, 2 or 3')
      IF ( msg==19 ) WRITE (nlg,99017)
99017 FORMAT (/T3,                                                      &
     &     'IUNITS must be 1, 2 or 3.  It has been reset to the default'&
     &     )
 
      IF ( msg==50 ) WRITE (nlg,99018)
 
!  Model Messages
 
99018 FORMAT (/T3,'QLYR must = 1 (invert on layer thickness)')
      IF ( msg==53 ) WRITE (nlg,99019)
99019 FORMAT (/T3,                                                      &
     &'A lithology must have a positive first component (resistivity) if&
     & it is to be'/T3,'applied to a layer.')
      IF ( msg==54 ) WRITE (nlg,99020)
99020 FORMAT (/T3,'LITH must be an integer between 1 & NLITH')
      IF ( msg==55 ) WRITE (nlg,99021)
99021 FORMAT (/T3,'Layer resistivities must be positive.')
 
      IF ( msg==201 ) WRITE (nlg,99022)
 
! Inversion messages
! -------------------
 
99022 FORMAT (/T3,                                                      &
     &'This version of Airbeo does not invert for AEM system parameters.&
     &')
      IF ( msg==202 ) WRITE (nlg,99023)
99023 FORMAT (/T3,                                                      &
     &'X component inversion was requested but Airbeo.inv contains only &
     &Z component data.')
      IF ( msg==203 ) WRITE (nlg,99024)
99024 FORMAT (/T3,                                                      &
     &'Z component inversion was requested but Airbeo.inv contains only &
     &X component data.')
      IF ( msg==204 ) WRITE (nlg,99025)
99025 FORMAT (/T3,                                                      &
     &        'The maximum number of iterations has been reduced to 20.'&
     &        )
      IF ( msg==205 ) WRITE (nlg,99026)
99026 FORMAT (/T3,'CNVRG must be 1 or 2.  It has been reset to 1.')
      IF ( msg==206 ) WRITE (nlg,99027)
99027 FORMAT (/T3,                                                      &
     &'For inversion, aircraft positions must be entered for every stati&
     &on',/T3,'The automatic course option, SURVEY = 1 is not allowed.',&
     &/T3,'SURVEY must be either 2, 3, -2, or -3.'/)
      IF ( msg==207 ) WRITE (nlg,99028)
99028 FORMAT (/T3,'KCMP must = 12 or 21 for Frequency-domain inversion.'&
     &        )
      IF ( msg==208 ) WRITE (nlg,99029)
99029 FORMAT (/T3,'ORDER must = 1 or 2.')
      IF ( msg==209 ) WRITE (nlg,99030)
99030 FORMAT (/T3,                                                      &
     &        'KCMP is restricted to the values: 1, 3, 13, 31, 123, 312'&
     &        )
      IF ( msg==210 ) WRITE (nlg,99031)
99031 FORMAT (/T3,                                                      &
     &'X & Z component inversion was requested but Airbeo.inv contains o&
     &nly X component data.')
      IF ( msg==211 ) WRITE (nlg,99032)
99032 FORMAT (/T3,                                                      &
     &'X & Z component inversion was requested but Airbeo.inv contains o&
     &nly Z component data.')
      IF ( msg==212 ) WRITE (nlg,99033)
99033 FORMAT (/T3,                                                      &
     &'3 component inversion was requested but Airbeo.inv contains only &
     &X component data.')
      IF ( msg==213 ) WRITE (nlg,99034)
99034 FORMAT (/T3,                                                      &
     &'3 component inversion was requested but Airbeo.inv contains only &
     &Z component data.')
      IF ( msg==214 ) WRITE (nlg,99035)
99035 FORMAT (/T3,                                                      &
     &'3 component inversion was requested but Airbeo.inv contains no Y &
     &component data.')
      IF ( msg==215 ) WRITE (nlg,99036)
99036 FORMAT (/T3,                                                      &
     & 'There is a component discrepency between the cfl and inv files.'&
     & )
      IF ( msg==220 ) WRITE (nlg,99037)
99037 FORMAT (/T3,                                                      &
     &'Order must be 1122 (for IIQQ),  1212 (for IQIQ),  2211 (for QQII)&
     &,  2121 (for QIQI)')
 
      END SUBROUTINE write_log_file
