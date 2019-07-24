      subroutine S32_interp_phiSS ( a_coeff, b_coeff, period, specT, a, b, iflag )
      
      implicit none
      
      real period (17), a_coeff(17), b_coeff(17), specT, a, b
      integer iflag, count1, count2, nPer, i

C     First check for the PGA
      if (specT .le. 0.0) then 
        a = a_coeff(1)
        b = b_coeff(1)
        iflag = 0
        return
      endif
      
      nPer = 17
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 100 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
100   call S24_interp ( period(count1), period(count2), a_coeff(count1), a_coeff(count2),
     &    specT, a, iflag )
      call S24_interp ( period(count1), period(count2), b_coeff(count1), b_coeff(count2),
     &    specT, b, iflag )

      return
      end

c  --------------------------------------

      subroutine S32_SWUS_PHISS_CA1 ( mag, specT, phiSS, iflag, iBranch ) 
      
      implicit none
      
      integer iflag, iBranch
      real period (17), a_high(17), b_high(17), mag, specT, phiSS
      real a_low(17), b_low(17), a_central(17), b_central(17), a, b

c     Coeff from SWUS report, Table 7.3.3-1
      data period / 0.01, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 
     1              0.5, 0.75, 1, 1.5, 2, 3, 5, 10 /
      data a_Central /  0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485,
     1                0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 
     2                0.485, 0.485, 0.485 /
      data b_central / 0.3087, 0.3183, 0.3278, 0.3363, 0.3425, 0.3512, 
     1                0.3571, 0.3648, 0.3699, 0.3736, 0.3796, 0.3835, 
     2                0.3883, 0.3913, 0.3951, 0.396, 0.396 /
      data a_High / 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 
     1              0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 
     2              0.5789, 0.5789, 0.5789, 0.5789, 0.5789 /
      data b_High / 0.3685, 0.3799, 0.3913, 0.4014, 0.4088, 0.4192, 
     1              0.4262, 0.4354, 0.4415, 0.4459, 0.4531, 0.4577, 
     2              0.4635, 0.4671, 0.4716, 0.4727, 0.4727 /
      data a_Low / 0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 
     1              0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 
     2              0.3882, 0.3882, 0.3882, 0.3882, 0.3882 /
      data b_Low / 0.2471, 0.2548, 0.2624, 0.2692, 0.2741, 0.2811, 
     1             0.2858, 0.292, 0.2961, 0.299, 0.3038, 0.3069, 
     2             0.3108, 0.3132, 0.3162, 0.3169, 0.3169 /

c     Set the branch to use for the PhiSS_CA1  model
      if ( iBranch .eq. 1 ) then
        call S32_interp_phiSS ( a_low, b_low, period, specT, a, b, iflag ) 
      elseif ( iBranch .eq. 1 ) then
        call S32_interp_phiSS ( a_central, b_central, period, specT, a, b, iflag ) 
      elseif ( iBranch .eq. 3 ) then
        call S32_interp_phiSS ( a_high, b_high, period, specT, a, b,iflag )      
      endif

c     From SWUS eq 7.3.3-1a
      if ( mag .le. 7. ) then
        phiSS = a + (mag-5)/2. * (b-a)
      else
        phiSS = b
      endif
      return
      end



c --------------------

      subroutine S32_SWUS_PHISS_CA2 ( mag, specT, phiSS, iflag, iBranch ) 
      
      implicit none
      
      integer iflag, iBranch
      real period (17), a_high(17), c_high(17), mag, specT, phiSS
      real a_low(17), c_low(17), a_central(17), c_central(17), a, c

c     Coeff from SWUS report, Table 7.3.3-2
      data period / 0.01, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 
     1              0.5, 0.75, 1, 1.5, 2, 3, 5, 10 /
      data a_Central / 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 
     1                 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 
     2                 0.485, 0.485, 0.485 /
      data c_Central / 0.3581, 0.3654, 0.3725, 0.379, 0.3837, 0.3902, 
     1                 0.3947, 0.4006, 0.4044, 0.4072, 0.4118, 0.4147, 
     2                 0.4183, 0.4206, 0.4234, 0.4265, 0.4297 /
      data a_High/ 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 0.5789,
     1             0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 0.5789, 
     2             0.5789, 0.5789, 0.5789 /
      data c_High / 0.4273, 0.4357, 0.4452, 0.4524, 0.4583, 0.4655, 
     1              0.4715, 0.4786, 0.4822, 0.4858, 0.4918, 0.4953, 
     2              0.4989, 0.5025, 0.5049, 0.5085, 0.5133 /
      data a_Low / 0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 
     1             0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 0.3882, 
     2             0.3882, 0.3882, 0.3882, 0.3882, 0.3882 /
      data c_Low / 0.2865, 0.2921, 0.2985, 0.3033, 0.3073, 0.3121, 
     1             0.3161, 0.3209, 0.3233, 0.3257, 0.3297, 0.3321, 
     2             0.3346, 0.337, 0.3386, 0.341, 0.3442 /

c     Set the branch to use for the PhiSS_CA2  model
      if ( iBranch .eq. 1 ) then
        call S32_interp_phiSS ( a_low, c_low, period, specT, a, c, iflag ) 
      elseif ( iBranch .eq. 2 ) then
        call S32_interp_phiSS ( a_central, c_central, period, specT, a, c, iflag ) 
      elseif ( iBranch .eq. 3 ) then
        call S32_interp_phiSS ( a_high, c_high, period, specT, a, c, iflag )      
      endif

c     From SWUS eq 7.3.3-1ab
      if ( mag .le. 5.5 ) then
        phiSS = a + (mag-5)/0.5 * (c-a)
      else
        phiSS = c
      endif
      return
      end

c --------------------

      subroutine S32_SWUS_PHISS_Global_R50 ( phiSS, iflag, iBranch ) 

      implicit none
      
      integer iflag, iBranch
      real phiSS

c     from SWUS table 7.3.2-1 
c     Set the branch to use for the PhiSS_CA2  model
      if ( iBranch .eq. 1 ) then
        phiSS = 0.350
      elseif ( iBranch .eq. 2 ) then
        phiSS = 0.437
      elseif ( iBranch .eq. 3 ) then
        phiSS = 0.522
      endif
      iflag = 0

      return
      end
     
c --------------------------------

c     This is from an earlier verision of the code used for testing 
c     during SWUS.  It is not a final model.

C     Main subroutine for Single Station Sigma (Phi) models
      subroutine S32_SSSPhiModel (ssscalc1, specT, mag, Rrup, phiSSS) 

      implicit none
      
      real specT, mag, Rrup, phiSSS
      integer ssscalc1, modelflag, zoneflag

C     Select the Single Station Sigma model based on the SSSCalc1 value. 
      if (ssscalc1 .eq. 1 ) then
         call S32_PhiSS_PRP_1Mean ( specT, PhiSSS )
      elseif (ssscalc1 .eq. 2) then
         call S32_PhiSS_PRP_1MeanMinus ( specT, PhiSSS )
      elseif (ssscalc1 .eq. 3) then
         call S32_PhiSS_PRP_1MeanPlus ( specT, PhiSSS )
      elseif (ssscalc1 .eq. 4) then
         call S32_PhiSS_PRP_2Mean ( specT, Rrup, PhiSSS )
      elseif (ssscalc1 .eq. 5) then
         call S32_PhiSS_PRP_3Mean ( specT, Rrup, Mag, PhiSSS )

C     Preliminary DCPP Models
C     California Constant - Base
      elseif (ssscalc1 .eq. 10) then
         call S32_DCPP_CAconst ( specT, Rrup, Mag, PhiSSS )
C     California Constant - Lower 
      elseif (ssscalc1 .eq. 11) then
         call S32_DCPP_CAconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     California Constant - Upper 
      elseif (ssscalc1 .eq. 12) then
         call S32_DCPP_CAconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05

C     Global Constant - Base
      elseif (ssscalc1 .eq. 13) then
         call S32_DCPP_GBconst ( specT, Rrup, Mag, PhiSSS )
C     Global Constant - Lower 
      elseif (ssscalc1 .eq. 14) then
         call S32_DCPP_GBconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     Global Constant - Upper 
      elseif (ssscalc1 .eq. 15) then
         call S32_DCPP_GBconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05

C     CA Mag Dep - Base
      elseif (ssscalc1 .eq. 16) then
         call S32_DCPP_CAMag ( specT, Rrup, Mag, PhiSSS )
C     CA Mag Dep - Lower 
      elseif (ssscalc1 .eq. 17) then
         call S32_DCPP_CAMag ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     CA Mag Dep - Upper 
      elseif (ssscalc1 .eq. 18) then
         call S32_DCPP_CAMag ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05



C     Preliminary PVNGS Models - Local Host source
C     Local Global Constant - Base
      elseif (ssscalc1 .eq. 20) then
         call S32_PVNGS_LGBconst ( specT, Rrup, Mag, PhiSSS )
C     Local Global Constant - Lower
      elseif (ssscalc1 .eq. 21) then
         call S32_PVNGS_LGBconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     Local Global Constant - Upper
      elseif (ssscalc1 .eq. 22) then
         call S32_PVNGS_LGBconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05
         
C     Preliminary PVNGS Models - Distant California Sources
C     Distant Global Constant - Base
      elseif (ssscalc1 .eq. 23) then
         call S32_PVNGS_DGBconst ( specT, Rrup, Mag, PhiSSS )
C     Distant Global Constant - Lower
      elseif (ssscalc1 .eq. 24) then
         call S32_PVNGS_DGBconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     Distant Global Constant - Upper
      elseif (ssscalc1 .eq. 25) then
         call S32_PVNGS_DGBconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05
         
C     Preliminary PVNGS Models - North Arizona Path Sources
C     Distant Global Constant - Base
      elseif (ssscalc1 .eq. 26) then
         call S32_PVNGS_NAZconst ( specT, Rrup, Mag, PhiSSS)
C     Distant Global Constant - Lower
      elseif (ssscalc1 .eq. 27) then
         call S32_PVNGS_NAZconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     Distant Global Constant - Upper
      elseif (ssscalc1 .eq. 28) then
         call S32_PVNGS_NAZconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05
         
C     Preliminary PVNGS Models - South Arizona Path Sources
C     Distant Global Constant - Base
      elseif (ssscalc1 .eq. 29) then
         call S32_PVNGS_SAZconst ( specT, Rrup, Mag, PhiSSS )
C     Distant Global Constant - Lower
      elseif (ssscalc1 .eq. 30) then
         call S32_PVNGS_SAZconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS - 0.05
C     Distant Global Constant - Upper
      elseif (ssscalc1 .eq. 31) then
         call S32_PVNGS_SAZconst ( specT, Rrup, Mag, PhiSSS )
         phiSSS = phiSSS + 0.05
         
C     Updated Single Station Sigma PhiSS Models (July, 2014)
C     DCPP Models
C     California Dataset Constant
C     Central Model
      elseif (ssscalc1 .eq. 40) then
         modelflag = 0
         call S32_DCPP_CAconst_2014 ( specT, modelflag, PhiSSS )
C     Low Model
      elseif (ssscalc1 .eq. 41) then
         modelflag = 1
         call S32_DCPP_CAconst_2014 ( specT, modelflag, PhiSSS )
C     High Model
      elseif (ssscalc1 .eq. 42) then
         modelflag = 2
         call S32_DCPP_CAconst_2014 ( specT, modelflag, PhiSSS )

C     Global Dataset Constant
C     Central Model
      elseif (ssscalc1 .eq. 43) then
         modelflag = 0
         call S32_DCPP_GBconst_2014 ( specT, modelflag, PhiSSS )
C     Low Model
      elseif (ssscalc1 .eq. 44) then
         modelflag = 1
         call S32_DCPP_GBconst_2014 ( specT, modelflag, PhiSSS )
C     High Model
      elseif (ssscalc1 .eq. 45) then
         modelflag = 2
         call S32_DCPP_GBconst_2014 ( specT, modelflag, PhiSSS )

C     California Dataset Magnitude dependent
C     Central Model
      elseif (ssscalc1 .eq. 46) then
         modelflag = 0
         call S32_DCPP_CAMag_2014 ( specT, mag, modelflag, PhiSSS )
C     Low Model
      elseif (ssscalc1 .eq. 47) then
         modelflag = 1
         call S32_DCPP_CAMag_2014 ( specT, mag, modelflag, PhiSSS )
C     High Model
      elseif (ssscalc1 .eq. 48) then
         modelflag = 2
         call S32_DCPP_CAMag_2014 ( specT, mag, modelflag, PhiSSS )

C     PVNGS Models
C     PhiSP-R, Arizona, Zone1: Central
      elseif (ssscalc1 .eq. 50) then
         modelflag = 0
         zoneflag = 1
         call S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSSS )
C     PhiSP-R, Arizona, Zone1: Low
      elseif (ssscalc1 .eq. 51) then
         modelflag = 1
         zoneflag = 1
         call S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSSS )
C     PhiSP-R, Arizona, Zone1: High
      elseif (ssscalc1 .eq. 52) then
         modelflag = 2
         zoneflag = 1
         call S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSSS )
C     PhiSP-R, Arizona, Zone2/3: Central
      elseif (ssscalc1 .eq. 53) then
         modelflag = 0
         zoneflag = 2
         call S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSSS )
C     PhiSP-R, Arizona, Zone2/3: Low
      elseif (ssscalc1 .eq. 54) then
         modelflag = 1
         zoneflag = 2
         call S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSSS )
C     PhiSP-R, Arizona, Zone2/3: High
      elseif (ssscalc1 .eq. 55) then
         modelflag = 2
         zoneflag = 2
         call S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSSS )

C     PhiSS, Global, Zone1,2,3: Central
      elseif (ssscalc1 .eq. 59) then
         modelflag = 0
         call S32_PVNGS_GB123_2014 ( specT, modelflag, PhiSSS )
C     PhiSS, Global, Zone1,2,3: Low
      elseif (ssscalc1 .eq. 60) then
         modelflag = 1
         call S32_PVNGS_GB123_2014 ( specT, modelflag, PhiSSS )
C     PhiSS, Global, Zone1,2,3: High
      elseif (ssscalc1 .eq. 61) then
         modelflag = 2
         call S32_PVNGS_GB123_2014 ( specT, modelflag, PhiSSS )

C     PhiSS, Global, Outside Zone1,2,3: Central
      elseif (ssscalc1 .eq. 62) then
         modelflag = 0
         call S32_PVNGS_PSSGB_2014 ( specT, modelflag, PhiSSS )
C     PhiSS, Global, Outside Zone1,2,3: Low
      elseif (ssscalc1 .eq. 63) then
         modelflag = 1
         call S32_PVNGS_PSSGB_2014 ( specT, modelflag, PhiSSS )
C     PhiSS, Global, Outside Zone1,2,3: High
      elseif (ssscalc1 .eq. 64) then
         modelflag = 2
         call S32_PVNGS_PSSGB_2014 ( specT, modelflag, PhiSSS )

C     PhiSS, European, Outside Zone1,2,3: Central
      elseif (ssscalc1 .eq. 65) then
         modelflag = 0
         call S32_PVNGS_PSSEuro_2014 ( specT, modelflag, PhiSSS )
C     PhiSS, European, Outside Zone1,2,3: Low
      elseif (ssscalc1 .eq. 66) then
         modelflag = 1
         call S32_PVNGS_PSSEuro_2014 ( specT, modelflag, PhiSSS )
C     PhiSS, European, Outside Zone1,2,3: High
      elseif (ssscalc1 .eq. 67) then
         modelflag = 2
         call S32_PVNGS_PSSEuro_2014 ( specT, modelflag, PhiSSS )

       
         
      else
         write (*,*) 'Invalid Single Station Sigma Model.'
         write (*,*) 'SSSCalc1 = ', ssscalc1
         write (*,*) 'Check input file.'
         stop 99
      endif

      return 
      end

C--------------------------------------------------------------------------------------   

c     SIGMA.f
c     File added by Linda on July 26, 2011
c     To add different single-station sigma models

      subroutine S32_PhiSS_PRP_1Mean ( specT, PhiSS )
      
c     This subroutine evaluates the PEGASOS constant PhiSS model (mean values)
c     By Linda
c     July 26, 2011

      implicit none 
      integer MAXPER
      parameter (MAXPER=7) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1 
      
      data period  / 0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.46, 0.45, 0.48, 0.48, 0.46, 0.45, 0.41 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 7
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Pegasos Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end subroutine S32_PhiSS_PRP_1Mean
C--------------------------------------------------------------------------------------   

      subroutine S32_PhiSS_PRP_1MeanMinus ( specT, PhiSS )
      
c     This subroutine evaluates the PEGASOS constant PhiSS model (mean PhiSS minus
c     1.6* stdDev values)
c     By Linda
c     July 26, 2011

      implicit none 
      integer MAXPER
      parameter (MAXPER=7) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, std_PhiSS1(MAXPER), std_PhiSS1T
      
      data period  / 0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.46, 0.45, 0.48, 0.48, 0.46, 0.45, 0.41 /
      data std_PhiSS1 / 0.08, 0.08, 0.11, 0.09, 0.08, 0.07, 0.07 /
          
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        std_PhiSS1T = std_PhiSS1(i1)
        goto 5
      endif
      
      nPer = 7
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Pegasos Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )
      call S24_interp ( period(count1), period(count2), std_PhiSS1(count1), std_PhiSS1(count2),
     &    specT, std_PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T - 1.6*std_PhiSS1T
      
      return
      end subroutine S32_PhiSS_PRP_1MeanMinus

C--------------------------------------------------------------------------------------   

      subroutine S32_PhiSS_PRP_1MeanPlus ( specT, PhiSS )
      
c     This subroutine evaluates the PEGASOS constant PhiSS model (mean PhiSS plus
c     1.6* stdDev values)
c     By Linda
c     July 26, 2011

      implicit none 
      integer MAXPER
      parameter (MAXPER=7) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, std_PhiSS1(MAXPER), std_PhiSS1T
      
      data period  / 0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.46, 0.45, 0.48, 0.48, 0.46, 0.45, 0.41 /
      data std_PhiSS1 / 0.08, 0.08, 0.11, 0.09, 0.08, 0.07, 0.07 /
          
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        std_PhiSS1T = std_PhiSS1(i1)
        goto 5
      endif
      
      nPer = 7
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Pegasos Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )
      call S24_interp ( period(count1), period(count2), std_PhiSS1(count1), std_PhiSS1(count2),
     &    specT, std_PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T + 1.6*std_PhiSS1T
      
      return
      end subroutine S32_PhiSS_PRP_1MeanPlus
C--------------------------------------------------------------------------------------    

      subroutine S32_PhiSS_PRP_2Mean ( specT, Rrup, PhiSS )
      
c     This subroutine evaluates the PEGASOS distance-dependent 
c     PhiSS model (mean values)
c     By Linda
c     July 26, 2011

      implicit none 
      integer MAXPER
      parameter (MAXPER=7) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, PhiSS2(MAXPER), Rc1(MAXPER), Rc2(MAXPER), PhiSS2T,
     &  Rc1T, Rc2T, Rrup
      
      data period  / 0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.56, 0.55, 0.62, 0.62, 0.58, 0.54, 0.53 /
      data PhiSS2 / 0.45, 0.44, 0.47, 0.47, 0.45, 0.44, 0.40 /
      data Rc1 / 16, 16, 16, 16, 16, 16, 16 /
      data Rc2 / 32, 32, 32, 32, 32, 32, 36 /


C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        PhiSS2T = PhiSS2(i1)
        Rc1T = Rc1(i1)
        Rc2T = Rc2(i1)
        goto 5
      endif
      
      nPer = 7
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Pegasos Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )
      call S24_interp ( period(count1), period(count2), PhiSS2(count1), PhiSS2(count2),
     &    specT, PhiSS2T, iflag )
      call S24_interp ( period(count1), period(count2), Rc1(count1), Rc1(count2),
     &    specT, Rc1T, iflag )
      call S24_interp ( period(count1), period(count2), Rc2(count1), Rc2(count2),
     &    specT, Rc2T, iflag )

5     period1 = specT
      if ( Rrup < Rc1T ) then
        PhiSS = PhiSS1T
      elseif ( Rrup <= Rc2T ) then
        PhiSS = PhiSS1T + (PhiSS2T-PhiSS1T) * (Rrup-Rc1T) / (Rc2T-Rc1T)
      else
        PhiSS = PhiSS2T
      endif        
      
      return
      end subroutine S32_PhiSS_PRP_2Mean
C--------------------------------------------------------------------------------------   

      subroutine S32_PhiSS_PRP_3Mean ( specT, Rrup, Mag, PhiSS )
      
c     This subroutine evaluates the PEGASOS distance- and magnitude-dependent 
c     PhiSS model (mean values)
c     By Linda
c     July 26, 2011

      implicit none 
      integer MAXPER
      parameter (MAXPER=7) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS11(MAXPER), PhiSS, PhiSS11T,
     &  period1, PhiSS21(MAXPER), Rc11(MAXPER), Rc21(MAXPER), PhiSS21T,
     &  Rc11T, Rc21T, Rrup, C2(MAXPER), C2T, Mc1(MAXPER), Mc2(MAXPER),
     &  Mc1T, Mc2T, C1T, mag
      
      data period  / 0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0 / 
      data PhiSS11 / 0.58, 0.54, 0.60, 0.63, 0.59, 0.54, 0.44 /
      data PhiSS21 / 0.47, 0.44, 0.49, 0.50, 0.48, 0.45, 0.37 /
      data Rc11 / 16, 16, 16, 16, 16, 16, 16 /
      data Rc21 / 36, 36, 36, 36, 36, 36, 36 /
      data C2 / 0.34, 0.43, 0.37, 0.36, 0.36, 0.37, 0.37 /
      data Mc1 / 5.2, 5.2, 5.2, 5.2, 5.2, 5.3, 5.5 /
      data Mc2 / 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0 /


C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS11T = PhiSS11(i1)
        PhiSS21T = PhiSS21(i1)
        Rc11T = Rc11(i1)
        Rc21T = Rc21(i1)
        C2T = C2(i1)
        Mc1T = Mc1(i1)
        Mc2T = Mc2(i1)
        goto 5
      endif
      
      nPer = 7
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'Pegasos Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS11(count1), PhiSS11(count2),
     &    specT, PhiSS11T, iflag )
      call S24_interp ( period(count1), period(count2), PhiSS21(count1), PhiSS21(count2),
     &    specT, PhiSS21T, iflag )
      call S24_interp ( period(count1), period(count2), Rc11(count1), Rc11(count2),
     &    specT, Rc11T, iflag )
      call S24_interp ( period(count1), period(count2), Rc21(count1), Rc21(count2),
     &    specT, Rc21T, iflag )
      call S24_interp ( period(count1), period(count2), c2(count1), c2(count2),
     &    specT, c2T, iflag )
      call S24_interp ( period(count1), period(count2), Mc1(count1), Mc1(count2),
     &    specT, Mc1T, iflag )
      call S24_interp ( period(count1), period(count2), Mc2(count1), Mc2(count2),
     &    specT, Mc2T, iflag )

5     period1 = specT
      if ( Rrup < Rc11T ) then
        C1T = PhiSS11T
      elseif ( Rrup <= Rc21T ) then
        C1T = PhiSS11T + (PhiSS21T-PhiSS11T) * (Rrup-Rc11T) / (Rc21T-Rc11T)
      else
        C1T = PhiSS21T
      endif    
      
      if ( Mag < Mc1T ) then
        PhiSS = C1T
      elseif ( Mag <= Mc2T ) then
        PhiSS = C1T + (C2T-C1T) * (Mag-Mc1T) / (Mc2T-Mc1T)
      else
        PhiSS = C2T
      endif    
               
      return
      end subroutine S32_PhiSS_PRP_3Mean


C--------------------------------------------------------------------------------------  
C     Preliminary PhiSS Model - DCPP California Constant
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_DCPP_CAConst ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=5) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, Rrup, Mag 
      
      data period  / 0.01, 0.1, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.3621, 0.3657, 0.4260, 0.4326, 0.4242 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'DCPP California Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     Preliminary PhiSS Model - DCPP Global Constant
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_DCPP_GBConst ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=5) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, Rrup, Mag 
      
      data period  / 0.01, 0.1, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.3889, 0.4164, 0.4344, 0.4593, 0.4671 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'DCPP Global Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     Preliminary PhiSS Model - DCPP California magnitude Dependent
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_DCPP_CAMag ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=5) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS, a(MAXPER), b(MAXPER), c(MAXPER),
     &  period1, Rrup, Mag, aT, bT, cT
      
      data period  / 0.01, 0.1, 0.5, 1.0, 3.0 / 
      data a / -0.06297, -0.06954, -0.03137, -0.03816, -0.1146 /
      data b / 0.7598, 0.7998, 0.6138, 0.66, 1.153 /
      data c / 0.31901, 0.31302, 0.39421, 0.39288, 0.3508 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        aT = a(i1)
        bT = b(i1)
        cT = c(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'DCPP California Magnitude PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), a(count1), a(count2),
     &    specT, aT, iflag )
      call S24_interp ( period(count1), period(count2), b(count1), b(count2),
     &    specT, bT, iflag )
      call S24_interp ( period(count1), period(count2), c(count1), c(count2),
     &    specT, cT, iflag )

5     period1 = specT

C     Compute magnitude dependent PhiSS 
      if (mag .lt. 7.0) then
         phiSS = aT*mag + bT
      else
         phiSS = cT
      endif
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     Preliminary PhiSS Model - PVNGS Local Normal Global Constant
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_PVNGS_LGBConst ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=5) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, Rrup, Mag 
      
      data period  / 0.01, 0.1, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.4066, 0.4378, 0.4137, 0.4381, 0.4442 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS Local Source Global Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     Preliminary PhiSS Model - PVNGS Local Normal Global Constant
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_PVNGS_DGBConst ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=5) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, Rrup, Mag 
      
      data period  / 0.01, 0.1, 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.4975, 0.5069, 0.5248, 0.4579, 0.4871 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS Distant Source Global Constant PhiSS Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     Preliminary PhiSP Model - PVNGS Region-North
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_PVNGS_NAZConst ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=3) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, Rrup, Mag 
      
      data period  / 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.4704, 0.5221, 0.6250 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS Single Path North Constant PhiSP Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     Preliminary PhiSP Model - PVNGS Region-South
C     Models from Linda Al-Atik (3/2014)

      subroutine S32_PVNGS_SAZConst ( specT, Rrup, Mag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=3) 
      integer nPer, count1, count2, i, iflag, i1
      real specT, period(MAXPER), PhiSS1(MAXPER), PhiSS, PhiSS1T,
     &  period1, Rrup, Mag 
      
      data period  / 0.5, 1.0, 3.0 / 
      data PhiSS1 / 0.3517, 0.4712, 0.4736 /
     
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSS1T = PhiSS1(i1)
        goto 5
      endif
      
      nPer = 5
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS Single Path South Constant PhiSP Model'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSS1(count1), PhiSS1(count2),
     &    specT, PhiSS1T, iflag )

5     period1 = specT
      PhiSS = PhiSS1T
      
      return
      end 





C--------------------------------------------------------------------------------------  
C     PhiSS Model - DCPP California Constant
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper

      subroutine S32_DCPP_CAConst_2014 ( specT, modelflag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER)
      real PhiSS, PhiSSCT, phiSSHT, PhiSSLT, period1
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data PhissC  /  0.358, 0.358, 0.358, 0.365, 0.373, 0.379, 0.384, 0.39, 0.395, 
     1          0.397, 0.401, 0.404, 0.407, 0.412, 0.415, 0.418, 0.421, 0.423, 0.425, 
     1          0.426, 0.428, 0.43  / 
      data PhiSSH  /  0.427, 0.427, 0.427, 0.436, 0.445, 0.452, 0.458, 0.466, 0.472, 
     1          0.475, 0.479, 0.482, 0.486, 0.492, 0.495, 0.499, 0.503, 0.505, 0.506, 
     1          0.509, 0.51, 0.513  / 
      data PhiSSL  /  0.287, 0.287, 0.287, 0.292, 0.299, 0.303, 0.307, 0.312, 0.316, 
     1          0.318, 0.321, 0.323, 0.326, 0.33, 0.332, 0.335, 0.337, 0.339, 0.34, 
     1          0.341, 0.342, 0.344  /       
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'DCPP California Constant PhiSS Model - July 2014'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT
      if (modelflag .eq. 0) then
         PhiSS = PhiSSCT
      elseif (modelflag .eq. 1) then
         PhiSS = PhiSSLT
      elseif (modelflag .eq. 2) then
         PhiSS = PhiSSHT
      endif 
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     PhiSS Model - DCPP Global Constant
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper

      subroutine S32_DCPP_GBConst_2014 ( specT, modelflag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER)
      real PhiSS, PhiSSCT, PhiSSHT, PhiSSLT, period1
      
      data Period  /  0.01, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
     1       0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data PhissC  /  0.417, 0.417, 0.417, 0.419, 0.422, 0.425, 0.427, 0.429, 0.431, 
     1       0.432, 0.433, 0.435, 0.436, 0.438, 0.439, 0.44, 0.441, 0.442, 0.443, 
     1       0.443, 0.444, 0.445  / 
      data PhiSSH  /  0.498, 0.498, 0.498, 0.5, 0.504, 0.507, 0.51, 0.512, 0.514, 
     1       0.515, 0.517, 0.519, 0.52, 0.523, 0.524, 0.525, 0.526, 0.528, 0.528, 
     1       0.529, 0.53, 0.531  / 
      data PhiSSL  /  0.334, 0.334, 0.334, 0.335, 0.338, 0.34, 0.342, 0.343, 0.345, 
     1       0.346, 0.347, 0.348, 0.349, 0.351, 0.351, 0.352, 0.353, 0.354, 0.354, 
     1       0.355, 0.355, 0.356  /      
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'DCPP Global Constant PhiSS Model - July 2014'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT
      if (modelflag .eq. 0) then
         PhiSS = PhiSSCT
      elseif (modelflag .eq. 1) then 
         PhiSS = PhiSSLT
      elseif (modelflag .eq. 2) then 
         PhiSS = PhiSSHT
      endif
  
      return
      end 


C--------------------------------------------------------------------------------------  
C     PhiSS Model - DCPP California Magnitude Dependent
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper

      subroutine S32_DCPP_CAMag_2014 ( specT, mag, modelflag, PhiSS )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag
      real specT, period(MAXPER), PhiSS,
     &  period1, aC(MAXPER), aH(MAXPER), aL(MAXPER), bC(MAXPER), bH(MAXPER), bL(MAXPER),
     &  aCT, aHT, aLT, bCT, bHT, bLT, mag, aT, bT
      
      data period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1      0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data aC  /  0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485,
     1      0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485, 0.485  / 
      data bC  /  0.309, 0.309, 0.309, 0.318, 0.328, 0.336, 0.343, 0.351, 0.357, 0.361, 0.365,
     1      0.37, 0.374, 0.38, 0.383, 0.388, 0.391, 0.395, 0.395, 0.396, 0.391, 0.385  / 
      data aH  /  0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579,
     1      0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579, 0.579  / 
      data bH  /  0.369, 0.369, 0.369, 0.38, 0.391, 0.401, 0.409, 0.419, 0.426, 0.43, 0.435, 0.442,
     1      0.446, 0.453, 0.458, 0.464, 0.467, 0.472, 0.472, 0.473, 0.467, 0.46  / 
      data aL  /  0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 
     1      0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388  / 
      data bL  /  0.247, 0.247, 0.247, 0.255, 0.262, 0.269, 0.274, 0.281, 0.286, 0.289, 0.292, 0.296, 
     1      0.299, 0.304, 0.307, 0.311, 0.313, 0.316, 0.317, 0.317, 0.313, 0.308  / 

C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        aCT = aC(1)
        bCT = bC(1)
        aHT = aH(1)
        bHT = bH(1)
        aLT = aL(1)
        bLT = bL(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'DCPP California Mag Dependent PhiSS Model - July 2014'
      write (*,*) 'PhiSS Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), aC(count1), aC(count2),
     &    specT, aCT, iflag )
      call S24_interp ( period(count1), period(count2), bC(count1), bC(count2),
     &    specT, bCT, iflag )

      call S24_interp ( period(count1), period(count2), aH(count1), aH(count2),
     &    specT, aHT, iflag )
      call S24_interp ( period(count1), period(count2), bH(count1), bH(count2),
     &    specT, bHT, iflag )

      call S24_interp ( period(count1), period(count2), aL(count1), aL(count2),
     &    specT, aLT, iflag )
      call S24_interp ( period(count1), period(count2), bL(count1), bL(count2),
     &    specT, bLT, iflag )

5     period1 = specT

C     Set coefficients based on selected model.
      if (modelflag .eq. 0) then
          aT = aCT
          bT = bCT
      elseif (modelflag .eq. 1) then
          aT = aLT
          bT = bLT
      elseif (modelflag .eq. 2) then
          aT = aHT
          bT = bHT
      endif
      
      if (mag .lt. 7.0) then
         PhiSS = aT + ((Mag-5.0)/2.0)*(bT-aT)
      else
         phiSS = bT
      endif
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     PhiSS Model - PVNGS Outside Zones 1,2,3: European Data
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper

      subroutine S32_PVNGS_PSSEuro_2014 ( specT, modelflag, PhiSS )

      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER),
     1     PhiSS, PhiSSCT, PhiSSHT, PhiSSLT, period1
      
      data period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1         0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data PhiSSC  /  0.46, 0.46, 0.46, 0.464, 0.468, 0.471, 0.474, 0.477, 0.48, 0.481, 0.483, 
     1         0.485, 0.487, 0.489, 0.491, 0.493, 0.494, 0.496, 0.496, 0.497, 0.498, 0.499  / 
      data PhiSSH  /  0.549, 0.549, 0.549, 0.554, 0.559, 0.562, 0.566, 0.569, 0.573, 0.575, 0.577,
     1         0.579, 0.581, 0.584, 0.586, 0.588, 0.59, 0.592, 0.593, 0.593, 0.594, 0.596  / 
      data PhiSSL  /  0.368, 0.368, 0.368, 0.371, 0.375, 0.377, 0.379, 0.382, 0.384, 0.385, 0.387,
     1         0.388, 0.39, 0.391, 0.393, 0.395, 0.395, 0.397, 0.397, 0.398, 0.398, 0.399  / 
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS PhiSS European, Outside Zones 1,2,3 Model - July 2014'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT

      if (modelflag .eq. 0) then
         PhiSS = PhiSSCT
      elseif (modelflag .eq. 1) then 
         PhiSS = PhiSSLT
      elseif (modelflag .eq. 2) then 
         PhiSS = PhiSSHT
      endif

      return
      end 

C--------------------------------------------------------------------------------------  
C     PhiSS Model - PVNGS Outside Zones 1,2,3: Global Data
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper

      subroutine S32_PVNGS_PSSGB_2014 ( specT, modelflag, PhiSS )

      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER), PhiSS, 
     1     PhiSSCT, PhiSSHT, PhiSSLT, period1
      
      data period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data PhiSSC  /  0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 
     1          0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437, 0.437  / 
      data PhiSSH  /  0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 
     1          0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522  / 
      data PhiSSL  /  0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35,
     1          0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35  / 
      
C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS PhiSS Global, Outside Zones 1,2,3 Model - July 2014'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT

      if (modelflag .eq. 0) then
         PhiSS = PhiSSCT
      elseif (modelflag .eq. 1) then 
         PhiSS = PhiSSLT
      elseif (modelflag .eq. 2) then 
         PhiSS = PhiSSHT
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     PhiSPR Model - PVNGS Zones 1,2,3: Arizona Data
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper
C     Zoneflag = 1, 2 (which applies for 3 also)

      subroutine S32_PVNGS_SPRAZ_2014 ( specT, modelflag, zoneflag, PhiSS )

      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag, zoneflag
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER), 
     1     PhiSS, PhiSSCT, PhiSSHT, PhiSSLT, period1
      
      data period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 
     1        1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
c      data PhiSSC  /  0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.327, 0.335, 0.34, 0.344, 0.348, 
c     1        0.354, 0.358, 0.366, 0.371, 0.378, 0.382, 0.388, 0.391, 0.395, 0.398, 0.402  / 
c      data PhiSSH  /  0.417, 0.417, 0.417, 0.417, 0.417, 0.417, 0.426, 0.436, 0.443, 0.448, 0.454, 
c     1        0.461, 0.467, 0.476, 0.483, 0.492, 0.498, 0.505, 0.509, 0.514, 0.518, 0.524  / 
c      data PhiSSL  /  0.219, 0.219, 0.219, 0.219, 0.219, 0.219, 0.223, 0.229, 0.233, 0.235, 0.238, 0.242, 
c     1        0.245, 0.25, 0.253, 0.258, 0.261, 0.265, 0.267, 0.27, 0.272, 0.275  / 
C Updated 8/26/14
      data PhiSSC /  0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304,
     1               0.304, 0.304, 0.339, 0.364, 0.399, 0.423, 0.423, 0.423, 0.423, 0.423, 0.423 /
      data PhiSSH /  0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388,
     1               0.388, 0.388, 0.432, 0.464, 0.508, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54 /
      data PhiSSL /  0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 0.217, 
     1               0.217, 0.217, 0.242, 0.26, 0.284, 0.302, 0.302, 0.302, 0.302, 0.302, 0.302 /

C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS PhiSP-R, Zones 1,2,3 Model - July 2014'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT

      if (modelflag .eq. 0) then
         PhiSS = PhiSSCT
      elseif (modelflag .eq. 1) then 
         PhiSS = PhiSSLT
      elseif (modelflag .eq. 2) then 
         PhiSS = PhiSSHT
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     PhiSS Model - PVNGS Zones 1,2,3: Global Data
C     July 2014
C     Modelflag: 0 = Central, 1 = Lower, 2 = Upper

      subroutine S32_PVNGS_GB123_2014 ( specT, modelflag, PhiSS )

      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag, i1, modelflag
      real specT, period(MAXPER), PhiSSC(MAXPER), PhiSSH(MAXPER), PhiSSL(MAXPER),
     1     PhiSS, PhiSSCT, PhiSSHT, PhiSSLT, period1

      data period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 
     1         0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data PhissC  /  0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.515, 0.508, 0.498, 
     1         0.49, 0.476, 0.466, 0.452, 0.443, 0.431, 0.417, 0.399, 0.517, 0.517, 0.517, 0.517 / 
      data PhissH  /  0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.617, 0.615, 0.606, 0.594, 0.585, 
     1         0.568, 0.556, 0.54, 0.529, 0.514, 0.498, 0.476, 0.617, 0.617, 0.617, 0.617  / 
      data PhissL  /  0.414, 0.414, 0.414, 0.414, 0.414, 0.414, 0.414, 0.412, 0.407, 0.399, 0.392,
     1         0.381, 0.373, 0.362, 0.355, 0.345, 0.334, 0.319, 0.414, 0.414, 0.414, 0.414  / 

C     First check for the PGA
      if (specT .le. 0.0) then 
        if ( specT .eq. 0.0 ) i1=1
        period1 = period(i1)
        PhiSSCT = PhiSSC(i1)
        PhiSSHT = PhiSSH(i1)
        PhiSSLT = PhiSSL(i1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the PhiSS Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'PVNGS PhiSS Global, Zones 1,2,3 Model - July 2014'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), PhiSSC(count1), PhiSSC(count2),
     &    specT, PhiSSCT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSH(count1), PhiSSH(count2),
     &    specT, PhiSSHT, iflag )
      call S24_interp ( period(count1), period(count2), PhiSSL(count1), PhiSSL(count2),
     &    specT, PhiSSLT, iflag )

5     period1 = specT

      if (modelflag .eq. 0) then
         PhiSS = PhiSSCT
      elseif (modelflag .eq. 1) then 
         PhiSS = PhiSSLT
      elseif (modelflag .eq. 2) then 
         PhiSS = PhiSSHT
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - DCPP Central 
C     Scalc = 13001

      subroutine S32_SWUS_Sigma_DCPP_Cen ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER), sig2(MAXPER)
      real sig1T, sig2T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.557, 0.557, 0.559, 0.560, 0.562, 0.564, 0.566, 0.568, 0.569, 0.571, 
     1            0.572, 0.573, 0.574, 0.575, 0.576, 0.578, 0.578, 0.579, 0.580, 0.581,
     1            0.581, 0.582  /
      data sig2 / 0.491, 0.491, 0.494, 0.496, 0.500, 0.505, 0.508, 0.512, 0.515, 0.517,
     1            0.519, 0.522, 0.524, 0.527, 0.529, 0.532, 0.533, 0.535, 0.536, 0.537, 
     1            0.537, 0.538  /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        sig2T = sig2(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: DCPP Central'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )
      call S24_interp ( period(count1), period(count2), sig2(count1), sig2(count2),
     &    specT, sig2T, iflag )

5     period1 = specT

C     Now compute the magnitude-dependent total sigma value. 
      if (mag .ge. 7.0) then
         sigma = sig2T
      else
         sigma = sig1T + ((mag-5.0)/2.0)*(sig2T-sig1T)
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - DCPP Low 
C     Scalc = 13002

      subroutine S32_SWUS_Sigma_DCPP_Low ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER), sig2(MAXPER)
      real sig1T, sig2T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.428, 0.428, 0.431, 0.432, 0.436, 0.440, 0.443, 0.446, 0.449, 0.450, 
     1            0.452, 0.454, 0.455, 0.457, 0.459, 0.461, 0.462, 0.463, 0.464, 0.464,
     1            0.465, 0.466 /
      data sig2 / 0.389, 0.389, 0.393, 0.395, 0.401, 0.406, 0.410, 0.416, 0.419, 0.422,
     1            0.424, 0.427, 0.429, 0.432, 0.435, 0.437, 0.439, 0.441, 0.442, 0.442,
     1            0.443, 0.443 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        sig2T = sig2(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: DCPP Low'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )
      call S24_interp ( period(count1), period(count2), sig2(count1), sig2(count2),
     &    specT, sig2T, iflag )

5     period1 = specT

C     Now compute the magnitude-dependent total sigma value. 
      if (mag .ge. 7.0) then
         sigma = sig2T
      else
         sigma = sig1T + ((mag-5.0)/2.0)*(sig2T-sig1T)
      endif
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - DCPP High 
C     Scalc = 13003

      subroutine S32_SWUS_Sigma_DCPP_High ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER), sig2(MAXPER)
      real sig1T, sig2T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.694, 0.694, 0.695, 0.695, 0.696, 0.697, 0.697, 0.698, 0.699, 0.700, 
     1            0.700, 0.701, 0.701, 0.702, 0.703, 0.703, 0.704, 0.705, 0.705, 0.705, 
     1            0.706, 0.706 /
      data sig2 / 0.600, 0.600, 0.602, 0.604, 0.607, 0.61, 0.613, 0.616, 0.618, 0.620, 
     1            0.622, 0.624, 0.625, 0.628, 0.63, 0.632, 0.633, 0.635, 0.636, 0.637,
     1            0.638, 0.638 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        sig2T = sig2(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: DCPP High'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )
      call S24_interp ( period(count1), period(count2), sig2(count1), sig2(count2),
     &    specT, sig2T, iflag )

5     period1 = specT

C     Now compute the magnitude-dependent total sigma value. 
      if (mag .ge. 7.0) then
         sigma = sig2T
      else
         sigma = sig1T + ((mag-5.0)/2.0)*(sig2T-sig1T)
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS Arizona Sources Central 
C     Scalc = 13004

      subroutine S32_SWUS_Sigma_PVNGS_AZ_Cen ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER), sig2(MAXPER)
      real sig1T, sig2T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.573, 0.573, 0.574, 0.574, 0.575, 0.576, 0.576, 0.577, 0.577, 0.578, 
     1            0.578, 0.578, 0.579, 0.579, 0.579, 0.580, 0.580, 0.580, 0.580, 0.580,
     1            0.581, 0.581 /
      data sig2 / 0.553, 0.553, 0.553, 0.553, 0.554, 0.555, 0.555, 0.556, 0.556, 0.557, 
     1            0.557, 0.557, 0.558, 0.558, 0.558, 0.559, 0.559, 0.559, 0.559, 0.559,
     1            0.560, 0.560 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        sig2T = sig2(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS Arizona Sources Central'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )
      call S24_interp ( period(count1), period(count2), sig2(count1), sig2(count2),
     &    specT, sig2T, iflag )

5     period1 = specT

C     Now compute the magnitude-dependent total sigma value. 
      if (mag .ge. 7.0) then
         sigma = sig2T
      else
         sigma = sig1T + ((mag-5.0)/2.0)*(sig2T-sig1T)
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS Arizona Sources Low 
C     Scalc = 13005

      subroutine S32_SWUS_Sigma_PVNGS_AZ_Low ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER), sig2(MAXPER)
      real sig1T, sig2T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.461, 0.461, 0.461, 0.461, 0.462, 0.462, 0.463, 0.463, 0.463, 0.463, 
     1            0.463, 0.464, 0.464, 0.464, 0.464, 0.464, 0.464, 0.464, 0.464, 0.464, 
     1            0.464, 0.464 /
      data sig2 / 0.459, 0.459, 0.459, 0.460, 0.460, 0.460, 0.461, 0.461, 0.461, 0.461, 
     1            0.461, 0.461, 0.461, 0.462, 0.462, 0.462, 0.462, 0.462, 0.462, 0.462,
     1            0.462, 0.462 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        sig2T = sig2(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS Arizona Sources Low'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )
      call S24_interp ( period(count1), period(count2), sig2(count1), sig2(count2),
     &    specT, sig2T, iflag )

5     period1 = specT

C     Now compute the magnitude-dependent total sigma value. 
      if (mag .ge. 7.0) then
         sigma = sig2T
      else
         sigma = sig1T + ((mag-5.0)/2.0)*(sig2T-sig1T)
      endif
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS Arizona Sources High 
C     Scalc = 13006

      subroutine S32_SWUS_Sigma_PVNGS_AZ_High ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER), sig2(MAXPER)
      real sig1T, sig2T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.694, 0.694, 0.695, 0.695, 0.696, 0.697, 0.698, 0.699, 0.700, 0.701, 
     1            0.701, 0.702, 0.702, 0.703, 0.704, 0.705, 0.705, 0.706, 0.706, 0.707, 
     1            0.707, 0.707 /
      data sig2 / 0.652, 0.652, 0.653, 0.653, 0.655, 0.656, 0.657, 0.658, 0.659, 0.660, 
     1            0.661, 0.662, 0.662, 0.664, 0.664, 0.665, 0.666, 0.667, 0.667, 0.667, 
     1            0.668, 0.668 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        sig2T = sig2(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS Arizona Sources High'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )
      call S24_interp ( period(count1), period(count2), sig2(count1), sig2(count2),
     &    specT, sig2T, iflag )

5     period1 = specT

C     Now compute the magnitude-dependent total sigma value. 
      if (mag .ge. 7.0) then
         sigma = sig2T
      else
         sigma = sig1T + ((mag-5.0)/2.0)*(sig2T-sig1T)
      endif
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS California Sources with Path Central 
C     Scalc = 13007

      subroutine S32_SWUS_Sigma_PVNGS_CAPath_Cen ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER)
      real sig1T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.449, 0.449, 0.449, 0.449, 0.449, 0.449, 0.449, 0.449, 0.449, 0.449,
     1            0.449, 0.449, 0.449, 0.473, 0.491, 0.517, 0.535, 0.535, 0.535, 0.535,
     1            0.535, 0.535 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS California Sources with Path Central'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )

5     period1 = specT

      sigma = sig1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS California Sources with Path Low 
C     Scalc = 13008

      subroutine S32_SWUS_Sigma_PVNGS_CAPath_Low ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER)
      real sig1T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.354, 0.354, 0.354, 0.354, 0.354, 0.354, 0.354, 0.354, 0.354, 0.354,
     1            0.354, 0.354, 0.354, 0.375, 0.390, 0.410, 0.425, 0.425, 0.425, 0.425,
     1            0.425, 0.425 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS California Sources with Path Low'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )

5     period1 = specT

      sigma = sig1T
      
      return
      end 


C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS California Sources with Path High 
C     Scalc = 13009

      subroutine S32_SWUS_Sigma_PVNGS_CAPath_High ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER)
      real sig1T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 0.552, 
     1            0.552, 0.552, 0.552, 0.579, 0.600, 0.631, 0.655, 0.655, 0.655, 0.655,
     1            0.655, 0.655 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS California Sources with Path High'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )

5     period1 = specT

      sigma = sig1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS California Sources no Path Central 
C     Scalc = 13010

      subroutine S32_SWUS_Sigma_PVNGS_CAnoPath_Cen ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER)
      real sig1T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.612, 0.609, 
     1            0.606, 0.598, 0.591, 0.579, 0.571, 0.561, 0.553, 0.544, 0.537, 0.532, 
     1            0.525, 0.519 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS California Sources no Path Central'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )

5     period1 = specT

      sigma = sig1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS California Sources no Path Low 
C     Scalc = 13011

      subroutine S32_SWUS_Sigma_PVNGS_CAnoPath_Low ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER)
      real sig1T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.512, 0.512, 0.512, 0.512, 0.512, 0.512, 0.512, 0.512, 0.511, 0.508, 
     1            0.506, 0.499, 0.493, 0.483, 0.476, 0.467, 0.460, 0.452, 0.446, 0.442,
     1            0.435, 0.430 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS California Sources no Path Low'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )

5     period1 = specT

      sigma = sig1T
      
      return
      end 

C--------------------------------------------------------------------------------------  
C     SWUS Total Sigma Model - PVNGS California Sources no Path High 
C     Scalc = 13012

      subroutine S32_SWUS_Sigma_PVNGS_CAnoPath_High ( mag, specT, sigma, iflag )
      
      implicit none 
      integer MAXPER
      parameter (MAXPER=22) 
      integer nPer, count1, count2, i, iflag
      real specT, period(MAXPER), sig1(MAXPER)
      real sig1T, period1, sigma, mag
      
      data Period  /  0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1          0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10  / 
      data sig1 / 0.720, 0.720, 0.720, 0.720, 0.720, 0.720, 0.720, 0.720, 0.718, 0.715, 
     1            0.712, 0.702, 0.694, 0.681, 0.672, 0.660, 0.652, 0.641, 0.634, 0.629, 
     1            0.620, 0.614 /

C     First check for the PGA
      if (specT .eq. 0.0) then 
        period1 = period(1)
        sig1T = sig1(1)
        goto 5
      endif
      
      nPer = 22
C     For other periods, loop over the spectral period range of the Sigma Model.
      do i = 2, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      
C     Selected spectral period is outside range defined by the model.
      write (*,*) 
      write (*,*) 'SWUS Total Sigma Model: PVNGS California Sources no Path High'
      write (*,*) 'Model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99 
      
C     Interpolate the coefficients for the requested spectral period.
1020  call S24_interp ( period(count1), period(count2), sig1(count1), sig1(count2),
     &    specT, sig1T, iflag )

5     period1 = specT

      sigma = sig1T
      
      return
      end 







