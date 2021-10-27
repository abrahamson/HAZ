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
      data sig1 / 0.576, 0.576, 0.577, 0.577, 0.578, 0.578, 0.579, 0.58, 0.581,
     1            0.581, 0.581, 0.582, 0.582, 0.583, 0.583, 0.584, 0.584, 0.585,
     2            0.585, 0.585, 0.585, 0.586 /
      data sig2 / 0.495, 0.495, 0.498, 0.499, 0.504, 0.507, 0.51, 0.514, 0.517,
     1            0.519, 0.52, 0.522, 0.524, 0.527, 0.529, 0.531, 0.532, 0.534,
     2            0.534, 0.535, 0.535, 0.536 /

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
      data sig1 / 0.456, 0.456, 0.457, 0.458, 0.46, 0.461, 0.462, 0.464, 0.465, 0.465,
     1            0.466, 0.466, 0.467, 0.468, 0.468, 0.469, 0.469, 0.47, 0.47, 0.47,
     2            0.471, 0.471 /
      data sig2 / 0.39, 0.39, 0.394, 0.396, 0.402, 0.407, 0.411, 0.416, 0.419, 0.422,
     1            0.424, 0.427, 0.429, 0.432, 0.434, 0.437, 0.439, 0.441, 0.441,
     2            0.441, 0.442, 0.442 /

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
      data sig1 / 0.699, 0.699, 0.699, 0.7, 0.7, 0.701, 0.702, 0.703, 0.703, 0.704,
     1            0.704, 0.704, 0.705, 0.705, 0.706, 0.706, 0.707, 0.707, 0.707,
     2            0.707, 0.708, 0.708 /
      data sig2 / 0.614, 0.614, 0.614, 0.615, 0.616, 0.617, 0.618, 0.62, 0.621,
     1            0.622, 0.623, 0.625, 0.626, 0.628, 0.629, 0.631, 0.632, 0.633,
     2            0.634, 0.634, 0.635, 0.635 /

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

c ---------------------------------------------------------------------
      subroutine S32_NGAEast_CompErgSig_Low ( m, specT, sigma, iflag )

c     model: NGA-East composite ergodic sigma, CENA, low
c     ref: Goulet et al., 2018 (PEER Report 2018/08), equation 11-14 for functional form,
c          Table 11-24 for coefficients

      implicit none

      integer MAX_FREQ
      parameter (MAX_FREQ=24)
      integer iflag, ifreq, c1, c2
      real m, specT, specF, freq(MAX_FREQ), sig1(MAX_FREQ), sig2(MAX_FREQ),
     1     sig3(MAX_FREQ), sig4(MAX_FREQ), sig1F, sig2F, sig3F, sig4F, sigma

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1           2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2           20.000, 25.000, 33.333, 50.000, 100.000, 100.000, -1.0 /

      data sig1 / 0.4826, 0.4871, 0.4983, 0.5076, 0.5229, 0.5497, 0.5703, 0.5977, 0.6145,
     1           0.6323, 0.6392, 0.6452, 0.6476, 0.6496, 0.6521, 0.6645, 0.6820,
     2           0.6874, 0.6758, 0.6610, 0.6533, 0.6539, 0.6539, 0.6034 /

      data sig2 / 0.4583, 0.4630, 0.4747, 0.4844, 0.5004, 0.5282, 0.5496, 0.5780, 0.5954,
     1           0.6138, 0.6210, 0.6271, 0.6296, 0.6317, 0.6343, 0.6471, 0.6650,
     2           0.6705, 0.6586, 0.6434, 0.6355, 0.6361, 0.6361, 0.5980 /

      data sig3 / 0.4207, 0.4256, 0.4379, 0.4480, 0.4648, 0.4945, 0.5156, 0.5430, 0.5586,
     1           0.5723, 0.5760, 0.5767, 0.5761, 0.5750, 0.5743, 0.5867, 0.6060,
     2           0.6118, 0.5987, 0.5819, 0.5731, 0.5736, 0.5736, 0.5587 /

      data sig4 / 0.4153, 0.4205, 0.4337, 0.4446, 0.4624, 0.4926, 0.5080, 0.5227, 0.5299,
     1           0.5293, 0.5244, 0.5119, 0.5040, 0.4958, 0.4876, 0.4982, 0.5200,
     2           0.5264, 0.5113, 0.4913, 0.4805, 0.4807, 0.4807, 0.4936 /

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(23)
        sig1F = sig1(23)
        sig2F = sig2(23)
        sig3F = sig3(23)
        sig4F = sig4(23)
        goto 1021
      elseif (specT .eq. -1.0) then
        specF = freq(24)
        sig1F = sig1(24)
        sig2F = sig2(24)
        sig3F = sig3(24)
        sig4F = sig4(24)
        goto 1021
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        c1 = 0
        c2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            c1 = ifreq
            c2 = ifreq+1
            goto 1020
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     interpolate between frequencies
1020  call S24_interp1 ( freq(c1), freq(c2), sig1(c1), sig1(c2), specF, sig1F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig2(c1), sig2(c2), specF, sig2F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig3(c1), sig3(c2), specF, sig3F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig4(c1), sig4(c2), specF, sig4F, iflag )

1021  if (m .le. 4.5) then
        sigma = sig1F
      elseif (m .gt. 4.5 .and. m .le. 5.0) then
        sigma = sig1F + (sig2F - sig1F) * ((m - 4.5)/0.5)
      elseif (m .gt. 5.0 .and. m .le. 5.5) then
        sigma = sig2F + (sig3F - sig2F) * ((m - 5.0)/0.5)
      elseif (m .gt. 5.5 .and. m .le. 6.5) then
        sigma = sig3F + (sig4F - sig3F) * ((m - 5.5)/1.0)
      else
        sigma = sig4F
      endif

      return
      end

c ---------------------------------------------------------------------
      subroutine S32_NGAEast_CompErgSig_Cen ( m, specT, sigma, iflag )

c     model: NGA-East composite ergodic sigma, CENA, central
c     ref: Goulet et al., 2018 (PEER Report 2018/08), equation 11-14 for functional form,
c          Table 11-24 for coefficients

      implicit none

      integer MAX_FREQ
      parameter (MAX_FREQ=24)
      integer iflag, ifreq, c1, c2
      real m, specT, specF, freq(MAX_FREQ), sig1(MAX_FREQ), sig2(MAX_FREQ),
     1     sig3(MAX_FREQ), sig4(MAX_FREQ), sig1F, sig2F, sig3F, sig4F, sigma

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1           2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2           20.000, 25.000, 33.333, 50.000, 100.000, 100.000, -1.0 /

      data sig1 / 0.6016, 0.6053, 0.6139, 0.6212, 0.6337, 0.6562, 0.6742, 0.6992, 0.7151,
     1           0.7335, 0.7412, 0.7488, 0.7523, 0.7556, 0.7597, 0.7726, 0.7881,
     2           0.7920, 0.7810, 0.7676, 0.7610, 0.7618, 0.7618, 0.7045 /

      data sig2 / 0.5829, 0.5867, 0.5956, 0.6032, 0.6160, 0.6392, 0.6576, 0.6832, 0.6995,
     1           0.7182, 0.7261, 0.7338, 0.7374, 0.7408, 0.7450, 0.7581, 0.7739,
     2           0.7779, 0.7667, 0.7530, 0.7463, 0.7471, 0.7471, 0.6995 /

      data sig3 / 0.5535, 0.5576, 0.5670, 0.5749, 0.5885, 0.6127, 0.6299, 0.6528, 0.6661,
     1           0.6786, 0.6825, 0.6842, 0.6844, 0.6844, 0.6850, 0.6969, 0.7133,
     2           0.7169, 0.7044, 0.6891, 0.6814, 0.6819, 0.6819, 0.6569 /

      data sig4 / 0.5285, 0.5327, 0.5426, 0.5509, 0.5651, 0.5907, 0.6050, 0.6215, 0.6283,
     1           0.6274, 0.6231, 0.6129, 0.6065, 0.6001, 0.5941, 0.6039, 0.6214,
     2           0.6239, 0.6090, 0.5905, 0.5808, 0.5809, 0.5809, 0.5825 /

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(23)
        sig1F = sig1(23)
        sig2F = sig2(23)
        sig3F = sig3(23)
        sig4F = sig4(23)
        goto 1021
      elseif (specT .eq. -1.0) then
        specF = freq(24)
        sig1F = sig1(24)
        sig2F = sig2(24)
        sig3F = sig3(24)
        sig4F = sig4(24)
        goto 1021
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        c1 = 0
        c2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            c1 = ifreq
            c2 = ifreq+1
            goto 1020
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     interpolate between frequencies
1020  call S24_interp1 ( freq(c1), freq(c2), sig1(c1), sig1(c2), specF, sig1F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig2(c1), sig2(c2), specF, sig2F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig3(c1), sig3(c2), specF, sig3F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig4(c1), sig4(c2), specF, sig4F, iflag )

1021  if (m .le. 4.5) then
        sigma = sig1F
      elseif (m .gt. 4.5 .and. m .le. 5.0) then
        sigma = sig1F + (sig2F - sig1F) * ((m - 4.5)/0.5)
      elseif (m .gt. 5.0 .and. m .le. 5.5) then
        sigma = sig2F + (sig3F - sig2F) * ((m - 5.0)/0.5)
      elseif (m .gt. 5.5 .and. m .le. 6.5) then
        sigma = sig3F + (sig4F - sig3F) * ((m - 5.5)/1.0)
      else
        sigma = sig4F
      endif

      return
      end

c ---------------------------------------------------------------------
      subroutine S32_NGAEast_CompErgSig_High ( m, specT, sigma, iflag )

c     model: NGA-East composite ergodic sigma, CENA, high
c     ref: Goulet et al., 2018 (PEER Report 2018/08), equation 11-14 for functional form,
c          Table 11-24 for coefficients

      implicit none

      integer MAX_FREQ
      parameter (MAX_FREQ=24)
      integer iflag, ifreq, c1, c2
      real m, specT, specF, freq(MAX_FREQ), sig1(MAX_FREQ), sig2(MAX_FREQ),
     1     sig3(MAX_FREQ), sig4(MAX_FREQ), sig1F, sig2F, sig3F, sig4F, sigma

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1           2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2           20.000, 25.000, 33.333, 50.000, 100.000, 100.000, -1.0 /

      data sig1 / 0.7519, 0.7545, 0.7591, 0.7627, 0.7687, 0.7807, 0.7916, 0.8096, 0.8228,
     1           0.8402, 0.8485, 0.8575, 0.8621, 0.8668, 0.8725, 0.8858, 0.8990,
     2           0.9011, 0.8909, 0.8791, 0.8737, 0.8747, 0.8747, 0.8134 /

      data sig2 / 0.7386, 0.7412, 0.7460, 0.7496, 0.7558, 0.7681, 0.7792, 0.7976, 0.8110,
     1           0.8286, 0.8370, 0.8461, 0.8508, 0.8555, 0.8613, 0.8747, 0.8880,
     2           0.8902, 0.8799, 0.8679, 0.8625, 0.8635, 0.8635, 0.8090 /

      data sig3 / 0.7160, 0.7188, 0.7238, 0.7276, 0.7343, 0.7471, 0.7568, 0.7717, 0.7814,
     1           0.7918, 0.7956, 0.7983, 0.7994, 0.8004, 0.8023, 0.8137, 0.8267,
     2           0.8275, 0.8158, 0.8022, 0.7958, 0.7963, 0.7963, 0.7661 /

      data sig4 / 0.6855, 0.6883, 0.6932, 0.6967, 0.7025, 0.7139, 0.7203, 0.7314, 0.7361,
     1           0.7360, 0.7340, 0.7295, 0.7270, 0.7246, 0.7233, 0.7328, 0.7457,
     2           0.7452, 0.7315, 0.7157, 0.7078, 0.7076, 0.7076, 0.7163 /

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(23)
        sig1F = sig1(23)
        sig2F = sig2(23)
        sig3F = sig3(23)
        sig4F = sig4(23)
        goto 1021
      elseif (specT .eq. -1.0) then
        specF = freq(24)
        sig1F = sig1(24)
        sig2F = sig2(24)
        sig3F = sig3(24)
        sig4F = sig4(24)
        goto 1021
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        c1 = 0
        c2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            c1 = ifreq
            c2 = ifreq+1
            goto 1020
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     interpolate between frequencies
1020  call S24_interp1 ( freq(c1), freq(c2), sig1(c1), sig1(c2), specF, sig1F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig2(c1), sig2(c2), specF, sig2F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig3(c1), sig3(c2), specF, sig3F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig4(c1), sig4(c2), specF, sig4F, iflag )

1021  if (m .le. 4.5) then
        sigma = sig1F
      elseif (m .gt. 4.5 .and. m .le. 5.0) then
        sigma = sig1F + (sig2F - sig1F) * ((m - 4.5)/0.5)
      elseif (m .gt. 5.0 .and. m .le. 5.5) then
        sigma = sig2F + (sig3F - sig2F) * ((m - 5.0)/0.5)
      elseif (m .gt. 5.5 .and. m .le. 6.5) then
        sigma = sig3F + (sig4F - sig3F) * ((m - 5.5)/1.0)
      else
        sigma = sig4F
      endif

      return
      end

c ---------------------------------------------------------------------
      subroutine S32_NGAEast_CompSSSig_Low ( m, specT, sigma, iflag )

c     model: NGA-East composite single-station sigma, CENA, low
c     ref: Goulet et al., 2018 (PEER Report 2018/08), equation 11-14 for functional form,
c          Table 11-18 for coefficients (also provided in Appendix H)

      implicit none

      integer MAX_FREQ
      parameter (MAX_FREQ=25)
      integer iflag, ifreq, c1, c2
      real m, specT, specF, freq(MAX_FREQ), sig1(MAX_FREQ), sig2(MAX_FREQ),
     1     sig3(MAX_FREQ), sig4(MAX_FREQ), sig1F, sig2F, sig3F, sig4F, sigma

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1           2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2           20.000, 25.000, 33.333, 40.000, 50.000, 100.000, 100.000, -1.0 /

      data sig1 / 0.4695, 0.4702, 0.4740, 0.4784, 0.4870, 0.5037, 0.5171, 0.5354, 0.5467,
     1            0.5595, 0.5650, 0.5707, 0.5736, 0.5765, 0.5795, 0.5825, 0.5840,
     2            0.5855, 0.5860, 0.5867, 0.5870, 0.5873, 0.5879, 0.5879, 0.5220 /

      data sig2 / 0.4444, 0.4451, 0.4491, 0.4538, 0.4628, 0.4803, 0.4944, 0.5135, 0.5253,
     1            0.5387, 0.5444, 0.5503, 0.5533, 0.5564, 0.5595, 0.5626, 0.5641,
     2            0.5657, 0.5662, 0.5669, 0.5672, 0.5676, 0.5682, 0.5682, 0.5159 /

      data sig3 / 0.4053, 0.4061, 0.4101, 0.4148, 0.4243, 0.4433, 0.4566, 0.4741, 0.4837,
     1            0.4913, 0.4931, 0.4928, 0.4924, 0.4921, 0.4917, 0.4933, 0.4946,
     2            0.4960, 0.4964, 0.4970, 0.4973, 0.4976, 0.4981, 0.4981, 0.4703 /

      data sig4 / 0.4017, 0.4025, 0.4071, 0.4124, 0.4226, 0.4416, 0.4482, 0.4508, 0.4504,
     1            0.4411, 0.4326, 0.4168, 0.4079, 0.3994, 0.3904, 0.3882, 0.3887,
     2            0.3893, 0.3895, 0.3897, 0.3898, 0.3899, 0.3901, 0.3901, 0.3928 /

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(24)
        sig1F = sig1(24)
        sig2F = sig2(24)
        sig3F = sig3(24)
        sig4F = sig4(24)
        goto 1021
      elseif (specT .eq. -1.0) then
        specF = freq(25)
        sig1F = sig1(25)
        sig2F = sig2(25)
        sig3F = sig3(25)
        sig4F = sig4(25)
        goto 1021
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        c1 = 0
        c2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            c1 = ifreq
            c2 = ifreq+1
            goto 1020
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     interpolate between frequencies
 1020 call S24_interp1 ( freq(c1), freq(c2), sig1(c1), sig1(c2), specF, sig1F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig2(c1), sig2(c2), specF, sig2F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig3(c1), sig3(c2), specF, sig3F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig4(c1), sig4(c2), specF, sig4F, iflag )

1021  if (m .le. 4.5) then
        sigma = sig1F
      elseif (m .gt. 4.5 .and. m .le. 5.0) then
        sigma = sig1F + (sig2F - sig1F) * ((m - 4.5)/0.5)
      elseif (m .gt. 5.0 .and. m .le. 5.5) then
        sigma = sig2F + (sig3F - sig2F) * ((m - 5.0)/0.5)
      elseif (m .gt. 5.5 .and. m .le. 6.5) then
        sigma = sig3F + (sig4F - sig3F) * ((m - 5.5)/1.0)
      else
        sigma = sig4F
      endif

      return
      end

c ---------------------------------------------------------------------
      subroutine S32_NGAEast_CompSSSig_Cen ( m, specT, sigma, iflag )

c     model: NGA-East composite single-station sigma, CENA, central
c     ref: Goulet et al., 2018 (PEER Report 2018/08), equation 11-14 for functional form,
c          Table 11-18 for coefficients (also provided in Appendix H)

      implicit none

      integer MAX_FREQ
      parameter (MAX_FREQ=25)
      integer iflag, ifreq, c1, c2
      real m, specT, specF, freq(MAX_FREQ), sig1(MAX_FREQ), sig2(MAX_FREQ),
     1     sig3(MAX_FREQ), sig4(MAX_FREQ), sig1F, sig2F, sig3F, sig4F, sigma

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1           2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2           20.000, 25.000, 33.333, 40.000, 50.000, 100.000, 100.000, -1.0 /

      data sig1 / 0.5853, 0.5858, 0.5889, 0.5925, 0.5998, 0.6143, 0.6262, 0.6432, 0.6541,
     1            0.6671, 0.6730, 0.6794, 0.6827, 0.6862, 0.6898, 0.6935, 0.6955,
     2            0.6974, 0.6981, 0.6990, 0.6994, 0.6998, 0.7006, 0.7006, 0.6303 /

      data sig2 / 0.5661, 0.5666, 0.5698, 0.5736, 0.5811, 0.5959, 0.6082, 0.6256, 0.6369,
     1            0.6502, 0.6563, 0.6628, 0.6662, 0.6698, 0.6735, 0.6773, 0.6793,
     2            0.6813, 0.6820, 0.6829, 0.6833, 0.6838, 0.6846, 0.6846, 0.6247 /

      data sig3 / 0.5357, 0.5364, 0.5397, 0.5438, 0.5517, 0.5674, 0.5781, 0.5921, 0.5999,
     1            0.6060, 0.6075, 0.6072, 0.6068, 0.6066, 0.6063, 0.6079, 0.6092,
     2            0.6105, 0.6110, 0.6116, 0.6119, 0.6122, 0.6127, 0.6127, 0.5765 /

      data sig4 / 0.5095, 0.5102, 0.5138, 0.5181, 0.5267, 0.5436, 0.5510, 0.5576, 0.5577,
     1           0.5482, 0.5399, 0.5254, 0.5171, 0.5092, 0.5007, 0.4978, 0.4977,
     2           0.4975, 0.4975, 0.4974, 0.4974, 0.4974, 0.4973, 0.4973, 0.4896 /

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(24)
        sig1F = sig1(24)
        sig2F = sig2(24)
        sig3F = sig3(24)
        sig4F = sig4(24)
        goto 1021
      elseif (specT .eq. -1.0) then
        specF = freq(25)
        sig1F = sig1(25)
        sig2F = sig2(25)
        sig3F = sig3(25)
        sig4F = sig4(25)
        goto 1021
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        c1 = 0
        c2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            c1 = ifreq
            c2 = ifreq+1
            goto 1020
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     interpolate between frequencies
1020  call S24_interp1 ( freq(c1), freq(c2), sig1(c1), sig1(c2), specF, sig1F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig2(c1), sig2(c2), specF, sig2F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig3(c1), sig3(c2), specF, sig3F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig4(c1), sig4(c2), specF, sig4F, iflag )

1021  if (m .le. 4.5) then
        sigma = sig1F
      elseif (m .gt. 4.5 .and. m .le. 5.0) then
        sigma = sig1F + (sig2F - sig1F) * ((m - 4.5)/0.5)
      elseif (m .gt. 5.0 .and. m .le. 5.5) then
        sigma = sig2F + (sig3F - sig2F) * ((m - 5.0)/0.5)
      elseif (m .gt. 5.5 .and. m .le. 6.5) then
        sigma = sig3F + (sig4F - sig3F) * ((m - 5.5)/1.0)
      else
        sigma = sig4F
      endif

      return
      end

c ---------------------------------------------------------------------
      subroutine S32_NGAEast_CompSSSig_High ( m, specT, sigma, iflag )

c     model: NGA-East composite single-station sigma, CENA, high
c     ref: Goulet et al., 2018 (PEER Report 2018/08), equation 11-14 for functional form,
c          Table 11-18 for coefficients (also provided in Appendix H)

      implicit none

      integer MAX_FREQ
      parameter (MAX_FREQ=25)
      integer iflag, ifreq, c1, c2
      real m, specT, specF, freq(MAX_FREQ), sig1(MAX_FREQ), sig2(MAX_FREQ),
     1     sig3(MAX_FREQ), sig4(MAX_FREQ), sig1F, sig2F, sig3F, sig4F, sigma

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1           2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2           20.000, 25.000, 33.333, 40.000, 50.000, 100.000, 100.000, -1.0 /

      data sig1 / 0.7357, 0.7358, 0.7364, 0.7372, 0.7391, 0.7444, 0.7502, 0.7611, 0.7697,
     1            0.7816, 0.7876, 0.7943, 0.7980, 0.8019, 0.8061, 0.8106, 0.8129,
     2            0.8153, 0.8162, 0.8173, 0.8178, 0.8183, 0.8193, 0.8193, 0.7480 /

      data sig2 / 0.7221, 0.7222, 0.7228, 0.7237, 0.7258, 0.7313, 0.7373, 0.7484, 0.7572,
     1            0.7693, 0.7753, 0.7822, 0.7859, 0.7899, 0.7941, 0.7986, 0.8010,
     2            0.8034, 0.8043, 0.8054, 0.8059, 0.8064, 0.8075, 0.8075, 0.7431 /

      data sig3 / 0.6989, 0.6990, 0.6999, 0.7009, 0.7033, 0.7093, 0.7137, 0.7209, 0.7256,
     1            0.7295, 0.7305, 0.7300, 0.7297, 0.7295, 0.7292, 0.7306, 0.7317,
     2            0.7329, 0.7333, 0.7338, 0.7340, 0.7343, 0.7348, 0.7348, 0.6960 /

      data sig4 / 0.6679, 0.6679, 0.6682, 0.6687, 0.6698, 0.6738, 0.6743, 0.6772, 0.6762,
     1            0.6681, 0.6623, 0.6537, 0.6492, 0.6453, 0.6413, 0.6395, 0.6390,
     2            0.6385, 0.6384, 0.6382, 0.6381, 0.6379, 0.6377, 0.6377, 0.6419 /

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(24)
        sig1F = sig1(24)
        sig2F = sig2(24)
        sig3F = sig3(24)
        sig4F = sig4(24)
        goto 1021
      elseif (specT .eq. -1.0) then
        specF = freq(25)
        sig1F = sig1(25)
        sig2F = sig2(25)
        sig3F = sig3(25)
        sig4F = sig4(25)
        goto 1021
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        c1 = 0
        c2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            c1 = ifreq
            c2 = ifreq+1
            goto 1020
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     interpolate between frequencies
1020  call S24_interp1 ( freq(c1), freq(c2), sig1(c1), sig1(c2), specF, sig1F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig2(c1), sig2(c2), specF, sig2F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig3(c1), sig3(c2), specF, sig3F, iflag )
      call S24_interp1 ( freq(c1), freq(c2), sig4(c1), sig4(c2), specF, sig4F, iflag )

1021  if (m .le. 4.5) then
        sigma = sig1F
      elseif (m .gt. 4.5 .and. m .le. 5.0) then
        sigma = sig1F + (sig2F - sig1F) * ((m - 4.5)/0.5)
      elseif (m .gt. 5.0 .and. m .le. 5.5) then
        sigma = sig2F + (sig3F - sig2F) * ((m - 5.0)/0.5)
      elseif (m .gt. 5.5 .and. m .le. 6.5) then
        sigma = sig3F + (sig4F - sig3F) * ((m - 5.5)/1.0)
      else
        sigma = sig4F
      endif

      return
      end
