
c ------------------------------------------------------------------            
C *** BCHydro Subduction (03/06/2009 - Model) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S05_BCHydroSub_V2 ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth )

      implicit none

      integer iflag, forearc     
      real mag, fType, rRup, vs30, pgaRock, lnSa, sigma, tau, period1, 
     1     sigma1, specT, depth, period0, vs30_rock, faba

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface
C      1       Intraslab
C
C     faba     Note
C     -------------------------
C      0       Non-Backarc site  
C      1       Backarc site  
C
c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1100.
      faba = real(forearc)

C     Compute regular ground motions. 
      call S05_BCHydroSub2008_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
c ----------------------------------------------------------------------

      subroutine S05_BCHydroSub2008_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth )
      
      implicit none
      
      integer MAXPER, count1, count2, iflag, i1, i, nPer
      parameter (MAXPER=22)
      real a1(MAXPER), a2(MAXPER), a3(MAXPER), a5(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, 
     1     mag, a4, a9, specT, tau 
      real a1T, a2T, a3T, a5T, a6T, a7T, a8T
      real a10T, a11T, a12T, a13T, a14T, sigsT, sigtT
      real vLinT, b_soilT, sum, Ftype, rhypo, faba
      real n, c, period1, R1, R2, a18, depth

      data vLin / 1000., 1000., 1000., 1000., 1000., 1000., 1000., 
     1            1000., 1000., 1000., 1000., 1000., 1000., 1000., 
     2            1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000./
      data b_soil /  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 0.000, 
     1               0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 0.000,
     2               0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 0.000 /
      data period /  0.00, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.6,  0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.00  /

      data a1 / 2.6644, 2.8484, 3.1579, 3.3794, 3.5007, 3.452, 3.2724, 3.1639,
     1          2.9919, 2.6301, 2.3927, 2.1078, 2.2869, 1.497, 1.476,  0.8737,
     2          0.9018, 0.3009, 0.0672, -0.1194, -0.3091, -1.3966 /
      data a2 / -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
     1          -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 /
      data a3 / 1.1521, 1.2787, 1.5491, 1.4819, 1.3539, 1.1972, 1.0734, 0.9436,
     1          0.7281, 0.5553, 0.4292, 0.2757, 0.0205,-0.2536,-0.4601,-0.5586,
     2         -0.6034,-0.6309,-0.5559,-0.4668,-0.3206,-0.2402 /

      data a5 / -0.0071, -0.0076, -0.0079, -0.0082, -0.0082, -0.008,  -0.0078,
     1          -0.0073, -0.0065, -0.0062, -0.0062, -0.0059, -0.0049, -0.0054,
     2          -0.0041, -0.0041, -0.0035, -0.004,  -0.0034, -0.0032, -0.0036, 0.0043 /
      data a6 / 1.6525, 1.993,  2.2121, 2.2901, 2.4929, 2.3676, 2.2501, 2.1631, 
     1          1.9101, 1.7412, 1.5027, 1.2842, 0.9797, 0.3446,-0.105, -0.5123,
     2         -0.7693,-1.1911,-1.6938,-2.1727,-2.7298,-3.5885 /
      data a7 / -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
     1          -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 /
      data a8 / -1.2978, -1.4614, -1.6816, -1.6447, -1.4951, -1.3442, -1.2068,
     1          -1.0841, -0.8765, -0.708,  -0.5685, -0.3989, -0.1883,  0.077, 
     2           0.2301,  0.3218, 0.3756,   0.4152,  0.3542,  0.24,    0.124, 0.0 /
      data a10 / -0.0071, -0.0071, -0.0066, -0.0066, -0.0071, -0.0073, -0.0069,
     1           -0.007,  -0.007,  -0.007,  -0.0069, -0.0072, -0.007,  -0.0061,
     2           -0.0053, -0.0045, -0.0042, -0.0041, -0.0032, -0.0017, -0.0007, 0.0014 /
      data a11 / 0.0092, 0.009,  0.0092, 0.0103, 0.0098, 0.0102, 0.0101, 0.01, 
     1           0.0085, 0.0104, 0.011,  0.0106, 0.0085, 0.0098, 0.0064, 0.0047,
     2           0.0023, 0.0007,-0.0004,-0.0035,-0.0048,-0.0086 /
      data a12 / -0.3887, -0.2674, -0.2095, -0.2302, -0.3215, -0.4166, -0.4824,
     1           -0.5463, -0.6388, -0.6778, -0.7366, -0.794,  -0.8085, -0.791, 
     2           -0.7707, -0.7469, -0.7054, -0.6446, -0.6081, -0.5817, -0.5324, -0.4974 /
      data a13 / -0.0076, -0.0043, -0.0114, -0.0074, -0.0156, -0.0078, -0.0034,
     1           -0.003,  -0.0071,  0.0124,  0.0222,  0.025,  -0.0322,  0.0031,
     2           -0.052,  -0.0143, -0.0518, -0.0263, -0.0472, -0.0688, -0.0828, -0.0908 /
      data a14 / 0.0074, 0.0147, 0.0234, 0.0364, 0.0153, 0.0196, 0.0043, -0.0073,
     1          -0.014, -0.0393,-0.0484,-0.0663,-0.0893,-0.1271,-0.1515, -0.1682, 
     2          -0.1828,-0.2032,-0.2066,-0.2036,-0.1925,-0.1537 /

      data sigs / 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61,
     1            0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61,
     2            0.61, 0.61 /
      data sigt /0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42,
     1           0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.42 /

C Constant parameters            
      n = 1.18
      c = 1.88
      a4 = 1.73
      a9 = 1.45
      a18 = 0.77
 
C Find the requested spectral period and corresponding coefficients
      nPer = 22

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a3T = a3(i1)
         a5T = a5(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1+3
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (3/2009 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +                   specT,a5T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)

 1011 period1 = specT                                                                                                              

c     Compute distances R1 or R2 depending on Event Type.
      R1 = rRup + 15.0*exp(0.7*(mag-6.0))

C Currently the rRup is set equal to Rhypo for Intraslab events. 
      rhypo = rRup
      R2 = rhypo + 15.0*exp(0.0*(mag-6.0))

C Compute the Base Model
C  Interface events
      if (fType .eq. 0) then
          if (mag .lt. 7.25) then
             sum = a1T + a4*(mag-7.25) + a13T*(10.0-mag)**2.0 + a2T * alog(R1) + 
     1                  a5T*rRup + a12T*alog(vs30/1000)
          else 
             sum = a1T + a18*(mag-7.25) + a13T*(10.0-mag)**2.0 + a2T * alog(R1) + 
     1                  a5T*rRup + a12T*alog(vs30/1000)
          endif

C  Intraslab events
      elseif (fType .eq. 1) then
          sum = a6T + a9*(mag-6.0) + a14T*(8.0-mag)**2.0 + a7T * alog(R2) + a10T*rhypo +
     1             a11T*(depth - 60.0) + faba*(a3T + a8T*alog(R2/40.0)) + a12T*alog(vs30/1000) 

      else
         write (*,*) 'Undefined Mechanism for BC Hydro Subduction Attenuation Model.'
         write (*,*) 'Ftype = ', ftype
         write (*,*) 'Check input fault file.'
         stop 99
      endif

C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = sum

      return
      end

c ----------------------------------------------------------------------

      subroutine S05_BCHydroSub2008_modelsm ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth )
      
      implicit none
      
      integer MAXPER, nPer, i1, count1, count2, iflag, i
      parameter (MAXPER=22)
      real a1(MAXPER), a2(MAXPER), a3(MAXPER), a4(MAXPER), a5(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, depth, faba, mag
      real a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real vLinT, b_soilT, sum, Ftype, rhypo
      real n, c, period1, c4, specT, tau

      data vLin / 1000., 1000., 1000., 1000., 1000., 1000., 1000., 
     1            1000., 1000., 1000., 1000., 1000., 1000., 1000., 
     2            1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000./
      data b_soil /  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 0.000, 
     1               0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 0.000,
     2               0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000, 0.000 /
      data period /  0.00, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.6,  0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.00  /

      data a1 / 2.549, 2.771, 2.96, 3.148, 3.222, 3.178, 3.074, 2.996, 2.779, 2.623, 2.449,
     1          2.163, 1.783, 1.302, 0.876, 0.719, 0.399, 0.095, -0.246, -0.583, -0.917, -1.291 /
      data a2 / 2.213, 2.58, 2.88, 3.013, 3.072, 2.962, 2.822, 2.747, 2.583, 2.375, 2.177, 
     1          1.907, 1.622, 0.882, 0.402, -0.006, -0.287, -0.819, -1.248, -1.621, -2.134, -2.669 /
      data a3 / -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, 
     1          -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54, -1.54 /
      data a4 / -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66,
     1          -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66, -0.66 /
      data a5 / -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, 	
     1          -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24, -2.24 /
      data a6 / -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, 
     1          -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4 /
      data a7 / 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56,
     1          0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56 /
      data a8 / -0.004, -0.004, -0.004, -0.004, -0.0037, -0.0032, -0.0028, -0.0024, -0.0019, 
     1         -0.0014, -0.0011, -0.0007, -0.0001, 0.0007, 0.0012, 0.0019, 0.0024, 0.0033, 
     2          0.0039, 0.0045, 0.0051, 0.006 /
      data a9 / -0.003, -0.003, -0.003, -0.003, -0.003, -0.003, -0.003, -0.003, -0.003, 
     1          -0.003, -0.003, -0.003, -0.003, -0.0016, -0.0006, 0.0002, 0.0008, 0.0018, 
     2          0.0026, 0.0032, 0.004, 0.005 /
      data a10 / 1.03, 0.971, 0.993, 0.972, 1.023, 1.033, 0.988, 0.922, 0.851, 0.797, 
     1          0.756, 0.707, 0.583, 0.459, 0.316, 0.221, 0.151, 0.065, -0.01, -0.093, -0.143, -0.238 /
      data a11 / -1.2, -1.2, -1.2, -1.2, -1.2, -1.2, -1.132, -1.076, -0.987, -0.919, -0.863, 
     1         -0.795, -0.706, -0.582, -0.494, -0.425, -0.369, -0.281, -0.213, -0.157, -0.088, 0.0 /
      data a12 / 0.0108, 0.0108, 0.0114, 0.0122, 0.0112, 0.0113, 0.0116, 0.0113, 0.0093, 0.0114, 
     1           0.0122, 0.0118, 0.0097, 0.0112, 0.0079, 0.0061, 0.0037, 0.0025, 0.0015, -0.0009, 
     2          -0.0013, -0.0031 /
      data a13 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.006,
     1       -0.0179, -0.0272, -0.0348, -0.0468, -0.0561, -0.0637, -0.073, -0.085 /
      data a14 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0145, -0.0373, -0.055, -0.0695, 
     1        -0.0872, -0.11, -0.1402, -0.1617, -0.1783, -0.1919, -0.2134, -0.23, -0.23,
     1          -0.23, -0.23 /
      data a15 / -0.3727, -0.2656, -0.2236, -0.2405, -0.3219, -0.4054, -0.4618, -0.5079, 
     1           -0.5921, -0.6325, -0.6964, -0.7709, -0.787, -0.7901, -0.7687, -0.7518, 
     2           -0.6965, -0.644, -0.606, -0.5716, -0.5295, -0.4905 /
      data a16 / -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, 
     1           -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4 /

      data sigs / 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
     1            0.6, 0.6, 0.6, 0.576, 0.54, 0.514, 0.493, 0.469, 0.44 /
      data sigt / 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 
     1            0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41 /

C Constant parameters            
      n = 1.18
      c = 1.88
 
C Find the requested spectral period and corresponding coefficients
      nPer = 22

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a3T = a3(i1)
         a4T = a4(i1)
         a5T = a5(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a9T = a9(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1+3
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (3/2009 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a3(count1),a3(count2),
     +                   specT,a3T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),a5(count1),a5(count2),
     +                   specT,a5T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a9(count1),a9(count2),
     +                   specT,a9T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)

 1011 period1 = specT                                                                                                              


C     Term needs to be set.
      c4 = 0.0

C Compute the Base Model
C  Interface events
      if (fType .eq. 0) then
          if (mag .lt. 6.5) then
             sum = a1T + a3T*(mag-6.5) + a13T*(10.0-mag)**2.0 +         
     1             (a6T + a7T*(mag-6.5)) * alog(Rrup + c4) + 
     1                  a8T*Rrup + a15T*alog(vs30/1000.0)
          elseif (mag .le. 7.25) then
             sum = a1T + a4T*(mag-6.5) + a13T*(10.0-mag)**2.0 +         
     1             (a6T + a7T*(mag-6.5)) * alog(Rrup + c4) + 
     1                  a8T*Rrup + a15T*alog(vs30/1000.0)
          else 
             sum = a1T + 0.75*a4T + a5T*(mag-7.25) + a13T*(10.0-mag)**2.0 +         
     1             (a6T + a7T*(mag-6.5)) * alog(Rrup + c4) + 
     1                  a8T*Rrup + a15T*alog(vs30/1000.0)
          endif

C  Intraslab events
      elseif (fType .eq. 1) then
          if (mag .lt. 6.5) then
             sum = a2T + a3T*(mag-6.5) + a14T*(8.0-mag)**2.0 +         
     1             a7T*(mag-6.0)*alog(100.0 + c4) + a16T*alog(rHypo + c4) + 
     1             a9T*rHypo + a12T*(depth - 60.0) + 
     1             faba*(a10T + a11T*alog((rHypo + c4)/40.0) ) + 
     1             a15T*alog(vs30/1000.0)
          elseif (mag .le. 7.25) then
             sum = a2T + a4T*(mag-6.5) + a14T*(8.0-mag)**2.0 +         
     1             a7T*(mag-6.0)*alog(100.0 + c4) + a16T*alog(rHypo + c4) + 
     1             a9T*rHypo + a12T*(depth - 60.0) + 
     1             faba*(a10T + a11T*alog((rHypo + c4)/40.0) ) + 
     1             a15T*alog(vs30/1000.0)
          else 
             sum = a2T + 0.75*a4T + a5T*(mag-7.25) + a14T*(8.0-mag)**2.0 +         
     1             a7T*(mag-6.0)*alog(100.0 + c4) + a16T*alog(rHypo + c4) + 
     1             a9T*rHypo + a12T*(depth - 60.0) + 
     1             faba*(a10T + a11T*alog((rHypo + c4)/40.0) ) + 
     1             a15T*alog(vs30/1000.0)
          endif
      else
         write (*,*) 'Undefined Mechanism for BC Hydro Subduction Attenuation Model.'
         write (*,*) 'Ftype = ', ftype
         write (*,*) 'Check input fault file.'
         stop 99
      endif

C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = sum

      return
      end

c ------------------------------------------------------------------            
C *** BCHydro Subduction (06/2010 - Model) Horizontal ***********
c ------------------------------------------------------------------            
      subroutine S05_BCHydroSub_V3 ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)

C     Compute Rock PGA
      call S05_BCHydroSub2010_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo, deltac1 )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S05_BCHydroSub2010_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
c ----------------------------------------------------------------------
      subroutine S05_BCHydroSub2010_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1 )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=23)
      real a1(MAXPER), a2(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a4, a5, a9 
      real a1T, a2T, a6T, a7T, a8T
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, deltac1, faba, R, testmag, VsStar, depth, specT

      data period /  0.00, 0.02, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.6,  0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.00  /
      data vLin / 865.1, 865.1, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1, 
     1            503.0,  456.6,  430.3,  410.5, 400.0, 400.0, 400.0, 400.0, 
     2            400.0,  400.0,  400.0,  400.0, 400.0, 400.0 / 
      data b_soil / -1.186,-1.186, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, -2.518,
     1               -2.657, -2.669, -2.599, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 0, 0 /
      data a1 / 4.2203,4.2203, 4.5371, 5.0733, 5.2892, 5.4563, 5.2684, 5.0594, 4.7945, 
     1          4.4644, 4.0181, 3.6055, 3.2174, 2.7981, 2.0123, 1.4128, 0.9976,
     2          0.6443, 0.0657, -0.4624, -0.9809, -1.6017, -2.2937 /
      data a2 / -1.35,-1.35, -1.4, -1.45, -1.45, -1.45, -1.4, -1.35, -1.28, -1.18, -1.08, -0.99,
     1          -0.91, -0.85, -0.77, -0.71, -0.67, -0.64, -0.58, -0.54, -0.5, -0.46, -0.4 /
      data a6 / -0.0012,-0.0012, -0.0012, -0.0012, -0.0012, -0.0014, -0.0018, -0.0023, -0.0027,
     1          -0.0035, -0.0044, -0.0050, -0.0058, -0.0062, -0.0064, -0.0064, -0.0064, 
     2          -0.0064, -0.0064, -0.0064, -0.0064, -0.0064, -0.0064 /
      data a7 / 1.0988,1.0988, 1.2536, 1.4175, 1.3997, 1.3582, 1.1648, 0.994, 0.8821,
     1          0.7046, 0.5799, 0.5021, 0.3687, 0.1746, -0.082, -0.2821, -0.4108,
     2         -0.4466, -0.4344, -0.4368, -0.4586, -0.4433, -0.4828 /
      data a8 / -1.42,-1.42, -1.65, -1.80, -1.80, -1.69, -1.49, -1.30, -1.18, -0.98, -0.82, 
     1          -0.70, -0.54, -0.34, -0.05, 0.12, 0.25, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 /
      data a10 / 3.12,3.12, 3.37, 3.37, 3.33, 3.25, 3.03, 2.8, 2.59, 2.2, 1.92, 
     1            1.7, 1.42, 1.1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 /
      data a11 / 0.0130,0.0130, 0.0130, 0.0130, 0.0130, 0.0130, 0.0129, 0.0129, 
     1           0.0128, 0.0127, 0.0125, 0.0124, 0.0120, 0.0114, 0.0100,
     2           0.0085, 0.0069, 0.0054, 0.0027, 0.0005, -0.0013, -0.0033, -0.0060 /
      data a12 / 0.980,0.980, 1.288, 1.483, 1.613, 1.882, 2.076, 2.248, 2.348, 2.427, 
     1           2.399, 2.273, 1.993, 1.470, 0.408, -0.401, -0.723, -0.673, -0.627,
     2          -0.596, -0.566, -0.528, -0.504 /
      data a13 / -0.0135,-0.0135, -0.0138, -0.0142, -0.0145, -0.0153, -0.0162, -0.0172, -0.0183,
     1           -0.0206, -0.0231, -0.0256, -0.0296, -0.0363, -0.0493, -0.0610, -0.0711, 
     2           -0.0798, -0.0935, -0.0980, -0.0980, -0.0980, -0.0980 /
      data a14 / -0.40,-0.40, -0.40, -0.40, -0.40, -0.40, -0.35, -0.31, -0.28, -0.23, -0.19, 
     1           -0.16, -0.12, -0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /
      data a15 / 0.9969,0.9969, 1.1030, 1.2732, 1.3042, 1.2600, 1.2230, 1.1600, 1.0500,
     1           0.8000, 0.6620, 0.5800, 0.4800, 0.3300, 0.3100, 0.3000, 0.3000, 
     2           0.3000, 0.3000, 0.3000, 0.3000, 0.3000, 0.3000 /
      data a16 / -1.00,-1.00, -1.18, -1.36, -1.36, -1.30, -1.25, -1.17, -1.06, -0.78, -0.62,
     1           -0.50, -0.34, -0.14, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /
C     Changed sigs=0.60 and tau=0.43 (January 17, 2014). Also changed SSS=0.60
C      data sigs /  0.603,0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 
c     1            0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603, 0.603 /
c      data sigt /  0.482,0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 
c     1            0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482, 0.482 /
      data sigs /  0.60,0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 
     1            0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60 /
      data sigt /  0.43,0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 
     1            0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43 /

C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a4 = 0.9
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8
 
C Find the requested spectral period and corresponding coefficients
      nPer = 23

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (6/2010 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)

 1011 period1 = specT                                                                                                              

C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if (ftype .eq. 0) then
         R = rRup + c4*exp( (mag-6.0)*a9 ) 
         sumgm = a1T + a4*deltaC1 + (a2T + a14T*ftype + a3*(mag - 7.8))*alog(R) + a6T*rRup + a10T*ftype
      elseif (ftype .eq. 1) then
         R = disthypo + c4*exp( (mag-6.0)*a9 ) 
         sumgm = a1T + a4*deltaC1 + (a2T + a14T*ftype + a3*(mag - 7.8))*alog(R) + a6T*disthypo + a10T*ftype
      else
         write (*,*) 'BC Hydro V3 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
      
C     Base model for Magnitude scaling.      
      testmag = (7.8 + deltaC1)
      if (mag .le. testmag ) then
         sumgm = sumgm + a4*(mag-testmag) + a13T*(10.0-mag)**2.0
      else
         sumgm = sumgm + a5*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif      
      
C     Depth Scaling
      if ( depth .gt. 120. ) then
        sumgm = sumgm + a11T*(120. - 60.0)*ftype
      else
        sumgm = sumgm + a11T*(depth - 60.0)*ftype
      endif

C     Forearc/Backarc scaling      
      if (ftype .eq. 1) then
         sumgm = sumgm + (a7T +a8T*alog(max(disthypo,85.0)/40.0))*faba
      elseif (ftype .eq. 0) then   
         sumgm = sumgm + (a15T +a16T*alog(max(rRup,100.0)/40.0))*faba
      endif 

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         sumgm = sumgm + a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         sumgm = sumgm + a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif

C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = sumgm

      return
      end

c ----------------------------------------------------------------------
      subroutine S05_BCHHR2Vs760 ( saRock, specT, lnSa )

      implicit none

      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=16)
      real b1(MAXPER), b2(MAXPER), b3(MAXPER)
      real c1(MAXPER), c2(MAXPER), period(MAXPER)
      real saRock, lnSa, specT, period1 
      real b1T, b2T, b3T, c1T, c2T, LsaRock, ampfac
      integer count1, count2, iflag

      data period /  0.00, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5,
     1               1.00, 2.00, 3.3,  5.0, 10.00  /
      data b1 / 0.044, 0.044, -0.219, -0.434, -0.424, -0.267, -0.102, 0.237, 0.29, 
     1          0.333, 0.379, 0.481, 0.418, 0.37, 0.307, 0.241 /
      data b2 / -0.16, -0.16, -0.354, -0.273, -0.126, -0.031, 0.015, -0.006, 0.0,
     1          -0.006, 0.019, 0.018, -0.002, 0.0, 0.0, 0.0 / 
      data b3 / 0.0072, 0.0072, 0.031, 0.0179, -0.0042, -0.0238, -0.0294, 
     1     -0.0333, -0.0384, -0.0289, 0.00 ,0.0174, 0.011, 0.0054, 0.0028, 0.0042 /
      data c1 / -3.0, -3.0, -2.5, -2.5, -2.5, -2.0, -2.0, -1.5, -1.0, -0.5, -2.0, -2.0, -2.5, -3.0, -4.0, -5.0 /
      data c2 / 1.0, 1.0, 1.8, 2.0, 1.9, 1.7, 1.6, 1.4, 1.3, 1.2, 0.9, 0.4, -0.2, -0.8, -1.4, -2.9 /
 
C Find the requested spectral period and corresponding coefficients
      nPer = 16

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         b1T = b1(i1)
         b2T = b2(i1)
         b3T = b3(i1)
         c1T = c1(i1)
         c2T = c2(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Hard Rock to Vs760 Amp Factors'
      write (*,*) 'model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),b1(count1),b1(count2),
     +                   specT,b1T,iflag)
            call S24_interp (period(count1),period(count2),b2(count1),b2(count2),
     +                   specT,b2T,iflag)
            call S24_interp (period(count1),period(count2),b3(count1),b3(count2),
     +                   specT,b3T,iflag)
            call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
            call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +                   specT,c2T,iflag)
 1011 period1 = specT                                                                                                              

C     Now compute the amp factor based on the given hard rock spectral acceleration.
C     Convert Hard Rock SA from Gals to G
      LsaRock = saRock - 6.89
C     Check for LsaRock value between bounding values C1 and C2 
      if (LsaRock .le. c1T) then
         ampfac = b1T
      elseif (lsaRock .gt. c2T) then
         ampfac = b1T + b2T*(c2T-c1T) + b3T*(c2T-c1T)*(c2T-c1T)
      else
         ampfac = b1T + b2T*(LsaRock-c1T) + b3T*(LsaRock-c1T)*(LsaRock-c1T)
      endif     

C     Adjust ground motions and convert back to gals
      lnSa = LsaRock + ampfac + 6.89

      return
      end

c ------------------------------------------------------------------            
C *** BCHydro Subduction update (2018 - Model) Horizontal ***********
c     version17, Central Model, const from all reg
c ------------------------------------------------------------------            
      subroutine S05_BCHydroSub_v4a ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)

C     Compute Rock PGA
      call S05_BCHydroSub2018_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo, deltac1 )
      pgaRock = exp(lnSa)
c      write (*,'( 2x,''rockPGA ='',f10.4)') pgaRock
 
C     Compute regular ground motions. 
      call S05_BCHydroSub2018_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
c ----------------------------------------------------------------------
C *** BCHydro Subduction update (2018 - Model) Horizontal ***********
c     Version17, high scaled backbone model
c ------------------------------------------------------------------            
      subroutine S05_BCHydroSub_V4b ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, epistemic
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)

C     Compute Rock PGA
      call S05_BCHydroSub2018_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo, deltac1 )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S05_BCHydroSub2018_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT
      
c     add epistemic term
      if (ftype .eq. 1 ) then
       if (specT .le. 0.25) then
         epistemic = 0.5
       elseif (specT .gt. 0.25 .and. specT .le. 1.) then
         epistemic = 0.5 + (-0.2)*(alog(specT)-alog(0.25)) / (alog(1.)-alog(0.25))
       else
         epistemic = 0.3
       endif
      else
        epistemic = 0.3
      endif     
      lnSA = lnSA + epistemic 

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
      
c ------------------------------------------------------------------            
C *** BCHydro Subduction update (2018 - Model) Horizontal ***********
c     Version17, low scaled backbone model
c ------------------------------------------------------------------            
      subroutine S05_BCHydroSub_V4c ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, epistemic
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)

C     Compute Rock PGA
      call S05_BCHydroSub2018_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo, deltac1 )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S05_BCHydroSub2018_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     add epistemic term
      if (ftype .eq. 1 ) then
       if (specT .le. 0.25) then
         epistemic = -0.5
       elseif (specT .gt. 0.25 .and. specT .le. 1.) then
         epistemic = -0.5 + (0.2)*(alog(specT)-alog(0.25)) / (alog(1.)-alog(0.25))
       else
         epistemic = -0.3
       endif
      else
         epistemic = -0.3
      endif     
      lnSA = lnSA + epistemic 


c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
      
c ----------------------------------------
      subroutine S05_BCHydroSub_Global ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, slab1
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)

c     temp set for south america      
      slab1 = 7.25

C     Compute Rock PGA
      call S05_BCHydroSub2018_Glob_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo, deltac1, slab1 )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S05_BCHydroSub2018_Glob_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1, slab1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end

c ----------------------------------------
      subroutine S05_BCHydroSub_southAm ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, slab1
      integer iflag, forearc

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)
      slab1 = 7.25

C     Compute Rock PGA
      call S05_BCHydroSub2018_SA_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau,
     2                     period0, Ftype, iflag, faba, depth, disthypo, deltac1, slab1 )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S05_BCHydroSub2018_SA_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1, slab1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
      
c ----------------------------------------------------------------------

      subroutine S05_BCHydroSub2018_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1T )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=24)
      real a1(MAXPER), a2(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER)
      real a45(MAXPER), a4(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a5, a9 
      real a1T, a2T, a6T, a7T, a8T, a4T, a45T
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, faba, R, testmag, VsStar, depth, specT
      real c1p, adjust_inter(MAXPER), adjust_slab(MAXPER)
      real adjust_slabT, adjust_interT, gm1, gm2, gm3, gm4, gm5
      real slabDC1(MAXPER), slabDC1T
      real deltaC1(MAXPER), deltaC1T

      data period /  0.00, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.60, 0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.01 /
      data vLin / 865.1, 865.1, 865.1, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1, 
     1            503.0,  456.6,  430.3,  410.5, 400.0, 400.0, 400.0, 400.0, 
     2            400.0,  400.0,  400.0,  400.0, 400.0, 400.0 / 
      data b_soil / -1.186,-1.219, -1.273, -1.346, -1.471, -1.624, -1.931,
     2              -2.188, -2.381, -2.518, -2.657, -2.669, -2.599, -2.401,
     3              -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 0, 0 /
      data a1 / 2.340, 2.360, 2.384, 2.446, 2.751, 3.019, 3.349, 3.284, 
     1        3.211, 3.145, 2.997, 2.839, 2.658, 2.346, 1.851, 1.216, 
     2        0.649, 0.082, -0.369, -1.034, -1.520, -1.810, -2.173, -2.712 /

      data a2 / -1.044, -1.044, -1.08, -1.11, -1.11, -1.11, -1.084, 
     1          -1.027, -0.983, -0.947, -0.89, -0.845, -0.809, -0.76,
     2 -0.698, -0.612, -0.55, -0.501, -0.46, -0.455, -0.45, -0.45, -0.45, -0.45 /

      data a4/0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.62, 0.64, 0.66,
     1        0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 
     2        0.73, 0.78, 0.84, 0.93 /
      data a6/ -0.00705, -0.00707, -0.00710, -0.00725, -0.00758, -0.00788, 
     1         -0.0082, -0.00835, -0.00835, -0.00828, -0.00797, -0.0077, 
     2         -0.0074, -0.00698, -0.00645, -0.0057, -0.0051, -0.00465, 
     3         -0.0043, -0.0039, -0.0037, -0.00357, -0.0034, -0.00327 /
      data a7 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a8 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      data a10/ 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     1          1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     2          1.73, 1.73, 1.73, 1.73 /
      data a11/ 0.017, 0.017, 0.017, 0.018, 0.018, 0.018, 0.0175, 0.017, 
     1          0.016, 0.0152, 0.014, 0.013, 0.0122, 0.0113, 0.01, 0.0082,
     2          0.007, 0.006, 0.0052, 0.004, 0.003, 0.0022, 0.0013, 0. /

      data a12/ 0.818, 0.857, 0.92, 1.007, 1.225, 1.457, 1.849, 2.082, 2.24,
     1          2.341, 2.415, 2.359, 2.227, 1.949, 1.402, 0.329, -0.487,
     2          -0.77, -0.7, -0.607, -0.54, -0.479, -0.393, -0.35 /
      data a13/ -0.0135, -0.0135, -0.0135, -0.0138, -0.0142, -0.0145,
     1       -0.0153, -0.0162, -0.0172, -0.0183, -0.0206, -0.0231, -0.0256,
     2       -0.0296, -0.0363, -0.0493, -0.061, -0.0711, -0.0798, -0.0935,
     3       -0.098, -0.098, -0.098, -0.098 /
      data a14/ -0.223, -0.196, -0.128, -0.13, -0.13, -0.13, -0.156, 
     1     -0.172, -0.184, -0.194, -0.21, -0.223, -0.233, -0.245, -0.261,
     2     -0.285, -0.301, -0.313, -0.323, -0.282, -0.25, -0.25, -0.25, -0.25 /
      data a15 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a16 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a45 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
     
      data adjust_inter / 1.04, 1.05, 1.23, 1.34, 1.32, 1.32, 1.21, 1.14,
     1       1.05, 0.95, 0.79, 0.66, 0.54, 0.36, 0.24, -0.08, -0.21, -0.21,
     2      -0.22, -0.06, 0.06, 0.09, 0.14, 0.31/
      data adjust_slab / 0.83, 0.79, 0.71, 0.98, 0.99, 1.00, 0.92, 0.88, 
     1      0.81, 0.75, 0.62, 0.54, 0.46, 0.33, 0.23, -0.01, -0.09, -0.05,
     2      -0.02, 0.06, 0.06, 0.10, 0.13, 0.24 /
     
      data slabDC1/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
     1     -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.6, -0.6, -0.6 / 

      data sigs / 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 
     1            0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62,
     2            0.62, 0.62, 0.62, 0.62, 0.62, 0.62 /
      data sigt / 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.56, 0.54, 0.52, 
     1           0.505, 0.48, 0.46, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45,
     2           0.45, 0.45, 0.45, 0.45, 0.45, 0.45 /
     
      data deltaC1 / 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
     1    0.4, 0.4, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0 /

C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8 
 
C Find the requested spectral period and corresponding coefficients
      nPer = 24

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a4T = a4(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         a45T = a45(i1)
         adjust_slabT = adjust_slab(i1)
         adjust_interT = adjust_inter(i1)
         slabDC1T = slabDC1(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         DeltaC1T = deltaC1(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (6/2018 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)
            call S24_interp (period(count1),period(count2),a45(count1),a45(count2),
     +                   specT,a45T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),adjust_slab(count1),adjust_slab(count2),
     +                   specT,adjust_slabT,iflag)
            call S24_interp (period(count1),period(count2),adjust_inter(count1),adjust_inter(count2),
     +                   specT,adjust_interT,iflag)
            call S24_interp (period(count1),period(count2),slabDC1(count1),slabDC1(count2),
     +                   specT,slabDC1T,iflag)
            call S24_interp (period(count1),period(count2),deltaC1(count1),deltaC1(count2),
     +                   specT,deltaC1T,iflag)
     
 1011 period1 = specT         

c      write (*,'( /,''period ='',f10.3)') period1                                                                                                     

C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if ( ftype .eq. 0 .or. ftype .eq. 1 ) then
        R = rRup + c4*exp( (mag-6.0)*a9 ) 
        gm1 = a1T + (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R)
     1          + a6T*rRup + a10T*ftype + ftype * (slabDC1T*a4T)
      else
         write (*,*) 'BC Hydro V3 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
c      write (*,'( 10f8.4)') a1T, a2T, a14T, a3, A6T, a10T, a4T
c      write (*,'( 10f10.4)') (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R),
c     1   ftype * (slabDC1T*a4T) 
c      write (*,'( f10.4)') gm1
c      if ( specT .gt. 8. ) pause
      
C     Base model for Magnitude scaling.      
      if ( ftype .eq. 0. ) then
        testMag = (7.8 + deltaC1T)
      else
        testMag = (7.8 + deltaC1T)+slabDC1T
      endif
      
      if (mag .le. testmag ) then
        if ( ftype .eq. 0 ) then
          gm2 = (a4T+a45T)*(mag-testmag) + a13T*(10.0-mag)**2.0
        else
          gm2 = a4T*(mag-testmag) + a13T*(10.0-mag)**2.0
        endif
      else
         gm2 = a5*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif            
c      write (*,'( 2x,''C1'',4f10.4)') testMag, a4T, a45T, a13T, gm2
      
C     Depth Scaling
      if ( depth .lt. 100. ) then
        gm3 = a11T*(depth - 60.0)*ftype
      else
        gm3 = a11T*(100. - 60.0)*ftype
      endif

C     Forearc/Backarc scaling      
c      if (ftype .eq. 1) then
c         sumgm = sumgm + (a7T +a8T*alog(max(disthypo,85.0)/40.0))*faba
c      elseif (ftype .eq. 0) then   
c         sumgm = sumgm + (a15T +a16T*alog(max(rRup,100.0)/40.0))*faba
c      endif 

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         gm4 = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         gm4 = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif
c      write (*,'( 10f10.4)') specT, a12T, VSstar,vLinT, b_soilT, n, c, gm4
c      pause
      
c     add the cascadia adjustment      
      if ( ftype .eq. 0. ) then
        gm5 = adjust_interT
      else
        gm5 = adjust_slabT
      endif
c      write (*,'( 2f10.4)') adjust_interT, adjust_slabT

C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = gm1 + gm2 + gm3 + gm4 + gm5

c      write (*,'( 10f10.4)') specT, vs30, gm1, gm2, gm3, gm4, gm5, exp(lnSA-gm5), lnSA, exp(lnSA)
c      if (specT .gt. 8. ) then
c        pause
c      endif
      return
      end
      
c ----------------------------------------------------------------------
      subroutine S05_BCHydroSub2018_glob_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1T, slab1 )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=24)
      real a1(MAXPER), a2(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER)
      real a45(MAXPER), a4(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a5, a9 
      real a1T, a2T, a6T, a7T, a8T, a4T, a45T
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, faba, R, testmag, VsStar, depth, specT
      real c1p, adjust_inter(MAXPER), adjust_slab(MAXPER)
      real adjust_slabT, adjust_interT, gm1, gm2, gm3, gm4, gm5
      real slabDC1(MAXPER), slabDC1T
      real deltaC1(MAXPER), deltaC1T, slab1

      data period /  0.00, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.60, 0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.01 /
      data vLin / 865.1, 865.1, 865.1, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1, 
     1            503.0,  456.6,  430.3,  410.5, 400.0, 400.0, 400.0, 400.0, 
     2            400.0,  400.0,  400.0,  400.0, 400.0, 400.0 / 
      data b_soil / -1.186,-1.219, -1.273, -1.346, -1.471, -1.624, -1.931,
     2              -2.188, -2.381, -2.518, -2.657, -2.669, -2.599, -2.401,
     3              -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 0, 0 /
      data a1 / 3.466, 3.570, 3.807, 4.135, 4.499, 4.615, 4.554, 
     1        4.326, 4.130, 3.965, 3.643, 3.332, 3.040, 2.585, 1.982, 
     2        1.100, 0.460, -0.042, -0.461, -0.955, -1.385, -1.669,
     3        -2.017, -2.456 /

      data a2 / -1.044, -1.044, -1.08, -1.11, -1.11, -1.11, -1.084, 
     1          -1.027, -0.983, -0.947, -0.89, -0.845, -0.809, -0.76,
     2 -0.698, -0.612, -0.55, -0.501, -0.46, -0.455, -0.45, -0.45, -0.45, -0.45 /

      data a4/0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.62, 0.64, 0.66,
     1        0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 
     2        0.73, 0.78, 0.84, 0.93 /
      data a6 / -0.0069, -0.0069, -0.0069, -0.0076, -0.0078, -0.0077, 
     1          -0.0074, -0.0072, -0.007, -0.0068, -0.0064, -0.0061, 
     2          -0.0058, -0.0054, -0.005, -0.0046, -0.0044, -0.0044, 
     3          -0.0044, -0.0044, -0.0044, -0.0044, -0.0044, -0.0044 /

      data a7 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a8 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      data a10/ 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     1          1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     2          1.73, 1.73, 1.73, 1.73 /
      data a11/ 0.017, 0.017, 0.017, 0.018, 0.018, 0.018, 0.0175, 0.017, 
     1          0.016, 0.0152, 0.014, 0.013, 0.0122, 0.0113, 0.01, 0.0082,
     2          0.007, 0.006, 0.0052, 0.004, 0.003, 0.0022, 0.0013, 0. /

      data a12 /0.983, 1.033, 1.126, 1.318, 1.536, 1.646, 1.826,
     1          1.998, 2.157, 2.266, 2.360, 2.334, 2.217, 1.941, 1.426,
     2          0.428, -0.367, -0.684, -0.650, -0.596, -0.560, -0.533,
     3          -0.505, -0.450 /
      data a13/ -0.0135, -0.0135, -0.0135, -0.0138, -0.0142, -0.0145,
     1       -0.0153, -0.0162, -0.0172, -0.0183, -0.0206, -0.0231, -0.0256,
     2       -0.0296, -0.0363, -0.0493, -0.061, -0.0711, -0.0798, -0.0935,
     3       -0.098, -0.098, -0.098, -0.098 /
      data a14/ -0.223, -0.196, -0.128, -0.13, -0.13, -0.13, -0.156, 
     1     -0.172, -0.184, -0.194, -0.21, -0.223, -0.233, -0.245, -0.261,
     2     -0.285, -0.301, -0.313, -0.323, -0.282, -0.25, -0.25, -0.25, -0.25 /
      data a15 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a16 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a45 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
     
      data slabDC1/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
     1     -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.6, -0.6, -0.6 / 

      data sigs / 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 
     1            0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62,
     2            0.62, 0.62, 0.62, 0.62, 0.62, 0.62 /
      data sigt / 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.56, 0.54, 0.52, 
     1           0.505, 0.48, 0.46, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45,
     2           0.45, 0.45, 0.45, 0.45, 0.45, 0.45 /
     
      data deltaC1 / 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
     1    0.4, 0.4, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0 /

C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8 
 
C Find the requested spectral period and corresponding coefficients
      nPer = 24

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a4T = a4(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         a45T = a45(i1)
         adjust_slabT = adjust_slab(i1)
         adjust_interT = adjust_inter(i1)
         slabDC1T = slabDC1(i1) + 0.05
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         DeltaC1T = deltaC1(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (6/2018 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)
            call S24_interp (period(count1),period(count2),a45(count1),a45(count2),
     +                   specT,a45T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),adjust_slab(count1),adjust_slab(count2),
     +                   specT,adjust_slabT,iflag)
            call S24_interp (period(count1),period(count2),adjust_inter(count1),adjust_inter(count2),
     +                   specT,adjust_interT,iflag)
            call S24_interp (period(count1),period(count2),slabDC1(count1),slabDC1(count2),
     +                   specT,slabDC1T,iflag)
            call S24_interp (period(count1),period(count2),deltaC1(count1),deltaC1(count2),
     +                   specT,deltaC1T,iflag)
     
 1011 period1 = specT         

c      write (*,'( /,''period ='',f10.3)') period1                                                                                                     

c     adjust mag break for slab
      slabDC1T = slabDC1T + slab1 - 7.2

C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if ( ftype .eq. 0 .or. ftype .eq. 1 ) then
        R = rRup + c4*exp( (mag-6.0)*a9 ) 
        gm1 = a1T + (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R)
     1          + a6T*rRup + a10T*ftype + ftype * (slabDC1T*a4T)
      else
         write (*,*) 'BC Hydro V3 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
c      write (*,'( 10f8.4)') a1T, a2T, a14T, a3, A6T, a10T, a4T
c      write (*,'( 10f10.4)') (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R),
c     1   ftype * (slabDC1T*a4T) 
c      write (*,'( f10.4)') gm1
c      if ( specT .gt. 8. ) pause
      
C     Base model for Magnitude scaling.      
      if ( ftype .eq. 0. ) then
        testMag = (7.8 + deltaC1T)
      else
        testMag = (7.8 + deltaC1T)+slabDC1T
      endif
      
      if (mag .le. testmag ) then
        if ( ftype .eq. 0 ) then
          gm2 = (a4T+a45T)*(mag-testmag) + a13T*(10.0-mag)**2.0
        else
          gm2 = a4T*(mag-testmag) + a13T*(10.0-mag)**2.0
        endif
      else
         gm2 = a5*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif            
c      write (*,'( 2x,''C1'',4f10.4)') testMag, a4T, a45T, a13T, gm2
      
C     Depth Scaling
      if ( depth .lt. 100. ) then
        gm3 = a11T*(depth - 60.0)*ftype
      else
        gm3 = a11T*(100. - 60.0)*ftype
      endif

C     Forearc/Backarc scaling      
c      if (ftype .eq. 1) then
c         sumgm = sumgm + (a7T +a8T*alog(max(disthypo,85.0)/40.0))*faba
c      elseif (ftype .eq. 0) then   
c         sumgm = sumgm + (a15T +a16T*alog(max(rRup,100.0)/40.0))*faba
c      endif 

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         gm4 = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         gm4 = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif
c      write (*,'( 10f10.4)') specT, a12T, VSstar,vLinT, b_soilT, n, c, gm4
c      pause
      
C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = gm1 + gm2 + gm3 + gm4 + gm5

c      write (*,'( 10f10.4)') specT, vs30, gm1, gm2, gm3, gm4, gm5, exp(lnSA-gm5), lnSA, exp(lnSA)
c      if (specT .gt. 8. ) then
c        pause
c      endif
      return
      end
      
c ----------------------------------------------------------------------
      subroutine S05_BCHydroSub2018_SA_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, Ftype, iflag, faba, depth, disthypo, deltac1T, slab1 )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=24)
      real a1(MAXPER), a2(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER)
      real a45(MAXPER), a4(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a5, a9 
      real a1T, a2T, a6T, a7T, a8T, a4T, a45T
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, faba, R, testmag, VsStar, depth, specT
      real c1p, adjust_inter(MAXPER), adjust_slab(MAXPER)
      real adjust_slabT, adjust_interT, gm1, gm2, gm3, gm4, gm5
      real slabDC1(MAXPER), slabDC1T
      real deltaC1(MAXPER), deltaC1T, slab1

      data period /  0.00, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.60, 0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.01 /
      data vLin / 865.1, 865.1, 865.1, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1, 
     1            503.0,  456.6,  430.3,  410.5, 400.0, 400.0, 400.0, 400.0, 
     2            400.0,  400.0,  400.0,  400.0, 400.0, 400.0 / 
      data b_soil / -1.186,-1.219, -1.273, -1.346, -1.471, -1.624, -1.931,
     2              -2.188, -2.381, -2.518, -2.657, -2.669, -2.599, -2.401,
     3              -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 0, 0 /
      data a1 / 3.410, 3.576, 3.904, 4.106, 4.371, 4.487, 4.564, 4.421, 
     1          4.191, 4.004, 3.602, 3.260, 2.977, 2.494, 1.802, 0.927,
     3          0.230, -0.254, -0.661, -1.132, -1.538, -1.780, -2.080, -2.479 /
      data a6 / -0.00616, -0.00632, -0.00647, -0.00694, -0.00733, -0.00732, 
     1          -0.00701, -0.00689, -0.00631, -0.00591, -0.00517, -0.00463,
     2          -0.00468, -0.00429, -0.00367, -0.00375, -0.00356, -0.00388,
     3          -0.00389, -0.00402, -0.0041, -0.00422, -0.00417, -0.00389 /
      data a12/ 0.977, 1.017, 1.097, 1.204, 1.313, 1.407, 1.745, 2.056, 
     1          2.323, 2.571, 2.791, 2.706, 2.644, 2.357, 1.867, 0.880, 
     2          0.137, -0.130, -0.157, -0.164, -0.159, -0.151, -0.179, -0.167 /

      data a2 / -1.044, -1.044, -1.08, -1.11, -1.11, -1.11, -1.084, 
     1          -1.027, -0.983, -0.947, -0.89, -0.845, -0.809, -0.76,
     2 -0.698, -0.612, -0.55, -0.501, -0.46, -0.455, -0.45, -0.45, -0.45, -0.45 /

      data a4/0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.62, 0.64, 0.66,
     1        0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 
     2        0.73, 0.78, 0.84, 0.93 /
      data a7 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a8 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      data a10/ 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     1          1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     2          1.73, 1.73, 1.73, 1.73 /
      data a11/ 0.017, 0.017, 0.017, 0.018, 0.018, 0.018, 0.0175, 0.017, 
     1          0.016, 0.0152, 0.014, 0.013, 0.0122, 0.0113, 0.01, 0.0082,
     2          0.007, 0.006, 0.0052, 0.004, 0.003, 0.0022, 0.0013, 0. /

      data a13/ -0.0135, -0.0135, -0.0135, -0.0138, -0.0142, -0.0145,
     1       -0.0153, -0.0162, -0.0172, -0.0183, -0.0206, -0.0231, -0.0256,
     2       -0.0296, -0.0363, -0.0493, -0.061, -0.0711, -0.0798, -0.0935,
     3       -0.098, -0.098, -0.098, -0.098 /
      data a14/ -0.223, -0.196, -0.128, -0.13, -0.13, -0.13, -0.156, 
     1     -0.172, -0.184, -0.194, -0.21, -0.223, -0.233, -0.245, -0.261,
     2     -0.285, -0.301, -0.313, -0.323, -0.282, -0.25, -0.25, -0.25, -0.25 /
      data a15 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a16 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a45 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
     
      data slabDC1/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
     1     -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.6, -0.6, -0.6 / 

      data sigs / 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 
     1            0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62,
     2            0.62, 0.62, 0.62, 0.62, 0.62, 0.62 /
      data sigt / 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.56, 0.54, 0.52, 
     1           0.505, 0.48, 0.46, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45,
     2           0.45, 0.45, 0.45, 0.45, 0.45, 0.45 /
     
      data deltaC1 / 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
     1    0.4, 0.4, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0 /

C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8 
 
C Find the requested spectral period and corresponding coefficients
      nPer = 24

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a4T = a4(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         a45T = a45(i1)
         adjust_slabT = adjust_slab(i1)
         adjust_interT = adjust_inter(i1)
         slabDC1T = slabDC1(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         DeltaC1T = deltaC1(i1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (6/2018 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)
            call S24_interp (period(count1),period(count2),a45(count1),a45(count2),
     +                   specT,a45T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),adjust_slab(count1),adjust_slab(count2),
     +                   specT,adjust_slabT,iflag)
            call S24_interp (period(count1),period(count2),adjust_inter(count1),adjust_inter(count2),
     +                   specT,adjust_interT,iflag)
            call S24_interp (period(count1),period(count2),slabDC1(count1),slabDC1(count2),
     +                   specT,slabDC1T,iflag)
            call S24_interp (period(count1),period(count2),deltaC1(count1),deltaC1(count2),
     +                   specT,deltaC1T,iflag)
     
     
 1011 period1 = specT         

c      write (*,'( /,''period ='',f10.3)') period1                                                                                                     

c     adjust mag break for slab
      slabDC1T = slabDC1T + slab1 - 7.2

C     Compute the R term and base model based on either Rupture Distance 
c         (Interface events) of Hypocentral distance (Intraslab events). 
      if ( ftype .eq. 0 .or. ftype .eq. 1 ) then
        R = rRup + c4*exp( (mag-6.0)*a9 ) 
        gm1 = a1T + (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R)
     1          + a6T*rRup + a10T*ftype + ftype * (slabDC1T*a4T)
      else
         write (*,*) 'BC Hydro V3 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
c      write (*,'( 10f8.4)') a1T, a2T, a14T, a3, A6T, a10T, a4T
c      write (*,'( 10f10.4)') (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R),
c     1   ftype * (slabDC1T*a4T) 
c      write (*,'( f10.4)') gm1
c      if ( specT .gt. 8. ) pause
      
C     Base model for Magnitude scaling.      
      if ( ftype .eq. 0. ) then
        testMag = (7.8 + deltaC1T)
      else
        testMag = (7.8 + deltaC1T)+slabDC1T
      endif
      
      if (mag .le. testmag ) then
        if ( ftype .eq. 0 ) then
          gm2 = (a4T+a45T)*(mag-testmag) + a13T*(10.0-mag)**2.0
        else
          gm2 = a4T*(mag-testmag) + a13T*(10.0-mag)**2.0
        endif
      else
         gm2 = a5*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif            
c      write (*,'( 2x,''C1'',4f10.4)') testMag, a4T, a45T, a13T, gm2
      
C     Depth Scaling
      if ( depth .lt. 100. ) then
        gm3 = a11T*(depth - 60.0)*ftype
      else
        gm3 = a11T*(100. - 60.0)*ftype
      endif

C     Forearc/Backarc scaling      
c      if (ftype .eq. 1) then
c         sumgm = sumgm + (a7T +a8T*alog(max(disthypo,85.0)/40.0))*faba
c      elseif (ftype .eq. 0) then   
c         sumgm = sumgm + (a15T +a16T*alog(max(rRup,100.0)/40.0))*faba
c      endif 

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         gm4 = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         gm4 = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif
c      write (*,'( 10f10.4)') specT, a12T, VSstar,vLinT, b_soilT, n, c, gm4
c      pause
      
C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = gm1 + gm2 + gm3 + gm4 + gm5

c      write (*,'( 10f10.4)') specT, vs30, gm1, gm2, gm3, gm4, gm5, exp(lnSA-gm5), lnSA, exp(lnSA)
c      if (specT .gt. 8. ) then
c        pause
c      endif
      return
      end
      
c ------------------------------------------------------------------            

      subroutine S05_NGASUB_GA2019 ( mag, fType, rRup, vs30, lnSa, sigma1, 
     2           specT, period1, iflag, forearc, depth, disthypo, deltac1,
     3           iRegion, z25 )

      implicit none
     
      real mag, fType, rRup, vs30, pgaRock, faba, vs30_rock, period0,
     1     lnSa, sigma, tau, period1, sigma1, disthypo, deltac1,
     2     depth, specT, slab1, z25
      integer iflag, forearc, iRegion
      real slabMag(7)
      
      data slabMag/ 7.95, 7.2, 7.6, 7.65, 7.6, 7.3, 7.7 /

c     Ftype defines an interface event or intraslab events      
C     fType    Event Type
C     -------------------
C      0       Interface  - use rupture distance
C      1       Intraslab  - use hypocentral distance
C
C     faba     Note
C     -------------------------
C      0       Forearc site  
C      1       Backarc site  
C

c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1000.
      faba = real(forearc)
      slab1 = slabMag(iRegion)

C     Compute Rock PGA
      call S05_AG2019_SUB_model ( mag, rRup, vs30_rock, pgaRock, lnSa, sigma, tau, 
     2        period0, Ftype, iflag, faba, depth, disthypo, z25, iRegion, slab1 )
      pgaRock = exp(lnSa)
 
C     Compute regular ground motions. 
      call S05_AG2019_SUB_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2        specT, Ftype, iflag, faba, depth, disthypo, z25, iRegion, slab1 )

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
      period1 = specT

c     Convert units spectral acceleration in gal                                
      lnSa = lnSa + 6.89                                                
      return
      end
  
c ------------------------------------------------------------------            

      subroutine S05_AG2019_SUB_model ( mag, rRup, vs30, pgaRock, lnSa, sigma, tau, 
     2        specT, Ftype, iflag, faba, depth, disthypo, z25, iRegion, slab1 )

      implicit none
      
      integer MAXPER, nPer, i1, i      
      parameter (MAXPER=24)
      real a1(MAXPER), a2(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a16(MAXPER)
      real a45(MAXPER), a4(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), sigs(MAXPER), sigt(MAXPER)
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo,
     1     mag, a3, a5, a9 
      real a1T, a2T, a6T, a7T, a8T, a4T, a45T
      real a10T, a11T, a12T, a13T, a14T, a15T, a16T, sigsT, sigtT
      real a17T, a18T, a19T, a20T, a21T, a22T, a23T, a24T, a25T, 
     1     a26T, a27T, a28T, a29T, a30T, a31T, a32T, a33T, a34T,
     1      a35T, a36T, a37T, a39T, a41T 
      real vLinT, b_soilT, sumgm, Ftype, tau, period1
      integer count1, count2, iflag
      real n, c, c4, c1, faba, R, testmag, VsStar, depth, specT
      real c1p, adjust_inter(MAXPER), adjust_slab(MAXPER)
      real adjust_slabT, adjust_interT, gm1, gm2, gm3, gm4, gm5
      real slabDC1(MAXPER), slabDC1T
      real deltaC1(MAXPER), deltaC1T
      real a17(MAXPER), a18(MAXPER), a19(MAXPER), a20(MAXPER), a21(MAXPER),
     1     a22(MAXPER), a23(MAXPER)
      real a24(MAXPER), a25(MAXPER), a26(MAXPER), a27(MAXPER), a28(MAXPER),
     1     a29(MAXPER), a30(MAXPER)
      real a31(MAXPER), a32(MAXPER), a33(MAXPER), a34(MAXPER), a35(MAXPER),
     1     a36(MAXPER), a37(MAXPER)
      real a39(MAXPER), a41(MAXPER)
      real CAS_adjust_slab(MAXPER), CAS_adjust_int(MAXPER)
      real z25, regTerm1, regTerm2, regTerm3, regTerm4, temp1 , z25Ref
      integer iRegion
      real slab1, z25_m

      data period /  0.00, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
     1               0.60, 0.75, 1.00,  1.5, 2.00, 2.5, 3.00, 4.0, 5.0, 
     2               6.0, 7.5, 10.01 /
      data vLin / 865.1, 865.1, 865.1, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1, 
     1            503.0,  456.6,  430.3,  410.5, 400.0, 400.0, 400.0, 400.0, 
     2            400.0,  400.0,  400.0,  400.0, 400.0, 400.0 / 
      data b_soil / -1.186,-1.219, -1.273, -1.346, -1.471, -1.624, -1.931,
     2              -2.188, -2.381, -2.518, -2.657, -2.669, -2.599, -2.401,
     3              -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 0, 0 /
      data a2 / -1.044, -1.044, -1.08, -1.11, -1.11, -1.11, -1.084, 
     1          -1.027, -0.983, -0.947, -0.89, -0.845, -0.809, -0.76,
     2 -0.698, -0.612, -0.55, -0.501, -0.46, -0.455, -0.45, -0.45, -0.45, -0.45 /

      data a4/0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.59, 0.62, 0.64, 0.66,
     1        0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 0.68, 
     2        0.73, 0.78, 0.84, 0.93 /
      data a6/ -0.0069, -0.0069, -0.0069, -0.0076, -0.0078, -0.0077, 
     1         -0.0074, -0.0072, -0.007, -0.0068, -0.0064, -0.0061, 
     2         -0.0058, -0.0054, -0.005, -0.0046, -0.0044, -0.0044, 
     2         -0.0044, -0.0044, -0.0044, -0.0044, -0.0044, -0.0044 /

      data a7 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a8 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      data a10/ 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     1          1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 1.73, 
     2          1.73, 1.73, 1.73, 1.73 /
      data a11/ 0.017, 0.017, 0.017, 0.018, 0.018, 0.018, 0.0175, 0.017, 
     1          0.016, 0.0152, 0.014, 0.013, 0.0122, 0.0113, 0.01, 0.0082,
     2          0.007, 0.006, 0.0052, 0.004, 0.003, 0.0022, 0.0013, 0. /

      data a12/ 0.983, 1.033, 1.126, 1.318, 1.536, 1.646, 1.826, 1.998,
     1          2.157, 2.266, 2.36, 2.334, 2.217, 1.941, 1.426, 0.428, 
     2         -0.367, -0.684, -0.65, -0.596, -0.56, -0.533, -0.505, -0.45 /
     
      data a13/ -0.0135, -0.0135, -0.0135, -0.0138, -0.0142, -0.0145,
     1       -0.0153, -0.0162, -0.0172, -0.0183, -0.0206, -0.0231, -0.0256,
     2       -0.0296, -0.0363, -0.0493, -0.061, -0.0711, -0.0798, -0.0935,
     3       -0.098, -0.098, -0.098, -0.098 /
     
      data a14/ -0.223, -0.196, -0.128, -0.13, -0.13, -0.13, -0.156, 
     1     -0.172, -0.184, -0.194, -0.21, -0.223, -0.233, -0.245, -0.261,
     2     -0.285, -0.301, -0.313, -0.323, -0.282, -0.25, -0.25, -0.25, -0.25 /
     
      data a15 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a16 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a45 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      data a17 / 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0.,
     1           0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  /


      data a18 / -0.165, -0.176, -0.205, -0.311, -0.311, -0.189, 0.023, 
     1           0.084, 0.083, 0.075, 0.055, 0.025, 0.01, 0.008, -0.024,
     2          -0.099, -0.12, -0.086, -0.05, -0.011, 0.02, 0.054, 0.112, 0.1 /

      data a19 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.01, -0.01, 0.008, 
     1       0.011, -0.001, -0.005, 0, 0.01, 0, 0, 0, 0.027, 0 /

      data a20 / -0.122, -0.116, -0.103, -0.091, -0.041, -0.019, -0.077,
     1          -0.116, -0.147, -0.165, -0.175, -0.177, -0.165, -0.127,
     2          -0.124, -0.155, -0.175, -0.176, -0.165, -0.146, -0.135,
     3          -0.106, -0.048, -0.06 /

      data a21 / 0.049, 0.047, 0.037, -0.041, -0.069, -0.015, 0.153, 
     1           0.204, 0.203, 0.205, 0.215, 0.205, 0.21, 0.228, 0.216,
     2           0.156, 0.13, 0.144, 0.15, 0.144, 0.13, 0.124, 0.131, 0.008 /

      data a22 / 0.049, 0.047, 0.037, -0.041, -0.069, -0.015, 0.153, 
     1           0.204, 0.203, 0.205, 0.215, 0.205, 0.21, 0.228, 0.216, 
     2           0.156, 0.13, 0.144, 0.15, 0.144, 0.13, 0.124, 0.131, 0.008 /

      data a23/ -0.046, -0.05, -0.064, -0.131, -0.175, -0.116, 0.038, 
     1           0.104, 0.145, 0.165, 0.185, 0.185, 0.18, 0.184, 0.146,
     2           0.046, -0.013, -0.016, -0.001, 0.014, 0.03, 0.054, 0.112, 0.1 /

      data a24 / 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0.,
     1           0., 0., 0., 0., 0., 0., 0., 0., 0., 0.  /

      data a25 / 0.0017, 0.0017, 0.0016, 0.0022, 0.0021, 0.0018, 0.0012,
     1           0.0009, 0.0007, 0.0005, 0.0002, 0, -0.0001, -0.0002, 
     2          -0.0003, 0, 0.0004, 0.0009, 0.0012, 0.0016, 0.0019, 
     3          0.002, 0.0022, 0.0024 /

      data a26 / 0.0002, -0.0003, -0.0011, -0.0012, -0.0014, -0.001, 
     1           -0.0005, 0.0002, 0.0003, 0.0003, 0.0001, 0, -0.0001, 
     2          -0.0003, -0.0002, -0.0002, -0.0002, 0, 0.0002, 0.0003,
     1           0.0004, 0.0005, 0.0006, 0.0007 /

      data a27 / -0.0036, -0.0036, -0.0037, -0.004, -0.0045, -0.0047, 
     1          -0.0045, -0.004, -0.0035, -0.0032, -0.0025, -0.002, 
     1          -0.0016, -0.0014, -0.0012, -0.0011, -0.0008, -0.0006, 
     1          -0.0006, -0.0004, -0.0003, -0.0002, -0.0001, 0. /

      data a28 / 0.0008, 0.0008, 0.0006, 0.0011, 0.0007, 0.0005, 0.0003,
     1           0.0004, 0.0004, 0.0003, 0.0004, 0.0004, 0.0004, 0.0002,
     1           0.0002, 0.0002, 0.0002, 0.0004, 0.0006, 0.0007, 0.0009,
     1           0.001, 0.0011, 0.0012 /

      data a29 / 0.0013, 0.0012, 0.001, 0.0013, 0.0012, 0.0008, 0.0006, 
     1          0.0007, 0.001, 0.0012, 0.0016, 0.0018, 0.0017, 0.0017, 
     1          0.0017, 0.0013, 0.0011, 0.0011, 0.0011, 0.0011, 0.0011,
     1           0.0011, 0.0011, 0.0011 /

      data a30 / 0.0064, 0.0064, 0.0064, 0.0071, 0.0073, 0.0072, 0.0069,
     1 0.0067, 0.0062, 0.0056, 0.0047, 0.004, 0.0035, 0.003, 0.0027, 
     2 0.0025, 0.0026, 0.0027, 0.0027, 0.0027, 0.0027, 0.0027, 0.0027, 0.0027 /
     
      data a31 / 3.466, 3.570, 3.807, 4.230, 4.499, 4.615, 4.554, 4.326, 
     1           4.130, 3.965, 3.643, 3.332, 3.040, 2.585, 1.982, 1.100,
     2           0.460, -0.042, -0.461, -0.955, -1.385, -1.669, -2.017, -2.456  /

      data a32 /1.813, 1.824, 1.909, 1.996, 2.207, 2.457, 2.797, 2.742,
     1   2.659, 2.600, 2.486, 2.399, 2.240, 1.930, 1.562, 0.925, 0.385,
     2   -0.227, -0.748, -1.460, -1.984, -2.257, -2.641, -3.218 /
      data a33 /3.014, 3.183, 3.539, 3.833, 4.106, 4.119, 4.010, 3.792,
     1   3.613, 3.478, 3.244, 3.002, 2.729, 2.299, 1.780, 0.984, 0.403, 
     2  -0.079, -0.509, -1.041, -1.440, -1.808, -2.378, -2.849 /
      data a34 /3.871, 3.999, 4.355, 4.787, 5.207, 5.387, 5.185, 4.814,
     1   4.512, 4.286, 3.834, 3.442, 3.094, 2.614, 2.015, 1.131, 0.454,
     2   -0.066, -0.488, -0.973, -1.403, -1.706, -2.093, -2.551 /
      data a35 /3.291, 3.351, 3.646, 3.790, 3.999, 4.187, 4.365, 4.270, 
     1  4.123, 3.976, 3.690, 3.379, 3.094, 2.672, 2.003, 1.067, 0.473, 
     2  -0.078, -0.456, -0.867, -1.221, -1.587, -1.951, -2.328 /
      data a36 /3.305, 3.473, 3.881, 4.147, 4.367, 4.569, 4.617, 4.381, 
     1   4.140, 3.931, 3.545, 3.222, 2.935, 2.449, 1.809, 0.934, 0.306, 
     2   -0.232, -0.681, -1.136, -1.526, -1.808, -2.114, -2.405 /
      data a37 /2.439, 2.490, 2.631, 2.813, 2.957, 3.158, 3.369, 3.282, 
     1   3.224, 3.193, 3.075, 2.907, 2.737, 2.394, 1.804, 1.059, 0.380,
     2   -0.130, -0.546, -1.097, -1.577, -1.832, -2.159, -2.681 /
     
      data a39 / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0.04,
     1           0.05, 0.06, 0.07, 0.07, 0.07, 0.07, 0.07 /

      data a41 / -0.034, -0.035, -0.04, -0.058, -0.064, -0.065, -0.062,
     1           -0.046, -0.029, -0.013, 0.012, 0.041, 0.067, 0.096, 
     2            0.133, 0.184, 0.21, 0.221, 0.229, 0.228, 0.218, 0.204,
     3             0.188, 0.143 /
     
      data a39 / 0.27, 0.26, 0.22, 0.19, 0.15, 0.15, 0.2, 0.26, 0.3, 
     1           0.34, 0.37, 0.4, 0.42, 0.44, 0.48, 0.54, 0.58, 0.61, 
     1           0.63, 0.65, 0.65, 0.65, 0.6, 0.55 /

      data a41 / -0.07, -0.07, -0.08, -0.11, -0.14, -0.12, -0.1, -0.08,
     1            -0.06, -0.04, 0.01, 0.07, 0.12, 0.17, 0.23, 0.32, 0.37, 
     1           0.38, 0.39, 0.4, 0.38, 0.35, 0.32, 0.27 /
 
      data adjust_inter /1.423, 1.422, 1.578, 1.810, 1.749, 1.759, 1.692,
     1 1.493, 1.394, 1.301, 1.120, 0.919, 0.804, 0.658, 0.408, 0.157, 
     2 0.028, 0.101, 0.192, 0.401, 0.503, 0.459, 0.446, 0.592 /
      data adjust_slab /1.058, 1.058, 1.210, 1.310, 1.280, 1.286, 1.240, 
     1   1.089, 1.024, 0.960, 0.831, 0.679, 0.592, 0.494, 0.306, 0.114, 
     2  0.023, 0.084, 0.160, 0.334, 0.419, 0.383, 0.371, 0.493 /
     
      data slabDC1/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
     1     -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.6, -0.6, -0.6 / 

      data sigs / 0.536, 0.538, 0.545, 0.576, 0.611, 0.617, 0.603, 0.587,
     1     0.583, 0.572, 0.567, 0.562, 0.56, 0.559, 0.555, 0.563, 0.576,
     2     0.569, 0.553, 0.516, 0.47, 0.447, 0.421, 0.388 /
      data sigt / 0.471, 0.472, 0.486, 0.532, 0.547, 0.528, 0.497, 0.473, 
     1     0.489, 0.477, 0.478, 0.459, 0.46, 0.463, 0.473, 0.505, 0.512,
     2     0.499, 0.486, 0.486, 0.449, 0.399, 0.399, 0.401 /

      data deltaC1 / 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 
     1    0.4, 0.4, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0 /

C Constant parameters            
      n = 1.18
      c = 1.88
      a3 = 0.1
      a5 = 0.0
      a9 = 0.4
      c4 = 10.0
      c1 = 7.8 
 
C Find the requested spectral period and corresponding coefficients
      nPer = 24

C First check for the PGA case 
      if (specT .eq. 0.0) then
         i1=1
         period1 = period(i1)
c         a1T = a1(i1)
         a2T = a2(i1)
         a4T = a4(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         a45T = a45(i1)
         adjust_slabT = adjust_slab(i1)
         adjust_interT = adjust_inter(i1)
         slabDC1T = slabDC1(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         sigtT = sigt(i1)
         sigsT = sigs(i1)
         DeltaC1T = deltaC1(i1)
         a17T = a17(i1)
         a18T = a18(i1)
         a19T = a19(i1)
         a20T = a20(i1)
         a21T = a21(i1)
         a22T = a22(i1)
         a23T = a23(i1)
         a24T = a24(i1)
         a25T = a25(i1)
         a26T = a26(i1)
         a27T = a27(i1)
         a28T = a28(i1)
         a29T = a29(i1)
         a30T = a30(i1)
         a31T = a31(i1)
         a32T = a32(i1)
         a33T = a33(i1)
         a34T = a34(i1)
         a35T = a35(i1)
         a36T = a36(i1)
         a37T = a37(i1)
         a39T = a39(i1)
         a41T = a41(i1)

         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'BCHydro Sub (6/2018 Model) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: ' 
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       continue
c           call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
c     +                   specT,a1T,iflag)
            call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
            call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,a6T,iflag)
            call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,a7T,iflag)
            call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,a8T,iflag)
            call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,a10T,iflag)
            call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,a11T,iflag)
            call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,a12T,iflag)
            call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,a13T,iflag)
            call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,a14T,iflag)
            call S24_interp (period(count1),period(count2),a15(count1),a15(count2),
     +                   specT,a15T,iflag)
            call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,a16T,iflag)
            call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soilT,iflag)
            call S24_interp (period(count1),period(count2),vLin(count1),vLin(count2),
     +                   specT,vLinT,iflag)
            call S24_interp (period(count1),period(count2),sigs(count1),sigs(count2),
     +                   specT,sigsT,iflag)
            call S24_interp (period(count1),period(count2),sigt(count1),sigt(count2),
     +                   specT,sigtT,iflag)
            call S24_interp (period(count1),period(count2),a45(count1),a45(count2),
     +                   specT,a45T,iflag)
            call S24_interp (period(count1),period(count2),a4(count1),a4(count2),
     +                   specT,a4T,iflag)
            call S24_interp (period(count1),period(count2),adjust_slab(count1),adjust_slab(count2),
     +                   specT,adjust_slabT,iflag)
            call S24_interp (period(count1),period(count2),adjust_inter(count1),adjust_inter(count2),
     +                   specT,adjust_interT,iflag)
            call S24_interp (period(count1),period(count2),slabDC1(count1),slabDC1(count2),
     +                   specT,slabDC1T,iflag)
            call S24_interp (period(count1),period(count2),deltaC1(count1),deltaC1(count2),
     +                   specT,deltaC1T,iflag)

            call S24_interp (period(count1),period(count2),a17(count1),a17(count2),
     +                   specT,a17T,iflag)
            call S24_interp (period(count1),period(count2),a18(count1),a18(count2),
     +                   specT,a18T,iflag)
            call S24_interp (period(count1),period(count2),a19(count1),a19(count2),
     +                   specT,a19T,iflag)
            call S24_interp (period(count1),period(count2),a20(count1),a20(count2),
     +                   specT,a20T,iflag)
            call S24_interp (period(count1),period(count2),a21(count1),a21(count2),
     +                   specT,a21T,iflag)
            call S24_interp (period(count1),period(count2),a22(count1),a22(count2),
     +                   specT,a22T,iflag)
            call S24_interp (period(count1),period(count2),a23(count1),a23(count2),
     +                   specT,a23T,iflag)
            call S24_interp (period(count1),period(count2),a24(count1),a24(count2),
     +                   specT,a24T,iflag)
            call S24_interp (period(count1),period(count2),a25(count1),a25(count2),
     +                   specT,a25T,iflag)
            call S24_interp (period(count1),period(count2),a26(count1),a26(count2),
     +                   specT,a26T,iflag)
            call S24_interp (period(count1),period(count2),a27(count1),a27(count2),
     +                   specT,a27T,iflag)
            call S24_interp (period(count1),period(count2),a28(count1),a28(count2),
     +                   specT,a28T,iflag)
            call S24_interp (period(count1),period(count2),a29(count1),a29(count2),
     +                   specT,a29T,iflag)
            call S24_interp (period(count1),period(count2),a30(count1),a30(count2),
     +                   specT,a30,iflag)
            call S24_interp (period(count1),period(count2),a31(count1),a31(count2),
     +                   specT,a31T,iflag)
            call S24_interp (period(count1),period(count2),a32(count1),a32(count2),
     +                   specT,a32T,iflag)
            call S24_interp (period(count1),period(count2),a33(count1),a33(count2),
     +                   specT,a33T,iflag)
            call S24_interp (period(count1),period(count2),a34(count1),a34(count2),
     +                   specT,a34T,iflag)
            call S24_interp (period(count1),period(count2),a35(count1),a35(count2),
     +                   specT,a35T,iflag)
            call S24_interp (period(count1),period(count2),a36(count1),a36(count2),
     +                   specT,a36T,iflag)
            call S24_interp (period(count1),period(count2),a37(count1),a37(count2),
     +                   specT,a37T,iflag)
            call S24_interp (period(count1),period(count2),a39(count1),a39(count2),
     +                   specT,a39T,iflag)
            call S24_interp (period(count1),period(count2),a41(count1),a41(count2),
     +                   specT,a41T,iflag)
     
 1011 period1 = specT         

c     adjust mag break for slab
      slabDC1T = slabDC1T + slab1 - 7.2

C     Compute the R term and base model 
      if ( ftype .eq. 0 .or. ftype .eq. 1 ) then
        R = rRup + c4*exp( (mag-6.0)*a9 ) 
        gm1 = (a2T + a14T*ftype + a3*(mag -7.8 ))*alog(R)
     1          + a6T*rRup + a10T*ftype + ftype * (slabDC1T*a4T)
      else
         write (*,*) 'BC Hydro V3 Model not defined for Ftype'
         write (*,*) 'other than 0 (interface) or 1 (intraslab)'
         stop 99
      endif
      
C     Base model for Magnitude scaling.      
      if ( ftype .eq. 0. ) then
        testMag = (7.8 + deltaC1T)
      else
        testMag = (7.8 + deltaC1T)+slabDC1T
      endif
      
      if (mag .le. testmag ) then
        if ( ftype .eq. 0 ) then
          gm2 = (a4T+a45T)*(mag-testmag) + a13T*(10.0-mag)**2.0
        else
          gm2 = a4T*(mag-testmag) + a13T*(10.0-mag)**2.0
        endif
      else
         gm2 = a5*(mag-testmag) + a13T*(10.0-mag)**2.0
      endif            
c      write (*,'( 2x,''C1'',4f10.4)') testMag, a4T, a45T, a13T, gm2
      
C     Depth Scaling
      if ( depth .lt. 100. ) then
        gm3 = a11T*(depth - 60.0)*ftype
      else
        gm3 = a11T*(100. - 60.0)*ftype
      endif

C     Forearc/Backarc scaling      
c      if (ftype .eq. 1) then
c         sumgm = sumgm + (a7T +a8T*alog(max(disthypo,85.0)/40.0))*faba
c      elseif (ftype .eq. 0) then   
c         sumgm = sumgm + (a15T +a16T*alog(max(rRup,100.0)/40.0))*faba
c      endif 

C     Site Response 
      if (vs30 .ge. 1000.0) then
          VsStar = 1000.0
      else
          VsStar = vs30
      endif
       
      if (vs30 .ge. VlinT) then
         gm4 = a12T*alog(VsStar/vLinT) + b_soilT*n*alog(VsStar/vLinT)
      else
         gm4 = a12T*alog(VsStar/vLinT) - b_soilT*alog(pgarock + c) +
     1          b_soilT*alog(pgarock + c*(VsStar/vlinT)**n)     
      endif
c      write (*,'( 10f10.4)') specT, pgarock, a12T, VSstar,vLinT, b_soilT, n, c, gm4
c      pause
      
c     add the cascadia adjustment      
      if ( iRegion .eq. 2 ) then
        if ( ftype .eq. 0. ) then
         gm5 = adjust_interT
        else
         gm5 = adjust_slabT
        endif
      else
        gm5 = 0.
      endif
c      write (*,'( 2f10.4)') adjust_interT, adjust_slabT

C     Set sigma values to return
      sigma = sigsT
      tau = sigtT

c     Set SA to return
      lnSa = gm1 + gm2 + gm3 + gm4 + gm5
      
      z25_m = z25 * 1000.
      
c     Add the regional terms
      if (iRegion .eq. 1 ) then
        regTerm1 = a17T*alog(VsStar/vLinT)
        regTerm2 = a24T * Rrup
        regTerm3 = a31T
        regTerm4 = 0.
      elseif (iRegion .eq. 2 ) then  
        regTerm1 = a18T*alog(VsStar/vLinT)
        regTerm2 = a25T * Rrup
        regTerm3 = a32T
        call S05_Calc_z25Ref ( vs30, 200., 570., 8.6, 6.9, z25ref) 
        temp1 = alog ( ( z25_m + 50.) /(z25ref+50.) )
        if ( temp1 .gt. 0.5 ) then
         regTerm4 = a39T * (temp1-0.5)
        else
          regTerm4 = 0.
        endif
      elseif (iRegion .eq. 3 ) then  
        regTerm1 = a19T*alog(VsStar/vLinT)
        regTerm2 = a26T * Rrup
        regTerm3 = a33T
        regTerm4 = 0.
      elseif (iRegion .eq. 4 ) then  
        regTerm1 = a20T*alog(VsStar/vLinT)
        regTerm2 = a27T * Rrup
        regTerm3 = a34T
        call S05_Calc_z25Ref ( vs30, 200., 800., 7.0, 4.1, z25ref) 
        temp1 = alog ( ( z25_m + 50.) /(z25ref+50.) )
        regTerm4 = a41T * temp1
      elseif (iRegion .eq. 5 ) then  
        regTerm1 = a21T*alog(VsStar/vLinT)
        regTerm2 = a28T * Rrup
        regTerm3 = a35T
        regTerm4 = 0.
      elseif (iRegion .eq. 6 ) then  
        regTerm1 = a22T*alog(VsStar/vLinT)
        regTerm2 = a29T * Rrup
        regTerm3 = a36T
        regTerm4 = 0.
      elseif (iRegion .eq. 7 ) then  
        regTerm1 = a23T*alog(VsStar/vLinT)
        regTerm2 = a30T * Rrup
        regTerm3 = a37T
        regTerm4 = 0.
      endif
      
      lnSa = lnSa + regTerm1 + regTerm2 + regTerm3 + regTerm4

c      write (*,'( 10f10.4)') specT, vs30, gm1, gm2, gm3, gm4, gm5, exp(lnSA-gm5), lnSA, exp(lnSA)
c      if (specT .gt. 8. ) then
c        pause
c      endif
      return
      end

c -------------------------------------------------------------------
	  
      subroutine S05_Calc_z25Ref ( vs30, x1, x2, y1, y2, y) 
      real vs30, x1, x2, y1, y2, y, x, x1Log, x2Log
      
      x = alog(vs30)
      x1Log = alog(x1)
      x2Log = alog(x2)
      if ( x .lt. x1Log ) then
        y = y1
      elseif ( x .gt. x2Log ) then
        y = y2
      else
        y = (x-x1Log)/(x2Log-x1Log) * (y2-y1) + y1
      endif
      y = exp(y)
      
      return
      end

