
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
      sumgm = sumgm + a11T*(depth - 60.0)*ftype

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

