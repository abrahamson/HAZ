!DEC$ FIXEDFORMLINESIZE:132
C Empirical Attenuation Models 
C Last Modified: 7/17/2015
C with haz43b
c ------------------------------------------------------------
C *** New Lin et al. (2011) Taiwan crustal hangingwall model SOIL (jClac=316)************
c ------------------------------------------------------------
      subroutine Lin_hw_soil(mag, rupDist, specT, period1, lnY, sigma, iflag )

      parameter (MAXPER=16)
      REAL Period(MAXPER), C1(MAXPER), C2(MAXPER), C3(MAXPER), C4(MAXPER), C5(MAXPER)
      real mag, rupDist, lnY, sigma, period1, sig(MAXPER)                                  
      real specT, c1T,  c2T,  c3T, c4T,  c5T
      integer nper, count1, count2,iflag
      
C.....MODEL COEFFICIENTS.....................
      data C1  / -3.248,-3.008,-1.994,-1.408, -1.508, -3.226, -4.050, -5.293, 
     +           -6.307,-7.209,-8.309,-9.868,-11.216,-12.806,-13.886,-14.606/

      data  C2 /  0.943,0.905,0.809,0.765,0.785,0.870,0.999,1.165, 
     +            1.291,1.395,1.509,1.691,1.798,2.005,2.099,2.160/ 

      data  C3 /  -1.471,-1.451,-1.500,-1.551,-1.551,-1.211,-1.205,-1.167, 
     +            -1.134,-1.099,-1.044,-1.004,-0.965,-0.975,-1.007,-1.114/

      data  C4 /  0.1000,0.1100,0.2510,0.2800,0.2800,0.0450,0.0300,0.0110, 
     +            0.0042,0.0016,0.0006,0.0004,0.0003,0.0005,0.0004,0.0004/  

      data  C5 /  0.648,0.638,0.518,0.510,0.500,0.708,0.788,0.958,                        
     +            1.118,1.258,1.408,1.485,1.522,1.528,1.548,1.562/                        

      data  sig / 0.627, 0.622, 0.686, 0.709, 0.713, 0.686, 0.656, 0.655,               
     +            0.653, 0.642, 0.650, 0.676, 0.721, 0.759, 0.787, 0.820/                        

      data  period/   .00, .01, .06,  .09,  .10,  .20,  .30,  .40,                  
     +                .50, .60, .75, 1.00, 1.50, 2.00, 3.00, 5.00/   
C Find the requested spectral period and corresponding coefficients
      nPer = 16

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

      write (*,*) 
      write (*,*) 'Lin (2011) crustal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c2(count1),c2(count2)
     +             ,specT,c2T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c5(count1),c5(count2)
     +             ,specT,c5T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

      lnY = C1T + C2T*mag + 
     1      C3T*alog(rupdist+C4T*exp(C5T*mag))
     
      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end                                                                       

c ------------------------------------------------------------
C *** New Lin et al. (2011) Taiwan crustal hangingwall model ROCK (jClac=317)************
c ------------------------------------------------------------
      subroutine Lin_hw_rock(mag, rupDist, specT, period1, lnY, sigma, iflag )

      parameter (MAXPER=16)
      REAL Period(MAXPER), C1(MAXPER), C2(MAXPER), C3(MAXPER), C4(MAXPER), C5(MAXPER)
      real mag, rupDist, lnY, sigma, period1, sig(MAXPER)                                  
      real specT, c1T,  c2T,  c3T, c4T,  c5T
      integer nper, count1, count2,iflag
      
C.....MODEL COEFFICIENTS.....................
      data C1  / -3.279,-3.253,-1.738,-1.237, -1.103, -2.767, -4.440, -5.630, 
     +           -6.746,-7.637,-8.641,-9.978,-11.617,-12.611,-13.303,-13.914/

      data  C2 /  1.035,1.018,0.908,0.841,0.841,0.980,1.186,1.335,
     +            1.456,1.557,1.653,1.800,1.976,2.058,2.036,1.958/

      data  C3 /  -1.651,-1.629,-1.769,-1.750,-1.765,-1.522,-1.438,-1.414, 
     +            -1.365,-1.348,-1.313,-1.286,-1.284,-1.261,-1.234,-1.156/

      data  C4 /  0.1520,0.1596,0.3270,0.4780,0.4550,0.0970,0.0275,0.0140, 
     +            0.0060,0.0033,0.0015,0.0008,0.0004,0.0005,0.0013,0.0012/  

      data  C5 /  0.623,0.612,0.502,0.402,0.417,0.627,0.823,0.932,                           
     +            1.057,1.147,1.257,1.377,1.508,1.497,1.302,1.241/                           

      data  sig / 0.651, 0.647, 0.701, 0.747, 0.750, 0.696, 0.685, 0.682,               
     +            0.678, 0.666, 0.652, 0.670, 0.682, 0.706, 0.702, 0.725/                        

      data  period/   .00, .01, .06,  .09,  .10,  .20,  .30,  .40,                  
     +                .50, .60, .75, 1.00, 1.50, 2.00, 3.00, 5.00/   
C Find the requested spectral period and corresponding coefficients
      nPer = 16

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

      write (*,*) 
      write (*,*) 'Lin (2011) crustal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c2(count1),c2(count2)
     +             ,specT,c2T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c5(count1),c5(count2)
     +             ,specT,c5T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

      lnY = C1T + C2T*mag + 
     1      C3T*alog(rupdist+C4T*exp(C5T*mag))
     
      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end                                                                       

c ------------------------------------------------------------
C *** New Lin et al. (2011) Taiwan crustal footwall model SOIL (jClac=316)************
c ------------------------------------------------------------
      subroutine Lin_fw_soil(mag, rupDist, specT, period1, lnY, sigma, iflag )

      parameter (MAXPER=16)
      REAL Period(MAXPER), C1(MAXPER), C2(MAXPER), C3(MAXPER), C4(MAXPER), C5(MAXPER)
      real mag, rupDist, lnY, sigma, period1, sig(MAXPER)                                  
      real specT, c1T,  c2T,  c3T, c4T,  c5T
      integer nper, count1, count2,iflag
      
C.....MODEL COEFFICIENTS.....................
      data C1  / -3.218,-3.306,-1.896,-1.256, -1.306, -3.310, -4.880, -5.628, 
     +           -6.284,-7.252,-8.355,-9.860,-11.750,-12.827,-13.795,-14.256/

      data  C2 /  0.935,0.937,0.977,0.907,0.907,0.957,1.219,1.239, 
     +            1.311,1.429,1.536,1.692,1.919,2.025,2.036,2.120/ 

      data  C3 /  -1.464,-1.454,-1.744,-1.754,-1.734,-1.291,-1.294,-1.181, 
     +            -1.160,-1.128,-1.065,-0.995,-0.997,-0.996,-0.989,-1.144/

      data  C4 /  0.1250,0.1000,0.1400,0.1510,0.1510,0.1000,0.0310,0.0122, 
     +            0.0057,0.0025,0.0008,0.0005,0.0005,0.0005,0.0005,0.0007/  

      data  C5 /  0.650,0.670,0.720,0.720,0.710,0.700,0.910,1.020,                        
     +            1.130,1.260,1.420,1.504,1.544,1.536,1.490,1.480/                        

      data  sig / 0.630, 0.624, 0.690, 0.723, 0.726, 0.690, 0.663, 0.654,               
     +            0.651, 0.640, 0.648, 0.673, 0.713, 0.756, 0.784, 0.822/                        

      data  period/   .00, .01, .06,  .09,  .10,  .20,  .30,  .40,                  
     +                .50, .60, .75, 1.00, 1.50, 2.00, 3.00, 5.00/   
C Find the requested spectral period and corresponding coefficients
      nPer = 16

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

      write (*,*) 
      write (*,*) 'Lin (2011) crustal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c2(count1),c2(count2)
     +             ,specT,c2T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c5(count1),c5(count2)
     +             ,specT,c5T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

      lnY = C1T + C2T*mag + 
     1      C3T*alog(rupdist+C4T*exp(C5T*mag))
     
      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end                                                                       

c ------------------------------------------------------------
C *** New Lin et al. (2011) Taiwan crustal footwall model ROCK (jClac=317)************
c ------------------------------------------------------------
      subroutine Lin_fw_rock(mag, rupDist, specT, period1, lnY, sigma, iflag )

      parameter (MAXPER=16)
      REAL Period(MAXPER), C1(MAXPER), C2(MAXPER), C3(MAXPER), C4(MAXPER), C5(MAXPER)
      real mag, rupDist, lnY, sigma, period1, sig(MAXPER)                                  
      real specT, c1T,  c2T,  c3T, c4T,  c5T
      integer nper, count1, count2,iflag
      
C.....MODEL COEFFICIENTS.....................
      data C1  / -3.232,-3.193,-2.643, -2.093, -1.993, -2.659, -4.387, -5.634, 
     +           -6.391,-7.217,-8.646,-10.031,-11.633,-12.599,-13.311,-13.957/

      data  C2 /  1.047,1.017,0.937,0.907,0.907,0.960,1.169,1.328, 
     +            1.410,1.507,1.684,1.777,1.930,1.989,1.971,1.981/ 

      data  C3 /  -1.662,-1.612,-1.602,-1.642,-1.652,-1.512,-1.422,-1.399, 
     +            -1.347,-1.315,-1.304,-1.240,-1.219,-1.174,-1.140,-1.183/

      data  C4 /  0.1950,0.2100,0.2300,0.2300,0.1900,0.1482,0.0440,0.0220, 
     +            0.0180,0.0080,0.0028,0.0007,0.0005,0.0005,0.0009,0.0049/  

      data  C5 /  0.630,0.590,0.550,0.550,0.590,0.610,0.790,0.900,                        
     +            0.950,1.081,1.237,1.416,1.463,1.464,1.306,1.011/                        

      data  sig / 0.654, 0.649, 0.710, 0.756, 0.756, 0.699, 0.686, 0.682,               
     +            0.734, 0.726, 0.707, 0.717, 0.678, 0.703, 0.700, 0.725/                        

      data  period/   .00, .01, .06,  .09,  .10,  .20,  .30,  .40,                  
     +                .50, .60, .75, 1.00, 1.50, 2.00, 3.00, 5.00/   
C Find the requested spectral period and corresponding coefficients
      nPer = 16

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

      write (*,*) 
      write (*,*) 'Lin (2011) crustal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c2(count1),c2(count2)
     +             ,specT,c2T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c5(count1),c5(count2)
     +             ,specT,c5T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

      lnY = C1T + C2T*mag + 
     1      C3T*alog(rupdist+C4T*exp(C5T*mag))
     
      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end                                                                       

c ------------------------------------------------------------
C *** Lin (2009) Doctoral thesis ****************** 
c ------------------------------------------------------------
      subroutine Lin2009(mag, rupDist, specT, period1, lnY, sigma, vs,
     + iflag, ftype )

      parameter (MAXPER=32)
      REAL Period(MAXPER), C1(MAXPER), C3(MAXPER), C4(MAXPER)
      REAL C6(MAXPER), C7(MAXPER), Phi1(MAXPER), sig(MAXPER)
      real mag, rupDist, lnY, sigma  , VS ,c2 ,c5 ,c8                                 
      real ftype , F_RV, F_NM, period1                                                            
      real specT, c1T,  c3T, c4T,  c6T, c7T, phi1T
      integer nper, count1, count2,iflag
    
      
C.....MODEL COEFFICIENTS.....................
      data Period / 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 
     +0.09 , 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 
     + 0.70, 0.80, 0.90, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 
     + 4.40, 5.00 /

      data c1  / 1.0109, 1.0209, 1.0416, 1.1961, 1.3834, 1.5612, 1.6907,
     + 1.7673, 1.8689, 1.9430, 2.0218, 2.0521, 2.0333, 1.9887, 1.8827, 
     + 1.7459, 1.6821, 1.6139, 1.5288, 1.3081, 1.1383, 1.0757, 0.9935, 
     + 0.8642, 0.3150, -0.1760, -0.4103, -0.5019, -0.7206, -0.9383, 
     + -1.0405, -1.3694  /

      data c3  / 0.0000, -0.0003, 0.0017, 0.0038, 0.0087, 0.0153, 0.0210
     + , 0.0261, 0.0273, 0.0276, 0.0254, 0.0100, -0.0091, -0.0293, 
     + -0.0459, -0.0600, -0.0737, -0.0861, -0.0960, -0.1133, -0.1292, 
     + -0.1442, -0.1577, -0.1687, -0.2006, -0.2190, -0.2319, -0.2431, 
     + -0.2479, -0.2493, -0.2559, -0.2535 /

      data c4  / -1.1634, -1.1633, -1.1668, -1.2028, -1.2499, -1.2957, 
     + -1.3218, -1.3336, -1.3440, -1.3435, -1.3409, -1.2578, -1.1769, 
     + -1.1153, -1.0726, -1.0307, -1.0116, -0.9939, -0.9755, -0.9407, 
     + -0.9193, -0.9167, -0.9104, -0.9001, -0.8696, -0.8328, -0.8415, 
     + -0.8684, -0.8689, -0.8618, -0.8472, -0.8287 /

      data c6  / -0.1907, -0.1922, -0.1942, -0.1990, -0.1959, -0.1922, 
     + -0.1984, -0.2011, -0.1947, -0.2011, -0.1817, -0.1851, -0.2265, 
     + -0.2355, -0.2163, -0.1949, -0.1955, -0.2011, -0.2089, -0.2212, 
     + -0.1900, -0.1865, -0.1643, -0.1505, -0.0377, 0.0780, 0.0907, 
     +  0.1195, 0.1206, 0.1267, 0.1655, 0.2208 /

      data c7  / 0.1322, 0.1314, 0.1311, 0.1314, 0.1362, 0.1417, 0.1500, 
     + 0.1557, 0.1627, 0.1589, 0.1607, 0.1212, 0.0999, 0.0994, 0.1036, 
     + 0.1029, 0.1099, 0.1178, 0.1142, 0.1016, 0.1036, 0.1058, 0.1165, 
     + 0.1372, 0.1572, 0.1660, 0.1648, 0.1790, 0.1629, 0.1262, 0.1486, 
     + 0.1648 /

      data phi1  / -0.4741, -0.4738, -0.4700, -0.4741, -0.4806, -0.4911, 
     + -0.4900, -0.4920, -0.4944, -0.4910, -0.4825, -0.4804, -0.4350, 
     + -0.4101, -0.4361, -0.4507, -0.4734, -0.4927, -0.5035, -0.5546, 
     + -0.6037, -0.6319, -0.6577, -0.6916, -0.7582, -0.7863, -0.7939, 
     + -0.7754, -0.7673, -0.7457, -0.7042, -0.6955 / !Phi1


      data sig / 0.5363, 0.5360, 0.5349, 0.5419, 0.5506, 0.5607, 0.5718, 
     + 0.5800, 0.5887, 0.5894, 0.5911, 0.5812, 0.5715, 0.5715, 0.5768, 
     + 0.5739, 0.5696, 0.5699, 0.5681, 0.5731, 0.5703, 0.5723, 0.5736, 
     + 0.5741, 0.5743, 0.5696, 0.5521, 0.5296, 0.5256, 0.5285, 0.5239, 
     + 0.5240 /

C Find the requested spectral period and corresponding coefficients
      nPer = 32

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         phi1T   = phi1(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

      write (*,*) 
      write (*,*) 'Lin (2009) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c6(count1),c6(count2)
     +             ,specT,c6T,iflag)
      call interp (period(count1),period(count2),c7(count1),c7(count2)
     +             ,specT,c7T,iflag)
      call interp (period(count1),period(count2),phi1(count1)
     +             ,phi1(count2),specT,phi1T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

c     Set the fault mechanism term.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
         if (ftype .eq. -1) then
            F_RV = 0.0
            F_NM = 1.0
         elseif (ftype .ge. 0.5) then
            F_RV = 1.0
            F_NM = 0.0
         else
            F_RV = 0.0
            F_NM = 0.0
         endif
    
C     Compute the ground motions.
      c2 = 0.3822
      c5 = 0.1722
      c8 = 1.5184 !H

c Magnitude scaling
        if (Mag .le. 6.3) then
           r1 = c2 * (Mag-6.3) 
         else
           r1 = (-c8*c5)*(Mag-6.3)
        endif

      lnY = r1 + c1T + c3T*((8.5-Mag)**2)  
     +     + (c4T+c5*(Mag-6.3))*alog(sqrt(Rupdist**2 + (exp(c8))**2))  
     +     +  c6T*F_NM  + c7T*F_RV  + phi1T*alog(Vs/1130)
     
      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end                                                                       

c ------------------------------------------------------------
C *** TG09221 (2012) Report: TNGA ATTEN ******************
c ------------------------------------------------------------
      subroutine TG09221_2012(mag, rupDist, specT, period1, lnY, sigma, 
     +vs, iflag, ftype )

      parameter (MAXPER=106)
      REAL Period(MAXPER), C1(MAXPER), C3(MAXPER), C4(MAXPER)
      REAL C6(MAXPER), C7(MAXPER), Phi1(MAXPER), sig(MAXPER)
      real mag, rupDist, lnY, sigma  , VS ,c2 ,c5 ,c8                                 
      real ftype , F_RV, F_NM                                                            
      real specT, c1T,  c3T, c4T,  c6T, c7T, phi1T
      integer nper, count1, count2,iflag
    
      

C.....MODEL COEFFICIENTS.....................
      data period / 0.0, 0.01, 0.02, 0.022, 0.025, 0.029, 0.03, 0.032, 
     +0.035, 0.036, 0.04, 0.042, 0.044, 0.045, 0.046, 0.048, 0.05, 0.055 
     +,0.06, 0.065, 0.067, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 
     +0.11, 0.12, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 
     +0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.3, 0.32, 0.34, 0.35, 0.36, 
     +0.38, 0.4, 0.42, 0.44, 0.45, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 
     +0.667, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5
     +, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.5, 2.6, 2.8, 3, 3.2, 3.4, 3.5
     +, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5,
     + 9, 9.5, 10 /

      data c1  /  1.3979, 1.4069, 1.4159, 1.4491, 1.4955, 1.5804, 1.6016
     +, 1.6423, 1.7017, 1.7272, 1.8259, 1.8686, 1.9111, 1.9345, 1.9566,  
     +1.9999,2.0390, 2.1336, 2.2167, 2.2760, 2.2949, 2.3191, 2.3688,  
     +2.4204,2.4539, 2.4822, 2.5116, 2.5394, 2.5675, 2.5749, 2.5519,  
     +2.5388,2.5060, 2.4755, 2.4498, 2.4090, 2.3590, 2.3267, 2.3015,  
     +2.2307,2.1619, 2.1322, 2.0926, 2.0320, 2.0119, 1.9988, 1.9627,  
     +1.9283,1.9069, 1.8812, 1.8405, 1.8006, 1.7530, 1.7102, 1.6983,  
     +1.6824,1.6374, 1.5897, 1.4532, 1.3448, 1.2526, 1.2301, 1.1925,  
     +1.1277,1.0771, 1.0178, 0.9502, 0.8649, 0.7693, 0.6114, 0.4751,  
     +0.3619,0.2645, 0.1667, 0.0782, -0.0145, -0.1075, -0.2007, -0.2795, 
     +-0.4395, -0.5726, -0.6316, -0.6953, -0.8150, -0.9282, -1.0357, 
     +-1.1346, -1.1780, -1.2245, -1.3238, -1.4226, -1.5270, -1.6251, 
     +-1.7127, -1.8017, -1.8939, -2.1013, -2.3053, -2.4776, -2.6386, 
     +-2.7980, -2.9512, -3.0910, -3.2297, -3.3597, -3.4859 /

      data c3  /  0.0, -0.0006, 0.0015, 0.0016, 0.0023, 0.0032, 0.0035, 
     +0.0042, 0.0053, 0.0056, 0.0074, 0.0083, 0.0094, 0.0099, 0.0104,  
     +0.0116,0.0128, 0.0155, 0.0175, 0.0197, 0.0207, 0.0220, 0.0229,  
     +0.0223,0.0221, 0.0215, 0.0203, 0.0190, 0.0162, 0.0136, 0.0104,  
     +0.0096,0.0079, 0.0044, 0.0014, -0.0008, -0.0032, -0.0068, -0.0105, 
     +-0.0176, -0.0237, -0.0266, -0.0295, -0.0355, -0.0390, -0.0424, 
     +-0.0485, -0.0546, -0.0571, -0.0592, -0.0634, -0.0676, -0.0717, 
     +-0.0756, -0.0777, -0.0797, -0.0831, -0.0862, -0.0937, -0.1012, 
     +-0.1083, -0.1107, -0.1150, -0.1214, -0.1272, -0.1322, -0.1374, 
     +-0.1417, -0.1453, -0.1528, -0.1586, -0.1640, -0.1687, -0.1723, 
     +-0.1759, -0.1787, -0.1808, -0.1833, -0.1861, -0.1894, -0.1914, 
     +-0.1922, -0.1927, -0.1932, -0.1940, -0.1945, -0.1950, -0.1952, 
     +-0.1953, -0.1948, -0.1942, -0.1936, -0.1934, -0.1933, -0.1930, 
     +-0.1923, -0.1914, -0.1902, -0.1891, -0.1876, -0.1864, -0.1849, 
     +-0.1838, -0.1829, -0.1822, -0.1811 /

      data c4  / -1.2273,-1.2261,-1.2265, -1.2340, -1.2448, -1.2648,  
     +-1.2698,-1.2796, -1.2938, -1.2998, -1.3230, -1.3325, -1.3431,  
     +-1.3486,-1.3539, -1.3647, -1.3743, -1.3935, -1.4078, -1.4180,  
     +-1.4208,-1.4237, -1.4299, -1.4326, -1.4315, -1.4270, -1.4219,  
     +-1.4169,-1.4039, -1.3861, -1.3626, -1.3552, -1.3379, -1.3149,  
     +-1.2949,-1.2752, -1.2545, -1.2364, -1.2201, -1.1845, -1.1554,  
     +-1.1419,-1.1276, -1.1038, -1.0955, -1.0885, -1.0739, -1.0605,  
     +-1.0543,-1.0482, -1.0390, -1.0312, -1.0232, -1.0161, -1.0138,  
     +-1.0110,-1.0030, -0.9958, -0.9767, -0.9667, -0.9595, -0.9591,  
     +-0.9593,-0.9560, -0.9561, -0.9546, -0.9499, -0.9426, -0.9347,  
     +-0.9239,-0.9200, -0.9183, -0.9197, -0.9204, -0.9218, -0.9229,  
     +-0.9231,-0.9223, -0.9216, -0.9212, -0.9283, -0.9326, -0.9361,  
     +-0.9449,-0.9507, -0.9559, -0.9608, -0.9644, -0.9667, -0.9693, 
     +-0.9703,-0.9706, -0.9708, -0.9718, -0.9710, -0.9701, -0.9712,  
     +-0.9698,-0.9717, -0.9733, -0.9723, -0.9712, -0.9702, -0.9674,  
     +-0.9647,-0.9634 /

      data c6  / -0.1934,-0.1921, -0.1782, -0.1773, -0.1767, -0.1684,  
     +-0.1714,-0.1701, -0.1539, -0.1474, -0.1347, -0.1382, -0.1354,  
     +-0.1350,-0.1361, -0.1351, -0.1347, -0.1282, -0.1351, -0.1263,  
     +-0.1275,-0.1270, -0.1330, -0.1321, -0.1338, -0.1407, -0.1445,  
     +-0.1479,-0.1382, -0.1483, -0.1692, -0.1757, -0.1847, -0.2013, 
     +-0.2270,-0.2293, -0.2238, -0.2241, -0.2180, -0.2181, -0.2132,  
     +-0.2203,-0.2245, -0.2371, -0.2378, -0.2422, -0.2473, -0.2513,  
     +-0.2548,-0.2567, -0.2564, -0.2484, -0.2437, -0.2478, -0.2448,  
     +-0.2410,-0.2307, -0.2215, -0.2286, -0.2231, -0.2189, -0.2139,  
     +-0.2001,-0.1919, -0.1789, -0.1703, -0.1666, -0.1705, -0.1746,  
     +-0.1839,-0.2018, -0.2127, -0.2170, -0.2209, -0.2225, -0.2233,  
     +-0.2190,-0.2098, -0.2003, -0.2123, -0.2298, -0.2352, -0.2368,  
     +-0.2322,-0.2395, -0.2443, -0.2344, -0.2297, -0.2253, -0.2187,  
     +-0.2138,-0.2100, -0.2080, -0.2106, -0.2173, -0.2269, -0.2400,  
     +-0.2517,-0.2558, -0.2511, -0.2471, -0.2378, -0.2231, -0.2134,  
     +-0.2082,-0.1998 /

      data c7  /  0.1122,0.1117, 0.1129, 0.1128, 0.1118, 0.1126, 0.1131, 
     +0.1148, 0.1183, 0.1200, 0.1233, 0.1251, 0.1287, 0.1299, 0.1311,  
     +0.1348,0.1372, 0.1421, 0.1449, 0.1513, 0.1547, 0.1564, 0.1666,  
     +0.1759,0.1815, 0.1824, 0.1832, 0.1856, 0.1887, 0.1764, 0.1641,  
     +0.1639,0.1605, 0.1492, 0.1391, 0.1345, 0.1326, 0.1351, 0.1367,  
     +0.1304,0.1281, 0.1217, 0.1157, 0.1015, 0.0975, 0.0939, 0.0824,  
     +0.0759,0.0751, 0.0728, 0.0670, 0.0620, 0.0636, 0.0654, 0.0655,  
     +0.0663,0.0614, 0.0540, 0.0411, 0.0336, 0.0255, 0.0255, 0.0275,  
     +0.0294,0.0335, 0.0373, 0.0383, 0.0379, 0.0403, 0.0350, 0.0320,  
     +0.0320,0.0337, 0.0315, 0.0325, 0.0329, 0.0299, 0.0305, 0.0270,  
     +0.0200,0.0149, 0.0120, 0.0092, 0.0040, 0.0005, -0.0070, -0.0175,  
     +-0.0208,-0.0239, -0.0356, -0.0477, -0.0506, -0.0533, -0.0599,  
     +-0.0694,-0.0772, -0.0830, -0.0848, -0.0864, -0.0852, -0.0819,  
     +-0.0786,-0.0770, -0.0747, -0.0717, -0.0665 /

      data phi1	 /  -0.4359,-0.4344, -0.4306, -0.4312, -0.4325, -0.4351, 
     +-0.4355, -0.4354, -0.4370, -0.4375, -0.4397, -0.4408, -0.4429,  
     +-0.4445,-0.4458, -0.4476, -0.4488, -0.4462, -0.4435, -0.4435,  
     +-0.4431,-0.4440, -0.4444, -0.4440, -0.4455, -0.4415, -0.4368,  
     +-0.4324,-0.4293, -0.4309, -0.4341, -0.4339, -0.4326, -0.4301,  
     +-0.4254,-0.4201, -0.4168, -0.4149, -0.4101, -0.4000, -0.3927,  
     +-0.3907,-0.3924, -0.3971, -0.4030, -0.4067, -0.4130, -0.4200, 
     +-0.4217, -0.4246, -0.4316, -0.4405, -0.4508, -0.4593, -0.4627,  
     +-0.4645,-0.4689, -0.4766, -0.5058, -0.5398, -0.5696, -0.5789,  
     +-0.5936,-0.6099, -0.6196, -0.6296, -0.6419, -0.6554, -0.6709,  
     +-0.6998,-0.7255, -0.7423, -0.7573, -0.7676, -0.7786, -0.7875,  
     +-0.7946,-0.8052, -0.8126, -0.8251, -0.8350, -0.8357, -0.8386,  
     +-0.8458,-0.8472, -0.8491, -0.8519, -0.8527, -0.8524, -0.8505,  
     +-0.8457,-0.8431, -0.8420, -0.8407, -0.8380, -0.8367, -0.8383,  
     +-0.8371,-0.8306, -0.8189, -0.8085, -0.7991, -0.7900, -0.7823,  
     +-0.7744,-0.7684 /

      data sig / 0.6740, 0.6734, 0.6729, 0.6753, 0.6790, 0.6856, 
     +0.6872, 0.6904, 0.6953, 0.6967, 0.7044, 0.7077, 0.7112, 0.7127, 
     +0.7145, 0.7178, 0.7218, 0.7302, 0.7378, 0.7433, 0.7459, 0.7488, 
     +0.7537, 0.7573, 0.7593, 0.7606, 0.7639, 0.7658, 0.7623, 0.7585, 
     +0.7504, 0.7475, 0.7432, 0.7369, 0.7315, 0.7237, 0.7163, 0.7105, 
     +0.7048, 0.6993, 0.6973, 0.6946, 0.6934, 0.6906, 0.6901, 0.6901, 
     +0.6906, 0.6916, 0.6918, 0.6917, 0.6905, 0.6886, 0.6890, 0.6900, 
     +0.6910, 0.6912, 0.6902, 0.6889, 0.6908, 0.6930, 0.6933, 0.6944, 
     +0.6960, 0.6985, 0.7011, 0.7026, 0.7050, 0.7071, 0.7090, 0.7164, 
     +0.7237, 0.7302, 0.7331, 0.7362, 0.7406, 0.7432, 0.7456, 0.7488, 
     +0.7517, 0.7564, 0.7629, 0.7664, 0.7704, 0.7775, 0.7841, 0.7917, 
     +0.7965, 0.7988, 0.8015, 0.8064, 0.8111, 0.8168, 0.8215, 0.8248, 
     +0.8280, 0.8308, 0.8343, 0.8345, 0.8330, 0.8301, 0.8264, 0.8219, 
     +0.8173, 0.8130, 0.8087, 0.8046     /

C Find the requested spectral period and corresponding coefficients
      nPer = 106

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         phi1T   = phi1(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
c            write(*,*) 'specT = ', specT

      write (*,*) 
      write (*,*) 'TG09221 (2012) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c6(count1),c6(count2)
     +             ,specT,c6T,iflag)
      call interp (period(count1),period(count2),c7(count1),c7(count2)
     +             ,specT,c7T,iflag)
      call interp (period(count1),period(count2),phi1(count1)
     +             ,phi1(count2),specT,phi1T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

c     Set the fault mechanism term.
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
         if (ftype .eq. -1) then
            F_RV = 0.0
            F_NM = 1.0
         elseif (ftype .ge. 0.5) then
            F_RV = 1.0
            F_NM = 0.0
         else
            F_RV = 0.0
            F_NM = 0.0
         endif
    
C     Compute the ground motions.
      c2 = 0.3700
      c5 = 0.2086
      c8 = 1.4877 !H

c Magnitude scaling
        if (Mag .le. 6.3) then
           r1 = c2 * (Mag-6.3) 
        else
           r1 = (-c8*c5)*(Mag-6.3)
        endif

      lnY = r1 + c1T + c3T*(8.5-Mag)**2 + (c4T+c5*(Mag-6.3))*
     +      alog(sqrt(Rupdist**2 + (exp(c8))**2)) + c6T*F_NM + c7T*F_RV + 
     +      phi1T*alog(Vs/1130)

      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end                                                                       

c ------------------------------------------------------------
C *** NCREE (2011)  Vs30¡Ù360m/sec  ******************
c ------------------------------------------------------------
      subroutine NCREE_2011(mag, rupDist, specT, period1, lnY, sigma )

      parameter (MAXPER=30)
      REAL Period(MAXPER), C1(MAXPER), C3(MAXPER), C4(MAXPER)
      REAL sig(MAXPER)
      real mag, rupDist, lnY, sigma   ,c2(MAXPER) ,c5(MAXPER)                                  
      real specT, c1T, c2T,  c3T, c4T,  c5T
      integer nper, count1, count2,iflag
    
      
C.....MODEL COEFFICIENTS.....................
      data Period / 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2,
     + 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
     + 1, 1.5, 2, 3, 4, 5, 6, 8, 10 /

      data C1      / 0.0052, 0.0053, 0.0055, 0.0058, 0.0062, 0.0077, 
     +0.0086, 0.0095, 0.0094, 0.0092, 0.0086, 0.0080, 0.0073, 0.0068, 
     +0.0063, 0.0058, 0.0054, 0.0050, 0.0046, 0.0044, 0.0041, 0.0033, 
     +0.0023, 0.0016, 0.0010, 0.0007, 0.0006, 0.0004, 0.0003, 0.0002 /

      data C2      / 1.7295, 1.7298, 1.7299, 1.7301, 1.7303, 1.7305, 
     +1.7308, 1.7313, 1.7308, 1.7303, 1.7302, 1.7300, 1.7299, 1.7298, 
     +1.7297, 1.7296, 1.7295, 1.7294, 1.7294, 1.7293, 1.7293, 1.7292, 
     +1.7292, 1.7293, 1.7293, 1.7293, 1.7293, 1.7293, 1.7293, 1.7293 /

      data C3      / 2.0655, 2.0653, 2.0653, 2.0652, 2.0651, 2.0649, 
     +2.0648, 2.0645, 2.0648, 2.0650, 2.0651, 2.0652, 2.0653, 2.0653, 
     +2.0654, 2.0655, 2.0655, 2.0656, 2.0656, 2.0656, 2.0657, 2.0657, 
     +2.0657, 2.0655, 2.0654, 2.0654, 2.0654, 2.0654, 2.0654, 2.0654 /

      data C4      / 0.1207, 0.1203, 0.1202, 0.1200, 0.1197, 0.1194, 
     +0.1185, 0.1179, 0.1187, 0.1194, 0.1195, 0.1197, 0.1198, 0.1200, 
     +0.1201, 0.1201, 0.1202, 0.1203, 0.1203, 0.1204, 0.1204, 0.1205, 
     +0.1205, 0.1209, 0.1210, 0.1211, 0.1211, 0.1211, 0.1210, 0.1210 /

      data C5      / 0.7968, 0.7960, 0.7958, 0.7954, 0.7949, 0.7944, 
     +0.7933, 0.7923, 0.7936, 0.7947, 0.7949, 0.7951, 0.7953, 0.7954, 
     +0.7955, 0.7956, 0.7956, 0.7957, 0.7957, 0.7958, 0.7958, 0.7959, 
     +0.7959, 0.7961, 0.7962, 0.7962, 0.7962, 0.7962, 0.7962, 0.7962 /

      data sig / 0.7477, 0.7604, 0.7745, 0.7962, 0.8280, 0.8828, 
     +0.9078, 0.8531, 0.7823, 0.7342, 0.7231, 0.7267, 0.7356, 0.7384, 
     +0.7437, 0.7528, 0.7553, 0.7620, 0.7740, 0.7906, 0.8123, 0.8654, 
     +0.9011, 0.9166, 0.9241, 0.9404, 0.9681, 0.9807, 0.9708, 0.9820 /   


C Find the requested spectral period and corresponding coefficients
      nPer = 30

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1  = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         sigT    = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) 
     +         then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
c      write(*,*) 'specT = ', specT

      write (*,*) 
      write (*,*) 'NCREE (2011) Horizontal atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call interp (period(count1),period(count2),c1(count1),c1(count2)
     +             ,specT,c1T,iflag)
      call interp (period(count1),period(count2),c2(count1),c2(count2)
     +             ,specT,c2T,iflag)
      call interp (period(count1),period(count2),c3(count1),c3(count2)
     +             ,specT,c3T,iflag)
      call interp (period(count1),period(count2),c4(count1),c4(count2)
     +             ,specT,c4T,iflag)
      call interp (period(count1),period(count2),c5(count1),c5(count2)
     +             ,specT,c5T,iflag)
      call interp (period(count1),period(count2),sig(count1),
     +             sig(count2), specT,sigT,iflag)

 1011 period1 = specT                                                                                                              

      lnY = alog(c1T) + c2T*Mag + (-c3T)*alog(Rupdist + c4T*exp(c5T*mag)  
     + )

      sigma = sigT
           
c     Convert units to spectral acceleration in gal                             
      lnY = lnY + 6.89                                                          

      return                                                                    
      end               