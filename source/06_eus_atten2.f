C     Various CEUS attenuation Models

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 2-Corner, Horizontal *************
c -------------------------------------------------------------------
                                                                               
      subroutine S06_PEA2C ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                               
      real lnY, m, dist, sigT, sig, period1
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName 
      
      parameter (MAXPER=26)
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)  
                                                                                                                                                                
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1              0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1              0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 3.54103, 5.06834, 5.03119, 4.86717, 4.69293, 4.0767, 3.99907, 3.92454, 3.7751,
     1          3.624, 3.30684, 2.60454, 2.14018, 1.99361, 1.42831, 0.86777, 0.1092, -0.96872,
     1         -1.95968, -3.77355, -5.47019, -7.68301, -8.7688, -11.04809, -13.88893, -17.74463 /
      data C2 / 0.18904, 0.14806, 0.15779, 0.17018, 0.18262, 0.22547, 0.23357, 0.24169, 0.25773,
     1          0.27369, 0.30373, 0.35667, 0.39715, 0.41219, 0.46988, 0.52085, 0.59537, 0.7137,
     1          0.8081, 0.98718, 1.1259, 1.34978, 1.452, 1.64665, 1.89859, 2.22485 /
      data C4 / 2.7, 2.8, 2.8, 2.8, 2.8, 2.7, 2.7, 2.7, 2.7,
     1          2.7, 2.7, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6,
     1          2.6, 2.6, 2.5, 2.5, 2.5, 2.4, 2.3, 2.1 /
      data C6 / -2.97418, -3.08409, -3.06134, -3.03252, -3.00672, -2.9004, -2.88734, -2.87495, -2.85226,
     1          -2.83163, -2.79751, -2.69927, -2.66676, -2.6551, -2.6138, -2.58506, -2.5412, -2.465,
     1          -2.41132, -2.28113, -2.13473, -1.94573, -1.86494, -1.7001, -1.54772, -1.40084 /
      data C7 / 0.19819, 0.19935, 0.19746, 0.1956, 0.19396, 0.1872, 0.18619, 0.18521, 0.18339,
     1          0.1817, 0.17893, 0.17238, 0.16973, 0.16868, 0.16486, 0.16235, 0.15808, 0.15003,
     1          0.14449, 0.13007, 0.1171, 0.09603, 0.08722, 0.07272, 0.06068, 0.05305 /
      data C10 / -0.05814, -0.05361, -0.05377, -0.05434, -0.0552, -0.05647, -0.05717, -0.05791, -0.05952,
     1          -0.06128, -0.06512, -0.06929, -0.07573, -0.07801, -0.08671, -0.09484, -0.10506, -0.11749,
     1          -0.12529, -0.13323, -0.1383, -0.16127, -0.18125, -0.22943, -0.2896, -0.31641 /
c      data sigma / 0.84, 0.8939, 0.8902, 0.8869, 0.8795, 0.8675, 0.8604, 0.8602, 0.8521,
c     1          0.8476, 0.8468, 0.8485, 0.8339, 0.829, 0.826, 0.8272, 0.8358, 0.832,
c     1          0.8426, 0.8815, 0.8739, 0.945, 1.0095, 1.0871, 1.2228, 1.3429 /
      data sigma / 0.7353, 0.7776, 0.7823, 0.7858, 0.7817, 0.7711, 0.7644, 0.7656, 0.7585,
     1       0.7534, 0.7507, 0.7503, 0.7328, 0.7271, 0.7247, 0.7274, 0.7395, 0.7396,
     1       0.7551, 0.805, 0.8021, 0.8874, 0.9591, 1.0462, 1.1933, 1.3243 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 2 Corner atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 2-Corner'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 2-Corner-Saturation, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA2CS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none      
                                                                                
      real lnY, m, dist, sigT, sig, period1
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName
      
      parameter (MAXPER=26) 
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)                                                                                   

      data period /0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 5.91196, 7.55648, 7.51145, 7.33736, 6.61204, 6.42423, 6.34238, 6.26384, 6.10708,
     1       5.94942, 5.13706, 4.81663, 4.34277, 4.19281, 3.61568, 3.04705, 2.27626, 1.17695,
     1       0.17104, -1.6801, -3.10841, -5.75016, -6.86049, -9.24347, -12.1791, -16.16329 /
      data C2 / -0.15727, -0.20898, -0.19862, -0.18563, -0.1337, -0.11726, -0.10886, -0.10044, -0.08387,
     1       -0.06741, -0.00173, 0.02793, 0.06911, 0.08441, 0.14311, 0.19471, 0.27031, 0.39078,
     1       0.48663, 0.66971, 0.79561, 1.05061, 1.15548, 1.36201, 1.62451, 1.96535 /
      data C4 / 2.9, 3, 3, 3, 2.9, 2.9, 2.9, 2.9, 2.9,
     1       2.9, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8,
     1       2.8, 2.8, 2.8, 2.7, 2.7, 2.6, 2.5, 2.3 /
      data C6 / -3.42401, -3.55306, -3.52888, -3.49824, -3.37593, -3.34604, -3.33222, -3.31911, -3.29509,
     1       -3.27328, -3.15185, -3.12224, -3.08805, -3.07579, -3.03239, -3.00223, -2.95623, -2.87626,
     1       -2.81997, -2.68318, -2.58562, -2.32003, -2.23472, -2.05193, -1.88291, -1.71374 /
      data C7 / 0.26564, 0.26853, 0.26652, 0.26456, 0.25613, 0.25401, 0.25295, 0.25192, 0.25,
     1       0.24822, 0.23929, 0.23686, 0.23409, 0.233, 0.229, 0.22639, 0.22193, 0.21352,
     1       0.20773, 0.19261, 0.18195, 0.1554, 0.1461, 0.12954, 0.11564, 0.10547 /
      data C10 / -0.07004, -0.06551, -0.06568, -0.06625, -0.06711, -0.06838, -0.06908, -0.06982, -0.07142,
     1       -0.07318, -0.07703, -0.08119, -0.08764, -0.08991, -0.09861, -0.10675, -0.11697, -0.1294,
     1       -0.13719, -0.14513, -0.1502, -0.17317, -0.19315, -0.24133, -0.3015, -0.32832 /
c      data sigma / 0.84, 0.8939, 0.8902, 0.8869, 0.8795, 0.8675, 0.8604, 0.8602, 0.8521,
c     1       0.8476, 0.8468, 0.8485, 0.8339, 0.829, 0.826, 0.8272, 0.8358, 0.832,
c     1       0.8426, 0.8815, 0.8739, 0.945, 1.0095, 1.0871, 1.2228, 1.3429 /
      data sigma / 0.7353, 0.7776, 0.7823, 0.7858, 0.7817, 0.7711, 0.7644, 0.7656, 0.7585,
     1       0.7534, 0.7507, 0.7503, 0.7328, 0.7271, 0.7247, 0.7274, 0.7395, 0.7396,
     1       0.7551, 0.805, 0.8021, 0.8874, 0.9591, 1.0462, 1.1933, 1.3243 /
                                                                                    
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 2 Corner-Sat atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 2-Corner-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Variable-High, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CVH ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none

      real lnY, m, dist, sigT, sig, period1
      real specT, c1T, c2T, c4T, c6T, c7T, c10T     
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName 
      
      parameter (MAXPER=26)                                                                                     
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)                                                                     

      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 5.19757, 7.3541, 6.75933, 6.5773, 6.39121, 6.20019, 6.11753, 6.03843, 5.41104,
     1       5.25826, 4.94207, 4.63032, 4.17095, 4.02502, 3.46126, 2.49641, 1.72806, 0.58917,
     1       -0.51056, -2.67167, -4.46472, -7.70148, -9.12315, -11.97362, -15.20886, -18.80138 /
      data C2 / 0.07129, -0.00721, 0.03739, 0.05102, 0.06416, 0.08062, 0.08905, 0.09752, 0.14236,
     1       0.15855, 0.18877, 0.21818, 0.25933, 0.27487, 0.33544, 0.41405, 0.49782, 0.63923,
     1       0.76645, 1.03112, 1.24156, 1.62326, 1.78482, 2.06358, 2.3499, 2.59958 /
      data C4 / 2.8, 3, 2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.8,
     1       2.8, 2.8, 2.8, 2.8, 2.8, 2.8, 2.7, 2.7, 2.7,
     1       2.7, 2.7, 2.7, 2.6, 2.6, 2.5, 2.4, 2.3 /
      data C6 / -3.13247, -3.35245, -3.23117, -3.1992, -3.17073, -3.13898, -3.12417, -3.11001, -3.0004,
     1       -2.97721, -2.93847, -2.90541, -2.86633, -2.85216, -2.80186, -2.69325, -2.64087, -2.5519,
     1       -2.48971, -2.34731, -2.25138, -2.01584, -1.94059, -1.77626, -1.62679, -1.51629 /
      data C7 / 0.20485, 0.21111, 0.203, 0.20091, 0.19905, 0.19672, 0.19554, 0.19437, 0.18726,
     1       0.18528, 0.18204, 0.17924, 0.17597, 0.17467, 0.16992, 0.16254, 0.1574, 0.14804,
     1       0.14173, 0.12639, 0.11635, 0.0946, 0.08684, 0.07329, 0.0623, 0.05717 /
      data C10 / -0.07375, -0.06367, -0.06354, -0.06426, -0.06525, -0.06658, -0.06731, -0.06812, -0.06993,
     1       -0.072, -0.07672, -0.08205, -0.09083, -0.09403, -0.10713, -0.12105, -0.1416, -0.17378,
     1       -0.20549, -0.26186, -0.30377, -0.357, -0.37239, -0.38117, -0.35359, -0.25763 /
      data sigma / 0.7405, 0.7829, 0.7873, 0.7907, 0.7866, 0.7758, 0.7691, 0.7703, 0.7629,
     1       0.7578, 0.7549, 0.7543, 0.7369, 0.7311, 0.7289, 0.7317, 0.7441, 0.7445,
     1       0.7601, 0.8085, 0.8032, 0.8847, 0.9557, 1.0436, 1.1928, 1.325 /
                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Var-High atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Var-High'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Variable-Med, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CVM ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=27)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)                               

      data period / 0.00, 0.01, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 4.0393, 4.0393, 6.07941, 6.02744, 5.37895, 5.2089, 5.03867, 4.96669, 4.89845, 4.33334,
     1       4.20238, 3.92885, 3.65946, 3.26041, 3.13556, 2.27495, 1.79388, 1.12628, 0.13162,
     1       -0.84738, -2.82095, -4.51914, -7.60788, -9.00041, -11.84462, -15.15004, -19.07223 /
      data C2 / 0.10412, 0.10412, 0.03289, 0.04417, 0.08559, 0.09698, 0.11102, 0.11815, 0.12529, 0.16542,
     1       0.17878, 0.20331, 0.22693, 0.25961, 0.2722, 0.344, 0.38804, 0.45746, 0.5789,
     1       0.6896, 0.93101, 1.1322, 1.50586, 1.66899, 1.96, 2.27308, 2.57205 /
      data C4 / 2.7, 2.7, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8, 2.8, 2.7,
     1       2.7, 2.7, 2.7, 2.7, 2.7, 2.6, 2.6, 2.6, 2.6,
     1       2.6, 2.6, 2.6, 2.5, 2.5, 2.4, 2.3, 2.1 /
      data C6 / -2.97465, -2.97465, -3.18403, -3.15877, -3.04366, -3.01742, -2.98849, -2.97508, -2.9623, -2.86188,
     1       -2.84105, -2.8063, -2.7766, -2.74131, -2.72838, -2.61448, -2.58195, -2.53338, -2.45001,
     1       -2.39187, -2.25774, -2.16445, -1.94031, -1.86794, -1.70638, -1.55609, -1.41166 /
      data C7 / 0.19631, 0.19631, 0.20265, 0.20038, 0.19337, 0.19172, 0.18968, 0.18865, 0.18763, 0.1811,
     1       0.17938, 0.17658, 0.17414, 0.17129, 0.17012, 0.16182, 0.15895, 0.1542, 0.14539,
     1       0.13949, 0.12494, 0.11502, 0.09384, 0.08623, 0.07232, 0.06043, 0.05292 /
      data C10 / -0.08874, -0.08874, -0.08044, -0.08027, -0.08079, -0.0815, -0.08242, -0.08293, -0.08349, -0.08477,
     1       -0.08624, -0.08961, -0.09345, -0.09985, -0.10222, -0.11211, -0.12283, -0.1393, -0.16638,
     1       -0.19435, -0.24823, -0.29235, -0.35415, -0.37576, -0.39806, -0.38898, -0.31205 /
      data sigma / 0.7353, 0.7353, 0.7776, 0.7823, 0.7858, 0.7817, 0.7711, 0.7644, 0.7656, 0.7585,
     1       0.7534, 0.7507, 0.7503, 0.7328, 0.7271, 0.7247, 0.7274, 0.7395, 0.7396,
     1       0.7551, 0.805, 0.8021, 0.8874, 0.9591, 1.0462, 1.1933, 1.3243 /
                                                                      
C Find the requested spectral period and corresponding coefficients
      nPer = 27

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Var-Med atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Var-Med'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Variable-Low, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CVL ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)                               

      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 3.42714, 5.43661, 5.39249, 4.76103, 4.60063, 4.44314, 4.37793, 4.31663, 3.7718,
     1       3.65547, 3.41148, 3.1727, 2.81889, 2.71047, 1.91753, 1.4958, 0.91358, 0.04301,
     1       -0.82196, -2.60545, -4.48436, -7.14902, -8.74448, -11.5582, -14.90966, -18.82818 /
      data C2 / 0.10323, 0.03559, 0.04598, 0.08583, 0.09599, 0.1084, 0.11466, 0.12089, 0.15884,
     1       0.17022, 0.19064, 0.2099, 0.23601, 0.24622, 0.30871, 0.34395, 0.40083, 0.50444,
     1       0.59874, 0.81461, 1.01787, 1.36498, 1.53854, 1.83734, 2.17243, 2.50853 /
      data C4 / 2.7, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8, 2.8, 2.7,
     1       2.7, 2.7, 2.7, 2.7, 2.7, 2.6, 2.6, 2.6, 2.6,
     1       2.6, 2.6, 2.5, 2.5, 2.4, 2.3, 2.2, 2.1 /
      data C6 / -2.91721, -3.12848, -3.10501, -2.99272, -2.96783, -2.94065, -2.92815, -2.91626, -2.81861,
     1       -2.79941, -2.76748, -2.74034, -2.70814, -2.69634, -2.58706, -2.55688, -2.51171, -2.43239,
     1       -2.37729, -2.24962, -2.10529, -1.93736, -1.82313, -1.66484, -1.51375, -1.39437 /
      data C7 / 0.19218, 0.19893, 0.19687, 0.19011, 0.18862, 0.18679, 0.18586, 0.18495, 0.17867,
     1       0.17715, 0.17469, 0.17258, 0.17011, 0.16908, 0.16126, 0.15868, 0.15434, 0.14596,
     1       0.1404, 0.12651, 0.11396, 0.09559, 0.08612, 0.07203, 0.05943, 0.05155 /
      data C10 / -0.09443, -0.0876, -0.08744, -0.08783, -0.08835, -0.08901, -0.08938, -0.08978, -0.0907,
     1       -0.09175, -0.09418, -0.09693, -0.10154, -0.10327, -0.11059, -0.11864, -0.13143, -0.15344,
     1       -0.17695, -0.22591, -0.26954, -0.33624, -0.36309, -0.39809, -0.40726, -0.35284 /
      data sigma / 0.7366, 0.779, 0.7837, 0.7873, 0.7833, 0.7729, 0.7662, 0.7676, 0.7605,
     1       0.7556, 0.7531, 0.7529, 0.7357, 0.73, 0.7278, 0.7305, 0.7425, 0.7423,
     1       0.7574, 0.8073, 0.8056, 0.8932, 0.9656, 1.0524, 1.1969, 1.3261 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Var-Low atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Var-Low'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-High, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CCH ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)                               

      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 3.52033, 5.57118, 5.51533, 4.86194, 4.68655, 4.50882, 4.43278, 4.36033, 4.21631,
     1       3.64678, 3.35414, 3.06506, 2.63659, 2.50153, 1.97802, 1.08397, 0.36441, -0.70305,
     1       -1.74607, -3.81429, -5.55701, -8.68001, -10.06892, -12.85572, -16.02752, -19.70343 /
      data C2 / 0.27213, 0.19995, 0.21155, 0.25301, 0.26499, 0.2799, 0.28752, 0.29517, 0.31022,
     1       0.35047, 0.37757, 0.40391, 0.44078, 0.4548, 0.5097, 0.58122, 0.65805, 0.78942,
     1       0.90876, 1.16106, 1.36536, 1.73516, 1.8928, 2.16545, 2.44418, 2.68814 /
      data C4 / 2.7, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8,
     1       2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.6, 2.6, 2.6,
     1       2.6, 2.6, 2.6, 2.5, 2.5, 2.4, 2.3, 2.1 /
      data C6 / -2.98288, -3.18948, -3.16325, -3.04739, -3.02045, -2.99056, -2.97664, -2.96337, -2.93884,
     1       -2.84039, -2.80403, -2.77296, -2.73607, -2.72265, -2.67484, -2.57366, -2.52362, -2.4384,
     1       -2.37871, -2.24164, -2.14802, -1.92392, -1.85191, -1.69535, -1.55334, -1.41959 /
      data C7 / 0.19476, 0.20062, 0.19826, 0.19126, 0.18954, 0.18741, 0.18632, 0.18525, 0.18324,
     1       0.17689, 0.17391, 0.17133, 0.16829, 0.16707, 0.16259, 0.15569, 0.1508, 0.14185,
     1       0.13582, 0.12109, 0.11128, 0.0906, 0.08324, 0.07035, 0.06006, 0.05456 /
      data C10 / -0.05886, -0.04972, -0.04959, -0.05023, -0.05109, -0.05224, -0.05287, -0.05357, -0.05514,
     1       -0.05692, -0.06102, -0.06565, -0.07333, -0.07615, -0.08779, -0.10026, -0.11897, -0.14883,
     1       -0.17883, -0.23358, -0.27558, -0.32977, -0.34592, -0.35582, -0.32899, -0.23554 /
      data sigma / 0.714, 0.7757, 0.7803, 0.7837, 0.7795, 0.7687, 0.7619, 0.763, 0.7556,
     1       0.7505, 0.7474, 0.7468, 0.7288, 0.7228, 0.72, 0.7223, 0.734, 0.7332,
     1       0.7482, 0.7972, 0.7931, 0.8786, 0.9514, 1.0412, 1.1907, 1.3212 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-High atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Const-High'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Medium, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CCM ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)    
                           
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1        0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1        0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 3.0073, 5.0311, 4.9836, 4.34502, 4.18017, 4.01659, 3.94814, 3.88344, 3.75552,
     1        3.20588, 2.9469, 2.69216, 2.31425, 2.19706, 1.74233, 0.91433, 0.2849, -0.65379,
     1        -1.58285, -3.4733, -5.12369, -8.14308, -9.51015, -12.32672, -15.60343, -19.48096 /
      data C2 / 0.25858, 0.19, 0.20066, 0.24062, 0.25125, 0.26433, 0.27096, 0.27758, 0.29046,
     1        0.32832, 0.35069, 0.37212, 0.40161, 0.41304, 0.45792, 0.51993, 0.58358, 0.69665,
     1        0.79993, 1.02939, 1.22405, 1.58833, 1.74832, 2.03581, 2.34394, 2.63369 /
      data C4 / 2.7, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8, 2.8, 2.8,
     1        2.7, 2.7, 2.7, 2.7, 2.7, 2.7, 2.6, 2.6, 2.6,
     1        2.6, 2.6, 2.6, 2.5, 2.5, 2.4, 2.3, 2.1 /
      data C6 / -2.94208, -3.15204, -3.12766, -3.01413, -2.98855, -2.96045, -2.94748, -2.93512, -2.91238,
     1        -2.81616, -2.78273, -2.75418, -2.72023, -2.70779, -2.66338, -2.56449, -2.5173, -2.4357,
     1        -2.37885, -2.24741, -2.15471, -1.93245, -1.86136, -1.70046, -1.55118, -1.40816 /
      data C7 / 0.19152, 0.1979, 0.19576, 0.18897, 0.18741, 0.18549, 0.18451, 0.18355, 0.18175,
     1        0.17563, 0.173, 0.17072, 0.16804, 0.16694, 0.16286, 0.15619, 0.15162, 0.14303,
     1        0.13729, 0.12308, 0.11324, 0.09238, 0.08496, 0.07122, 0.0596, 0.05251 /
      data C10 / -0.05571, -0.04818, -0.04804, -0.0485, -0.04913, -0.04994, -0.05038, -0.05087, -0.05199,
     1        -0.05326, -0.05619, -0.05952, -0.06508, -0.06715, -0.07586, -0.08536, -0.10016, -0.12494,
     1        -0.1509, -0.20234, -0.24573, -0.30768, -0.33001, -0.35378, -0.3457, -0.27037 /
      data sigma / 0.7329, 0.7754, 0.7802, 0.7837, 0.7797, 0.7691, 0.7624, 0.7636, 0.7565,
     1        0.7515, 0.7487, 0.7483, 0.7308, 0.7249, 0.7224, 0.7248, 0.7365, 0.7357,
     1        0.7506, 0.8, 0.7972, 0.8842, 0.957, 1.0451, 1.1924, 1.3226 /
                                                                   
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Med atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Const-Med'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Low, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CCL ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    
    
      implicit none 
                                                                               
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER) 
                              
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 2.57877, 4.54101, 4.50286, 3.90146, 3.74605, 3.59389, 3.53102, 3.47181, 2.95159,
     1       2.83795, 2.59748, 2.36108, 2.00886, 1.90139, 1.13055, 0.71018, 0.13269, -0.72368,
     1       -1.58102, -3.33687, -4.89906, -7.79761, -9.12644, -11.90955, -15.22296, -19.09237 /
      data C2 / 0.20187, 0.13536, 0.14513, 0.18325, 0.19303, 0.205, 0.21106, 0.21711, 0.25367,
     1       0.26482, 0.28497, 0.30408, 0.33012, 0.34023, 0.40133, 0.43636, 0.49251, 0.59402,
     1       0.68709, 0.89916, 1.0839, 1.44039, 1.60015, 1.89523, 2.22463, 2.55127 /
      data C4 / 2.6, 2.8, 2.8, 2.7, 2.7, 2.7, 2.7, 2.7, 2.6,
     1       2.6, 2.6, 2.6, 2.6, 2.6, 2.5, 2.5, 2.5, 2.5,
     1       2.5, 2.5, 2.5, 2.4, 2.4, 2.3, 2.2, 2.1 /
      data C6 / -2.92333, -3.12969, -3.10735, -2.99999, -2.97578, -2.94929, -2.93706, -2.92538, -2.83154,
     1       -2.81238, -2.78005, -2.75216, -2.71828, -2.70593, -2.59747, -2.56429, -2.51516, -2.43076,
     1       -2.36933, -2.23329, -2.13693, -1.91307, -1.8433, -1.68527, -1.53646, -1.41992 /
      data C7 / 0.1963, 0.20328, 0.20134, 0.19487, 0.19341, 0.19162, 0.19071, 0.1898, 0.18373,
     1       0.18218, 0.17962, 0.17736, 0.17461, 0.17349, 0.16549, 0.16246, 0.15754, 0.1484,
     1       0.14193, 0.12682, 0.11627, 0.0945, 0.08701, 0.07303, 0.06087, 0.05359 /
      data C10 / -0.04493, -0.03854, -0.03842, -0.03882, -0.03934, -0.03999, -0.04035, -0.04075, -0.04164,
     1       -0.04266, -0.04496, -0.04754, -0.05184, -0.05345, -0.06026, -0.0677, -0.07959, -0.10021,
     1       -0.12239, -0.1695, -0.21228, -0.27892, -0.30641, -0.34299, -0.35343, -0.29994 /
      data sigma / 0.738, 0.7809, 0.7856, 0.7891, 0.785, 0.7745, 0.7678, 0.7691, 0.7619,
     1       0.757, 0.7541, 0.7536, 0.7359, 0.7301, 0.7272, 0.7292, 0.7405, 0.7391,
     1       0.7532, 0.8022, 0.8, 0.8885, 0.9618, 1.0496, 1.1946, 1.3239 /
                                                                              
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Low atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Const-Low'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-High-Sat, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CCHS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)   
                            
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 5.87466, 8.16024, 7.49863, 7.31723, 7.13246, 6.94435, 6.86357, 6.25069, 6.10266,
     1       5.95548, 5.65135, 5.35255, 4.91274, 4.77362, 3.76921, 3.241, 2.50733, 1.41567,
     1       0.35595, -1.75176, -3.52171, -6.43297, -8.18384, -11.07071, -14.33292, -17.92122 /
      data C2 / -0.06918, -0.16443, -0.11031, -0.09725, -0.08472, -0.06911, -0.06114, -0.01554, -0.00018,
     1       0.01492, 0.04288, 0.06996, 0.10766, 0.12203, 0.21066, 0.26094, 0.33905, 0.47284,
     1       0.59376, 0.85016, 1.0572, 1.42055, 1.60124, 1.88506, 2.17359, 2.41939 /
      data C4 / 2.9, 3.1, 3, 3, 3, 3, 3, 2.9, 2.9,
     1       2.9, 2.9, 2.9, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8,
     1       2.8, 2.8, 2.8, 2.8, 2.7, 2.6, 2.5, 2.4 /
      data C6 / -3.42987, -3.67439, -3.5427, -3.51057, -3.48193, -3.45017, -3.4354, -3.32741, -3.30216,
     1       -3.27916, -3.24075, -3.20795, -3.16905, -3.15491, -3.02221, -2.98696, -2.93443, -2.84495,
     1       -2.78233, -2.63826, -2.53979, -2.35414, -2.21757, -2.04376, -1.88598, -1.76932 /
      data C7 / 0.26131, 0.27088, 0.26124, 0.25921, 0.2574, 0.25514, 0.25399, 0.24639, 0.24432,
     1       0.24241, 0.23928, 0.23656, 0.23338, 0.23211, 0.22178, 0.21867, 0.21356, 0.20419,
     1       0.19789, 0.18244, 0.17215, 0.15261, 0.14122, 0.12639, 0.11439, 0.10861 /
      data C10 / -0.07032, -0.06117, -0.06105, -0.06168, -0.06255, -0.06369, -0.06433, -0.06502, -0.06659,
     1       -0.06838, -0.07247, -0.07711, -0.08479, -0.08761, -0.09924, -0.11171, -0.13042, -0.16028,
     1       -0.19028, -0.24504, -0.28704, -0.34123, -0.35737, -0.36727, -0.34044, -0.247 /
      data sigma / 0.7433, 0.7859, 0.7906, 0.7938, 0.7897, 0.779, 0.7723, 0.7734, 0.7659,
     1       0.7606, 0.7575, 0.7567, 0.7389, 0.7331, 0.7302, 0.7324, 0.7439, 0.7431,
     1       0.7578, 0.806, 0.8016, 0.8859, 0.9578, 1.0463, 1.1943, 1.3238 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-High-Sat atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Const-High-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Med-Sat, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CCMS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER) 
                              
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 5.35011, 7.60902, 6.96175, 6.79091, 6.61705, 6.44374, 6.37083, 6.30181, 5.63835,
     1       5.50813, 5.23836, 4.97456, 4.58602, 4.46498, 3.53193, 3.06862, 2.42583, 1.46377,
     1       0.51857, -1.41, -3.08744, -6.23481, -7.62375, -10.54155, -13.9107, -17.91423 /
      data C2 / -0.08193, -0.17372, -0.1209, -0.10914, -0.09799, -0.08429, -0.07735, -0.0704, -0.01981,
     1       -0.00703, 0.01612, 0.03821, 0.06845, 0.08019, 0.1588, 0.19946, 0.2643, 0.37971,
     1       0.48452, 0.71797, 0.91539, 1.29417, 1.45642, 1.75532, 2.07364, 2.37754 /
      data C4 / 2.9, 3.1, 3, 3, 3, 3, 3, 3, 2.9,
     1       2.9, 2.9, 2.9, 2.9, 2.9, 2.8, 2.8, 2.8, 2.8,
     1       2.8, 2.8, 2.8, 2.7, 2.7, 2.6, 2.5, 2.3 /
      data C6 / -3.38707, -3.63508, -3.50625, -3.47572, -3.44851, -3.41866, -3.40488, -3.39174, -3.27511,
     1       -3.25386, -3.21848, -3.1883, -3.15246, -3.13934, -3.01048, -2.9773, -2.92775, -2.84201,
     1       -2.78232, -2.6441, -2.54658, -2.30231, -2.22717, -2.04882, -1.8834, -1.71861 /
      data C7 / 0.25794, 0.26806, 0.2587, 0.25685, 0.2552, 0.25316, 0.25213, 0.25111, 0.24282,
     1       0.24112, 0.23835, 0.23596, 0.23315, 0.232, 0.22207, 0.2192, 0.21443, 0.20543,
     1       0.19943, 0.18452, 0.17419, 0.15083, 0.14299, 0.12727, 0.11388, 0.10433 /
      data C10 / -0.06717, -0.05964, -0.05949, -0.05996, -0.06059, -0.06139, -0.06184, -0.06233, -0.06344,
     1       -0.06471, -0.06764, -0.07097, -0.07654, -0.07861, -0.08732, -0.09681, -0.11161, -0.13639,
     1       -0.16236, -0.21379, -0.25719, -0.31913, -0.34147, -0.36524, -0.35716, -0.28182 /
      data sigma / 0.7427, 0.7853, 0.7901, 0.7934, 0.7894, 0.779, 0.7724, 0.7736, 0.7664,
     1       0.7612, 0.7585, 0.758, 0.7405, 0.7348, 0.7322, 0.7345, 0.7461, 0.7453,
     1       0.7599, 0.8087, 0.8057, 0.8917, 0.9637, 1.0506, 1.1963, 1.3253 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Med-Sat atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Const-Med-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Low-Sat, Horizontal *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEA1CCLS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none 
                                                                               
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)   
                            
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1        0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1        0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 4.83071, 7.01484, 6.4123, 6.25258, 6.08901, 5.92808, 5.35709, 5.29605, 5.17534,
     1        5.05581, 4.8055, 4.56076, 3.75047, 3.64126, 3.21855, 2.78943, 2.19908, 1.32039,
     1        0.44711, -1.3454, -2.93381, -5.63265, -7.30278, -9.90982, -13.58249, -17.52612 /
      data C2 / -0.12898, -0.21767, -0.16752, -0.15669, -0.14646, -0.13394, -0.09141, -0.08523, -0.07324,
     1        -0.06166, -0.04079, -0.02108, 0.03785, 0.04811, 0.08824, 0.12398, 0.1813, 0.2851,
     1        0.37973, 0.59571, 0.78321, 1.13348, 1.31553, 1.60565, 1.96072, 2.29529 /
      data C4 / 2.8, 3, 2.9, 2.9, 2.9, 2.9, 2.8, 2.8, 2.8,
     1        2.8, 2.8, 2.8, 2.7, 2.7, 2.7, 2.7, 2.7, 2.7,
     1        2.7, 2.7, 2.7, 2.7, 2.6, 2.6, 2.4, 2.3 /
      data C6 / -3.35311, -3.59567, -3.47458, -3.44589, -3.4202, -3.39211, -3.2904, -3.2784, -3.2562,
     1        -3.23598, -3.20187, -3.17249, -3.05739, -3.04474, -2.99921, -2.96449, -2.91311, -2.82477,
     1        -2.76056, -2.61803, -2.51701, -2.32961, -2.19863, -2.0728, -1.85975, -1.73008 /
      data C7 / 0.26108, 0.27165, 0.26279, 0.26109, 0.25956, 0.25767, 0.25047, 0.24954, 0.24779,
     1        0.24617, 0.24348, 0.24111, 0.23266, 0.23151, 0.22724, 0.22409, 0.21897, 0.20944,
     1        0.20271, 0.18691, 0.17588, 0.15521, 0.14381, 0.13069, 0.11403, 0.10536 /
      data C10 / -0.05639, -0.04999, -0.04987, -0.05027, -0.0508, -0.05145, -0.05181, -0.0522, -0.0531,
     1        -0.05411, -0.05641, -0.059, -0.0633, -0.0649, -0.07171, -0.07915, -0.09104, -0.11167,
     1        -0.13385, -0.18095, -0.22374, -0.29038, -0.31786, -0.35444, -0.36489, -0.31139 /
      data sigma / 0.7477, 0.7908, 0.7956, 0.7989, 0.7949, 0.7844, 0.7779, 0.7791, 0.7718,
     1        0.7668, 0.7638, 0.7632, 0.7456, 0.7398, 0.7369, 0.7389, 0.7499, 0.7485,
     1        0.7625, 0.8109, 0.8086, 0.8961, 0.9688, 1.0556, 1.199, 1.3269 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Low-Sat atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. (2002) 1-Corner-Const-Low-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 2-Corner, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG2C ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    
 
      implicit none
                                                                               
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER) 
                              
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1        0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1        0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 9.90148, 14.33926, 14.14384, 13.81046, 12.42738, 12.02874, 11.85646, 11.69837, 10.51625,
     1        10.26367, 9.7961, 9.36118, 7.97786, 7.79199, 7.09082, 5.74732, 4.8783, 2.9958,
     1        1.93674, -0.12511, -2.15831, -5.10947, -6.28665, -8.83853, -12.20468, -16.41379 /
      data C2 / -0.12757, -0.28175, -0.25488, -0.22653, -0.15043, -0.12544, -0.11325, -0.10103, -0.03055,
     1        -0.00857, 0.03112, 0.06628, 0.15425, 0.16914, 0.2271, 0.31977, 0.39924, 0.55002,
     1        0.64984, 0.8388, 0.99778, 1.24004, 1.33826, 1.54958, 1.83553, 2.20767 /
      data C4 / 3.4, 3.7, 3.7, 3.7, 3.6, 3.6, 3.6, 3.6, 3.5,
     1        3.5, 3.5, 3.5, 3.4, 3.4, 3.4, 3.3, 3.3, 3.2,
     1        3.2, 3.2, 3.1, 3, 3, 2.9, 2.7, 2.5 /
      data C6 / -4.30771, -4.93249, -4.87587, -4.80279, -4.56505, -4.49489, -4.46406, -4.43516, -4.23161,
     1        -4.18673, -4.11197, -4.0525, -3.85005, -3.8308, -3.76365, -3.58708, -3.51506, -3.30621,
     1        -3.22152, -3.0503, -2.82918, -2.53021, -2.43987, -2.21747, -1.95111, -1.74567 /
      data C7 / 0.25806, 0.27571, 0.26995, 0.26447, 0.25176, 0.24783, 0.24596, 0.2441, 0.23289,
     1        0.22992, 0.22526, 0.22188, 0.21146, 0.21046, 0.20689, 0.19746, 0.19281, 0.17993,
     1        0.17312, 0.15823, 0.14191, 0.11775, 0.10973, 0.09356, 0.07737, 0.06829 /
      data C10 / -0.05882, -0.05331, -0.05318, -0.0537, -0.05445, -0.05578, -0.0565, -0.05724, -0.05885,
     1        -0.06062, -0.06459, -0.06891, -0.07567, -0.07794, -0.0868, -0.09504, -0.10565, -0.11838,
     1        -0.12608, -0.13454, -0.13999, -0.16323, -0.18303, -0.23217, -0.29415, -0.33131 /
      data sigma / 0.9031, 0.9685, 0.9722, 0.9821, 0.9842, 0.9799, 0.9733, 0.973, 0.9664,
     1        0.9628, 0.9596, 0.9587, 0.9369, 0.9306, 0.9177, 0.9111, 0.9109, 0.8991,
     1        0.904, 0.9349, 0.9271, 0.9956, 1.0612, 1.1393, 1.2665, 1.3791 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 2 Corner atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 2-Corner'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 2-Corner-Saturation, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG2CS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)   
                            
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 13.52127, 17.40998, 17.2009, 16.84962, 16.45698, 16.02362, 15.83611, 14.53975, 14.23988,
     1       13.96635, 13.46434, 12.01553, 11.39716, 11.20299, 9.5919, 8.93501, 8.03762, 5.96206,
     1       4.87153, 2.74561, 0.53026, -2.0663, -3.8314, -6.54516, -9.76826, -14.54243 /
      data C2 / -0.58872, -0.72596, -0.69766, -0.66794, -0.64342, -0.61652, -0.60342, -0.52162, -0.49678,
     1       -0.47347, -0.43179, -0.33472, -0.28656, -0.27131, -0.15735, -0.10061, -0.01963, 0.14984,
     1       0.25189, 0.44592, 0.62351, 0.85356, 0.9871, 1.21492, 1.50496, 1.91703 /
      data C4 / 3.6, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7,
     1       3.7, 3.7, 3.6, 3.6, 3.6, 3.5, 3.5, 3.5, 3.4,
     1       3.4, 3.4, 3.3, 3.3, 3.2, 3.1, 3, 2.7 /
      data C6 / -4.95888, -5.48554, -5.42663, -5.35052, -5.27833, -5.20228, -5.16886, -4.95186, -4.89689,
     1       -4.84844, -4.76782, -4.53955, -4.46693, -4.44627, -4.22645, -4.1674, -4.09059, -3.85124,
     1       -3.76122, -3.57911, -3.32873, -3.0922, -2.90058, -2.65122, -2.41187, -2.10859 /
      data C7 / 0.34391, 0.3583, 0.35229, 0.34658, 0.34234, 0.33809, 0.33606, 0.32294, 0.31927,
     1       0.31607, 0.31108, 0.29762, 0.29383, 0.29278, 0.28002, 0.27663, 0.27175, 0.25592,
     1       0.24875, 0.23302, 0.2137, 0.19167, 0.17769, 0.15876, 0.14194, 0.12608 /
      data C10 / -0.07073, -0.06521, -0.06508, -0.06561, -0.06636, -0.06769, -0.0684, -0.06914, -0.07075,
     1       -0.07252, -0.07649, -0.08082, -0.08758, -0.08984, -0.09871, -0.10695, -0.11756, -0.13028,
     1       -0.13799, -0.14644, -0.1519, -0.17514, -0.19493, -0.24408, -0.30605, -0.34322 /
      data sigma / 0.9031, 0.9685, 0.9722, 0.9821, 0.9842, 0.9799, 0.9733, 0.973, 0.9664,
     1       0.9628, 0.9596, 0.9587, 0.9369, 0.9306, 0.9177, 0.9111, 0.9109, 0.8991,
     1       0.904, 0.9349, 0.9271, 0.9956, 1.0612, 1.1393, 1.2665, 1.3791 /
                                                                                    
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 2 Corner-Sat atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 2-Corner-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Variable-High, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CVH ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)  
                             
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 17.56501, 22.54184, 22.40623, 22.02865, 21.61053, 21.11875, 20.91108, 20.72257, 20.37544,
     1       20.04863, 19.43335, 18.86865, 16.81254, 16.59083, 15.77669, 15.04806, 12.94874, 11.55566,
     1       10.33931, 6.90409, 4.94057, 0.72588, -1.56307, -5.06249, -9.02136, -13.91863 /
      data C2 / -0.73081, -0.94055, -0.91637, -0.88044, -0.8526, -0.81789, -0.80174, -0.78585, -0.75402,
     1       -0.72294, -0.66605, -0.61723, -0.48116, -0.46325, -0.39607, -0.33259, -0.17886, -0.03985,
     1       0.08773, 0.41108, 0.63072, 1.0561, 1.26035, 1.59204, 1.95631, 2.33294 /
      data C4 / 3.7, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9,
     1       3.9, 3.9, 3.9, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7,
     1       3.7, 3.6, 3.6, 3.5, 3.4, 3.3, 3.2, 3 /
      data C6 / -5.65962, -6.41953, -6.37768, -6.29329, -6.21225, -6.11772, -6.07691, -6.039, -5.96788,
     1       -5.9014, -5.78203, -5.68075, -5.34165, -5.30738, -5.1878, -5.08734, -4.77368, -4.60821,
     1       -4.47467, -4.06719, -3.88811, -3.47121, -3.23534, -2.93546, -2.64856, -2.33817 /
      data C7 / 0.36566, 0.39634, 0.39115, 0.38353, 0.37806, 0.37124, 0.36809, 0.36501, 0.35891,
     1       0.35305, 0.34275, 0.33451, 0.31279, 0.31025, 0.30141, 0.29385, 0.27355, 0.2601,
     1       0.24889, 0.21878, 0.20282, 0.17081, 0.15492, 0.13367, 0.11437, 0.09908 /
      data C10 / -0.07661, -0.06885, -0.06821, -0.06773, -0.06793, -0.06833, -0.06854, -0.06877, -0.06938,
     1       -0.07023, -0.07256, -0.07559, -0.08098, -0.08296, -0.09155, -0.10095, -0.11589, -0.14124,
     1       -0.16562, -0.21837, -0.25891, -0.32004, -0.34181, -0.365, -0.3607, -0.29154 /
      data sigma / 0.91, 0.9758, 0.9794, 0.9892, 0.9912, 0.9866, 0.98, 0.9795, 0.9726,
     1       0.9687, 0.9649, 0.9635, 0.9412, 0.9348, 0.9217, 0.9152, 0.915, 0.9039,
     1       0.9093, 0.9403, 0.9321, 0.9978, 1.0621, 1.139, 1.2674, 1.3821 /
                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Var-High atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Var-High'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Variable-Med, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CVM ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)    
                           
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 15.27441, 19.88182, 19.78142, 20.8987, 20.50874, 18.64419, 18.46167, 18.29779, 17.99875,
     1       17.71756, 17.18425, 16.69027, 16.0265, 15.82366, 13.93331, 13.29372, 12.43075, 10.18655,
     1       9.10831, 6.01523, 4.21778, 0.25617, -1.92096, -5.3148, -9.25169, -14.54986 /
      data C2 / -0.61726, -0.80582, -0.78593, -0.84219, -0.81779, -0.703, -0.68951, -0.67629, -0.64981,
     1       -0.62387, -0.57629, -0.53555, -0.48449, -0.4692, -0.34606, -0.2936, -0.21755, -0.04192,
     1       0.0681, 0.35576, 0.55654, 0.96349, 1.16422, 1.49937, 1.88136, 2.30998 /
      data C4 / 3.6, 3.8, 3.8, 3.9, 3.9, 3.8, 3.8, 3.8, 3.8,
     1       3.8, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7, 3.7, 3.6,
     1       3.6, 3.5, 3.5, 3.4, 3.3, 3.2, 3.1, 2.8 /
      data C6 / -5.30301, -6.01166, -5.97732, -6.14187, -6.06641, -5.74551, -5.70923, -5.67576, -5.61318,
     1       -5.55455, -5.44841, -5.35752, -5.24545, -5.21303, -4.9083, -4.81633, -4.69883, -4.37417,
     1       -4.24952, -3.87179, -3.70387, -3.30857, -3.08551, -2.79932, -2.51971, -2.15716 /
      data C7 / 0.34239, 0.37073, 0.36646, 0.37424, 0.36944, 0.34952, 0.34685, 0.34425, 0.33907,
     1       0.33408, 0.32521, 0.31804, 0.30976, 0.30744, 0.28841, 0.28167, 0.27269, 0.25066,
     1       0.24024, 0.21222, 0.19724, 0.16627, 0.15079, 0.12984, 0.1101, 0.09152 /
      data C10 / -0.09155, -0.08525, -0.08469, -0.08421, -0.0843, -0.08451, -0.08461, -0.08473, -0.08507,
     1       -0.08558, -0.08711, -0.08918, -0.09294, -0.09434, -0.10055, -0.10754, -0.119, -0.13932,
     1       -0.15976, -0.20693, -0.24619, -0.31189, -0.33869, -0.37439, -0.38752, -0.34105 /
      data sigma / 0.9031, 0.9685, 0.9722, 0.9821, 0.9842, 0.9799, 0.9733, 0.973, 0.9664,
     1       0.9628, 0.9596, 0.9587, 0.9369, 0.9306, 0.9177, 0.9111, 0.9109, 0.8991,
     1       0.904, 0.9349, 0.9271, 0.9956, 1.0612, 1.1393, 1.2665, 1.3791 /
                                                                      
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Var-Med atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Var-Med'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Variable-Low, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CVL ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)  
                             
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 14.35825, 18.89574, 18.81796, 19.93414, 19.56869, 19.13672, 18.95886, 17.43775, 17.16559,
     1       16.90986, 16.42234, 15.96772, 15.35608, 15.16957, 13.35601, 12.77568, 11.99662, 9.88541,
     1       8.91163, 6.03088, 4.3827, 0.60771, -1.52175, -4.86624, -9.24105, -14.25163 /
      data C2 / -0.58678, -0.76779, -0.75053, -0.8076, -0.78618, -0.75882, -0.74635, -0.65344, -0.63027,
     1       -0.60751, -0.56572, -0.53011, -0.48604, -0.47298, -0.35932, -0.31521, -0.25104, -0.09326,
     1       0.00198, 0.26122, 0.44249, 0.82916, 1.02753, 1.36569, 1.78174, 2.2322 /
      data C4 / 3.6, 3.8, 3.8, 3.9, 3.9, 3.9, 3.9, 3.8, 3.8,
     1       3.8, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7, 3.7, 3.6,
     1       3.6, 3.5, 3.5, 3.4, 3.3, 3.2, 3, 2.8 /
      data C6 / -5.18268, -5.88521, -5.85612, -6.0219, -5.95141, -5.86837, -5.83325, -5.57562, -5.51821,
     1       -5.46428, -5.36596, -5.28098, -5.17562, -5.14508, -4.8492, -4.76255, -4.65152, -4.33601,
     1       -4.21624, -3.84799, -3.68677, -3.29808, -3.07493, -2.78892, -2.4321, -2.1362 /
      data C7 / 0.33157, 0.35942, 0.35579, 0.36387, 0.35968, 0.35426, 0.3518, 0.33613, 0.33158,
     1       0.32715, 0.31923, 0.31278, 0.30531, 0.30322, 0.28512, 0.27903, 0.27084, 0.24975,
     1       0.23987, 0.21274, 0.19846, 0.16792, 0.15219, 0.13082, 0.10715, 0.09013 /
      data C10 / -0.09714, -0.09203, -0.09155, -0.09109, -0.09113, -0.09122, -0.09125, -0.0913, -0.09148,
     1       -0.09179, -0.09279, -0.09419, -0.0968, -0.09778, -0.10218, -0.10726, -0.11582, -0.13159,
     1       -0.14814, -0.1886, -0.22467, -0.29094, -0.32102, -0.36709, -0.39711, -0.37542 /
      data sigma / 0.9003, 0.9651, 0.9688, 0.9788, 0.9811, 0.977, 0.9705, 0.9704, 0.9642,
     1       0.9609, 0.9583, 0.9578, 0.9367, 0.9305, 0.9182, 0.9119, 0.9119, 0.9001,
     1       0.9045, 0.9346, 0.9268, 0.9968, 1.0631, 1.1426, 1.269, 1.379 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Var-Low atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Var-Low'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-High, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CCH ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)  
                             
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 16.08302, 19.51657, 19.40352, 20.51285, 20.10895, 19.63311, 19.43335, 19.25288, 18.92158,
     1       17.25027, 16.69021, 16.17236, 15.47703, 15.26405, 13.33561, 12.66006, 11.74743, 9.4332,
     1       8.29484, 5.08922, 3.21915, -0.80232, -2.37232, -6.35091, -10.18213, -14.92343 /
      data C2 / -0.53416, -0.64891, -0.62778, -0.68238, -0.6565, -0.624, -0.60894, -0.59416, -0.56453,
     1       -0.45606, -0.40538, -0.36176, -0.30669, -0.29006, -0.16266, -0.10531, -0.02246, 0.162,
     1       0.28003, 0.58119, 0.78973, 1.19788, 1.36543, 1.71518, 2.06757, 2.43027 /
      data C4 / 3.7, 3.8, 3.8, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9,
     1       3.8, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7, 3.7, 3.6,
     1       3.6, 3.5, 3.5, 3.4, 3.4, 3.2, 3.1, 2.9 /
      data C6 / -5.54089, -6.0488, -6.01142, -6.17374, -6.09544, -6.00383, -5.96446, -5.92801, -5.85972,
     1       -5.57097, -5.46052, -5.36632, -5.25065, -5.2172, -4.90916, -4.81457, -4.69411, -4.36652,
     1       -4.23968, -3.85783, -3.68683, -3.29152, -3.17502, -2.78676, -2.51603, -2.22382 /
      data C7 / 0.35572, 0.37128, 0.3667, 0.37408, 0.36898, 0.36256, 0.35961, 0.35673, 0.351,
     1       0.33242, 0.3231, 0.31558, 0.30691, 0.30448, 0.28523, 0.2782, 0.26892, 0.24675,
     1       0.23616, 0.20794, 0.19276, 0.16227, 0.15221, 0.12722, 0.10917, 0.09485 /
      data C10 / -0.0616, -0.05461, -0.05403, -0.05357, -0.05373, -0.05406, -0.05422, -0.05441, -0.05491,
     1       -0.05561, -0.05758, -0.06017, -0.06481, -0.06652, -0.07402, -0.08232, -0.09568, -0.11874,
     1       -0.14131, -0.19127, -0.23069, -0.29176, -0.3141, -0.33853, -0.33532, -0.2685 /
      data sigma / 0.9054, 0.9715, 0.9752, 0.9849, 0.987, 0.9824, 0.9757, 0.9753, 0.9685,
     1       0.9646, 0.9608, 0.9592, 0.9369, 0.9303, 0.917, 0.9098, 0.909, 0.8966,
     1       0.901, 0.9311, 0.9225, 0.9905, 1.0564, 1.1358, 1.2655, 1.3793 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
       
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-High atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Const-High'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Medium, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CCM ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    
 
      implicit none
                                                                               
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)  
                             
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 14.09083, 18.66451, 18.5759, 19.69186, 19.31416, 18.86815, 18.68342, 17.14437, 16.8585,
     1       16.58985, 16.07898, 15.60406, 14.96522, 14.76997, 12.91716, 12.30467, 11.47972, 9.29703,
     1       8.2654, 5.26938, 3.53746, -0.33936, -2.48893, -5.83777, -9.7274, -14.97295 /
      data C2 / -0.44176, -0.62536, -0.60707, -0.66304, -0.64033, -0.61146, -0.59824, -0.50458, -0.47991,
     1       -0.45572, -0.41128, -0.37326, -0.32579, -0.3116, -0.19398, -0.14553, -0.07511, 0.09136,
     1       0.19447, 0.46791, 0.65933, 1.05439, 1.25139, 1.58067, 1.95572, 2.37485 /
      data C4 / 3.6, 3.8, 3.8, 3.9, 3.9, 3.9, 3.9, 3.8, 3.8,
     1       3.8, 3.8, 3.8, 3.8, 3.8, 3.7, 3.7, 3.7, 3.6,
     1       3.6, 3.5, 3.5, 3.4, 3.3, 3.2, 3.1, 2.8 /
      data C6 / -5.23954, -5.94544, -5.91373, -6.07864, -6.00557, -5.91969, -5.88319, -5.62237, -5.56225,
     1       -5.50586, -5.40344, -5.31533, -5.20638, -5.17479, -4.87438, -4.78472, -4.66998, -4.34955,
     1       -4.22696, -3.85384, -3.68869, -3.29591, -3.07316, -2.78845, -2.51098, -2.14967 /
      data C7 / 0.33332, 0.36127, 0.35737, 0.36517, 0.3607, 0.35496, 0.35235, 0.33654, 0.33169,
     1       0.327, 0.31863, 0.31184, 0.30397, 0.30176, 0.28333, 0.2769, 0.26828, 0.24683,
     1       0.23667, 0.20926, 0.19464, 0.16405, 0.14866, 0.12797, 0.10865, 0.09041 /
      data C10 / -0.05839, -0.05272, -0.05222, -0.05176, -0.05184, -0.05201, -0.05209, -0.05218, -0.05245,
     1       -0.05288, -0.05418, -0.05594, -0.05917, -0.06037, -0.06574, -0.07185, -0.08198, -0.10023,
     1       -0.11889, -0.1629, -0.20041, -0.26496, -0.29199, -0.32881, -0.34332, -0.29863 /
      data sigma / 0.9, 0.9655, 0.9692, 0.9791, 0.9814, 0.9771, 0.9706, 0.9704, 0.9639,
     1       0.9604, 0.9574, 0.9565, 0.9348, 0.9285, 0.9156, 0.909, 0.9085, 0.8962,
     1       0.9004, 0.9302, 0.9219, 0.9913, 1.0576, 1.1373, 1.2658, 1.3782 /
                                                                   
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Med atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Const-Med'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Low, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CCL ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                               
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)  
                             
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 13.95924, 17.18549, 18.4581, 18.16727, 17.82567, 17.4211, 17.25532, 17.1086, 16.84332,
     1       16.59372, 16.11475, 14.4586, 13.87411, 13.69482, 13.03902, 11.41205, 10.65378, 8.63053,
     1       7.67899, 5.67074, 3.27417, -0.39847, -1.89095, -5.22266, -9.56051, -14.51574 /
      data C2 / -0.51587, -0.6188, -0.68363, -0.65766, -0.6382, -0.61316, -0.60178, -0.59071, -0.56847,
     1       -0.54655, -0.50595, -0.40049, -0.35857, -0.34602, -0.29926, -0.19463, -0.13196, 0.02055,
     1       0.11381, 0.32399, 0.54493, 0.92377, 1.08997, 1.4263, 1.83677, 2.27743 /
      data C4 / 3.6, 3.7, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8,
     1       3.8, 3.8, 3.7, 3.7, 3.7, 3.7, 3.6, 3.6, 3.5,
     1       3.5, 3.5, 3.4, 3.3, 3.3, 3.2, 3, 2.8 /
      data C6 / -5.26211, -5.74801, -5.94442, -5.87782, -5.81156, -5.73333, -5.70035, -5.67026, -5.61422,
     1       -5.56147, -5.46456, -5.17926, -5.07762, -5.04794, -4.9434, -4.67925, -4.56857, -4.26127,
     1       -4.13972, -3.90708, -3.60867, -3.22226, -3.10319, -2.81135, -2.45407, -2.16004 /
      data C7 / 0.33959, 0.3548, 0.36483, 0.35931, 0.35548, 0.35047, 0.34821, 0.34603, 0.34167,
     1       0.3374, 0.32967, 0.31164, 0.3044, 0.30235, 0.29508, 0.27864, 0.27022, 0.24915,
     1       0.23879, 0.21814, 0.19579, 0.16458, 0.15363, 0.1316, 0.10808, 0.09156 /
      data C10 / -0.04706, -0.04232, -0.04192, -0.04152, -0.04159, -0.04171, -0.04176, -0.04182, -0.04202,
     1       -0.04233, -0.04331, -0.04464, -0.04708, -0.04799, -0.05208, -0.05678, -0.06472, -0.0794,
     1       -0.09489, -0.13325, -0.16807, -0.23333, -0.26352, -0.31057, -0.34227, -0.3221 /
      data sigma / 0.8999, 0.965, 0.9686, 0.9786, 0.981, 0.977, 0.9705, 0.9705, 0.9642,
     1       0.9611, 0.9586, 0.9582, 0.9368, 0.9306, 0.9181, 0.9114, 0.9108, 0.8982,
     1       0.9018, 0.9306, 0.9217, 0.9915, 1.0581, 1.1387, 1.2666, 1.3771 /
                                                                              
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Low atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf (2002) 1-Corner-Const-Low'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-High-Sat, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CCHS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)   
                            
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 19.2806, 24.64848, 24.51845, 24.13701, 23.71273, 23.21284, 23.00275, 22.8127, 22.46355,
     1       22.13508, 21.51625, 20.94827, 20.19141, 18.50771, 17.69717, 16.97448, 14.74399, 13.37181,
     1       12.17467, 8.6523, 6.70599, 2.38331, 0.76119, -2.92833, -7.63385, -12.65135 /
      data C2 / -0.99348, -1.24267, -1.21926, -1.18374, -1.15652, -1.12231, -1.10647, -1.09092, -1.05978,
     1       -1.02936, -0.97369, -0.92609, -0.86646, -0.7544, -0.69006, -0.62936, -0.46111, -0.32798,
     1       -0.20531, 0.12586, 0.34088, 0.77804, 0.94998, 1.29933, 1.70983, 2.09962 /
      data C4 / 3.8, 4, 4, 4, 4, 4, 4, 4, 4,
     1       4, 4, 4, 4, 3.9, 3.9, 3.9, 3.8, 3.8,
     1       3.8, 3.7, 3.7, 3.6, 3.6, 3.5, 3.3, 3.1 /
      data C6 / -6.11546, -6.93799, -6.89791, -6.81358, -6.73191, -6.63632, -6.59524, -6.55721, -6.48596,
     1       -6.41929, -6.29904, -6.19648, -6.07056, -5.79682, -5.67505, -5.57259, -5.23503, -5.06579,
     1       -4.92914, -4.49751, -4.31367, -3.87013, -3.74466, -3.40943, -2.99118, -2.65372 /
      data C7 / 0.44085, 0.47804, 0.47308, 0.46562, 0.4603, 0.45359, 0.45051, 0.4475, 0.44153,
     1       0.43578, 0.42562, 0.41744, 0.40801, 0.39016, 0.38135, 0.37377, 0.35058, 0.33699,
     1       0.32564, 0.29269, 0.27643, 0.24129, 0.23048, 0.2058, 0.17807, 0.15931 /
      data C10 / -0.07305, -0.06606, -0.06549, -0.06503, -0.06519, -0.06551, -0.06568, -0.06586, -0.06636,
     1       	-0.06706, -0.06904, -0.07162, -0.07627, -0.07798, -0.08547, -0.09378, -0.10713, -0.13019,
     1       -0.15276, -0.20273, -0.24215, -0.30322, -0.32555, -0.34999, -0.34677, -0.27995 /
      data sigma / 0.9129, 0.9786, 0.9822, 0.9922, 0.9941, 0.9896, 0.9831, 0.9827, 0.976,
     1       0.9721, 0.9683, 0.9669, 0.9448, 0.9383, 0.925, 0.9181, 0.9173, 0.905,
     1       0.9094, 0.9392, 0.9306, 0.998, 1.0634, 1.1422, 1.271, 1.3837 /

C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-High-Sat atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf (2002) 1-Corner-Const-High-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Med-Sat, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CCMS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)                               

      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 18.4435, 23.74361, 23.64102, 23.29188, 22.89529, 22.42674, 22.23245, 22.05839, 21.74098,
     1       21.44261, 20.87722, 20.35514, 18.20699, 18.00368, 17.26232, 16.60473, 14.47052, 13.22694,
     1       12.13819, 8.82898, 7.023, 2.84583, 1.25366, -3.08593, -7.18334, -12.36761 /
      data C2 / -0.97271, -1.21411, -1.19389, -1.16219, -1.13831, -1.10792, -1.09401, -1.08042, -1.05317,
     1       -1.02649, -0.97756, -0.9359, -0.7901, -0.77536, -0.7206, -0.66904, -0.51361, -0.39862,
     1       -0.29103, 0.01214, 0.20984, 0.63403, 0.80656, 1.20331, 1.59845, 2.03153 /
      data C4 / 3.8, 4, 4, 4, 4, 4, 4, 4, 4,
     1       4, 4, 4, 3.9, 3.9, 3.9, 3.9, 3.8, 3.8,
     1       3.8, 3.7, 3.7, 3.6, 3.6, 3.4, 3.3, 3.1 /
      data C6 / -6.00839, -6.82593, -6.79201, -6.71452, -6.63833, -6.54872, -6.51064, -6.47563, -6.41027,
     1       -6.34895, -6.23744, -6.14146, -5.78572, -5.75279, -5.63758, -5.54036, -5.20994, -5.04733,
     1       -4.91519, -4.49285, -4.3152, -3.87433, -3.74757, -3.29709, -2.98527, -2.63877 /
      data C7 / 0.43044, 0.46721, 0.46299, 0.45635, 0.45168, 0.44569, 0.44297, 0.44033, 0.43505,
     1       0.42994, 0.42083, 0.41343, 0.38965, 0.38735, 0.37933, 0.37238, 0.34992, 0.33707,
     1       0.32617, 0.29408, 0.27841, 0.24315, 0.23195, 0.2001, 0.17746, 0.15707 /
      data C10 / -0.06984, -0.06417, -0.06368, -0.06322, -0.0633, -0.06346, -0.06354, -0.06363, -0.06391,
     1       -0.06434, -0.06563, -0.0674, -0.07062, -0.07183, -0.0772, -0.08331, -0.09344, -0.11168,
     1       -0.13035, -0.17435, -0.21186, -0.27642, -0.30344, -0.34026, -0.35477, -0.31008 /
      data sigma / 0.9073, 0.9726, 0.9761, 0.9861, 0.9883, 0.9841, 0.9777, 0.9775, 0.9711,
     1       0.9676, 0.9646, 0.9639, 0.9425, 0.9362, 0.9235, 0.9169, 0.9165, 0.9044,
     1       0.9086, 0.9382, 0.93, 0.9987, 1.0648, 1.1441, 1.2715, 1.3828 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Med-Sat atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Const-Med-Sat'                      
                                                                                
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end                                                                       

c -------------------------------------------------------------------           
C *** Silva et al. (2002) 1-Corner-Constant-Low-Sat, Horizontal Gulf *************
c -------------------------------------------------------------------
                                                                                
      subroutine S06_PEAG1CCLS ( m, dist, lnY, sig, specT,                      
     1                  attenName, period1,iflag )                                    

      implicit none
                                                                                
      real lnY, m, dist, sigT, sig, period1                              
      real specT, c1T, c2T, c4T, c6T, c7T, c10T
      integer nper, count1, count2, iflag, MAXPER, i
      character*80 attenName                                                    
                                                                              
      parameter (MAXPER=26)                                                      
      real c1(MAXPER), c2(MAXPER), c4(MAXPER), sigma(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c10(MAXPER), period(MAXPER)
                               
      data period / 0.00, 0.02, 0.025, 0.0323, 0.04, 0.05, 0.055, 0.06, 0.07,
     1       0.08, 0.1, 0.12, 0.15, 0.16, 0.2, 0.24, 0.3, 0.4,
     1       0.5, 0.75, 1, 1.6, 2, 3.000, 5, 10 /
      data C1 / 16.96725, 21.95228, 21.87723, 21.57041, 21.21243, 20.78819, 20.61413, 20.45986, 20.18058,
     1       19.91774, 19.41416, 18.94247, 18.30605, 18.11166, 16.08791, 15.48361, 14.67303, 12.35996,
     1       11.3548, 9.24442, 6.58784, 2.63452, 1.0916, -2.46635, -7.12762, -12.34486 /
      data C2 / -0.95374, -1.17915, -1.16223, -1.13489, -1.11449, -1.08818, -1.07622, -1.0646, -1.04126,
     1       -1.01826, -0.9757, -0.93912, -0.89356, -0.87999, -0.74305, -0.69729, -0.63086, -0.45077,
     1       -0.35323, -0.13451, 0.1121, 0.51855, 0.68925, 1.04807, 1.49148, 1.95845 /
      data C4 / 3.7, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9,
     1       3.9, 3.9, 3.9, 3.9, 3.9, 3.8, 3.8, 3.8, 3.7,
     1       3.7, 3.7, 3.6, 3.5, 3.5, 3.4, 3.2, 3 /
      data C6 / -5.807, -6.58099, -6.55296, -6.48374, -6.41476, -6.33324, -6.29889, -6.26753, -6.20917,
     1       -6.15421, -6.05319, -5.96502, -5.85465, -5.82243, -5.49275, -5.40014, -5.28067, -4.9283,
     1       -4.79775, -4.54798, -4.20864, -3.777, -3.64922, -3.32053, -2.91057, -2.57317 /
      data C7 / 0.42132, 0.45634, 0.45285, 0.4471, 0.44311, 0.43789, 0.43554, 0.43327, 0.42872,
     1       0.42427, 0.41621, 0.40954, 0.4017, 0.39947, 0.37755, 0.3709, 0.36186, 0.33648,
     1       0.32541, 0.30336, 0.27695, 0.24128, 0.22956, 0.20386, 0.17493, 0.15404 /
      data C10 / -0.05852, -0.05377, -0.05337, -0.05297, -0.05305, -0.05317, -0.05322, -0.05328, -0.05347,
     1       -0.05379, -0.05476, -0.0561, -0.05854, -0.05945, -0.06353, -0.06824, -0.07617, -0.09086,
     1       -0.10635, -0.14471, -0.17953, -0.24479, -0.27498, -0.32203, -0.35373, -0.33355 /
      data sigma / 0.9075, 0.9723, 0.9759, 0.9861, 0.9884, 0.9844, 0.9781, 0.978, 0.9719,
     1       0.9687, 0.9662, 0.9657, 0.9446, 0.9385, 0.9261, 0.9195, 0.919, 0.9065,
     1       0.9101, 0.9387, 0.93, 0.9992, 1.0655, 1.1457, 1.2728, 1.3821 /
                                                                                
C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c4T     = c4(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c10T     = c10(1)
         sigT     = sigma(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Silva et al. (2002) 1 Corner-Const-Low-Sat atttenuation Gulf model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),sigma(count1),sigma(count2),
     +             specT,sigT,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'Silva et al. Gulf(2002) 1-Corner-Const-Low-Sat'                      
                                                                               
      lnY = c1T + c2T*M + (c6T+c7T*M)*alog(dist+exp(c4T)) + c10t*(M-6)*(M-6)                                                       
                                                                                
c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      sig = sigT                                            
                                                                                
      return                                                                    
      end            
c ----------------------------

c ------------------------------------------------------------------           
C *** Campbell (2003) Horizontal for Hybrid CEUS ****
c ------------------------------------------------------------------

      subroutine S06_CHY03 ( mag, rupdist, lnY, sigma, specT, period, iflag )                                                

c    Campbell (2003) Horizontal Hybrid CEUS
C    Hard Rock Site Conditions

      implicit none

      real mag, rupDist, lnY, sigma, period                                    
      real R, term1, term2, term3                                                           
      real c1(18), c2(18), c3(18), c4(18), c5(18), c6(18)
      real c7(18), c8(18), c9(18), c10(18), c11(18)
      real c12(18), c13(18)
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T
      real c8T, c9T, c10T, c11T, c12T, c13T
      integer nper, count1, count2, iflag, i
      real period1(18), r1, r2

      data period1 / 0.00, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20,   
     &               0.30, 0.50, 0.75, 1.00, 1.50,  2.00, 3.00, 4.00, 5.00 /
      data c1  / 0.0305,  0.0305,  1.3535,  1.186,  0.3736, -0.0395, -0.1475,      
     &          -0.1901, -0.4328, -0.6906, -0.5907, -0.5429, -0.6104,   
     &          -0.9666, -1.4306, -2.2331, -2.7975, -3.2365 /
      data c2  / 0.633,  0.633, 0.63, 0.622, 0.616, 0.615, 0.613, 0.616, 0.617, 
     &           0.609, 0.534, 0.48, 0.451, 0.441, 0.459, 0.492, 0.507, 0.525  /
      data c3  / -0.0427, -0.0427, -0.0404, -0.0362, -0.0353, -0.0353, -0.0353,  
     &           -0.0478, -0.0586, -0.0786, -0.1379, -0.1806, -0.209,   
     &           -0.2405, -0.2552, -0.2646, -0.2738, -0.2792  /
      data c4  / -1.591, -1.591, -1.787, -1.691, -1.469, -1.383, -1.369, -1.368, 
     &           -1.32, -1.28, -1.216, -1.184, -1.158, -1.135, -1.124,  
     &           -1.121, -1.119, -1.117  /
      data c5  / -0.00428, -0.00428, -0.00388, -0.00367, -0.00378, -0.00421, -0.00454, 
     &           -0.00473, -0.0046, -0.00414, -0.00341, -0.00288, -0.00255, 
     &           -0.00213, -0.00187, -0.00154, -0.00135, -0.00117  /
      data c6  / 0.000483, 0.000483, 0.000497, 0.000501, 0.0005, 0.000486, 0.00046, 
     &           0.000393, 0.000337, 0.000263, 0.000194, 0.00016, 0.000141, 
     &           0.000119, 0.000103, 0.000084, 0.000074, 0.000064  /
      data c7  / 0.683, 0.683, 1.02, 0.922, 0.63, 0.491, 0.484, 0.461, 0.399, 0.349, 
     &           0.318, 0.304, 0.299, 0.304, 0.31, 0.31, 0.294, 0.293  /
      data c8  / 0.416, 0.416, 0.363, 0.376, 0.423, 0.463, 0.467, 0.478, 0.493, 0.502, 
     &           0.503, 0.504, 0.503, 0.5, 0.499, 0.499, 0.506, 0.507  /
      data c9  / 1.14, 1.14, 0.851, 0.759, 0.771, 0.955, 1.096, 1.239, 1.25, 1.241, 
     &           1.166, 1.11, 1.067, 1.029, 1.015, 1.014, 1.018, 1.018 /
      data c10 / -0.873, -0.873, -0.715, -0.922, -1.239, -1.349, -1.284, -1.079, -0.928, 
     &           -0.753, -0.606, -0.526, -0.482, -0.438, -0.417, -0.393, -0.386, -0.374 /
      data c11 / 1.03, 1.03, 1.03, 1.03, 1.042, 1.052, 1.059, 1.068, 1.077, 1.081, 
     &           1.098, 1.105, 1.11, 1.099, 1.093, 1.09, 1.092, 1.092  /
      data c12 / -0.086, -0.086, -0.086, -0.086, -0.0838, -0.0838, -0.0838, -0.0838, 
     &           -0.0838, -0.0838, -0.0824, -0.0806, -0.0793, -0.0771, -0.0758, 
     &           -0.0737, -0.0722, -0.0722 /
      data c13 / 0.414, 0.414, 0.414, 0.414, 0.443, 0.453, 0.46, 0.469, 0.478, 
     &           0.482, 0.508, 0.528, 0.543, 0.547, 0.551, 0.562, 0.575, 0.575 /

C Find the requested spectral period and corresponding coefficients
      nPer = 18

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         c3T      = c3(1)
         c4T      = c4(1)
         c5T      = c5(1)
         c6T      = c6(1)
         c7T      = c7(1)
         c8T      = c8(1)
         c9T      = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Campbell (2003) Hor. Hybrid-CEUS atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period1(count1),period1(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period1(count1),period1(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period1(count1),period1(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period1(count1),period1(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period1(count1),period1(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period1(count1),period1(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period1(count1),period1(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period1(count1),period1(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period1(count1),period1(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period1(count1),period1(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period1(count1),period1(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period1(count1),period1(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)

 1011 period = specT                                                                                                              

C     Compute the ground motions.
      r1 = 70.0
      r2 = 130.0
      R = sqrt(rupdist*rupdist + (c7T*exp(c8T*mag))**2.0)
      term1 = c2T*mag + c3T*(8.5-mag)**2.0
      term2 = c4T*alog(R) + (c5T + c6T*mag)*rupdist

      if (rupdist .le. r1) then
         term3 = 0.0
      elseif (rupdist .le. r2) then
         term3 = c9T*(alog(rupdist) - alog(r1) )
      else
         term3 = c9T*(alog(rupdist) - alog(r1) ) + c10T*(alog(rupdist) - alog(r2) )
      endif

      lnY = c1T + term1 + term2 + term3

C     Now compute the sigma
      if (mag .lt. 7.16) then
         sigma = c11T + c12T*mag
      else
         sigma = c13T
      endif

C     Now convert to Ln Units in gals.
      lnY = lnY + 6.89

      return                                                                    
      end                                                                                                                                                   

c ------------------------------------------------------------------           
C *** Campbell (2003) Epistemic Model from B. Youngs ****
c ------------------------------------------------------------------

      subroutine S06_CHY03Eps ( mag, rupdist, sigmaeps, specT, period, iflag )                                                

c    Campbell (2003) Horizontal Hybrid CEUS
C    Hard Rock Site Conditions
C    Functional fit to Epistemic Values

      implicit none

      real mag, rupDist, sigmaeps, period                                    
      real EU1(19), EU2(19), EU3(19), EU4(19) 
      real EU5(19), EU6(19), EU7(19), EU8(19) 
      real EU1T, EU2T, EU3T, EU4T, EU5T, EU6T, EU7T, EU8T
      real period1(19), specT
      integer nper, count1, count2, iflag, i

      data Period1 / 0.00, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1,
     1               0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 
     2               3.0, 4.0 /
      data EU1 / -0.033125, 0.51389, 0.4732, 0.43996, 0.35982, 0.29766, 
     1            0.14199, 0.070928,0.090445,0.094299,0.228,0.071481,
     2           -0.049923,-0.024928,-0.31915,0.30874,0.38988,0.030017,
     3            0.17177 /
      data EU2 / 0.03861, -0.012873, -0.023817, -0.032759, -0.022174,
     1          -0.013963, 0.0035871, 0.010224, 0.0016472, -0.0018457,
     2          -0.024652, -0.0067245, 0.0071813, -0.00096697, 0.059038,
     3          -0.045457, -0.048256, 0.018689, 0.0040366 /
      data EU3 / 0.0087249, 0.0024367, 0.0025327, 0.0026111, 0.0038935,
     1           0.0048883, 0.007836, 0.0092802, 0.0070123, 0.0058854,
     2          -0.00036376, 0.0027743, 0.0052083, 0.003623, 0.0072105,
     3          -0.0018214, -0.0026968, -0.00083408, -0.0029759 /
      data EU4 / -0.0011583, -0.0003466, -0.0001655, -0.000017529,
     1           -0.00022502, -0.00038596, -0.00082995, -0.0010257,
     2           -0.00062561, -0.00036639, 0.00058705, 0.00017773,
     3           -0.00013977, 0.00010171, -0.00068873, 0.00074708,
     4            0.00079911,0.00020333,0.00040235 /
      data EU5 / 0.48433, 0.077638, -0.0016085, -0.066357, -0.084339,
     1          -0.098287, -0.21126, -0.18045, -0.25327, -0.19499,
     2           0.0026986, 0.22268, 0.39332, 0.41276, 0.80731,
     3           0.26528, 0.22295, 0.83934, 1.3955 /
      data EU6 / -0.054313, -0.00050347, 0.018919, 0.034789, 0.032453,
     1            0.030642, 0.055054, 0.054054, 0.072855, 0.066763,
     2            0.028926, -0.0035868, -0.028806, -0.029292, -0.10475,
     3           -0.010409, -0.012206, -0.11457, -0.19531 /
      data EU7 / 0.10412, -0.28007, -0.15434, -0.051614, 0.0030899,
     1           0.045522, 0.089641, -0.012585, -0.16112, -0.22373,
     2          -0.31465, -0.21393, -0.13581, -0.13727, 0.069083,
     3          -0.1091, -0.07587, 0.24009, 0.53063 /
      data EU8 / -0.029092, 0.032305, 0.016081, 0.0028249, -0.0091361,
     1           -0.018414, -0.032056, -0.018555, 0.0042663, 0.013275,
     2            0.030403, 0.020432, 0.012697, 0.013937, -0.027262,
     3            0.01131, 0.0059025, -0.057348, -0.11271 /

C Find the requested spectral period and corresponding coefficients
      nPer = 19

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         EU1T      = EU1(1)
         EU2T      = EU2(1)
         EU3T      = EU3(1)
         EU4T      = EU4(1)
         EU5T      = EU5(1)
         EU6T      = EU6(1)
         EU7T      = EU7(1)
         EU8T      = EU8(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Campbell (2003) Hor. Hybrid-CEUS Epistemic model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),EU1(count1),EU1(count2),
     +             specT,EU1T,iflag)
      call S24_interp (period1(count1),period1(count2),EU2(count1),EU2(count2),
     +             specT,EU2T,iflag)
      call S24_interp (period1(count1),period1(count2),EU3(count1),EU3(count2),
     +             specT,EU3T,iflag)
      call S24_interp (period1(count1),period1(count2),EU4(count1),EU4(count2),
     +             specT,EU4T,iflag)
      call S24_interp (period1(count1),period1(count2),EU5(count1),EU5(count2),
     +             specT,EU5T,iflag)
      call S24_interp (period1(count1),period1(count2),EU6(count1),EU6(count2),
     +             specT,EU6T,iflag)
      call S24_interp (period1(count1),period1(count2),EU7(count1),EU7(count2),
     +             specT,EU7T,iflag)
      call S24_interp (period1(count1),period1(count2),EU8(count1),EU8(count2),
     +             specT,EU8T,iflag)

 1011 period = specT                                                                                                              

C     Compute the epistemic sigma value.
      sigmaeps = EU1T + EU2T*mag + (EU3T + EU4T*mag)*rupDist + 
     1          (EU5T + EU6T*mag)*exp( (EU7T + EU8T*mag)*rupDist )

      return                                                                    
      end                                                                                                                                                   

c ------------------------------------------------------------------           
C *** Atkinson and Boore (2006) Horizontal for CEUS Hard Rock ****
c ------------------------------------------------------------------

      subroutine S06_AB06 ( mag, rupdist, lnY, sigma, specT, period, iflag )                                                

c    Atkinson and Boore (2006) Horizontal CEUS Hard Rock
C    Hard Rock Site Conditions

      implicit none

      real mag, rupDist, lnY, sigma, period                                    
      real term1, term2                                                           
      real c1(25), c2(25), c3(25), c4(25), c5(25), c6(25)
      real c7(25), c8(25), c9(25), c10(25)
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T
      real c8T, c9T, c10T, r0, r1, r2
      integer nper, count1, count2, iflag, i
      real period1(25), term0

      data Period1 / 0.00, 0.025, 0.031, 0.04, 0.05, 0.063, 0.079,      
     &               0.1, 0.125, 0.158, 0.199, 0.251, 0.315, 0.397,     
     &               0.5, 0.629, 0.794, 1.00, 1.25, 1.59, 2.00, 2.5,    
     &               3.13, 4.00, 5.00 / 
      data c1 / 0.907, 1.52, 1.44, 1.26, 1.11, 0.911, 0.691, 0.48, 0.214, 
     &         -0.146, -0.615, -1.12, -1.72, -2.44, -3.22, -3.92, -4.6, 
     &         -5.27, -5.72, -6.04, -6.18, -6.17, -6.04, -5.79, -5.41 /
      data c2 / 0.983, 0.96, 0.959, 0.968, 0.972, 0.98, 0.997, 1.02, 1.05, 
     &          1.12, 1.23, 1.34, 1.48, 1.65, 1.83, 1.99, 2.13, 2.26, 2.32, 
     &          2.34, 2.3, 2.21, 2.08, 1.92, 1.71 /
      data c3 / -0.066, -0.0635, -0.0628, -0.0623, -0.062, -0.0621, -0.0628, 
     &          -0.064, -0.0666, -0.0714, -0.0789, -0.0872, -0.0974, -0.108, 
     &          -0.12, -0.131, -0.141, -0.148, -0.151, -0.15, -0.144, -0.135,
     &          -0.122, -0.107, -0.0901 /
      data c4 / -2.7, -2.81, -2.71, -2.58, -2.47, -2.36, -2.26, -2.2, -2.15, 
     &          -2.12, -2.09, -2.08, -2.08, -2.05, -2.02, -2.05, -2.06, -2.07, 
     &          -2.1, -2.16, -2.22, -2.3, -2.37, -2.44, -2.54 /
      data c5 / 0.159, 0.146, 0.14, 0.132, 0.128, 0.126, 0.125, 0.127, 0.13, 0.13, 
     &          0.131, 0.135, 0.138, 0.136, 0.134, 0.142, 0.147, 0.15, 0.157, 0.166, 
     &          0.177, 0.19, 0.2, 0.211, 0.227 /
      data c6 / -2.8, -3.65, -3.73, -3.64, -3.39, -2.97, -2.49, -2.01, -1.61, -1.3, 
     &          -1.12, -0.971, -0.889, -0.843, -0.813, -0.782, -0.797, -0.813, -0.82, 
     &          -0.87, -0.937, -0.986, -1.07, -1.16, -1.27 /
      data c7 / 0.212, 0.236, 0.234, 0.228, 0.214, 0.191, 0.164, 0.133, 0.105, 0.0831, 
     &          0.0679, 0.0563, 0.0487, 0.0448, 0.0444, 0.043, 0.0435, 0.0467, 0.0519, 
     &          0.0605, 0.0707, 0.0786, 0.0895, 0.102, 0.116 /
      data c8 / -0.301, -0.654, -0.543, -0.351, -0.139, 0.107, 0.214, 0.337, 0.427, 
     &           0.562, 0.606, 0.614, 0.61, 0.739, 0.884, 0.788, 0.775, 0.826, 0.856, 
     &           0.921, 0.952, 0.968, 1, 1.01, 0.979 /
      data c9 / -0.0653, -0.055, -0.0645, -0.0813, -0.0984, -0.117, -0.121, -0.127, 
     &          -0.13, -0.144, -0.146, -0.143, -0.139, -0.156, -0.175, -0.159, -0.156, 
     &          -0.162, -0.166, -0.173, -0.177, -0.177, -0.18, -0.182, -0.177 /
      data c10 / -0.000448, -0.0000485, -0.0000323, -0.000123, -0.000317, -0.000579, 
     &           -0.000847, -0.00105, -0.00115, -0.00118, -0.00113, -0.00106, 
     &           -0.000954, -0.000851, -0.00077, -0.000695, -0.000579, -0.000386, 
     &           -0.000433, -0.000375, -0.000322, -0.000282, -0.000231, -0.000201, 
     &           -0.000176 /
 
C Find the requested spectral period and corresponding coefficients
      nPer = 25

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         c3T      = c3(1)
         c4T      = c4(1)
         c5T      = c5(1)
         c6T      = c6(1)
         c7T      = c7(1)
         c8T      = c8(1)
         c9T      = c9(1)
         c10T     = c10(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Atkinson and Boore (2006) Hor. CEUS Hard Rock atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period1(count1),period1(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period1(count1),period1(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period1(count1),period1(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period1(count1),period1(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period1(count1),period1(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period1(count1),period1(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period1(count1),period1(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period1(count1),period1(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period1(count1),period1(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)

 1011 period = specT                                                                                                              

C     Compute the ground motions.
      r0 = 10.0
      r1 = 70.0 
      r2 = 140.0
      term0 = max(alog10(r0/rupdist),0.0)
      term1 = min(alog10(rupdist),alog10(r1))
      term2 = max(alog10(rupdist/r2),0.0)

      lnY = c1T + c2T*mag + c3T*mag*mag + (c4T+c5T*mag)*term1 + 
     1        (c6T + c7T*mag)*term2 + (c8T + c9T*mag)*term0 + c10T*rupdist

C     Now convert to Ln Units.
      lnY = alog(10.0)*lnY
c      lnY = lnY + 6.89
      sigma = alog(10.0)*0.30

      return                                                                    
      end                                                                                                                                                   


c ------------------------------------------------------------------           
C *** Atkinson and Boore (2006) Horizontal for CEUS Hard Rock ****
C        **** Alternative Stress Drop Scale Factors ****
c ------------------------------------------------------------------

      subroutine S06_AB06SF2 ( mag, specT, period, SF2 )                                                

c     Atkinson and Boore (2006) Horizontal CEUS Hard Rock
C     Hard Rock Site Conditions
C     Alternative Stress Drop values (see Eq. 6 and Table 7)

      implicit none

      real mag, period, SF2, specT, term1, term2
      real period1(26)
      real Mh(26), Ml(26), Delta, MhT, MlT
      integer nper, count1, count2, iFlag, i

      data Period1 / 0.00, 0.01, 0.025, 0.031, 0.04, 0.05, 0.063, 0.079,      
     &               0.1, 0.125, 0.158, 0.199, 0.251, 0.315, 0.397,     
     &               0.5, 0.629, 0.794, 1.00, 1.25, 1.59, 2.00, 2.5,    
     &               3.13, 4.00, 5.00 / 
      data Mh / 5.5, 5.5, 5.0, 5.0, 5.0, 5.0, 5.17, 5.34, 5.5, 5.67, 5.84, 
     1          6.0, 6.12, 6.25, 6.37, 6.5, 6.7, 6.95, 7.2, 7.45, 7.7,
     2          8.0 ,8.12, 8.25 ,8.37 ,8.5 /

      data Ml / 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.17, 0.34, 0.5, 1.15, 1.85, 
     1          2.5, 2.9, 3.3, 3.65, 4.0, 4.17, 4.34, 4.5, 4.67, 4.84, 
     2          5.0, 5.25, 5.5, 5.75, 6.0 /

C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         MlT      = Ml(1)
         MhT      = Mh(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Atkinson and Boore (2006) Hor. CEUS Hard Rock atttenuation model'
      write (*,*) 'variable stress drop factors'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),Ml(count1),Ml(count2),
     +             specT,MlT,iflag)
      call S24_interp (period1(count1),period1(count2),Mh(count1),Mh(count2),
     +             specT,MhT,iflag)

 1011 period = specT                                                                                                              

C     Compute the ground motions.
      delta = 0.15
      term1 = max( (mag - mlT), 0.0)
      term2 = 0.05 + delta*(term1/(MhT-MlT))  
      SF2 = min(delta+0.05, term2)

      return                                                                    
      end                                                                                                                                                   

c ------------------------------------------------------------------           
C *** Atkinson and Boore (2006) Horizontal for CEUS Vs=760m/sec ****
c ------------------------------------------------------------------

      subroutine S06_AB06vs760 ( mag, rupdist, lnY, sigma, specT, period, iflag )                                                

c    Atkinson and Boore (2006) Horizontal CEUS Vs=760m/sec

      implicit none

      real mag, rupDist, lnY, sigma, period                                    
      real term1, term2                                                         
      real c1(26), c2(26), c3(26), c4(26), c5(26), c6(26)
      real c7(26), c8(26), c9(26), c10(26)
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T
      real c8T, c9T, c10T
      integer nper, count1, count2, iflag, i
      real period1(26), r0, r1, r2, term0

      data Period1 / 0.00, 0.0100, 0.025, 0.031, 0.04, 0.05, 0.063, 0.079,      
     &               0.1, 0.125, 0.158, 0.199, 0.251, 0.315, 0.397,     
     &               0.5, 0.629, 0.794, 1.00, 1.25, 1.59, 2.00, 2.5,    
     &               3.13, 4.00, 5.00 / 
      data c1 / 0.523, 0.523, 1.05, 1.19, 1.26, 1.21, 1.11, 0.967, 0.782, 0.536, 
     &          0.119, -0.306, -0.876, -1.56, -2.28, -3.01, -3.75, -4.45, 
     &         -5.06, -5.49, -5.75, -5.85, -5.8, -5.59, -5.26, -4.85 /
      data c2 / 0.969, 0.969, 0.903, 0.888, 0.879, 0.883, 0.888, 0.903, 0.924, 
     &          0.965, 1.06, 1.16, 1.29, 1.46, 1.63, 1.8, 1.97, 2.12, 2.23,
     &          2.29, 2.29, 2.23, 2.13, 1.97, 1.79, 1.58 /
      data c3 / -0.062, -0.062, -0.0577, -0.0564, -0.0552, -0.0544, -0.0539,    
     &          -0.0548, -0.0556, -0.0584, -0.0647, -0.0721, -0.0819,   
     &          -0.0931, -0.105, -0.118, -0.129, -0.139, -0.145, -0.148, 
     &          -0.145, -0.139, -0.128, -0.114, -0.0979, -0.0807 /
      data c4 / -2.44, -2.44, -2.57, -2.58, -2.54, -2.44, -2.33, -2.25, -2.17, 
     &          -2.11, -2.05, -2.04, -2.01, -1.98, -1.97, -1.98, -2.00, 
     &          -2.01, -2.03, -2.08, -2.13, -2.2, -2.26, -2.33, -2.44, -2.53 /
      data c5 / 0.147, 0.147, 0.148, 0.145, 0.139, 0.13, 0.123, 0.122, 0.119,  
     &          0.121, 0.119, 0.122, 0.123, 0.121, 0.123, 0.127, 0.131, 
     &          0.136, 0.141, 0.15, 0.158, 0.169, 0.179, 0.191, 0.207, 0.222 /
      data c6 / -2.34, -2.34, -2.65, -2.84, -2.99, -3.04, -2.88, -2.53, -2.1,  
     &          -1.67, -1.36, -1.15, -1.03, -0.947, -0.888, -0.847, -0.842, 
     &          -0.858, -0.874, -0.9, -0.957, -1.04, -1.12, -1.2, -1.31, -1.43 /
      data c7 / 0.191, 0.191, 0.207, 0.212, 0.216, 0.213, 0.201, 0.178, 0.148, 
     &          0.116, 0.0916, 0.0738, 0.0634, 0.0558, 0.0503, 0.047, 0.0482,
     &          0.0498, 0.0541, 0.0579, 0.0676, 0.08, 0.0954, 0.11, 0.121, 0.136 /
      data c8 / -0.087, -0.087, -0.408, -0.437, -0.391, -0.21, -0.0319, 0.1, 0.285, 
     &           0.343, 0.516, 0.508, 0.581, 0.65, 0.684, 0.667, 0.677, 0.708, 
     &           0.792, 0.821, 0.867, 0.867, 0.891, 0.845, 0.734, 0.634 /
      data c9 / -0.0829, -0.0829, -0.0577, -0.0587, -0.0675, -0.09, -0.107, -0.115, 
     &          -0.132, -0.132, -0.15, -0.143, -0.149, -0.156, -0.158,  
     &          -0.155, -0.156, -0.159, -0.17, -0.172, -0.179, -0.179,  
     &          -0.18, -0.172, -0.156, -0.141 /
      data c10 / -0.00063, -0.00063, -0.000512, -0.000433, -0.000388, -0.000415,  
     &           -0.000548, -0.000772, -0.00099, -0.00113, -0.00118,    
     &           -0.00114, -0.00105, -0.000955, -0.000859, -0.000768,   
     &           -0.000676, -0.000575, -0.000489, -0.000407, -0.000343, 
     &           -0.000286, -0.00026, -0.000245, -0.000196, -0.000161 /

C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         c3T      = c3(1)
         c4T      = c4(1)
         c5T      = c5(1)
         c6T      = c6(1)
         c7T      = c7(1)
         c8T      = c8(1)
         c9T      = c9(1)
         c10T     = c10(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'Atkinson and Boore (2006) Hor. CEUS Vs=760m/s atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period1(count1),period1(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period1(count1),period1(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period1(count1),period1(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period1(count1),period1(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period1(count1),period1(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period1(count1),period1(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period1(count1),period1(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period1(count1),period1(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period1(count1),period1(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)

 1011 period = specT                                                                                                              

C     Compute the ground motions.
      r0 = 10.0
      r1 = 70.0 
      r2 = 140.0
      term0 = max(alog10(r0/rupdist),0.0)
      term1 = min(alog10(rupdist),alog10(r1))
      term2 = max(alog10(rupdist/r2),0.0)

      lnY = c1T + c2T*mag + c3T*mag*mag + (c4T+c5T*mag)*term1 + 
     1        (c6T + c7T*mag)*term2 + (c8T + c9T*mag)*term0 + c10T*rupdist

C     Now convert to Ln Units.
      lnY = alog(10.0)*lnY
c      lnY = lnY + 6.89
      sigma = alog(10.0)*0.30

      return                                                                    
      end                                                                                                                                                   

c ------------------------------------------------------------------           
C *** Atkinson (2008) Horizontal for CEUS-NGA BA08 Vs=760m/sec ****
c ------------------------------------------------------------------

      subroutine S06_A08vs760 ( mag, jbdist, specT, BA08lnY, period, LnY, sigma, iflag )                                                

c    Atkinson (2008) Horizontal CEUS-NGA BA08 Vs=760m/sec

      implicit none

      real mag, jbDist, lnY, sigma, period, BA08LnY                                    
      real c1(7), c2(7), c0(7)
      real specT, c1T, c2T, c0T
      integer nper, count1, count2, iflag, i
      real per1(7), factor

      data Per1 / 0.00, 0.10, 0.20, 0.50, 1.00, 2.00, 5.00 / 
      data c1 / 0.00120, 0.00124, 0.00144, 0.00113, 0.000556, 0.000520, -0.00107 /
      data c2 / 0.00000230, 0.00000199, 0.00000127, 0.000000698, 
     1          0.000000744, 0.000000376, 0.00000149 /
      data c0 / 0.163, 0.093, -0.155, -0.356, -0.404, -0.379, -0.319 /

C Find the requested spectral period and corresponding coefficients
      nPer = 7
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = per1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         c0T      = c0(1)
         goto 1011
      elseif (specT .ne. 0.0) then

C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. per1(i) .and. specT .le. per1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
     
      write (*,*) 
      write (*,*) 'Atkinson (2008) Hor. CEUS-NGA BA08 Vs760 atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (per1(count1),per1(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (per1(count1),per1(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (per1(count1),per1(count2),c0(count1),c0(count2),
     +             specT,c0T,iflag)

 1011 period = specT                                                                                                              

c     Compute the ground motions adjustment factor.
      factor = c0T + c1T*jbdist + c2T*jbdist*jbdist
      factor = factor*alog(10.0)

      lnY = BA08LnY + factor

      return                                                                    
      end                                                                                                                                                   
     

c ------------------------------------------------------------------           
C *** Atkinson (2008) Horizontal for CEUS-NGA BA08 Vs=760m/sec   ****
C      *** Alternative C0 values                                 ****
c ------------------------------------------------------------------

      subroutine S06_A08vs760C0 ( mag, jbdist, specT, BA08lnY, period, LnY, sigma, iflag )                                                

c    Atkinson (2008) Horizontal CEUS-NGA BA08 Vs=760m/sec
C     Alternative C0 values

      implicit none

      real mag, jbDist, lnY, sigma, period, BA08LnY                                    
      real c1(7), c2(7), c0(7)
      real specT, c1T, c2T, c0T
      integer nper, count1, count2, iflag, i
      real period1(7), factor

      data Period1 / 0.00, 0.10, 0.20, 0.50, 1.00, 2.00, 5.00 / 
      data c1 / 0.00120, 0.00124, 0.00144, 0.00113, 0.000556, 0.000520, -0.00107 /
      data c2 / 0.00000230, 0.00000199, 0.00000127, 0.000000698, 
     1          0.000000744, 0.000000376, 0.00000149 /
      data c0 / 0.287, 0.143, -0.102, -0.364, -0.376, -0.419, -0.271 /

C Find the requested spectral period and corresponding coefficients
      nPer = 7

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         c1T      = c1(1)
         c2T      = c2(1)
         c0T      = c0(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
     
      write (*,*) 
      write (*,*) 'Atkinson (2008) Hor. CEUS-NGA BA08 Vs760 atttenuation model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period1(count1),period1(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period1(count1),period1(count2),c0(count1),c0(count2),
     +             specT,c0T,iflag)

 1011 period = specT                                                                                                              

c     Compute the ground motions adjustment factor.
      factor = c0T + c1T*jbdist + c2T*jbdist*jbdist
      factor = factor*alog(10.0)

      lnY = BA08LnY + factor

      return                                                                    
      end                                                                                                                                                   
     
c ------------------------------------------------------------------           
C *** NGA Sigma (Average from NGA Models for M5.5) ****
c ------------------------------------------------------------------

      subroutine S06_NGASigma ( specT, period, sigma )                                                

C     Average sigma from NGA models for M5.5.

      implicit none

      real period, sigma, specT
      real period1(21)
      real sig(21), sigT
      integer nper, count1, count2, iFlag, i

      data Period1 / 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 
     1               0.2, 0.25, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0, 3.0, 
     2               4.0, 5.0, 7.5, 10.0 / 

      data sig / 0.684, 0.684, 0.684, 0.684, 0.684, 0.684 ,0.684, 
     1           0.684, 0.684, 0.684, 0.684, 0.682, 0.686, 0.701, 
     2           0.712, 0.721, 0.727, 0.732, 0.755, 0.758, 0.809 /

 
C Find the requested spectral period and corresponding coefficients
      nPer = 21

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = period1(1)
         sigT      = sig(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period1(i) .and. specT .le. period1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'NGA sigma model (M5.5)'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period1(count1),period1(count2),sig(count1),sig(count2),
     +             specT,sigT,iflag)
 
 1011 period = specT                                                                                                              

C     Set the sigma value.
      sigma = sigT

      return                                                                    
      end                                                                                                                                                   

c ------------------------------------------------------------------           
C *** Atkinson (2010) Fena adjustment factors ****
c ------------------------------------------------------------------

      subroutine S06_Fena ( jbdist, specT, Factor )                                                

c    Atkinson (2008) Horizontal CEUS-NGA BA08 Vs=760m/sec

      implicit none

      real jbDist, period, Factor                              
      real c(12), d(12)
      real specT, cT, dT
      integer nper, count1, count2, iflag, i
      real per1(12)

      data Per1 / 0.00, 0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00, 
     1            2.00, 3.03, 5.00, 10.00 / 
      data c / 0.419, 0.417, 0.417, 0.245, 0.042, -0.078, -0.180, 
     1        -0.248, -0.214, -0.084, 0.0, 0.0 /
      data d / 0.00211, 0.00192, 0.00192, 0.00273, 0.00232, 0.00190, 
     1         0.00180, 0.00153, 0.00117, 0.00091, 0.00, 0.00  /

C Find the requested spectral period and corresponding coefficients
      nPer = 12
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period   = per1(1)
         cT      = c(1)
         dT      = d(1)
         goto 1011
      elseif (specT .ne. 0.0) then

C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. per1(i) .and. specT .le. per1(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
     
      write (*,*) 
      write (*,*) 'Atkinson (2010) Hor. Fena adustment model'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (per1(count1),per1(count2),c(count1),c(count2),
     +             specT,cT,iflag)
      call S24_interp (per1(count1),per1(count2),d(count1),d(count2),
     +             specT,dT,iflag)

 1011 period = specT                                                                                                              

c     Compute the ground motions adjustment factor.
      factor = cT + dT*jbdist
      factor = factor*alog(10.0)

      return                                                                    
      end                                                                                                                                                   

     
c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster01-Low: Mid-Continent, Functional Model1&3, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C1Low ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR, R1, R2, R3
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T, c15T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / -4.848, -4.848, -2.088, -4.659, -6.172, -9.963, -17.78, -25.86 /
      Data C2 / 1.621, 1.621, 1.266, 1.645, 2.026, 2.924, 4.751, 6.414 /
      Data C3 / -0.09426, -0.09426, -0.07828, -0.09066, -0.1225, -0.1829, -0.2995, -0.3878 /
      Data C4 / -1.801, -1.801, -2.145, -1.635, -1.804, -1.783, -1.752, -1.271 /
      Data C5 / 0.1038, 0.1038, 0.1427, 0.08719, 0.1179, 0.1187, 0.1182, 0.06127 /
      Data C6 / -2.658, -2.658, -2.327, -2.496, -2.993, -3.206, -3.106, -2.64 /
      Data C7 / 0.2127, 0.2127, 0.1083, 0.1976, 0.2776, 0.3074, 0.2928, 0.2503 /
      Data C8 / -2.531, -2.531, -2.188, -1.75, -2.138, -2.374, -2.226, -2.372 /
      Data C9 / 0.1731, 0.1731, 0.06719, 0.03638, 0.1328, 0.1777, 0.1522, 0.1943 /
      Data C10 / -0.0009935, -0.0009935, -0.001056, -0.002131, -0.000478, -0.00006192, -0.0002885, 0.0007184 /
      Data C11 / -0.00001544, -0.00001544, 0.0000762, 0.0001211, -0.0001307, -0.000144, -0.00004353, -0.0001961 /
      Data C12 / 1.737, 1.737, 2.378, 1.679, 1.829, 1.904, 2.028, 1.743 /
      Data C13 / 0.04032, 0.04032, -0.03753, 0.0453, 0.01974, 0.009728, -0.006279, 0.01898 /
      Data C14 / 70, 70, 70, 70, 70, 70, 70, 70 /
      Data C15 / 130, 130, 130, 130, 130, 130, 130, 130 /
      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         c15T     = c15(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster01-Low, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +             specT,c15T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster01-Low, MidC, Hor'                      

      RR = sqrt(dist*dist + (exp(c12T+c13T*M))**2.0 )
      R1 = min(alog(RR),alog(C14T))
      R2 = max( min(alog(RR/C14T),alog(c15T/c14T)),0.0)
      R3 = max(alog(RR/c15T),0.0)

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*m**2.0 + (c4T+c5T*M)*R1 + (c6T+c7T*M)*R2 + (c8T+c9T*M)*R3 + (c10T+c11T*M)*RR
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigM5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster01-Med: Mid-Continent, Functional Model1&3, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C1Med ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR, R1, R2, R3
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T, c15T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / -3.319, -3.319, -2.235, -3.135, -4.677, -8.615, -16.93, -24.33 /
      Data C2 / 1.268, 1.268, 1.212, 1.306, 1.557, 2.453, 4.368, 6.022 /
      Data C3 / -0.06382, -0.06382, -0.05816, -0.06238, -0.07678, -0.1346, -0.2544, -0.3514 /
      Data C4 / -1.892, -1.892, -1.88, -1.754, -1.658, -1.551, -1.378, -1.258 /
      Data C5 / 0.1107, 0.1107, 0.1037, 0.09942, 0.09415, 0.08433, 0.06594, 0.05294 /
      Data C6 / -2.389, -2.389, -2.384, -2.139, -1.985, -1.859, -1.62, -1.385 /
      Data C7 / 0.1557, 0.1557, 0.1443, 0.1316, 0.1268, 0.1176, 0.09272, 0.06748 /
      Data C8 / -2.497, -2.497, -2.642, -2.4, -2.071, -1.877, -1.642, -1.433 /
      Data C9 / 0.1647, 0.1647, 0.1607, 0.1523, 0.1347, 0.1184, 0.0951, 0.07481 /
      Data C10 / -0.0009028, -0.0009028, -0.0008213, -0.001312, -0.001459, -0.001222, -0.0009006, -0.0006909 /
      Data C11 / 0.00002638, 0.00002638, 0.0000274, 0.00002413, 0.00002092, 0.00001907, 0.00001591, 0.00001218 /
      Data C12 / 1.781, 1.781, 1.845, 1.783, 1.753, 1.749, 1.754, 1.714 /
      Data C13 / 0.03712, 0.03712, 0.03645, 0.03298, 0.02906, 0.02812, 0.02622, 0.02746 /
      Data C14 / 70, 70, 70, 70, 70, 70, 70, 70 /
      Data C15 / 130, 130, 130, 130, 130, 130, 130, 130 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         c15T     = c15(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster01-Med, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +             specT,c15T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster01-Med, MidC, Hor'                      

      RR = sqrt(dist*dist + (exp(c12T+c13T*M))**2.0 )
      R1 = min(alog(RR),alog(C14T))
      R2 = max( min(alog(RR/C14T),alog(c15T/c14T)),0.0)
      R3 = max(alog(RR/c15T),0.0)

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*m**2.0 + (c4T+c5T*M)*R1 + (c6T+c7T*M)*R2 + (c8T+c9T*M)*R3 + (c10T+c11T*M)*RR
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster01-High: Mid-Continent, Functional Model1&3, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C1High ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR, R1, R2, R3
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T, c15T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /


      Data C1 / -1.741, -1.741, -2.065, -1.502, -3.787, -8.305, -17.05, -22.98 /
      Data C2 / 0.9078, 0.9078, 1.094, 0.947, 1.183, 2.132, 4.106, 5.659 /
      Data C3 / -0.03331, -0.03331, -0.03549, -0.03383, -0.03079, -0.08381, -0.2037, -0.315 /
      Data C4 / -1.995, -1.995, -1.67, -1.901, -1.343, -1.005, -0.6753, -1.192 /
      Data C5 / 0.1193, 0.1193, 0.07244, 0.1158, 0.04361, -0.0003668, -0.03943, 0.03643 /
      Data C6 / -2.112, -2.112, -2.404, -1.758, -1.148, -0.8359, -0.4802, -0.1836 /
      Data C7 / 0.09747, 0.09747, 0.175, 0.06194, 0.002991, -0.0211, -0.05237, -0.1069 /
      Data C8 / -2.468, -2.468, -3.118, -3.063, -1.912, -1.206, -0.8714, -0.4637 /
      Data C9 / 0.157, 0.157, 0.2572, 0.2703, 0.1223, 0.03171, 0.008667, -0.04932 /
      Data C10 / -0.0007989, -0.0007989, -0.0005394, -0.0004617, -0.002657, -0.002793, -0.001954, -0.002172 /
      Data C11 / 0.00006614, 0.00006614, -0.00002835, -0.00007789, 0.0002064, 0.0002469, 0.000145, 0.0002316 /
      Data C12 / 1.841, 1.841, 1.356, 1.923, 1.333, 0.8951, 0.6201, 1.575 /
      Data C13 / 0.03107, 0.03107, 0.1037, 0.01455, 0.0931, 0.1573, 0.1938, 0.0529 /
      Data C14 / 70, 70, 70, 70, 70, 70, 70, 70 /
      Data C15 / 130, 130, 130, 130, 130, 130, 130, 130 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         c15T     = c15(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster01-High, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +             specT,c15T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster01-High, MidC, Hor'                      

      RR = sqrt(dist*dist + (exp(c12T+c13T*M))**2.0 )
      R1 = min(alog(RR),alog(C14T))
      R2 = max( min(alog(RR/C14T),alog(c15T/c14T)),0.0)
      R3 = max(alog(RR/c15T),0.0)

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*m**2.0 + (c4T+c5T*M)*R1 + (c6T+c7T*M)*R2 + (c8T+c9T*M)*R3 + (c10T+c11T*M)*RR
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

    
c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster02-Low: Mid-Continent, Functional Model2, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C2Low ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / 9.375, 9.375, 8.657, -4.012, -24.31, -20.56, -34.54, -40.21 /
      Data C2 / -3.378, -3.378, -2.236, 2.86, 10.7, 7.546, 13.34, 14.72 /
      Data C3 / 0.6098, 0.6098, 0.4056, -0.3493, -1.393, -0.8576, -1.674, -1.802 /
      Data C4 / -0.03561, -0.03561, -0.02389, 0.01202, 0.05803, 0.02973, 0.06694, 0.07154 /
      Data C5 / -3.081, -3.081, -3.485, -2.979, -2.23, -1.942, -1.97, -1.712 /
      Data C6 / 0.2714, 0.2714, 0.3281, 0.2882, 0.1765, 0.1713, 0.1351, 0.09196 /
      Data C7 / 0.2179, 0.2179, 0.1994, 0.2132, 0.1822, 0.1637, 0.2772, 0.2523 /
      Data C8 / -0.008817, -0.008817, -0.008201, -0.006515, -0.005683, -0.004415, -0.002111, -0.002285 /
      Data C9 / 0.0008775, 0.0008775, 0.0007111, 0.0006072, 0.0005755, 0.0004946, 0.0003024, 0.0003816 /
      Data C10 / -0.0002916, -0.0002916, -0.0002179, -0.0002223, -0.0001698, -0.000057, -0.0001211, -0.0001491 /
      Data C11 / 2.025, 2.025, 1.096, 0.6107, -0.4214, -1.526, 0.1708, 0.8019 /
      Data C12 / 0.09247, 0.09247, 0.3018, 0.3378, 0.5205, 0.6114, 0.3537, 0.2328 /
      Data C13 / 0.1355, 0.1355, 0.2184, 0.07885, -0.1417, -0.2739, -0.1513, -0.1221 /
      Data C14 / 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster02-Low, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster02-Low, MidC, Hor'                      

      CC = exp(c11T + c12T*min(M,c14T) + c13T*max((M-c14T),0.0))
      RR = dist + CC

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*M**2.0 + c4T*M**3.0 + 
     1     (c5T + c6T*min(M,c14T) +  c7T*max((M-c14T),0.0))*alog(RR) + 
     1     (c8T + c9T*min(M,c14T) + c10T*max((M-c14T),0.0))*RR    
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster02-Med: Mid-Continent, Functional Model2, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C2Med ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / -3.738, -3.738, -3.297, -1.548, -11.82, -15.27, -20.3, -26.51 /
      Data C2 / 2.318, 2.318, 2.551, 1.581, 5.062, 5.412, 6.868, 8.708 /
      Data C3 / -0.2634, -0.2634, -0.3267, -0.1876, -0.5837, -0.5534, -0.7453, -0.9685 /
      Data C4 / 0.008203, 0.008203, 0.01261, 0.006226, 0.02059, 0.01602, 0.02524, 0.0356 /
      Data C5 / -2.925, -2.925, -2.951, -2.905, -2.424, -2.246, -2.53, -2.496 /
      Data C6 / 0.3044, 0.3044, 0.3136, 0.3141, 0.2343, 0.2193, 0.2623, 0.2609 /
      Data C7 / 0.206, 0.206, 0.2043, 0.1875, 0.1604, 0.147, 0.1747, 0.155 /
      Data C8 / -0.005209, -0.005209, -0.005333, -0.004107, -0.003891, -0.003463, -0.002268, -0.001852 /
      Data C9 / 0.0002607, 0.0002607, 0.0002453, 0.0002144, 0.0002842, 0.0003825, 0.000357, 0.0003285 /
      Data C10 / -0.00004905, -0.00004905, -0.00004371, -0.00004316, -0.00006888, -0.00009736, -0.000095, -0.00009018 /
      Data C11 / 0.005719, 0.005719, 0.06974, 0.3146, -0.2284, -1.516, -1.321, -1.083 /
      Data C12 / 0.3629, 0.3629, 0.3513, 0.3066, 0.4259, 0.6191, 0.5725, 0.5213 /
      Data C13 / -0.07923, -0.07923, -0.0464, -0.06699, -0.1991, -0.2994, -0.2209, -0.212 /
      Data C14 / 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster02-Med, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster02-Med, MidC, Hor'                      

      CC = exp(c11T + c12T*min(M,c14T) + c13T*max((M-c14T),0.0))
      RR = dist + CC

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*M**2.0 + c4T*M**3.0 + 
     1     (c5T + c6T*min(M,c14T) +  c7T*max((M-c14T),0.0))*alog(RR) + 
     1     (c8T + c9T*min(M,c14T) + c10T*max((M-c14T),0.0))*RR    
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       
        
c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster02-High: Mid-Continent, Functional Model2, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C2High ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / -7.709, -7.709, 0.3463, -7.616, -22.82, -13.7, -4.622, -6.92 /
      Data C2 / 3.787, 3.787, 0.2886, 4.28, 10.95, 5.552, 0.01532, 0.01103 /
      Data C3 / -0.5, -0.5, 0.004643, -0.6293, -1.625, -0.7189, 0.1854, 0.2609 /
      Data C4 / 0.0217, 0.0217, -0.001578, 0.03088, 0.08046, 0.03353, -0.01315, -0.01948 /
      Data C5 / -2.693, -2.693, -2.5, -2.651, -2.544, -2.614, -3.128, -3.291 /
      Data C6 / 0.3243, 0.3243, 0.3002, 0.3058, 0.2771, 0.2899, 0.3995, 0.4352 /
      Data C7 / 0.1614, 0.1614, 0.1519, 0.1674, 0.1768, 0.09939, 0.05889, 0.05389 /
      Data C8 / -0.003865, -0.003865, -0.004288, -0.003348, -0.002797, -0.002955, -0.003016, -0.002403 /
      Data C9 / -0.00002854, -0.00002854, 0.00002666, 0.00006314, 0.00008068, 0.000295, 0.0004848, 0.0004176 /
      Data C10 / 0.0001343, 0.0001343, 0.0001158, 0.00004686, -0.0000288, -0.00005291, -0.00005643, -0.00006663 /
      Data C11 / -3.671, -3.671, -0.668, -0.4119, -0.1841, -1.3, -2.07, -1.904 /
      Data C12 / 0.8593, 0.8593, 0.2715, 0.3345, 0.363, 0.5627, 0.6302, 0.5819 /
      Data C13 / -0.3748, -0.3748, -0.3441, -0.3332, -0.5366, -0.2887, -0.1715, -0.1641 /
      Data C14 / 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.8 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster02-High, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster02-High, MidC, Hor'                      

      CC = exp(c11T + c12T*min(M,c14T) + c13T*max((M-c14T),0.0))
      RR = dist + CC

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*M**2.0 + c4T*M**3.0 + 
     1     (c5T + c6T*min(M,c14T) +  c7T*max((M-c14T),0.0))*alog(RR) + 
     1     (c8T + c9T*min(M,c14T) + c10T*max((M-c14T),0.0))*RR    
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster03-Low: Mid-Continent, Functional Model1&3, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C3Low ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR, R1, R2, R3
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T, c15T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / -3.387, -3.387, -2.462, -4.478, -7.819, -12.56, -20.1, -23.41 /
      Data C2 / 1.892, 1.892, 1.837, 2.011, 2.732, 3.837, 5.511, 5.931 /
      Data C3 / -0.1469, -0.1469, -0.1377, -0.1461, -0.1941, -0.2678, -0.3759, -0.3897 /
      Data C4 / -2.544, -2.544, -2.45, -2.074, -1.898, -1.844, -1.865, -1.996 /
      Data C5 / 0.1571, 0.1571, 0.1345, 0.1276, 0.1171, 0.1165, 0.1269, 0.1511 /
      Data C6 / 0.3513, 0.3513, 0.3929, 0.4439, 0.3467, 0.3279, 0.3541, 0.3506 /
      Data C7 / -0.0488, -0.0488, -0.06395, -0.04426, -0.02948, -0.02606, -0.02821, -0.02736 /
      Data C8 / -2.695, -2.695, -3.581, -2.023, -1.137, -0.8706, -0.8142, -0.8959 /
      Data C9 / 0.2068, 0.2068, 0.232, 0.1365, 0.07169, 0.04971, 0.04863, 0.0667 /
      Data C10 / -0.001827, -0.001827, -0.0009334, -0.003142, -0.003267, -0.002626, -0.001796, -0.00139 /
      Data C11 / 0.00005335, 0.00005335, 0.00003054, 0.00005509, 0.00005484, 0.00005281, 0.00004872, 0.00003937 /
      Data C12 / 3.038, 3.038, 2.92, 2.973, 2.97, 2.955, 2.912, 2.847 /
      Data C13 / -0.1709, -0.1709, -0.1392, -0.1729, -0.1822, -0.1846, -0.1788, -0.1698 /
      Data C14 / 70, 70, 70, 70, 70, 70, 70, 70 /
      Data C15 / 140, 140, 140, 140, 140, 140, 140, 140 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         c15T     = c15(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster03-Low, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +             specT,c15T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster03-Low, MidC, Hor'                      

      RR = sqrt(dist*dist + (exp(c12T+c13T*M))**2.0 )
      R1 = min(alog(RR),alog(C14T))
      R2 = max( min(alog(RR/C14T),alog(c15T/c14T)),0.0)
      R3 = max(alog(RR/c15T),0.0)

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*m**2.0 + (c4T+c5T*M)*R1 + (c6T+c7T*M)*R2 + (c8T+c9T*M)*R3 + (c10T+c11T*M)*RR
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster03-Med: Mid-Continent, Functional Model1&3, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C3Med ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR, R1, R2, R3
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T, c15T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / -2.27, -2.27, -1.35, -3.44, -6.752, -11.5, -19.01, -22.22 /
      Data C2 / 1.679, 1.679, 1.624, 1.818, 2.532, 3.635, 5.311, 5.726 /
      Data C3 / -0.1277, -0.1277, -0.1184, -0.1282, -0.1758, -0.2495, -0.3578, -0.372 /
      Data C4 / -2.649, -2.649, -2.554, -2.176, -2.015, -1.961, -1.982, -2.133 /
      Data C5 / 0.1655, 0.1655, 0.1425, 0.1357, 0.1269, 0.1264, 0.1366, 0.1624 /
      Data C6 / -0.04737, -0.04737, -0.007665, 0.1136, 0.06214, 0.04284, 0.027, 0.03288 /
      Data C7 / -0.009272, -0.009272, -0.02411, -0.01666, -0.003979, -0.0004893, 0.0008474, -0.0004211 /
      Data C8 / -2.61, -2.61, -3.471, -1.999, -1.122, -0.8568, -0.8034, -0.8856 /
      Data C9 / 0.1876, 0.1876, 0.2084, 0.1279, 0.06538, 0.04361, 0.04278, 0.06077 /
      Data C10 / -0.001255, -0.001255, -0.0004193, -0.002493, -0.002706, -0.002063, -0.001216, -0.0007902 /
      Data C11 / 0.00003347, 0.00003347, 0.00002056, 0.00002176, 0.00002053, 0.000018, 0.00001282, 0.000004342 /
      Data C12 / 3.049, 3.049, 2.934, 2.98, 2.979, 2.965, 2.921, 2.859 /
      Data C13 / -0.1689, -0.1689, -0.1382, -0.1694, -0.1778, -0.1801, -0.1738, -0.1638 /
      Data C14 / 70, 70, 70, 70, 70, 70, 70, 70 /
      Data C15 / 140, 140, 140, 140, 140, 140, 140, 140 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         c15T     = c15(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster03-Med, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +             specT,c15T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster03-Med, MidC, Hor'                      

      RR = sqrt(dist*dist + (exp(c12T+c13T*M))**2.0 )
      R1 = min(alog(RR),alog(C14T))
      R2 = max( min(alog(RR/C14T),alog(c15T/c14T)),0.0)
      R3 = max(alog(RR/c15T),0.0)

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*m**2.0 + (c4T+c5T*M)*R1 + (c6T+c7T*M)*R2 + (c8T+c9T*M)*R3 + (c10T+c11T*M)*RR

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c -----------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster03-High: Mid-Continent, Functional Model1&3, Horizontal ***
C -----------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C3High ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                
      real lnY, m, dist, period1, RR, R1, R2, R3
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, c8T, c9T, C10T
      real c11T, c12T, c13T, c14T, c15T
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), c8(MAXPER), c9(MAXPER), c10(MAXPER)
      real c11(MAXPER), c12(MAXPER), c13(MAXPER), c14(MAXPER), c15(MAXPER)
      real period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data Period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /
      Data C1 / -1.152, -1.152, -0.2207, -2.407, -5.689, -10.43, -17.92, -21.03 /
      Data C2 / 1.466, 1.466, 1.409, 1.627, 2.332, 3.435, 5.11, 5.522 /
      Data C3 / -0.1084, -0.1084, -0.09901, -0.1104, -0.1575, -0.2312, -0.3396, -0.3542 /
      Data C4 / -2.755, -2.755, -2.661, -2.276, -2.131, -2.078, -2.1, -2.271 /
      Data C5 / 0.1739, 0.1739, 0.151, 0.1435, 0.1365, 0.1362, 0.1465, 0.1738 /
      Data C6 / -0.4456, -0.4456, -0.4064, -0.2181, -0.226, -0.2418, -0.2985, -0.285 /
      Data C7 / 0.03018, 0.03018, 0.01547, 0.01117, 0.02205, 0.02506, 0.02971, 0.02658 /
      Data C8 / -2.525, -2.525, -3.362, -1.974, -1.107, -0.8422, -0.7946, -0.8754 /
      Data C9 / 0.1684, 0.1684, 0.1852, 0.1192, 0.05907, 0.03742, 0.0372, 0.05484 /
      Data C10 / -0.0006815, -0.0006815, 0.0001003, -0.001846, -0.002142, -0.001503, -0.0006309, -0.0001901 /
      Data C11 / 0.00001348, 0.00001348, 0.000009821, -0.00001132, -0.00001408, -0.00001647, -0.00002377, -0.00003068 /
      Data C12 / 3.059, 3.059, 2.95, 2.984, 2.984, 2.972, 2.931, 2.87 /
      Data C13 / -0.167, -0.167, -0.1377, -0.1656, -0.1734, -0.1756, -0.1693, -0.1583 /
      Data C14 / 70, 70, 70, 70, 70, 70, 70, 70 /
      Data C15 / 140, 140, 140, 140, 140, 140, 140, 140 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         c8T     = c8(1)
         c9T     = c9(1)
         c10T     = c10(1)
         c11T     = c11(1)
         c12T     = c12(1)
         c13T     = c13(1)
         c14T     = c14(1)
         c15T     = c15(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster03-High, Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),c8(count1),c8(count2),
     +             specT,c8T,iflag)
      call S24_interp (period(count1),period(count2),c9(count1),c9(count2),
     +             specT,c9T,iflag)
      call S24_interp (period(count1),period(count2),c10(count1),c10(count2),
     +             specT,c10T,iflag)
      call S24_interp (period(count1),period(count2),c11(count1),c11(count2),
     +             specT,c11T,iflag)
      call S24_interp (period(count1),period(count2),c12(count1),c12(count2),
     +             specT,c12T,iflag)
      call S24_interp (period(count1),period(count2),c13(count1),c13(count2),
     +             specT,c13T,iflag)
      call S24_interp (period(count1),period(count2),c14(count1),c14(count2),
     +             specT,c14T,iflag)
      call S24_interp (period(count1),period(count2),c15(count1),c15(count2),
     +             specT,c15T,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)

 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster03-High, MidC, Hor'                      

      RR = sqrt(dist*dist + (exp(c12T+c13T*M))**2.0 )
      R1 = min(alog(RR),alog(C14T))
      R2 = max( min(alog(RR/C14T),alog(c15T/c14T)),0.0)
      R3 = max(alog(RR/c15T),0.0)

C     Compute median ground motions
      lnY = c1T + c2T*M + c3T*m**2.0 + (c4T+c5T*M)*R1 + (c6T+c7T*M)*R2 + (c8T+c9T*M)*R3 + (c10T+c11T*M)*RR
      
C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       


c ---------------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster04-Low (Rift): Mid-Continent, Functional Model4, Horizontal ***
C ---------------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C4RLow ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                                                                                                
      real lnY, m, dist, period1, r1, m1, m2, d, d1
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, hT
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), h(MAXPER), period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / -0.2099, -0.2099, 0.4771, 0.4391, 0.3459, 0.1753, -0.7732, -1.657 /
      Data C2 / 0.6754, 0.6754, 0.6754, 0.6754, 0.6743, 0.6744, 0.6782, 0.6843 /
      Data C3 / -0.6303, -0.6303, -0.6303, -0.6303, -0.6251, -0.6103, -0.6406, -0.6616 /
      Data C4 / 0.08405, 0.08405, 0.08405, 0.08405, 0.08303, 0.083, 0.08261, 0.08232 /
      Data C5 / -0.005586, -0.005586, -0.005586, -0.005586, -0.005461, -0.005159, -0.004133, -0.002758 /
      Data C6 / -0.3669, -0.3669, -0.3669, -0.3669, -0.3867, -0.4668, -0.6525, -0.8382 /
      Data C7 / -0.01757, -0.01757, -0.01757, -0.01757, -0.01793, -0.06973, -0.1198, -0.1574 /
      Data h / 5.711, 5.711, 5.711, 5.711, 5.662, 5.657, 5.658, 5.595 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         hT      = h(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster04-Low (Rift), Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster04-Low-Rift, MidC, Hor'                      
      r1 = 50.0
      m1 = 6.4
      m2 = 8.5

      d = sqrt (dist*dist + hT*hT )
      d1 = sqrt (r1*r1 + hT*hT )

      if (dist .lt. r1) then
          lnY = c1T + c2T*(M-m1) + c3T*alog(d)+ c4T*(M-m1)*alog(d) + c5T*dist + c7T*(m2-M)**2.0
      else
          lnY = c1T + c2T*(M-m1) + c3T*alog(d1)+ c4T*(M-m1)*alog(d) + c5T*dist + c6T*(alog(d)-alog(d1)) + c7T*(m2-M)**2.0
      endif

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       
        


c ---------------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster04-Med (Rift): Mid-Continent, Functional Model4, Horizontal ***
C ---------------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C4RMed ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                                                                                                
      real lnY, m, dist, period1, r1, m1, m2, d, d1
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, hT
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), h(MAXPER), period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / 0.239, 0.239, 0.926, 0.888, 0.793, 0.622, -0.307, -1.132 /
      Data C2 / 0.805, 0.805, 0.805, 0.805, 0.805, 0.805, 0.805, 0.805 /
      Data C3 / -0.679, -0.679, -0.679, -0.679, -0.679, -0.664, -0.696, -0.728 /
      Data C4 / 0.0861, 0.0861, 0.0861, 0.0861, 0.0861, 0.0861, 0.0861, 0.0861 /
      Data C5 / -0.00498, -0.00498, -0.00498, -0.00498, -0.00498, -0.00468, -0.00362, -0.00221 /
      Data C6 / -0.477, -0.477, -0.477, -0.477, -0.477, -0.557, -0.755, -0.946 /
      Data C7 / 0, 0, 0, 0, 0, -0.0518, -0.102, -0.14 /
      Data h / 6, 6, 6, 6, 6, 6, 6, 6 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         hT      = h(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster04-Med (Rift), Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster04-Med-Rift, MidC, Hor'                      
      r1 = 50.0
      m1 = 6.4
      m2 = 8.5

      d = sqrt (dist*dist + hT*hT )
      d1 = sqrt (r1*r1 + hT*hT )

      if (dist .lt. r1) then
          lnY = c1T + c2T*(M-m1) + c3T*alog(d)+ c4T*(M-m1)*alog(d) + c5T*dist + c7T*(m2-M)**2.0
      else
          lnY = c1T + c2T*(M-m1) + c3T*alog(d1)+ c4T*(M-m1)*alog(d) + c5T*dist + c6T*(alog(d)-alog(d1)) + c7T*(m2-M)**2.0
      endif

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       
 
c ---------------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster04-High (Rift): Mid-Continent, Functional Model4, Horizontal ***
C ---------------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C4RHigh ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                                                                                                
      real lnY, m, dist, period1, r1, m1, m2, d, d1
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, hT
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), h(MAXPER), period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / 0.6912, 0.6912, 1.378, 1.34, 1.244, 1.073, 0.1632, -0.6014 /
      Data C2 / 0.9348, 0.9348, 0.9348, 0.9348, 0.9359, 0.9358, 0.932, 0.9259 /
      Data C3 / -0.7286, -0.7286, -0.7286, -0.7286, -0.734, -0.7188, -0.7525, -0.796 /
      Data C4 / 0.0881, 0.0881, 0.0881, 0.0881, 0.08912, 0.08915, 0.08956, 0.08983 /
      Data C5 / -0.004375, -0.004375, -0.004375, -0.004375, -0.0045, -0.004201, -0.003107, -0.001662 /
      Data C6 / -0.587, -0.587, -0.587, -0.587, -0.5671, -0.647, -0.8574, -1.054 /
      Data C7 / 0.01757, 0.01757, 0.01757, 0.01757, 0.01793, -0.03387, -0.08425, -0.1226 /
      Data h / 6.275, 6.275, 6.275, 6.275, 6.319, 6.324, 6.322, 6.38 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         hT      = h(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster04-High (Rift), Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster04-High-Rift, MidC, Hor'                      
      r1 = 50.0
      m1 = 6.4
      m2 = 8.5

      d = sqrt (dist*dist + hT*hT )
      d1 = sqrt (r1*r1 + hT*hT )

      if (dist .lt. r1) then
          lnY = c1T + c2T*(M-m1) + c3T*alog(d)+ c4T*(M-m1)*alog(d) + c5T*dist + c7T*(m2-M)**2.0
      else
          lnY = c1T + c2T*(M-m1) + c3T*alog(d1)+ c4T*(M-m1)*alog(d) + c5T*dist + c6T*(alog(d)-alog(d1)) + c7T*(m2-M)**2.0
      endif

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       
        
 
c ------------------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster04-Low (NonRift): Mid-Continent, Functional Model4, Horizontal ***
C ------------------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C4NRLow ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                                                                                                
      real lnY, m, dist, period1, r1, m1, m2, d, d1
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, hT
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), h(MAXPER), period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / -0.03467, -0.03467, 0.6464, 0.6183, 0.5269, 0.4, -0.6087, -1.46 /
      Data C2 / 0.6775, 0.6775, 0.6775, 0.6775, 0.6762, 0.6763, 0.6801, 0.6861 /
      Data C3 / -0.6782, -0.6782, -0.6783, -0.6782, -0.673, -0.6731, -0.6826, -0.6866 /
      Data C4 / 0.06322, 0.06322, 0.06322, 0.06322, 0.06224, 0.06222, 0.06181, 0.06156 /
      Data C5 / -0.006614, -0.006614, -0.006614, -0.006614, -0.006489, -0.005858, -0.004492, -0.003727 /
      Data C6 / -0.1914, -0.1914, -0.1914, -0.1914, -0.2113, -0.3334, -0.557, -0.5948 /
      Data C7 / -0.01757, -0.01757, -0.01757, -0.01757, -0.01793, -0.06973, -0.1198, -0.1574 /
      Data h / 5.71, 5.71, 5.711, 5.71, 5.664, 5.665, 5.66, 5.591 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         hT      = h(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster04-Low (NonRift), Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster04-Low-NonRift, MidC, Hor'                      
      r1 = 50.0
      m1 = 6.4
      m2 = 8.5

      d = sqrt (dist*dist + hT*hT )
      d1 = sqrt (r1*r1 + hT*hT )

      if (dist .lt. r1) then
          lnY = c1T + c2T*(M-m1) + c3T*alog(d)+ c4T*(M-m1)*alog(d) + c5T*dist + c7T*(m2-M)**2.0
      else
          lnY = c1T + c2T*(M-m1) + c3T*alog(d1)+ c4T*(M-m1)*alog(d) + c5T*dist + c6T*(alog(d)-alog(d1)) + c7T*(m2-M)**2.0
      endif

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       
        
c ------------------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster04-Med (NonRift): Mid-Continent, Functional Model4, Horizontal ***
C ------------------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C4NRMed ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                                                                                                
      real lnY, m, dist, period1, r1, m1, m2, d, d1
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, hT
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), h(MAXPER), period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / 0.418, 0.418, 1.099, 1.071, 0.978, 0.851, -0.139, -0.932 /
      Data C2 / 0.808, 0.808, 0.808, 0.808, 0.808, 0.808, 0.808, 0.808 /
      Data C3 / -0.728, -0.728, -0.728, -0.728, -0.728, -0.728, -0.739, -0.754 /
      Data C4 / 0.0651, 0.0651, 0.0651, 0.0651, 0.0651, 0.0651, 0.0651, 0.0651 /
      Data C5 / -0.00601, -0.00601, -0.00601, -0.00601, -0.00601, -0.00538, -0.00398, -0.00318 /
      Data C6 / -0.301, -0.301, -0.301, -0.301, -0.301, -0.423, -0.659, -0.702 /
      Data C7 / 0, 0, 0, 0, 0, -0.0518, -0.102, -0.14 /
      Data h / 6, 6, 6, 6, 6, 6, 6, 6 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         hT      = h(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster04-Med (NonRift), Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster04-Med-NonRift, MidC, Hor'                      
      r1 = 50.0
      m1 = 6.4
      m2 = 8.5

      d = sqrt (dist*dist + hT*hT )
      d1 = sqrt (r1*r1 + hT*hT )

      if (dist .lt. r1) then
          lnY = c1T + c2T*(M-m1) + c3T*alog(d)+ c4T*(M-m1)*alog(d) + c5T*dist + c7T*(m2-M)**2.0
      else
          lnY = c1T + c2T*(M-m1) + c3T*alog(d1)+ c4T*(M-m1)*alog(d) + c5T*dist + c6T*(alog(d)-alog(d1)) + c7T*(m2-M)**2.0
      endif

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

c ------------------------------------------------------------------------------------------------
C *** EPRI Update (2013) Cluster04-High (NonRift): Mid-Continent, Functional Model4, Horizontal ***
C ------------------------------------------------------------------------------------------------
                                                                               
      subroutine S06_EPRI13C4NRHigh ( m, dist, lnY, specT,                      
     1                  attenName, period1,iflag, sig )                                    
                                                                                                                                                                
      real lnY, m, dist, period1, r1, m1, m2, d, d1
      real specT, c1T, c2T, c3T, c4T, c5T, c6T, c7T, hT
      integer nper, count1, count2, iflag
      real sigM5T, sigM6T, sigM7T, sig
      character*80 attenName                                                    
                                                                                
      parameter (MAXPER=8)                                                      
      real c1(MAXPER), c2(MAXPER), c3(MAXPER), c4(MAXPER), c5(MAXPER)         
      real c6(MAXPER), c7(MAXPER), h(MAXPER), period(MAXPER)
      real sigM5(MAXPER), sigM6(MAXPER), sigM7(MAXPER)

      Data period / 0.0, 0.01, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0 /

      Data C1 / 0.8738, 0.8738, 1.555, 1.527, 1.433, 1.306, 0.3346, -0.398 /
      Data C2 / 0.9387, 0.9387, 0.9387, 0.9387, 0.9399, 0.9398, 0.936, 0.93 /
      Data C3 / -0.7786, -0.7786, -0.7786, -0.7786, -0.784, -0.784, -0.7965, -0.8229 /
      Data C4 / 0.06695, 0.06695, 0.06695, 0.06695, 0.06793, 0.06795, 0.06837, 0.06861 /
      Data C5 / -0.005406, -0.005406, -0.005406, -0.005406, -0.005531, -0.004903, -0.003468, -0.002633 /
      Data C6 / -0.4105, -0.4105, -0.4105, -0.4105, -0.3906, -0.5125, -0.7609, -0.809 /
      Data C7 / 0.01757, 0.01757, 0.01757, 0.01757, 0.01793, -0.03387, -0.08425, -0.1226 /
      Data h / 6.276, 6.276, 6.276, 6.276, 6.318, 6.318, 6.321, 6.384 /

      Data sigM5 / 0.68, 0.68, 0.74, 0.74, 0.72, 0.72, 0.75, 0.78 /
      Data sigM6 / 0.63, 0.63, 0.69, 0.69, 0.68, 0.69, 0.74, 0.78 /
      Data sigM7 / 0.60, 0.60, 0.67, 0.67, 0.66, 0.67, 0.73, 0.77 /

C Find the requested spectral period and corresponding coefficients
      nPer = 8

C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         period1 = period(1)
         c1T     = c1(1)
         c2T     = c2(1)
         c3T     = c3(1)
         c4T     = c4(1)
         c5T     = c5(1)
         c6T     = c6(1)
         c7T     = c7(1)
         hT      = h(1)
         sigM5T   = sigM5(1)
         sigM6T   = sigM6(1)
         sigM7T   = sigM7(1)
         goto 1011
      elseif (specT .ne. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif
        
      write (*,*) 
      write (*,*) 'EPRI Update (2013) Cluster04-High (NonRift), Mid-C, Horizontal'
      write (*,*) 'is not defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +             specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c2(count1),c2(count2),
     +             specT,c2T,iflag)
      call S24_interp (period(count1),period(count2),c3(count1),c3(count2),
     +             specT,c3T,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +             specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +             specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +             specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c7(count1),c7(count2),
     +             specT,c7T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +             specT,hT,iflag)
      call S24_interp (period(count1),period(count2),sigM5(count1),sigM5(count2),
     +             specT,sigM5T,iflag)
      call S24_interp (period(count1),period(count2),sigM6(count1),sigM6(count2),
     +             specT,sigM6T,iflag)
      call S24_interp (period(count1),period(count2),sigM7(count1),sigM7(count2),
     +             specT,sigM7T,iflag)
 1011 period1 = specT                                                                                                              
                                                                                
c     Set atten name                                                            
      attenName = 'EPRI Update(2013), Cluster04-High-NonRift, MidC, Hor'                      
      r1 = 50.0
      m1 = 6.4
      m2 = 8.5

      d = sqrt (dist*dist + hT*hT )
      d1 = sqrt (r1*r1 + hT*hT )

      if (dist .lt. r1) then
          lnY = c1T + c2T*(M-m1) + c3T*alog(d)+ c4T*(M-m1)*alog(d) + c5T*dist + c7T*(m2-M)**2.0
      else
          lnY = c1T + c2T*(M-m1) + c3T*alog(d1)+ c4T*(M-m1)*alog(d) + c5T*dist + c6T*(alog(d)-alog(d1)) + c7T*(m2-M)**2.0
      endif

C     Compute the Sigma value.
      if (m .le. 5.0) then
         sig = sigm5T
      elseif (m .ge. 7.0) then
         sig = sigM7T
      elseif (m .gt. 5.0 .and. m .lt. 6.0) then
         sig = sigM5T + (sigM6T - sigM5T)*(m-5.0)/(6.0-5.0)
      elseif (m .ge. 6.0 .and. m .lt. 7.0) then
         sig = sigM6T + (sigM7T - sigM6T)*(m-6.0)/(7.0-6.0)
      endif

c     Convert to spectral acceleration in gal                                   
      lnY = lnY + 6.89                                                          
                                                                                
      return                                                                    
      end                                                                       

                                      
     
