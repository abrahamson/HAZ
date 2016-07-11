
c -------------------------------------------------------------------
      subroutine DirJWL_V3Uni (specT, Rrup, Rx1, Ry, Ruplen, Mag, mech, 
     1        RupWidth, dip, HWFlag, medadj, sigadj )

      real Rx, Ry, Ruplen, Mag, medadj, sigadj, mech, per, Rrup,
     1     taperDist, taperMag, Cos2TheaAvg, Rymin, specT, RyRatio, RyRatioBar
      real period (22), c8b(22), c8bT
      real c8rev(22), c8revT, c8org(22), c8orgT
      real b0m, b1m, b2m, b3m, b4m, b5m, b6m, b1am, bMm, sigm
      real b0s, b1s, b2s, b3s, b4s, b5s, b6s, b1as, bMs, sigs
      real MDF, SDF, MDA, SDA

      real b0_ssm, b1_ssm, b2_ssm, b3_ssm, bM_ssm
      real b0_sss, b1_sss, b2_sss, b3_sss, bM_sss

      real b0_rvm, b1_rvm, b2_rvm, b3_rvm, b4_rvm, b5_rvm, b6_rvm, b7_rvm, b8_rvm, bM_rvm
      real b0_rvs, b1_rvs, b2_rvs, b3_rvs, b4_rvs, b5_rvs, bM_rvs      

      real RupWidth, dip, cos2phiAvg, RyRatioRV
      real thick, RxThick, RyThick, RyminThick 
      integer iflag, count1, count2, HWflag   

      data period     /
     1              0.0000, 0.0100, 0.0200,  0.0300, 0.0500,
     1              0.0750, 0.1000, 0.1500,
     1              0.2000, 0.2500, 0.3000,  0.4000, 0.5000,
     1              0.7500, 1.0000, 1.5000,  2.0000, 3.0000,
     1              4.0000, 5.0000, 7.5000, 10.0000 /
      data c8rev / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1       0.0991, 0.1982, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 / 
      data c8b / 0.4833, 0.4833, 1.2144, 1.6421, 2.181, 2.6087, 2.9122, 3.3399, 3.6434, 3.8787, 4.0711, 
     1       4.3745, 4.6099, 5.0376, 5.3411, 5.7688, 6.0723, 6.5, 6.8035, 7.0389, 7.4666, 7.77 / 
      data c8org / 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 
     1         0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 /

c     Flip sign of Rx (JWL model used different convention
c     THis will be corrected in the final documentation)
      Rx = -Rx1


C     Mean Model Coefficients - StrikeSlip
      b0_ssm = -0.10489
      b1_ssm = 0.0034569
      b2_ssm = 0.32012
      b3_ssm = 0.053815
      bM_ssm = -0.29074

C     Sigma Model Coefficients - StrikeSlip
      b0_sss = 0.046052
      b1_sss = 0.21626
      b2_sss = 0.42046
      b3_sss = -0.55467
      bM_sss = -0.15878

C     Mean Model Coefficients - Reverse
      b0_rvm = -0.084732
      b1_rvm = -0.024127
      b2_rvm = 0.032329
      b3_rvm = 0.070904
      b4_rvm = 0.083838
      b5_rvm = 0.014869
      b6_rvm = -0.11502
      b7_rvm = 0.010686
      b8_rvm = 0.081255
      bM_rvm = -0.26578

C     Sigma Model Coefficients - Reverse
      b0_rvs = 0.051969
      b1_rvs = -0.016977
      b2_rvs = -0.029564
      b3_rvs = 0.011029
      b4_rvs = -0.037247
      b5_rvs = 0.073575
      bM_rvs = 0.0


C     First interpolate terms for given spectral period. 
C Find the requested spectral period and corresponding coefficients
      nPer = 22
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         c8revT  = c8rev(1)
         c8bT  = c8b(1)
         c8orgT  = c8org(1)

         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Interpolate the coefficients for the requested spectral period.
 1020       call interp (period(count1),period(count2),c8rev(count1),c8rev(count2),
     +                   specT,c8revT,iflag)
            call interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
            call interp (period(count1),period(count2),c8org(count1),c8org(count2),
     +                   specT,c8orgT,iflag)

 1011 period1 = specT                                                                                                              

C     Compute the mag and distance tapers. 
      taperdist = max(1.0-max(Rrup-40.0,0.0)/30.0,0.0)
      taperMag = min(max(mag-5.5,0.0)/0.8,1.0)
      RyRatio = min(abs(Ry)/(Ruplen/2.0),1.0)

C     Now calculate the parameters for the model given input values.
C     Classification is strike-slip or normal --> Strike-slip coefficients
C                       reverse or oblique    --> Reverse coefficients
      if (mech .le. 0.0) then

         cos2thetaAvg = ( ( (Ry+RupLen/2.0)-2.0*abs(Rx)*atan2((Ry+Ruplen/2.0),abs(Rx)) ) -
     1          ( (Ry-Ruplen/2.0)-2.0*abs(Rx)*atan2((Ry-Ruplen/2.0),abs(Rx)) ) ) /Ruplen

         MDF = b0_ssm + b1_ssm*(max(RyRatio*cos2ThetaAvg,-0.5)) + b2_ssm*(max(RyRatio*cos2ThetaAvg,-0.5))**2.0 +
     1            b3_ssm*(max(RyRatio*cos2ThetaAvg,-0.5))**3.0    
         SDF = b0_sss + b1_sss*(max(RyRatio*cos2ThetaAvg,-0.5)) + b2_sss*(max(RyRatio*cos2ThetaAvg,-0.5))**2.0 +
     1            b3_sss*(max(RyRatio*cos2ThetaAvg,-0.5))**3.0         

         Medadj = MDF*exp(bM_ssm*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag
         SDA =    SDF*exp(bM_sss*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag

      else
         cosdip = cos(dip*3.14159/180.0)
         sindip = sin(dip*3.14159/180.0)
         cos2thetaAvg = ( ( (Ry+RupLen/2.0)-2.0*abs(Rx*sindip)*atan2((Ry+Ruplen/2.0),abs(Rx*sindip)) ) -
     1          ( (Ry-Ruplen/2.0)-2.0*abs(Rx*sindip)*atan2((Ry-Ruplen/2.0),abs(Rx*sindip)) ) ) /Ruplen

         cos2phiAvg = ( ( (Rx*sindip)-2.0*abs(Rx*cosdip)*atan2(Rx*sindip,abs(Rx*cosdip))) -
     1         ( (Rx*sindip-RupWidth)-2.0*abs(Rx*cosdip)*atan2((Rx*sindip-RupWidth),abs(Rx*cosdip))) ) / RupWidth

         MDF = b0_rvm + b1_rvm*cos2phiAvg + b2_rvm*cos2phiAvg**2.0 + b3_rvm*cos2phiAvg**3.0 +
     1           b4_rvm*cos2thetaAvg + b5_rvm*cos2thetaAvg**2.0 + b6_rvm*cos2thetaAvg**3.0 + 
     1           b7_rvm*RyRatio + b8_rvm*RyRatio**2.0

         SDF = b0_rvs+ b1_rvs*cos2phiAvg + b2_rvs*cos2phiAvg**2.0 + b3_rvs*cos2phiAvg**3.0 +
     1           b4_rvs*RyRatio + b5_rvs*RyRatio**2.0
         Medadj = MDF*exp(bM_rvs*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag
         SDA =    SDF*exp(bM_rvs*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag

      endif

C     Check for sigma adjustment values less than 0.0, if so set equal to 0.0
      if (sda .lt. 0.0) SDA=0.0
      sigadj = SDA

      return
      end

c -------------------------------------------------------------------
      subroutine DirJWL_V3Pre (specT, Rrup, Rx, Ry, Ruplen, Mag, mech, 
     1        RupWidth, dip, HWFlag, medadj, sigadj )

      real Rx, Ry, Ruplen, Mag, medadj, sigadj, mech, per, Rrup,
     1     taperDist, taperMag, Cos2TheaAvg, Rymin, specT, RyRatio, RyRatioBar
      real period (22), c8b(22), c8bT
      real c8rev(22), c8revT, c8org(22), c8orgT
      real b0m, b1m, b2m, b3m, b4m, b5m, b6m, b1am, bMm, sigm
      real b0s, b1s, b2s, b3s, b4s, b5s, b6s, b1as, bMs, sigs
      real MDF, SDF, MDA, SDA

      real b0_ssm, b1_ssm, b2_ssm, b3_ssm, bM_ssm
      real b0_sss, b1_sss, b2_sss, b3_sss, bM_sss

      real b0_rvm, b1_rvm, b2_rvm, b3_rvm, b4_rvm, b5_rvm, b6_rvm, b7_rvm, b8_rvm, bM_rvm
      real b0_rvs, b1_rvs, b2_rvs, b3_rvs, b4_rvs, b5_rvs, bM_rvs      

      real RupWidth, dip, cos2phiAvg, RyRatioRV
      real thick, RxThick, RyThick, RyminThick 
      integer iflag, count1, count2, HWflag   

      data period     /
     1              0.0000, 0.0100, 0.0200,  0.0300, 0.0500,
     1              0.0750, 0.1000, 0.1500,
     1              0.2000, 0.2500, 0.3000,  0.4000, 0.5000,
     1              0.7500, 1.0000, 1.5000,  2.0000, 3.0000,
     1              4.0000, 5.0000, 7.5000, 10.0000 /
      data c8rev / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1       0.0991, 0.1982, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 / 
      data c8b / 0.4833, 0.4833, 1.2144, 1.6421, 2.181, 2.6087, 2.9122, 3.3399, 3.6434, 3.8787, 4.0711, 
     1       4.3745, 4.6099, 5.0376, 5.3411, 5.7688, 6.0723, 6.5, 6.8035, 7.0389, 7.4666, 7.77 / 
      data c8org / 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 
     1         0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 /

C     Mean Model Coefficients - StrikeSlip
      b0_ssm = -0.078101
      b1_ssm = -0.033923
      b2_ssm = 0.20067
      b3_ssm = 0.149361
      bM_ssm = -0.3004

C     Sigma Model Coefficients - StrikeSlip
      b0_sss = 0.029001
      b1_sss = 0.21112
      b2_sss = 0.53116
      b3_sss = -0.60145
      bM_sss = -0.16960

C     Mean Model Coefficients - Reverse
      b0_rvm = -0.13414
      b1_rvm = 0.022371
      b2_rvm = 0.047363
      b3_rvm = 0.02078
      b4_rvm = 0.1042
      b5_rvm = 0.00716
      b6_rvm = -0.1239
      b7_rvm = 0.069512
      b8_rvm = 0.076094
      bM_rvm = -0.26717

C     Sigma Model Coefficients - Reverse
      b0_rvs = 0.044478
      b1_rvs = -0.021982
      b2_rvs = -0.025129
      b3_rvs = 0.020343
      b4_rvs = -0.022130
      b5_rvs = 0.030626
      bM_rvs = 0.0


C     First interpolate terms for given spectral period. 
C Find the requested spectral period and corresponding coefficients
      nPer = 22
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         c8revT  = c8rev(1)
         c8bT  = c8b(1)
         c8orgT  = c8org(1)

         goto 1011
      elseif (specT .gt. 0.0) then
C Now loop over the spectral period range of the attenuation relationship.
         do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
               count1 = i
               count2 = i+1
               goto 1020 
            endif
         enddo
      endif

C Interpolate the coefficients for the requested spectral period.
 1020       call interp (period(count1),period(count2),c8rev(count1),c8rev(count2),
     +                   specT,c8revT,iflag)
            call interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
            call interp (period(count1),period(count2),c8org(count1),c8org(count2),
     +                   specT,c8orgT,iflag)

 1011 period1 = specT                                                                                                              

C     Compute the mag and distance tapers. 
      taperdist = max(1.0-max(Rrup-40.0,0.0)/30.0,0.0)
      taperMag = min(max(mag-5.5,0.0)/0.8,1.0)
      RyRatio = min(abs(Ry)/(Ruplen/2.0),1.0)

C     Now calculate the parameters for the model given input values.
C     Classification is strike-slip or normal --> Strike-slip coefficients
C                       reverse or oblique    --> Reverse coefficients
      if (mech .le. 0.0) then

         cos2thetaAvg = ( ( (Ry+RupLen/2.0)-2.0*abs(Rx)*atan2((Ry+Ruplen/2.0),abs(Rx)) ) -
     1          ( (Ry-Ruplen/2.0)-2.0*abs(Rx)*atan2((Ry-Ruplen/2.0),abs(Rx)) ) ) /Ruplen

         MDF = b0_ssm + b1_ssm*(max(RyRatio*cos2ThetaAvg,-0.5)) + b2_ssm*(max(RyRatio*cos2ThetaAvg,-0.5))**2.0 +
     1            b3_ssm*(max(RyRatio*cos2ThetaAvg,-0.5))**3.0    
         SDF = b0_sss + b1_sss*(max(RyRatio*cos2ThetaAvg,-0.5)) + b2_sss*(max(RyRatio*cos2ThetaAvg,-0.5))**2.0 +
     1            b3_sss*(max(RyRatio*cos2ThetaAvg,-0.5))**3.0         

         Medadj = MDF*exp(bM_ssm*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag
         SDA =    SDF*exp(bM_sss*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag

      else
         cosdip = cos(dip*3.14159/180.0)
         sindip = sin(dip*3.14159/180.0)
         cos2thetaAvg = ( ( (Ry+RupLen/2.0)-2.0*abs(Rx*sindip)*atan2((Ry+Ruplen/2.0),abs(Rx*sindip)) ) -
     1          ( (Ry-Ruplen/2.0)-2.0*abs(Rx*sindip)*atan2((Ry-Ruplen/2.0),abs(Rx*sindip)) ) ) /Ruplen

         cos2phiAvg = ( ( (Rx*sindip)-2.0*abs(Rx*cosdip)*atan2(Rx*sindip,abs(Rx*cosdip))) -
     1         ( (Rx*sindip-RupWidth)-2.0*abs(Rx*cosdip)*atan2((Rx*sindip-RupWidth),abs(Rx*cosdip))) ) / RupWidth

         MDF = b0_rvm + b1_rvm*cos2phiAvg + b2_rvm*cos2phiAvg**2.0 + b3_rvm*cos2phiAvg**3.0 +
     1           b4_rvm*cos2thetaAvg + b5_rvm*cos2thetaAvg**2.0 + b6_rvm*cos2thetaAvg**3.0 + 
     1           b7_rvm*RyRatio + b8_rvm*RyRatio**2.0

         SDF = b0_rvs+ b1_rvs*cos2phiAvg + b2_rvs*cos2phiAvg**2.0 + b3_rvs*cos2phiAvg**3.0 +
     1           b4_rvs*RyRatio + b5_rvs*RyRatio**2.0
         Medadj = MDF*exp(bM_rvs*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag
         SDA =    SDF*exp(bM_rvs*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag

      endif

C     Check for sigma adjustment values less than 0.0, if so set equal to 0.0
      if (sda .lt. 0.0) SDA=0.0
      sigadj = SDA

      return
      end
