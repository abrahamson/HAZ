      
      subroutine ruptdirct2012_jrb ( specT, distrup, mag, rupLength, width, ftype, 
     &                                theta, rake, Rx, X, Y, az, lnfd, lnfn, lnfp ) 

       implicit none
            
       real rake, Q1Rake, DipWt, StrikeWt, az     
       real specT
       real lnfd, lnfn, lnfp,ftype, width
       real slnDird, dlnDird, slnDirn, dlnDirn, slnDirp, dlnDirp 
       real X, Y, theta, mag, distrup, Rx, ruplength, pi
       integer nper          
                             
       
C***************************************************************************************************************************      
C      REQUIRED FUNCTION INPUT
C      
C      X, Y        = X (for strike slip) and Y (for dip slip) directivity parameter, SSGA97
C      theta, phi  = theta (for strike slip) and phi (for dip slip) directivity parameter, SSGA97. **Must be in DEGREES**
C                    theta = [0,90] and phi = [0,180]
C      M           = moment magnitude
C      Rrup        = closest distance to the fault rupture, km
C      L           = fault length, km
C      W           = down-dip fault width, km
C      Az          = source-site azimuth as defined by PEER NGA. **Must be in DEGREES**
C      Dip         = fault dipping angle.  **Must be POSITIVE and in DEGREES**
C      rake        = rake angle, [-180 180]  **Must be in DEGREES**
C      T           = Period of vibration, sec. Note: directivity effects only observed for T greater than about 0.5 sec
C      
C      FUNCTION OUTPUT
C      
C      fD.FN       = directivity correction in natural log units, FN component
C      fD.FP       = directivity correction in natural log units, FP component
C      fD.RotD50   = directivity correction in natural log units, RotD50 component
C      
C      terms.T_AZ  = azimuth taper term
C      terms.T_M   = magnitude taper term
C      terms.T_CD  = distance taper term
C      terms.fgeom = geometry directivity predictor term
C      terms.C     = period dependent constant coefficients for each component
C***************************************************************************************************************************             

C      Test input parameters for Owens Valley fault
c       rake = 0.
c       ftype = 1
c       mag = 7.5
c       distrup=5.1
c       theta=0.047
c       Rx=5.1
c       ruplength=108.2
c       X=1.         

C      Test input parameters for Sierra Frontal 
c       rake = -90.
c       ftype = 0
c       mag = 7.75
c       distrup=0.6
c       Rx=2.0
c       ruplength=112.5
c       Y=1. 
c       width=22.8  
c       az=90.*pi/180.
       
       pi = acos(-1.0)
       theta=theta*pi/180.
C       az=az*pi/180.     
         
       if (ftype. eq. 0.) then
           StrikeWt=1.
           DipWt=0.
       elseif (ftype. eq. 1. or. ftype. eq. -1.) then
           StrikeWt=0.
           DipWt=1.
       else
           if (abs(rake).gt.90.) then
             Q1Rake=180-abs(rake)
           else
             Q1Rake=abs(rake)
           endif
           DipWt=Q1Rake/90.
           StrikeWt=1-DipWt     
       endif

     
       if (specT.lt.0.1) specT=0.1 ! No Directivity effect for PGA, PGV etc
       if (specT.gt.10.) specT=10. ! Model capped at 10s 
     
       if (ftype. eq. 0.) then 
         call rupdirct_strike2012 ( mag, 
     1         distrup, theta, Rx, ruplength, specT, X,
     2         slnDird, slnDirn, slnDirp )
       
       else 
       
         call rupdirct_dip2012 ( mag, width, az, 
     1         distrup, Rx, ruplength, specT, Y,
     2         dlnDird, dlnDirn, dlnDirp )
       
       endif 
       
       lnfd = StrikeWt * slnDird + DipWt * dlnDird
       lnfn = StrikeWt * slnDirn + DipWt * dlnDirn
       lnfp = StrikeWt * slnDirp + DipWt * dlnDirp
             
       return
       END


C*********************************************************************************
      subroutine rupdirct_strike2012 ( mag, 
     1           distrup, theta, Rx, ruplength, specT, X,
     2           lnDird, lnDirn, lnDirp )      


      implicit none 

C     -- input parameters
      real mag
      real distrup, theta, Rx, ruplength, specT, X    
C     -- output parameters        
      real lnDird, lnDirn, lnDirp 
      

C     -- internal paramaters
      real CDL,TCD,TMW,TAZ,fgeom,fd,fn,fp
      real C,C1,Cn,C1n,Cp,C1p
      real cav_ss(16),c1av_ss(16),cfn_ss(16),c1fn_ss(16),cfp_ss(16),c1fp_ss(16)
      real period(16),Sf,exp1
      real pi, pio2
      integer nper, count1, count2, iflag, i

C     Apply directivity and fault normal/fault parallel effects to hazard.
C       Dirflag        Effect
C       -------        ------
C          0           No directivity
C         11           Directivity average component
C         12           Directivity fault normal component
C         13           Directivity fault parallel component
 
      data period/ 0.1  , 0.15 , 0.2  , 0.25 , 0.3  , 0.4  , 0.5  , 0.75, 
     1             1.   , 1.5  , 2.   , 3.   , 4.   , 5.   , 7.5  , 10. /
      data cav_ss/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.060, 
     1             -0.120,-0.175,-0.210,-0.235,-0.255,-0.275,-0.290,-0.300/
      data c1av_ss/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.038, 
     1             0.075, 0.090, 0.095, 0.099, 0.103, 0.108, 0.112, 0.115/
      data cfn_ss/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.080, 
     1            -0.225,-0.300,-0.325,-0.365,-0.390,-0.410,-0.420,-0.425/
      data c1fn_ss/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.055, 
     1             0.110, 0.135, 0.160, 0.185, 0.205, 0.215, 0.220, 0.225/
      data cfp_ss/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     1             0.015, 0.030, 0.050, 0.070, 0.080, 0.090, 0.100, 0.108/
      data c1fp_ss/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     1             0.000,-0.025,-0.040,-0.045,-0.050,-0.060,-0.070,-0.071/

      pi = acos(-1.0)
      pio2 = pi / 2.
      nper = 16

C     -- Find correct coefficient
      C=0.0
      C1=0.0
      do i=1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1001
         endif
      enddo
      if(specT .gt. period(nper) ) then
            C=cav_ss(nper)
            C1=c1av_ss(nper)
            Cn=cfn_ss(nper)
            C1n=c1fn_ss(nper)
            Cp=cfp_ss(nper)
            C1p=c1fp_ss(nper)
      endif
      goto 1002
1001     call interp (period(count1),period(count2),cav_ss(count1),
     1                cav_ss(count2),specT,C,iflag)
         call interp (period(count1),period(count2),c1av_ss(count1),
     1                c1av_ss(count2),specT,C1,iflag)
         call interp (period(count1),period(count2),cfn_ss(count1),
     1                cfn_ss(count2),specT,Cn,iflag)
         call interp (period(count1),period(count2),c1fn_ss(count1),
     1                c1fn_ss(count2),specT,C1n,iflag)
         call interp (period(count1),period(count2),cfp_ss(count1),
     1                cfp_ss(count2),specT,Cp,iflag)
         call interp (period(count1),period(count2),c1fp_ss(count1),
     1                c1fp_ss(count2),specT,C1p,iflag)   
      
1002  continue
C     -- Calculate Sf
      exp1 = exp(1.)
      Sf = X*ruplength
      if(Sf .lt. exp1) Sf = exp1
      
C     -- Main term
      fgeom = alog(Sf)*(0.5*cos(2*theta)+0.5)
c      write (*,'( 2x,''mag, x, ...''5f10.3)') mag, X, rupLength, sf, theta, fgeom
c      pause

C     -- Distance taper
      CDL = distrup/rupLength
      if(CDL .lt. 0.5) then
         TCD = 1.
      elseif(CDL .ge. .5 .and. CDL .lt. 1.0) then
         TCD=1.-(CDL-0.5)/0.5
      else
         TCD=0.0
      endif

C     -- Magnitude taper
      if (mag .lt. 5.0) then
         TMW=0.0
      elseif ( mag .lt. 6.5 ) then
         TMW=1.-(6.5-mag)/1.5
      else
         TMW=1.
      endif

C     -- Azimuth taper
      TAZ=1.
      
c      print*,'          X,        theta,        C,        fgeom,        TCD,        TMW,        TAZ'
c      print*,X,theta,C,fgeom,TCD,TMW,TAZ      

       
      fd = (C+C1*fgeom)*TCD*TMW*TAZ
      fn = (Cn+C1n*fgeom)*TCD*TMW*TAZ 
      fp = (Cp+C1p*fgeom)*TCD*TMW*TAZ 
      
      lnDird = fd
      lnDirn = fn
      lnDirp = fp
      
      return
      end

C-----------------------------------------------------------------------
      subroutine rupdirct_dip2012 ( mag, W, az,
     1           distrup, Rx, ruplength, specT, Y,
     2           lnDird, lnDirn, lnDirp )

      implicit none


      real mag, W, az
      real distrup, Rx, ruplength, specT, Y    
C     -- output parameters        
      real lnDird, lnDirn, lnDirp 

C     -- internal paramaters
      real CDW,TCD,TMW,TAZ,fgeom,fd,fn,fp,RxW,D
      real C,C1,Cn,C1n,Cp,C1p
      real cav_ds(16),c1av_ds(16),cfn_ds(16),c1fn_ds(16),cfp_ds(16),c1fp_ds(16)
      real period(16),pio2,Djeff
      real z1, pi, dist2, npio2, tpio3
      real strikeX, strikeY, astrike, dx, dy
      integer nper, count1, count2, iflag, i, iCellRupStrike, iCellRupDip
  
               
C     Apply directivity and fault normal/fault parallel effects to hazard.
C       Dirflag        Effect
C       -------        ------
C         11           Directivity average component
C         12           Directivity fault normal component
C         13           Directivity fault parallel component
 
      data period/ 0.1  , 0.15 , 0.2  , 0.25 , 0.3  , 0.4  , 0.5  , 0.75,
     1             1.   , 1.5  , 2.   , 3.   , 4.   , 5.   , 7.5  , 10./
      data cav_ds/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     1             0.000, 0.000, 0.000,-0.033,-0.089,-0.133,-0.160,-0.176/
      data c1av_ds/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     1             0.000, 0.000, 0.034, 0.093, 0.128, 0.150, 0.165, 0.179/
      data cfn_ds/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     1             0.000, 0.000, 0.000,-0.034,-0.092,-0.115,-0.122,-0.125/
      data c1fn_ds/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     1             0.000, 0.000, 0.056, 0.120, 0.142, 0.160, 0.165, 0.170/
      data cfp_ds/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     1             0.000, 0.000, 0.000,-0.034,-0.110,-0.175,-0.195,-0.200/
      data c1fp_ds/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     1             0.000, 0.000, 0.030, 0.080, 0.120, 0.150, 0.170, 0.175/

      pi = acos(-1.0)
      nper = 16
      npio2 = -pi / 2.
      tpio3 = 2. * pi / 3.


C     -- Find correct coefficient
      C=0.0
      C1=0.0
      do i=1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1001
         endif
      enddo
      if(specT .gt. period(nper) ) then
            C=cav_ds(nper)
            C1=c1av_ds(nper)
            Cn=cfn_ds(nper)
            C1n=c1fn_ds(nper)
            Cp=cfp_ds(nper)
            C1p=c1fp_ds(nper)
      endif
      goto 1002
1001     call interp (period(count1),period(count2),cav_ds(count1),
     1                cav_ds(count2),specT,C,iflag)
         call interp (period(count1),period(count2),c1av_ds(count1),
     1                c1av_ds(count2),specT,C1,iflag)
         call interp (period(count1),period(count2),cfn_ds(count1),
     1                cfn_ds(count2),specT,Cn,iflag)
         call interp (period(count1),period(count2),c1fn_ds(count1),
     1                c1fn_ds(count2),specT,C1n,iflag)
         call interp (period(count1),period(count2),cfp_ds(count1),
     1                cfp_ds(count2),specT,Cp,iflag)
         call interp (period(count1),period(count2),c1fp_ds(count1),
     1                c1fp_ds(count2),specT,C1p,iflag)
1002  continue

C     -- Calculate d
      Djeff=Y*W
      if(Djeff .lt. 1.) then
         Djeff = 1.
      else
         Djeff = Y*W
      endif
C     -- Calculate RxW
      RxW=Rx/W
      if(RxW .lt. npio2) then
         RxW = npio2
      else if(RxW .gt. tpio3) then
         RxW = tpio3
      endif
            
C     -- Main term
      fgeom = alog(Djeff)*cos(RxW)

C     -- Distance taper
      CDW=distrup/W
      if(CDW .lt. 1.5) then
         TCD = 1.
      elseif(CDW .ge. 1.5 .and. CDW .lt. 2.0) then
         TCD=1.-(CDW-1.5)/0.5
      else
         TCD=0.0
      endif

C     -- Magnitude taper
      if (mag .lt. 5.0) then
         TMW=0.0
      elseif ( mag .lt. 6.5 ) then
         TMW=1.-(6.5-mag)/1.5
      else
         TMW=1.
      endif

C     -- Azimuth taper
      TAZ=sin(az)*sin(az)  
     
      fd = (C+C1*fgeom)*TCD*TMW*TAZ
      fn = (Cn+C1n*fgeom)*TCD*TMW*TAZ 
      fp = (Cp+C1p*fgeom)*TCD*TMW*TAZ 
      
      lnDird = fd
      lnDirn = fn
      lnDirp = fp

c      print*,'          Y,          C,        fgeom,        TCD,        TMW,        TAZ    '
c      print*,Y,C,fgeom,TCD,TMW,TAZ

      return
      end
