c      subroutine rupdirct_strike (dirflag,iAtten, sigmaY, 
c     1           dirindex, Lginten, ftype,iFlt, diptmp, mag, 
c     2           xRup, yRup, zRup, x0, y0, z0, distrup, fwflag, 
c     3           hwflag,theta, Xcost, YcosP, nSegRup, tapflag,
c     4           nhypo_strike, istrike, specT )

      subroutine rupdirct_strike (period, rrup, mag, rupLen, rupwidth, ftype, tprime,  
     1                            uprime, s, lnDir, sigdirY, dflag, X, Y )

C This subroutine will apply when requested the fault directivity model
C from Sommerville et al., (1997), SRL, Vol. 68, No. 1. The model has
C been modified to include a saturation effect (Abrahamson, 2000).
C This subroutine will work for strike-slip faults.

c      include 'pfrisk.h'

      real c1ss(10), c2ss(10), c1dp(10), c2dp(10)
      real c1fn(10), c2fn(10), c3fn(10), per(10)
      real c1ssT, c2ssT, c1dpT, c2dpT, c1fnT, c2fnT, c3fnT
      real rrup, mag, uprime, tprime, theta, rupLen, s, sigdirY, lnDir
      real pi, pio2, period, yCosP, XCosT, slope1, dirfact, fnfact, rupWidth
      real ftype, X, Y
      integer iflag, dflag, count1, count2

c      real ftype(MAX_FLT), diptmp, mag , xRup(4,MAX_SEG)
c      real yRup(4,MAX_SEG), zRup(4,MAX_SEG), slope1
c      real x0, y0, z0, distrup, theta, Xcost, YcosP
c      real deltax, deltaY, strike1,ruplength1, hypox, hypoy
c      real az1, pi, dist2, pio2, xx, projdist, deltax1, deltay1
c      real specT, c1ssT, c2ssT, c1fnT, c2fnT, c3fnT, period(8)
c      integer nper, count1, count2, iflag
c      integer hwflag, dirflag(MAX_PROB), tapflag(MAX_PROB)
c      integer dirindex(MAX_PROB), fwflag, iAtten, iFlt
c      integer nhypo_strike, istrike

C
C Apply directivity and fault normal/fault parallel effects to hazard.
C       Dirflag        Effect
C       -------        ------
C          0           No directivity
C          1           Directivity average component
C          2           Directivity fault normal component
C          3           Directivity fault parallel component
C          4           Directivity fault normal component w/o average horizontal effect
C          5           Directivity fault parallel component w/o average horizontal effect
 
      data per /0.60, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00, 5.00, 6.00, 10.00 /
      data c1ss /0.000, -0.084, -0.192, -0.344, -0.452, -0.605, -0.713,
     1          -0.797, -0.797, 0.0 /
      data c2ss /0.000, 0.185, 0.423, 0.759, 0.998, 1.333, 1.571, 1.757, 1.757, 0.0 /
      data c1dp /0.000, -0.045, -0.104, -0.186, -0.245, -0.327, -0.386,
     1          -0.431, -0.431, -0.431 /
      data c2dp /0.000, 0.008, 0.178, 0.318, 0.418, 0.559, 0.659, 0.737, 0.737, 0.737 /
      data c1fn /0.0, 0.061, 0.104, 0.164, 0.207, 0.353, 0.456, 0.450, 0.424, 0.424 /
      data c2fn /0.0, -0.0155, -0.0255, -0.0490, -0.0613, -0.1007, 
     1          -0.1282, -0.1269, -0.1223, -0.1223 / 
      data c3fn /0.0, 0.000, 0.000, 0.034, 0.059, 0.093, 0.118, 0.137, 0.152, 0.152 /
      pi = 3.1415926

c      deltax = xRup(2,nSegRup) - xRup(1,1)  
c      deltaY = yRup(2,nSegRup) - yRup(1,1)
c      strike1 = atan2(deltay,deltax)
c      ruplength1 = sqrt ( deltax**2 + deltay**2 )
c      hypox = float(istrike)/(nhypo_strike+1) * 
c     1        deltax + xRup(1,1)
c      hypoy = float(istrike)/(nhypo_strike+1) * 
c     1        deltay + yRup(1,1)
c      deltax1 = x0 - hypox
c      deltay1 = y0 - hypoy 
c      az1 = atan2(deltay1, deltax1 )
c      theta = az1 - strike1
c      if ( theta .lt. -pi ) theta = theta + 2*pi
c      if ( theta .gt. pi ) theta = theta - 2.*pi
c      dist2 = sqrt( deltax1**2 + deltay1**2)

      theta = atan2(tprime, uprime)
      if ( theta .lt. -pi ) theta = theta + 2*pi
      if ( theta .gt. pi ) theta = theta - 2.*pi
      pio2 = 3.1415926 / 2.
      X = s/ruplen
      Y = litd/rupwidth

C Check to see if theta is between -90 to 90 degrees.
c
c      if ( theta .lt. pio2 .and. theta .gt. -pio2 ) then
c         x = 1. - float(istrike)/(nhypo_strike+1)
c         xx = x
c         projdist = cos(theta) * dist2
c
c         if ( projdist .lt. x*ruplength1 ) then
c            x = projdist/ruplength1
c         endif
c      else
c         x = float(istrike)/(nhypo_strike+1)
c         xx = x
c         projdist = -cos(theta) * dist2
c         if ( projdist .lt. x*ruplength1 ) then
c            x = projdist/ruplength1
c         endif
c      endif

C Apply Bay Bridge taper adjustment to the directivity model (Abrahamson, 2000).
c Interpolate the slope factor.

      nper = 10
      do i=1, nper-1
         if (period .ge. per(i) .and. period .le. per(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1001
         endif
      enddo

      if (period .gt. per(nper) ) then
         c1ssT = c1ss(nper)
         c2ssT = c2ss(nper)
         c1fnT = c1fn(nper)
         c2fnT = c2fn(nper)
         c3fnT = c3fn(nper)
         goto 1002
       endif

C Interpolate the coefficients for the requested spectral period.
 1001    call interp (per(count1),per(count2),c2ss(count1),c2ss(count2),
     +                period,c2ssT,iflag)
         call interp (per(count1),per(count2),c1ss(count1),c1ss(count2),
     +                period,c1ssT,iflag)
         call interp (per(count1),per(count2),c1fn(count1),c1fn(count2),
     +                period,c1fnT,iflag)
         call interp (per(count1),per(count2),c2fn(count1),c2fn(count2),
     +                period,c2fnT,iflag)
         call interp (per(count1),per(count2),c3fn(count1),c3fn(count2),
     +                period,c3fnT,iflag)
                                              
c 1002 slope1 = c2ssT * 1.96
 1002 slope1 = c2ssT * 1.88
      if ( x .gt. 0.4 ) x = 0.4
      dirfact = c1ssT + slope1*X*abs(cos(theta))

C Taper Dirfact based on rupdist.

      if (rrup .gt. 30.0) then
         if (rrup .lt. 60.0) then
            dirfact = dirfact-dirfact*(rrup-30.0)*0.0333
         else
            dirfact = 0.0
            endif
      endif

      Xcost = X*abs(cos(theta))
  
C  Only add for cases where M>6.5, taper for m=6 to 6.5
C  and for cases in which it is requested (i.e., not DIRFLAG=4 and 5).
      if (dflag .lt. 4) then
         if (mag .lt. 6.0) then
            LnDir = 0.0
         elseif ( mag .lt. 6.5 ) then
            LnDir = dirfact*(mag-6.0)/0.5
         else
            LnDir = dirfact 
         endif
      endif

C Add Fault Normal Effect if requested or subtract for Fault 
C Parallel effect.
      if (dflag .eq. 2 .and. period .gt. 0.60 ) then
          fnfact = c1fnT+
     1             c2fnT*alog(rrup+1)+
     2             c3fnT*(mag-6)
          if (theta .lt. 0.785398) then
             LnDir = LnDir + (fnfact*cos(2.0*theta)) 
          endif
      elseif (dflag .eq. 4 .and. period .gt. 0.60 ) then
          fnfact = c1fnT+
     1             c2fnT*alog(rrup+1)+
     2             c3fnT*(mag-6)
          if (theta .lt. 0.785398) then
             LnDir = LnDir + (fnfact*cos(2.0*theta)) 
          endif
      elseif (dflag .eq. 3 .and. period .gt. 0.60 ) then
         fnfact = c1fnT+
     1            c2fnT*alog(rrup+1)+
     2            c3fnT*(mag-6)
         if (theta .lt. 0.785398) then
            LnDir = LnDir - (fnfact*cos(2.0*theta)) 
         endif
      elseif (dflag .eq. 5 .and. period .gt. 0.60 ) then
         fnfact = c1fnT+
     1            c2fnT*alog(rrup+1)+
     2            c3fnT*(mag-6)
         if (theta .lt. 0.785398) then
            LnDir = LnDir - (fnfact*cos(2.0*theta)) 
         endif
      endif

C     Reduce sigma based on a reduction of 0.05 ln units at T=3.0sec.
C     Sigma will be reduced (i.e., subtracted) back in main code.
      sigDirY = 0.05*(c2ssT/c2ss(6))



      YcosP = 2.0

      return
      end

C-----------------------------------------------------------------------
      subroutine rupdirct_dip (dirflag,iAtten, sigmaY,dirindex,
     1           lginten, ftype,iFlt, diptmp, mag, xRup, yRup, zRup,
     2           x0, y0, z0, distrup, fwflag, hwflag, theta,
     3           Xcost, YcosP, nSegRup, tapflag, nhypo_dip, idip, specT,
     4           iSegClose )

C This subroutine will apply when requested the fault directivity model
C from Sommerville et al., (1997), SRL, Vol. 68, No. 1. The model has not
C been modified to include a saturation effect like the strike-slip case.

      include 'pfrisk.h'

      real c1ss(8), c2ss(8), c1dp(8),c2dp(8), dirfact, Y
      real c1fn(8), c2fn(8), c3fn(8), fnfact, sigmaY, Lginten
      real ftype(MAX_FLT), diptmp, mag , xRup(4,MAX_SEG)
      real yRup(4,MAX_SEG), zRup(4,MAX_SEG)
      real x0, y0, z0, distrup, theta, Xcost, YcosP
      real deltax, deltaY, strike1, pi, d1a, d3, zeta
      real zz, rr, rr1, phi, L1, d1, d2, unitX, unitY
      real specT, c1dpT, c2dpT, c1fnT, c2fnT, c3fnT, period(8)
      integer nper, count1, count2, iflag
      integer hwflag, dirflag(MAX_PROB), tapflag(MAX_PROB)
      integer dirindex(MAX_PROB), fwflag, iAtten, iFlt
      integer nhypo_dip, idip, iSegClose

C Sommerville et al. (1997) Model coefficient values:

      data period /0.60, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00, 5.00/
      data c1ss /0.000, -0.084, -0.192, -0.344, -0.452, -0.605, -0.713,
     1          -0.797 /
      data c2ss /0.000, 0.185, 0.423, 0.759, 0.998, 1.333, 1.571,1.757/
      data c1dp /0.000, -0.045, -0.104, -0.186, -0.245, -0.327, -0.386,
     1          -0.431 /
      data c2dp /0.000, 0.008, 0.178, 0.318, 0.418, 0.559, 0.659,0.737/
      data c1fn /0.0, 0.061, 0.104, 0.164, 0.207, 0.353, 0.456, 0.450/
      data c2fn /0.0, -0.0155, -0.0255, -0.0490, -0.0613, -0.1007, 
     1          -0.1282, -0.1269/ 
      data c3fn /0.0, 0.000, 0.000, 0.034, 0.059, 0.093, 0.118, 0.137/

      pi = 3.1415926
               
c The directivity indexes given below are no longer used in the
c hazard code. The coefficients for the model are now determined
c based in the specT value entered in the input file. For cases
c between the defined spectral period range of 0.6-5.0 seconds 
c the model will be interpolated for periods other than those listed
c below. For any hazard run which has a period greater than 5.0 seconds
c the program will terminate. 

C Sommerville et al. (1997) Model coefficient values:
C              DirIndex      Period (sec)
C              --------      ------------
C                 1              0.60
C                 2              0.75
C                 3              1.00
C                 4              1.50
C                 5              2.00
C                 6              3.00
C                 7              4.00
C                 8              5.00
C
C Apply directivity and fault normal/fault parallel effects to hazard.
C       Dirflag        Effect
C       -------        ------
C          0           No directivity
C          1           Directivity average component
C          2           Directivity fault normal component
C          3           Directivity fault parallel component
C          4           Directivity fault normal component w/o average horizontal effect
C          5           Directivity fault parallel component w/o average horizontal effect


C Only apply the directivity factor to sites located on the hanging
C wall or foot wall. 
      if ( fwflag .eq. 1 .or. hwflag .eq. 1 ) then

C Geometry problem of finding the angle phi.

           y = 1. - (float(idip)/(nhypo_dip+1))/2.
           deltax = xRup(2,iSegClose) - xRup(1,iSegClose)  
           deltaY = yRup(2,iSegClose) - yRup(1,iSegClose)
           deltaz =  zRup(3,iSegClose) - zRup(1,iSegClose)  
           L1 = sqrt( (yRup(1,iSegClose)-y0)**2 + 
     1          (xRup(1,iSegClose)-x0)**2 ) 
           zz = - (1.0 - y) * deltaz + zRup(3,iSegClose) - zRup(1,iSegClose)

C Set up the unit vector along the fault strike.
           unitX = deltaX/sqrt(deltaX**2 + deltaY**2)
           unitY = deltaY/sqrt(deltaX**2 + deltaY**2)
           d1a = (xRup(1,iSegClose)-x0)*unitX + 
     1             (yRup(1,iSegClose)-y0)*unitY

           if (abs(d1a) .gt. abs(L1)) then
              d1 = 0.0
           else
              if (fwflag .eq. 1) then
                 d1 = -1.0*sqrt( L1*L1 - d1a*d1a)
              else
                 d1 = sqrt( L1*L1 - d1a*d1a)
              endif
           endif

           d2 = zz/abs(tan(diptmp*3.14159/180.0))
           d3 = d2 - d1

           zeta = atan2(d3,zz+zRup(1,1))
           phi = 1.5708 - diptmp*3.14159/180.0 - zeta

C Now compute the directivity factor.
c Interpolate the slope factor.

      nper = 8

      do i=1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1001
         endif
      enddo

      if (specT .gt. period(nper) ) then
         c2dpT = c2dp(nper)
         c1dpT = c1dp(nper)
         c1fnT = c1fn(nper)
         c2fnT = c2fn(nper)
         c3fnT = c3fn(nper)
         goto 1002
       endif

C Interpolate the coefficients for the requested spectral period.
 1001    call interp (period(count1),period(count2),c2dp(count1),c2dp(count2),
     +                specT,c2dpT,iflag)
         call interp (period(count1),period(count2),c1dp(count1),c1dp(count2),
     +                specT,c1dpT,iflag)
         call interp (period(count1),period(count2),c1fn(count1),c1fn(count2),
     +                specT,c1fnT,iflag)
         call interp (period(count1),period(count2),c2fn(count1),c2fn(count2),
     +                specT,c2fnT,iflag)
         call interp (period(count1),period(count2),c3fn(count1),c3fn(count2),
     +                specT,c3fnT,iflag)
                                              
 1002    dirfact = c1dpT + 
     1               c2dpT*Y*abs(cos(phi))

         Yphi = Y*abs(cos(phi))

C  only add for cases where M>6.5
C  and for cases in which it is requested (i.e., not DIRFLAG=4 and 5).
        gmtest = Lginten
        if (dirflag(iAtten) .lt. 4) then
           if (mag .lt. 6.5) then
               Lginten = Lginten
           else
               Lginten = Lginten + dirfact 
           endif
        endif

c   Add Fault Normal Effect if requested or subtract for Fault 
c   Parallel effect.
           if (dirflag(iAtten) .eq. 2 .and. 
     1                 specT .gt. 0.60 ) then
               fnfact = c1fnT +
     1                  c2fnT*alog(distrup+1) +
     2                  c3fnT*(mag-6)
               if (phi .lt. 0.785398) then
                   Lginten = Lginten + (fnfact*cos(2.0*phi)) 
               endif
           elseif (dirflag(iAtten) .eq. 4 .and. 
     1                 specT .gt. 0.60 ) then
               fnfact = c1fnT +
     1                  c2fnT*alog(distrup+1) +
     2                  c3fnT*(mag-6)
               if (phi .lt. 0.785398) then
                   Lginten = Lginten + (fnfact*cos(2.0*phi)) 
               endif
           elseif (dirflag(iAtten) .eq. 3 .and. 
     1                     specT .gt. 0.60 ) then
               fnfact = c1fnT+
     1                  c2fnT*alog(distrup+1)+
     2                  c3fnT*(mag-6)
               if (phi .lt. 0.785398) then
                   Lginten = Lginten - (fnfact*cos(2.0*phi)) 
               endif
           elseif (dirflag(iAtten) .eq. 5 .and. 
     1                     specT .gt. 0.60 ) then
               fnfact = c1fnT+
     1                  c2fnT*alog(distrup+1)+
     2                  c3fnT*(mag-6)
               if (phi .lt. 0.785398) then
                   Lginten = Lginten - (fnfact*cos(2.0*phi)) 
               endif
           endif
        endif

c       write (*,*) (phi*180)/3.14159, dirfact, Y, fnfact,
c     1      exp(gmtest),exp(lginten)
    
c        pause

C Set the strike slip factor to 2 for deaggregation bins.
       Xcost = 2.0

       return
       end
       
C-----------------------------------------------------------
      Subroutine sigmaP_strike ( sigmaY, dirindex, iAtten, specT ) 

      real c1ss(8), c2ss(8), c1dp(8),c2dp(8)
      real c1fn(8), c2fn(8), c3fn(8),sigmaY
      real c2ssT, period(8), specT
      integer nper, count1, count2
      integer iAtten, dirindex(1)
      
C Sommerville et al. (1997) Model coefficient values:

      data period /0.60, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00, 5.00/
      data c1ss /0.000, -0.084, -0.192, -0.344, -0.452, -0.605, -0.713,
     1          -0.797 /
      data c2ss /0.000, 0.185, 0.423, 0.759, 0.998, 1.333, 1.571,1.757/
      data c1dp /0.000, -0.045, -0.104, -0.186, -0.245, -0.327, -0.386,
     1          -0.431 /
      data c2dp /0.000, 0.008, 0.178, 0.318, 0.418, 0.559, 0.659,0.737/
      data c1fn /0.0, 0.061, 0.104, 0.164, 0.207, 0.353, 0.456, 0.450/
      data c2fn /0.0, -0.0155, -0.0255, -0.0490, -0.0613, -0.1007, 
     1          -0.1282, -0.1269/ 
      data c3fn /0.0, 0.000, 0.000, 0.034, 0.059, 0.093, 0.118, 0.137/

C Directivity Model takes 0.05 ln units off of the sigma 
c (at T=3 sec)

c Interpolate the c2ss term.

      nper = 8

      do i=1, nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1001
         endif
      enddo

      if (specT .gt. period(nper) ) then
         c2ssT = c2ss(nper)
         goto 1002
       endif

C Interpolate the coefficients for the requested spectral period.
 1001 call interp (period(count1),period(count2),c2ss(count1),c2ss(count2),
     +                specT,c2ssT,iflag)
                                              
 1002 sigmaY = sigmaY - 0.05 * (c2ssT/c2ss(6))

      return
      end


C---------------------------------------------------------------------
      subroutine NGArupdirct (jcalc1, specT, lgInten,
     1                sigmaY, distRup, distJB, distHypo, mag, attenName1,
     2                period1, hypodepth, fd, vrup )

C This subroutine will apply the NGA model specific directity factors 
C to a given ground motion based on the Spudich and Chiou (2007) paper. 
C The subroutine will pass the jcalc number to determine which NGA 
C attenuation relationship is being modified. 

      include 'pfrisk.h'

      real period(8), sigAS(8), sigBA(8), sigCB(8), sigCY(8)
      real sigAST, sigBAT, sigCBT, sigCYT, pi
      real a0(4), a1(4), a2(4), b0(4), b1(4), b2(4)
      real specT, period1, mag, lgInten, sigmaY, distRup, DistJB, distHypo
      real hypodepth, hypodepth1, vrup
      real IDP, D, C, s, Rri, fd, fm, fr, aterm, bterm, ctildaprime
      integer NGAmodel, jcalc1, nPer, count1, count2
      character*80 attenName1

c The directivity indexes given below are no longer used in the
c hazard code. The coefficients for the model are now determined
c based in the specT value entered in the input file. For cases
c between the defined spectral period range of 0.6-5.0 seconds 
c the model will be interpolated for periods other than those listed
c below. For any hazard run which has a period greater than 5.0 seconds
c the program will terminate. 
C
C Apply directivity and fault normal/fault parallel effects to hazard.
C       Dirflag        Effect
C       -------        ------
C          0           No directivity
C          1           Directivity average component
C          2           Directivity fault normal component
C          3           Directivity fault parallel component
C          4           Directivity fault normal component w/o average horizontal effect
C          5           Directivity fault parallel component w/o average horizontal effect
C          6           NGA Model Specific Spudich and Chiou Model 
 
      data period / 1.00, 1.50, 2.00, 3.00, 4.00, 5.00, 7.50, 10.00 /

      data a0 / 0.06368428, 0.1461653, 0.09952188, 0.1516156 /
      data a1 / -0.01698134, -0.05711179, -0.03828315, -0.08022648 /
      data a2 / -0.00697153, -0.004596812, -0.006164576, 0.0007366842 /

      data b0 / -0.3911646, -0.07668303, -0.07454529, -0.1169328 /
      data b1 / 0.04640031, 0.07837993, 0.05155013, 0.06438914 /
      data b2 / 0.0, -0.002863621, 0.0, -0.001839794 /

      data sigAS / 0.00, 0.00,  0.00,  -0.015, -0.027, -0.041, -0.076, -0.179 / 
      data sigBA / 0.00, 0.00, -0.001, -0.021, -0.041, -0.063, -0.106, -0.213 / 
      data sigCB / 0.00, 0.00,  0.00,  -0.017, -0.031, -0.054, -0.092, -0.201 / 
      data sigCY / 0.00, 0.00,  0.00,  -0.010, -0.018, -0.029, -0.048, -0.122 / 

      pi = 3.1415926
      nPer = 8
      NGAmodel = 0

C     Check for given attenuation model is one of the four specific NGA models.
      if (jcalc1 .eq. 783) then
          NGAmodel = 1      
          attenname1 = 'A&S(072007 Model)-Hor,Est. Vs30m, NGA Direct'
      elseif (jcalc1 .eq. 784) then
          NGAmodel = 1      
          attenname1 = 'A&S(072007 Model)-Hor,Meas. Vs30m, NGA Direct'
      elseif (jcalc1 .eq. 785) then
          NGAmodel = 1      
          attenname1 = 'A&S(072007 Model)-Hor,No Soil Depth, Est. Vs30m, NGA Direct'
      elseif (jcalc1 .eq. 786) then
          NGAmodel = 1      
          attenname1 = 'A&S(072007 Model)-Hor,No Soil Depth, Meas. Vs30m, NGA Direct'
      elseif (jcalc1 .eq. 921 ) then          
          NGAmodel = 2
          attenname1 = 'Boore&Atkinson(07-07)-Hor, NGA Direct'
      elseif (jcalc1 .eq. 835 ) then          
          NGAmodel = 3
          attenname1 = 'Campbell&Bozorgnia(05-07)-Horizontal, NGA Direct'
      elseif (jcalc1 .eq. 794 ) then          
          NGAmodel = 4
          attenname1 = 'Chiou&Youngs(08-07)-Hor, NGA Direct'
      else
          write (*,*) 'NGA Directity Model only allowed for'
          write (*,*) 'NGA Ground Motion Models!!!'
          write (*,*) 'Current Model in code being called: '
          write (*,*) attenName1
          write (*,*) '*** Check input file ***'
          stop 99                   
      endif

C   Check for the need for interpolation of the sigma adjustments based 
c   on the requested spectral period.
      do i=1,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'NGA model specific directity model is not'
      write (*,*) 'defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call interp (period(count1),period(count2),sigAS(count1),sigAS(count2),
     +                   specT,sigAST,iflag)
            call interp (period(count1),period(count2),sigBA(count1),sigBA(count2),
     +                   specT,sigBAT,iflag)
            call interp (period(count1),period(count2),sigCB(count1),sigCB(count2),
     +                   specT,sigCBT,iflag)
            call interp (period(count1),period(count2),sigCY(count1),sigCY(count2),
     +                   specT,sigCYT,iflag)
       
C Compute the distance taper function.
      if (distrup .gt. 70.0) then
         fr = 0.0
      elseif (distrup .le. 40.0) then
         fr = 1.0
      else
         fr = 1.0 - (distrup-40.0)/30.0
      endif
       
C Compute the magnitude taper function.
      if (mag .gt. 6.4) then
         fm = 1.0
      elseif (mag .le. 6.0) then
         fm = 0.0
      else
         fm = (mag-6)/0.4
      endif

C Compute the a and b terms.        
       aterm = a0(NGAmodel) + a1(NGAmodel)*specT + a2(NGAmodel)*specT*specT
       if (aterm .gt. 0.0) then
           aterm = 0.0
       endif
       bterm = b0(NGAmodel) + b1(NGAmodel)*specT + b2(NGAmodel)*specT*specT
       if (bterm .lt. 0.0) then
           bterm = 0.0
       endif

C Compute distance term s (horizontal distance along fault from hypo to closest point on 
C                          fault plane to station)
       s = 1.0

C Compute distance term D (distance along fault from hypo to closest point on fault plane to station)
       D = 1.0

C Compute ctildaprime.
       if (D .eq. 0.0) then
          ctidlaprime = vrup
       else
          ctildaprime = (1.0/vrup) + ((distrup - distHypo)/D)
          ctildaprime = 1.0 / ctildaprime
       endif
       
C Compute C term. 
       if (ctildaprime .gt. 2.0) then
           C = (2.0 - 0.8)/1.2
       else
           C = (ctildaprime - 0.8)/1.2
       endif
         
C Compute Radiation Pattern.  
       Rri = 1.0

C First compute the angle Phi for multiple segment faults. 

       r = disthypo
       dh = depthhypo
       Rcap = sqrt(disthypo*disthypo-dh*dh)
       uprime = us - uh
       zprime = -dh
       d0 = depthtop
       tprime = -ws/sin(dip*pi/180.0) - (-d0/sin(dip*pi/180.0)+vh)*cos(dip*pi/180.0)

C Set minimum Radiation Pattern Term
       if (Rri .le. 0.2) then
          Rri = 0.2
       endif

C Compute IDP term.
       if (s .lt. 1.0) then
           s = 1.0
       endif 
       if (C*alog(s) .gt. 4.0) then
          IDP = Rri*4.0
       else
          IDP = Rri*C*alog(s)
       endif 

C Compute the directivity adjustment. 
       fd = fr*fm*(aterm+bterm*IDP)

C Adjust the sigma value.
       if (NGAmodel .eq. 1) then
          sigmaY = sigmaY + sigAST
       elseif (NGAModel .eq. 2) then
          sigmaY = sigmaY + sigBAT
       elseif (NGAModel .eq. 3) then
          sigmaY = sigmaY + sigCBT
       elseif (NGAModel .eq. 4) then
          sigmaY = sigmaY + sigCYT
       endif

c       write (*,*) 'Coming out of NGA Rupt'

       return 
       end
          
          
