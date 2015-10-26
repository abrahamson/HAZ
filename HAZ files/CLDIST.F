c closestDist.f

      subroutine CalcDist (sourceType, pscorflag, hypoDepth, RupWidth, RupLen, 
     1             r_horiz, seismoDepth, nFltGrid, n1AS, iLocX, iLocY, n2AS,
     2             fltgrid_x, fltgrid_y, fltgrid_z, fltgrid_x1, fltgrid_y1, 
     3             fltgrid_z1, fltgrid_x2, fltgrid_y2, fltgrid_x3, fltgrid_y3,
     4             fltgrid_x4, fltgrid_y4, fltgrid_z4, fltgrid_Rrup, fltgrid_Rjb,
     5             distJB, distRup, ZTOR, distSeismo, distepi, disthypo, HWFlag,
     6             dipavg, n1, n2, Rx, Ry, Ry0, icellRupstrike, icellRupdip, dip, iFltWidth, iFlt)
           
c     This subroutine computes the distance parameters for the GMPEs
c     (e.g. Rjb, Rrup, ZTOR, Rx, Ry, Ry0)

      implicit none
      include 'pfrisk.h'

c     declarations passed in 
      integer sourceType, pscorflag, nfltGrid(2), n1AS(MAXFLT_AS), iLocX, 
     1        iLocY, n2AS(MAXFLT_AS), iFltWidth, iFlt
      real hypoDepth, RupWidth, RupLen, r_horiz, seismoDepth,
     1     fltgrid_x(MAXFLT_DD,MAXFLT_AS), fltgrid_y(MAXFLT_DD,MAXFLT_AS),
     2     fltgrid_z(MAXFLT_DD,MAXFLT_AS), fltgrid_x1(MAXFLT_DD,MAXFLT_AS), 
     3     fltgrid_y1(MAXFLT_DD,MAXFLT_AS), fltgrid_z1(MAXFLT_DD,MAXFLT_AS), 
     4     fltgrid_x2(MAXFLT_DD,MAXFLT_AS), fltgrid_y2(MAXFLT_DD,MAXFLT_AS), 
     5     fltgrid_x3(MAXFLT_DD,MAXFLT_AS), fltgrid_y3(MAXFLT_DD,MAXFLT_AS), 
     6     fltgrid_x4(MAXFLT_DD,MAXFLT_AS), fltgrid_y4(MAXFLT_DD,MAXFLT_AS), 
     7     fltgrid_z4(MAXFLT_DD,MAXFLT_AS), fltgrid_Rrup(MAXFLT_DD,MAXFLT_AS),
     8     fltgrid_Rjb(MAXFLT_DD,MAXFLT_AS)
       real dip(MAX_FLT,MAX_WIDTH, MAX_SEG)

c     declarations passed out
      real distJB, distRup, ZTOR, distSeismo, distepi, disthypo,
     1     dipavg, Rx, Ry, Ry0
      integer HWFlag, n1, n2, icellRupstrike, icellRupdip  

c     declarations only used within subroutine
      integer i, j
      real hypo1, distRx, rupLength2, r_horiz1 
      
c     Initialize closest distance      
      distRup = 1.0e30
      distJB = 1.0e30
      distSeismo = 1.0e30
      disthypo = 1.0e30
      distepi = 1.0e30
      distRx = 1.0e30
      
c     Compute distances for areal sources      
        if ( sourceType .eq. 2 .or. sourceType .eq. 3 .or. sourceType .eq. 4 ) then
 
c         Approximate correction for using point source with extended source model
C         added limit on only shallow areal sources (Hypo>30.0km)
          if ( psCorFlag .eq. 1  .and. hypodepth .le. 30.0) then
            hypo1=hypoDepth-rupWidth/2.   
          if (hypo1 .lt. 0.0) hypo1=0.0
          rupLength2 = rupLen/2.
          if ( r_horiz .gt. rupLength2 ) then
            r_horiz1 = r_horiz - 0.625*rupLength2*(1.-1./(1.5*(r_horiz/rupLength2)**(1.5)+1))
          else
            r_horiz1 = r_horiz/1.57
          endif
          else
            r_horiz1 = r_horiz
            hypo1 = hypoDepth         
          endif
          distJB = r_horiz1       
          distRup = sqrt( r_horiz1**2 + hypo1**2 )                  
          ZTOR = hypo1
          if ( hypo1 .lt. seismoDepth) then
            distSeismo = sqrt(distJB**2 + seismoDepth**2)
          else
            distSeismo = distRup
          endif 
       
c         Set other distance values equal to corresponding distance measures.
          distepi = distJB
          distRx = distJB
          disthypo = distRup
C         Set HWflag = 0 for areal or grid sources.
          HWFlag = 0
C         Set Dipaverage equal to 90 for gridded sources.
c          dipavg = 3.14159/2.0      
          dipavg = dip(iFlt,ifltWidth,1) * 3.1415926/180.
          return

c     Otherwise use fault sources (i.e., Sourcetype = 1, 5, or 6)
        else 
     
c         n1 is the last cell that makes up the rupture plane down dip
c         n2 is the last cell that makes up the rupture plane along strike       
          n1 = nfltgrid(1)-n1AS(iLocX)+iLocY
          n2 = n2AS(iLocX)
       
c         Calculate Rx, Ry, Ry0, dipavg, and HWFlag
c         Average Strike method
c          call AveStrike (fltgrid_x, fltgrid_y, fltgrid_z, fltgrid_x1,
c     1                    fltgrid_y1, fltgrid_x2, fltgrid_y2, fltgrid_x4,
c     2                    fltgrid_y4, rupLen, distRx, iLocX, iLocY, n2, 
c     3                    n1, Rx, Ry, Ry0, HWFlag, dipavg)

c         OR Global Coordinate System 2 method
          call GC2 (iLocX, iLocY, n2, n1, fltgrid_x1, fltgrid_y1, fltgrid_z1, 
     1              fltgrid_x2, fltgrid_y2, fltgrid_x3, fltgrid_y3,
     2              fltgrid_x4, fltgrid_y4, fltgrid_z4, Rx, Ry, Ry0, 
     3              HWFlag, dipavg)          

c         Calculate ZTOR
          ZTOR = fltgrid_z1(iLocY,iLocX)             

c         Calculate Rrup and Rjb
C         Keep track of the fault grid cell for the closest rupture distance
C         for this rupture area since it is needed for the NGA directivity models.     
     
          do j=iLocX,n2
            do i=iLocY,n1
              if (distRup .gt. fltgrid_Rrup(i,j)) then
                distRup = fltgrid_rRup(i,j)
                icellRupstrike = j
                icellRupdip = i
              endif

              if (distJB .gt. fltgrid_Rjb(i,j)) then
                distJB = fltgrid_Rjb(i,j)
              endif

C         Compute the DistSeismo value based on min depth being equal to seismoDepth parameter.
              if (fltgrid_z(i,j) .gt. seismoDepth) then                 
                if (distSeismo .gt. fltgrid_Rrup(i,j)) then
                  distSeismo = fltgrid_Rrup(i,j)
                endif
              endif 
            enddo
          enddo

C         Set the JB and Hypo distances for faults.
          distepi = fltgrid_Rjb(iLocY,iLocX)
          disthypo = fltgrid_Rrup(iLocY,iLocX)
      
        endif

      return
      end

c ----------------------------------------------------------------------

      Subroutine SetnRupLoc ( n1, n2, nHypoX, pHypoX, nHypoXStep, 
     1                        nHypoZ, pHypoZ, nHypoZstep ) 

      integer n1, n2, nHypoX, nHypoXstep, nHypoZ, nHypoZstep
      real pHypoX, pHypoZ

C     First set up the number of hypocenter locations for a given fault rupture area
C     If there are less than 10 cells in either along strike or along dip direction
C     just use each cell. Otherwise take 10 locations along strike and dip

      if (n2 .lt. 10) then
         nHypoX = n2
         phypoX = real(1.0/nHypoX)
         nHypoXstep = 1
      else
         nHypoX = int(n2/10)*10
         phypoX = real(1.0/10.0)
         nHypoXstep = int(n2/10)
      endif

c     Compute the step sizes down dip 
      if (n1 .lt. 10) then
         nHypoZ = n1
         phypoZ = real(1.0/nHypoZ)
         nHypoZstep = 1
      else
         nHypoZ = int(n1/10)*10
         phypoZ = real(1.0/10.0)
         nHypoZstep = int(n1/10)
      endif
                
      return
      end

c -------------------------------------------------------------------
      subroutine DetermDist (hAScell, hDDcell, icellRupStrike, icellRupdip, 
     1                      fltgrid_x, fltgrid_y, fltgrid_z, n2, n1, dipavg,
     2                      iLocY, iLocX, rupLen, RupWidth, x0, y0, z0,
     3                      edist, hdist, slit, azp1p2, step, dlit, phiang, FltGrid_rRup,
     4                      s2site, Rx, astrike )

      include 'pfrisk.h'
      include 'declare1.h'
      integer hAScell, hDDcell, icellRupStrike, icellRupDip
      integer n2, n1, iLocx, iLocY
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS),
     1      fltGrid_y(MAXFLT_DD,MAXFLT_AS), dx, dy, dz, edist, hdist
      real dist, step, phiang, dlit, strikeX, strikeY, astrike
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), Rx

C     Variable: hDDcell, hAScell --> Location of hypocenter locations down-dip and along strike.
C               iLocY, iLocX --> number of rupture areas to complete fill fault plane.
C               n2, n1 --> number of cells down-dip and along strike.
C               iCellRupdip, iCellRupStirke --> cell location down-dip and along strike for closest point to station.

C      Compute average strike for given fault plane.
       strikeX = fltgrid_x(1,n2) - fltgrid_x(1,1)
       strikeY = fltgrid_y(1,n2) - fltgrid_y(1,1)
       if (strikeX .eq. 0.0) then
          astrike = 0.0
       else
          astrike = atan2(strikeX,strikeY)
       endif

C     Compute distance along strike between closest cell on fault plane and hypocenter location. 
      dx = (fltgrid_x(hDDcell,hAScell) - fltgrid_x(hDDcell,icellRupStrike))
      dy = (fltgrid_y(hDDcell,hAScell) - fltgrid_y(hDDcell,icellRupStrike))
      slit = sqrt( dx**2.0 + dy**2.0)

C     Compute the downdip distance dlit.
      dlit = real(step*(abs(icellRupDip-hDDcell)))

C     Compute the angle Phi between Hypocenter and Station location. 
      phiang = (180.0/3.14159)*atan2(Rx,dlit) - (90.0 - dipavg)

C     Compute Hypocentral and Epicentral Distances.
      dx = fltGrid_x(hDDcell,hAScell) - x0
      dy = fltGrid_y(hDDcell,hAScell) - y0
      dz = fltGrid_z(hDDcell,hAScell) - z0
      hdist = sqrt ( dx*dx + dy*dy + dz*dz)
      edist = sqrt ( dx*dx + dy*dy )

C     Compute azimuth between epicenter and station location.
      dx = x0 - fltGrid_x(hDDcell,hAScell) 
      dy = y0 - fltGrid_y(hDDcell,hAScell) 
      azp1p2 =  atan2(dx,dy)*(180.0/3.14159) - astrike*180/3.14159

C     Compute azimuth between closest point and site.
      dx = x0 - fltGrid_x(iCellRupdip,iCellRupstrike) 
      dy = y0 - fltGrid_y(iCellRupdip,iCellRupstrike)
      s2site =  atan2(dx,dy)*(180.0/3.14159) - astrike*180/3.14159


C     Check to see if the slit value is greater than the limited c1 value. If so recompute values. 
C     Case where hypocenter is further down strike than closest point cell.
      c1temp = 0.0
      if (hAScell .gt. icellRupstrike) then
         do icell=icellrupstrike, hAScell-1, 1
            dx = (fltgrid_x(icellRupdip,icell) - fltgrid_x(icellRupdip,icell+1))
            dy = (fltgrid_y(icellRupdip,icell) - fltgrid_y(icellRupdip,icell+1))
            dist = sqrt( dx**2.0 + dy**2.0)
            c1temp = c1temp + dist
         enddo
C     Case where closest point cell is further down strike than hypocenter.
      else
         do icell=hAScell,icellrupstrike-1, 1
            dx = (fltgrid_x(icellRupdip,icell) - fltgrid_x(icellRupdip,icell+1))
            dy = (fltgrid_y(icellRupdip,icell) - fltgrid_y(icellRupdip,icell+1))
            dist = sqrt( dx**2.0 + dy**2.0)
            c1temp = c1temp + dist
         enddo

      endif
      
      return
      end

c -------------------------------------------------------------------
      subroutine ngaRDirmodel (iHypoX, iHypoZ, icellRupStrike, icellRupdip, 
     1                      fltgrid_x, fltgrid_y, fltgrid_z, n2, n1, dipavg,
     2                      iLocY, iLocX, fltgrid_rrup, Rx, HWFlag, flttype,
     3                      dflag, lnDir, period, mag, sigDirY, rupLen, RupWidth )

      include 'pfrisk.h'
      include 'declare1.h'
c      integer iHypoX, iHypoZ, icellRupStrike, icellRupDip, HWFlag, dflag
      integer iHypoX, iHypoZ, icellRupStrike, icellRupDip, dflag
      integer n2, n1, iLocx, iLocY
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS),
     1      fltGrid_y(MAXFLT_DD,MAXFLT_AS), dipavg, fltgrid_rrup(MAXFLT_DD,MAXFLT_AS),
     2      Rx, dprod1, distcorn, ctildap, BigD, BigS, h, litd, velrat, BigC,
     3      Ts, Us, Th, Uh, Radri, IDP, flttype, BigRu, BigRt, radterm1, 
     2      radterm2, radterm3, radterm4, radterm5, radterm6, zhyp, lnDir,
c     3      period, mag, sigDirY, rupLen, RupWidth
     3      period, sigDirY, Y, X


C     Hypocenter is in (iHypoX, iHypoZ) location for given (iLocX,iLocY) rupture area.
C     Compute the distance down the fault plane for hypocenter
cnjg      h = abs((fltgrid_z(iLocY,iHypoX) - fltgrid_z(iHypoY,iHypoX)) / sin(dipavg) )
      h = abs((fltgrid_z(iLocY,iHypoX) - fltgrid_z(iHypoZ,iHypoX)) / sin(dipavg) )

C     Compute distances from hypocenter cell to closest distance cell.
      BigD = sqrt( (fltgrid_x(iHypoZ,iHypoX) - fltgrid_x(icellRupdip,icellRupStrike))**2.0 +
     1             (fltgrid_y(iHypoZ,iHypoX) - fltgrid_y(icellRupdip,icellRupStrike))**2.0 + 
     2             (fltgrid_z(iHypoZ,iHypoX) - fltgrid_z(icellRupdip,icellRupStrike))**2.0 )        

      litd = abs((fltgrid_z(iHypoZ,icellRupStrike) - fltgrid_z(icellRupDip,icellRupStrike)) /sin(dipavg) )

      if (litd .gt. BigD) then
         s = sqrt (litd*litd - BigD*BigD)
      else
         s = sqrt (BigD*BigD - litd*litd)
      endif
C     Compute c tilda prime (Eq. 4 of Spudich and Chiou).
      velrat = 0.80
      rhypo = fltgrid_rrup(iHypoZ,iHypoX)
      rrup  = fltgrid_rrup(icellRupdip,icellRupStrike)

      if (BigD .eq. 0.0) then
         ctildap = velrat
      else
         ctildap = ((1.0/velrat) - ((rhypo-rrup)/BigD))
         ctildap = 1.0/ctildap
      endif

C     Now compute BigC term. 
      if (ctildap .lt. 2.45) then 
         BigC = (ctildap - 0.80) / (2.45 - 0.80)
      else
         BigC = (2.45 - 0.80) / (2.45 - 0.80)
      endif   
C     Now compute BigS term.
      if (s .gt. 0.0 .or. h .gt. 0.0) then
         BigS = alog( min(75.0,max(s,h)) )
      else
         BigS = 0.0
      endif
C     Compute the Radiation term.
C     Set Ts (perpendicular distance from station to top of fault trace)
C         distance equal to Rx distance. Note that Rx value is based on 
C         average strike for given fault rupture segment. Hanging wall sites 
C         are given positive values and footwall sites are negative.
      if (HWFlag .eq. 1) then
         Ts = Rx
      else
         Ts = -1.0*Rx
      endif
C     Now compute the Us value which is the distance along the strike at
C         the top of the fault rupture area to the closest perpendicular
C         distance to the site. To be consistent with the calculation of Rx
C         the calculation of Us is based on the average strike for a given
C         rupture area. 
      Us = 0.0
C     Compute distance from first cell of rupture area and site location
C         which is assumed to be at (0.0, 0.0, 0.0).
      Distcorn = sqrt ( (fltgrid_x(iLocY,iLocX)*fltgrid_x(iLocY,iLocX)) +
     1                  (fltgrid_y(iLocY,iLocX)*fltgrid_y(iLocY,iLocX)) )
      if (distcorn*distcorn .gt. Rx*Rx) then
         Us = sqrt ( Distcorn*Distcorn - Rx*Rx )
      else
         Us = 0.0
      endif
C     Now check for sign on Us value. Take dot product and if positive then
C         Us --> positive and if negative then Us --> negative      
      dprod1 = ( (fltgrid_x(iLocY,n2)-fltgrid_x(iLocY,iLocX)) * (0.0 - fltgrid_x(iLocY,iLocX) ) +
     1          (fltgrid_y(iLocY,n2)-fltgrid_y(iLocY,iLocX)) * (0.0 - fltgrid_y(iLocY,iLocX) ) ) 
      if (dprod1 .lt. 0.0) then
         Us = -1.0*Us
      endif

C     Now compute the Th and Uh distances for the given rupture area and hypocenter location.
      Uh = 0.0
      do i=iLocX,iHypoX
         Uh = Uh + sqrt( (fltgrid_x(iLocY,i)*fltgrid_x(iLocY,i)) +      
     1                   (fltgrid_y(iLocY,i)*fltgrid_y(iLocY,i)) )  
      enddo

C     Compute Th distance checking for near vertical dipping faults.       
C         If Dip > 80 then assume vertical dip for depth calculation.
      Th = 0.0
      if (abs(dipavg) .gt. 1.396) then
cnjg         Th = (fltgrid_z(iLocY,iHypoX) - fltgrid_z(iHypoY,iHypoX) )
         Th = (fltgrid_z(iLocY,iHypoX) - fltgrid_z(iHypoZ,iHypoX) )
      else
cnjg         Th = (fltgrid_z(iLocY,iHypoX) - fltgrid_z(iHypoY,iHypoX) ) / tan(dipavg)
         Th = (fltgrid_z(iLocY,iHypoX) - fltgrid_z(iHypoZ,iHypoX) ) / tan(dipavg)
      endif

C     Now compute the additional distance terms for the radiation pattern. 
      uprime = Us - Uh
      tprime = Ts - Th
      BigR = sqrt ( uprime*uprime + tprime*tprime )
cnjg      zhyp = fltgrid_z(iHypoY,iHypoX)
      zhyp = fltgrid_z(iHypoZ,iHypoX)
      rh = sqrt ( zhyp*zhyp + BigR*BigR ) 

C     Now assign rake angles given fault mechanism type. 
C          Mech Type    Mech Flag     Rake Angle
C          -------------------------------------
C          Strike-slip     0.0            0.0
C          Reverse         1.0           60.0
C          Reverse/Obl     0.5           60.0
C          Normal         -1.0          -50.0
C          Normal/Obl     -0.5          -50.0
      if (flttype .eq. 0.0) then
          rake = 0.0*3.14159/180.0
      elseif (flttype .eq. 0.5) then
          rake = 60.0*3.14159/180.0
      elseif (flttype .eq. 1.0) then
          rake = 60.0*3.14159/180.0
      elseif (flttype .eq. -0.5) then
          rake = -50.0*3.14159/180.0
      elseif (flttype .eq. -1.0) then
          rake = -50.0*3.14159/180.0
      endif

      radterm1 = tprime*sin(dipavg)/rh + zhyp*cos(dipavg)/rh
      radterm2 = (zhyp*tprime*sin(dipavg))/(rh*BigR) - (BigR*cos(dipavg)/rh)
      radterm3 = uprime*sin(dipavg)/BigR
      radterm4 = -cos(rake)*uprime/rh + sin(rake)*cos(dipavg)*tprime/rh -
     1           sin(rake)*sin(dipavg)*zhyp/rh
      radterm5 = (-cos(rake)*zhyp*uprime)/(rh*BigR) + (sin(rake)*cos(dipavg)*zhyp*tprime)/(rh*BigR) +
     1           sin(rake)*sin(dipavg)*BigR/rh
      radterm6 = cos(rake)*tprime/BigR + sin(rake)*cos(dipavg)*uprime/BigR
      BigRu = (-uprime/BigR)*(radterm1*radterm5+radterm2*radterm4) +    
     1            (tprime/BigR)*(radterm1*radterm6+radterm3*radterm4)

      BigRt = (-tprime/BigR)*(radterm1*radterm5+radterm2*radterm4) -    
     1            (uprime/BigR)*(radterm1*radterm6+radterm3*radterm4)

      Radri = sqrt (BigRu*BigRu + BigRt*BigRt)
      Radri = max(Radri,0.2)

C     Compute the IDP term.
      IDP = BigC*BigS*Radri

C     Now call the specified directivity model based on dflag value:
c        dflag    Directivity Model
c          1      Average Horizontal Directivity (Abrahamson, 2000)
c          2      Fault Normal with Directivity (Abrahamson, 2000)
c          3      Fault Parallel with Directivity (Abrahamson, 2000)
c          4      Fault Normal without Directivity (Abrahamson, 2000)
c          5      Fault Parallel without Directivity (Abrahamson, 2000)
c          6      Abrahamson and Silva NGA (2008)
c          7      Boore and Atkinson NGA (2008)
c          8      Campbell and Bozorgnia NGA (2008)
c          9      Chiou and Youngs NGA (2008)
c         10      Abrahamson NGA (2010) - strike slip events only

      if (dflag .ge. 1 .and. dflag .lt. 6) then
c      if (dflag .eq. 1 ) then

          if (abs(flttype) .ne. 1) then

          call rupdirct_strike (period, rrup, mag, rupLen, rupwidth, flttype, tprime, 
     1                         uprime, s, lnDir, sigdirY, dflag, X, Y )
          endif

      endif

c     NGA Specific Directivity Models
      if (dflag .ge. 6 .and. dflag .le. 9) then
         call ngadircalc (dflag, period, IDP, rrup, mag, lnDir, sigdirY)
      endif

C     Abrahamson (2010) NGA Directivity Model. Note only applies for strike-slip sources. 
      if (dflag .eq. 10 .and. period .gt. 1.0) then
         if (flttype .eq. 0) then
            call dirmodelA2010 (dflag, period, rrup, mag, lnDir, sigdirY,
     1                          iHypoX, iHypoZ, icellRupStrike, icellRupdip, 
     2                          fltgrid_x, fltgrid_y, fltgrid_z, n2, n1 )
         endif
      endif


      return
      end

c -------------------------------------------------------------------
      subroutine dirmodelA2010 (dflag, per, rrup, mag, lnDir, sigDirY,
     1                          iHypoX, iHypoZ, icellRupStrike, icellRupdip, 
     2                          fltgrid_x, fltgrid_y, fltgrid_z, n2, n1 )

C     Subroutine to compute the directivity factors for the Abrahamson 2010 model. 
C     Model is only applicable to Strike-slip events. 

      include 'pfrisk.h'
c      include 'declare1.h'

      real c1T, b2T, a1T, a2T, mag, somag
      real per, rrup, sigDirY, lnDir, taper1, taper2, lnDirT3
      real so, sdist, p1p2ang, dx, dy, s
      real Period(7), c1(7), b2(7), a1(7), a2(7)
      integer nper, dflag, count1, count2
      integer iHypoX, iHypoZ, icellRupstrike, icellRupdip, n2, n1
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS)
      real  fltGrid_y(MAXFLT_DD,MAXFLT_AS)

      data period / 1.0, 2.00, 3.00, 4.00, 5.00, 7.00, 10.00 /
      data c1 / 10.0, 20.0, 30.0, 40.0, 50.0, 50.0, 50.0 /
      data b2 / 0.018, 0.041, 0.044, 0.037, 0.034, 0.028, 0.023 /
      data a1 / -2.07, -0.27, -0.32, 0.43, 0.42, 0.22, 0.04 /
      data a2 / -0.061, -0.201, -0.303, -0.371, -0.391, -0.380, -0.392 /

      nper = 7
C   Check for the need for interpolation of the sigma adjustments based 
c   on the requested spectral period.
      do i=1,nper-1
         if (per .ge. period(i) .and. per .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'Abrahamson 2010 directity model is not'
      write (*,*) 'defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',per
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   per,c1T,iflag)
            call interp (period(count1),period(count2),b2(count1),b2(count2),
     +                   per,b2T,iflag)
            call interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   per,a1T,iflag)
            call interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   per,a2T,iflag)

C     Compute the Taper 1 model based on rupture distance
      if (rrup .le. 10.0) then
         taper1 = 1.0
      elseif (rrup .lt. 40.0) then
         taper1 = (rrup - 10.0)/30.0
      else
         taper1 = 0.0
      endif
C     Compute taper 2 model based on magnitude
      if (mag .le. 6.0) then
         taper2 = 0.0
      elseif (mag .lt. 6.5) then
         taper2 = 2.0*(mag-6.0)
      else
         taper2 = 1.0
      endif

C     Now compute the length term So
      somag = 10.0 + 20.0*(mag-6.0)
      so = min(c1T,somag)

C     Now find the location of the cell for the So distance.
      sdist = 0.0
      if (icellRupstrike .gt. iHypoX) then
         do icell=icellRupstrike,iHypoX+1,-1
            dx = fltGrid_x(iHypoZ,icell) - fltGrid_x(iHypoZ,icell-1) 
            dy = fltGrid_y(iHypoZ,icell) - fltGrid_y(iHypoZ,icell-1) 
            sdist = sdist + sqrt(dx*dx + dy*dy)
            if (sdist .le. so) then
               iAS = icell
               s = sdist
            endif
         enddo
      elseif (icellRupstrike .lt. iHypoX) then
         do icell=iHypoX,icellRupstrike+1,-1
            dx = fltGrid_x(iHypoZ,icell) - fltGrid_x(iHypoZ,icell-1) 
            dy = fltGrid_y(iHypoZ,icell) - fltGrid_y(iHypoZ,icell-1) 
            sdist = sdist + sqrt(dx*dx + dy*dy)
            if (sdist .le. so) then
               iAS = icell
               s = sdist
            endif
         enddo
      elseif (icellRupstrike .eq. iHypoX) then
         iAS = icell
         s = 0.0
      endif

C     Now compute the angle between the station and the cell for So distance.
C     Assumes that the station is at location (0.0,0.0,0.0).
      dx = -fltGrid_x(iHypoZ,iAS)  
      dy = -fltGrid_y(iHypoZ,iAS) 
      p1p2ang = atan2(dx,dy)

C     Now compute directivity factor/adjustment.
      lnDir = (b2T*min(so,s)*cos(p1p2ang) - exp(a1T-a2T*(8.5-mag)**2.0))*taper1*taper2

C     Now adjust the sigma values. Adjustment is based on a reduction of 0.05 Ln Units for T=3.0sec.
      lnDirT3 = (b2(3)*min(so,s)*cos(p1p2ang) - exp(a1(3)-a2(3)*(8.5-mag)**2.0))*taper1*taper2

      sigDirY = 0.05*(lnDir/lnDirT3)

      write (94,'(5i4,14f10.3)') iHypoZ, iHypoX, icellRupDip, iCellRupStrike, iAS, mag, rrup,
     1       somag, c1T, so, s, p1p2ang*180.0/3.14159, taper1, taper2, lnDir, lnDirT3, sigDirY, 
     1       fltGrid_x(iHypoZ,iAS), fltgrid_y(iHypoZ,iAS)

      return
      end
c -------------------------------------------------------------------
      subroutine ngadirCalc (dflag, per, IDP, rrup, mag, lnDir, sigDirY)

C     This subroutine computes the NGA based Directivity factors given IDP and spectral period. 
    
      real IDP, lnDir, per, rrup, mag, fr, fm, sigDirY
      real Period(10), aAS(10), aBA(10), aCB(10), aCY(10)
      real bAS(10), bBA(10), bCB(10), bCY(10)
      real sigAS(10), sigBA(10), sigCB(10), sigCY(10)
      real aAST, aBAT, aCBT, aCYT, bAST, bBAT, bCBT, bCYT
      real sigAST, sigBAT, sigCBT, sigCYT
      integer dflag, iflag, count1, count2
      character*80 attenname1

      data period / 0.5, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00, 5.00, 7.50, 10.00 /
      data aAS / 0.0000, -0.0447, -0.0765, -0.1213, -0.1531, -0.1979,   
     1          -0.2296, -0.2542, -0.3636, -0.5755 /
      data aBA / 0.0000, -0.0532, -0.0910, -0.1443, -0.1821, -0.2353,   
     1          -0.2731, -0.3021, -0.4627, -0.8285 /
      data aCB / 0.0000, 0.0000, -0.0329, -0.0795, -0.1125, -0.1590,   
     1          -0.1921, -0.2172, -0.3227, -0.6419 /
      data aCY / 0.0000, 0.0000, -0.0260, -0.0627, -0.0887, -0.1254,   
     1          -0.1514, -0.1715, -0.2797, -0.4847 /

      data bAS / 0.0000, 0.0298, 0.0510, 0.0809, 0.1020, 0.1319,   
     1           0.1530, 0.1695, 0.2411, 0.3489 /
      data bBA / 0.0000, 0.0355, 0.0607, 0.0962, 0.1214, 0.1569,   
     1           0.1821, 0.2015, 0.2727, 0.4141 /
      data bCB / 0.0000, 0.0000, 0.0220, 0.0530, 0.0750, 0.1060,   
     1           0.1280, 0.1450, 0.2147, 0.3522 /
      data bCY / 0.0000, 0.0000, 0.0200, 0.0482, 0.0682, 0.0965,   
     1           0.1165, 0.1320, 0.1865, 0.2933 /

      data sigAS / 0.000, 0.002, 0.008, 0.010, 0.019, 0.040, 0.039,      
     1             0.041, 0.074, 0.160 / 
      data sigBA / 0.000, 0.001, 0.009, 0.008, 0.021, 0.052, 0.056,      
     1             0.077, 0.111, 0.205 / 
      data sigCB / 0.000, 0.000, 0.002, 0.002, 0.010, 0.032, 0.033,      
     1             0.050, 0.082, 0.174 / 
      data sigCY / 0.000, 0.000, 0.002, 0.003, 0.008, 0.022, 0.022,      
     1             0.023, 0.041, 0.098 / 

      nPer = 10

C     Check for given attenuation model is one of the four specific NGA models.
      if (dflag .eq. 6) then
          attenname1 = 'A&S, NGA Directivity Model'
      elseif (dflag .eq. 7 ) then          
          attenname1 = 'Boore&Atkinson, NGA Directivity Model'
      elseif (dflag .eq. 8 ) then          
          attenname1 = 'Campbell&Bozorgnia, NGA Directivity Model'
      elseif (dflag .eq. 9 ) then          
          attenname1 = 'Chiou&Youngs, NGA Directivity Model'
      endif

C   Check for the need for interpolation of the sigma adjustments based 
c   on the requested spectral period.
      do i=1,nper-1
         if (per .ge. period(i) .and. per .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020 
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*) 
      write (*,*) 'NGA model specific directity model is not'
      write (*,*) 'defined for a spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',per
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*) 
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020       call interp (period(count1),period(count2),aAS(count1),aAS(count2),
     +                   per,aAST,iflag)
            call interp (period(count1),period(count2),aBA(count1),aBA(count2),
     +                   per,aBAT,iflag)
            call interp (period(count1),period(count2),aCB(count1),aCB(count2),
     +                   per,aCBT,iflag)
            call interp (period(count1),period(count2),aCY(count1),aCY(count2),
     +                   per,aCYT,iflag)
            call interp (period(count1),period(count2),bAS(count1),bAS(count2),
     +                   per,bAST,iflag)
            call interp (period(count1),period(count2),bBA(count1),bBA(count2),
     +                   per,bBAT,iflag)
            call interp (period(count1),period(count2),bCB(count1),bCB(count2),
     +                   per,bCBT,iflag)
            call interp (period(count1),period(count2),bCY(count1),bCY(count2),
     +                   per,bCYT,iflag)
            call interp (period(count1),period(count2),sigAS(count1),sigAS(count2),
     +                   per,sigAST,iflag)
            call interp (period(count1),period(count2),sigBA(count1),sigBA(count2),
     +                   per,sigBAT,iflag)
            call interp (period(count1),period(count2),sigCB(count1),sigCB(count2),
     +                   per,sigCBT,iflag)
            call interp (period(count1),period(count2),sigCY(count1),sigCY(count2),
     +                   per,sigCYT,iflag)
       
C Compute the distance taper function.
      if (rrup .gt. 70.0) then
         fr = 0.0
      elseif (rrup .le. 40.0) then
         fr = 1.0
      else
         fr = 1.0 - (rrup-40.0)/30.0
      endif
       
C Compute the magnitude taper function.
      if (mag .gt. 6.0) then
         fm = 1.0
      elseif (mag .le. 5.6) then
         fm = 0.0
      else
         fm = (mag-5.6)/0.4
      endif

C Compute the directivity adjustment. 
C Adjust the sigma value.

       if (dflag .eq. 6) then
          lndir = fr*fm*(aAST+bAST*IDP)
          sigdirY = sigAST
       elseif (dflag .eq. 7) then
          lndir = fr*fm*(aBAT+bBAT*IDP)
          sigdirY = sigBAT
       elseif (dflag .eq. 8) then
          lndir = fr*fm*(aCBT+bCBT*IDP)
          sigdirY = sigCBT
       elseif (dflag .eq. 9) then
          lndir = fr*fm*(aCYT+bCYT*IDP)
          sigdirY = sigCYT
       endif

      return
      end

c -------------------------------------------------------------------

      subroutine Get_plane_dist (x, y, z, x0, y0, z0, insideFlag, dist)
      
      implicit none
      
c     declarations passed in      
      real x(5), y(5), z(5), x0, y0, z0
      
c     declarations passed out  
      integer insideFlag
      real dist

c     declarations only used in this subroutine    
      integer i, i1, nSeg
      real dist1

c     Determine if the site is inside the surface projection of the cell
c     boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, x, y, x0, y0, insideFlag )

c     Compute closest distance for inside
      if ( insideFlag .eq. 1 ) then
         call CalcPlaneDist ( x0, y0, z0, x, y, z, dist )
         return
      endif
      
c     Compute the shortest dist to each edge
      dist = 1.0e30
      do i=1,4
         i1 = i+1
         if ( i1 .gt. 4 ) i1=1
         call Calc_LineSeg_dist ( x(i), y(i), z(i), x(i1), 
     1       y(i1), z(i1), x0, y0, z0, dist1 )         
         if ( dist1 .lt. dist ) then
            dist = dist1
         endif
      enddo
      return
      end

c ---------------------------------------------------------------

      subroutine Calc_LineSeg_dist ( x1, y1, z1, x2, y2, z2, x0, y0,
     1           z0, dist )
      real x0, x1, x2, y0, y1, y2, z0, z1, z2, dist
      real t1, t2, x, y, z, L, L1, L2, d1, d2
 
c     Find shortest distance to line (without ends)
c     Interesection at (x,y,z)     
      if ( z1 .ne. z2 ) then
         t1 = (x2-x1)/(z2-z1) 
         t2 = (y2-y1)/(z2-z1)
         z =  (z0 - (-z1*t1 + x1 - x0)*t1  - (-z1*t2 + y1 - y0)*t2 )
     1     / ( t1**2 + t2**2 + 1 )
         x = t1 * (z-z1) + x1
         y = t2 * (z-z1) + y1
      elseif ( y1 .ne. y2 ) then
         z = z1
         t1 = (x2-x1)/(y2-y1)
         y = (y0 - (-y1*t1 + x1 - x0)*t1) / (t1**2 + 1)      
         x = t1 * (y-y1) + x1
      else
         z = z1
         y = y1
         x = x0
      endif
      dist = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )

c     Check if intersection is outside of edge
      L = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
      L1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 )
      L2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2 )
      if ( L1 .le. L .and. L2 .le. L ) then
         return
      endif
      
c     Intersection does not hit segment
      d1 = sqrt( (x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2 )
      d2 = sqrt( (x0-x2)**2 + (y0-y2)**2 + (z0-z2)**2 )
      dist = min ( d1, d2 )

      return
      end

c -----------------------------

      subroutine Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, 
     1           insideFlag )
     
      implicit none

c     declarations passed in      
      integer nSeg
      real xSeg(1), ySeg(1), x0, y0

c     declarations passed out
      integer insideFlag      

c     declarations only used within subroutine    
      integer i  
      real twoPi, pi, theta1, theta2, dTheta, dx1, dx2, dy1, dy2,
     1     sumTheta, test, tol

c     This subroutine determines if a point (x0,y0) is inside of
c     the polygon given by xSeg, ySeg.

      pi = 3.1415926
      twoPi = 2. * pi
      sumTheta = 0.
      
      do i=1,nSeg

c       Compute Azimuth to ends of segments
        dy1 = ySeg(i) - y0
        dy2 = ySeg(i+1) - y0
        dx1 = xSeg(i) - x0
        dx2 = xSeg(i+1) - x0
        theta1 = atan2 ( dy1, dx1 )
        theta2 = atan2 ( dy2, dx2 )
        dTheta = theta2 - theta1

c       Check if theta range is greater than pi (wrap around)
        if ( dTheta .gt. pi ) then
           dTheta = dTheta - twoPi
        elseif ( dTheta .lt. -pi ) then
           dTheta = dTheta + twoPi
        endif

c       Compute sum of azimuth ranges
        sumTheta = sumTheta + dTheta
      enddo

c     Determine if point is inside the polygon
c     If sumTheta = +-2 pi , then point is inside
      test = abs ( abs(sumTheta) - twoPi ) 
      tol = 0.01
      
      if ( test .lt. tol ) then
         insideFlag = 1
      else
         insideFlag = 0
      endif
      return
      end

c -------------------------------------------------------------

      subroutine CalcPlaneDist ( x0, y0, z0, x, y, z, dist )
      parameter (MAXPTS=3, MAXTERM=3)
c     implicit none 
      integer arow, acol, i, nterm, npts
      real*8 A(MAXPTS,MAXTERM), b(MAXPTS,1)
      real*8 xhat(MAXPTS,1), eps, deter, work(1000)
      real x(1), y(1), z(1), x0, y0, z0, dist, cx, cy, cz, c
      eps = 1.0e-11
      arow = MAXPTS
      acol = MAXTERM

c     Compute the equation for the plane
      if ( y(1) .eq. y(2) .and. z(1) .eq. z(2) ) then
         if ( z(1) .eq. z(2) .and. z(2) .eq. z(3) ) then
             cx = 0.
             cy = 0.
             cz = 1.
             c = -z(1)
         else
             do i=1,3
                A(i,1) = x(i)
                A(i,2) = z(i)
                A(i,3) = 1.
                b(i,1) = -y(i)
            enddo
            nTerm = 3
            nPts = 3
            call simul ( nterm, A, work, eps, -1, acol, deter )
            call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1         acol )
            cy = 1.
            cx = xhat(1,1)
            cz = xhat(2,1)
            c = xhat(3,1)
         endif
      else
         do i=1,3
            A(i,1) = y(i)
            A(i,2) = z(i)
            A(i,3) = 1.
            b(i,1) = -x(i)
         enddo
      nTerm = 3
      nPts = 3
      call simul ( nterm, A, work, eps, -1, acol, deter )
      call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1       acol )
      cx = 1.
      cy = xhat(1,1)
      cz = xhat(2,1)
      c = xhat(3,1)
      endif
      
c     Compute distance from point to plane
      dist = abs( cx*x0 + cy*y0 + cz*z0 + c) /
     1       sqrt( cx**2 + cy**2 + cz**2 )

      return
      end
      
      
c -----------------------------------------------------------------

      subroutine CalcBottomPt ( x1, y1, z1, strike, dip0, width,
     1           x3, y3, z3 )
      real width, dip0
      real L1, L2, Lz
      pi = 3.1415926

c     Points 1 and 2 are at the top of the plane      
      dip = dip0/180 * pi
      theta1 = strike + pi/2.
      
      if ( dip .lt. 0 ) then
         theta1 = theta1 + pi
      endif

c     Set coordinates of lower corners of plane
c     point 3 is below point 1
      L1 = cos(dip)*width
      L2 = L1
      Lz = abs(sin(dip))*width
      cos1 = cos(theta1)
      sin1 = sin(theta1) 
      x3 = x1 + cos1 * L1
      y3 = y1 + sin1 * L1
      z3 = z1 + Lz

      return
      end

c ------------------------------------------

      subroutine transpose ( A, arow, m, n, At, atrow )
      integer arow, atrow
      real*8 A(arow,1), At(atrow,1)

      do i=1,m
         do j=1,n
            At(j,i) = A(i,j)
         enddo
      enddo
      return
      end

c --------------------------------------------------------------------

      subroutine CalcDistDensity (nPts, xFlt2, yFlt2, distDensity,
     1                            dr, nr, x0, y0, step, minDist )
      include 'pfrisk.h'
      real xFlt2(MAX_DD,MAX_SEG), yFlt2(MAX_DD,MAX_SEG), minDist, maxdist
      real xFlt(MAX_SEG), yFlt(MAX_SEG)
      real distDensity(MAX_DIST1)
      
      nnGrid = 1

c     Copy to 1-D array
      iz = 1
      do i=1,nPts
        xFlt(i) = xFlt2(iz,i)
        yFlt(i) = yFlt2(iz,i)
      enddo

c     Set bounding rectangle for source
      xMin = 1.0e30
      xMax = -1.0e30
      yMin = 1.0e30
      yMax = -1.0e30
      iz = 1
       do i=1,nPts
        if ( xFlt(i) .lt. xMin ) xMin = xFlt(i)
        if ( xFlt(i) .gt. xMax ) xMax = xFlt(i)
        if ( yFlt(i) .lt. yMin ) yMin = yFlt(i)
        if ( yFlt(i) .gt. yMax ) yMax = yFlt(i)
      enddo

c     Initialize distDensity
      do i=1,MAX_DIST1
        distDensity(i) = 0.
      enddo
      minDist = 1.e10
      maxdist = 0.
      
c     Loop over bounding rectangle to estimate density function
      sum = 0.
      nx = (xMax-xMin)/step + 1
      ny = (yMax-yMin)/step + 1
      x = xMin + step/2.

      do ix=1,nx
        y = yMin + step/2.
        do iy=1,ny
          call Inside_OutSide ( nPts-1, xFlt, yFlt, x, y, insideFlag)
          if ( insideFlag .eq. 1 ) then
             horDist = sqrt( (x-x0)**2 + (y-y0)**2 )
             iBin = int( horDist / dr ) + 1             
             call CheckDim ( iBin, MAX_DIST1, 'MAX_DIST1  ' )

C     Reset Max and Min distances if needed
             if ( horDist .gt. maxDist ) then
               maxDist = horDist
               iBinMax = iBin               
             endif
             if (horDist .le. mindist) then
                minDist = horDist
             endif

             distDensity(iBin) = distDensity(iBin) + 1.
             sum = sum + 1.

          endif         

          y = y + step
        enddo
        x = x + step
      enddo

c     Normalize density function
      do iBin=1,iBinMax
        distDensity(iBin) = distDensity(iBin) / sum
      enddo
      nr = iBinMax
      return
      end
 
c --------------------------------------------------------------------

      subroutine CalcDistDensity1 ( iFlt, grid_a, grid_x, grid_y,
     1           grid_dx, grid_dy,
     1           grid_n, distDensity, dr, nr, x0, y0, step, minDist )
     
      include 'pfrisk.h'
c     implicit none

      real grid_a(MAX_FLT,MAX_GRID)
      real grid_x(MAX_GRID), grid_y(MAX_GRID)
      real grid_dx, grid_dy, minDist
      integer grid_n(MAX_FLT)
      real distDensity(MAX_DIST1)
      
      nnGrid = 1
    
      nx = grid_dx/step + 1
      ny = grid_dy/step + 1
      sum = 0.0
      minDist = 1.e10

C     Reset distdensity array.
      do i=1,max_dist1
         distdensity(i) = 0.0
      enddo
      
      do i=1,grid_n(iFlt)
        sum = sum + grid_a(iFlt,i)
        rate1 = grid_a(iFlt,i)/(nx*ny)
        x = step/2. + grid_x(i)
        do ix=1,nx
          y = step/2. + grid_y(i)
          do iy=1,ny
             horDist = sqrt( (x-x0)**2 + (y-y0)**2 )
             if (horDist .le. mindist) then
                minDist = horDist
             endif
             iBin = int( horDist / dr ) + 1
             call CheckDim ( iBin, MAX_DIST1, 'MAX_DIST1  ' )
             distDensity(iBin) = distDensity(iBin) + rate1
             y = y + step
          enddo
          x = x + step
        enddo
      enddo
      
c     Normalize density function
      do iBin=1,max_dist1
        distDensity(iBin) = distDensity(iBin) / sum
        if ( distDensity(iBin) .ne. 0 ) nr = iBin
      enddo
      return
      end
      
c --------------------------------------------------------------------

      subroutine CalcDistDensity2 ( iFlt, grid_a, grid_n, distDensity2 )
     
      include 'pfrisk.h'
c     implicit none

      real grid_a(MAX_FLT,MAX_GRID)
      integer grid_n(MAX_FLT)
      real distDensity2(MAX_GRID)
      
      sum = 0.
      do i=1,grid_n(iFlt)
        sum = sum + grid_a(iFlt,i)
      enddo
      
c     Normalize density function
      do i=1,grid_n(iFlt)
        distDensity2(i) = grid_a(iFlt,i) / sum
      enddo
      return
      end

c -------------------------------------------------------------------

      subroutine Get_plane_dist2 ( xRup, yRup, zRup, iSeg, 
     1           x0, y0, z0, dip0, dist, xclp, yclp, zclp )
c     implicit none
      real xRup(4,1), yRup(4,1), zRup(4,1)
      real x(4), y(4), z(4), xSeg(5), ySeg(5)
      real dist, dip0
      real x0, y0, z0, dist1
      real pi, dip, theta, theta1, sin1, cos1, 
     1     tan1, d1, d2
      real xclp, yclp, zclp, xclp2, yclp2, zclp2
      integer i, i1, nSeg, insideFlag, iSeg
      pi = 3.1415926

      do i=1,4
         x(i) = xRup(i,iSeg)
         y(i) = yRup(i,iSeg)
         z(i) = zRup(i,iSeg)
      enddo
      
c     Set angles
      dip = dip0/180 * pi
      theta = atan2 ( y(2)-y(1), x(2)-x(1) )
      theta1 = theta + pi/2.
      if ( dip .lt. 0 ) then
         theta1 = theta1 + pi
      endif
      cos1 = cos(theta1)
      sin1 = sin(theta1) 
      tan1 = tan(dip)
      
c     Set boundary for points for which the closest point is to 
C     the plane

C Note: Reorders points in clockwise rotation.
      d1 = z(1) * abs(tan1)
      d2 = z(3) * abs(tan1)
      xSeg(1) = x(1) + d1*cos1
      xSeg(2) = x(2) + d1*cos1
      xSeg(3) = x(4) + d2*cos1
      xSeg(4) = x(3) + d2*cos1
      xSeg(5) = xSeg(1)
      ySeg(1) = y(1) + d1*sin1
      ySeg(2) = y(2) + d1*sin1
      ySeg(3) = y(4) + d2*sin1
      ySeg(4) = y(3) + d2*sin1
      ySeg(5) = ySeg(1)

c     Determine if the site is inside this boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, insideFlag )
      
c     Compute closest distance for inside
      if ( insideFlag .eq. 1 ) then
         call CalcPlaneDist2 ( x0, y0, z0, x, y, z, dist, xclp, yclp, zclp )
         return
      endif
      
c     Compute the shortest dist to each edge
      dist = 1.0e30
      xclp = 1.0e30
      yclp = 1.0e30
      zclp = 1.0e30
      do i=1,4
         i1 = i+1
         if ( i1 .gt. 4 ) i1=1
         call Calc_LineSeg_dist2 ( x(i), y(i), z(i), x(i1), 
     1       y(i1), z(i1), x0, y0, z0, dist1, xclp2, yclp2, zclp2 )

         if ( dist1 .lt. dist ) then
            dist = dist1
            xclp = xclp2
            yclp = yclp2
            zclp = zclp2
         endif
      enddo
      return
      end

c ---------------------------------------------------------------

      subroutine Calc_LineSeg_dist2 ( x1, y1, z1, x2, y2, z2, x0, y0,
     1           z0, dist, x, y, z )
      real x0, x1, x2, y0, y1, y2, z0, z1, z2, dist
      real t1, t2, x, y, z, L, L1, L2, d1, d2
 
c     Find shortest distance to line (without ends)
c     Interesection at (x,y,z)     
      if ( z1 .ne. z2 ) then
         t1 = (x2-x1)/(z2-z1) 
         t2 = (y2-y1)/(z2-z1)
         z =  (z0 - (-z1*t1 + x1 - x0)*t1  - (-z1*t2 + y1 - y0)*t2 )
     1     / ( t1**2 + t2**2 + 1 )
         x = t1 * (z-z1) + x1
         y = t2 * (z-z1) + y1
      elseif ( y1 .ne. y2 ) then
         z = z1
         t1 = (x2-x1)/(y2-y1)
         y = (y0 - (-y1*t1 + x1 - x0)*t1) / (t1**2 + 1)      
         x = t1 * (y-y1) + x1
      else
         z = z1
         y = y1
         x = x0
      endif
      dist = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )

c     Check if intersection is outside of edge
      L = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
      L1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 )
      L2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2 )
      if ( L1 .le. L .and. L2 .le. L ) then
         return
      endif
      
c     Intersection does not hit segment
      d1 = sqrt( (x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2 )
      d2 = sqrt( (x0-x2)**2 + (y0-y2)**2 + (z0-z2)**2 )
      dist = min ( d1, d2 )

      return
      end

c -------------------------------------------------------------

      subroutine CalcPlaneDist2 ( x0, y0, z0, x, y, z, dist, cx, cy, cz )
      parameter (MAXPTS=3, MAXTERM=3)
c     implicit none 
      integer arow, acol, i, nterm, npts
      real*8 A(MAXPTS,MAXTERM), b(MAXPTS,1)
      real*8 xhat(MAXPTS,1), eps, deter, work(1000)
      real x(1), y(1), z(1), x0, y0, z0, dist, cx, cy, cz, c
      eps = 1.0e-11
      arow = MAXPTS
      acol = MAXTERM

c     Compute the equation for the plane
      if ( y(1) .eq. y(2) .and. z(1) .eq. z(2) ) then
         if ( z(1) .eq. z(2) .and. z(2) .eq. z(3) ) then
             cx = 0.
             cy = 0.
             cz = 1.
             c = -z(1)
         else
             do i=1,3
                A(i,1) = x(i)
                A(i,2) = z(i)
                A(i,3) = 1.
                b(i,1) = -y(i)
            enddo
            nTerm = 3
            nPts = 3
            call simul ( nterm, A, work, eps, -1, acol, deter )
            call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1         acol )
            cy = 1.
            cx = xhat(1,1)
            cz = xhat(2,1)
            c = xhat(3,1)
         endif
      else
         do i=1,3
            A(i,1) = y(i)
            A(i,2) = z(i)
            A(i,3) = 1.
            b(i,1) = -x(i)
         enddo
      nTerm = 3
      nPts = 3
      call simul ( nterm, A, work, eps, -1, acol, deter )
      call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1       acol )
      cx = 1.
      cy = xhat(1,1)
      cz = xhat(2,1)
      c = xhat(3,1)
      endif
      
c     Compute distance from point to plane
      dist = abs( cx*x0 + cy*y0 + cz*z0 + c) /
     1       sqrt( cx**2 + cy**2 + cz**2 )

      return
      end
      
      

