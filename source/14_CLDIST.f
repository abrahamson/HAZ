
      subroutine S14_CalcDist (sourceType, pscorflag, nFltGrid, n1AS, iLocX, iLocY, n2AS, 
     1             iFltWidth, iFlt, iMag, ystep, grid_top, RupWidth, RupLen, r_horiz, seismoDepth, 
     2             fltgrid_x, fltgrid_y, fltgrid_z, fltgrid_x1, fltgrid_y1, fltgrid_z1, 
     3             fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3, fltgrid_z3, 
     4             fltgrid_x4, fltgrid_y4, fltgrid_z4, fltgrid_Rrup, fltgrid_Rjb, dip, dipS7,
     5             distS7, HWFlag, n1, n2, icellRupstrike, icellRupdip, hypoDepth, distJB, 
     6             distRup, ZTOR, distSeismo, distepi, disthypo, dipavgd, Rx, Ry, Ry0,
     7             h_listric, dMag1_listric, logA_shal, logA_deep,
     9             iDD_shal_min, iDD_shal_max, iDD_deep_min, iDD_deep_max)
           
c     This subroutine computes the distance parameters for the GMPEs
c     (e.g. Rjb, Rrup, ZTOR, Rx, Ry, Ry0)

      implicit none
      include 'pfrisk.h'

c     declarations passed in 
      integer sourceType, pscorflag, nfltGrid(2), n1AS(MAXFLT_AS), iLocX, 
     1        iLocY, n2AS(MAXFLT_AS), iFltWidth, iFlt, iMag
      real ystep, grid_top(MAX_FLT,MAX_GRID), RupWidth, RupLen, r_horiz, seismoDepth,
     1     fltgrid_x(MAXFLT_DD,MAXFLT_AS), fltgrid_y(MAXFLT_DD,MAXFLT_AS),
     2     fltgrid_z(MAXFLT_DD,MAXFLT_AS), fltgrid_x1(MAXFLT_DD,MAXFLT_AS), 
     3     fltgrid_y1(MAXFLT_DD,MAXFLT_AS), fltgrid_z1(MAXFLT_DD,MAXFLT_AS), 
     4     fltgrid_x2(MAXFLT_DD,MAXFLT_AS), fltgrid_y2(MAXFLT_DD,MAXFLT_AS), 
     5     fltgrid_z2(MAXFLT_DD,MAXFLT_AS), fltgrid_x3(MAXFLT_DD,MAXFLT_AS), 
     6     fltgrid_y3(MAXFLT_DD,MAXFLT_AS), fltgrid_z3(MAXFLT_DD,MAXFLT_AS),
     6     fltgrid_x4(MAXFLT_DD,MAXFLT_AS), fltgrid_y4(MAXFLT_DD,MAXFLT_AS), 
     7     fltgrid_z4(MAXFLT_DD,MAXFLT_AS), fltgrid_Rrup(MAXFLT_DD,MAXFLT_AS),
     8     fltgrid_Rjb(MAXFLT_DD,MAXFLT_AS), dip(MAX_FLT,MAX_WIDTH,MAX_SEG),
     9     distS7(MAX_FLT,MAX_S7), dipS7(MAX_FLT,MAX_S7)
      real h_listric, dMag1_listric
      real logA_shal, logA_deep
      integer iDD_shal_min, iDD_shal_max, iDD_deep_min, iDD_deep_max

c     declarations passed out
      integer HWFlag, n1, n2, icellRupstrike, icellRupdip  
      real hypoDepth, distJB, distRup, ZTOR, distSeismo, distepi, disthypo,
     1     dipavgd, Rx, Ry, Ry0

c     declarations only used within subroutine
      integer i, j, celly, cellx
      real top_grid, hypo1, rupLength2, r_horiz1, testy, testx, hypox, hypoy
      real magUpper, magLower, h_dist, v_dist
      
c     Initialize closest distance      
      distRup = 1.0e30
      distJB = 1.0e30
      distSeismo = 1.0e30
      disthypo = 1.0e30
      distepi = 1.0e30
      
c     Set top_grid for sourceTypes 2/3 and 4 (different)
      if (sourceType .eq. 2 .or. sourceType .eq. 3) then
        top_grid = grid_top(iFlt,1)
      else if (sourceType .eq. 4) then
        top_grid = grid_top(iFlt,iLocX)
      endif
      
c     Compute distances for areal sources      
      if ( sourceType .eq. 2 .or. sourceType .eq. 3 .or. sourceType .eq. 4 ) then

c      Set hypoDepth for areal sources
       hypoDepth = (iLocY-0.5)*ystep + top_grid
 
c      Approximate correction for using point source with extended source model
C      added limit on only shallow areal sources (Hypo>30.0km)
       if ( psCorFlag .eq. 1  .and. hypoDepth .le. 30.0) then
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
       
c      Set other distance values equal to corresponding distance measures.
       distepi = distJB
       disthypo = distRup
       Rx = distJB
       Ry = 0.0
       Ry0 = 0.0

C      Set HWflag = 0 for areal or grid sources.
       HWFlag = 0    
       dipavgd = dip(iFlt,ifltWidth,1)
       return
          
c     Compute distances for SourceType 7
      elseif (sourceType .eq. 7) then
        distrup = distS7(iFlt,iMag)
        distJB  = distS7(iFlt,iMag)
        distseismo = distS7(iFlt,iMag)
        dipavgd = dipS7(iFlt,iMag)
        disthypo = distS7(iFlt,iMag)

c       Fixed Parameters
        ZTOR = 0.0
        HWFlag = 0
        hypoDepth = 8.0
        Rx = distrup

c     Otherwise use fault sources (i.e., Sourcetype = 1, 5, or 6)
      else 
     
c       n1 is the last cell that makes up the rupture plane down dip
c       n2 is the last cell that makes up the rupture plane along strike       
        n1 = nfltgrid(1)-n1AS(iLocX)+iLocY
        n2 = n2AS(iLocX)
          
c       Calculate hypoDepth (assumes hypocenter is in center of rupture plane)
        if (sourceType .eq. 1 .or. sourceType .eq. 5 .or. sourceType .eq. 6) then       
           testy = ((n1-iLocY)+1.)/2.
           celly = int(iLocY+testy) 
           testx = ((n2-iLocX)+1.)/2.
           cellx = int(iLocX+testx)           
           if (int(testy) .eq. testy .and. int(testx) .eq. testx) then                    
             hypoDepth = fltgrid_z1(celly,cellx)
             hypox = fltgrid_x1(celly,cellx)
             hypoy = fltgrid_y1(celly,cellx) 
           else if (int(testy) .eq. testy .and. int(testx) .ne. testx) then
             hypoDepth = (fltgrid_z2(celly,cellx) + fltgrid_z1(celly,cellx))/2.
             hypox = (fltgrid_x2(celly,cellx) + fltgrid_x1(celly,cellx))/2.
             hypoy = (fltgrid_y2(celly,cellx) + fltgrid_y1(celly,cellx))/2.
           else if (int(testy) .ne. testy .and. int(testx) .eq. testx) then
             hypoDepth = (fltgrid_z4(celly,cellx) + fltgrid_z1(celly,cellx))/2.  
             hypox = (fltgrid_x4(celly,cellx) + fltgrid_x1(celly,cellx))/2.   
             hypoy = (fltgrid_y4(celly,cellx) + fltgrid_y1(celly,cellx))/2.             
           else 
             hypoDepth = fltgrid_z(celly,cellx)    
             hypox = fltgrid_x(celly,cellx)
             hypoy = fltgrid_y(celly,cellx)                 
           endif
        endif

c       Calculate Rx, Ry, Ry0, dipavgd, and HWFlag        
c       Global Coordinate System 2 method
        call S22_GC2 (iLocX, iLocY, n2, n1, fltgrid_x1, fltgrid_y1, fltgrid_z1, 
     1              fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3,
     2              fltgrid_z3, fltgrid_x4, fltgrid_y4, fltgrid_z4, Rx, Ry, Ry0, 
     3              HWFlag, dipavgd) 

c       NAA, added for listric for source type 5  
c       apply the Taiwan model for setting the average dip for listric faults     
        if (sourceType .eq. 5 .and. h_listric .gt. 0. ) then
          magUpper = logA_shal + 4.0
          magLower = logA_deep + 4.0

c         logArea1 is the log10 of the area of the rupture shallower than h_listic            
c         logArea1 is the log10 of the area of the rupture deepter than h_listic            
          if ( magUpper - MagLower .gt. dMag1_listric ) then
c           compute dip of the shallow part
            h_dist = sqrt( (fltgrid_x4(iDD_shal_max,iLocY) - fltgrid_x1(iDD_shal_min,iLocY))**2 +
     1               (fltgrid_y4(iDD_shal_max,iLocY) - fltgrid_y1(iDD_shal_min,iLocY))**2 )
            v_dist = abs( fltgrid_z4(iDD_shal_max,iLocY) - fltgrid_z1(iDD_shal_min,iLocY) )
            dipavgd = atan2(v_dist,h_dist)*180/3.1415926
c            write (*,'( 2x,''shallow'',f10.3)') dipavgd
            
          elseif ( magUpper - MagLower .lt. 0. ) then 

c           compute dip of the deep part
            h_dist = sqrt( (fltgrid_x4(iDD_deep_max,iLocY) - fltgrid_x1(iDD_deep_min,iLocY))**2 +
     1               (fltgrid_y4(iDD_deep_max,iLocY) - fltgrid_y1(iDD_deep_min,iLocY))**2 )
            v_dist = abs( fltgrid_z4(iDD_deep_max,iLocY) - fltgrid_z1(iDD_deep_min,iLocY) )
            dipavgd = atan2(v_dist,h_dist)*180/3.1415926
c            write (*,'( 2x,''deep'',f10.3)') dipavgd
          endif
        endif

c       Calculate ZTOR
        ZTOR = fltgrid_z1(iLocY,iLocX)             

c       Calculate Rrup and Rjb
C       Keep track of the fault grid cell for the closest rupture distance
C       for this rupture area since it is needed for the NGA directivity models.     
     
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

C           Compute the DistSeismo value based on min depth being equal to seismoDepth parameter.
            if (fltgrid_z(i,j) .gt. seismoDepth) then                 
              if (distSeismo .gt. fltgrid_Rrup(i,j)) then
                distSeismo = fltgrid_Rrup(i,j)
              endif
            endif 
          enddo
        enddo
c        write (*,'(4i5,f10.3)') iLocX, n2, iLocY, n1, distRup

C       Set the Epi and Hypo distances for faults.
        distepi = sqrt(hypox**2.+hypoy**2.)
        disthypo = sqrt(hypox**2.+hypoy**2.+hypoDepth**2.)
      
      endif

      return
      end

c ----------------------------------------------------------------------

      subroutine S14_Set_MinDist (sourceType, iFlt, iFltWidth, distRup, distJB, distSeismo, 
     1                        SourceDist, MinRrup_temp, MinRjb_temp, MinSeismo_temp)

      implicit none
      include 'pfrisk.h'
      
      integer sourceType, iFlt, iFltWidth
      real distRup, distJB, distSeismo, SourceDist(MAX_FLT,MAX_WIDTH,3),
     1     MinRrup_temp, MinRjb_temp, MinSeismo_temp
      
       if ( distRup .lt. SourceDist(iFlt,iFltWidth,1) ) then
         SourceDist(iFlt,iFltWidth,1)=distRup
       endif
       if ( distJB .lt. SourceDist(iFlt,iFltWidth,2) ) then
         SourceDist(iFlt,iFltWidth,2)=distJB
       endif
       if ( distSeismo .lt. SourceDist(iFlt,iFltWidth,3) ) then
         SourceDist(iFlt,iFltWidth,3)=distSeismo
       endif
       MinRrup_temp = SourceDist(iFlt,iFltWidth,1)
       MinRjb_temp = SourceDist(iFlt,iFltWidth,2)
       MinSeismo_temp = SourceDist(iFlt,iFltWidth,3)

      return
      end

c ----------------------------------------------------------------------

      subroutine S14_SetnRupLoc ( n1, n2, nHypoX, pHypoX, nHypoXStep, 
     1                        nHypoZ, pHypoZ, nHypoZstep ) 

      implicit none

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
      subroutine S14_DetermDist (hAScell, hDDcell, icellRupStrike, icellRupdip, 
     1                      fltgrid_x, fltgrid_y, fltgrid_z, n2, n1, dipavgd,
     2                      iLocY, iLocX, rupLen, rupWidth, x0, y0, z0,
     3                      edist, hdist, slit, azp1p2, step, dlit, phiang, FltGrid_rRup,
     4                      s2site, Rx, astrike)

      implicit none
      include 'pfrisk.h'

      integer hAScell, hDDcell, icellRupStrike, icellRupDip
      integer n2, n1, iLocx, iLocY, icell
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS),
     1      fltGrid_y(MAXFLT_DD,MAXFLT_AS), dx, dy, dz, edist, hdist
      real dist, step, phiang, dlit, strikeX, strikeY, astrike
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), Rx, dipavgd, x0, y0, z0, rupLen,
     1     rupWidth, azp1p2, slit, s2site, c1temp

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
      phiang = (180.0/3.14159)*atan2(Rx,dlit) - (90.0 - dipavgd)

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

      subroutine S14_Get_plane_dist (x, y, z, x0, y0, z0, insideFlag, dist)
      
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
      call S14_Inside_OutSide ( nSeg, x, y, x0, y0, insideFlag )
      
c     Compute the shortest dist to each edge
      dist = 1.0e30
      do i=1,4
         i1 = i+1
         if ( i1 .gt. 4 ) i1=1
         call S14_Calc_LineSeg_dist ( x(i), y(i), z(i), x(i1), 
     1       y(i1), z(i1), x0, y0, z0, dist1 )         
         if ( dist1 .lt. dist ) then
            dist = dist1
         endif
      enddo
      return
      end

c ---------------------------------------------------------------

      subroutine S14_Calc_LineSeg_dist ( x1, y1, z1, x2, y2, z2, x0, y0,
     1           z0, dist )
     
      implicit none
     
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

      subroutine S14_Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, 
     1           insideFlag )
     
      implicit none
      include 'pfrisk.h'

c     declarations passed in      
      integer nSeg
      real xSeg(MAX_SEG), ySeg(MAX_SEG), x0, y0

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

      subroutine S14_CalcPlaneDist ( x0, y0, z0, x, y, z, dist )

      implicit none

      integer MAXPTS, MAXTERM, arow, acol, i, nterm, npts
      parameter (MAXPTS=3, MAXTERM=3)       
      real*8 A(MAXPTS,MAXTERM), b(MAXPTS,1)
      real*8 xhat(MAXPTS,1), eps, deter, work(1000)
      real x(5), y(5), z(5), x0, y0, z0, dist, cx, cy, cz, c
      
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
            call S26_simul ( nterm, A, work, eps, -1, acol, deter )
            call S26_mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
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
      call S26_simul ( nterm, A, work, eps, -1, acol, deter )
      call S26_mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
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

c --------------------------------------------------------------------

      subroutine S14_CalcDistDensity (nPts, xFlt2, yFlt2, distDensity,
     1                            dr, nr, x0, y0, step)

      implicit none  
      include 'pfrisk.h'

      integer i, ix, iy, nx, ny, nPts, nr, insideFlag, iBin, iz, iBinMax
      real xFlt2(MAX_DD,MAX_SEG), yFlt2(MAX_DD,MAX_SEG), maxdist
      real xFlt(MAX_SEG), yFlt(MAX_SEG), xMin, xMax, yMin, yMax, sum
      real distDensity(MAX_DIST1), dr, step, x, y, x0, y0, horDist

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
      maxdist = 0.
      
c     Loop over bounding rectangle to estimate density function
      sum = 0.
      nx = (xMax-xMin)/step + 1
      ny = (yMax-yMin)/step + 1
      x = xMin + step/2.
      
      do ix=1,nx
        y = yMin + step/2.
        do iy=1,ny
          call S14_Inside_OutSide ( nPts-1, xFlt, yFlt, x, y, insideFlag)
          if ( insideFlag .eq. 1 ) then
             horDist = sqrt( (x-x0)**2 + (y-y0)**2 )
             iBin = int( horDist / dr ) + 1             
             call S21_CheckDim ( iBin, MAX_DIST1, 'MAX_DIST1  ' )

C     Reset Max and Min distances if needed
             if ( horDist .gt. maxDist ) then
               maxDist = horDist
               iBinMax = iBin               
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

      subroutine S14_CalcDistDensity1 ( iFlt, grid_a, grid_x, grid_y,
     1           grid_dx, grid_dy, grid_n, distDensity, dr, nr, 
     2           x0, y0, step)
     
      implicit none
      include 'pfrisk.h'
 
      integer grid_n(MAX_FLT), iFlt, nr, i, nx, ny, ix, iy, iBin
      real grid_a(MAX_FLT,MAX_GRID), grid_x(MAX_GRID), grid_y(MAX_GRID),
     1     grid_dx, grid_dy, distDensity(MAX_DIST1), dr, x0, y0,
     2     step, x, y, horDist, sum, rate1
    
      nx = grid_dx/step + 1
      ny = grid_dy/step + 1
      write (*,'( 3f10.4)')grid_dx, grid_dy, step
      sum = 0.0

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
             iBin = int( horDist / dr ) + 1
             call S21_CheckDim ( iBin, MAX_DIST1, 'MAX_DIST1  ' )
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

      subroutine S14_CalcDistDensity2 ( iFlt, grid_a, grid_n, distDensity2 )
     
      
      implicit none
      include 'pfrisk.h'

      integer grid_n(MAX_FLT), iFlt, i
      real grid_a(MAX_FLT,MAX_GRID), distDensity2(MAX_GRID), sum
      
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

      subroutine S14_Get_plane_dist2 ( xRup, yRup, zRup, iSeg, 
     1           x0, y0, z0, dip0, dist, xclp, yclp, zclp )

      implicit none
      
      integer i, i1, nSeg, insideFlag, iSeg
      real xRup(4,1), yRup(4,1), zRup(4,1), x(4), y(4), z(4), xSeg(5), 
     1     ySeg(5), dist, dip0, x0, y0, z0, dist1, pi, dip, theta, 
     2     theta1, sin1, cos1, tan1, d1, d2, xclp, yclp, zclp, xclp2, 
     3     yclp2, zclp2

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
      call S14_Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, insideFlag )
      
c     Compute closest distance for inside
      if ( insideFlag .eq. 1 ) then
         call S14_CalcPlaneDist2 ( x0, y0, z0, x, y, z, dist, xclp, yclp, zclp )
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
         call S14_Calc_LineSeg_dist2 ( x(i), y(i), z(i), x(i1), 
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

      subroutine S14_Calc_LineSeg_dist2 ( x1, y1, z1, x2, y2, z2, x0, y0,
     1           z0, dist, x, y, z )
     
      implicit none
     
      real x0, x1, x2, y0, y1, y2, z0, z1, z2, dist, t1, t2, x, y, z, 
     1     L, L1, L2, d1, d2
 
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

      subroutine S14_CalcPlaneDist2 ( x0, y0, z0, x, y, z, dist, cx, cy, cz )

      implicit none

      integer MAXPTS, MAXTERM, arow, acol, i, nterm, npts
      parameter (MAXPTS=3, MAXTERM=3)
      real x(4), y(4), z(4), x0, y0, z0, dist, cx, cy, cz, c      
      real*8 A(MAXPTS,MAXTERM), b(MAXPTS,1), xhat(MAXPTS,1), eps, 
     1       deter, work(1000)

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
            call S26_simul ( nterm, A, work, eps, -1, acol, deter )
            call S26_mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
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
      call S26_simul ( nterm, A, work, eps, -1, acol, deter )
      call S26_mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
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
      
      

