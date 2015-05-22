c closestDist.f

      subroutine CalcDist ( iLocX, iLocY, RupLen, RupWidth, fltGrid_rRup,
     1             fltGrid_RJB, fltGrid_Rx, fltGrid_HW, fltGrid_Rseis, 
     2             nFltGrid, distRup, distJB, distSeismo, Rx, hwFlag, zTOR,
     3             r_horiz, hypoDepth, fltGRid_z, sourceType, iiflag,   
     4             fltgrid_x, fltgrid_y, iCellRupStrike, iCellRupDip, n2, n1, dipavg, pscorflag,
     5             seismoDepth, distepi, disthypo, stepX, stepY, runflag, iSite,
     6             hDDkm, hASkm, fltgrid_a, n1AS, n2AS, Ry, Ry0  )        
     
c     This subroutine computes the rupture distance and JB distance
c     (e.g. closest distance from the rupture plane to the site)

      include 'pfrisk.h'
      integer nfltGrid(2), sourceType, HWFLAG, FWFlag, runflag, iSite
      real fltGrid_w(MAXFLT_DD,MAXFLT_AS),  fltGrid_a(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_RJB(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_Rseis(MAXFLT_DD,MAXFLT_AS), fltGrid_Rx(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_HW(MAXFLT_DD,MAXFLT_AS)
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS),
     1      fltGrid_y(MAXFLT_DD,MAXFLT_AS)   
      real xtest(5), ytest(5), xtestfw(5), ytestfw(5), Rx, dipavg, Ruplen, seismoDepth
      real distepi, disthypo, hDDkm, hASkm, Ry, MidDist
      real stepX, stepY, deltaX, deltaY, deltaZ
      integer pscorflag, n1var(MAXFLT_DD,MAXFLT_AS), n1AS(MAXFLT_AS), n2AS(MAXFLT_AS)
      integer mAS, mDD
      real rup_E(MAXFLT_AS), rup_N(MAXFLT_AS), Seg_length(MAXFLT_AS), Strike_slope(MAXFLT_AS),
     1     Normal_slope(MAXFLT_AS), a, P90_E(MAXFLT_AS), P90_N(MAXFLT_AS), Site_E, Site_N,
     2     t_local(MAXFLT_AS), u_local(MAXFLT_AS), Seg_weight(MAXFLT_AS), 
     3     Seg_weight_t(MAXFLT_AS), sum_Weight, rec_Weight,
     4     Seg_x(MAXFLT_AS), Seg_wxu(MAXFLT_AS), sum_Swt, sum_Swxu,
     5     Global_T, Global_U
      integer iGC2 
         

c     Initialize closest distance      
      distRup = 1.0e30
      distJB = 1.0e30
      distSeismo = 1.0e30
      disthypo = 1.0e30
      distepi = 1.0e30
      distRx = 1.0e30
      
c     Check for areal sources      
      if ( sourceType .eq. 2 .or. sourceType .eq. 3 .or. sourceType .eq. 4 ) then
 
c      Approximate corrction for using point source with extended source model
C     added limit on only shallow areal sources (Hypo>30.0km)
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
c     Set other distance values equal to corresponding distance measures.
       distepi = distJB
       distRx = distJB
       disthypo = distRup
C     Set HWflag = 0 for areal or grid sources.
       hwflag = 0
C     Set Dipaverage equal to 90 for gridded sources.
       dipavg = 3.14159/2.0 
       return

c     Otherwise use fault sources (i.e., Sourcetype = 1, 5, or 6)
      else
cnjg       n2 = int(RupLen) + iLocX
       n2 = int(RupLen/StepX) + iLocX
       if ( n2 .eq. 0 ) n2 = 1
       if ( n2 .gt. nFltGrid(2)) n2=nFltGrid(2)
cnjg       n1 = int(RupWidth) + iLocY
       n1 = int(RupWidth/StepY) + iLocY
       if ( n1 .eq. 0 ) n1 = 1
       if ( n1 .gt. nFltGrid(1)) n1=nFltGrid(1)
                     
C     New method for the summation of fault areas to get the right number of cells 
C     given the expected rupture area for a magnitude.
      sumA = 0.0
      RupArea = RupLen*RupWidth
C     Start with loop over given along strike location then loop down dip     
        do in2=iLocX, nFltgrid(2)
c     Option of limiting rupture area no greater than Rupwidth. 
c         do in1=iLocY, n1var(iLocY,in2)

         
c CH: ------------- start        
c     save the rupture cell locations along strike in two arrays 
c     with their easting and northing coordinates
            
         rup_E(in2) = fltgrid_x(iLocY,in2)
         rup_N(in2) = fltgrid_y(iLocY,in2)
c CH: ------------- end
         
         do in1=iLocY, nfltgrid(1)-n1AS(iLocX)+iLocY
            sumA = sumA + fltgrid_a(in1,in2)         
         enddo
C     Check to see that summed rupture area is equal to requested rupture area. 
C     If not sum up another column of cells only downdip the preset number of cells for this rupture.
         if (sumA .ge. RupArea) then
c            n1 = in1
            n2 = in2           
            goto 2345
         endif       
      enddo                  
2345  continue
        

c CH: -------------- start
c     calculate the length of each segment of the rupture (a segment of the
c     rupture is from the center of one rupture cell to the center of the 
c     next rupture cell, so the segment length should be approximately equal
c     to the step size along strike)         

        Site_E = 0.0
        Site_N = 0.0

        do iGC2=1, n2-1
          Seg_length(iGC2) = sqrt(((rup_E(iGC2+1)-rup_E(iGC2))**2)+((
     1                    rup_N(iGC2+1)-rup_N(iGC2))**2))
          Strike_slope(iGC2) = (rup_N(iGC2+1)-rup_N(iGC2))/(rup_E(iGC2+1)-
     1                    rup_E(iGC2))
          Normal_slope(iGC2) = (-1)/Strike_slope(iGC2) 
        enddo      

c       calculate Easting and Northing coordinate of point where local 
c       u and t form 90 degree angle 

          do iGC2=1, n2-1
            a = rup_E(iGC2+1) - rup_E(iGC2)  
              if (a.EQ.0) then
                P90_E(iGC2) = rup_E(iGC2)
                P90_N(iGC2) = Site_N
              else              
                P90_E(iGC2) = (((-1)*Strike_slope(iGC2)*rup_E(iGC2)) +
     1                     rup_N(iGC2) + (Normal_slope(i)*Site_E) -
     2                     Site_N)/(Normal_slope(iGC2) - 
     3                     Strike_slope(iGC2))     
                P90_N(iGC2) = (Strike_slope(iGC2) * (P90_E(iGC2) - 
     1                     rup_E(iGC2))) + rup_N(iGC2)  
              endif
          enddo

c       calculate local strike-normal coordinate t  

          do iGC2=1, n2-1
              if (Site_E.GT.P90_E(iGC2)) then
                t_local(iGC2) = -sqrt(((P90_E(iGC2) - Site_E)**2) +
     1                         ((P90_N(iGC2) - Site_N)**2)) 
              else              
                t_local(iGC2) = sqrt(((P90_E(iGC2) - Site_E)**2) +
     1                         ((P90_N(iGC2) - Site_N)**2)) 
              endif
          enddo

c       calculate local strike-parallel coordinate u  

          do iGC2=1, n2-1
              if (P90_N(iGC2).GT.rup_N(iGC2)) then
                u_local(iGC2) = -sqrt(((P90_E(iGC2) - rup_E(iGC2))**2) +
     1                         ((P90_N(iGC2) - rup_N(iGC2))**2)) 
              else              
                u_local(iGC2) = sqrt(((P90_E(iGC2) - rup_E(iGC2))**2) +
     1                         ((P90_N(iGC2) - rup_N(iGC2))**2)) 
              endif
          enddo

c        check for special case t=0, on extension and assign 
c          alternative segment weight
c        check for special case t=0, on segment and assign dummy
c          segment weight = 0, which will not be used 
c        if no special case, assign Segment weight as analytical 
c          solution of 1/(r^2) evaluated from 0 to Seg_length

          do iGC2=1, n2-1
            if (t_local(iGC2).EQ.0) then
              if (u_local(iGC2).LT.0) then
                Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) - 
     1                            (1/u_local(iGC2))
              else if (u_local(iGC2).GT.Seg_length(iGC2)) then
                Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) - 
     1                            (1/u_local(iGC2))
              else
                Seg_weight(iGC2) = 0.0
              endif
            else
              Seg_weight(iGC2) = (1/t_local(iGC2))*((ATAN(
     1                          (Seg_length(iGC2) - u_local(iGC2))/
     2                          t_local(iGC2))) - (ATAN(((-1)*
     3                          u_local(iGC2))/t_local(iGC2))))
            endif
              Seg_weight_t(iGC2) = Seg_weight(iGC2) * t_local(iGC2)
          enddo 

c       calculate the reciprocal of the sum of the segment weights, 
c       rec_Weight

          sum_Weight = 0.0
          do iGC2=1, n2-1
            sum_Weight = sum_Weight + Seg_weight(iGC2)
          enddo
            rec_Weight = 1/sum_Weight

c       calculate where you are on the strike of the fault for each
c       segment and segment weight * Seg_x + u_local

        Seg_x(1)=0.0
        do iGC2=2, n2-1
          Seg_x(iGC2) = Seg_x(iGC2-1) + Seg_length(iGC2-1)
        enddo
        
        do iGC2=1, n2-1
          Seg_wxu(iGC2) = Seg_weight(iGC2)*(Seg_x(iGC2)+u_local(iGC2))
        enddo
        
c       calculate Global Coordinate T
c       check for special case t=0 on segment

          sum_Swt = 0.0
          do iGC2=1, n2-1
            sum_Swt = sum_Swt + Seg_weight_t(iGC2)
            Global_T = rec_Weight*sum_Swt
          enddo  
          do iGC2=1, n2-1
            if (t_local(iGC2).EQ.0) then
              if (u_local(iGC2).GE.0) then
                if (u_local(iGC2).LE.Seg_length(iGC2)) then 
                Global_T = 0.0
                endif
              endif
            endif              
          enddo      

c       calculate Global Coordinate U
c       check for special case t=0 on segment

          sum_Swxu = 0.0
          do iGC2=1, n2-1
            sum_Swxu = sum_Swxu + Seg_wxu(iGC2)
            Global_U = rec_Weight*sum_Swxu
          enddo
          do iGC2=1, n2-1  
            if (t_local(iGC2).EQ.0) then
              if (u_local(iGC2).GE.0) then
                if (u_local(iGC2).LE.Seg_length(iGC2)) then 
                Global_U = Seg_x(iGC2)+u_local(iGC2)
                endif
              endif
            endif              
          enddo      
  
c CH: -------------- end
                                      
      n1 = nfltgrid(1)-n1AS(iLocX)+iLocY
      n2 = n2AS(iLocX)

C     For distance case set the hypocenter location based on user defined values. 
c      if (runflag .eq. 4) then
c         iLocY = int(nDDkm/stepY)
c         iLocX = int(nASkm/stepX)
c         if (iLocY .eq. 0) iLocY = 1
c         if (iLocX .eq. 0) iLocX = 1
c      endif

       ZTOR = fltGRid_z(iLocY,iLocX)

C      Compute average strike for given rupture area.
       strikeX = fltgrid_x(iLocY,n2) - fltgrid_x(iLocY,iLocX)
       strikeY = fltgrid_y(iLocY,n2) - fltgrid_y(iLocY,iLocX)
c       if (strikeY .eq. 0.0) then
c          strike = 0.0
c       else
c          strike = atan2(strikeY,strikeX)
c       endif
cnew
       if (strikeX .eq. 0.0) then
          strike = 0.0
       else
          strike = atan2(strikeX,strikeY)
       endif
      
C      Compute the Rx distance for each rupture area. 
C      First extend the two end points by 1000 km along the average strike for this rupture area
C      Site is assumed to be location: (0.0, 0.0, 0.0)
c       xtemp = 1000.0*cos(strike)
c       ytemp = 1000.0*sin(strike)
cnew
       xtemp = 1000.0*sin(strike)
       ytemp = 1000.0*cos(strike)
       xfltend1 = fltgrid_x(iLocY,iLocx) - xtemp
       xfltend2 = fltgrid_x(iLocY,iLocx) + xtemp
       yfltend1 = fltgrid_y(iLocY,iLocx) - ytemp
       yfltend2 = fltgrid_y(iLocY,iLocx) + ytemp

       call Calc_LineSeg_Dist (xfltend1,yfltend1, 0.0,
     1             xfltend2,yfltend2, 0.0, 0.0, 0.0, 0.0, distRx) 
       Rx = distRx

C     Now compute the Ry distance metric which is defined as the distance
C         along the strike mesured from the center of the rupture plane. 
C         For this calculation the site is assumed to be at a location of 0.0, 0.0.
       mDD = iLocY
       mAS = int((n2-iLocX)/2) + iLocX
       midDist = sqrt (fltgrid_x(iLocY,mAS)**2.0 + fltgrid_y(iLocY,mAS)**2.0 )
       if (Rx .lt. middist) then
          Ry = sqrt ( midDist**2.0 - Rx**2.0 ) 
       else
c  naa :  this hould never happen
          Ry = sqrt ( Rx**2.0 - midDist**2.0 ) 
       endif

c  CH: Mar 2015
c      compute the Ry0 (for a straight fault only)  Need to change this
       if (Ry.le.rupLen/2.) then
         Ry0 = 0.
       else
         Ry0 = Ry - rupLen/2.  
       endif         

C      Now determine if station is located on Hanging wall side or footwall side of fault rupture area. 
C      Reset HW or FW Flag for each rupture area.
       HWFlag = 0
       FWFlag = 0

C      Extend downdip point to check for site being on HW site.
C      First compute the distances and angles between upper and lower points on the rupture area. 
       adistX1 = fltgrid_x(n1,iLocX) - fltgrid_x(iLocY,iLocX)
       adistY1 = fltgrid_y(n1,iLocX) - fltgrid_y(iLocY,iLocX)
       adist1 = sqrt (adistX1*adistX1 + adistY1*adistY1)

C       if (astrike1 .eq. 0.0) then
C          astrike1 = 0.0
C       else
C          astrike1 = atan2(adistY1,adistX1)
C       endif 
cnew
       astrike1 = atan2(adistY1,adistX1)
    
c      Extend top of rupture end point locations by 100 additional km down dip
c             and end points by 1000 km along strike.
c      Set up testing points: 1-Upper left, 2-Upper right, 3-Lower left, 4-Lower right
       xtest(1) = xfltend1 
       ytest(1) = yfltend1 
       xtest(2) = xfltend2 
       ytest(2) = yfltend2 

c       xtest(3) = xfltend1 + cos(astrike1)*(adist1+100.0)
c       ytest(3) = yfltend1 + sin(astrike1)*(adist1+100.0)
c       xtest(4) = xfltend2 + cos(astrike1)*(adist1+100.0)
c       ytest(4) = yfltend2 + sin(astrike1)*(adist1+100.0)
cnew
       xtest(3) = xfltend1 + 100.0*cos(astrike1)*(adist1)
       ytest(3) = yfltend1 + 100.0*sin(astrike1)*(adist1)
       xtest(4) = xfltend2 + 100.0*cos(astrike1)*(adist1)
       ytest(4) = yfltend2 + 100.0*sin(astrike1)*(adist1)

       xtest(5) = xtest(1)
       ytest(5) = ytest(1)

c      Extend top of rupture end point locations by 100 additional km away from the down dip direction
c             and end points by 1000 km along strike
c      Set up testing points: 1-Upper left, 2-Upper right, 3-Lower left, 4-Lower right
       xtestfw(1) = xfltend1 
       ytestfw(1) = yfltend1 
       xtestfw(2) = xfltend2 
       ytestfw(2) = yfltend2 

       xtestfw(3) = xtestfw(1) - (xtest(3) - xtest(1) )
       ytestfw(3) = ytestfw(1) - (ytest(3) - ytest(1) )
       xtestfw(4) = xtestfw(2) - (xtest(4) - xtest(2) )
       ytestfw(4) = ytestfw(2) - (ytest(4) - ytest(2) )
       xtestfw(5) = xtestfw(1)
       ytestfw(5) = ytestfw(1)

C      Compute the average dip angle for the given rupture area. 
       if ( (fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)) .eq. 0.0) then
          dipavg = 3.14159/2.0
       else
          dipavg = atan2((fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)),sqrt(adistX1*adistX1+adistY1*adistY1))
       endif

C      Check to see if site is located over extended rupture area (HW). If so then Rx=Rrup.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtest, ytest, 0.0, 0.0, hwflag)

C      Check to see if site is located over extended rupture area (FW). If so then Rx=Rrup.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestfw, ytestfw, 0.0, 0.0, fwflag)

C      Loop over grid cell for given rupture area to find distances
C      Keep track of the fault grid cell for the closest rupture distance
C      for this rupture area since it is needed for the NGA directivity models.     
     
       do j=iLocX,n2
         do i=iLocY,n1
           if (distRup .gt. fltGrid_Rrup(i,j)) then
              distRup = fltGRid_rRup(i,j)
              icellRupstrike = j
              icellRupdip = i
           endif

           if (distJB .gt. fltGrid_Rjb(i,j)) then
             distJB = fltGRid_Rjb(i,j)
           endif

C     Compute the DistSeismo value based on min depth being equal to seismoDepth parameter.
           if (fltGrid_z(i,j) .gt. seismoDepth) then                 
              if (distSeismo .gt. fltGrid_Rrup(i,j)) then
                distSeismo = fltGRid_Rrup(i,j)
              endif
           endif 
         enddo

       enddo

C     Set the JB and Hypo distances for faults.
      distepi = fltGrid_Rjb(iLocY,iLocX)
      disthypo = fltGrid_Rrup(iLocY,iLocX)
      
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

      subroutine hanging_wall ( xRup, yRup, iSeg, x0, y0, dip0, hwflag)
c     implicit none
      real xRup(4,1), yRup(4,1) 
      real x(4), y(4), xSeg(5), ySeg(5)
      real dip0
      real x0, y0
      real pi, dip, theta, theta1, sin1, cos1 
      integer i, iSeg, nSeg, hwflag
      pi = 3.1415926

c       write (*,*) 'Dip = ', dip0

c     for (near) vertical fault, there is no hanging wall
      if (abs(dip0) .ge. 80.) return

      do i=1,4
         x(i) = xRup(i,iSeg)
         y(i) = yRup(i,iSeg)
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
      
c     Set boundary for points for which the closest point is to 
C     the plane
      xSeg(1) = x(1) 
      xSeg(2) = x(2) 
      xSeg(3) = x(3) + 100.*cos1
      xSeg(4) = x(4) + 100.*cos1
      xSeg(5) = xSeg(1)
      ySeg(1) = y(1) 
      ySeg(2) = y(2) 
      ySeg(3) = y(3) + 100.*sin1
      ySeg(4) = y(4) + 100.*sin1
      ySeg(5) = ySeg(1)

c     Determine if the site is inside this boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, hwflag)
      return
      end

c -------------------------------------------------------------------

      subroutine foot_wall ( xRup, yRup, iSeg, x0, y0, dip0, fwflag)
c     implicit none
      real xRup(4,1), yRup(4,1) 
      real x(4), y(4), xSeg(5), ySeg(5)
      real dip0
      real x0, y0
      real pi, dip, theta, theta1, sin1, cos1 
      integer i, iSeg, nSeg, fwflag
      pi = 3.1415926

c     for (near) vertical fault, there is no foot wall
      if (abs(dip0) .ge. 80.) return

      do i=1,4
         x(i) = xRup(i,iSeg)
         y(i) = yRup(i,iSeg)
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
      
c     Set boundary for points for which the closest point is to 
C     the plane
      xSeg(1) = x(1) 
      xSeg(2) = x(2) 
      xSeg(3) = x(3) - 100.*cos1
      xSeg(4) = x(4) - 100.*cos1
      xSeg(5) = xSeg(1)
      ySeg(1) = y(1) 
      ySeg(2) = y(2) 
      ySeg(3) = y(3) - 100.*sin1
      ySeg(4) = y(4) - 100.*sin1
      ySeg(5) = ySeg(1)

c     Determine if the site is inside this boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, fwflag)
      
      return
      end

c -------------------------------------------------------------------

      subroutine Get_plane_dist ( xRup, yRup, zRup, iSeg, 
     1           x0, y0, z0, dip0, dist )
c     implicit none
      real xRup(4,1), yRup(4,1), zRup(4,1)
      real x(4), y(4), z(4), xSeg(5), ySeg(5)
      real dist, dip0
      real x0, y0, z0, dist1
      real pi, dip, theta, theta1, sin1, cos1, 
     1     tan1, d1, d2
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

c ----------------------------------------------------

      subroutine CalcJB_dist ( xRup, yRup, iSeg, x0, y0, dip,
     1           distJB )
      real xRup(4,1), yRup(4,1), dip
      real x(4), y(4), xSeg(5), ySeg(5)
      
      do i=1,4
         x(i) = xRup(i,iSeg)
         y(i) = yRup(i,iSeg)
      enddo
      
c     Check for vertical dip
      if ( dip .eq. 90. ) then
         call Calc_LineSeg_dist ( x(1), y(1), 0., x(2), 
     1              y(2), 0., x0, y0, 0., distJB )
         return
      endif

      xSeg(1) = x(1)
      xSeg(2) = x(2)
      xSeg(3) = x(4)
      xSeg(4) = x(3)
      xSeg(5) = xSeg(1)
      ySeg(1) = y(1)
      ySeg(2) = y(2)
      ySeg(3) = y(4)
      ySeg(4) = y(3)
      ySeg(5) = ySeg(1)

c     Determine if the site is inside this boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, insideFlag )

c     Compute closest distance for inside
      if ( insideFlag .eq. 1 ) then
         distJB = 0.
         return
      endif
      
c     Compute the shortest dist to each edge
      distJB = 1.0e30
      do i=1,4
         i1 = i+1
         if ( i1 .gt. 4 ) i1=1
         call Calc_LineSeg_dist ( x(i), y(i), 0., x(i1), 
     1              y(i1), 0., x0, y0, 0., dist1 )
         if ( dist1 .lt. distJB ) then
            distJB = dist1
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
      integer nSeg, insideFlag
      real xSeg(1), ySeg(1), x0, y0
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
      real xFlt2(MAX_DD,1), yFlt2(MAX_DD,1), minDist, maxdist
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
      
      
c -------------------------------------------------------------------
      subroutine DirJWL (specT, Rrup, Rx, Ry, Ruplen, Mag, mech, 
     1        RupWidth, dip, HWFlag, medadj, sigadj )

      real Rx, Ry, Ruplen, Mag, medadj, sigadj, mech, per, Rrup,
     1     taperDist, taperMag, CosTheaAvg, Rymin, specT, RyRatio, RyRatioBar
      real period (25), c8b(25), c8bT
      real b0m, b1m, b2m, b3m, b4m, b5m, bMm, b6m, b1am      
      real b0s, b1s, b2s, b3s, b4s, b5s, bMs, b6s, b1as
      real rvb0m, rvb1m, rvb2m, rvb3m, rvb4m, rvb5m, rvbMm, rvb6m, rvb1am      
      real rvb0s, rvb1s, rvb2s, rvb3s, rvb4s, rvb5s, rvbMs, rvb6s, rvb1as
      real RupWidth, dip, cosphiAvg, RyRatioRV
      real thick, RxThick, RyThick, RyminThick 
      integer iflag, count1, count2, HWflag   

      data period     /
     1              0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500,
     1              0.0750, 0.1000, 0.1200, 0.1500, 0.1700,
     1              0.2000, 0.2500, 0.3000, 0.4000, 0.5000,
     1              0.7500, 1.0000, 1.5000, 2.0000, 3.0000,
     1              4.0000, 5.0000, 7.5000,10.0000/
      data c8b     /
     1             0.4833,  0.4833,  1.2144,  1.6421,  1.9456,  2.1810,
     1             2.6087,  2.9122,  3.1045,  3.3399,  3.4719,
     1             3.6434,  3.8787,  4.0711,  4.3745,  4.6099,
     1             5.0376,  5.3411,  5.7688,  6.0723,  6.5000,
     1             6.8035,  7.0389,  7.4666,  7.7700/

C     Mean Model Coefficients - StrikeSlip
      b0m = -0.142176807
      b1m = 0.227619199
      b2m = -0.056291119
      b3m = 0.993687135
      b4m = -0.234500608
      b5m = 0.63407793
      bMm = -0.28886311
      b6m = 0.020796267
      b1am = 0.347351648

C     Sigma Model Coefficients - StrikeSlip
      b0s = 0.114652061
      b1s = -0.079941683
      b2s = 0.215778271
      b3s = 0.612074481
      b4s = 1.365351827
      b5s = -4.669684236
      bMs = -0.254587556
      b6s = 0.020881559
      b1as = -0.03262467

C     Mean Model Coefficients - StrikeSlip
      rvb0m = -0.004584596
      rvb1m = 0.1600731578
      rvb2m = 0.1205763622
      rvb3m = -0.233024195
      rvb4m = 0.1153573318
      rvb5m = -0.149307015
      rvbMm = -0.271711183
      rvb6m = -0.000303303
      rvb1am = 0.1392414442

C     Sigma Model Coefficients - StrikeSlip
      rvb0s = 0.1697209219
      rvb1s = -0.024765316
      rvb2s = -0.021401123
      rvb3s = -0.324332252
      rvb4s = -0.929915749
      rvb5s = 2.2941530314
      rvbMs = -0.268399438
      rvb6s = 0.0139274036
      rvb1as = -0.073533106

C     First interpolate C8b term for given spectral period. 
C Find the requested spectral period and corresponding coefficients
      nPer = 25
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         c8bT  = c8b(1)
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
 1020       call interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
      
 1011 period1 = specT                                                                                                              

C     Now calculate the parameters for the model given input values.
      Rymin = min(Ry,Ruplen/2.0)
      costhetaAvg = (max((sqrt(Rx**2.0+Ry**2.0) - sqrt(Rx**2.0+(Ry-Rymin)**2.0))/Rymin,0.3)-0.5)

C     Compute for Reverse Case
C     Might need to be updated and confirmed with current model.
      Thick = Rupwidth/(sin(dip*3.14159/180.0))
      RyThick = Rx*sin(dip*3.14159/180.0)
      if (HWFlag .eq. 0) then
         RxThick = Rx*cos(dip*3.14159/180.0)
      else
         RxThick = -1.0*Rx*cos(dip*3.14159/180.0)
      endif
      RyminThick = min(RyThick,Thick/2.0)
      cosPhiAvg = (max((sqrt(RxThick**2.0+RyThick**2.0) - 
     1      sqrt(RxThick**2.0+(RyThick-RyminThick)**2.0))/RyminThick,0.3)-0.5)

      RyRatio = min(Ry/(Ruplen/2.0),1.0)
      RyRatioBar = RyRatio - (1.8709-0.1712*mag)      
      RyRatioBarRV = RyRatio - (2.1213-0.204*mag)

C     Compute the mag and distance tapers. 
      taperdist = max(1.0-max(Rrup-40.0,0.0)/30.0,0.0)
      taperMag = min(max(mag-5.5,0.0)/0.8,1.0)
 
C     Apply to Strike-Slip and Oblique slip cases
      if (mech .eq. 0.0 .or. mech .eq. 0.5) then      
C     Compute the mean adjustment
         medadj = taperDist*taperMag*((b0m + b1m*RyRatioBar + b1am*RyratioBar**2.0 + b2m*cosThetaAvg +
     1            b3m*costhetaAvg**2.0 + b4m*costhetaAvg**3.0 + b5m*costhetaAvg**4.0 + b6m*alog(Rrup/10.0) ) *
     1            exp(bMm*(mag-c8bT)**2.0) )
      
C     Compute the sigma adjustment
         sigadj = taperDist*taperMag*((b0s + b1s*RyRatioBar + b1as*RyratioBar**2.0 + b2s*cosThetaAvg +
     1            b3s*costhetaAvg**2.0 + b4s*costhetaAvg**3.0 + b5s*costhetaAvg**4.0 + b6s*alog(Rrup/10.0) ) *
     1            exp(bMs*(mag-c8bT)**2.0) )

C     Apply to Reverse cases
      elseif ( mech .eq. 1.0 ) then

C     Compute the mean adjustment
         medadj = taperDist*taperMag*((rvb0m + rvb1m*RyRatioBarRV + rvb1am*RyratioBarRV**2.0 + rvb2m*cosPhiAvg +
     1            rvb3m*cosPhiAvg**2.0 + rvb4m*cosPhiAvg**3.0 + rvb5m*cosPhiAvg**4.0 + rvb6m*alog(Rrup/10.0) ) *
     1            exp(rvbMm*(mag-c8bT)**2.0) )
      
C     Compute the sigma adjustment
         sigadj = taperDist*taperMag*((rvb0s + rvb1s*RyRatioBarRV + rvb1as*RyratioBarRV**2.0 + rvb2s*cosPhiAvg +
     1            rvb3s*cosPhiAvg**2.0 + rvb4s*cosPhiAvg**3.0 + rvb5s*cosPhiAvg**4.0 + rvb6s*alog(Rrup/10.0) ) *
     1            exp(rvbMs*(mag-c8bT)**2.0) )

c      write (*,*) Rx, Ry, ruplen, mag, Rrup, Rupwidth, Dip, Hwflag, medadj, Sigadj
c      write (*,*) cosPhiAvg, RyRatio, RyRAtioBarRV, RxThick, RyThick, Thick
c      pause


      endif

      return
      end
      


c -------------------------------------------------------------------
      subroutine DirJWL_V2 (specT, Rrup, Rx, Ry, Ruplen, Mag, mech, 
     1        RupWidth, dip, HWFlag, medadj, sigadj )

      real Rx, Ry, Ruplen, Mag, medadj, sigadj, mech, per, Rrup,
     1     taperDist, taperMag, CosTheaAvg, Rymin, specT, RyRatio, RyRatioBar
      real period (22), c8b(22), c8bT
      real c8rev(22), c8revT, c8org(22), c8orgT
      real b0m, b1m, b2m, b3m, b4m, b5m, b6m, b1am, bMm, sigm
      real b0s, b1s, b2s, b3s, b4s, b5s, b6s, b1as, bMs, sigs
      real MDF, SDF, MDA, SDA

      real b0_ssm(22), b1_ssm(22), b2_ssm(22), b3_ssm(22), b4_ssm(22), b5_ssm(22), bM_ssm(22)
      real b6_ssm(22), b1a_ssm(22), sigma_ssm(22)
      real b0_ssmT, b1_ssmT, b2_ssmT, b3_ssmT, b4_ssmT, b5_ssmT, bM_ssmT
      real b6_ssmT, b1a_ssmT, sigma_ssmT

      real b0_sss(22), b1_sss(22), b2_sss(22), b3_sss(22), b4_sss(22), b5_sss(22), bM_sss(22)
      real b6_sss(22), b1a_sss(22), sigma_sss(22)
      real b0_sssT, b1_sssT, b2_sssT, b3_sssT, b4_sssT, b5_sssT, bM_sssT
      real b6_sssT, b1a_sssT, sigma_sssT

      real b0_rvm(22), b1_rvm(22), b2_rvm(22), b3_rvm(22), b4_rvm(22), b5_rvm(22), bM_rvm(22)
      real b6_rvm(22), b1a_rvm(22), sigma_rvm(22)
      real b0_rvmT, b1_rvmT, b2_rvmT, b3_rvmT, b4_rvmT, b5_rvmT, bM_rvmT
      real b6_rvmT, b1a_rvmT, sigma_rvmT

      real b0_rvs(22), b1_rvs(22), b2_rvs(22), b3_rvs(22), b4_rvs(22), b5_rvs(22), bM_rvs(22)
      real b6_rvs(22), b1a_rvs(22), sigma_rvs(22)
      real b0_rvsT, b1_rvsT, b2_rvsT, b3_rvsT, b4_rvsT, b5_rvsT, bM_rvsT
      real b6_rvsT, b1a_rvsT, sigma_rvsT
      
      real sigchange(22), sigchangeT

      real RupWidth, dip, cosphiAvg, RyRatioRV
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

C     Mean Model Coefficients - StrikeSlip
      data b0_ssm  /  -0.142176807, -0.142176807, -0.142176807, -0.142176807, -0.142176807, 
     1        -0.142176807, -0.142176807, -0.142176807, -0.142176807, -0.142176807, -0.142176807,
     1        -0.142176807, -0.142176807, -0.142176807, -0.142176807, -0.142176807, -0.142176807, 
     1        -0.142176807, -0.142176807, -0.142176807, -0.142176807, -0.142176807  / 
      data b1_ssm  /  0.227619199, 0.227619199, 0.227619199, 0.227619199, 0.227619199, 0.227619199,
     1         0.227619199, 0.227619199, 0.227619199, 0.227619199, 0.227619199, 0.227619199, 
     1         0.227619199, 0.227619199, 0.227619199, 0.227619199, 0.227619199, 0.227619199, 
     1         0.227619199, 0.227619199, 0.227619199, 0.227619199  / 
      data b2_ssm  /  -0.056291119, -0.056291119, -0.056291119, -0.056291119, -0.056291119, 
     1        -0.056291119, -0.056291119, -0.056291119, -0.056291119, -0.056291119, -0.056291119, 
     1        -0.056291119, -0.056291119, -0.056291119, -0.056291119, -0.056291119, -0.056291119, 
     1        -0.056291119, -0.056291119, -0.056291119, -0.056291119, -0.056291119  / 
      data b3_ssm  /  0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135,
     1         0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135,
     1         0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135, 0.993687135,
     1         0.993687135, 0.993687135  / 
      data b4_ssm  /  -0.234500608, -0.234500608, -0.234500608, -0.234500608, -0.234500608, 
     1        -0.234500608, -0.234500608, -0.234500608, -0.234500608, -0.234500608, -0.234500608,
     1        -0.234500608, -0.234500608, -0.234500608, -0.234500608, -0.234500608, -0.234500608,
     1        -0.234500608, -0.234500608, -0.234500608, -0.234500608, -0.234500608  / 
      data b5_ssm  /  0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793, 
     1         0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793,
     1         0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793, 0.63407793,
     1         0.63407793, 0.63407793  / 
      data bM_ssm  /  -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311,
     1        -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311,
     1        -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311, -0.28886311,
     1        -0.28886311, -0.28886311  / 
      data b6_ssm /  0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267, 
     1         0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267,
     1         0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267, 0.020796267,
     1         0.020796267, 0.020796267  / 
      data b1a_ssm  /  0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648,
     1         0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648,
     1         0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648, 0.347351648,
     1         0.347351648, 0.347351648  / 
      data Sigma_ssm  /  0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429,
     1         0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429,
     1         0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429, 0.0446429  / 

C     Sigma Model Coefficients - StrikeSlip
      data b0_sss  /  0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061,
     1         0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061,
     1         0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061, 0.114652061, 
     1         0.114652061, 0.114652061  / 
      data b1_sss  /  -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, 
     1        -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, 
     1        -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, -0.079941683, 
     1        -0.079941683, -0.079941683  / 
      data b2_sss  /  0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271,
     1         0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271,
     1         0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271, 0.215778271  / 
      data b3_sss  /  0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481,
     1         0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481,
     1         0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481, 0.612074481  / 
      data b4_sss  /  1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827,
     1         1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827,
     1         1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827, 1.365351827 / 
      data b5_sss  /  -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236,
     1        -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, 
     1        -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, -4.669684236, 
     1        -4.669684236  / 
      data bM_sss  /  -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556, 
     1        -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556,
     1        -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556, -0.254587556,
     1        -0.254587556, -0.254587556  / 
      data b6_sss  /  0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559, 
     1         0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559,
     1         0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559, 0.020881559,
     1         0.020881559, 0.020881559  / 
      data b1a_sss  /  -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467,
     1        -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467,
     1        -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467, -0.03262467,
     1        -0.03262467, -0.03262467  / 
      data Sigma_sss  /  0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 
     1         0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027, 0.031027,
     1         0.031027, 0.031027, 0.031027, 0.031027, 0.031027  / 

C     Mean Model Coefficients - Reverse
      data b0_rvm  /  -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596,
     1        -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596,
     1        -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596, -0.004584596,
     1        -0.004584596, -0.004584596  / 
      data b1_rvm  /  0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158, 
     1         0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158,
     1         0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158, 0.160073158,
     1         0.160073158, 0.160073158 / 
      data b2_rvm  /  0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362, 
     1         0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362,
     1         0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362, 0.120576362,
     1         0.120576362, 0.120576362  / 
      data b3_rvm  /  -0.233024195, -0.233024195, -0.233024195, -0.233024195, -0.233024195, -0.233024195,
     1        -0.233024195, -0.233024195, -0.233024195, -0.233024195, -0.233024195, -0.233024195, 
     1        -0.233024195, -0.233024195, -0.233024195, -0.233024195, -0.233024195, -0.233024195, 
     1        -0.233024195, -0.233024195, -0.233024195, -0.233024195  / 
      data b4_rvm  /  0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332,
     1         0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332,
     1         0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332, 0.115357332,
     1         0.115357332, 0.115357332  / 
      data b5_rvm  /  -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015,
     1        -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015,
     1        -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015, -0.149307015,
     1        -0.149307015, -0.149307015  / 
      data bM_rvm  /  -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183,
     1        -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183,
     1        -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183, -0.271711183,
     1        -0.271711183, -0.271711183  / 
      data b6_rvm  /  -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303, 
     1        -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303,
     1        -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303, -0.000303303,
     1        -0.000303303, -0.000303303  / 
      data b1a_rvm  /  0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444,
     1         0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444,
     1         0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444, 0.139241444,
     1         0.139241444, 0.139241444  / 
      data Sigma_rvm  /  0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932,
     1         0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932,
     1         0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932, 0.0302932  / 

C     Sigma Model Coefficients - Reverse
      data b0_rvs  /  0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 
     1         0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 
     1         0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922, 0.169720922,
     1         0.169720922, 0.169720922  / 
      data b1_rvs  /  -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316,
     1        -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316,
     1        -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316, -0.024765316,
     1        -0.024765316, -0.024765316  / 
      data b2_rvs  /  -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123,
     1        -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123,
     1        -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123, -0.021401123,
     1        -0.021401123, -0.021401123  / 
      data b3_rvs  /  -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252, 
     1        -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252,
     1        -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252, -0.324332252,
     1        -0.324332252, -0.324332252  / 
      data b4_rvs  /  -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749, 
     1        -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749,
     1        -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749, -0.929915749,
     1        -0.929915749, -0.929915749  / 
      data b5_rvs  /  2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031,
     1         2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031,
     1         2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031, 2.294153031, 
     1         2.294153031, 2.294153031  / 
      data bM_rvs  /  -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438,
     1        -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438,
     1        -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438, -0.268399438,
     1        -0.268399438, -0.268399438  / 
      data b6_rvs  /  0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404,
     1         0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404,
     1         0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404, 0.013927404,
     1         0.013927404, 0.013927404  / 
      data b1a_rvs  /  -0.073533106, -0.073533106, -0.073533106, -0.073533106, -0.073533106,
     1        -0.073533106, -0.073533106, -0.073533106, -0.073533106, -0.073533106, -0.073533106,
     1        -0.073533106, -0.073533106, -0.073533106, -0.073533106, -0.073533106, -0.073533106,
     1        -0.073533106, -0.073533106, -0.073533106, -0.073533106, -0.073533106  / 
      data Sigma_rvs  /  0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974,
     1         0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974,
     1         0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974, 0.0163974  / 

      data sigchange / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1         0.0, 0.0, 0.1205, 0.1535, 0.1913, 0.2285, 0.2478 /

      data c8org / 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 
     1         0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 /

C     First interpolate terms for given spectral period. 
C Find the requested spectral period and corresponding coefficients
      nPer = 22
C First check for the PGA case (i.e., specT=0.0) 
      if (specT .eq. 0.0) then
         c8revT  = c8rev(1)
         c8bT  = c8b(1)
         c8orgT  = c8org(1)
         sigChangeT = sigchange(1)
         b0_ssmT = b0_ssm(1)
         b1_ssmT = b1_ssm(1)
         b2_ssmT = b2_ssm(1)
         b3_ssmT = b3_ssm(1)
         b4_ssmT = b4_ssm(1)
         b5_ssmT = b5_ssm(1)
         bM_ssmT = bM_ssm(1)
         b6_ssmT = b6_ssm(1)
         b1a_ssmT = b1a_ssm(1)
         sigma_ssmT = sigma_ssm(1)

         b0_sssT = b0_sss(1)
         b1_sssT = b1_sss(1)
         b2_sssT = b2_sss(1)
         b3_sssT = b3_sss(1)
         b4_sssT = b4_sss(1)
         b5_sssT = b5_sss(1)
         bM_sssT = bM_sss(1)
         b6_sssT = b6_sss(1)
         b1a_sssT = b1a_sss(1)
         sigma_sssT = sigma_sss(1)

         b0_rvmT = b0_rvm(1)
         b1_rvmT = b1_rvm(1)
         b2_rvmT = b2_rvm(1)
         b3_rvmT = b3_rvm(1)
         b4_rvmT = b4_rvm(1)
         b5_rvmT = b5_rvm(1)
         bM_rvmT = bM_rvm(1)
         b6_rvmT = b6_rvm(1)
         b1a_rvmT = b1a_rvm(1)
         sigma_rvmT = sigma_rvm(1)

         b0_rvsT = b0_rvs(1)
         b1_rvsT = b1_rvs(1)
         b2_rvsT = b2_rvs(1)
         b3_rvsT = b3_rvs(1)
         b4_rvsT = b4_rvs(1)
         b5_rvsT = b5_rvs(1)
         bM_rvsT = bM_rvs(1)
         b6_rvsT = b6_rvs(1)
         b1a_rvsT = b1a_rvs(1)
         sigma_rvsT = sigma_rvs(1)

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
            call interp (period(count1),period(count2),sigchange(count1),sigchange(count2),
     +                   specT,sigchangeT,iflag)
     
            call interp (period(count1),period(count2),b0_ssm(count1),b0_ssm(count2),
     +                   specT,b0_ssmT,iflag)
            call interp (period(count1),period(count2),b1_ssm(count1),b1_ssm(count2),
     +                   specT,b1_ssmT,iflag)
            call interp (period(count1),period(count2),b2_ssm(count1),b2_ssm(count2),
     +                   specT,b2_ssmT,iflag)
            call interp (period(count1),period(count2),b3_ssm(count1),b3_ssm(count2),
     +                   specT,b3_ssmT,iflag)
            call interp (period(count1),period(count2),b4_ssm(count1),b4_ssm(count2),
     +                   specT,b4_ssmT,iflag)
            call interp (period(count1),period(count2),b5_ssm(count1),b5_ssm(count2),
     +                   specT,b5_ssmT,iflag)
            call interp (period(count1),period(count2),bM_ssm(count1),bM_ssm(count2),
     +                   specT,bM_ssmT,iflag)
            call interp (period(count1),period(count2),b6_ssm(count1),b6_ssm(count2),
     +                   specT,b6_ssmT,iflag)
            call interp (period(count1),period(count2),b1a_ssm(count1),b1a_ssm(count2),
     +                   specT,b1a_ssmT,iflag)
            call interp (period(count1),period(count2),sigma_ssm(count1),sigma_ssm(count2),
     +                   specT,sigma_ssmT,iflag)

            call interp (period(count1),period(count2),b0_sss(count1),b0_sss(count2),
     +                   specT,b0_sssT,iflag)
            call interp (period(count1),period(count2),b1_sss(count1),b1_sss(count2),
     +                   specT,b1_sssT,iflag)
            call interp (period(count1),period(count2),b2_sss(count1),b2_sss(count2),
     +                   specT,b2_sssT,iflag)
            call interp (period(count1),period(count2),b3_sss(count1),b3_sss(count2),
     +                   specT,b3_sssT,iflag)
            call interp (period(count1),period(count2),b4_sss(count1),b4_sss(count2),
     +                   specT,b4_sssT,iflag)
            call interp (period(count1),period(count2),b5_sss(count1),b5_sss(count2),
     +                   specT,b5_sssT,iflag)
            call interp (period(count1),period(count2),bM_sss(count1),bM_sss(count2),
     +                   specT,bM_sssT,iflag)
            call interp (period(count1),period(count2),b6_sss(count1),b6_sss(count2),
     +                   specT,b6_sssT,iflag)
            call interp (period(count1),period(count2),b1a_sss(count1),b1a_sss(count2),
     +                   specT,b1a_sssT,iflag)
            call interp (period(count1),period(count2),sigma_sss(count1),sigma_sss(count2),
     +                   specT,sigma_sssT,iflag)

            call interp (period(count1),period(count2),b0_rvm(count1),b0_rvm(count2),
     +                   specT,b0_rvmT,iflag)
            call interp (period(count1),period(count2),b1_rvm(count1),b1_rvm(count2),
     +                   specT,b1_rvmT,iflag)
            call interp (period(count1),period(count2),b2_rvm(count1),b2_rvm(count2),
     +                   specT,b2_rvmT,iflag)
            call interp (period(count1),period(count2),b3_rvm(count1),b3_rvm(count2),
     +                   specT,b3_rvmT,iflag)
            call interp (period(count1),period(count2),b4_rvm(count1),b4_rvm(count2),
     +                   specT,b4_rvmT,iflag)
            call interp (period(count1),period(count2),b5_rvm(count1),b5_rvm(count2),
     +                   specT,b5_rvmT,iflag)
            call interp (period(count1),period(count2),bM_rvm(count1),bM_rvm(count2),
     +                   specT,bM_rvmT,iflag)
            call interp (period(count1),period(count2),b6_rvm(count1),b6_rvm(count2),
     +                   specT,b6_rvmT,iflag)
            call interp (period(count1),period(count2),b1a_rvm(count1),b1a_rvm(count2),
     +                   specT,b1a_rvmT,iflag)
            call interp (period(count1),period(count2),sigma_rvm(count1),sigma_rvm(count2),
     +                   specT,sigma_rvmT,iflag)

            call interp (period(count1),period(count2),b0_rvs(count1),b0_rvs(count2),
     +                   specT,b0_rvsT,iflag)
            call interp (period(count1),period(count2),b1_rvs(count1),b1_rvs(count2),
     +                   specT,b1_rvsT,iflag)
            call interp (period(count1),period(count2),b2_rvs(count1),b2_rvs(count2),
     +                   specT,b2_rvsT,iflag)
            call interp (period(count1),period(count2),b3_rvs(count1),b3_rvs(count2),
     +                   specT,b3_rvsT,iflag)
            call interp (period(count1),period(count2),b4_rvs(count1),b4_rvs(count2),
     +                   specT,b4_rvsT,iflag)
            call interp (period(count1),period(count2),b5_rvs(count1),b5_rvs(count2),
     +                   specT,b5_rvsT,iflag)
            call interp (period(count1),period(count2),bM_rvs(count1),bM_rvs(count2),
     +                   specT,bM_rvsT,iflag)
            call interp (period(count1),period(count2),b6_rvs(count1),b6_rvs(count2),
     +                   specT,b6_rvsT,iflag)
            call interp (period(count1),period(count2),b1a_rvs(count1),b1a_rvs(count2),
     +                   specT,b1a_rvsT,iflag)
            call interp (period(count1),period(count2),sigma_rvs(count1),sigma_rvs(count2),
     +                   specT,sigma_rvsT,iflag)

 1011 period1 = specT                                                                                                              

C     Now calculate the parameters for the model given input values.

C     First set coefficients b0-b6 based on source mechanism
C     Classification is strike-slip or normal --> Strike-slip coefficients
C                       reverse or oblique    --> Reverse coefficients
      if (mech .le. 0.0) then
          b0m = b0_ssmT
          b1m = b1_ssmT
          b2m = b2_ssmT
          b3m = b3_ssmT
          b4m = b4_ssmT
          b5m = b5_ssmT
          b6m = b6_ssmT
          b1am = b1a_ssmT
          bMm = bM_ssmT
          sigm = sigma_ssmT

          b0s = b0_sssT
          b1s = b1_sssT
          b2s = b2_sssT
          b3s = b3_sssT
          b4s = b4_sssT
          b5s = b5_sssT
          b6s = b6_sssT
          b1as = b1a_sssT
          bMs = bM_sssT
          sigs = sigma_sssT
      else
          b0m = b0_rvmT
          b1m = b1_rvmT
          b2m = b2_rvmT
          b3m = b3_rvmT
          b4m = b4_rvmT
          b5m = b5_rvmT
          b6m = b6_rvmT
          b1am = b1a_rvmT
          bMm = bM_rvmT
          sigm = sigma_rvmT

          b0s = b0_rvsT
          b1s = b1_rvsT
          b2s = b2_rvsT
          b3s = b3_rvsT
          b4s = b4_rvsT
          b5s = b5_rvsT
          b6s = b6_rvsT
          b1as = b1a_rvsT
          bMs = bM_rvsT
          sigs = sigma_rvsT
      endif

C     Now compute the Directivity factors and adjustments.

      RyRatio = min(Ry/(Ruplen/2.0),1.0)
      RyRatioBar = RyRatio - (1.8709-0.1712*mag)
      costhetaAvg = max( ((sqrt(Rx**2.0+Ry**2.0) - sqrt(Rx**2.0 + (Ry-min(Ry,RupLen/2.0))**2.0) )/min(Ry,Ruplen/2.0)),0.3) - 0.5

      cosdip = cos(dip*3.14159/180.0)
      sindip = sin(dip*3.14159/180.0)

      term1 = sqrt ( (Rx*sindip)**2.0 + abs((-1.0*Rx*cosdip+Rupwidth/2.0)**2.0)) 
      term2 = -1.0*sqrt ( (Rx*sindip)**2.0 + (abs(-1.0*Rx*cosdip+RupWidth/2.0) - 
     1            min(abs(-1.0*Rx*cosdip+RupWidth/2.0),RupWidth/2.0))**2.0  ) 
      term3 = min(abs(-1.0*Rx*cosdip+RupWidth/2.0),RupWidth/2.0)
      cosPhiAvg = max( (term1+term2)/term3,0.3) - 0.5

C     Compute the mag and distance tapers. 
      taperdist = max(1.0-max(Rrup-40.0,0.0)/30.0,0.0)
      taperMag = min(max(mag-5.5,0.0)/0.8,1.0)
 

      if (mech .le. 0.0) then
         MDF = b0m + b1m*RyRatioBar + b1am*RyratioBar*RyRatioBar + b2m*costhetaAvg + b3m*costhetaAvg**2.0 +
     1         b4m*costhetaAvg**3.0 + b5m*costhetaAvg**4.0 + b6m*alog(Rrup/10.0)
         SDF = b0s + b1s*RyRatioBar + b1as*RyratioBar*RyRatioBar + b2s*costhetaAvg + b3s*costhetaAvg**2.0 +
     1         b4s*costhetaAvg**3.0 + b5s*costhetaAvg**4.0 + b6s*alog(Rrup/10.0)     
      else
         MDF = b0m + b1m*RyRatioBar + b1am*RyratioBar*RyRatioBar + b2m*cosPhiAvg + b3m*cosPhiAvg**2.0 +
     1         b4m*cosPhiAvg**3.0 + b5m*costhetaAvg**4.0 + b6m*alog(Rrup/10.0)
         SDF = b0s + b1s*RyRatioBar + b1as*RyratioBar*RyRatioBar + b2s*cosPhiAvg + b3s*cosPhiAvg**2.0 +
     1         b4s*cosPhiAvg**3.0 + b5s*cosPhiAvg**4.0 + b6s*alog(Rrup/10.0)           
      endif

C     Now compute the median and sigma adjustments

      Medadj = MDF*exp(bMm*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag
      SDA = SDF*exp(bMs*(mag-c8bT)**2.0)*(c8revT/C8orgT)*taperdist*tapermag

      if (sigchangeT .le. SDA) then
         Sigadj = sqrt (SDA*SDA - sigchangeT*sigchangeT)
      else
         Sigadj = 0.0
      endif

      return
      end
      






c -------------------------------------------------------------------
      subroutine DirJWL_V3 (specT, Rrup, Rx, Ry, Ruplen, Mag, mech, 
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
      b0_ssm = -0.14556
      b1_ssm = 0.036602
      b2_ssm = 0.55010
      b3_ssm = -0.14676
      bM_ssm = -0.27145

C     Sigma Model Coefficients - StrikeSlip
      b0_sss = 0.045097
      b1_sss = 0.16267
      b2_sss = 0.20112
      b3_sss = -0.37626
      bM_sss = -0.10699

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

      return
      end
      

