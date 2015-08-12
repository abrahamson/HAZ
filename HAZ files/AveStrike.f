       subroutine AveStrike (fltgrid_x, fltgrid_y, fltgrid_z, fltgrid_x1,
     1                       fltgrid_y1, fltgrid_x2, fltgrid_y2, fltgrid_x4,
     2                       fltgrid_y4, rupLen, distRx, iLocX, iLocY, n2, 
     3                       n1, Rx, Ry, Ry0, HWFlag, dipavg)
      
       implicit none
       include 'pfrisk.h'
       
c      declarations passed in 
       integer iLocX, iLocY, n2, n1
       real fltgrid_x(MAXFLT_DD,MAXFLT_AS), fltgrid_y(MAXFLT_DD,MAXFLT_AS), 
     1      fltgrid_z(MAXFLT_DD,MAXFLT_AS), fltgrid_x1(MAXFLT_DD,MAXFLT_AS),
     2      fltgrid_y1(MAXFLT_DD,MAXFLT_AS), fltgrid_x2(MAXFLT_DD,MAXFLT_AS),
     3      fltgrid_y2(MAXFLT_DD,MAXFLT_AS), fltgrid_x4(MAXFLT_DD,MAXFLT_AS),
     4      fltgrid_y4(MAXFLT_DD,MAXFLT_AS), rupLen, distRx        
       
c      declarations passed out
       integer HWFlag
       real Rx, Ry, Ry0, dipavg

c      declarations only used within subroutine
       integer mDD, mAS, FWFlag 
       real StrikeX, StrikeY, strike, xtemp, ytemp, xfltend1, xfltend2,
     1      yfltend1, yfltend2, midDist, adistX1, adistY1, adist1, astrike1,
     2      xtest(5), ytest(5), xtestfw(5), ytestfw(5) 

C      Compute average strike for given rupture area.
       strikeX = fltgrid_x2(iLocY,n2) - fltgrid_x1(iLocY,iLocX)
       strikeY = fltgrid_y2(iLocY,n2) - fltgrid_y1(iLocY,iLocX)
       if (strikeX .eq. 0.0) then
          strike = 0.0
       else
          strike = atan2(strikeX,strikeY)
       endif
      
C      Compute the Rx distance for each rupture area. 
C      First extend the two end points by 1000 km along the average strike for this rupture area
C      Site is assumed to be location: (0.0, 0.0, 0.0)

       xtemp = 1000.0*sin(strike)
       ytemp = 1000.0*cos(strike)
       xfltend1 = fltgrid_x1(iLocY,iLocx) - xtemp
       xfltend2 = fltgrid_x1(iLocY,iLocx) + xtemp
       yfltend1 = fltgrid_y1(iLocY,iLocx) - ytemp
       yfltend2 = fltgrid_y1(iLocY,iLocx) + ytemp     

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
c  naa :  this should never happen
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

C      Extend downdip point to check for site being on HW side.
C      First compute the distances and angles between upper and lower points on the rupture area. 
       adistX1 = fltgrid_x4(n1,iLocX) - fltgrid_x1(iLocY,iLocX)
       adistY1 = fltgrid_y4(n1,iLocX) - fltgrid_y1(iLocY,iLocX)
       adist1 = sqrt (adistX1*adistX1 + adistY1*adistY1)

       astrike1 = atan2(adistY1,adistX1)     
    
c      Extend top of rupture end point locations by 100 additional km down dip
c             and end points by 1000 km along strike.
c      Set up testing points: 1-Upper left, 2-Upper right, 3-Lower left, 4-Lower right
       xtest(1) = xfltend1 
       ytest(1) = yfltend1 
       xtest(2) = xfltend2 
       ytest(2) = yfltend2 

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
       call Inside_OutSide ( 4, xtest, ytest, 0.0, 0.0, HWFlag)

C      Check to see if site is located over extended rupture area (FW). If so then Rx=Rrup.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestfw, ytestfw, 0.0, 0.0, FWFlag)
       
        return
       end
