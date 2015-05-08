      
c -------------------------------------------------------------

      subroutine SetRupLoc ( hx, hy, rupArea, rupWidth0,
     1           faultWidth, faultLength, cumLength, strike,
     2           nPts, xFlt, yFlt, zFlt, dip, xRup, yRup, zRup,
     3           dipRup, nSegRup, minDepth, hypodepth )

      real xFlt(1), yFlt(1), zFlt(1), dip(1)
      real xRup(4,1), yRup(4,1), zRup(4,1), dipRup(1)
  
      return
      end

c ----------------------------------------------------

      subroutine ConvertCoordinates ( fLat, fLong, fZ, siteX, siteY,
     1         xFlt, yFlt, zFlt, nfp, iFlt, nPts1, x0, y0,
     2         z0, iCoor, grid_lat, grid_long, grid_x, grid_y,
     3   grid_dlat, grid_dlong, grid_dx, grid_dy, grid_n, sourceType, nDD )
     
      include 'pfrisk.h'

      real fLat(MAX_FLT,MAX_DD,MAX_SEG), fLong(MAX_FLT,MAX_DD,MAX_SEG), fZ(MAX_FLT,MAX_DD,MAX_SEG),
     1     dip(MAX_FLT,1)
      real xFlt(MAX_DD,1), yFlt(MAX_DD,1), zFlt(MAX_DD,1), dip1(1)
      integer nfp(1), sourceType(1), grid_n(MAX_FLT), nDD(1)
      real grid_lat(MAX_FLT,MAX_GRID)
      real grid_long(MAX_FLT,MAX_GRID)
      real grid_dlong(MAX_FLT), grid_dlat(MAX_FLT)
      real grid_x(MAX_GRID), grid_y(MAX_GRID), grid_dx, grid_dy, conlat
      
            
c     Load npts in fault and dip into 1-D arrays
      npts1 = nfp(iFlt)
      
      if ( iCoor .eq. 1. ) then
c     Convert Lat, long to km (using site coordinates as ref)
      refLat = siteY
      refLong = siteX
C     Use average latitude of fault to scale flat earth conversion. 
C     Note only use for souretype = 1 or 5. For other sourcetypes 
C     (i.e., areal and gridded sources) use site latitude for conversion. 
      if (sourcetype(iFlt) .eq. 1 ) then
         conlat = (fLat(iFlt,1,nfp(iFlt)) + fLat(iFlt,1,1) ) / 2.0
      elseif (sourcetype(iFlt) .eq. 5 ) then
         conlat = (fLat(iFlt,1,nfp(iFlt)) + fLat(iFlt,1,1) ) / 2.0
      elseif (sourcetype(iFlt) .eq. 6 ) then
         conlat = (fLat(iFlt,1,nfp(iFlt)) + fLat(iFlt,1,1) ) / 2.0
      else
         conlat = refLat
      endif

      xScale = 111.12 * cos(conlat/180*3.1415926)
      yScale = 111.12

      else
c     Units are in km already 
      reflat = 0.
      refLong = 0.
      xScale = 1.
      yScale = 1.
      endif

c     Check for grid source
      if ( sourceType(iFlt) .eq. 3 .or. sourceType(iFlt) .eq. 4 ) then
        do i=1,grid_n(iFlt)
          grid_y(i) = (grid_lat(iFlt,i) - reflat ) * yScale
          grid_x(i) = (grid_long(iFlt,i) - reflong ) * xScale
        enddo
        grid_dx = grid_dlong(iFlt)*xscale 
        grid_dy = grid_dlat(iFlt)*yscale
      else

        do iz=1,nDD(iFlt)
          do iPt=1,nfp(iFlt)
            xFlt(iz,iPt) = ( fLong(iFlt,iz,iPt) - reflong ) * xScale
            yFlt(iz,iPt) = ( fLat(iFlt,iz,iPt) - reflat ) * yScale
            zFlt(iz,iPt) = fZ(iFlt,iz,iPt)            
          enddo
        enddo
      endif
      x0 = ( siteX - reflong ) * xScale
      y0 = ( siteY - reflat ) * yScale
      z0 = 0.

      return 
      end

c ---------------------------------------------------------------------

      subroutine calcFltGrid ( xFlt, yFlt, zFlt, npts, nDD, fltGrid_X, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGRid_Rrup, FltGRid_Rjb, fltGrid_Rseis, fltGrid_Rx, FltGrid_HW,
     3               faultArea, faultLen, faultW, distMin1, step, FltGrid_fLen )  
     
      include 'pfrisk.h'
      include 'declare1.h'
cnjg      integer nDD, n2(MAX_DD), nfltGrid(2), n1(MAX_SEG)
      integer nDD, n2(MAXFLT_DD), nfltGrid(2), n1(MAXFLT_AS)
      real fltGrid_X(MAXFLT_DD,MAXFLT_AS), fltGrid_y(MAXFLT_DD,MAXFLT_AS), 
     1     fltGrid_z(MAXFLT_DD,MAXFLT_AS),fltGrid_fLen(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_w(MAXFLT_DD,MAXFLT_AS), fltGrid_a(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_RJB(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_Rseis(MAXFLT_DD,MAXFLT_AS), fltGrid_Rx(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_HW(MAXFLT_DD,MAXFLT_AS)

c     set grid step in km
c      step = 0.25
cnjg      step = 0.25  

      distMin1 = 1.E10

c     Find the widest part of the fault and set the number of grid points down dip
      sum1 = 0.
      nn = 0
      do i=1,nDD-1
        width1 = 0.
        do j=1,npts
          dx = xFlt(i+1,j) - xFlt(i,j)
          dy = yFlt(i+1,j) - yFlt(i,j)
          dz = zFlt(i+1,j) - zFlt(i,j)
          w = sqrt( dx**2 + dy**2 + dz**2 ) 
          if ( w .gt. width1 ) width1 = w
cnjg          sum1 = sum1 + w
          sum1 = sum1 + real(w/npts)
          nn = nn + 1
        enddo
        n2(i) = Int( width1 / step ) + 1
      enddo
cnjg      faultW = sum1 / nn
      faultW = sum1 
      nx1 = 0
      ny1 = 0
      ii = 0
      i0 = 0
      j0 = 0
      faultArea = 0
      faultLen = 0
      do i=1,nDD-1
        ny1 = ny1 + n2(i)

        do j=2,npts
          if ( i .eq. 1 ) then
            xLtop =  xFlt(i,j) - xFlt(i,j-1)
            yLtop =  yFlt(i,j) - yFlt(i,j-1)
            n1(j) = int(sqrt(xLtop**2 + yLtop**2)/step) + 1
            nx1 = nx1 + n1(j)
            faultLen = faultLen + sqrt(xLtop**2 + yLtop**2)
          endif
          
          xLtop =  xFlt(i,j) - xFlt(i,j-1)
          yLtop =  yFlt(i,j) - yFlt(i,j-1)
          zLtop =  zFlt(i,j) - zFlt(i,j-1)
          xLbot =  xFlt(i+1,j) - xFlt(i+1,j-1)
          yLbot =  yFlt(i+1,j) - yFlt(i+1,j-1)
          zLbot =  zFlt(i+1,j) - zFlt(i+1,j-1)
          
          do j1=1,n1(j)
            j2 = j1 + j0
            x1 = xFlt(i,j-1) + (XLtop/n1(j))*(j1-1+0.5)
            y1 = yFlt(i,j-1) + (yLtop/n1(j))*(j1-1+0.5)
            z1 = zFlt(i,j-1) + (zLtop/n1(j))*(j1-1+0.5)
            x2 = xFlt(i+1,j-1) + (XLbot/n1(j))*(j1-1+0.5)
            y2 = yFlt(i+1,j-1) + (yLbot/n1(j))*(j1-1+0.5)
            z2 = zFlt(i+1,j-1) + (zLbot/n1(j))*(j1-1+0.5)
            
            dx = (x2 - x1) / n2(i)
            dy = (y2 - y1) / n2(i)
            dz = (z2 - z1) / n2(i)

c           Set coordinates of center of grid and width and area of grid
            do i1=1,n2(i)
              i2 = i1 + i0
              fltGrid_x(i2,j2) = x1 + dx*(i1-1+0.5)
              fltGrid_y(i2,j2) = y1 + dy*(i1-1+0.5)
              fltGrid_z(i2,j2) = z1 + dz*(i1-1+0.5)
              

c   *** this needs to be fixed.  It now assumes right angles ***
              fltGrid_w(i2,j2) = sqrt(dx**2+dy**2+dz**2)
              fltGrid_a(i2,j2) = fltGrid_w(i2,j2) * sqrt( (xLtop/n1(j))**2 + (ylTop/n1(j))**2 )
              faultArea = faultArea + fltGrid_a(i2,j2)  

            enddo
            i0 = ii
          enddo
          j0 = j0 + n1(j)
        enddo
        ii = ii + n2(i)
        i0 = ii
        j0 = 0
      enddo   
      nfltGrid(1) = ny1
      nfltGrid(2) = nx1
      
c     Find the distance to each grid point. 
      do i=1,ny1
        do j=1,nx1
          dx = fltGrid_x(i,j) - x0
          dy = fltGrid_y(i,j) - y0
          dz = fltGrid_z(i,j) - z0
          FltGrid_Rrup(i,j) = sqrt( dx**2 + dy**2 + dz**2 ) 
          FltGrid_Rjb(i,j) = sqrt( dx**2 + dy**2 )
          fltGrid_Rx(i,j) = 0
          fltGrid_HW(i,j) = 0
          if (fltGrid_Rrup(i,j) .lt. distMin1) distMin1=fltGrid_Rrup(i,j)

        enddo
      enddo

C     Compute the fault lengths for each cell relative to start of fault. 
      do i=1,ny1
        dx = fltGrid_x(i,1) - fltGrid_x(i,2)
        dy = fltGrid_y(i,1) - fltGrid_y(i,2)
        FltGrid_fLen(i,1) = sqrt( dx**2 + dy**2 ) 
        do j=2,nx1
          dx = fltGrid_x(i,j-1) - fltGrid_x(i,j)
          dy = fltGrid_y(i,j-1) - fltGrid_y(i,j)
          FltGrid_fLen(i,j) = FltGrid_fLen(i,j-1) + sqrt( dx**2 + dy**2 ) 
        enddo
      enddo
      
      return
      end

c -------------------------------------------------------------------

      subroutine SetFltBottom ( dip2, nDD, faultThick, flong, flat, fZ, nfp, iFlt,
     1           iCoor )
      include 'pfrisk.h'
      real lat1
      real fLong(MAX_FLT,MAX_DD,MAX_SEG), fLat(MAX_FLT,MAX_DD,MAX_SEG),
     1     fZ(MAX_FLT,MAX_DD,MAX_SEG)
      integer nfp(MAX_FLT), nDD(MAX_FLT) 


C     Set bottom for case in which fault is defined in latitude and longitude. 
      if (iCoor .eq. 1) then     
         sin_theta = sin( abs(dip2) *3.1415926/180.)
         nDD(iFlt) = 2.
         top = fZ(iFlt,1,1)
         bottom = top + faultThick
         lat1 = ( flat(iFlt,1,nfp(iflt)) + flat(iFlt,1,1) )/2. * 3.14159/180.
         dx = ( flong(iFlt,1,nfp(iflt)) - flong(iFlt,1,1) ) * cos(lat1) * 111.12
         dy = (flat(iFlt,1,nfp(iflt)) - flat(iFlt,1,1)) * 111.12
         strike = atan2(dx,dy)
         if ( dip2 .gt. 0 ) then
           az1 = strike + 3.1415926/2.
         else
           az1 = strike - 3.1415926/2.
         endif
         R1 = faultThick / tan(abs(dip2)*3.1415926/180.)
         do ipt=1,nfp(iflt)
           x1 = R1 * sin(az1)
           y1 = R1 * cos(az1)
           flong(iFlt,2,ipt) = flong(iFlt,1,ipt) + x1/(111.12*cos(lat1))
           flat(iFlt,2,ipt)  = flat(iFlt,1,ipt) + y1/111.12
           fZ(iFlt,2,ipt) = bottom
         enddo
C     Case where fault is defined in terms of X,Y coordinates in km. 
      elseif (iCoor .eq. 0) then          
         sin_theta = sin( abs(dip2) *3.1415926/180.)
         nDD(iFlt) = 2.
         top = fZ(iFlt,1,1)
         bottom = top + faultThick
         lat1 = 0.0
         dx = ( flong(iFlt,1,nfp(iflt)) - flong(iFlt,1,1) ) 
         dy = (flat(iFlt,1,nfp(iflt)) - flat(iFlt,1,1)) 
         strike = atan2(dx,dy)
         if ( dip2 .gt. 0 ) then
           az1 = strike + 3.1415926/2.
         else
           az1 = strike - 3.1415926/2.
         endif
         R1 = faultThick / tan(abs(dip2)*3.1415926/180.)
         do ipt=1,nfp(iflt)
           x1 = R1 * sin(az1)
           y1 = R1 * cos(az1)
           flong(iFlt,2,ipt) = flong(iFlt,1,ipt) + x1
           flat(iFlt,2,ipt)  = flat(iFlt,1,ipt) + y1
           fZ(iFlt,2,ipt) = bottom
c        write (*,'( i5,10f10.3)') ipt, flong(iFlt,2,ipt),flat(iFlt,2,ipt) 
c     &, fZ(iFlt,2,ipt)
         enddo

      endif

      return
      end
      
c -----------------------------------------------------------------


        
