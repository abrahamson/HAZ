
      subroutine ConvertCoordinates2 (MAX_FLT, MAX_GRID, MAX_DD, MAX_SEG, nfp, iFlt, 
     1           iCoor, grid_n, sourceType, nDD, siteX, siteY, fLat, fLong, fZ, 
     2           grid_lat, grid_long, grid_dlat, grid_dlong, nPts1, xFlt, yFlt, 
     3           zFlt, grid_x, grid_y, grid_dx, grid_dy, x0, y0, z0)
     
c     This subroutine uses the haversine formula to calculate distance, bearing, and
c     convert latitude/longitude to x/y coordinates on the basis of a spherical
c     earth (ignoring ellipsoidal effects).

c     Reference: http://www.movable-type.co.uk/scripts/latlong.html
     
      implicit none

      integer MAX_FLT, MAX_GRID, MAX_DD, MAX_SEG, nfp, iFlt, iCoor, grid_n, sourceType, 
     1        nDD, nPts1
      real siteX, siteY, fLat(MAX_FLT,MAX_DD,MAX_SEG), fLong(MAX_FLT,MAX_DD,MAX_SEG), 
     1     fZ(MAX_FLT,MAX_DD,MAX_SEG), grid_lat(MAX_FLT,MAX_GRID), 
     2     grid_long(MAX_FLT,MAX_GRID), grid_dlat(MAX_FLT), grid_dlong(MAX_FLT), 
     3     xFlt(MAX_DD,1), yFlt(MAX_DD,1), zFlt(MAX_DD,1), grid_x(MAX_GRID), 
     4     grid_y(MAX_GRID), grid_dx, grid_dy, x0, y0, z0
     
      integer i, iz, iPt
      real Rearth, refLat, refLong, phi1, phi2, deltaphi, lambda1, lambda2, 
     1     deltalambda, a, c, dist, bearing, xscale, yscale
      
      Rearth = 6371.        
     
c     Load npts in fault and dip into 1-D arrays
      npts1 = nfp
      
      if (iCoor .eq. 1.) then
c     Convert Lat, long to km (using site coordinates as ref)
        refLat = siteY
        refLong = siteX

c       Check for grid source
        if ( sourceType .eq. 3 .or. sourceType .eq. 4 ) then       
          do i=1,grid_n
              phi1 = refLat*3.1415926/180.
              phi2 = grid_lat(iFlt,i)*3.1415926/180.
              deltaphi = phi2 - phi1
              lambda1 = refLong*3.1415926/180.
              lambda2 = grid_long(iFlt,i)*3.1415926/180.
              deltalambda = lambda2 - lambda1
              a = (sin(deltaphi/2.)*sin(deltaphi/2.))+(cos(phi1)*cos(phi2)
     1            *sin(deltalambda/2.)*sin(deltalambda/2.))
              c = 2*atan2(sqrt(a),sqrt(1.-a))
              dist = Rearth*c
              if (dist .eq. 0.) then
                bearing = 0.
              else
                bearing = atan2(sin(lambda2-lambda1)*cos(phi2),cos(phi1)*
     1                    sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda2-lambda1))
              endif
              grid_x(i) = dist*sin(bearing)
              grid_y(i) = dist*cos(bearing)
          enddo
c         Use flat earth scaling for grid_dx and grid_dy only
          xscale = 111.12 * cos(reflat/180*3.1415926)
          yscale = 111.12
          grid_dx = grid_dlong(iFlt)*xscale 
          grid_dy = grid_dlat(iFlt)*yscale
        else

          do iz=1,nDD
            do iPt=1,nfp
              phi1 = refLat*3.1415926/180.
              phi2 = fLat(iFlt,iz,iPt)*3.1415926/180.
              deltaphi = phi2 - phi1
              lambda1 = refLong*3.1415926/180.
              lambda2 = fLong(iFlt,iz,iPt)*3.1415926/180.
              deltalambda = lambda2 - lambda1
              a = (sin(deltaphi/2.)*sin(deltaphi/2.))+(cos(phi1)*cos(phi2)
     1            *sin(deltalambda/2.)*sin(deltalambda/2.))
              c = 2*atan2(sqrt(a),sqrt(1.-a))
              dist = Rearth*c
              if (dist .eq. 0.) then
                bearing = 0.
              else
                bearing = atan2(sin(lambda2-lambda1)*cos(phi2),cos(phi1)*
     1                    sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda2-lambda1))
              endif
              xFlt(iz,iPt) = dist*sin(bearing)
              yFlt(iz,iPt) = dist*cos(bearing)
              zFlt(iz,iPt) = fZ(iFlt,iz,iPt)            
            enddo
          enddo
        endif
        
      else
c     Units are in km already 
      x0 = 0.
      y0 = 0.
      z0 = 0.
      endif

      return 
      end


c ---------------------------------------------------------------------

      subroutine calcFltGrid ( xFlt, yFlt, zFlt, npts, nDD, fltGrid_X, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGRid_Rrup, FltGRid_Rjb, fltGrid_Rseis, fltGrid_Rx, FltGrid_HW,
     3               faultArea, faultLen, faultW, distMin1, step, FltGrid_fLen )  
     
      include 'pfrisk.h'
      include 'declare1.h'
      integer nDD, n2(MAXFLT_DD), nfltGrid(2), n1(MAXFLT_AS)
      real fltGrid_X(MAXFLT_DD,MAXFLT_AS), fltGrid_y(MAXFLT_DD,MAXFLT_AS), 
     1     fltGrid_z(MAXFLT_DD,MAXFLT_AS),fltGrid_fLen(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_w(MAXFLT_DD,MAXFLT_AS), fltGrid_a(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_RJB(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_Rseis(MAXFLT_DD,MAXFLT_AS), fltGrid_Rx(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_HW(MAXFLT_DD,MAXFLT_AS), faultArea
      real*8 dfaultArea

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
          sum1 = sum1 + real(w/npts)
          nn = nn + 1
        enddo
        n2(i) = nint( width1 / step )
        if (n2(i).eq.0) then
          n2(i) = 1
        endif
      enddo
      
      faultW = sum1 
      nx1 = 0
      ny1 = 0
      ii = 0
      i0 = 0
      j0 = 0
      dfaultArea = 0
      faultLen = 0
      do i=1,nDD-1
        ny1 = ny1 + n2(i)

        do j=2,npts
          if ( i .eq. 1 ) then
            xLtop =  xFlt(i,j) - xFlt(i,j-1)
            yLtop =  yFlt(i,j) - yFlt(i,j-1)
            n1(j) = nint(sqrt(xLtop**2 + yLtop**2)/step)
            if (n1(j) .eq. 1) then
              n1(j) = 1
            endif
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
              dfaultArea = dfaultArea + dble(fltGrid_a(i2,j2))  

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
      
      faultArea = real(dfaultArea)
      
      return
      end

c -------------------------------------------------------------------

      subroutine SetFltBottom (MAX_FLT, MAX_DD, MAX_SEG, iCoor, iFlt, 
     1                         nfp, dip2, faultThick, fZ, flat, flong, nDD)
     
      implicit none
     
c      declarations passed in
       integer MAX_FLT, MAX_DD, MAX_SEG, iCoor, iFlt, nfp(MAX_FLT)
       real dip2, faultThick        

c      declarations passed out
       integer nDD(MAX_FLT)
       
c      declarations passed in and out (changed) 
       real fZ(MAX_FLT,MAX_DD,MAX_SEG), fLat(MAX_FLT,MAX_DD,MAX_SEG),
     1      fLong(MAX_FLT,MAX_DD,MAX_SEG)

c      declarations only used within subroutine   
       integer ipt
       real sin_theta, top, bottom, lat1, dx, dy, strike, az1, R1, x1, y1         

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
         enddo

      endif

      return
      end
      
c -----------------------------------------------------------------


        
