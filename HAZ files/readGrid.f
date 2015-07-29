c  --------------------------------------------------------------------

       subroutine RdGrid1 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale )
      include 'pfrisk.h'
         
      character*80 filein
      real*8 sum, sum1
      real minLat, minLong, maxLat, maxLong
      real grid_long(MAX_FLT,MAX_GRID), grid_lat(MAX_FLT,MAX_GRID),
     1     grid_a(MAX_FLT,MAX_GRID)
      real grid_dlong(MAX_FLT), grid_dlat(MAX_FLT)
      integer grid_n(MAX_FLT)
                

      read (10,'( a80)') filein
      open (11,file=filein,status='old')
          
c     Read header
      read (11,*) nHead
      do i=1,nHead
        read (11,'( a1)') dummy
      enddo
      read (11,*) n, grid_dlong(iFlt), grid_dlat(iFlt)
      if ( n .gt. MAX_GRID ) then
        write (*,'( 2x,''Increase MAX_GRID to '',i7)') n
        stop 99
      endif
          
c     read grid data keeping only grid points with non-zero rates
c     and between min and max long and lat
      j = 1
      sum = 0.
      sum1 = 0.
      do i=1,n
        read (11,*) grid_long(iFlt,j), grid_lat(iFlt,j),grid_a(iFlt,j)
        sum = sum + grid_a(iFlt,j)
        if ( grid_a(iFlt,j) .gt. 0. 
     1       .and. grid_long(iFlt,j) .ge. minLong  
     1       .and. grid_long(iFlt,j) .le. maxLong 
     2       .and. grid_lat(iFlt,j) .ge. minLat 
     2       .and. grid_lat(iFlt,j).le. maxLat ) then
          sum1 = sum1 + grid_a(iFlt,j)
          j = j + 1
        endif
      enddo
      close (11)
      grid_n(iFlt) = j - 1
      rate_scale = sum1/sum
      
c     Note: Rate_scale keeps track of the activity rate that is removed because
c     it is at to large a distance (but is still in the activity rate of the input)            

      return
      end

c ----------------------------------------------------------------------
     
      subroutine RdGrid2 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale,
     2           grid_top )
      include 'pfrisk.h'
         
      character*80 filein
      real*8 sum, sum1
      real minLat, minLong, maxLat, maxLong
      real grid_long(MAX_FLT,1), grid_lat(MAX_FLT,1), grid_a(MAX_FLT,1)
      real grid_dlong(1), grid_dlat(1), grid_top(MAX_FLT,1)
      integer grid_n(1)
                      

      read (10,'( a80)') filein
      open (11,file=filein,status='old')
          
c     Read header
      read (11,*) nHead
      do i=1,nHead
        read (11,'( a1)') dummy
      enddo
      read (11,*) n, grid_dlong(iFlt), grid_dlat(iFlt)
      if ( n .gt. MAX_GRID ) then
        write (*,'( 2x,''Increase MAX_GRID to '',i7)') n
        stop 99
      endif
          
c     read grid data keeping only grid points with non-zero rates
c     and between min and max long and lat
      j = 1
      sum = 0.
      sum1 = 0.
      do i=1,n
        read (11,*) grid_long(iFlt,j), grid_lat(iFlt,j),grid_top(iFlt,j), grid_a(iFlt,j)
        sum = sum + grid_a(iFlt,j)
        if ( grid_a(iFlt,j) .gt. 0. 
     1       .and. grid_long(iFlt,j) .ge. minLong  
     1       .and. grid_long(iFlt,j) .le. maxLong 
     2       .and. grid_lat(iFlt,j) .ge. minLat 
     2       .and. grid_lat(iFlt,j).le. maxLat ) then
          sum1 = sum1 + grid_a(iFlt,j)
          j = j + 1
        endif
      enddo
      close (11)
      grid_n(iFlt) = j - 1
      rate_scale = sum1/sum
      
c     Note: Rate_scale keeps track of the activity rate that is removed because
c     it is at to large a distance (but is still in the activity rate of the input)            

      return
      end


     
