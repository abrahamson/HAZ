c  --------------------------------------------------------------------

       subroutine RdGrid1 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale )
      include 'pfrisk.h'
         
      character*80 filein, dummy
      real*8 sum, sum1
      real minLat, minLong, maxLat, maxLong
      real grid_long(MAX_FLT,MAX_GRID), grid_lat(MAX_FLT,MAX_GRID),
     1     grid_a(MAX_FLT,MAX_GRID)
      real grid_dlong(MAX_FLT), grid_dlat(MAX_FLT)
      integer grid_n(MAX_FLT)
                

      read (10,'( a80)') filein
      open (11,file=filein,status='old',err=2100)
          
c     Read header
      read (11,*,err=2001) nHead
      write (*,'( i5)') nHead
      
      do i=1,nHead
        read (11,'( a80)',err=2002) dummy
        write (*,'( a80)') dummy
      enddo
      read (11,*,err=2003) n, grid_dlong(iFlt), grid_dlat(iFlt)
c      write (*,'( i5, 2f10.4)') n, grid_dlong(iFlt), grid_dlat(iFlt)
      if ( n .gt. MAX_GRID ) then
        write (*,'( 2x,''Increase MAX_GRID to '',i7)') n
        stop 99
      endif
      write (*,'( i5)') n
          
c     read grid data keeping only grid points with non-zero rates
c     and between min and max long and lat
      j = 1
      sum = 0.
      sum1 = 0.
      do i=1,n
        read (11,*,err=2004) grid_long(iFlt,j), grid_lat(iFlt,j),grid_a(iFlt,j)
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
      
 2100 write (*,'( 2x,''bad gridded seismicity file'')')
      write (*,'( a80)') filein
      stop 99
 2001 write (*,'( 2x,''gridded seismicity file error: nHeader'')')
      stop 99
 2002 write (*,'( 2x,''gridded seismicity file error: header'')')
      stop 99
 2003 write (*,'( 2x,''gridded seismicity file error: n, grid_dLong, grid_dlat'')')
      stop 99
 2004 write (*,'( 2x,''gridded seismicity file error: long, lat, a-value'')')
      write (*,'( 2x,''entry:'',i5)') i
      stop 99
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


     
