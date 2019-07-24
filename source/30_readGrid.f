c  --------------------------------------------------------------------

       subroutine S30_RdGrid1 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale )

      implicit none
      include 'pfrisk.h'
         
      integer grid_n(MAX_FLT), iFlt, nHead, i, n, j
      real grid_long(MAX_FLT,MAX_GRID), grid_lat(MAX_FLT,MAX_GRID),
     1     grid_a(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT), 
     2     grid_dlat(MAX_FLT), minLat, minLong, maxLat, maxLong, rate_scale     
      real*8 sum, sum1  
      character*80 filein, dummy  
                
      read (10,'( a80)') filein
      open (11,file=filein,status='old',err=2100)
          
c     Read header
      read (11,*,err=2001) nHead
      
      do i=1,nHead
        read (11,'( a80)',err=2002) dummy
c        write (*,'( a80)') dummy
      enddo
      read (11,*,err=2003) n, grid_dlong(iFlt), grid_dlat(iFlt)
      if ( n .gt. MAX_GRID ) then
        write (*,'( 2x,''Increase MAX_GRID to '',i7)') n
        stop 99
      endif
c      write (*,'( i5)') n
          
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
c     it is at too large a distance (but is still in the activity rate of the input)            

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
     
      subroutine S30_RdGrid2 ( iFlt, grid_a,grid_dlong,grid_dlat,grid_n,
     1           grid_long, grid_lat, minLat, minLong, maxLat, maxLong, rate_scale,
     2           grid_top )

      implicit none
      include 'pfrisk.h'

      integer grid_n(1), iFlt, nHead, i, n, j
      real grid_long(MAX_FLT,1), grid_lat(MAX_FLT,1), grid_a(MAX_FLT,1),
     1     grid_dlong(1), grid_dlat(1), grid_top(MAX_FLT,1), minLat, 
     2     minLong, maxLat, maxLong, dummy, rate_scale
      real*8 sum, sum1    
      character*80 filein                     

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
c     it is at too large a distance (but is still in the activity rate of the input)            

      return
      end

c  --------------------------------------------------------------------

      subroutine S30_RdSource7 ( iFlt, mag, rate, dist, dip, mech, ncount )
 
      implicit none
      include 'pfrisk.h'
      
      integer ncount(MAX_FLT), srflag, rupid, nearID, i, iFlt         
      real rate(MAX_FLT,MAX_S7), mag(MAX_FLT,MAX_S7), dist(MAX_FLT,MAX_S7), 
     1     lat, long, strike, Dip(MAX_FLT,MAX_S7), rake, mech(MAX_FLT,MAX_S7)
      character*80 filein
                
      read (10,'( a80)') filein
      open (11,file=filein,status='old')
          
      ncount(iflt) = 0
      do i=1,MAX_S7
         read (11,*,END=777) srflag,rupID,mag(iFlt,i),rate(iFlt,i),nearID,dist(iFlt,i),
     1             lat,long,strike,Dip(iFlt,i),rake
     
         if (rake .ge. 30.0 .and. rake .le. 150.0) then    
           mech(iFlt,i) = 1.0
         elseif (rake .ge. -120.0 .and. rake .le. -60.0) then
           mech(iFlt,i) = -1.0
         else
           mech(iFLt,i) = 0.0
         endif
      
         ncount(iFlt) = ncount(iFlt) + 1
      enddo

      write (*,*) 'More than 70,000 source in data file.'
      write (*,*) 'Need to separate into two files.'
      stop 99
      
 777  continue

C     Classify the Rake angle with Mechanims choices of SS(0), RV(1), or NM(-1)
C     The following Rake Angles are used to classify Fault Mechanims:
C     Reverse: 30<=Rake<=150  (includes RV/OB as RV)
C     Normal:  -120<=Rake<=-60
C     Strike-Slip: -180<=Rake<-120
C                  -60<Rake<=30
C                  150<Rake<=180
C        (includes NM/OB as SS)

      return
      end

     
