       subroutine nLocXYdims (sourceType, rupWidth, aveWidth, rupArea, faultLen, 
     1                       xStep, yStep, Grid_n, faultWidth, nLocX, nLocY, 
     2                       nLocXST1, nLocYST1, nLocXAS, nLocYAS, rupLen)
       
       implicit none
       integer sourceType, Grid_n
       real rupWidth, aveWidth, rupArea, faultLen, xStep, yStep, faultWidth
       
       integer nLocX, nLocY, nLocXST1, nLocYST1, nLocXAS, nLocYAS
       real rupLen
       
c           Set rupture locations and probabilities for faults
            if (sourceType .eq. 1 ) then                     
              if (rupWidth .gt. aveWidth ) then
                rupWidth = aveWidth
              endif
              rupLen = rupArea / rupWidth
              if (rupLen .gt. faultLen) then
                rupLen = faultLen
              endif
              nLocX = (faultLen - RupLen)/xStep + 1       
              nLocY = (aveWidth - rupWidth)/yStep + 1
              nLocXST1 = nLocX
              nLocYST1 = nLocY           
            elseif (sourceType .eq. 5 .or. sourceType .eq. 6) then
              rupLen = rupArea / rupWidth
              if (rupLen .gt. faultLen) rupLen = faultLen
              nLocX = (faultLen - RupLen)/xStep + 1
            elseif (sourceType .eq. 4 ) then
              nLocXAS = Grid_n 
              nLocYAS = faultWidth / yStep + 1
              if (rupWidth .gt. faultWidth) then                                        
                 rupWidth = faultWidth                                                          
              endif                                                                     
              rupLen = rupArea / rupWidth
              
c           Otherwise source is area  
            else                  
              nLocYAS = faultWidth / yStep
              if (nLocYAS.eq.0) then
                nLocYAS = 1
              endif    
              if (rupWidth .gt. faultWidth) then                                        
                rupWidth = faultWidth                                                          
              endif                                                                     
              rupLen = rupArea / rupWidth
            endif

        return
       end

c ----------------------------------------------------------------------
       subroutine nLocXcells (MAXFLT_AS, MAXFLT_DD, nfltgrid, sourceType, nLocYST1, 
     1                        rupWidth, fltgrid_w, fltgrid_a, ruparea, nLocXAS, 
     2                        nLocX, n1AS, n2AS)
     
       implicit none        
       integer MAXFLT_AS, MAXFLT_DD, nfltgrid(2), sourceType, nLocYST1, nLocXAS, 
     1         nLocX, n1AS(MAXFLT_AS), n2AS(MAXFLT_AS) 
       real rupWidth, fltgrid_w(MAXFLT_DD,MAXFLT_AS), fltgrid_a(MAXFLT_DD,MAXFLT_AS), 
     1      ruparea 
        
       integer countnLocX, in2m, in2, iw, in1, iii, iiii, jjj, in1m, 
     1         n1var(MAXFLT_DD,MAXFLT_AS) 
       real colarea, colarew, awid, cwid(MAXFLT_AS), asum, arwid(MAXFLT_AS), rwidth, 
     1      rcolwid(MAXFLT_AS), colarel
     
C     Need to do loop along strike adding up area for all cells with width 
C     less than RupWidth

       countnLocX = 0       
       do 1234 in2m=1,nfltgrid(2) 
          colarea = 0.0
          
C     Check for case in which total area is less than RupArea.
          if (in2m .gt. 1 .and. countnLocX .eq. 0) then   
             goto 1233
          endif
          do 1235 in2=in2m,nfltgrid(2)  
             colarew = 0.0
             if (in2 .eq. in2m) then            
                do iw=1,nfltgrid(1)              
                   colarew = colarew + fltgrid_w(iw,in2)

C      Find the number of downdip cells for this along strike location to 
C      get the correct rupture width.
                   if (colarew .gt. rupWidth) then                   
                      colarew = 0.0
                      goto 2341
                   endif
                enddo
 2341           continue
                if (iw .gt. nfltgrid(1) ) iw = nfltgrid(1)                
                endif
             do 1236 in1=1,iw             
                   colarea = colarea + fltgrid_a(in1,in2)           
                   
 1236        continue
                if (colarea .ge. ruparea) then              
                   countnLocX = countnLocX + 1
                   n1AS(in2m) = nfltgrid(1) - iw + 1
                   n2AS(in2m) = in2                  
                   goto 1234
c                endif
                elseif (in2.eq.nfltgrid(2) .and. colarea.lt.ruparea) then
                  goto 1233
                endif

 1235     continue
 1234  continue

 1233  continue
      if (countnLocX .eq. 0) then
         nLocX = 1
         if (sourceType .eq. 1) then
            n1AS(1) = nLocYST1
         else
            n1AS(1) = 1
         endif
         n2AS(1) = nfltgrid(2)
      else 
         nLocX = countnLocX
      endif   

C     Now compute the unequal wts for given number of ruptures along strike
      if (sourceType .eq. 6) then
C     First compute the average width of the entire fault.
       awid = 0.0
       do iii=1,nfltgrid(2)
          cwid(iii) = 0.0
          do jjj=1,nfltgrid(1)
             cwid(iii) = cwid(iii) + fltgrid_w(jjj,iii) 
          enddo
          awid = awid + cwid(iii)/nfltgrid(2)
       enddo

C     Next compute the average rupture width for a given rupture area. 
       asum = 0.0
       do iii=1,nLocX
          arwid(iii) = 0.0
          do iiii=iii,n2AS(iii)
             arwid(iii) = arwid(iii) + cwid(iiii)/(n2AS(iii)-iii)
          enddo
          arwid(iii) = arwid(iii)/awid
          asum = asum + arwid(iii)
       enddo

C     Now normalize the arwid values.
       do iii=1,nLocX
          arwid(iii) = arwid(iii)/asum
       enddo

C     First sum up all of the column widths for a given along strike locations 
        rwidth = 0.0
         do iii=1,nLocX
            rcolwid(iii) = 0.0
            do jjj=1,nfltgrid(1)
               rwidth = rwidth + fltgrid_w(jjj,iii)
               rcolwid(iii) = rcolwid(iii) + fltgrid_w(jjj,iii)
            enddo
         enddo      
C     Next normalize each column location based on the column width divided by total width. 
         do iii=1,nLocX
               rcolwid(iii) = rcolwid(iii)/rwidth
         enddo      
      endif
 
C     Loop along each column along the strike
       do in2=1,nfltgrid(2)
c     Loop at each down dip cell for a given strike location (i.e., column). 
          do in1m=1,nfltGrid(1)
             n1var(in1m,in2) = 0
             colarea = 0.0
             colarew = 0.0
             colarel = 0.0
c     Loop down the given column and find the number of cells to get the expected rupture width
               do in1=in1m,nfltgrid(1)
                  colarea = colarea + fltgrid_a(in1,in2)
                  colarew = colarew + fltgrid_w(in1,in2)
c             Set the number of cells downdip for a given along strike column that is equal to the rupWidth
                  if (colarew .ge. rupWidth) then
                     if (n1var(in1m,in2) .eq. 0) then
                         n1var(in1m,in2) = in1
                     endif
                  endif
                  if (in1 .eq. nfltgrid(1) .and. n1var(in1m,in2) .eq. 0) then
                     n1var(in1m,in2) = nfltgrid(1)
                  endif
               enddo
          enddo
          colarel = colarel + fltgrid_a(1,in2)/fltgrid_w(1,in2)
       enddo
       
C     Reset nLocX for areal and grid sources
      if (sourceType .eq. 2 .or. sourceType .eq. 3) then
         nLocX = nLocXAS         
      endif 
        
c CH            if (sourcetype(iFlt) .eq. 1) nLocX = nLocXST1

       return
       end

c ----------------------------------------------------------------------
       subroutine nLocYcells (MAXFLT_AS, MAXFLT_DD, MAX_DIST1, MAX_GRID, 
     1                       nfltgrid, fltgrid_w, iLocX, rupWidth, n1AS, 
     2                       sourceType, nLocX, xStep, nLocYAS, distDensity,
     3                       distDensity2, grid_x, x0, grid_y, y0, nLocY, pLocX,
     4                       r_horiz)
     
       implicit none        
       integer MAXFLT_AS, MAXFLT_DD, MAX_DIST1, MAX_GRID, nfltgrid(2), 
     1         iLocX, n1AS(MAXFLT_AS), sourceType, nLocX, nLocYAS, 
     2         grid_x(MAX_GRID), grid_y(MAX_GRID), nLocY    
       real fltgrid_w(MAXFLT_DD,MAXFLT_AS), rupWidth, distDensity(MAX_DIST1), 
     1      xStep, distDensity2(MAX_GRID), x0, y0, pLocX
       
       integer irw
       real swidth, testw, r_horiz
     
C      New code for Varibable fault plane case
C      Compute the number of rupture locations down dip for given location along strike and rupwidth
         swidth = 0.0
           do irw=1,nfltgrid(1)
             swidth = swidth + fltgrid_w(irw,iLocX)               
           enddo
             testw = 0.0
             
C     If rupture width at given along strike location is less than requested rupture width
C        set number of down dip rupture locations equal to 1 (i.e., entire width ruptures).
C        For this case extra columns of rupture will be needed to make the areas needed and 
C        will be adjusted in the CALDIST subroutine. 
              if (swidth .le. rupWidth) then
                 nLocY = 1
                 goto 2345
              endif
              do irw=nfltgrid(1),1,-1
                 testw = testw + fltgrid_w(irw,iLocX)
                 if (testw .gt.rupWidth) then
                    nLocY = irw
                    goto 2345
                 endif
                 nLocY = irw
              enddo
 2345         continue

C     Set nLocY equal to n1AS values for given column along strike location.
           nLocY = n1AS(iLocX)

              if (sourceType .eq. 1 ) then
                pLocX = 1./nLocX                
              elseif (sourceType .eq. 5 ) then
                pLocX = 1./nLocX
              elseif (sourceType .eq. 6 ) then
                pLocX = 1./nLocX
              elseif ( sourceType .eq. 2 .or. sourceType .eq. 3 ) then
                pLocX = distDensity(iLocX)               
                if ( pLocX .ne. 0. ) then
                  r_horiz = xStep * (iLocX-0.5)                
                  nLocY = nLocYAS  
                endif             
              elseif ( sourceType .eq. 4 ) then
                pLocX = distDensity2(iLocX)
                r_horiz = sqrt( (grid_x(iLocX)-x0)**2 + (grid_y(iLocX)-y0)**2 )
              endif
              
c CH             if (sourcetype(iFlt) .eq. 1) nLocY = nLocYST1

       return
       end
