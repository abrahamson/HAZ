      subroutine CalcDepthProb (iDepthModel, depthParam, iFlt, pLocY,
     1              sourceType, nLocY, yStep, top, width, rupWidth, dip )

      implicit none
      include 'pfrisk.h'

c     declarations passed in 
      integer iDepthModel, iFlt, sourceType, nLocY 
      real depthParam(MAX_FLT,3), yStep, top, width, rupWidth, dip

c     declarations passed out
      real pLocY(MAXFLT_AS)
      
c     declarations only used within subroutine 
      integer iLocY     
      real meanDepth, z1, z2, z3, depth1, sigma, depth0,
     1     depth2, z_half, p1, p2, zstep, d, sum
      
c     Compute depth probabilities for areal sources (point sources)
      if (sourceType .eq. 2 .or. sourceType .eq. 3 .or. sourceType .eq. 4 ) then
        if (iDepthModel .eq. 0 ) then
c         Uniform distribution
          do iLocY=1,nLocY
            pLocY(iLocY) = 1./float(nLocY)
          enddo

        elseif (iDepthModel .eq. 2 ) then
c         Triangle model
          z1 = depthParam(iFlt,1)
          z2 = depthParam(iFlt,2)
          z3 = depthParam(iFlt,3)
          do iLocY=1,nLocY
            depth1 = top + yStep*(iLocy-1) + yStep/2.              
            if ( depth1 .lt. z1 ) then
              pLocY(iLocY) = 0.
            elseif ( depth1 .lt. z2 ) then
              pLocY(iLocY) = (depth1-z1)/(z2-z1)
            elseif ( depth1 .lt. z3 ) then
              pLocY(iLocY) = 1. - (depth1-z2)/(z3-z2)
            elseif ( depth1 .ge. z3 ) then
              pLocY(iLocY) = 0.
            endif
          enddo

        elseif (iDepthModel .eq. 1 ) then
c         normal distribution
          meanDepth = depthParam(iFlt,1)                          
          sigma = depthParam(iFlt,2)  

          do iLocY=1,nLocY
            depth0 = top + yStep*(iLocy-1) + yStep/2.   
            depth1 = depth0 - yStep/2.         
            depth2 = depth0 + yStep/2. 

c           COMPUTE PROBABILITY OF depth BETWEEN depth1 and depth2                
            z1 = (depth1-meanDepth)/sigma    
            z2 = (depth2-meanDepth)/sigma    
            call NDTR(z1,p1,d)        
            call NDTR(z2,p2,d)        
            pLocY(iLocY) = (p1 - p2)               
          enddo

        else
          write (*,'( 2x,''invalid depth model'')')
          stop 99
        endif
        
c     Set depth probability for SourceType 7      
      elseif (sourceType .eq. 7) then
        do iLocY=1,nLocY
          pLocY(iLocY) = 1.0 
        enddo
      
c     Otherwise use fault sources (i.e., Sourcetype = 1, 5, or 6)
      else

c      Uniform distribution for faults
        if (iDepthModel .eq. 0 ) then
c         Uniform distribution 
          do iLocY=1,nLocY
            pLocY(iLocY) = 1./float(nLocY)
          enddo

        elseif (iDepthModel .eq. 1) then
          write (*,'( 2x,''normal depth model not working for faults'')')
          stop 99

        elseif (iDepthModel .eq. 2 ) then
c         Triangle distribution (using the hypocenter at the center of the rupture)
          z_half = rupWidth * sin(dip*3.1415926/180.) / 2.

          do iLocY=1,nLocY
            z1 = depthParam(iFlt,1) 
            z2 = depthParam(iFlt,2) 
            z3 = depthParam(iFlt,3) 
            zStep = yStep * sin(dip*3.1415926/180.)

            depth1 = top + zStep*(iLocy-1) + z_half              
            if ( depth1 .lt. z1 ) then
              pLocY(iLocY) = 0.
            elseif ( depth1 .lt. z2 ) then
              pLocY(iLocY) = (depth1-z1)/(z2-z1)
            elseif ( depth1 .lt. z3 ) then
              pLocY(iLocY) = 1. - (depth1-z2)/(z3-z2)
            elseif ( depth1 .ge. z3 ) then
              pLocY(iLocY) = 0.
            endif

          enddo
        endif
       endif

c      normalize prob
       sum = 0.
       do iLocY=1,nLocY
         sum = sum + pLocY(iLocY)
       enddo
       do iLocY=1,nLocY
         pLocY(iLocY) = pLocY(iLocY) / sum
       enddo


      return
      end
        
