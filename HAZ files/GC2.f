
       subroutine GC2 (n3, MAXFLT_AS, rup_x, rup_y, rup_xb, rup_yb, Global_T, Global_U)
      
       implicit none
       integer n3, MAXFLT_AS
       real rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), rup_xb(MAXFLT_AS), rup_yb(MAXFLT_AS), 
     1      Global_T, Global_U

       integer iGC2, tflag, uflag       
       real Seg_length(MAXFLT_AS), Strike_slope(MAXFLT_AS), 
     1      Normal_slope(MAXFLT_AS), a, P90_x(MAXFLT_AS), P90_y(MAXFLT_AS), 
     2      Site_x, Site_y, t_local(MAXFLT_AS), u_local(MAXFLT_AS), 
     3      Seg_weight(MAXFLT_AS), Seg_weight_t(MAXFLT_AS), sum_Weight, 
     4      rec_Weight, Seg_x(MAXFLT_AS), Seg_wxu(MAXFLT_AS), sum_Swt, sum_Swxu

c     calculate the length of each segment of the rupture (a segment of the
c     rupture is from the center of one rupture cell to the center of the 
c     next rupture cell, so the segment length should be approximately equal
c     to the step size along strike)         

        Site_x = 0.0
        Site_y = 0.0

        do iGC2=1, n3-1
          Seg_length(iGC2) = sqrt(((rup_x(iGC2+1)-rup_x(iGC2))**2)+((
     1                    rup_y(iGC2+1)-rup_y(iGC2))**2))
          Strike_slope(iGC2) = (rup_y(iGC2+1)-rup_y(iGC2))/(rup_x(iGC2+1)-
     1                    rup_x(iGC2))
          Normal_slope(iGC2) = (-1)/Strike_slope(iGC2) 
        enddo      

c       calculate x and y coordinate of point where local 
c       u and t form 90 degree angle 

          do iGC2=1, n3-1
            a = rup_x(iGC2+1) - rup_x(iGC2)  
              if (a.EQ.0) then
                P90_x(iGC2) = rup_x(iGC2)
                P90_y(iGC2) = Site_y
              else              
                P90_x(iGC2) = (((-1)*Strike_slope(iGC2)*rup_x(iGC2)) +
     1                     rup_y(iGC2) + (Normal_slope(iGC2)*Site_x) -
     2                     Site_y)/(Normal_slope(iGC2) - 
     3                     Strike_slope(iGC2))     
                P90_y(iGC2) = (Strike_slope(iGC2) * (P90_x(iGC2) - 
     1                     rup_x(iGC2))) + rup_y(iGC2)  
              endif
          enddo

c       calculate local strike-normal coordinate t  

          do iGC2=1, n3-1
              call GC2_tsign (iGC2, MAXFLT_AS, rup_x, rup_y, rup_xb, rup_yb, tflag) 
              t_local(iGC2) = tflag*sqrt(((P90_x(iGC2) - Site_x)**2) +
     1                        ((P90_y(iGC2) - Site_y)**2))               
          enddo

c       calculate local strike-parallel coordinate u  

          do iGC2=1, n3-1
              call GC2_usign (iGC2, MAXFLT_AS, rup_x, rup_y, uflag)             
              u_local(iGC2) = uflag*sqrt(((P90_x(iGC2) - rup_x(iGC2))**2) +
     1                         ((P90_y(iGC2) - rup_y(iGC2))**2))                
          enddo

c        check for special case t=0, on extension and assign 
c          alternative segment weight
c        check for special case t=0, on segment and assign dummy
c          segment weight = 0, which will not be used 
c        if no special case, assign Segment weight as analytical 
c          solution of 1/(r^2) evaluated from 0 to Seg_length

          do iGC2=1, n3-1
            if (t_local(iGC2).EQ.0 .and. u_local(iGC2).LT.0) then
              Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) - 
     1                           (1/u_local(iGC2))
            else if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GT.Seg_length(iGC2)) then
              Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) - 
     1                           (1/u_local(iGC2))
            else if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GE.0 .and. 
     1        u_local(iGC2).LE.Seg_length(iGC2)) then
              Seg_weight(iGC2) = 0.0
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
          do iGC2=1, n3-1
            sum_Weight = sum_Weight + Seg_weight(iGC2)
          enddo
            rec_Weight = 1/sum_Weight

c       calculate where you are on the strike of the fault for each
c       segment and segment weight * Seg_x + u_local

        Seg_x(1)=0.0
        do iGC2=2, n3-1
          Seg_x(iGC2) = Seg_x(iGC2-1) + Seg_length(iGC2-1)
        enddo
        
        do iGC2=1, n3-1
          Seg_wxu(iGC2) = Seg_weight(iGC2)*(Seg_x(iGC2)+u_local(iGC2))
        enddo
        
c       calculate Global Coordinate T
c       check for special case t=0 on segment

          sum_Swt = 0.0
          do iGC2=1, n3-1
            sum_Swt = sum_Swt + Seg_weight_t(iGC2)
            Global_T = rec_Weight*sum_Swt
          enddo  
          do iGC2=1, n3-1
            if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GE.0 .and. 
     1        u_local(iGC2).LE.Seg_length(iGC2)) then 
              Global_T = 0.0
            endif              
          enddo      

c       calculate Global Coordinate U
c       check for special case t=0 on segment

          sum_Swxu = 0.0
          do iGC2=1, n3-1
            sum_Swxu = sum_Swxu + Seg_wxu(iGC2)
            Global_U = rec_Weight*sum_Swxu
          enddo
          do iGC2=1, n3-1  
            if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GE.0 .and. 
     1        u_local(iGC2).LE.Seg_length(iGC2)) then 
              Global_U = Seg_x(iGC2)+u_local(iGC2)
            endif              
          enddo    
  
        return
       end

c ----------------------------------------------------------------------

       subroutine GC2_tsign (iGC2, MAXFLT_AS, rup_x, rup_y, rup_xb, rup_yb, tflag)
       
       implicit none
       
       integer iGC2, MAXFLT_AS, tflag 
       real rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), rup_xb(MAXFLT_AS), rup_yb(MAXFLT_AS)
        
       integer tpositive, tnegative 
       real strikeX, strikeY, strike, xtemp, ytemp, xrupend1, xrupend2, yrupend1, 
     1      yrupend2, adistX1, adistY1, adist1, astrike1, xtest(5), ytest(5), 
     2      xtestfw(5), ytestfw(5)

C      Compute average strike for segment.
       strikeX = rup_x(iGC2+1) - rup_x(iGC2)
       strikeY = rup_y(iGC2+1) - rup_y(iGC2)
       if (strikeX .eq. 0.0) then
          strike = 0.0
       else
          strike = atan2(strikeX,strikeY)
       endif

C      Extend the two end points by 1000 km along the strike for this segment
C      Site is assumed to be location: (0.0, 0.0, 0.0)      
       xtemp = 1000.0*sin(strike)
       ytemp = 1000.0*cos(strike)
       xrupend1 = rup_x(iGC2) - xtemp
       xrupend2 = rup_x(iGC2) + xtemp
       yrupend1 = rup_y(iGC2) - ytemp
       yrupend2 = rup_y(iGC2) + ytemp
       
C      Now determine if station is located on Hanging wall side or footwall side of fault rupture area. 
C      Reset HW or FW Flag for each rupture area.
       tpositive = 0
       tnegative = 0
       
C      Extend downdip point to check for site being on HW side.
C      First compute the distances and angles between upper and lower points on the rupture area. 
       adistX1 = rup_xb(iGC2)- rup_x(iGC2)
       adistY1 = rup_yb(iGC2) - rup_y(iGC2)
       adist1 = sqrt (adistX1*adistX1 + adistY1*adistY1)
       astrike1 = atan2(adistY1,adistX1)     

c      Extend segment end point locations by 100 additional km down dip
c      and end points by 1000 km along strike.
c      Set up testing points: 1-Upper left, 2-Upper right, 3-Lower left, 4-Lower right
       xtest(1) = xrupend1 
       ytest(1) = yrupend1 
       xtest(2) = xrupend2 
       ytest(2) = yrupend2 
       xtest(3) = xrupend1 + 100.0*cos(astrike1)*(adist1)
       ytest(3) = yrupend1 + 100.0*sin(astrike1)*(adist1)
       xtest(4) = xrupend2 + 100.0*cos(astrike1)*(adist1)
       ytest(4) = yrupend2 + 100.0*sin(astrike1)*(adist1)
       xtest(5) = xtest(1)
       ytest(5) = ytest(1)

c      Extend top of rupture end point locations by 100 additional km away from the 
c      down dip direction and end points by 1000 km along strike
c      Set up testing points: 1-Upper left, 2-Upper right, 3-Lower left, 4-Lower right
       xtestfw(1) = xrupend1 
       ytestfw(1) = yrupend1 
       xtestfw(2) = xrupend2 
       ytestfw(2) = yrupend2 
       xtestfw(3) = xtestfw(1) - (xtest(3) - xtest(1) )
       ytestfw(3) = ytestfw(1) - (ytest(3) - ytest(1) )
       xtestfw(4) = xtestfw(2) - (xtest(4) - xtest(2) )
       ytestfw(4) = ytestfw(2) - (ytest(4) - ytest(2) )
       xtestfw(5) = xtestfw(1)
       ytestfw(5) = ytestfw(1)

C      Check to see if site is located over extended rupture area (HW).
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtest, ytest, 0.0, 0.0, tpositive)

C      Check to see if site is located over extended rupture area (FW).
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestfw, ytestfw, 0.0, 0.0, tnegative)       

       if (tpositive.eq.1.) then
         tflag = 1.
       elseif (tnegative.eq.1.) then
         tflag = -1.
       endif
       
        return
       end
       
c ----------------------------------------------------------------------

       subroutine GC2_usign (iGC2, MAXFLT_AS, rup_x, rup_y, uflag)
       
       implicit none
       
       integer iGC2, MAXFLT_AS, uflag 
       real rup_x(MAXFLT_AS), rup_y(MAXFLT_AS)
        
       integer upositive, unegative 
       real strikeX, strikeY, strike, normal, xtemp, ytemp, xtendp, ytendp,
     1      xtendn, ytendn, xtempp, ytempp, xtestp(5), ytestp(5), xtempn, 
     2      ytempn, xtestn(5), ytestn(5)

C      Compute strike and strike-normal for segment
       strikeX = rup_x(iGC2+1) - rup_x(iGC2)
       strikeY = rup_y(iGC2+1) - rup_y(iGC2)
       strike = atan2(strikeX,strikeY)   
       normal = strike - 3.141592653590/2.0  
       
       write (*,*) 'iGC2 ', iGC2
       write (*,*) 'rup_x(iGC2) ', rup_x(iGC2)
       write (*,*) 'rup_y(iGC2) ', rup_y(iGC2)
       write (*,*) 'rup_x(iGC2+1) ', rup_x(iGC2+1)
       write (*,*) 'rup_y(iGC2+1) ', rup_y(iGC2+1)
       write (*,*) 'strikeX ', strikeX
       write (*,*) 'strikeY ', strikeY
       write (*,*) 'strike ', strike
       write (*,*) 'normal ', normal
       pause

C      Extend the first point of the segment by 1000 km in each 
C      of the normal directions (t positive and t negative)
C      Site is assumed to be location: (0.0, 0.0, 0.0)      
       xtemp = 1000.0*cos(normal)
       ytemp = 1000.0*sin(normal)
       xtendp = rup_x(iGC2) + xtemp
       ytendp = rup_y(iGC2) + ytemp       
       xtendn = rup_x(iGC2) - xtemp
       ytendn = rup_y(iGC2) - ytemp
       
       write (*,*) 'xtemp ', xtemp
       write (*,*) 'ytemp ', ytemp
       write (*,*) 'xtendp ', xtendp
       write (*,*) 'ytendp ', ytendp
       write (*,*) 'xtendn ', xtendn
       write (*,*) 'ytendn ', ytendn
       pause
       
C      Now determine if site is located in the direction of rupture 
C      (u+) or in the opposite direction of rupture (u-) relative to 
C      the first point of this segment. Reset uflag for each segment.
       upositive = 0
       unegative = 0     

c      Extend segment end point locations by 1000 km in the direction 
C      of rupture. Set up testing points for u+
       
       xtempp = 1000*cos(strike)
       ytempp = 1000*sin(strike)
       
       write (*,*) 'xtempp ', xtempp
       write (*,*) 'ytempp ', ytempp
       pause
       
       xtestp(1) = xtendp 
       ytestp(1) = ytendp 
       xtestp(2) = xtendn 
       ytestp(2) = ytendn 
       xtestp(3) = xtestp(2) + xtempp
       ytestp(3) = ytestp(2) + ytempp
       xtestp(4) = xtestp(1) + xtempp
       ytestp(4) = ytestp(1) + ytempp
       xtestp(5) = xtestp(1)
       ytestp(5) = ytestp(1)
       
       write (*,*) 'xtestp(1) ', xtestp(1)
       write (*,*) 'ytestp(1) ', ytestp(1)
       write (*,*) 'xtestp(2) ', xtestp(2)
       write (*,*) 'ytestp(2) ', ytestp(2)
       write (*,*) 'xtestp(3) ', xtestp(3)
       write (*,*) 'ytestp(3) ', ytestp(3)
       write (*,*) 'xtestp(4) ', xtestp(4)
       write (*,*) 'ytestp(4) ', ytestp(4)
       write (*,*) 'xtestp(5) ', xtestp(5)
       write (*,*) 'ytestp(5) ', ytestp(5)
       pause
       
c      Extend segment end point locations by 1000 km in the direction 
C      opposite of rupture. Set up testing points for u-

       xtempn = 1000*cos(-strike)
       ytempn = 1000*sin(-strike)
       
       write (*,*) 'xtempn ', xtempn
       write (*,*) 'ytempn ', ytempn
       pause
       
       xtestn(1) = xtendp 
       ytestn(1) = ytendp 
       xtestn(2) = xtendn 
       ytestn(2) = ytendn 
       xtestn(3) = xtestn(2) + xtempn
       ytestn(3) = ytestn(2) + ytempn
       xtestn(4) = xtestn(1) + xtempn
       ytestn(4) = ytestn(1) + ytempn
       xtestn(5) = xtestn(1)
       ytestn(5) = ytestn(1)
       
       write (*,*) 'xtestn(1) ', xtestn(1)
       write (*,*) 'ytestn(1) ', ytestn(1)
       write (*,*) 'xtestn(2) ', xtestn(2)
       write (*,*) 'ytestn(2) ', ytestn(2)
       write (*,*) 'xtestn(3) ', xtestn(3)
       write (*,*) 'ytestn(3) ', ytestn(3)
       write (*,*) 'xtestn(4) ', xtestn(4)
       write (*,*) 'ytestn(4) ', ytestn(4)
       write (*,*) 'xtestn(5) ', xtestn(5)
       write (*,*) 'ytestn(5) ', ytestn(5)
       pause

C      Check to see if site is located in the u+ testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestp, ytestp, 0.0, 0.0, upositive)

C      Check to see if site is located in the u- testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestn, ytestn, 0.0, 0.0, unegative)       

       if (upositive.eq.1.) then
         uflag = 1.
       elseif (unegative.eq.1.) then
         uflag = -1.
       endif
       
       write (*,*) 'iGC2 ', iGC2
       write (*,*) 'uflag ', uflag
       pause

        return
       end       
