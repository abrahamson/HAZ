     
       subroutine GC2 (MAXFLT_AS, MAXFLT_DD, iLocX, iLocY, n2, n1, fltgrid_x, 
     1                 fltgrid_y, fltgrid_z, Rx, Ry, Ry0, HWFlag, dipavg)
      
       implicit none
       
c      declarations passed in        
       integer MAXFLT_AS, MAXFLT_DD, iLocX, iLocY, n2, n1 
       real fltgrid_x(MAXFLT_DD,MAXFLT_AS), fltgrid_y(MAXFLT_DD,MAXFLT_AS),
     1      fltgrid_z(MAXFLT_DD,MAXFLT_AS)  

c      declarations passed out
       integer HWFlag
       real Rx, Ry, Ry0, dipavg

c      declarations only used within subroutine
       integer inorm, irup, n3, iGC2, tflag, uflag       
       real rup_xt, rup_yt, rup_xbt, rup_ybt, rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), 
     1      rup_xb(MAXFLT_AS), rup_yb(MAXFLT_AS), Seg_length(MAXFLT_AS), 
     2      Strike_slope(MAXFLT_AS), Normal_slope(MAXFLT_AS), GC2_ruplength, a, b,
     3      P90_x(MAXFLT_AS), P90_y(MAXFLT_AS), Site_x, Site_y, t_local(MAXFLT_AS), 
     4      u_local(MAXFLT_AS), Seg_weight(MAXFLT_AS), Seg_weight_t(MAXFLT_AS), 
     5      sum_Weight, rec_Weight, Seg_x(MAXFLT_AS), Seg_wxu(MAXFLT_AS), sum_Swt, 
     6      sum_Swxu, Global_T, Global_U, adistX1, adistY1
         
c      save rupture grid cell locations in new arrays
       inorm = 0
       do irup=iLocX, n2
         rup_xt = fltgrid_x(iLocY,irup)
         rup_yt = fltgrid_y(iLocY,irup)
         rup_xbt = fltgrid_x(n1,irup)
         rup_ybt = fltgrid_y(n1,irup)
         inorm = inorm + 1
         rup_x(inorm) = rup_xt
         rup_y(inorm) = rup_yt
         rup_xb(inorm) = rup_xbt
         rup_yb(inorm) = rup_ybt
       enddo
       n3 = inorm   

c       calculate the length of each segment of the rupture (a segment of the
c       rupture is from the center of one rupture cell to the center of the 
c       next rupture cell, so the segment length should be approximately equal
c       to the step size along strike)         

        Site_x = 0.0
        Site_y = 0.0

        do iGC2=1, n3-1
          Seg_length(iGC2) = sqrt(((rup_x(iGC2+1)-rup_x(iGC2))**2)+((
     1                    rup_y(iGC2+1)-rup_y(iGC2))**2))
          Strike_slope(iGC2) = (rup_y(iGC2+1)-rup_y(iGC2))/(rup_x(iGC2+1)-
     1                    rup_x(iGC2))
          Normal_slope(iGC2) = (-1)/Strike_slope(iGC2)          
        enddo 

c       calculate total segment length (rupture length) for Ry0 adjustment
c       later

        GC2_ruplength = 0.0
        do iGC2=1, n3-1
          GC2_ruplength = GC2_ruplength + Seg_length(iGC2)
        enddo  

c       calculate x and y coordinate of point where local 
c       u and t form 90 degree angle 

          do iGC2=1, n3-1
            a = rup_x(iGC2+1) - rup_x(iGC2)  
            b = rup_y(iGC2+1) - rup_y(iGC2)
              if (a.EQ.0) then
                P90_x(iGC2) = rup_x(iGC2)
                P90_y(iGC2) = Site_y
              else if (b.EQ.0) then
                P90_x(iGC2) = Site_x
                P90_y(iGC2) = rup_y(iGC2)
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
              call GC2_tsign (iGC2, MAXFLT_AS, MAXFLT_DD, n1, iLocX, iLocY, 
     1                        fltgrid_z, rup_x, rup_y, rup_xb, rup_yb, tflag) 
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
          
c       calculate Rx and assign HWFlag from Global Coordinate T  
         HWFlag = 0
         if ((fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)).eq. 0.0) then
           Rx = (-1)*(abs(Global_T))
           HWFlag = 0
         else
           Rx = Global_T  
           if (Global_T .LE. 0.0) then
             HWFlag = 0
           elseif (Global_T .GT. 0.0) then
             HWFlag = 1
           endif    
         endif             

c       calculate Ry from Global Coordinate U
         if (Global_U .LT. 0) then
           Ry = abs(Global_U) + (0.5*GC2_ruplength)
         elseif (Global_U .GE. 0 .and. Global_U .LE. GC2_ruplength) then
           Ry = abs((0.5*GC2_ruplength)-Global_U)
         elseif (Global_U .GT. GC2_ruplength) then
           Ry = Global_U - (0.5*GC2_ruplength)
         endif  

c       calculate Ry0 from Global Coordinate U
         if (Global_U .LT. 0) then
           Ry0 = abs(Global_U)
         elseif (Global_U .GE. 0 .and. Global_U .LE. GC2_ruplength) then
           Ry0 = 0.0
         elseif (Global_U .GT. GC2_ruplength) then
           Ry0 = Global_U - GC2_ruplength
         endif  

c       calculate dipavg
         if ((fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)).eq. 0.0) then
           dipavg = 3.141592653590/2.0
         else
           adistX1 = fltgrid_x(n1,iLocX) - fltgrid_x(iLocY,iLocX)
           adistY1 = fltgrid_y(n1,iLocX) - fltgrid_y(iLocY,iLocX)
           dipavg = atan2((fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)),sqrt(adistX1*adistX1+adistY1*adistY1))
         endif
  
        return
       end
       
c ----------------------------------------------------------------------

       subroutine GC2_tsign (iGC2, MAXFLT_AS, MAXFLT_DD, n1, iLocX, iLocY, fltgrid_z, 
     1                       rup_x, rup_y, rup_xb, rup_yb, tflag)
       
       implicit none
       
c      declarations passed in        
       integer iGC2, MAXFLT_AS, MAXFLT_DD, n1, iLocX, iLocY  
       real fltgrid_z(MAXFLT_DD,MAXFLT_AS), rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), rup_xb(MAXFLT_AS), 
     1      rup_yb(MAXFLT_AS)

c      declarations passed out
       integer tflag
 
c      declarations only used within subroutine       
       integer tpositive, tnegative 
       real strikeX, strikeY, dipX, dipY, ddip, strike, xtemp, ytemp, xtendp, ytendp,
     1      xtendn, ytendn, xtempp, ytempp, xtestp(5), ytestp(5), xtempn, 
     2      ytempn, xtestn(5), ytestn(5)

c      Compute direction of strike and direction of dip (ddip)
c      For pure strike slip fault, use right hand rule for dummy dip direction
       strikeX = rup_x(iGC2+1) - rup_x(iGC2)
       strikeY = rup_y(iGC2+1) - rup_y(iGC2)
       strike = atan2(strikeY,strikeX) 
       if ((fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)).eq. 0.0) then 
         ddip = strike - 3.141592653590/2.0                
       else
         dipX = rup_xb(iGC2) - rup_x(iGC2)
         dipY = rup_yb(iGC2) - rup_y(iGC2)                  
         ddip = atan2(dipY,dipX)    
       endif         

c      Extend the first point of the segment by 1000 km in each 
c      of the strike directions
c      Site is assumed to be location: (0.0, 0.0, 0.0)      
       xtemp = 1000.0*cos(strike)
       ytemp = 1000.0*sin(strike)
       xtendp = rup_x(iGC2) + xtemp
       ytendp = rup_y(iGC2) + ytemp       
       xtendn = rup_x(iGC2) - xtemp
       ytendn = rup_y(iGC2) - ytemp
       
c      Now determine if site is located on the hanging wall (t+) or on
c      the footwall (t-) relative to the segment. 
c      Reset tflag for each segment.
       tpositive = 0
       tnegative = 0     

c      Extend the first point of the segment by 1000 km in the direction 
c      of dip. Set up testing points for t+
       
       xtempp = 1000*cos(ddip)
       ytempp = 1000*sin(ddip)
       
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
       
c      Extend the first point of the segment by 1000 km in the direction 
c      opposite of dip. Set up testing points for t-

       xtempn = 1000*cos(ddip - 3.141592653590)
       ytempn = 1000*sin(ddip - 3.141592653590)
       
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

c      Check to see if site is located in the t+ testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestp, ytestp, 0.0, 0.0, tpositive)

c      Check to see if site is located in the t- testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestn, ytestn, 0.0, 0.0, tnegative)       

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

c      Compute strike and strike-normal for segment
       strikeX = rup_x(iGC2+1) - rup_x(iGC2)
       strikeY = rup_y(iGC2+1) - rup_y(iGC2)
       strike = atan2(strikeY,strikeX)   
       normal = strike - 3.141592653590/2.0 

c      Extend the first point of the segment by 1000 km in each 
c      of the strike normal directions
c      Site is assumed to be location: (0.0, 0.0, 0.0)      
       xtemp = 1000.0*cos(normal)
       ytemp = 1000.0*sin(normal)
       xtendp = rup_x(iGC2) + xtemp
       ytendp = rup_y(iGC2) + ytemp       
       xtendn = rup_x(iGC2) - xtemp
       ytendn = rup_y(iGC2) - ytemp
       
c      Now determine if site is located in the direction of rupture 
c      (u+) or in the opposite direction of rupture (u-) relative to 
c      the first point of this segment. Reset uflag for each segment.
       upositive = 0
       unegative = 0     

c      Extend segment end point locations by 1000 km in the direction 
c      of rupture. Set up testing points for u+
       
       xtempp = 1000*cos(strike)
       ytempp = 1000*sin(strike)
       
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
       
c      Extend segment end point locations by 1000 km in the direction 
c      opposite of rupture. Set up testing points for u-

       xtempn = 1000*cos(strike - 3.141592653590)
       ytempn = 1000*sin(strike - 3.141592653590)
       
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

c      Check to see if site is located in the u+ testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestp, ytestp, 0.0, 0.0, upositive)

c      Check to see if site is located in the u- testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestn, ytestn, 0.0, 0.0, unegative)       

       if (upositive.eq.1.) then
         uflag = 1.
       elseif (unegative.eq.1.) then
         uflag = -1.
       endif

        return
       end       
