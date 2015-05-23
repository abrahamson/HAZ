
      subroutine GC2 (n2, MAXFLT_AS, rup_E, rup_N, Global_T, Global_U)
      
      implicit none
      real rup_E(MAXFLT_AS), rup_N(MAXFLT_AS), Global_T, Global_U
      integer n2, MAXFLT_AS
      real Seg_length(MAXFLT_AS), Strike_slope(MAXFLT_AS), 
     1     Normal_slope(MAXFLT_AS), a, P90_E(MAXFLT_AS), P90_N(MAXFLT_AS), 
     2     Site_E, Site_N, t_local(MAXFLT_AS), u_local(MAXFLT_AS), 
     3     Seg_weight(MAXFLT_AS), Seg_weight_t(MAXFLT_AS), sum_Weight, 
     4     rec_Weight, Seg_x(MAXFLT_AS), Seg_wxu(MAXFLT_AS), sum_Swt, sum_Swxu
      integer iGC2 

c     calculate the length of each segment of the rupture (a segment of the
c     rupture is from the center of one rupture cell to the center of the 
c     next rupture cell, so the segment length should be approximately equal
c     to the step size along strike)         

        Site_E = 0.0
        Site_N = 0.0

        do iGC2=1, n2-1
          Seg_length(iGC2) = sqrt(((rup_E(iGC2+1)-rup_E(iGC2))**2)+((
     1                    rup_N(iGC2+1)-rup_N(iGC2))**2))
          Strike_slope(iGC2) = (rup_N(iGC2+1)-rup_N(iGC2))/(rup_E(iGC2+1)-
     1                    rup_E(iGC2))
          Normal_slope(iGC2) = (-1)/Strike_slope(iGC2) 
        enddo      

c       calculate Easting and Northing coordinate of point where local 
c       u and t form 90 degree angle 

          do iGC2=1, n2-1
            a = rup_E(iGC2+1) - rup_E(iGC2)  
              if (a.EQ.0) then
                P90_E(iGC2) = rup_E(iGC2)
                P90_N(iGC2) = Site_N
              else              
                P90_E(iGC2) = (((-1)*Strike_slope(iGC2)*rup_E(iGC2)) +
     1                     rup_N(iGC2) + (Normal_slope(iGC2)*Site_E) -
     2                     Site_N)/(Normal_slope(iGC2) - 
     3                     Strike_slope(iGC2))     
                P90_N(iGC2) = (Strike_slope(iGC2) * (P90_E(iGC2) - 
     1                     rup_E(iGC2))) + rup_N(iGC2)  
              endif
          enddo

c       calculate local strike-normal coordinate t  

          do iGC2=1, n2-1
              if (Site_E.GT.P90_E(iGC2)) then
                t_local(iGC2) = -sqrt(((P90_E(iGC2) - Site_E)**2) +
     1                         ((P90_N(iGC2) - Site_N)**2)) 
              else              
                t_local(iGC2) = sqrt(((P90_E(iGC2) - Site_E)**2) +
     1                         ((P90_N(iGC2) - Site_N)**2)) 
              endif
          enddo

c       calculate local strike-parallel coordinate u  

          do iGC2=1, n2-1
              if (P90_N(iGC2).GT.rup_N(iGC2)) then
                u_local(iGC2) = -sqrt(((P90_E(iGC2) - rup_E(iGC2))**2) +
     1                         ((P90_N(iGC2) - rup_N(iGC2))**2)) 
              else              
                u_local(iGC2) = sqrt(((P90_E(iGC2) - rup_E(iGC2))**2) +
     1                         ((P90_N(iGC2) - rup_N(iGC2))**2)) 
              endif
          enddo

c        check for special case t=0, on extension and assign 
c          alternative segment weight
c        check for special case t=0, on segment and assign dummy
c          segment weight = 0, which will not be used 
c        if no special case, assign Segment weight as analytical 
c          solution of 1/(r^2) evaluated from 0 to Seg_length

          do iGC2=1, n2-1
            if (t_local(iGC2).EQ.0) then
              if (u_local(iGC2).LT.0) then
                Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) - 
     1                            (1/u_local(iGC2))
              else if (u_local(iGC2).GT.Seg_length(iGC2)) then
                Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) - 
     1                            (1/u_local(iGC2))
              else
                Seg_weight(iGC2) = 0.0
              endif
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
          do iGC2=1, n2-1
            sum_Weight = sum_Weight + Seg_weight(iGC2)
          enddo
            rec_Weight = 1/sum_Weight

c       calculate where you are on the strike of the fault for each
c       segment and segment weight * Seg_x + u_local

        Seg_x(1)=0.0
        do iGC2=2, n2-1
          Seg_x(iGC2) = Seg_x(iGC2-1) + Seg_length(iGC2-1)
        enddo
        
        do iGC2=1, n2-1
          Seg_wxu(iGC2) = Seg_weight(iGC2)*(Seg_x(iGC2)+u_local(iGC2))
        enddo
        
c       calculate Global Coordinate T
c       check for special case t=0 on segment

          sum_Swt = 0.0
          do iGC2=1, n2-1
            sum_Swt = sum_Swt + Seg_weight_t(iGC2)
            Global_T = rec_Weight*sum_Swt
          enddo  
          do iGC2=1, n2-1
            if (t_local(iGC2).EQ.0) then
              if (u_local(iGC2).GE.0) then
                if (u_local(iGC2).LE.Seg_length(iGC2)) then 
                Global_T = 0.0
                endif
              endif
            endif              
          enddo      

c       calculate Global Coordinate U
c       check for special case t=0 on segment

          sum_Swxu = 0.0
          do iGC2=1, n2-1
            sum_Swxu = sum_Swxu + Seg_wxu(iGC2)
            Global_U = rec_Weight*sum_Swxu
          enddo
          do iGC2=1, n2-1  
            if (t_local(iGC2).EQ.0) then
              if (u_local(iGC2).GE.0) then
                if (u_local(iGC2).LE.Seg_length(iGC2)) then 
                Global_U = Seg_x(iGC2)+u_local(iGC2)
                endif
              endif
            endif              
          enddo      
  
      return
      end
