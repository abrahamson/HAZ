       subroutine DPP (xPrup, yPrup, zPrup, xPd, yPd, zPd)

       implicit none
     
c      declarations passed in        
       real xPrup(5), yPrup(5), zPrup(5) 

c      declarations passed out
       real xPd, yPd, zPd
     
c      declarations only used within subroutine
       integer Ppflag
       real xPs, yPs, zPs, xu, yu, zu, xv, yv, zv, xn, yn, zn, t, xPp, yPp, zPp     

c      site is assumed to be location 0,0,0
       xPs = 0.0
       yPs = 0.0
       zPs = 0.0

c      three non-collinear points are needed to describe the rupture plane
c      use the first three points
       
c      these three points can be described by two vectors, u and v
       xu = xPrup(2) - xPrup(1)
       yu = yPrup(2) - yPrup(1)
       zu = zPrup(2) - zPrup(1)
       xv = xPrup(3) - xPrup(1)
       yv = yPrup(3) - yPrup(1)
       zv = zPrup(3) - zPrup(1)

c      the cross product of vectors u and v gives us a vector normal to the plane
       xn = (yu*zv)-(zu*yv)
       yn = (-1)*((xu*zv)-(zu*xv))
       zn = (xu*yv)-(yu*xv)
       
c      calculate t
       t = ((xn*xPrup(1))-(xn*xPs)+(yn*yPrup(1))-(yn*yPs)+
     1     (zn*zPrup(1))-(zn*zPs))/(xn**2.+yn**2.+zn**2.)
       
c      calculate the coordinates for Point p
       xPp = xPs+(t*xn)
       yPp = yPs+(t*yn)
       zPp = zPs+(t*zn)        

c      check for special case Point p inside Gamma       
       call Inside_OutSide ( 4, xPrup, yPrup, xPp, yPp, Ppflag)
 
       if (Ppflag.eq.1.) then
         xPd = xPp
         yPd = yPp
         zPd = zPp
       else
         xPd = 1.
         yPd = 1.
         zPd = 1.
       endif

        return
       end
       
c ----------------------------------------------------------------------

       subroutine Gamma_4points (MAXFLT_AS, MAXFLT_DD, n3, iLocX, iLocY,
     1                           dip, rup_x, rup_y, Seg_length, 
     2                           GC2_ruplength, Global_T, Global_U,
     3                           fltgrid_z, xPrup, yPrup, zPrup)  

       implicit none

c      declarations passed in 
       integer MAXFLT_AS, MAXFLT_DD, n1, n2, n3, iLocX, iLocY
       real dip, rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), Seg_length(MAXFLT_AS), 
     1      GC2_ruplength, Global_T, Global_U, fltgrid_z(MAXFLT_DD,MAXFLT_AS)
       
c      declarations passed out   
       real xPrup(5), yPrup(5), zPrup(5)    

c      declarations only used within subroutine
       integer iGC2
       real dipr, xPs, yPs, zPs, sum1, segX, segY, theta(MAXFLT_AS), Ave_theta,
     1      J, t, dx1, dy1, dx2, dy2, z_rup, cdip, dxy3, dx3, dy3

       dipr = dip*3.141592653590/180.

       xPs = 0.0
       yPs = 0.0
       zPs = 0.0

c      1 = first rupture cell, top
c      2 = last rupture cell, top
c      3 = first rupture cell, bottom
c      4 = last rupture cell, bottom
       
       sum1 = 0.0
       do iGC2=1, n3-1
         segX = rup_x(iGC2+1) - rup_x(iGC2)
         segY = rup_y(iGC2+1) - rup_y(iGC2)                  
         theta(iGC2) = atan2(segY,segX) 
         sum1 = sum1 + (Seg_length(iGC2) / GC2_ruplength * theta(iGC2)) 
       enddo
       Ave_theta = sum1
       
       J = sqrt(Global_T**2. + Global_U**2.)
       t = tan(abs(Global_T)/abs(Global_U))
       dx1 = (cos(t+Ave_theta))*J       
       dy1 = (sin(t+Ave_theta))*J       
       
       xPrup(1) = xPs - dx1      
       yPrup(1) = yPs - dy1
       zPrup(1) = fltgrid_z(iLocY,iLocX)
       
       dx2 = cos(Ave_theta)*GC2_ruplength
       dy2 = sin(Ave_theta)*GC2_ruplength
       
       xPrup(2) = xPrup(1) + dx2
       yPrup(2) = yPrup(1) + dy2
       zPrup(2) = fltgrid_z(iLocY,n2)
       
       z_rup = fltgrid_z(n1,iLocX)-fltgrid_z(iLocY,iLocX)
       cdip = 3.141592653590/2. - dipr
       dxy3 = tan(cdip)*z_rup
       dx3 = cos(Ave_theta - 3.141592653590/2.)*dxy3
       dy3 = sin(Ave_theta - 3.141592653590/2.)*dxy3
       
       xPrup(3) = xPrup(2) + dx3
       yPrup(3) = yPrup(2) + dy3
       zPrup(3) = fltgrid_z(n1,n2)
       
       xPrup(4) = xPrup(1) + dx3
       yPrup(4) = yPrup(1) + dy3
       zPrup(4) = fltgrid_z(n1,iLocX)
       
       xPrup(5) = xPrup(1) 
       yPrup(5) = yPrup(1)
       zPrup(5) = zPrup(1)
       
        return
       end             
