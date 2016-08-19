c     calc_DPP

c      Subroutine created to calculate the DPP parameter defined in Spudich  al. (2013)

c      Inputs
c      fltgrid_x, fltgrid_y, fltgrid_z  - matrices of x, y and z-axis locations in km of 
c       one grid square of the gridded rupture.
c      x0, y0, z0 - x, y and z-axis locations in km of the site.
c      iHypoX, iHypoZ - down-dip and along-strike grid square where the hypocenter is located.
c      icellRupStrike, icellRupDip - down-dip and along-strike grid square where the closest grid square to the site.
c      period - period in sec for which DDP should be calculated.
c      ftype - fault type, SS=0
c      n1, n2 - number of downdip and along strike grid cells
c      rupLen, rupWidth - rupture length and rupture width

c      Outputs
c      iDirectStrike, iDirectDip - down-dip and along-strike grid square of direct point.
c      DPP

c      Checks needed
c      ftype - Is SS ftype 0?

c      April 10, 2016


c -------------------------------------------------------------------      subroutine calc_DPP (iHypoX, iHypoZ, icellRupStrike, icellRupdip, 

      subroutine calc_DPP (hypoX, hypoY, hypoZ, 
     1           fltgrid_x, fltgrid_y, fltgrid_z, x0, y0, z0, n2, n1, specT, 
     2           RupLength, rupWidth, ftype,DPP, iLocAS, iLocDD)  

      implicit none

      include 'pfrisk.h'
      
      integer iHypoX, iHypoZ
      integer n2, n1
      integer nStrike
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS),
     1      fltGrid_y(MAXFLT_DD,MAXFLT_AS)
      real  ctildap, Rhypo, Rd, f, FS_bar, velrat, DPP, period, ftype
      real Ix, In, Iphi, L2, E, cosPhi, sinPhi, zs
      real hypoX, hypoY, hypoZ, x0, y0, z0, specT
      real rupLength, rupWidth
      integer iLocAS, iLocDD
      real pp_x, pp_y, pp_z, pd_x, pd_y, pd_z
      real xleft, xright
      
c     For now, simplfy to a straight fault with strike along the x-axis and a 90 degree dip

c     Set coordinates of point Pp
      Pp_x = 0.
      Pp_y = 0.
      Pp_z = 0.
            
c     Set coordinates of the direct point

c     First, check if the rupture is past the site 
      xLeft = fltgrid_x(1,iLocAS) 
      xright = fltgrid_x(1,n2) 
      if (xLeft .le. x0 .and. xRight .ge. x0 ) then
         pd_x = 0.
         pd_y = 0.
         pd_z = fltgrid_z(iLocDD,1)
      elseif (xRight .le. x0 .and. xLeft .ge. x0 ) then
         pd_x = 0.
         pd_y = 0.
         pd_z = fltgrid_z(iLocDD,1)
      else
        if ( xRight .lt. x0 ) then
          pd_x = fltgrid_x(1,n2)
          pd_y = 0.
          pd_z = hypoZ - hypoZ*abs(hypoX-fltgrid_x(iLocDD,n2)) /
     1         abs( hypoX-x0)
        else
          pd_x = fltgrid_x(1,iLocAS)
          pd_y = 0.
          pd_z = hypoZ - hypoZ*abs(hypoX-fltgrid_x(iLocDD,iLocAS)) /
     1         abs( hypoX-x0)
        endif
      endif

c     Find the E
      E = sqrt( (hypoX - pd_x)**2 + (hypoY - pd_y)**2 + (hypoZ - pd_z)**2 )
c     add 0.5 km for the grid spacing to correct for use of centers
      E = E + 0.5      

C     find the FS (simplifed geometry used here)
      zs = fltgrid_y(iLocDD,iLocAS)
      L2 = sqrt( (hypoX - pp_x)**2 + (hypoY - pp_y)**2 + (hypoZ - pp_z)**2 )
      Rd = sqrt( (L2-E)**2 + zs**2)
      Rhypo = sqrt( L2**2 + zs**2)
      cosPhi = (pd_x-hypoX) / E
      sinPhi = (pd_z-hypoZ) / E

c     Equation 6.7
      Ix = cosPhi * ( 2*Zs * (L2/Rhypo - (L2-E)/Rd) 
     1      - zs * alog ( (L2 + Rhypo) / (L2 - E + Rd) ) )
      
c     Equation 6.8
      In = cosPhi * (-2 * zs**2 * ( 1/RHypo - 1/Rd ) - ( Rhypo - Rd ) )
      
c     Equation 6.9
      if (cosphi .ge. 0.999999) then
       Iphi = 0
      else
       Iphi = sinPhi * ( zs * alog ( (L2 + Rhypo) / (L2 - E + Rd) ) )
      endif

c     Equation 6.11
      FS_bar = sqrt ( Ix**2 + In**2 + Iphi**2 ) / E

C     Compute c tilda prime (eq 4 of spudich and chiou).
      velrat = 0.80

      if (E .eq. 0.0) then
         ctildap = velrat
         E = 0.5
      else
         ctildap = ((1.0/velrat) - ((RHypo-Rd)/E))
         ctildap = 1.0/ctildap
      endif      
      
c      f = max(rupLength, rupWidth)
      
C     Compute DPP prime (Eq. 6.12).
      DPP = alog ( ctildap * max(E,0.1*rupLength) * max(FS_bar, 0.2) )
      
c      write (*,'( 2x,''E, FS, DPP'',3f10.3)') E, FS_bar, DPP

      return
      end

   
c--------------------------------------------------------------------
