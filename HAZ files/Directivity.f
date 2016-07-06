      subroutine Directivity ( dirFlag, specT, Rrup, ZTOR, 
     1                 x0, y0, z0,
     1                 Rx, Ry, Ry0, mag, ftype, 
     2                 RupWidth, RupLength, Dipaverage, HWflag, dirMed, dirSigma, 
     3                 fltgrid_x, fltgrid_y, fltgrid_z,
     6                 n1, n2, iRupAS, iRupDD, dip, fs, fd, dpp_flag,
     7                 iLocAS, iLocDD)

c     3                 fltgrid_x1, fltgrid_y1, fltgrid_z1,
c     4                 fltgrid_x2, fltgrid_y2, fltgrid_z2,
c     5                 fltgrid_x3, fltgrid_y3,fltgrid_z3,
c     5                 fltgrid_x4, fltgrid_y4, fltgrid_z4, 

c      implicit none
      include 'pfrisk.h'
      
      integer MAXPERCY
      parameter (MAXPERCY=25)
      REAL Period(MAXPERCY), c8b(MAXPERCY), c8, c8a, c8bT
      real phi2_red(MAXPERCY)
      real fltGrid_x(MAXFLT_DD,MAXFLT_AS), fltGrid_y(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_z(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_x1(MAXFLT_DD,MAXFLT_AS), fltGrid_y1(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_z1(MAXFLT_DD,MAXFLT_AS), fltGrid_x2(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_y2(MAXFLT_DD,MAXFLT_AS), fltGrid_z2(MAXFLT_DD,MAXFLT_AS),
     3     fltGrid_x3(MAXFLT_DD,MAXFLT_AS), fltGrid_y3(MAXFLT_DD,MAXFLT_AS),
     4     fltGrid_z3(MAXFLT_DD,MAXFLT_AS), fltGrid_x4(MAXFLT_DD,MAXFLT_AS),
     5     fltGrid_y4(MAXFLT_DD,MAXFLT_AS), fltGrid_z4(MAXFLT_DD,MAXFLT_AS)
      real x0, y0, z0 
      real specT, Rrup, zTOR, Rx, Ry, Ry0, Mag, ftype
      real RupWidth, RupLength, Dipaverage, dirMed, dirSigma
      integer HWflag, dpp_flag
      real dip, len1, wid1, sum, sum1, sum2
      real fs, fd, aveDPP_est
      real hypoX, hypoY, hypoZ      
      real medadj, sigadj
      integer iRupAS, iRupDD, n1, n2
      integer dirFlag, iDD, iAS, iHypoAS, iHypoDD
      integer iDirectStrike, iDirectDip
      real cell_Width, cell_length, fractL, fractW
      real DPP, lnfd, lnfn, lnfp
      integer count1, count2
      
      data period     /
     1              0.0000, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500,
     1              0.0750, 0.1000, 0.1200, 0.1500, 0.1700,
     1              0.2000, 0.2500, 0.3000, 0.4000, 0.5000,
     1              0.7500, 1.0000, 1.5000, 2.0000, 3.0000,
     1              4.0000, 5.0000, 7.5000,10.0000/
      
      data c8b     /
     1             0.4833,  0.4833,  1.2144,  1.6421,  1.9456,  2.1810,
     1             2.6087,  2.9122,  3.1045,  3.3399,  3.4719,
     1             3.6434,  3.8787,  4.0711,  4.3745,  4.6099,
     1             5.0376,  5.3411,  5.7688,  6.0723,  6.5000,
     1             6.8035,  7.0389,  7.4666,  7.7700/
      data phi2_red /
     1             0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
     2             0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     3             -0.12, -0.15, -0.19, -0.23, -0.25 /
     
      c8 = 0.2154
      c8a = 0.2695
      nper = MAXPERCY

c     iRupY is the first cell in the rupture along strike
c     iRupZ is the first cell in the rupture down dip
c     n1 is the last cell that makes up the rupture plane down dip
c     n2 is the last cell that makes up the rupture plane along strike       
c      write (*,'( 2f10.2, 2x,''fs, fd'')') fs, fd
c      write (*,'( 3f10.3, 2x, '' mag, Rrup, ztor'' )') mag, Rrup, ztor
c      write (*,'( 6i5)') iRupAS, iRupDD, n1, n2
c      pause 'test iRupAS, iRupDD, n1, n2'
c      write (*,'( 3f10.3)') 
c     1       fltGrid_x1(iRupDD,iRupAS), fltGrid_y1(iRupDD,iRupAS), fltGrid_z1(iRupDD,iRupAS),
c     1       fltGrid_x2(iRupDD,iRupAS), fltGrid_y2(iRupDD,iRupAS), fltGrid_z2(iRupDD,iRupAS),
c     1       fltGrid_x3(iRupDD,iRupAS), fltGrid_y3(iRupDD,iRupAS), fltGrid_z3(iRupDD,iRupAS),
c     1       fltGrid_x4(iRupDD,iRupAS), fltGrid_y4(iRupDD,iRupAS), fltGrid_z4(iRupDD,iRupAS)
c
c      pause
      
c     Check the mag and period range for applying directivity      
c     **** Later, make these input parameters ****
      if (mag .lt. 5.6 .or. specT .lt. 0.50 ) return

c     Bayless and Somerville model, DIRFLAG = 30
      if (dirflag .eq. 30 ) then 
        len1 = RupLength * fs
      
c       Set X, Y, theta
        Y = fd

c       is the site along the rupture?
        if ( ry0 .eq. 0 ) then
          s = abs( Ry - (fs-0.5)*RupLength)
          x =  s / RupLength
          theta = atan(Rx/s)
        else
          s = len1
          x =  s / RupLength
          theta = atan( Rx/abs(s+Ry0) )
        endif

c       Change to degrees
        theta = 180./3.14 * theta
        
c       set rake
        if ( ftype .eq. 0. ) then
          rake = 0.
        elseif ( ftype .eq. 0.5 ) then
          rake = 45.
        elseif ( ftype .eq. 1. ) then
          rake = 90
        elseif ( ftype .eq. -0.5 ) then
          rake = -45.
        elseif ( ftype .eq. -1 ) then
          rake = -90.
        endif
        
c       Compute the source-to-site azimuth term (as defined by PEER)
        if ( Ry0 .eq. 0. ) then
          if ( Rx .ge. 0. ) then
            az = 90.
          else
            az = -90.
          endif
        else
          az = atan( Rx/Ry0) * 180./3.14 
        endif

c       compute Bayless model
        call ruptdirct2012_jrb ( specT, Rrup, mag, rupLength, rupWidth, ftype, 
     1                           theta, rake, Rx, X, Y, az, lnfd, lnfn, lnfp ) 
        medadj = lnfd
        
c       set sigma reduction based on JWL model (not used, but keep in case we want to add it)
        if (specT .eq. 0.0) then
          count1=1
          count2=1
        elseif (specT .gt. 0.0) then
          do i=2,nper-1
            if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
              count1 = i
              count2 = i+1
            endif
          enddo
        endif
        call interp (period(count1),period(count2),phi2_red(count1),phi2_red(count2),
     +                   specT,sigadj,iflag)

c      zero sigma adj for now    
       sigadj = 0.

c        write (*,'( 20f8.3)') mag, rupLength, rupWidth, ftype, theta, rake, Rx, X, Y, az, 
c     1        specT, fs, fd, lnfd, dirSigma, sigadj
      endif
       
c     JWL 2015 model, uniform Hypocenter model for SS and RV faults
c     Dirflag = 105
      if (dirflag .eq. 105) then
        call DirJWL_V3Uni (specT, rRup, Rx, Ry, RupLength, Mag, ftype, 
     1                    RupWidth, Dipaverage, HWflag, medadj, sigadj )
      endif

c     JWL 2015 model, preferred model for SS (SWUS, Appendix D) and RV (Chapter 3) faults
c     Dirflag = 106
      if (dirflag .eq. 106) then
        call DirJWL_V3Pre (specT, rRup, Rx, Ry, RupLength, Mag, ftype, 
     1                    RupWidth, Dipaverage, HWflag, medadj, sigadj )
c        write (*,'( 2f10.3)') medadj, sigadj
c        pause
      endif

c     Choiu and Spudich DPP model
c     Dirflag = 20
      if (dirflag .eq. 20) then

c       Compute average DPP
        call calc_aveDPP ( dpp_flag, mag, Rrup, Ztor, fs, fd, aveDPP_est ) 
c        write (*,'( 2x,''aveDPP ='',5f10.3)') fs, fd, mag, Rrup, aveDPP_est
c        pause

c       Check for bad value of aveDPP
        if ( aveDPP_est .lt. -900. ) then
          dirMed = 0
          dirSigma = 0
          return
        endif
        
c       Compute the coordinates of the hypocenter 
        len1 = RupLength * fs
        Wid1 = RupWidth * fd
        hypoX = len1 + fltGrid_x(iLocDD,iLocAS)
        hypoY = 0.
        hypoZ = Wid1 + fltGrid_z(iLocDD,iLocAS)

        
c        write (*,'( 2x,''Hypo:'',5f10.3)') hypox, hypoY, hypoZ, fs, rupLength
c        write (*,'( 4i5)') iLocAS, iLocDD, n1, n2
c        write (*,'( 2f10.4)') fltgrid_x(1,iLocAS), fltgrid_x(1,n2)

c       Calculate the DPP value
        call calc_DPP (hypoX, hypoY, hypoZ, 
     1           fltgrid_x, fltgrid_y, fltgrid_z, x0, y0, z0, n2, n1, specT, 
     2           RupLength, rupWidth, ftype,DPP, iLocAS, iLocDD)  
c        write (*,'( 2x,''DPP='',e12.4)') DPP
c        pause 'DPP'     
        
      if (specT .eq. 0.0) then
       count1=1
       count2=1
      elseif (specT .gt. 0.0) then
       do i=2,nper-1
        if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
         count1 = i
         count2 = i+1
        endif
       enddo
      endif
c      write (*,'( f10.4,3i5)') specT, nper, count1, count2
c      pause 'counts'
      
      call interp (period(count1),period(count2),c8b(count1),c8b(count2),
     +                   specT,c8bT,iflag)
      call interp (period(count1),period(count2),phi2_red(count1),phi2_red(count2),
     +                   specT,sigadj,iflag)
 
        cDPP = (DPP - aveDPP_est)
        medadj = c8 * exp(-c8a * (Mag-c8bT)**2) *
     1       max(0.0, 1.0-max(0.0,Rrup-40.0)/30.0) *
     1       min(max(0.0,Mag-5.5)/0.8, 1.0) * cDPP
     
c       write (*,'( 5f10.4)') c8, c8a, c8bT, mag, rrup 
c       write (*,'( 10f10.4)') cDPP, DPP, aveDPP_est, medadj, sigadj
c       pause



      endif

      dirMed = medadj
      dirSigma = sigadj


      return
      end

c —————————————————————————————
      subroutine calc_aveDPP ( dpp_flag, mag, Rrup, Ztor, fs, fd, 
     1              aveDPP_est ) 
      implicit none
      integer dpp_flag, nMag(10)
      integer i1, i2, i3, i4, i5, i, k
      integer n_fs, n_fd, nZTOR, ndpp_R
      integer iMag
      real aveDPP(10,15,10,10,31), Mag1(15), fd1(10)
      real fs1(10), R1(31), ztor1(10)
      real Len, Width, rake,dip, mag, rRup, zTor, fs, fd, aveDPP_est
      real fs2
      character*80 dummy
      common /Save_AVEDPP/ n_fs, n_fd, nZtor, nMag, ndpp_R,
     1       aveDPP, Mag1, fd1,
     2       fs1, R1, ztor1
     
c      write (*,'( 5f10.3)') mag, Rrup, Ztor, fs, fd
c      pause 'into aveDPP'

c     Check if the aveDPP file has been read into memory      
      if ( dpp_flag .eq. 0 ) then
        open (44,file='SS13_aveDPP1.txt',status='old')
        read (44,*) n_fd, n_fs, nZtor
        read (44,*) (nMag(k),k=1,nZtor)
        read (44,*) ndpp_R, (R1(k),k=1,ndpp_R)
        read (44,'( a80)') dummy
        do i1=1,nZtor
         do i2=1,nMag(i1)
          do i3=1,n_fs
           do i4=1,n_fd
            read (44,*) i, mag1(i2),Len, Width, ztor1(i1), rake,
     1          dip, fs1(i3), fd1(i4), 
     2          (aveDPP(i1,i2,i3,i4,i5),i5=1,ndpp_r) 
           enddo
          enddo
         enddo
        enddo

c       reset flag
        dpp_flag = 1
      endif
        
c     Interpolate to get DPP

c     Find the ztor index
      if ( ztor .eq. 0 ) then
        i1 = 1
      elseif ( ztor .gt. ztor1(nZtor) ) then
        i1 = nZtor
      else
        do i1=2,nZtor
          if ( zTOR .lt. ztor1(i1) ) goto 10
        enddo
      endif
 10   continue

c     Find the mag index
      imag = int((mag-5.6)/0.2 ) + 1
      if (imag .gt. nMag(i1)) then
        i2 = nMag(i1)
      else
        i2 = imag
      endif

c     Find the fs index
      fs2 = fs
      if (fs .gt. 0.5) fs2 = 1.-fs
      i3 = int(fs/0.1 +0.5)
      if ( i3 .gt. n_fs ) i3=n_fs

c     Find the fd index
      i4 = int(fd/0.1 +0.5)
      if ( i4 .gt. n_fd ) i4=n_fd
c      write (*,'( 2x,''nfd, nfs'',2i5)') n_fd, n_fs

c     find the distance index
c      write (*,'( 5f10.3)') (R1(k),k=1,ndpp_r)
c      pause 'dist check'     
      if ( rRUp .lt. R1(1) ) then
        i5 = 1
      elseif ( rRup .gt. R1(ndpp_r) ) then
        i5 = ndpp_r
      else
        do i5=1,ndpp_r
          if ( rRup .lt. r1(i5) ) goto 14
        enddo
      endif
 14   continue 

      aveDPP_est = aveDPP(i1,i2,i3,i4,i5) 
c      write (*,'( 5i5,6f10.3)') i1,i2,i3,i4,i5, mag, rrup, ztor, fs, fd
c      write (*,'( 2x,''aveDPP ='',f10.3)') aveDPP(i1,i2,i3,i4,i5)  
        
      return
      end               

c ------------------------------------------------------------------
