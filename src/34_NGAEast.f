c  ***** PEER NGAEast Models (2018) **********

c ---------------------------------------------------------------------
c ** Goulet et al., 2018 (PEER Report 2018/08) **
c    Central and Eastern North America Ground-Motion Characterization
c    NGA-East Final Report
c
c    Applicable Range:
c    M 4.0 - 8.2
c    Rrup < 1500 km
c
c    Reference Site Conditions:
c    VS30 = 3000 m/s, kappa = 0.006 s
c    
c    Appendix H (17 median models)
c ---------------------------------------------------------------------
      subroutine S34_NGAEast_Med ( m, Rrup, specT, imod, period2, lnY, iflag )

      implicit none
      include 'pfrisk.h'

      integer MAX_FREQ
      parameter (MAX_FREQ=25)
      integer iflag, iflagRead(17), imod, isc, nsc, ifreq, f1, f2, m1d1, m1d2, m2d1, m2d2
      real m, Rrup, specT, period2, lnY, freq(MAX_FREQ), table_m(MAX_EAST),
     1     table_r(MAX_EAST), table_lnSA(17,MAX_FREQ,MAX_EAST), dum, specF,
     2     lnSAf(4), lnSAfm(2), period1
      character*80 table, dummy, freq_ch(MAX_FREQ)
      character*10 number

      common/GMTable/ iflagRead, table_m, table_r, table_lnSA, nsc

      data freq / 0.100, 0.133, 0.200, 0.250, 0.333, 0.500, 0.667, 1.000, 1.333,
     1            2.000, 2.500, 3.333, 4.000, 5.000, 6.667, 10.000, 13.333,
     2            20.000, 25.000, 33.333, 40.000, 50.000, 100.000, 100.000, -1.0 /

      data freq_ch / 'F0.100', 'F0.133', 'F0.200', 'F0.250', 'F0.333', 'F0.500', 'F0.667',
     1            'F1.000', 'F1.333', 'F2.000', 'F2.500', 'F3.333', 'F4.000', 'F5.000',
     2            'F6.667', 'F10.000', 'F13.333', 'F20.000', 'F25.000', 'F33.333', 'F40.000',
     3            'F50.000', 'F100.000', 'PGA', 'PGV' /

c     number of scenarios in the NGA-East Tables (rows)
      nsc = 374

c     write the model number integer into a string
      write(number, '(i0)') imod

c     read in the NGA-East Tables
      if (iflagRead(imod) .eq. 0) then
        do ifreq=1,MAX_FREQ
          table = 'H.7.1_Median_Models/'//TRIM(freq_ch(ifreq))//'/NGA_East_Model_'//trim(adjustl(number))//'.csv'
          open (404,file=table,status='old')
          read (404,*) dummy
          do isc=1,nsc
            if (ifreq .eq. 1) then
              read (404,*) table_m(isc), table_r(isc), table_lnSA(imod,ifreq,isc)
            else
              read (404,*) dum, dum, table_lnSA(imod,ifreq,isc)
            endif
          enddo
          close (404)
        enddo
        iflagRead(imod) = 1
      endif

c     m or Rrup is outside range defined by ground motion model
      if (m .lt. 4. .or. m .gt. 8.2) then
        write (*,*)
        write (*,*) 'The NGA-East ground motion models'
        write (*,*) '(Goulet et al., 2018) are not defined'
        write (*,*) 'for a magnitude less than 4 or greater than 8.2 '
        write (*,'(a27,f10.5,a11)') 'The specific magnitude of: ', m, 'is invalid.'
        write (*,*) 'Please check the input file.'
        write (*,*)
        stop 99
      elseif (Rrup .gt. 1500.) then
        write (*,*)
        write (*,*) 'The NGA-East ground motion models'
        write (*,*) '(Goulet et al., 2018) are not defined'
        write (*,*) 'for Rrup greater than 1500 km'
        write (*,'(a22,f10.5,a14)') 'The specific Rrup of: ', Rrup, 'km is invalid.'
        write (*,*) 'Please check the input file.'
        write (*,*)
        stop 99
      endif

c     Check for PGA, PGV
      if (specT .eq. 0.0) then
        specF = freq(24)
        goto 1019
      elseif (specT .eq. -1.0) then
        specF = freq(25)
        goto 1019
      endif

c     find frequencies for interpolation
      if (specT .gt. 0.0) then
        specF = 1./specT
        f1 = 0
        f2 = 0
        do ifreq=1,MAX_FREQ
          if (specF .ge. freq(ifreq) .and. specF .le. freq(ifreq+1) ) then
            f1 = ifreq
            f2 = ifreq+1
            goto 1019
          endif
        enddo
      endif

c     specified spectral period is outside range defined by ground motion model
      write (*,*)
      write (*,*) 'The NGA-East ground motion models'
      write (*,*) '(Goulet et al., 2018) are not defined'
      write (*,'(a26,f10.5,a2)') ' for a spectral period of: ', specT, ' s'
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

c     find magnitudes and distances for interpolation
 1019 m1d1 = 0
      m1d2 = 0
      m2d1 = 0
      m2d2 = 0
      do isc=1,nsc
        if (m .ge. table_m(isc) .and. m .le. table_m(isc+34) ) then
          if (Rrup .ge. table_r(isc) .and. Rrup .le. table_r(isc+1)) then
            m1d1 = isc
            m1d2 = isc + 1
            m2d1 = isc + 34
            m2d2 = isc + 35
            goto 1020
          endif
        endif
      enddo

c     Set lnSAf for PGA, PGV
 1020 if (specT .eq. 0.0) then
        lnSAf(1) = table_lnSA(imod,24,m1d1)
        lnSAf(2) = table_lnSA(imod,24,m1d2)
        lnSAf(3) = table_lnSA(imod,24,m2d1)
        lnSAf(4) = table_lnSA(imod,24,m2d2)
        goto 1021
      elseif (specT .eq. -1.0) then
        lnSAf(1) = table_lnSA(imod,25,m1d1)
        lnSAf(2) = table_lnSA(imod,25,m1d2)
        lnSAf(3) = table_lnSA(imod,25,m2d1)
        lnSAf(4) = table_lnSA(imod,25,m2d2)
        goto 1021
      endif

c     interpolate between frequencies
      lnSAf(1) = ((((table_lnSA(imod,f1,m1d1)-table_lnSA(imod,f2,m1d1))/(alog(freq(f1))-alog(freq(f2))))
     1                *(alog(freq(f1))-alog(specF)))-table_lnSA(imod,f1,m1d1))*(-1)
      lnSAf(2) = ((((table_lnSA(imod,f1,m1d2)-table_lnSA(imod,f2,m1d2))/(alog(freq(f1))-alog(freq(f2))))
     1                *(alog(freq(f1))-alog(specF)))-table_lnSA(imod,f1,m1d2))*(-1)
      lnSAf(3) = ((((table_lnSA(imod,f1,m2d1)-table_lnSA(imod,f2,m2d1))/(alog(freq(f1))-alog(freq(f2))))
     1                *(alog(freq(f1))-alog(specF)))-table_lnSA(imod,f1,m2d1))*(-1)
      lnSAf(4) = ((((table_lnSA(imod,f1,m2d2)-table_lnSA(imod,f2,m2d2))/(alog(freq(f1))-alog(freq(f2))))
     1                *(alog(freq(f1))-alog(specF)))-table_lnSA(imod,f1,m2d2))*(-1)

c     interpolate between magnitudes
 1021 lnSAfm(1) = ((((lnSAf(1)-lnSAf(3))/(table_m(m1d1)-table_m(m2d1)))
     1                *(table_m(m1d1)-m))-lnSAf(1))*(-1)
      lnSAfm(2) = ((((lnSAf(2)-lnSAf(4))/(table_m(m1d2)-table_m(m2d2)))
     1                *(table_m(m1d2)-m))-lnSAf(2))*(-1)

c     interpolate between distances
      if (Rrup .le. 1.) then
        lnY = ((((lnSAfm(1)-lnSAfm(2))/(table_r(m1d1)-table_r(m1d2)))
     1                *(table_r(m1d1)-Rrup))-lnSAfm(1))*(-1)
      else
        lnY = ((((lnSAfm(1)-lnSAfm(2))/(alog(table_r(m1d1))-alog(table_r(m1d2))))
     1                *(alog(table_r(m1d1))-alog(Rrup)))-lnSAfm(1))*(-1)
      endif

c     convert ground motion to units of gals (ok even for PGV, gets converted back)
      lnY = lnY + 6.89

      period1 = specT
      period2 = period1

      return
      end
