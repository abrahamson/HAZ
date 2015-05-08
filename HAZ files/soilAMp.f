      subroutine ApplySoilAmp (mag, lnY, specT, 
     1      refPeriod, nRefPer, RefGM, nRefGM,
     2      RefGM_mag, nRefMag, amp, lgSoilSa )

      include 'pfrisk.h'

      integer nRefPer, nRefGM, nRefMag, iflag
      real mag, lnY, specT, lgSoilSa
      real refPeriod(1), RefGM_Mag(1)
      real RefGM(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM)
      real amp(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM)
      real amp3 
      integer iMag, igm, iper
      real amp1(MAX_AMPPER,MAX_AMPGM), amp2(MAX_AMPGM),
     1     gm1(MAX_AMPPER,MAX_AMPGM), gm2(MAX_AMPGM)
      
c     Note: amp factor and ref GM are in log units. 
      gm = exp(lnY)   
      

c      do iMag=1,nRefMag
c        do igm=1,nRefGM
c          do iper=1,nRefPer
c            write (*,'( 3i5,4f10.4)') iMag, igm, iper, 
c     1          refGM(iMag,iper,igm), amp(iMag,iper,igm)
c          enddo
c        enddo
c      enddo
c      pause '13'
                                                   
c     Interpolate Amp factor and ref GM on magnitude
      if ( mag .lt. refGM_mag(1)) then
        iMag = 1
         do iper=1,nper
           do igm=1,nRefGM
             amp1(iper,igm) = amp(iMag,iper,iGM)
             gm1(iper,igm) = refgm(iMag,iper,iGM)
           enddo
         enddo
      elseif ( mag .gt. refGM_mag(nRefMag) ) then
        iMag = nRefMag
         do iper=1,nRefPer
           do igm=1,nRefGM
             amp1(iper,igm) = amp(iMag,iper,iGM)
             gm1(iper,igm) = refGM(iMag,iper,iGM)
           enddo
         enddo
      else
        iMag = 2
        do while (refGM_mag(iMag) .lt. mag ) 
          iMag = iMag + 1
        enddo
        do iper=1,nRefPer
           do igm=1,nRefGM
             amp1(iper,igm) = amp(iMag-1,iper,iGM) +
     1            (amp(iMag,iper,iGM)-amp(iMag-1,iper,iGM)) /
     2            (refGM_Mag(iMag) - refGM_Mag(iMag-1))
             gm1(iper,igm) = refGM(iMag-1,iper,iGM) +
     1            (refGM(iMag,iper,iGM)-refGM(iMag-1,iper,iGM)) /
     2            (refGM_Mag(iMag) - refGM_Mag(iMag-1))
           enddo
         enddo
      endif
      

c     Interpolate on log period                                                                                                                                                   
C     Find index of the requested spectral period.
C     First Check for the PGA case (i.e., specT=0.0)
      if (specT .eq. 0.0) then
        do iGM=1,nRefGM
          iper = 1
          amp2(iGM) = amp1(iper,iGM)
          gm2(iGM) = gm1(iper,iGM)
        enddo 
      elseif ( specT .gt. refPeriod(nRefPer) ) then
        do iGM=1,nRefGM
          iper = nRefPer
          amp2(iGM) = amp1(iper,iGM)
          gm2(iGM) = gm1(iper,iGM)
        enddo 
      else 
        do i=1, nRefPer-1
          if (specT .ge. refPeriod(i) .and. specT .le. refPeriod(i+1) ) then
            do iGM=1,nRefGM
              call interp ( refPeriod(i), refPeriod(i+1), 
     1                           amp1(i,iGM), amp1(i+1,iGM), 
     2                           specT, amp2(iGM), iflag )
c              write (*,'(10 f10.4)') refPeriod(i), refPeriod(i+1), 
c     1                           amp1(i,iGM), amp1(i+1,iGM), 
c     2                           specT, amp2(iGM)
              call interp ( refPeriod(i), refPeriod(i+1), 
     1                           gm1(i,iGM), gm1(i+1,iGM), 
     2                           specT, gm2(iGM), iflag )
            enddo 
          endif
         enddo
       endif
       
c       do igm=1,nRefGM
c         write (*,'( i5,2f10.4)') igm, amp2(iGm), gm2(iGM)
c       enddo
       
c       write (*,'( f10.4)') gm

C     Interpolate the amp for the specified GM 
c     level.
      if ( gm .le. gm2(1) ) then
c        write (*,'( 2x,'' small'',2f10.4)') gm, gm2(1)
        amp3 = amp2(1)
      elseif ( gm .ge. gm2(nRefGM)) then
c         write (*,'( 2x,'' large'',2f10.4)') gm, gm2(nRefGM)
       amp3 = amp2(nRefGM)
      else
        do iGM=1, nRefGM-1
          if ( gm .ge. gm2(iGM) .and. gm .le. gm2(iGM+1)) then
c            write (*,'( 2x,''Call interp'')')
            call interp ( gm2(iGM), gm2(iGM+1), amp2(iGM), 
     1           amp2(iGM+1), gm, amp3, iflag ) 
          endif
        enddo
      endif
c      write (*,'( 5f10.4)') mag, gm, specT, amp3
c      pause 'amp3'

C     Now apply the rock to soil amp factor.
      lgSoilSa = amp3 + LnY

      return
      end
      
c -------------------------------------------------

      subroutine RdSoilAmpModel ( refPeriod, nRefPer, RefGM, nRefGM,
     2      RefGM_mag, nRefMag, amp )    

      include 'pfrisk.h'

      integer nRefPer, nRefGM, nRefMag
c      real specT 
      real refPeriod(1), RefGM_Mag(1)
      real RefGM(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM)
      real amp(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM)
c      real amp1(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM)
      real pga(MAX_AMPGM), refShape(MAX_AMPMAG,MAX_AMPPER)
      character*80 filein

c     Open file      
      read (13,'( a80)') filein
c      write (*,'( a80)') filein
      open (31,file=filein,status='old')
      rewind (31)
      
c     Read soil amp file
      read (31,*) nRefMag, (RefGM_mag(k),k=1,nRefMag)
      read (31,*) nRefPer 
      read (31,*) nRefGM, (pga(k),k=1,nRefGM)
c      write (*,'( 3i5)') nRefMag, nRefPer, nRefGM
c      pause '10'

      do iMag =1,nRefMag
        do iper = 1,nRefPer
          read (31,*) refPeriod(iper), refShape(iMag,iper), 
     1    (amp(iMag,iper,igm),igm=1,nRefGM)
        enddo
      enddo
c      pause '11'
      
c     Convert shape to Sa and amp factor to log (amp)
      do iMag=1,nRefMag
        do igm=1,nRefGM
          do iper=1,nRefPer
            refGM(iMag,iper,igm) = refShape(iMag,iper)*pga(igm)
            amp(iMag,iper,igm) = alog( amp(iMag,iper,igm) )
c            write (*,'( 3i5,4f10.4)') iMag, igm, iper, 
c     1          refGM(iMag,iper,igm), amp(iMag,iper,igm)
          enddo
        enddo
      enddo
c      pause '12'

      close (31)
      
      return
      end             
      
