
c --------------------------------------------------------------

      subroutine S23_InitMinDis (iFlt, nWidth, MinRrup_temp, MinRjb_temp,
     1                       MinSeismo_temp, SourceDist)

      implicit none
      include 'pfrisk.h'

      integer iFlt, nWidth, iWidth, j
      real MinRrup_temp, MinRjb_temp, MinSeismo_temp,
     1     SourceDist(MAX_FLT,MAX_WIDTH,3)

        MinRrup_temp = 99999
        MinRjb_temp = 1.e10
        MinSeismo_temp = 1.e10
        do iWidth=1,nWidth
          do j=1,3
            SourceDist(iFlt,iWidth,j) = 1.e10
          enddo
        enddo

      return
      end

c --------------------------------------------------------------

      subroutine S23_InitHaz ( Haz )

      implicit none
      include 'pfrisk.h'

      real*8 Haz( MAX_INTEN, MAX_PROB, MAX_FLT  )
      integer i2, i3, i4

      do 90 i2 = 1, MAX_PROB
        do 80 i3 = 1, MAX_FLT
          do 70 i4 = 1, MAX_INTEN
                Haz( i4, i2, i3) = 0.0
  70      continue
  80    continue
  90  continue
      return
      end

c --------------------------------------------------------------

      subroutine S23_InitHaztree ( Haz )

      implicit none
      include 'pfrisk.h'

      real*8 Haz( MAX_INTEN, MAX_PROB, MAX_ATTEN  )
      integer i2, i3, i4

      do 90 i2 = 1, MAX_PROB
        do 80 i3 = 1, MAX_ATTEN
          do 70 i4 = 1, MAX_INTEN
                Haz( i4, i2, i3) = 0.0
  70      continue
  80    continue
  90  continue
      return
      end

c --------------------------------------------------------------

      subroutine S23_InitHazBins ( HazBins )

      implicit none
      include 'pfrisk.h'

      real HazBins( MAX_MAG, MAX_DIST, MAX_EPS,MAX_PROB, MAX_INTEN )
      integer i1, i2, i2b, i3, i4

      do 90 i1 = 1, MAX_MAG
        do 85 i2 = 1, MAX_DIST
         do 80 i2b = 1, MAX_EPS
          do 70 i3 = 1, MAX_PROB
            do i4 = 1,MAX_INTEN
              HazBins( i1, i2, i2b, i3, i4) = 0.0
            enddo
  70      continue
  80     continue
  85    continue
  90  continue
      return
      end

c --------------------------------------------------------------

      subroutine S23_InitHazBinsX ( HazBinsX )

      implicit none
      include 'pfrisk.h'

      real HazBinsX(MAX_XCost,MAX_PROB, MAX_INTEN )
      integer i2b, i3, i4

         do 80 i2b = 1, MAX_XCost
          do 70 i3 = 1, MAX_PROB
            do i4 = 1,  MAX_INTEN
              HazBinsX(i2b, i3, i4) = 0.0
            enddo
  70      continue
  80     continue
      return
      end

c --------------------------------------------------------------

      subroutine S23_initRup ( sigma, nValue, sigmaMax, step, iFlt)

      implicit none

      real sigma(1), step, sigmaMax
      integer nValue(1), iFlt

c     For Small Sigma, only one rupture length is needed
      IF ( sigma(iFlt) .LE. 0.02 ) nValue(iFlt) = 1

c     Set Rupture Length Integration parameters
c     (maximum number of standard deviations and step size)
      sigmaMax = 2.0
      step = 2.*sigmaMax / nValue(iFlt)

      return
      end

c --------------------------------------------------------------

      subroutine S23_SetBin ( nBins, bins, x, iBin )

      implicit none

      integer nBins, iBin, i
      real bins(1), x

      ibin = 0

      do i=2,nBins
        if ( x .le. bins(i) .and. x .ge. bins(i-1)) then
          iBin = i-1
          return
        endif
      enddo

c     Reset Bin Range to include testing value.
c     First check for value less than lowest bin value.
      if (ibin .eq. 0 .and. x .lt. bins(1) ) then
         bins(1) = x
         iBin = 1
      endif

c     Now check for value greater than highest bin value.
      if (ibin .eq. 0 .and. x .gt. bins(nBins) ) then
         bins(nBins) = x
         ibin = nBins - 1
      endif

      return
      end

c --------------------------------------------------------------

      subroutine S23_SetBin_Eps ( nBins, bins, x, iBin )

      implicit none

      integer nBins, iBin, i
      real bins(1), x

      ibin = 0

      do i=2,nBins
        if ( x .le. bins(i) .and. x .ge. bins(i-1)) then
          iBin = i-1
          return
        endif
      enddo

c     Reset Bin Range to include testing value.
c     First check for value less than lowest bin value.
      if (ibin .eq. 0 .and. x .lt. bins(1) ) then
         bins(1) = x
         iBin =1
      endif

c     Now check for value greater than highest bin value.
      if (ibin .eq. 0 .and. x .gt. bins(nBins) ) then
         bins(nBins) = x
         iBin = nBins - 1
      endif

      return
      end

c -------------------------------------------------------------

      subroutine S23_Init_tempHaz ( tempHaz )

      implicit none
      include 'pfrisk.h'

      real*8 tempHaz(MAXPARAM,MAX_INTEN, MAX_PROB, MAX_ATTEN,MAX_FTYPE)
      integer iParam, iInten, iProb, iAtten, iFtype

      do iParam=1,MAXPARAM
          do iInten=1,MAX_INTEN
            do iProb=1,MAX_PROB
              do iAtten=1,MAX_ATTEN
                do iFtype=1,MAX_FTYPE
                  tempHaz(iParam,iInten,iProb,iAtten,iFtype) = 0.0
                enddo
              enddo
            enddo
          enddo
      enddo
      return
      end

c -------------------------------------------------------------

      subroutine S23_Init_tempHaz1 ( tempHaz1 )

      implicit none
      include 'pfrisk.h'

      real*8 tempHaz1(MAXPARAM,MAX_INTEN, MAX_PROB,MAX_FTYPE)
      integer iParam, iInten, iProb, iFtype

      do iParam=1,MAXPARAM
          do iInten=1,MAX_INTEN
            do iProb=1,MAX_PROB
                do iFtype=1,MAX_FTYPE
                  tempHaz1(iParam,iInten,iProb,iFtype) = 0.0
                enddo
            enddo
          enddo
      enddo
      return
      end

c -------------------------------------------------------------

      subroutine S23_Init_tempHaz2 ( tempHaz2 )

      implicit none
      include 'pfrisk.h'

      real*8 tempHaz2(4, MAX_INTEN, MAX_PROB, MAX_ATTEN)
      integer i, iInten, iProb, iAtten

      do i=1,4
       do iInten=1,MAX_INTEN
         do iProb=1,MAX_PROB
           do iAtten=1,MAX_ATTEN
             tempHaz2(i, iInten,iProb,iAtten) = 0.0
           enddo
         enddo
       enddo
      enddo
      return
      end

c -------------------------------------------------------------

      subroutine S23_InitDeagg ( m_bar, d_bar, e_bar, Xcost_bar )

      implicit none
      include 'pfrisk.h'

      real*8 m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN)
      real*8 Xcost_bar(MAX_PROB,MAX_INTEN)
      real*8 e_bar(MAX_PROB,MAX_INTEN)
      integer i, ii

      do i=1,MAX_PROB
         do ii=1,MAX_INTEN
           m_bar(i,ii) = 0.
           d_bar(i,ii) = 0.
           e_bar(i,ii) = 0.
           Xcost_bar(i,ii) = 0.
         enddo
      enddo
      return
      end

c -------------------------------------------------------------

      subroutine S23_InitDeagg2 ( m_bar_s, rrup_bar_s, rjb_bar_s, rx_bar_s, e_bar_s,
     1                        rSeismo_bar_s, ry0_bar_s, mag_bar_s, ftype_bar_s,
     2                        hypodepth_bar_s, dipavgd_bar_s, ztor_bar_s,
     3                        thetasite_bar_s, rupwidth_bar_s, rHypo_bar_s )

      implicit none
      include 'pfrisk.h'

      real*8 m_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN), rrup_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     1       rjb_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN), rx_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN),
     2       e_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), rSeismo_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     3       ry0_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), mag_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     4       ftype_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), hypodepth_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     5       dipavgd_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), ztor_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     6       thetasite_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), rupwidth_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     7       rHypo_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN)


      integer iFlt, i, ii

      do iFlt=1,MAX_FLT
        do i=1,MAX_PROB
          do ii=1,MAX_INTEN
            m_bar_s(iFlt,i,ii) = 0.
            rrup_bar_s(iFlt,i,ii) = 0.
            rjb_bar_s(iFlt,i,ii) = 0.
            rx_bar_s(iFlt,i,ii) = 0.
            e_bar_s(iFlt,i,ii) = 0.
            rSeismo_bar_s(iFlt,i,ii) = 0.
            ry0_bar_s(iFlt,i,ii) = 0.
            mag_bar_s(iFlt,i,ii) = 0.
            ftype_bar_s(iFlt,i,ii) = 0.
            hypodepth_bar_s(iFlt,i,ii) = 0.
            dipavgd_bar_s(iFlt,i,ii) = 0.
            ztor_bar_s(iFlt,i,ii) = 0.
            thetasite_bar_s(iFlt,i,ii) = 0.
            rupwidth_bar_s(iFlt,i,ii) = 0.
            rHypo_bar_s(iFlt,i,ii) = 0.
          enddo
        enddo
      enddo
      return
      end

c -------------------------------------------------------------
