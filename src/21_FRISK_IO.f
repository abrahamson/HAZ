
      subroutine S21_RdInput ( nProb, nAttenType, nAtten, jcalc, specT, sigTrunc,
     1               gmScale, dirFlag, PCflag, nInten, testInten, lgTestInten,
     2               psCorFlag, minlat, maxlat, minlong, maxlong, distmax,
     3               nMagBins, magBins, nDistBins, distBins, nepsBins, epsBins,
     4               nXcostBins, xcostBins, soilAmpFlag, gm_wt, sigvaradd,
     5               sCalc, sigfix, bnumflag, cfcoefrrup, cfcoefrjb,
     6               coefcountrrup, coefcountRjb, iMixture, version, starttime )

      implicit none
      include 'pfrisk.h'

      integer nInten(MAX_PROB), runFlag, bnumflag, nMagBins, nDistBins,
     1        nEpsBins, nXcostBins, nProb, nwr, psCorFlag, ichart,
     2        dirflag(MAX_PROB), jCalc(MAX_PROB,4,MAX_ATTEN), i, j,
     3        sCalc(MAX_PROB,4,MAX_ATTEN), nAttenType, j1, k, jj,
     4        nAtten(MAX_PROB,4), soilAmpFlag, charcount
      integer coefcountrrup, coefcountrjb, PCflag(MAX_PROB),
     1        iMixture(MAX_PROB,4,MAX_ATTEN)
      integer*4 starttime(3)
      real testInten(MAX_PROB,MAX_INTEN), distmax, minlat, maxlat,
     1     minlong, maxlong, lgTestInten(MAX_PROB,MAX_INTEN),
     2     distBins(MAX_DIST), epsBins(MAX_EPS), XcostBins(MAX_Xcost),
     3     sigTrunc(MAX_PROB), sigfix(MAX_PROB,4,MAX_ATTEN),
     4     magBins(MAX_MAG), specT(MAX_PROB), gm1, gm2
      real gmScale(MAX_PROB,4,MAX_ATTEN), gm_wt(MAX_PROB,4,MAX_ATTEN),
     1     sigvaradd(MAX_PROB,4,MAX_ATTEN), cfcoefrrup(MAX_Atten,11),
     2     cfcoefrjb(MAX_Atten,11), checkattenwt, version
      character*80 filein, title, file1, filelog

c     Set Data file units
      nwr = 11

C     Run Ground Motion Case or PSHA case.
 123  write (*,*) 'Enter 0 to run hazard'
      write (*,*) 'Enter 1 to run attenuation models'
      write (*,*) 'Enter 2 to run spectra models'

      if (bnumflag .eq. 0) then
         read (*,*) runFlag
      else
         runflag = 0
      endif

c     Open Input and Output Files
      if (runFlag .eq. 0) then
         write (*,*) 'Enter the input file name for PSHA runs.'
         if (bnumflag .eq. 0 ) then
           read (*,'( A80)') file1
         else
           read (77, '(a80)') file1
         endif

         open (13,file=file1,status='old',err=2100)

c        Create log filename for PHSA run.
         charcount = 80
         do ichart=80,1,-1
           if (file1(ichart:ichart) .ne. ' ') then
             charcount = ichart
             goto 111
           endif
         enddo
 111     continue
         do ichart=1,charcount,1
           filelog(ichart:ichart) = file1(ichart:ichart)
         enddo
         filelog(charcount+1:charcount+1) = '.'
         filelog(charcount+2:charcount+2) = 'l'
         filelog(charcount+3:charcount+3) = 'o'
         filelog(charcount+4:charcount+4) = 'g'
         open (18,file=filelog(1:charcount+4))
      elseif (runFlag .eq. 1) then
         write (*,*) 'Enter the input file name for attenuation runs.'
         read (*,'( A80)') file1
         open (13,file=file1,status='old')
         call S15_determ (runFlag)
         stop
      elseif (runFlag .eq. 2) then
         write (*,*) 'Enter the input file name for spectra runs.'
         read (*,'( A80)') file1
         open (13,file=file1,status='old')
         call S15_determ (runFlag)
         stop
      else
         write (*,*) 'Invalid entry!!!'
         write (*,*) 'Try again.'
         goto 123
      endif

c     Read fault file name and version format
      read (13,'( a80)') filein
      open (10,file=filein,status='old',err=2101)
      read (13,*,err=2020) version

c     Read min max long and lat for grid sources
      read (13,*,err=2001) minlat, maxlat, minlong, maxlong
      read (13,*,err=2002) distmax

c     Input options
      read(13,'( a80)') title
      read(13,*,err=2003) nProb, nAttenType

      call S21_CheckDim ( nProb, MAX_Prob, 'MAX_PROB ' )
      call S21_CheckDim ( nProb, MAX_Prob, 'MAX_ATTEN ' )
      write (18,*) '45.3 Haz45.3 log file'
      write (18,'( a80)') title
      write (18,*)
      write (18,*) 'Summary of Attenuation relationship used in the analysis:'
      write (18,*)
      write (18,'(a20,2x,i4)') 'Number of Problems: ', nProb
      write (18,'(a29,2x,i4)') 'Number of Attenuation Types: ', nAttentype
      write (18,*)

c     Enter GM models for each problem
      do i=1,nProb
        read (13,*,err=2004) specT(i), sigTrunc(i), dirflag(i), PCflag(i)
        call S16_CheckDir (dirflag(i))
        call S21_CheckPC (PCflag(i))
        read (13,*,err=2005) nInten(i), (testInten(i,j),j=1,nInten(i))
        call S21_CheckDim ( nInten(i), MAX_INTEN, 'MAX_INTEN ' )
        do j1=1,nInten(i)
          lgTestInten(i,j1) = alog(testInten(i,j1))
        enddo

        write (18,'(a19,2x,f8.4)') 'Spectral Period  = ', specT(i)
        write (18,'(a19,2x,f8.4)') 'Sigma Trunc      = ', sigTrunc(i)
        write (18,'(a19,2x,i8)')   'Directivity Flag = ', dirflag(i)
        write (18,'(a19,2x,i8)')    'Poly Chaos Flag = ', PCflag(i)
        write (18,'(a21)')         'Ground Motion Values:'
        do j1=1,nInten(i)
          write (18,'(f12.4)') testInten(i,j1)
        enddo
        do j=1,nAttenType
          write (18,*)
          write (18,'(a19,2x,i4)') 'Attenuation Type = ', j
          write (18,*) ' Jcalc    Const1    Const2    Weight    SigVaradd    iMix    j_sigma    fix_sigma'

          checkattenwt = 0.0
          read (13,*,err=2006) nAtten(i,j)
          call S21_CheckDim ( nAtten(i,j), MAX_ATTEN, 'MAX_ATTEN ' )
          do k=1,nAtten(i,j)
            read (13,*,err=2007) jCalc(i,j,k), gm1, gm2, gm_wt(i,j,k), sigvaradd(i,j,k), iMixture(i,j,k)

C           Check for either fixed sigma or different sigma model (i.e., jCalc<0)
            if (jcalc(i,j,k) .lt. 0) then
              backspace (13)
              read (13,*,err=2006) jCalc(i,j,k), gm1, gm2, gm_wt(i,j,k),
     1              sigvaradd(i,j,k), iMixture(i,j,k), sCalc(i,j,k), sigfix(i,j,k)
              write (18,'(i6,3x,f8.4,2x,f8.4,1x,f9.5,f11.4,2x,i6,4x,i6,6x,f8.4)')
     1              jCalc(i,j,k), gm1, gm2, gm_wt(i,j,k),
     1              sigvaradd(i,j,k), iMixture(i,j,k), sCalc(i,j,k), sigfix(i,j,k)
            else
              write (18,'(i6,3x,f8.4,2x,f8.4,2x,f9.5,f11.4,2x,i6)') jCalc(i,j,k), gm1, gm2, gm_wt(i,j,k),
     1              sigvaradd(i,j,k), iMixture(i,j,k)
            endif

C     Check for Common Functional Form with Rrup Distance (10000<jcalc<11000) selected and if so read in coefficients.
            if (abs(jcalc(i,j,k)) .gt. 10000 .and. abs(jcalc(i,j,k)) .lt. 11000) then
                coefcountrrup = abs(jcalc(i,j,k)) - 10000
                read (13,*,err=2008) (cfcoefrrup(coefcountrrup,jj),jj=1,11)
                write (18,'(a16,2xi8,11(2x,f12.6))')
     1             'CF(RRup) Coefs: ', coefcountrrup,(cfcoefrrup(coefcountrrup,jj),jj=1,11)
            endif
C     Check for Common Functional Form with RJB Distance (11000<jcalc<12000) selected and if so read in coefficients.
            if (abs(jcalc(i,j,k)) .gt. 11000 .and. abs(jcalc(i,j,k)) .lt. 12000) then
                coefcountrjb = abs(jcalc(i,j,k)) - 11000
                read (13,*,err=2008) (cfcoefrjb(coefcountrjb,jj),jj=1,11)
                write (18,'(a16,2x,i8,11(2x,f12.6))')
     1              'CF(RJB) Coefs:  ', coefcountrjb, (cfcoefrrup(coefcountrjb,jj),jj=1,11)
            endif
C     Check for DCPP Common Functional Form with Rrup Distance (12000<jcalc<13000) selected and if so read in coefficients.
            if (abs(jcalc(i,j,k)) .gt. 12000 .and. abs(jcalc(i,j,k)) .lt. 13000) then
                coefcountrrup = abs(jcalc(i,j,k)) - 12000
                read (13,*,err=2008) (cfcoefrrup(coefcountrrup,jj),jj=1,11)
                write (18,'(a16,2x,i8,11(2x,f12.6))')
     1              'CF(RRup) Coefs: ', coefcountrrup,(cfcoefrrup(coefcountrrup,jj),jj=1,11)
            endif

            gmScale(i,j,k) = gm1 + gm2
            checkattenwt = checkattenwt + gm_wt(i,j,k)
          enddo
          write (18,*)
C     Check that attenuation model weights sum to unity.
          if (abs(checkattenwt-1.0) .gt. 0.001 ) then
             write (*,*) 'Attenuation Model Weights do not sum to 1.0!!!'
             write (*,*) 'iProb =', i
             write (*,*) 'iAttenType = ', j
             write (*,*) 'Total Weights = ', checkattenwt
             write (*,*) 'Check input run file.'
             stop 99
          endif
        enddo
      enddo
      read (13,*,err=2009) psCorFlag
      write (18,*)
      write (18,*) 'psCorFlag = ', psCorflag
      write (18,*)

c     Read bins for de-aggregating the hazard
      read (13,*,err=2010) nMagBins
      call S21_CheckDim ( nMagBins, MAX_MAG,   'MAX_MAG   ' )
      read (13,*,err=2011) (magBins(k),k=1,nMagBins)
      read (13,*,err=2012) nDistBins
      call S21_CheckDim ( nDistBins, MAX_DIST, 'MAX_DIST  ' )
      read (13,*,err=2013) (distBins(k),k=1,nDistBins)
      read (13,*,err=2014) nEpsBins
      call S21_CheckDim ( nEpsBins, MAX_EPS,   'MAX_EPS   ' )
      read (13,*,err=2015) (epsBins(k),k=1,nEpsBins)

      read (13,*,err=2016) nXcostBins
      call S21_CheckDim ( nXcostBins, MAX_Xcost, 'MAX_Xcost ' )
      read (13,*,err=2017) (XcostBins(k),k=1,nXcostBins)

      read (13,*,err=2018) soilAmpFlag

      call itime(starttime)

      return
 2001 write (*,'( 2x,''input file error: minlat, maxlat line'')')
      stop 99
 2002 write (*,'( 2x,''input file error: maxdist line'')')
      stop 99
 2003 write (*,'( 2x,''input file error: nProb, nAttenType'')')
      stop 99
 2004 write (*,'( 2x,''input file error: spectral period line, iProb='',i5)') i
      stop 99
 2005 write (*,'( 2x,''input file error: Nz & Zvalues line, iProb='', i5)') i
      stop 99
 2006 write (*,'( 2x,''input file error: Natten line for iProb and jAttenType='',2i5)') i, j
      stop 99
 2007 write (*,'( 2x,''input file error: jcalc line for iProb, jAttenType, and iAtten='')')
      write (*,'( 3i5)') i,j,k
      stop 99
 2008 write (*,'( 2x,''input file error: common form coeff for iProb, jAttenType, and iAtten='',3i5)') i,j,k
      stop 99
 2009 write (*,'( 2x,''input file error: pscorflag'')')
      stop 99
 2010 write (*,'( 2x,''input file error: deagg nMagbins'')')
      stop 99
 2011 write (*,'( 2x,''input file error: deagg magbins'')')
      stop 99
 2012 write (*,'( 2x,''input file error: deagg nDistbins'')')
      stop 99
 2013 write (*,'( 2x,''input file error: deagg distbins'')')
      stop 99
 2014 write (*,'( 2x,''input file error: deagg nEpsbins'')')
      stop 99
 2015 write (*,'( 2x,''input file error: deagg epsbins'')')
      stop 99
 2016 write (*,'( 2x,''input file error: deagg nXcosTbins'')')
      stop 99
 2017 write (*,'( 2x,''input file error: deagg xcosT bins'')')
      stop 99
 2018 write (*,'( 2x,''input file error: soil amp flag'')')
      stop 99
 2020 write (*,'( 2x,''input file error: fault file version number'')')
      stop 99
 2100 write (*,'( 2x,''run file does not exist'')')
      write (*,'( 2x,a80)') file1
      stop 99
 2101 write (*,'( 2x,''fault file does not exist'')')
      write (*,'( 2x,a80)') filein
      stop 99

      end

c  --------------------------
      subroutine S21_CheckDim ( n, nMax, name )

      implicit none

      character(len=*) :: name
      integer n, nMax

      if ( n .gt. nMax ) then
        write (*,*)
        write (*,'( 2x,''Array Dimension Too Small'')')
        write (*,'( 2x,''Increase '',a10,''to '',i5)') name, n
        write (*,*)
        stop 99
      endif
      return
      end

c  --------------------------

      subroutine S21_CheckWt ( x, n, fName, name )

      implicit none
      include 'pfrisk.h'

      integer i, n
      real x(*), sum
      character*20 name, fName

      sum = 0.
      do i=1,n
        sum = sum + x(i)
      enddo
      sum = sum - 1.
      if ( abs(sum) .gt. 0.001 ) then
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 2x,a11)') name
        write (*,'( i5, 10f10.4)') n, (x(i),i=1,n)
        write (*,'( 2x,''sum ='',f12.7)') sum
        write (*,'( 2x,a80)') fName
        stop 99
      endif
      return
      end

c  --------------------------

      subroutine S21_CheckWt1 ( x, n, j, n1, fName, name  )

      implicit none

      integer n, j, n1, i
      real x(n1,1), sum
      character*20 fName, name

      sum = 0.
      do i=1,n
        sum = sum + x(j,i)
      enddo
      if ( abs(sum-1.)  .gt. 0.001 ) then
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 10f10.4)') (x(j,i),i=1,n),sum
        write (*,'( 2x,a20)') name
        write (*,'( 2x,a80)') fName
        stop 99
      endif
      return
      end
c  --------------------------

      subroutine S21_CheckWt1a ( x, n, j, n1, fName, name  )

      implicit none
      include 'pfrisk.h'

      integer n, j, n1, i
      real x(MAX_FLT,MAX_N1), sum
      character*20 fName, name

      sum = 0.
      do i=1,n
        sum = sum + x(j,i)
      enddo
      if ( abs(sum-1.)  .gt. 0.001 ) then
        write (*,'( 2x,''Error -- Weights do not sum to unity'')')
        write (*,'( 10f10.4)') (x(j,i),i=1,n),sum
        write (*,'( 2x,a20)') name
        write (*,'( 2x,a80)') fName
        stop 99
      endif
      return
      end

c  -------------------------------------------------------------------
      subroutine S21_WriteTempHaz ( PCflag, tempHaz, tempPC_D, nPC, nParamVar, nInten,
     1                nProb, nAtten, iFlt, jtype, nfType, iFltWidth, nWidth )

      implicit none
      include 'pfrisk.h'

      integer nProb, jType, iProb, iFtype, j, nParamVar(MAX_FLT,MAX_WIDTH),
     1        nInten(MAX_PROB), nAtten(MAX_PROB,4), iflt, nfType(MAX_FLT),
     2        iFltWidth, nWidth(MAX_FLT), iAtten, iParam, nPC, iPC, PCflag(MAX_PROB)
      real*8 tempHaz(MAXPARAM,MAX_INTEN,MAX_PROB,MAX_ATTEN,MAX_FTYPE),
     1       tempPC_D(7,MAXPARAM,MAX_INTEN,MAX_PROB,MAX_FTYPE)

      do iProb=1,nProb

        write (11,'( 20i5)') iFlt, iFltWidth, iProb, nAtten(iProb,jType),
     1     nWidth(iFlt), nFtype(iFlt),nParamVar(iFlt,iFltWidth), nInten(iProb)

        do iAtten = 1,nAtten(iProb,jType)
          do iFtype=1,nFtype(iFlt)
            do iParam=1,nParamVar(iFlt,iFltWidth)
              if (PCflag(iProb) .eq. 0) then
                write (11,'( 20e15.6 )')  (tempHaz(iParam,j,iProb,iAtten,iFtype),j=1,nInten(iProb))
              elseif (PCflag(iProb) .eq. 1) then
                do iPC=1,nPC
                  write (11,'( 20e15.6 )') (tempPC_D(iPC,iParam,j,iProb,iFtype),j=1,nInten(iProb))
                enddo
              endif
            enddo
          enddo
        enddo
      enddo

      return
      end

c  -------------------------------------------------------------------

      subroutine S21_WriteTempHaz1 ( tempHaz1,nParamVar, nInten,
     1                nProb, nAtten, iFlt, jtype, nfType, iFltWidth, nWidth )

      implicit none
      include 'pfrisk.h'

      integer nProb, jType, nwr, iProb, i, iFtype, j, nParamVar(MAX_FLT,MAX_WIDTH),
     1        nInten(MAX_PROB), nAtten(MAX_PROB,4), iFlt, nfType(MAX_FLT),
     2        iFltWidth, nWidth(MAX_FLT), iParam
      real*8 tempHaz1(MAXPARAM,MAX_INTEN, MAX_PROB,MAX_FTYPE)

      nwr = 27

      do iProb=1,nProb

      write (27,'( 20i5)') iFlt, iProb, nAtten(iProb,jType), nWidth(iFlt),
     1         nFtype(iFlt),(nParamVar(iFlt,i),i=1,nWidth(iFlt)),
     2         nInten(iProb)

          do iFtype=1,nFtype(iFlt)
           do iParam=1,nParamVar(iFlt,iFltWidth)
             write (27, '( 15e12.4)')  (tempHaz1(iParam,j,iProb,iFtype),j=1,nInten(iProb))
           enddo
          enddo
      enddo

      return
      end

c  -------------------------------------------------------------------

      subroutine S21_WriteTempHaz2 ( tempHaz2, nInten, nProb, nAtten, nattenType )

      implicit none
      include 'pfrisk.h'

      integer nProb, nwr, iProb, jType, j, nInten(MAX_PROB),
     1        nAtten(MAX_PROB,4), nAttenType, iAtten
      real*8 tempHaz2(4, MAX_INTEN, MAX_PROB,MAX_ATTEN)

      nwr = 28

      do iProb=1,nProb
        do jType=1,nattenType
           do iAtten = 1,nAtten(iProb,jType)
               write (28,'( 3i5, 100e12.4 )')  iProb, jType, iAtten,
     1             (tempHaz2(jType,j,iProb,iAtten),j=1,nInten(iProb))
           enddo
         enddo
      enddo

      return
      end

c  -------------------------------------------------------------------

      subroutine S21_output_TotalHaz ( isite, sitex, sitey, ati, nInten,
     1           nFlt, nAtten, risk, fName, jCalc, sigTrunc, csrflag,
     2           attenName, period1, probAct, nWidth,
     4           m_bar, d_bar, e_bar, riskBins, nMagBins, nDistBins,
     5           nEpsBins, magBins, distBins, epsBins,
     6           al_segWt, MinRrup, nAttenType, attenType, pdfsum, segWt1,
     7           dirflag,tapflag, intflag, fsys, SourceDist, mMagout,
     8           hwflagout, ftype, vs, nmaxmag2, mMagoutwt, specT)

      implicit none
      include 'pfrisk.h'

      real siteX, siteY, ati(MAX_PROB, MAX_INTEN), pdfsum(MAX_FLT)
      real*8 risk(MAX_INTEN,MAX_PROB,MAX_FLT)
      real al_segWt(MAX_FLT), segWt1(MAX_FLT), magbins(MAX_MAG), distbins(MAX_DIST),
     1     epsBins(MAX_EPS)
      real*8 m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN),
     3     e_bar(MAX_PROB,MAX_INTEN)
      real riskBins(MAX_MAG,MAX_DIST,MAX_EPS,MAX_PROB,MAX_INTEN),
     5     outrisk(MAX_MAG,MAX_DIST,MAX_PROB,MAX_INTEN),
     9     sumbar, xx
      real mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real MinRrup(MAX_FLT), mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real*8 risk2(MAX_INTEN), sum, risk2t(MAX_INTEN)
      real  period1(4,MAX_PROB),
     1     sigTrunc(MAX_PROB), probAct(MAX_FLT), SourceDist(MAX_FLT,MAX_WIDTH,3)
      integer csrflag(MAX_PROB), attenType(MAX_FLT), nAttenType, dirflag(MAX_PROB)
      integer tapflag(MAX_PROB), intflag, nAtten, nwr, ii, k, l
      integer nMagBins, nDistBins, nEpsBins, fsys(MAX_FLT), kk
      integer isite, nInten(MAX_PROB), nFlt, jCalc(MAX_PROB,4,MAX_ATTEN), nWidth(MAX_FLT)
      integer hwflagout(MAX_FLT), nMaxmag2(MAX_FLT), casecount, iProb
      integer J2, iInten, iFlt, M, iMagBin, iDistBin, iEpsBin
      real ftype(1), vs, specT(MAX_PROB)
      character*80 fName(MAX_FLT), file1, attenName

      nwr = 12
c     Loop over different maximum magnitude values to get
c     a total number of fault listings in the output file.
      casecount = 0
      do ii=1,nFlt
         do k=1,nWidth(ii)
            do l=1,nMaxMag2(ii)
               casecount = casecount + 1
            enddo
         enddo
      enddo

      nwr = 12
      write (nwr,'(2i8,10x,''Number of faults, fault cases'')')
     1       nFlt,casecount

      do ii=1,nFlt
         do k=1,nWidth(ii)
            do l=1,nMaxMag2(ii)

               write (nwr,1111) fname(ii), ii, fsys(ii),
     1         (SourceDist(ii,k,kk),kk=1,3), k,l, mMagout(ii,k,l),
     2         mMagoutWt(ii,k,l), ftype(ii), hwflagout(ii), vs

            enddo
         enddo
      enddo

 1111 format (a30,i6,4x,i5,7x,f8.1,7x,f6.1,4x,f8.1,4x,2i6,6x,f4.2,
     1        f7.3,5x,f4.1,2x,i8,3x,f8.1)

      write (18,*)
      write (18,'(1x,a34,2x,a80)') 'Output filename for hazard curves:', file1
      write (18,*)
      write (18,*) 'Summary of Faults used in the analysis:'
      write (18,*)
      write (18,'(i8,10x,''Number of faults'')') nFlt
      write (18,*)
      write (18,'(a50)') 'Fault Name                     nFlt  nSys  p1sum '
      write (18,'(1x,a50)') '--------------------------------------------------'
      do ii=1,nFlt
         write (18,'(a30,i5,i5,4x,f7.5)') fname(ii), ii, fsys(ii), pdfsum(ii)
      enddo

c     Write Site Coordinates
      write(NWR,910) iSite, siteX, siteY
      write (nwr,'( i5, 2x,''nAtten'')') nAtten

c     Loop over each attenuation model
      do 900 iProb=1,nAtten

c       Label Attenuation
        write (nwr,'( 2x,''Attenuation: '',i5,2x,f10.4)') iProb, specT(iProb)

c       Write Test Intensities
        write (nwr,'( i5,2x,''nAmp'')') nInten(iProb)
        write(NWR,920) (ati(iProb,J2),J2=1,nInten(iProb))

c       Add up Risks (number of events)
        do iInten=1,nInten(iProb)
          sum = 0.
          do iFlt=1,nFlt
            sum = sum + risk(iInten,iProb,iFlt)
          enddo
          risk2( iInten) = sum
        enddo

c       Write Results
        do iFlt=1,nFlt

          write (nwr,'( 2x,a38,2f6.3,f8.1,x,30e12.4)') fName(iFlt),
     1          segwt1(iFlt), al_segWt(iFlt),
     2          MinRrup(iFLt),
     3          (risk(iInten,iProb,iFlt),iInten=1,nInten(iProb))
        enddo

c       Write Total Number of Events
        xx = 0.
        write(NWR,'( 2x,''Wt_Total_Events/yr'',20x,2f6.3,f8.1,1x,30e12.4)' )
     1        xx, xx, xx, (risk2( M ), M=1,nInten(iProb))

c       Compute probabilitiy of exceedane (Poisson)
        do iInten=1,nInten(iProb)
          risk2t(iInten) = risk2(iInten)
          if (risk2(iInten) .gt. 0.0001) then
            risk2(iInten) = (1. - exp(-risk2(iInten)))
          endif
        enddo

c       Write hazard
        write(NWR,'( 2x,''Poisson_Prob:  '',23x,2f6.3,f8.1,1x,30e12.4)')
     1        xx, xx, xx, ( risk2(M), M=1,nInten(iProb))

c       Write out M,D and epsilon Bar values.
c       Normalize for each intensity
        do iInten=1,nInten(iProb)
          sumbar = 0.
          do iMagBin=1,nMagBins
            do iDistBin=1,nDistBins
              do iEpsBin=1,nEpsBins
                sumbar = sumbar + riskBins(iMagBin,iDistBin,iEpsBin,
     1                   iProb,iInten)
              enddo
            enddo
          enddo
          if (sumbar.gt.0.0) then
            do iMagBin=1,nMagBins
              do iDistBin=1,nDistBins
                do iEpsBin=1,nEpsBins
                  riskBins(iMagBin,iDistBin,iEpsBin,
     1                       iProb,iInten) =
     1            riskBins(iMagBin,iDistBin,iEpsBin,
     2                       iProb,iInten)/sumbar
                  outrisk(iMagBin,iDistBin,iProb,iInten) =
     1              outrisk(iMagBin,iDistBin,
     1                       iProb,iInten) +
     2              riskBins(iMagBin,iDistBin,
     2                       iEpsBin,iProb,iInten)
                enddo
              enddo
            enddo
            m_bar(iProb,iInten) = m_bar(iProb,iInten)/risk2t(iInten)
            d_bar(iProb,iInten) = d_bar(iProb,iInten)/risk2t(iInten)
            e_bar(iProb,iInten) = e_bar(iProb,iInten)/risk2t(iInten)
          endif
        enddo

c       Write MBar, DBar and EpsBar Results

        write (nwr,'(2x,''M_bar'',33x,2f6.3,f8.1,2x,30(e10.3,2x))')  xx, xx, xx,
     1  (m_bar(iProb,iInten),iInten=1,nInten(iProb))
            write (nwr,'(2x,''D_bar'',33x,2f6.3,f8.1,2x,30(e10.3,2x))') xx, xx, xx,
     1  (d_bar(iProb,iInten),iInten=1,nInten(iProb))
            write (nwr,'(2x,''Eps_bar'',31x,2f6.3,f8.1,2x,30(e10.3,2x))') xx, xx, xx,
     1  (e_bar(iProb,iInten),iInten=1,nInten(iProb))
        write (nwr,'( /,''------------------'',/)')

 900  continue
      close (nwr)
      return

 910  FORMAT(//,115('-'),//' SITE',I2,' COORDINATES: ',2F9.3)
 920  FORMAT(1X,'AMP:',36x,' pSeg al_Wt  MinRrup',30F12.6)
      end

c  ------------------------------------------------------------------

      subroutine S21_output_hazBins ( isite, sitex, sitey, testInten, nInten,
     1           nProb, riskBins, nMagBins, nDistBins, nEpsBins, magBins,
     2           distBins, epsBins, m_bar, d_bar, e_bar, Xcost_bar,
     3           nXcostBins, XcostBins, RiskBinsX )

      implicit none
      include 'pfrisk.h'

      integer nProb, i, isite, nInten(MAX_PROB), nMagBins, nDistBins,
     1        nEpsBins, nXcostBins, J2, iEpsBin, iDistBin, iMagBin,
     2        iXcostBin, iInten, iMagbins, iDistbins, iProb, nwr
      real siteX, siteY, testInten(MAX_PROB,MAX_INTEN), distBins(MAX_DIST),
     1     riskBins(MAX_MAG,MAX_DIST,MAX_EPS,MAX_PROB,MAX_INTEN),
     2     epsBins(MAX_EPS), XcostBins(MAX_Xcost), magBins(MAX_MAG),
     3     outrisk(MAX_MAG,MAX_DIST,MAX_PROB,MAX_INTEN),
     4     RiskBinsX(MAX_Xcost,MAX_PROB,MAX_INTEN)
      real*8 sumX, sum, m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN),
     1     e_bar(MAX_PROB,MAX_INTEN), Xcost_bar(MAX_PROB,MAX_INTEN)
      character*80 file1

      write (18,'(1x,a50)') '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
      write (18,*)
      write (18,'(1x,a34,2xa80)') 'Filename for deaggregation output:', file1
      write (18,*)

      nwr = 14

c     Write number of bins and bin values.
      write (nwr,'(i5,5x,a25)') nMagbins,  'Number of Magnitude Bins'
      write (nwr,'(20f10.3)') (magBins(i),i=1,nMagbins)
      write (nwr,'(i5,5x,a25)') nDistbins, 'Number of Distance Bins '
      write (nwr,'(20f10.3)') (distBins(i),i=1,nDistbins)
      write (nwr,'(i5,5x,a25)') nEpsbins,  'Number of Epsilon Bins  '
      write (nwr,'(20f10.3)') (epsBins(i),i=1,nepsbins)
      write (nwr,'(i5,5x,a25)') nXcostbins,'Number of Xcos(T) Bins  '
      write (nwr,'(20f10.3)') (xCostBins(i),i=1,nxCostbins)

      write (18,'(i5,5x,a25)') nMagbins,  'Number of Magnitude Bins'
      write (18,'(20f10.3)') (magBins(i),i=1,nMagbins)
      write (18,'(i5,5x,a25)') nDistbins, 'Number of Distance Bins '
      write (18,'(20f10.3)') (distBins(i),i=1,nDistbins)
      write (18,'(i5,5x,a25)') nEpsbins,  'Number of Epsilon Bins  '
      write (18,'(20f10.3)') (epsBins(i),i=1,nepsbins)
      write (18,'(i5,5x,a25)') nXcostbins,'Number of Xcos(T) Bins  '
      write (18,'(20f10.3)') (xCostBins(i),i=1,nxCostbins)

      write (18,*)
      write (18,'(1x,a50)') '--------------------------------------------------'
      write (18,*)

c     Write Site Coordinates
      write(NWR,910) iSite, siteX, siteY
      write (nwr,'( i5, 2x,''nProb'')') nProb

C     Initialize Outrisk Bins for each site location.
      do iMagbins=1,MAX_MAG
         do iDistbins=1,MAX_DIST
            do iProb=1,MAX_PROB
               do iInten=1,MAX_INTEN
                  outrisk(iMagBins, iDistbins, iProb, iInten) = 0.0
               enddo
            enddo
         enddo
      enddo

c     Loop over each attenuation relationship.

      do 900 iProb=1,nProb

c       Write Test Intensities
        write (nwr,'( 2i5,''nAmp'')') nInten(iProb)
        write(NWR,920) (testInten(iProb,J2),J2=1,nInten(iProb))

c       Normalize for each intensity
        do iInten=1,nInten(iProb)
          sum = 0.
          sumX = 0.
          do iXcostBin=1,nXcostBins
             sumX = sumX + RiskBinsX(iXcostBin,iProb,iInten)
          enddo

          if (sumX .gt. 0.0) then
             do iXcostBin=1,nXcostBins
                 RiskBinsX(iXcostBin,iProb,iInten) =
     1               RiskBinsX(iXcostBin,iProb,iInten)/sumX
             enddo
          endif

          do iMagBin=1,nMagBins
            do iDistBin=1,nDistBins
              do iEpsBin=1,nEpsBins
                sum = sum + riskBins(iMagBin,iDistBin,iEpsBin,
     1               iProb,iInten)
               enddo
             enddo
           enddo

           if (sum.gt.0.0) then
            do iMagBin=1,nMagBins
              do iDistBin=1,nDistBins
                do iEpsBin=1,nEpsBins
                    riskBins(iMagBin,iDistBin,iEpsBin,
     1                       iProb,iInten) =
     1              riskBins(iMagBin,iDistBin,iEpsBin,
     2                       iProb,iInten)/sum

                    outrisk(iMagBin,iDistBin,iProb,iInten) =
     1                       outrisk(iMagBin,iDistBin,
     1                       iProb,iInten) +
     2                       riskBins(iMagBin,iDistBin,
     2                       iEpsBin,iProb,iInten)
                enddo
              enddo
            enddo
            m_bar(iProb,iInten) = m_bar(iProb,iInten)/sum
            d_bar(iProb,iInten) = d_bar(iProb,iInten)/sum
            e_bar(iProb,iInten) = e_bar(iProb,iInten)/sum
            Xcost_bar(iProb,iInten) = Xcost_bar(iProb,iInten)/sumX

          endif
        enddo

c       Write Results

        write (nwr,'( /,''------------------'',/)')

        write (nwr,'(2x,''M_bar'',39x,30(e10.3,2x))')
     1  (m_bar(iProb,iInten),iInten=1,nInten(iProb))
        write (nwr,'(2x,''D_bar'',39x,30(e10.3,2x))')
     1  (d_bar(iProb,iInten),iInten=1,nInten(iProb))
           write (nwr,'(2x,''Eps_bar'',37x,30(e10.3,2x))')
     1  (e_bar(iProb,iInten),iInten=1,nInten(iProb))

           write (nwr,'(2x,''XCosT_bar'',35x,30(e10.3,2x))')
     1  (Xcost_bar(iProb,iInten),iInten=1,nInten(iProb))
        write (nwr,'( /,''------------------'',/)')

        write (nwr,'('' XCos(Theta) Bins'')')
        do iXcostBin=1,nXcostBins-1
           write (nwr,'(2f8.4,2x,30(e10.3,2x))') XcostBins(iXcostBin),
     1           XcostBins(iXcostBin+1),
     1      (RiskBinsX(iXcostBin,iProb,iInten),iInten=1,nInten(iProb))
        enddo

        write (nwr,'( /,''------------------'',/)')
        write (nwr,'(/,2x,3a14,/)') 'Eps. Range ', 'Mag. Range',
     1  'Dist. Range'

        do iMagBin=1,nMagBins-1
          do iDistBin=1,nDistBins-1
            write (nwr,'(2x,2f7.3,2f7.2,2f7.1,2x,30(e10.3,2x))') epsBins(1),
     1             epsBins(nEpsBins),
     1             magBins(iMagBin),
     1             magBins(iMagBin+1),
     1             distBins(iDistBin), distBins(iDistBin+1),
     2             (outrisk(iMagBin, iDistBin, iProb, iInten),
     2             iInten=1,nInten(iProb))
          enddo
        enddo

        write (nwr,'( /,''------------------'',/)')
        write (nwr,'(/,2x,3a14,/)') 'Eps. Range ', 'Mag. Range',
     1          'Dist. Range'

        do iEpsBin=1,nEpsBins-1
         do iMagBin=1,nMagBins-1
          do iDistBin=1,nDistBins-1
           write (nwr,'(2x,2f7.3,2f7.2,2f7.1,2x,30(e10.3,2x))')epsBins(iEpsBin),
     1             epsBins(iEpsBin+1),
     1             magBins(iMagBin),
     1             magBins(iMagBin+1),
     1             distBins(iDistBin), distBins(iDistBin+1),
     2             (riskBins(iMagBin,iDistBin,iEpsBin,iProb,iInten),
     2             iInten=1,nInten(iProb))
          enddo
         enddo
        enddo

        write (nwr,'( /,''-------------------------------------
     1--------------------------------------'',/)')

 900  continue
      close (nwr)
      return
 910  FORMAT(//,115('-'),//' SITE',I2,' COORDINATES: ',2F9.3)
 920  FORMAT(1X,'AMP:',41x,30(F10.3,2x))
      end

c ----------------------------------------------------------------------
      real function f_magModel(c,x, i)

      implicit none

      real c(4), x
      integer i

      if ( x .lt. c(3) .or. c(3) .eq. 0.) then
        f_magModel = c(1) + c(2)*alog10(x)
      else
        f_magModel = c(1) + c(2)*alog10(c(3)) + c(4)*(alog10(x)-alog10(c(3)))
      endif
      return
      end

c  -------------------------------------------------------------------
      subroutine S21_writedisthead (ioflag)

      implicit none

      integer ioflag
      if (ioflag .eq. 1) then
         write (19,*) 'Column Heading Information:   '
         write (19,*) 'HW = Hanging Wall (1) or Footwall (0) flag'
         write (19,*) 'S2Site = Azimuth between closest point and station location'
         write (19,*) 'Flen = Fault length in km'
         write (19,*) 'FWid = Fault width in km'
         write (19,*) 'Strike = Strike of fault plane'
         write (19,*) 'Dip = Dip of fault plane'
         write (19,*) 'Step = Step size in km for fault plane'
         write (19,*) 'cDD = Cell number downdip for closest point'
         write (19,*) 'nDD = Number of downdip cells for fault plane'
         write (19,*) 'cAS = Cell number along strike for closest point'
         write (19,*) 'nAS = Number of along strike cells for fault plane'
         write (19,*) 'cFlen = Distance along strike from Hypocenter to closest point'
         write (19,*) 'hDDkm = Hypocenter location (km) downdip'
         write (19,*) 'hASkm = Hypocenter location (km) along strike'
         write (19,*) 'hDD = Hypocenter cell location downdip'
         write (19,*) 'hAS = Hypocenter location along strike'
         write (19,*) 'eDist = Epicentral distance'
         write (19,*) 'hDist = Hypocentral distance'
         write (19,*) 'sLit = Distance along strike from Hypocenter to closest point'
         write (19,*) 'AzP1P2 = Angle between hypocenter and station location'
         write (19,*) 'X = Ratio of sLit/Fault length'
         write (19,*) 'dLit = Distance down dip from Hypocenter to closest point'
         write (19,*) 'PhiAng = Angle between hypocenter and station location'
         write (19,*) 'Y = Ratio of dLit/Fault width'
         write (19,*) 'Hypo(x,y,z) = Hypocenter location in km with station location at (0,0,0)'
         write (19,*) 'Closest(x,y,z) = Closest point location in km with station location at (0,0,0)'
         write (19,*)
      elseif (ioflag .eq. 2) then
         write (19,*)
         write (19,'(4a61,a42)') 'Site   Fault Name                    RupDist  JBDist SeismoD ',
     1                           '    Rx    HW  S2Site    FLen    FWid  Strike    Dip     Step ',
     2                           ' cDD  nDD  cAS  nAS   cFLen   hDDkm   hASkm  hDD  hAS   eDist',
     3                           '   hDist    sLit  AzP1P2     X      dlit  phiang     Y       ',
     4                           '   Hypo(x,y,z)              Closest(x,y,z)'
      endif

      return
      end
c  -------------------------------------------------------------------

      subroutine S21_output_Sourcedeagg ( isite, sitex, sitey, testInten, nInten,
     1           nFlt, haz, fName, m_bar_s, rrup_bar_s, rjb_bar_s,
     2           rx_bar_s, e_bar_s, rSeismo_bar_s, ry0_bar_s, mag_bar_s,
     3           ftype_bar_s, hypodepth_bar_s, dipavgd_bar_s, ztor_bar_s,
     4           thetasite_bar_s, rupwidth_bar_s, rHypo_bar_s, specT, nProb)

      implicit none
      include 'pfrisk.h'

      integer isite, nInten(MAX_PROB), iFlt, nFlt, iProb, J2, nProb, iInten, nwr
      real siteX, siteY, testInten(MAX_PROB, MAX_INTEN), specT(MAX_PROB)
      real*8 haz(MAX_INTEN,MAX_PROB,MAX_FLT), m_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN),
     1       rrup_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN), rjb_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN),
     2       rx_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN), e_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN),
     3       rSeismo_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), ry0_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     4       mag_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), ftype_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     5       hypodepth_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), dipavgd_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     6       ztor_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), thetasite_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     7       rupwidth_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN), rHypo_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN)
      character*80 fName(MAX_FLT)

      nwr = 29

c     Write Site Coordinates
      write(nwr,910) iSite, siteX, siteY


        do iProb=1,nProb
c       Write Test Intensities
         write (nwr,'( 2x,''iProb, specT: '',i5,2x,f10.4)') iProb, specT(iProb)
         write (nwr,'( i5,2x,''nAmp'')') nInten(iProb)
         write(nwr,920) (testInten(iProb,J2),J2=1,nInten(iProb))

c        Write MBar, RrupBar, Rjbbar, Rxbar and EpsBar Results
         do iFlt=1,nFlt
          write (nwr,'( 2x,a30,''MBar      '',30f10.2)') fName(iFlt),
     1          (m_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''RrupBar   '',30f10.2)') fName(iFlt),
     1          (rrup_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''RjbBar    '',30f10.2)') fName(iFlt),
     1          (rjb_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''RxBar     '',30f10.2)') fName(iFlt),
     1          (rx_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''EBar      '',30f10.2)') fName(iFlt),
     1          (e_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))

          write (nwr,'( 2x,a30,''rSeismoBar'',30f10.2)') fName(iFlt),
     1          (rSeismo_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''ry0Bar    '',30f10.2)') fName(iFlt),
     1          (ry0_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''rHypoBar  '',30f10.2)') fName(iFlt),
     1          (rHypo_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''DipAvgBar '',30f10.2)') fName(iFlt),
     1          (dipavgd_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''ZtorBar   '',30f10.2)') fName(iFlt),
     1          (Ztor_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''thetasiteBar '',1f7.2, 29f10.2)') fName(iFlt),
     1          (thetasite_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''RupWidthBar '',1f8.2,29f10.2)') fName(iFlt),
     1          (rupwidth_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))
          write (nwr,'( 2x,a30,''HypoDepthBar '',1f7.2,29f10.2)') fName(iFlt),
     1          (hypodepth_bar_s(iFlt,iProb,iInten)/haz(iInten,iProb,iFlt),iInten=1,nInten(iProb))


         enddo
        write (nwr,'( /,''------------------'',/)')
        enddo


 900  continue
      close (nwr)
      return

 910  FORMAT(//,115('-'),//' SITE',I2,' COORDINATES: ',2F9.3)
 920  FORMAT(1X,'AMP:',39x,30(F8.3,2x))
      end

c  -------------------------------------------------------------------

      subroutine S21_CheckPC (PCflag)

      implicit none

      integer PCflag

        if (PCflag .eq. 0) then
          goto 100
        elseif (PCflag .eq. 1) then
          goto 100
        else
          write (*,*) 'Invalid PC flag in run file.'
          stop 99
        endif

 100  continue

      return
      end
