
      subroutine RdInput ( nProb, nAttenType, nAtten, jcalc, specT, sigTrunc,
     1               gmScale, dirFlag, nInten, testInten, lgTestInten, 
     2               psCorFlag, minlat, maxlat, minlong, maxlong, distmax,
     3               nMagBins, magBins, nDistBins, distBins, nepsBins, epsBins,
     4               nXcostBins, xcostBins, soilAmpFlag, gm_wt, runflag, sigvaradd,
     5               sCalc, sigfix, sssCalc, bnumflag, cfcoefrrup, cfcoefrjb, coefcountrrup, coefcountRjb,
     6               iMixture )

      include 'pfrisk.h'

      real testInten(MAX_PROB, MAX_INTEN), lgTestInten(MAX_PROB, MAX_INTEN)
      real magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS),  XcostBins(MAX_Xcost)
      real sigTrunc(MAX_PROB) 
      integer attenType(MAX_FLT), nInten(MAX_PROB), runFlag, bnumflag
      integer nMagBins, nDistBins, nEpsBins, nXcostBins
      integer psCorFlag,dirflag(MAX_PROB)
      integer jCalc(MAX_PROB,4,MAX_ATTEN)
      integer sCalc(MAX_PROB,4,MAX_ATTEN),sssCalc(MAX_PROB,4,MAX_ATTEN)
      real sigfix(MAX_PROB,4,MAX_ATTEN)
      character*80 filein, title, file1, filelog, filehyp
      real minlat, maxlat, minlong, maxlong, vs
      integer nProb, nAttenType, nAtten(MAX_PROB,4), soilAmpFlag
      real specT(MAX_PROB)
      real gmScale(MAX_PROB,4,MAX_ATTEN), gm_wt(MAX_PROB,4,MAX_ATTEN), 
     1     sigvaradd(MAX_PROB,4,MAX_ATTEN)
      integer ichart, charcount, coefcountrrup, coefcountrjb, nDD(MAX_FLT)
      integer iMixture(4, MAX_PROB, MAX_ATTEN)
      real cfcoefrrup(MAX_Atten,11), cfcoefrjb(MAX_Atten,11)
      real syn_dip(MAX_FLT,MAX_SYN), syn_zTOR(MAX_FLT,MAX_SYN)
      real syn_RupWidth(MAX_FLT,MAX_SYN), syn_RX(MAX_FLT,MAX_SYN), syn_Ry0(MAX_FLT,MAX_SYN)

c     Set Data file units
      nwr = 11

C     Run Ground Motion Case or PSHA case.
 123  write (*,*) 'Enter 0 to run hazard'
      write (*,*) 'Enter 1 to run attenuation models'
      write (*,*) 'Enter 2 to run spectra models'
      write (*,*) 'Enter 3 for Fault File testing'
      write (*,*) 'Enter 4 for Distance calculation testing'

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
c     Create log filename for PHSA run.
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
         call determ (runFlag)
         stop
      elseif (runFlag .eq. 2) then
         write (*,*) 'Enter the input file name for spectra runs.'
         read (*,'( A80)') file1
         open (13,file=file1,status='old')
         call determ (runFlag)
         stop
      elseif (runFlag .eq. 3) then
         write (*,*) 'Enter the input fault file name for Fault File testing.'
         read (*,'( A80)') file1
         open (10,file=file1,status='old')
         write (*,*) 'Enter the ouput file name for Fault File testing.'
         read (*,'( A80)') filelog
         open (17,file=filelog,status='new')        
         write (17,*) 'Fault File Testing for fault file: ', file1
         write (17,*) 

c     read fault File
      call Rd_Fault_Data ( nFlt, fName, minMag, magStep, xStep,
     1     yStep, segModelWt, rateParam, rateParamWt, beta,
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, 
     3     maxMag,  maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, 
     4     ftype, sourceType, nRupArea, coeff_area, sigArea, nRupWidth, 
     5     coeff_width, sigWidth, nParamVar, iCoor,minDepth,
     6     fIndex, probAct, nWidth, mpdf_param, 
     7     al_segWt, attenType, sampleStep,
     8     grid_a,grid_dlong,grid_dlat,grid_n,grid_long, grid_lat,
     9     grid_top, minlat, maxlat, minlong, maxlong, scaleRate, fsys,
     1     mMagout, mMagoutWt, fltDirect, synchron, nsyn_Case, synatten,
     1     synmag, syndistRup, syndistJB, synDistSeismo, synHypo,
     2     synftype, synhwflag, synwt, RateType, iDepthModel, depthParam, 
     3     nMaxmag2, segWt1, faultFlag, nDD, nFtype, ftype_wt, 
     4     br_index, br_wt, segModelFlag, nSegModel, segModelWt1, runflag,
     7     syn_dip, syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0 )

           write (*,*)
           write (*,*) 'Normal termination of fault file for testing.'
           write (17,*)
           write (17,*) '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
           write (17,*)
           write (17,*) 'Normal termination of fault file for testing.'

           close (10)
           close (17)

         stop
      elseif (runFlag .eq. 4) then
         write (*,*) 'Enter the input file name for distance checking.'
         read (*,'( A80)') file1
         open (13,file=file1,status='old')
         write (*,*) 'Enter the output file name for distance checking.'
         read (*,'( A80)') filelog
         open (19,file=filelog,status='new')
         write (*,*) 'Enter the hypocenter location file name for distance checking.'
         read (*,'( A80)') filehyp
         open (29,file=filehyp,status='old')
c         call determ (runFlag)
c         stop
      else
         write (*,*) 'Invalid entry!!!'
         write (*,*) 'Try again.'
         goto 123
      endif
 
c     Read fault file name      
      read (13,'( a80)') filein
      open (10,file=filein,status='old',err=2101)

      if (runflag .eq. 4) then
         write (19,*) 'Distance Checking Results:'
         write (19,'(a28,a80)') 'Input Run filename       :  ', file1
         write (19,'(a28,a80)') 'Output distance filename :  ', filelog
         write (19,'(a28,a80)') 'Fault File filename      :  ', filein
         write (19,'(a28,a80)') 'Hypocenter File filename :  ', filehyp
         write (19,*) 
      endif

c     Read min max long and lat for grid sources
      read (13,*,err=2001) minlat, maxlat, minlong, maxlong
      read (13,*,err=2002) distmax

c     Input options 
      read(13,'( a80)') title
      read(13,*,err=2003) nProb, nAttenType
      
      call CheckDim ( nProb, MAX_Prob, 'MAX_PROB ' )        
      call CheckDim ( nProb, MAX_Prob, 'MAX_ATTEN ' )
      write (18,'( a80)') title
      write (18,*) 
      write (18,*) 'Summary of Attenuation relationship used in the analysis:'
      write (18,*) 
      write (18,'(a20,2x,i4)') 'Number of Problems: ', nProb
      write (18,'(a29,2x,i4)') 'Number of Attenuation Types: ', nAttentype
      write (18,*)

c     Enter GM models for each problem
      do i=1,nProb
        read (13,*,err=2004) specT(i), sigTrunc(i), dirflag(i)
        read (13,*,err=2005) nInten(i), (testInten(i,j),j=1,nInten(i))
        call CheckDim ( nInten(i), MAX_INTEN, 'MAX_INTEN ' )
        do j1=1,nInten(i)
           lgTestInten(i,j1) = alog(testInten(i,j1))
        enddo

        write (18,'(a19,2x,f8.4)') 'Spectral Period  = ', specT(i)
        write (18,'(a19,2x,f8.4)') 'Sigma Trunc      = ', sigTrunc(i)
        write (18,'(a19,2x,i8)')   'Directivity Flag = ', dirflag(i)
        write (18,'(a21)')         'Ground Motion Values:'
        do j1=1,nInten(i)
           write (18,'(f12.4)') testInten(i,j1)
        enddo
        do j=1,nAttenType
           write (18,*) 
           write (18,'(a19,2x,i4)') 'Attenuation Type = ', j
           write (18,*) ' Jcalc    Const1   Const2   Weight   SigVaradd  Scalc    SigFix   SSSCalc    iMix'

          checkattenwt = 0.0
          read (13,*,err=2006) nAtten(i,j)
          call CheckDim ( nAtten(i,j), MAX_ATTEN, 'MAX_ATTEN ' )
          do k=1,nAtten(i,j)
            read (13,*,err=2007) jCalc(i,j,k), gm1,gm2, gm_wt(i,j,k), sigvaradd(i,j,k), iMixture(j, i, k)  

C     Check for either fixed sigma or different sigma model (i.e., jCalc<0)
            if (jcalc(i,j,k) .lt. 0) then
               backspace (13)
               read (13,*,err=2006) jCalc(i,j,k), gm1,gm2, gm_wt(i,j,k), 
     1              sigvaradd(i,j,k), iMixture(j, i, k), sCalc(i,j,k), sigfix(i,j,k)
               write (18,'(i6,3x,f8.4,1x,f8.4,f9.5,f11.4,2x,i6,2x,f8.4,2x,2i6)')
     1              jCalc(i,j,k), gm1,gm2, gm_wt(i,j,k), 
     1             sigvaradd(i,j,k), sCalc(i,j,k), sigfix(i,j,k), ssscalc(i,j,k), iMixture(j, i, k)  
            else
            write (18,'(i6,3x,f8.4,1x,f8.4,f9.5,f11.4,31x,i6)') jCalc(i,j,k), gm1,gm2, gm_wt(i,j,k), 
     1           sigvaradd(i,j,k), iMixture(j, i, k)
            endif

C     Check for Common Functional Form with Rrup Distance (10000<jcalc<11000) selected and if so read in coefficients.
            if (abs(jcalc(i,j,k)) .gt. 10000 .and. abs(jcalc(i,j,k)) .lt. 11000) then
                coefcountrrup = abs(jcalc(i,j,k)) - 10000
                read (13,*,err=2008) (cfcoefrrup(coefcountrrup,jj),jj=1,11)
                write (18,'(a16,2xi8,10(2x,f12.6))') 
     1             'CF(RRup) Coefs: ', coefcountrrup,(cfcoefrrup(coefcountrrup,jj),jj=1,11)
            endif
C     Check for Common Functional Form with RJB Distance (11000<jcalc<12000) selected and if so read in coefficients.
            if (abs(jcalc(i,j,k)) .gt. 11000 .and. abs(jcalc(i,j,k)) .lt. 12000) then
                coefcountrjb = abs(jcalc(i,j,k)) - 11000
                read (13,*,err=2008) (cfcoefrjb(coefcountrjb,jj),jj=1,11)
                write (18,'(a16,2x,i8,10(2x,f12.6))') 
     1              'CF(RJB) Coefs:  ', coefcountrjb, (cfcoefrrup(coefcountrjb,jj),jj=1,11)
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
      call CheckDim ( nMagBins, MAX_MAG,   'MAX_MAG   ' )
      read (13,*,err=2011) (magBins(k),k=1,nMagBins)
      read (13,*,err=2012) nDistBins
      call CheckDim ( nDistBins, MAX_DIST, 'MAX_DIST  ' )
      read (13,*,err=2013) (distBins(k),k=1,nDistBins)
      read (13,*,err=2014) nEpsBins
      call CheckDim ( nEpsBins, MAX_EPS,   'MAX_EPS   ' )
      read (13,*,err=2015) (epsBins(k),k=1,nEpsBins)

      read (13,*,err=2016) nXcostBins
      call CheckDim ( nXcostBins, MAX_Xcost, 'MAX_Xcost ' )
      read (13,*,err=2017) (XcostBins(k),k=1,nXcostBins)

      read (13,*,err=2018) soilAmpFlag

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
 2006 write (*,'( 2x,''input file error: Natten line for iProb and jATtenType='',2i5)') i, j 
      stop 99
 2007 write (*,'( 2x,''input file error: jclac line for iProb, jATtenType, and iAtten='')') 
      write (*,'( 3i5)') i,j,k
      pause 'test 2'
      stop 99
 2008 write (*,'( 2x,''input file error: common for coeff for iProb, jATtenType, and iAtten='',3i5)') i,j,k 
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
 2100 write (*,'( 2x,''run file does not exist'')')
      write (*,'( 2x,a80)') file1
      stop 99
 2101 write (*,'( 2x,''fault file does not exist'')')
      write (*,'( 2x,a80)') filein
      stop 99

      end



c  --------------------------
      subroutine CheckDim ( n, nMax, name )
      character*10 name
      integer n, nMax
      
      if ( n .gt. nMax ) then
        write (*,'( 2x,''Array Dimension Too Small'')')
        write (*,'( 2x,''Increase '',a10,''to '',i5)') name, n
        pause
        stop 99
      endif
      return
      end

c  --------------------------

      subroutine CheckWt ( x, n, fName, name )
      real x(1)
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

      subroutine CheckWt1 ( x, n, j, n1, fName, name  )
      real x(n1,1)
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

      subroutine CheckWt1a ( x, n, j, n1, fName, name  )
      include 'pfrisk.h'
      real x(MAX_FLT,MAX_N1)
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
      subroutine WriteTempHaz ( tempHaz,nParamVar, nInten, 
     1                nProb, nAtten, iFlt, jtype, nftype, iFltWidth, nWidth )
 
      include 'pfrisk.h'
      include 'declare1.h'

      do iProb=1,nProb

      if ( iFltWidth .eq. 1 ) write (11,'( 20i5)') iFlt, iProb, nAtten(iProb,jType),
     1         nWidth(iFlt), nFtype(iFlt),(nParamVar(iFlt,i),i=1,nWidth(iFlt)), 
     1         nInten(iProb)

       do iAtten = 1,nAtten(iProb,jType)
          do iFtype=1,nFtype(iFlt)
           do iParam=1,nParamVar(iFlt,iFltWidth)
             write (11,'( 20e15.6 )')  (tempHaz(iParam,j,iProb,iAtten,iFtype),j=1,nInten(iProb))
           enddo
          enddo
       enddo
      enddo
    
      return
      end
c  -------------------------------------------------------------------
      subroutine WriteTempHaz1 ( tempHaz1,nParamVar, nInten, 
     1                nProb, nAtten, iFlt, jtype, nftype, iFltWidth, nWidth )
 
      include 'pfrisk.h'
      include 'declare1.h'

      real*8 tempHaz1(MAXPARAM,MAX_INTEN, MAX_PROB,MAX_FTYPE)
c      integer nInten(1), nWidth(1), nParamVar(MAX_FLT,1)
c      integer nAtten(MAX_PROB,4)

      nwr = 27

      do iProb=1,nProb

      write (27,'( 20i5)') iFlt, iProb, nAtten(iProb,jType), nWidth(iFlt), 
     1         nFtype(iFlt),(nParamVar(iFlt,i),i=1,nWidth(iFlt)), 
     1         nInten(iProb)

c       do iAtten = 1,nAtten(iProb,jType)
          do iFtype=1,nFtype(iFlt)
           do iParam=1,nParamVar(iFlt,iFltWidth)
             write (27, '( 15e12.4)')  (tempHaz1(iParam,j,iProb,iFtype),j=1,nInten(iProb))
           enddo
          enddo
c       enddo
      enddo
    
      return
      end

c  -------------------------------------------------------------------
      subroutine WriteTempHaz2 ( tempHaz2, nInten, nProb, nAtten, nattenType )

      include 'pfrisk.h'
      include 'declare1.h'

      real*8 tempHaz2(4, MAX_INTEN, MAX_PROB,MAX_ATTEN)
c      integer nInten(1), nWidth(1), nParamVar(MAX_FLT,1)
c      integer nAtten(MAX_PROB,4)

      nwr = 28

      write (*,'( i5)') nProb
      do iProb=1,nProb
        do jType=1,nattenType
           do iAtten = 1,nAtten(iProb,jType)
               write (28,'( 3i5, 100e12.4 )')  iProb, jType, iAtten, 
     1             (tempHaz2(jType, j,iProb,iAtten),j=1,nInten(iProb))
           enddo
         enddo
      enddo
    
      return
      end

c  -------------------------------------------------------------------
 
      subroutine output_TotalHaz ( isite, sitex, sitey, ati, nInten, 
     1           nFlt, nAtten, risk, fName, jCalc, sigTrunc, csrflag,
     2           attenName, period1, probAct, nWidth, 
     4           m_bar, d_bar, e_bar, riskBins, nMagBins, nDistBins,
     5           nEpsBins, magBins, distBins, epsBins,
     6           al_segWt, minDist, nAttenType, attenType, segWt1,
     7           dirflag,tapflag, intflag, fsys, faultdist,mMagout,
     8           hwflagout, ftype, vs, nmaxmag2, mMagoutwt, specT)

      include 'pfrisk.h'

      real siteX, siteY, ati(MAX_PROB, MAX_INTEN) 
      real*8     risk(MAX_INTEN,MAX_PROB,MAX_FLT)
      real al_segWt(MAX_FLT), segWt1(MAX_FLT), magbins(MAX_MAG), distbins(MAX_DIST), 
     1     epsBins(MAX_EPS)
      real*8 m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN),
     3     e_bar(MAX_PROB,MAX_INTEN)
      real    riskBins(MAX_MAG,MAX_DIST,MAX_EPS,MAX_PROB,MAX_INTEN),
     5     outrisk(MAX_MAG,MAX_DIST,MAX_PROB,MAX_INTEN),
     9     sumbar
      real mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real minDist(MAX_FLT), mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real*8 risk2(MAX_INTEN), sum, risk2t(MAX_INTEN)
      real  period1(4,MAX_PROB),
     1     sigTrunc(MAX_PROB), probAct(MAX_FLT), faultdist(MAX_FLT,MAX_WIDTH,3)
      integer csrflag(MAX_PROB), attenType(MAX_FLT), nAttenType, dirflag(MAX_PROB)
      integer tapflag(MAX_PROB), intflag(4,MAX_PROB)
      integer nMagBins, nDistBins, nEpsBins, fsys(MAX_FLT)
      integer isite, nInten(MAX_PROB), nFlt, jCalc(MAX_PROB,4,MAX_ATTEN), nWidth(MAX_FLT) 
      integer hwflagout(MAX_FLT), nMaxmag2(MAX_FLT), casecount
      real ftype(1), vs, specT(MAX_PROB)
      character*80 fName(MAX_FLT), file1, attenName(4,MAX_ATTEN)

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
      write (nwr,'(a28,5x,a19,a61,a33)') 
     1 'Fault Name                  ',
     2 'nFlt System Number ',
     3 'Rupture Dist   BJF Dist  Campbell Dist nWidth nMaxmag Max Mag',
     4 '  Wts.  Fault Type   HWFlag    Vs'

      do ii=1,nFlt
         do k=1,nWidth(ii)
            do l=1,nMaxMag2(ii)

               write (nwr,1111) fname(ii), ii, fsys(ii),
     1         (faultdist(ii,k,kk),kk=1,3), k,l, mMagout(ii,k,l),
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
      write (18,'(a50)') 'Fault Name                  nFlt System Number'
      write (18,'(1x,a50)') '--------------------------------------------------'
      do ii=1,nFlt
         write (18,'(a30,i5,i5)') fname(ii), ii, fsys(ii)
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
     2          minDist(iFLt),
     3          (risk(iInten,iProb,iFlt),iInten=1,nInten(iProb))
        enddo

c       Write Total Number of Events
        xx = 0.
        write(NWR,'( 2x,''Wt_Total_Events/yr'',20x,2f6.3,f8.1,1x,30e12.4)' )
     1        xx, xx, xx, (risk2( M ), M=1,nInten(iProb))

c       Compute Probability of Exceedance  (Poisson)
        do iInten=1,nInten(iProb)
          risk2t( iInten ) = risk2 (iInten)
          if ( risk2( iInten ) .GT. 0.0001 ) then
            risk2 ( iInten ) = (1. - EXP( -risk2( iInten)))
          endif
        enddo
        
c       Write hazard        
        write(NWR,'( 2x,''Poisson_Prob:  '',23x,2f6.3,f8.1,1x,30e12.4)')
     1        xx, xx, xx, ( risk2( M ), M=1,nInten(iProb))

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
 919  FORMAT(' LOG (AMPLITUDE):',7X,10F10.2)
 920  FORMAT(1X,'AMP:',36x,' pSeg al_Wt  MinDis',30F12.6)
 950  FORMAT('     FAULT',I2,' E(NO/YR)',5X,10E10.3)
 960  FORMAT('     TOTAL   E(NO/YR)',5X,10E10.3)
 975  FORMAT('     TOTAL risk',11X,10E10.3)
      end
      
c  ------------------------------------------------------------------

      subroutine output_hazBins ( isite, sitex, sitey, testInten,  
     1           nInten, nProb, riskBins, jCalc, sigTrunc,csrflag,
     1           nMagBins,
     2           nDistBins, nEpsBins, magBins, distBins, epsBins, 
     3           attenName, period1, m_bar, d_bar, e_bar,
     4           nAttenType, attenType, Xcost_bar, nXcostBins, 
     5           XcostBins, RiskBinsX )

      include 'pfrisk.h'

      real siteX, siteY,  
     1     riskBins(MAX_MAG,MAX_DIST,MAX_EPS,MAX_PROB,MAX_INTEN),
     2     magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS),
     2     testInten(MAX_PROB, MAX_INTEN),
     3     period1(4,MAX_PROB), sigTrunc(MAX_PROB), 
     6     outrisk(MAX_MAG,MAX_DIST,MAX_PROB,MAX_INTEN),
     7      XcostBins(MAX_Xcost),
     8     RiskBinsX(MAX_Xcost,MAX_PROB,MAX_INTEN)
     
      real*8 sumX, sum
      real*8  m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN),
     5     e_bar(MAX_PROB,MAX_INTEN),Xcost_bar(MAX_PROB,MAX_INTEN)


      integer csrflag(MAX_PROB), attenType(MAX_FLT)
      integer isite, nInten(MAX_PROB), jCalc(MAX_PROB,4,MAX_ATTEN) 
      integer nMagBins, nDistBins, nEpsBins, nXcostBins, iXcostBins
      character*80 file1, attenName(4,MAX_ATTEN) 

      write (18,'(1x,a50)') '-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
      write (18,*)
      write (18,'(1x,a34,2xa80)') 'Filename for deaggregation output:', file1
      write (18,*) 
      
      nwr = 14

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
            write (nwr,'(2x,6f7.1,2x,30(e10.3,2x))') epsBins(1),
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
           write (nwr,'(2x,6f7.1,2x,30(e10.3,2x))')epsBins(iEpsBin),
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
 919  FORMAT(' LOG (AMPLITUDE):',7X,10F10.2)
 920  FORMAT(1X,'AMP:',41x,30(F10.3,2x))
 950  FORMAT('     FAULT',I2,' E(NO/YR)',5X,10E10.3)
 960  FORMAT('     TOTAL   E(NO/YR)',5X,10E10.3)
 975  FORMAT('     TOTAL risk',11X,10E10.3)
      end
      
c ----------------------------------------------------------------------
      real function f_magModel(c,x, i) 
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
 
      subroutine Write_BR_Haz ( iFlt, nInten, nBr, nNode, Br_haz, Br_wt1, nProb )

      include 'pfrisk.h'
      integer iFlt, nNode(MAX_NODE), nInten(MAX_PROB), nBr, nProb
      real*8  BR_haz(MAX_INTEN, MAX_PROB,MAX_BRANCH,MAX_NODE)
      real    Br_wt1(MAX_BRANCH,MAX_NODE)


      do iProb=1,nProb
        do jInten = 1, nInten(iProb)
          do iBR=1,nBR
            write (23,'(5i5,100e12.4)') iFlt, iProb, jInten, iBR, 
     1             nNode(iBR), (BR_haz(JInten,iProb,iBr,iNode),iNode=1,nNode(iBR) ),
     2                         (br_wt1(iBr,iNode),iNode=1,nNode(iBR) )
          enddo
        enddo
      enddo

      return
      end
      
c  -------------------------------------------------------------------
 
      subroutine Set_Br_Haz ( nBr, Br_Index, Br_wt, Br_Haz, Br_wt1,     
     1            iFtype, ftype_Wt, nSegModel, segModelWt1, iflt, ifltwidth, 
     2            iParam, nNode, jInten, iProb, iSeg )

      include 'pfrisk.h'
      integer nBR, iFtype, nSegModel(MAX_FLT), iFlt, iFltwidth, iParam
      integer jInten, iProb, iSeg
      integer BR_index(MAX_FLT,20,MAX_WIDTH,MAXPARAM), nNode(MAX_NODE)
      real    BR_wt(MAX_FLT,20,MAX_WIDTH,MAXPARAM), br_wt1(MAX_BRANCH,MAX_NODE)
      real    ftype_wt(MAX_FLT,MAX_N1), segModelWt1(MAX_FLT,100)
      real*8  BR_haz(MAX_INTEN, MAX_PROB,MAX_BRANCH,MAX_NODE)


C     Program is currently set for 13 main branches. 
       nBr = 13
       do iBR=1,nBR
        if ( iBR .le. 11 ) then
          iNode = BR_index(iFlt,iBr,iFltWidth, iParam)
          wt5 = br_wt(iFlt,iBr,iFltWidth, iParam)
          if (iNode .gt. nNode(iBR) ) nNode(iBR) = iNode
          if (iNode .ne. 0 .and. wt5 .ne. 0. ) then
              BR_haz(JInten,iProb,iBr,iNode) = BR_haz(JInten,iProb,iBr,iNode) + wt*mHaz / wt5
              br_wt1(iBr,iNode) = wt5
          endif
        else
          if ( iBR .eq. 12 ) then
            iNode = iFtype
            wt5 = ftype_wt(iFlt,iFtype)
            if (iNode .gt. nNode(iBR) ) nNode(iBR) = iNode
            if (iNode .ne. 0 .and. wt5 .ne. 0. ) then
              BR_haz(JInten,iProb,iBr,iNode) = BR_haz(JInten,iProb,iBr,iNode) + wt*mHaz / wt5
              br_wt1(iBr,iNode) = wt5
            endif
          elseif ( iBR .eq. 13 ) then
             nNode(iBR) = nSegModel(iFlt)
            do iSeg=1,nSegModel(iFlt)
             iNode = iSeg
             wt5 = segModelWt1(iFlt,iSeg)  
             if (iNode .ne. 0 .and. wt5 .ne. 0. ) then
               BR_haz(JInten,iProb,iBr,iNode) = BR_haz(JInten,iProb,iBr,iNode) + wt*mHaz / wt5
              br_wt1(iBr,iNode) = wt5
             endif
            enddo
          endif
        endif    
      enddo  


      return
      end
      
c  -------------------------------------------------------------------
      subroutine writedisthead (ioflag)
      
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
