      program haz45a

c     Probabilisitic Seismic Hazard Program (PSHA)

      include 'pfrisk.h'
      include 'declare1.h'
      integer faultFlag(MAX_FLT,100,MAX_FLT), nDD(MAX_FLT)
      real segWt1(MAX_FLT)
      real fltGrid_X(MAXFLT_DD,MAXFLT_AS), fltGrid_y(MAXFLT_DD,MAXFLT_AS), 
     1     fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_fLen(MAXFLT_DD,MAXFLT_AS)
      real rupGrid_X(MAXFLT_DD,MAXFLT_AS), rupGrid_y(MAXFLT_DD,MAXFLT_AS), 
     1     rupGrid_z(MAXFLT_DD,MAXFLT_AS), hDD(MAX_FLT), hAS(MAX_FLT)
      integer nfltGrid(2), nRupGrid(2), hDDcell, hAScell
      real fltGrid_w(MAXFLT_DD,MAXFLT_AS),  fltGrid_a(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_RJB(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_Rseis(MAXFLT_DD,MAXFLT_AS), fltGrid_Rx(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_HW(MAXFLT_DD,MAXFLT_AS)
      real testsum(1000), sum1(1000,10), dipaverage(1), Rx, Ry, Ry0
      real*8 p1_sum, wt, p1
      real*8 BR_haz(MAX_INTEN, MAX_PROB,MAX_BRANCH,MAX_NODE)
      integer BR_index(MAX_FLT,20,MAX_WIDTH,MAXPARAM), nNode(MAX_NODE)
      integer segModelFlag(MAX_FLT,100), nSegModel(MAX_FLT)
      integer icellRupStrike, icellRupDip, runflag
      real    BR_wt(MAX_FLT,20,MAX_WIDTH,MAXPARAM), br_wt1(MAX_BRANCH,MAX_NODE)
      real segModelWt1(MAX_FLT,100), distDensity2(MAX_GRID), lnDir, lgIo, lgIntenscl
      real sigDirY, sigtemp, lg1, sig1, lat1
      integer n1AS(MAXFLT_AS), n2AS(MAXFLT_AS)
      real phi, tau, medadj, sigadj, phiSSS
      character*80 filebmode
      integer bnum, bnumflag, coefcountRrup, coefcountRjb, iMixture
      real pLocY(MAXFLT_AS), sigmaTotal, sigma1, sigma2
      real*8 prock1, prock2
   
c     Write Program information to the screen.
      write (*,*) '*********************************'
      write (*,*) '*   Hazard Code: Version 45i    *'
      write (*,*) '*          May, 2015            *'
      write (*,*) '*********************************'
      write (*,*)

      write (*,*) 'Enter number of cases to run in batch mode.'
      write (*,*) '         (For a single run enter 0)'
      read (*,*) bnum
      if (bnum .eq. 0) then
          bnumflag = 0 
          bnum = 1
      else
          bnumflag = 1
      endif

      if ( bnumflag .eq. 1) then
         write (*,*) 'Enter the batch mode filename.'
         read (*,'(a80)') filebmode
         open (77,file=filebmode,status='old')
      endif
      
c     Start loop over number of batch mode runs
      do 2000 ibnum=1,bnum
         if (bnumflag .eq. 1) then
            write (*,*) 'Looping over Number of Batch Mode Cases: ', ibnum, bnum
         endif

c     Read Run File
      call RdInput ( nProb, nAttenType, nAtten, jcalc, specT, sigTrunc,
     1               gmScale, dirFlag, nInten, testInten, lgTestInten, 
     2               psCorFlag, minlat, maxlat, minlong, maxlong, distmax,
     3               nMagBins, magBins, nDistBins, distBins, nepsBins, epsBins,
     4               nXcostBins, xcostBins, soilAmpFlag, gm_wt, runflag, sigvaradd,
     5               sCalc, sigfix, ssscalc, bnumflag, cfcoefRrup, cfcoefRjb, 
     6               coefcountRrup, coefcountRjb, iMixture )
 
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
     4     br_index, br_wt, segModelFlag, nSegModel, segModelWt1, runflag )

C     Write out Column header explanation to distance file if runflag=4.
      if (runflag .eq. 4) then
          call writedisthead (1)          
      endif

C     For Runflag=5 case the program will run in batch mode with multiple
C         input files being from batch input file. 
      if (runflag .eq. 5) then
         
      endif
      
c     Loop Over Number of Sites

      read (13,*) nSite
      
      do 1000 iSite = 1, nSite      
            
c      Read site coordinates and properties
       read (13,*) SiteX, SiteY, vs, depthvs10, depthvs15, D25, vrup, forearc

       if (runflag .eq. 4) then
          call writedisthead (2)
       endif

c      Read site-specific site amplification
       if ( soilAmpFlag .eq. 1 ) then
        call RdSoilAmpModel ( refPeriod, nRefPer, RefGM, nRefGM, RefGM_mag, nRefMag, amp )
       endif

c      Output1 file which will contain the individual hazard curves. Note output files are not created for 
C      distance run case (i.e., runflag=4).
          read (13,'( a80)') file1
       if (runflag .ne. 4) then
          open (11,file=file1,status='unknown')
       endif
c      Open Output file2 for Probability of Magnitude Density for each parameter combination.
          read (13,'( a80)') file2
       if (runflag .ne. 4) then
          open (17,file=file2,status='unknown')
          write (17,'(a25,i10)') 'Total Number of Faults = ', nFlt
          write (17,*)
       endif

c      Open Output file3 for tornado plot data. Note output files are not created for 
C      distance run case (i.e., runflag=4).
          read (13,'( a80)') file2
       if (runflag .ne. 4) then
          open (23,file=file2,status='unknown')
          write (23,'(i10,2x,''nFlt'')') nFlt
          write (23,'(i10,2x,''nProb'')') nProb
          write (23,'( 100i5)') ( nInten(iProb), iProb=1,nProb )
       endif

c      Initialize Haz Arrays
       call InitHaz ( Haz )
       call InitHaz ( magbar1 )
       call InitHazBins ( HazBins )
       call InitHazBinsX ( HazBinsX )

c      Initialize Mean Deagg values for this site
       call InitDeagg ( m_bar, d_bar, e_bar, Xcost_bar )

c      Sum Over Number of Faults
       do 900 iFlt = 1, nFlt
        write (*,'(2x,''Site = '',i4,'', '',''iFlt ='',4i5)') iSite, iFlt, nFlt, sourceType(iFlt), ibnum

C     For distance case, read in hypocenter location in km down-dip and along strike from separate file. 
C     Only read the locations in for the first site which will be used for all additional sites in distance case. 
         if (runflag .eq. 4) then   
            if (iSite .eq. 1) then
                read (29,*) hDD(iFlt), hAS(iFlt)
             endif
         endif

c       Initialize branch hazard sums
        call InitBrhaz ( nProb, nInten, nNode, br_haz )

        p1_sum = 0.

c       Initialize closest distance for each distance metric
        minDist = 1.e10
        do iWidth=1,nWidth(iFlt)
          do j=1,3
             FaultDist(iFlt,iWidth,j) = 1.e10
          enddo
        enddo

C       Initialize temp hazard array for this source
        call Init_tempHaz ( tempHaz )

c       Loop over alternative Fault Widths (epistemic)
c       (This changes the geometry of the fault)
        do 860 iFltWidth=1,nWidth(iFlt)
        	
c        Set bottom of fault for standard faults
          if ( sourceType(iFlt) .eq. 1. ) then
            call SetFltBottom (iCoor, iFlt, nfp, dip(iFlt,iFltWidth,1), 
     1                         faultWidth(iFlt,iFltWidth), fZ, flat, flong, nDD)
          endif

c        Convert Long, Lat to x,y in km and put into new array (1-D)
         call ConvertCoordinates2 (nfp(iFlt), iFlt, iCoor, grid_n(iFlt), 
     1           sourceType(iFlt), nDD(iFlt), siteX, siteY, fLat, fLong, fZ, 
     2           grid_lat, grid_long, grid_dlat, grid_dlong, nPts, xFlt, yFlt, 
     3           zFlt, grid_x, grid_y, grid_dx, grid_dy, x0, y0, z0) 

c        Turn fault into a grid 
         if ( sourceType(iFlt) .eq. 1 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_X, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGRid_Rrup, FltGRid_Rjb, fltGrid_Rseis, fltGrid_Rx, FltGrid_HW,  
     3               faultArea, faultLen, aveWidth, minDist, xStep(iFlt), fltGrid_fLen )  
         elseif ( sourceType(iFlt) .eq. 5 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_X, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGRid_Rrup, FltGRid_Rjb, fltGrid_Rseis, fltGrid_Rx, FltGrid_HW,  
     3               faultArea, faultLen, aveWidth, minDist, xStep(iFlt), fltGrid_fLen )  
         elseif ( sourceType(iFlt) .eq. 6 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_X, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGRid_Rrup, FltGRid_Rjb, fltGrid_Rseis, fltGrid_Rx, FltGrid_HW,  
     3               faultArea, faultLen, aveWidth, minDist, xStep(iFlt), fltGrid_fLen )  
         endif

c        Initialize Deterministic Values for this Fault
         call InitFltMax ( maxmag1, minDist, maxInten )

c        Set Sampling of Rupture Area and Rupture Width Distributions
         call initRup ( sigArea, nRupArea, sigMaxArea, areaStep, iFlt)
         call initRup ( sigWidth, nRupWidth, sigMaxWidth, widthStep, iFlt)

c        Compute horizontal distance density function for areal sources (polygons or gridded seismicity)
         if ( sourceType(iFlt) .eq. 2 ) then        
           call CalcDistDensity (nPts, xFlt, yFlt, distDensity,
     1         xStep(iFlt), nLocXAS, x0, y0, sampleStep(iFlt), minDist )
         elseif ( sourceType(iFlt) .eq. 3 ) then
           call CalcDistDensity1 ( iFlt, grid_a, grid_x, grid_y, grid_dx,
     1             grid_dy, grid_n, distDensity, xStep(iFlt), nLocXAS, 
     2             x0, y0, sampleStep(iFlt), minDist )         
         elseif ( sourceType(iFlt) .eq. 4 ) then
           call CalcDistDensity2 ( iFlt, grid_a, grid_n, distDensity2 )
         endif  
          
c        Compute activity rate: N(Mmin)
         call Set_Rates ( nParamVar, MagRecur, rate, beta, minMag,
     1         maxMag, iFlt, iFltWidth, faultArea, 
     2         RateParam, mpdf_param, magStep, RateType, 
     1         charMeanMo, expMeanMo )

c        Intergrate Over Magnitude (from minMag to maxMag) (Aleatory)
         do iParam=1,nparamVar(iFlt,iFltWidth)
           sum1(iParam,iFltWidth) = 0.
         enddo
         
c     Set max value for distance checking option (i.e. Runflag=4)
         if (runflag .eq. 4) then
            nMag(iFlt) = 1
            minMag(iFlt) = 8.5
         endif

         do 800 iMag = 1, nMag(iFlt)
          mag = minMag(iFlt) + (iMag-0.5) * magStep(iFlt)
          magTotal = mag

c         Set the magnitude bin for deagregating
          call SetBin ( nMagBins, magBins, mag, iMagBin )

c         Compute Probability of mag between mag-magStep/2 and mag+magStep/2 
c         using the magnitude pdf for each parameter variation
          call magProb ( mag, maxMag, minMag, magStep, beta, iFlt, 
     1           pMag, nParamVar, nWidth, MagRecur, 
     2           mpdf_param, ExpMeanMo, CharMeanMo )
         do iParam=1,nparamVar(iFlt,iFltWidth)
           sum1(iParam,iFltWidth) = sum1(iParam,iFltWidth) + pmag(iParam,iFltWidth)
         enddo

c         Echo magnitude integration step over magnitude to the screen 
c         as a check of the programs progress.
          if (runflag .ne. 4) then
             write (*,'( 2x,2I5,f10.3)') IFLT, ifltWidth,MAG
          endif
          
c         Intergrate Over Rupture Area for this mag (aleatory)
          do 750 iArea = 1, nRupArea(iFlt)

c          Compute Rupture Area and Probability of Rupture Area
           call rupDimProb ( mag, coeff_area, sigArea, areaStep, 
     1             sigMaxArea, rupArea, pArea, iFlt, iArea )

c          Intergrate Over Rupture Width for this mag (aleatory)
           do 700 iWidth = 1, nRupWidth(iFlt)

c           Compute Rupture Width and Probability of Rupture Width
            call rupDimProb ( mag, coeff_width, sigWidth, 
     1         widthStep, sigMaxWidth, rupWidth, pWidth, iFlt, iWidth)

c-----------temporary code for Test 3---------
c-----------Christie Hale-------
c-----------written April 16 2015, don't forget to remove--------

c            rupWidth = sqrt(rupArea/2.)    

c------------end temporary code

        call RupDims (sourcetype(iFlt), rupWidth, aveWidth, rupArea, faultLen,
     1                faultWidth(iFlt,iFltWidth), nLocYST1, yStep(iFlt), rupLen)      
 
        call nLocXcells (sourceType(iFlt), nLocXAS, grid_n(iFlt), nfltgrid, fltgrid_w,
     1                   rupWidth, fltgrid_a, ruparea, nLocYST1, nLocX, n1AS, n2AS)
                        
c           Integrate Over Rupture Location - along strike (aleatory)
c           This is along strike for faults and epicentral distance for source zones
            iDepthFlag = 0
            do 650 iLocX = 1, nLocX
            
            call nLocYcells (iLocX, n1AS, sourceType(iFlt), nLocX, distDensity, xStep(iFlt),
     1                       faultWidth(iFlt,iFltWidth), yStep(iFlt), distDensity2, grid_x,
     2                       grid_y, x0, y0, nLocY, pLocX, r_horiz)

             if ( pLocX .eq. 0. ) then
               goto 650
             endif

c            set the probabilities for the depths 
             if ( iDepthFlag .eq. 0 ) then
               call CalcDepthProb ( iDepthModel(iFlt), depthParam, iFlt, pLocY,
     1              sourceType(iFlt), nLocY, yStep(iFlt), zFlt(1,1),
     2              faultWidth(iFlt,iFltWidth), rupWidth, dip(iFlt,iWidth,1) )
               if (sourceType(iFlt).eq.1 .or. sourceType(iFlt).eq.2) then
                 iDepthFlag = 1
               endif
             endif

c            Integrate Over Rupture Location - Down Dip (aleatory)
             do 600 iLocY = 1, nLocY

c             SourceType 1 fixed, assumes hypocenter is in middle of rupture
c             Set the hypocentral depth (is this really ztor??)
              if (sourceType(iFlt) .eq. 1 ) then
                hypoDepth = (iLocY-1.)*ystep(iFlt)*sin(abs(dip(iFlt,iWidth,1))*3.14159/180.0) + zFlt(1,1)
     1          + ((0.5*rupWidth)*sin(abs(dip(iFlt,iWidth,1))*3.14159/180.0))
              elseif (sourceType(iFlt) .eq. 5 ) then
                hypoDepth = fltgrid_Z(iLocY,iLocX)
              elseif ( sourceType(iFlt) .eq. 2 .or. sourceType(iFlt) .eq. 3 ) then
                hypoDepth = (iLocY-0.5)*ystep(iFlt) + grid_top(iFlt,1)
              elseif ( sourceType(iFlt) .eq. 4 ) then
                hypoDepth = (iLocY-0.5)*ystep(iFlt) + grid_top(iFlt,1)
              endif  

c            Find the Closest Distances for this rupture
C            Pass along fault grid locations for calculation of HW and Rx values within CalcDist subroutine.     
             call CalcDist (sourceType(iflt), pscorflag, hypoDepth, RupWidth, RupLen, 
     1             r_horiz, mindepth(iflt), nFltGrid, n1AS, iLocX, iLocY, n2AS,
     2             fltGrid_x, fltGrid_y, fltGrid_z, fltGrid_Rrup, fltGrid_Rjb,
     3             distJB, distRup, ZTOR, distSeismo, distepi, disthypo, HWFlag,
     4             dipavg, n1, n2, Rx, Ry, Ry0, icellRupstrike, icellRupdip)          

c             Set minimum distances for output files.
              if ( distRup .lt. FaultDist(iFlt,iFltWidth,1) ) then
                FaultDist(iFlt,iFltWidth,1)=distRup
              endif
              if ( distJB .lt. FaultDist(iFlt,iFltWidth,2) ) then
                FaultDist(iFlt,iFltWidth,2)=distJB
              endif
              if ( distSeismo .lt. FaultDist(iFlt,iFltWidth,3) ) then
                FaultDist(iFlt,iFltWidth,3)=distSeismo
              endif

              if (sourceType(iFlt) .eq. 1 ) then
                 MinDist = FaultDist(iFlt,iFltWidth,1)
              elseif (sourceType(iFlt) .eq. 5 ) then
                 MinDist = FaultDist(iFlt,iFltWidth,1)
              elseif (sourceType(iFlt) .eq. 6 ) then
                 MinDist = FaultDist(iFlt,iFltWidth,1)
              endif
c             Check if the distance is greater than Max dist for hazard     
              if ( distRup .gt. distmax .and. distJB .gt. distmax) goto 600

c             Set the distance bin for de-aggregation (using rupture distance)                     
              call SetBin ( nDistBins, distBins, distRup,iDistBin)

c             Loop over ftypes
              if (runflag .eq. 4) nFtype(iflt) = 1
              do 561 iFtype=1,nFtype(iFlt)

c             Loop Over Number of Problems (e.g. spectral periods)
              if (runflag .eq. 4) nProb=1
              do 560 iProb=1,nProb
               jType = attenType(iFlt)
               if (runflag .eq. 4) nAtten(iProb,jtype)=1
               do 550 iAtten = 1,nAtten(iProb,jType)


c               Check for negative jcalc values which will set the corresponding sigma to 
C               either a fixed value or sigma from another model.
                sigflag = 0
                if (jcalc(iProb,jType,iAtten) .lt. 0) then
                   jcalc1 = abs(jcalc(iProb,jType,iAtten) )
                   scalc1 = scalc(iProb,jtype,iAtten) 
                   sigfix1 = sigfix(iProb,jType,iAtten)
                   ssscalc1 = ssscalc(iProb,jType,iAtten)
C                Check for either fixed sigma value (scalc1<0) or other sigma model
                   if (scalc1 .lt. 0) then
                      sigflag = 2
                   else
                      sigflag = 1
                   endif
                else
                   jcalc1 = jcalc(iProb,jType,iAtten) 
                endif

               dipaverage(1) = dipavg*180.0/3.14159  

c              Call for median ground motions
               call meanInten ( distRup, distJB, distSeismo,
     1               HWFlag, mag, jcalc1, specT(iProb),  
     2               lgInten,sigmaY, ftype(iFlt,iFtype), attenName, period1, 
     3               iAtten, jType, vs, hypodepth, intflag, AR, dipaverage(1),
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefrrup, cfcoefrjb, Ry0 )

                lgIntenscl = lgInten + gmScale(iProb,jType,iAtten)

C               Second call for different sigma model 
                if (sigflag .eq. 1) then
               call meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, scalc1, specT(iProb),  
     2               temp, sigmaY, ftype(iFlt,iFtype), sigmaName, period1, 
     3               iAtten, jType, vs, hypodepth, intflag, AR, dipaverage(1),
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi, 
     6               cfcoefrrup, cfcoefrjb, Ry0 )

C      Call for Single Station Sigma Phi if requested (i.e., SssCalc1>0)
C             Phi only Models: 0 < ssscalc1 < 100
                     if (ssscalc1 .gt. 0 .and. ssscalc1 .lt. 100) then
                        call sssphimodel (ssscalc1, specT(iProb), mag, Rrup, phiSSS )
                        sigmaY = sqrt (tau*tau + phiSSS*phiSSS)
                        phi = phiSSS

C             Tau only Models: 100 < ssscalc1 < 200
C              Tau Model - base case (July 2014)
                     elseif (ssscalc1 .eq. 100 ) then
                        if (mag .lt. 7.0) then
                            tauSSS = 0.386 + (mag-5.0)*(0.338-0.386)/(2.0)
                        else
                            tauSSS = 0.338
                        endif
                        sigmaY = sqrt (tauSSS*tauSSS + phi*phi)
                        tau = tauSSS
C              Tau Model - Lower Eps Case (July 2014)
                     elseif (ssscalc1 .eq. 101 ) then
                        if (mag .lt. 7.0) then
                            tauSSS = 0.226 + (mag-5.0)*(0.226-0.226)/(2.0) 
                        else
                            tauSSS = 0.226
                        endif
                        sigmaY = sqrt (tauSSS*tauSSS + phi*phi)
                        tau = tauSSS
C              Tau Model - Upper Eps Case (July 2014)
                     elseif (ssscalc1 .eq. 102 ) then
                        if (mag .lt. 7.0) then
                            tauSSS = 0.539 + (mag-5.0)*(0.443-0.539)/(2.0) 
                        else
                            tauSSS = 0.443 
                        endif
                        sigmaY = sqrt (tauSSS*tauSSS + phi*phi)
                        tau = tauSSS

C        Combined Phi and Tau Models
C             Phi and Tau (Central) Models: 200 < ssscalc1 < 300
                     elseif (ssscalc1 .gt. 200 .and. ssscalc1 .lt. 300) then
                        ssscalctemp = ssscalc1 - 200
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.386 + (mag-5.0)*(0.338-0.386)/(2.0)
                        else
                            tauSSS = 0.338
                        endif
                        sigmaY = sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS
                        phi = phiSSS

C             Phi and Tau (Low) Models: 300 < ssscalc1 < 400
                     elseif (ssscalc1 .gt. 300 .and. ssscalc1 .lt. 400) then
                        ssscalctemp = ssscalc1 - 300
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.226 + (mag-5.0)*(0.226-0.226)/(2.0) 
                        else
                            tauSSS = 0.226
                        endif
                        sigmaY = sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS
                        phi = phiSSS

C             Phi and Tau (High) Models: 400 < ssscalc1 < 500
                     elseif (ssscalc1 .gt. 400 .and. ssscalc1 .lt. 500) then
                        ssscalctemp = ssscalc1 - 400
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.539 + (mag-5.0)*(0.443-0.539)/(2.0) 
                        else
                            tauSSS = 0.443 
                        endif
                        sigmaY = sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS
                        phi = phiSSS


C        Combined Phi and Tau Models scaled by 0.8
C             Phi and Tau (Central) Models: 500 < ssscalc1 < 600
                     elseif (ssscalc1 .gt. 500 .and. ssscalc1 .lt. 600) then
                        ssscalctemp = ssscalc1 - 500
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.386 + (mag-5.0)*(0.338-0.386)/(2.0)
                        else
                            tauSSS = 0.338
                        endif
                        sigmaY = 0.8*sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS*0.8
                        phi = phiSSS*0.8

C             Phi and Tau (Low) Models: 600 < ssscalc1 < 700
                     elseif (ssscalc1 .gt. 600 .and. ssscalc1 .lt. 700) then
                        ssscalctemp = ssscalc1 - 600
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.226 + (mag-5.0)*(0.226-0.226)/(2.0) 
                        else
                            tauSSS = 0.226
                        endif
                        sigmaY = 0.8*sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS*0.8
                        phi = phiSSS*0.8

C             Phi and Tau (High) Models: 700 < ssscalc1 < 800
                     elseif (ssscalc1 .gt. 700 .and. ssscalc1 .lt. 800) then
                        ssscalctemp = ssscalc1 - 700
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.539 + (mag-5.0)*(0.443-0.539)/(2.0) 
                        else
                            tauSSS = 0.443 
                        endif
                        sigmaY = 0.8*sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS*0.8
                        phi = phiSSS*0.8



C        Combined Phi and Tau Models scaled by 1.2
C             Phi and Tau (Central) Models: 800 < ssscalc1 < 900
                     elseif (ssscalc1 .gt. 800 .and. ssscalc1 .lt. 900) then
                        ssscalctemp = ssscalc1 - 800
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.386 + (mag-5.0)*(0.338-0.386)/(2.0)
                        else
                            tauSSS = 0.338
                        endif
                        sigmaY = 1.2*sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS*1.2
                        phi = phiSSS*1.2

C             Phi and Tau (Low) Models: 900 < ssscalc1 < 1000
                     elseif (ssscalc1 .gt. 900 .and. ssscalc1 .lt. 1000) then
                        ssscalctemp = ssscalc1 - 900
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.226 + (mag-5.0)*(0.226-0.226)/(2.0) 
                        else
                            tauSSS = 0.226
                        endif
                        sigmaY = 1.2*sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS*1.2
                        phi = phiSSS*1.2

C             Phi and Tau (High) Models: 1000 < ssscalc1 < 1100
                     elseif (ssscalc1 .gt. 1000 .and. ssscalc1 .lt. 1100) then
                        ssscalctemp = ssscalc1 - 1000
                        call sssphimodel (ssscalctemp, specT(iProb), mag, Rrup, phiSSS )
                        if (mag .lt. 7.0) then
                            tauSSS = 0.539 + (mag-5.0)*(0.443-0.539)/(2.0) 
                        else
                            tauSSS = 0.443 
                        endif
                        sigmaY = 1.2*sqrt (tauSSS*tauSSS + phiSSS*phiSSS)
                        tau = tauSSS*1.2
                        phi = phiSSS*1.2

                     endif


                elseif (sigflag .eq. 2) then
                    sigmaY = sigfix1
                endif

C               Adjust sigma value if sigvaradd .ne. 0.0
                if (sigvaradd(iProb,jType,iAtten) .lt. 0.0) then
c               Check for reduction of sigma is greater than 0.0.
cnjg                   if ( (sigmaY+sigvaradd(iProb,jType,iAtten)) .le. 0.0) then
                   if ( (sigmaY*sigmaY+sigvaradd(iProb,jType,iAtten)) .le. 0.0) then
                      sigmaY = 0.0
                   else
                      sigmaY = sqrt(sigmaY*sigmaY + sigvaradd(iProb,jType,iAtten) )
                   endif
                elseif (sigvaradd(iProb,jType,iAtten) .gt. 0.0) then
                   sigmaY = sqrt(sigmaY*sigmaY + sigvaradd(iProb,jType,iAtten) )
                endif

C               Set Sigma = 0 if requested by a negative jcalc value. 
c                if (sigflag .eq. 1) then
c                  sigmaY = 0.0001
c                endif
 
                sigmaTotal = sigmaY

C      Application of Directivity model if requested. 
                if ( fltDirect(iFlt) .eq. 1) then
                   if (dirflag(iProb) .ge. 1 .and. dirflag(iProb) .lt. 100 .and. mag .gt. 5.6          
     1                      .and. specT(iProb) .ge. 0.50 ) then

C      First set up the number of hypocenter locations for a given fault rupture area
C      If there are less than 10 cells in either along strike or along dip direction
C      just use each cell. Otherwise take 10 locations along strike and dip
                   call SetnRupLoc ( n1, n2, nHypoX, pHypoX, nHypoXStep, 
     1                               nHypoZ, pHypoZ, nHypoZstep ) 

C      Directivity not needed based on fault flag, magnitude, spectral period range or Directivity model (i.e., JWL not hypo dependent). 
C      Set -->  nHypoX = nHypoZ = 1
                   else
                      nHypoX = 1
                      pHypoX = 1.
                      nHypoXStep = 1

                      nHypoZ = 1
                      pHypoZ = 1.
                      nHypoZStep = 1
                   endif
                else
                   nHypoX = 1
                   pHypoX = 1.
                   nHypoXStep = 1

                   nHypoZ = 1
                   pHypoZ = 1.
                   nHypoZStep = 1
                endif
                
c               Loop over hypocenter location along strike (aleatory)
                do 540 iHypoX=1,nHypoX,nHypoXstep

c                Loop over hypocenter location down dip (aleatory)
                 do 530 iHypoZ=1,nHypoZ,nHypoZstep

C     Call the DirModelDist Subroutine to compute necessary directivity Distance values. 
                if (runflag .eq. 4) then 
                   if (sourceType(iFlt) .eq. 1 .or. sourceType(iFlt) .eq. 5) then
                      hDDcell = int(hDD(iFlt)/xstep(iFlt))
                      hAScell = int(hAS(iFlt)/xstep(iFlt))
                      if (hDDcell .eq. 0) hDDcell=1
                      if (hAScell .eq. 0) hAScell=1
                      if (hDDcell .gt. n1) then
                         hDDcell=n1
                         write (*,*) 'Caution - Hypocenter location downdip greater than fault width downdip!!!'
                         write (*,*) 'Hypocenter placed at bottom of downdip fault plane.'
                      endif   
                      if (hAScell .gt. n2) then
                         hAScell=n2
                         write (*,*) 'Caution - Hypocenter location along strike greater than fault length!!!'
                         write (*,*) 'Hypocenter place at end of fault plane.'
                      endif
                   else
                      hDDcell = 1
                      hAScell = 1
                   endif
                   call DetermDist (hAScell, hDDcell, icellRupStrike, icellRupdip, 
     1                      fltgrid_x, fltgrid_y, fltgrid_z, n2, n1, dipaverage(1),
     2                      iLocY, iLocX, RupLen, RupWidth, x0, y0, z0, 
     3                      edist, hdist, slit, azp1p2, ystep(iFlt), dlit, phiang, fltGrid_rRup,
     4                      s2Site, Rx, astrike )
                endif

C     Now make the call to the NGA rupture directivity Subroutine if applicable
                if ( fltDirect(iFlt) .eq. 1) then
                  if (dirflag(iProb) .ge. 1 .and. mag .ge. 5.6          
     1                      .and. specT(iProb) .ge. 0.50 ) then

                    if (dirflag(iProb) .ge. 1 .and. dirflag(iProb) .lt. 100 ) then 

                       call ngaRDirmodel (iHypoX, iHypoZ, icellRupStrike, icellRupdip, 
     1                      fltgrid_x, fltgrid_y, fltgrid_z, n2, n1, dipavg,
     2                      iLocY, iLocX, fltgrid_rrup, Rx, HWFlag, ftype(iFlt,iFtype),
     3                      dirflag(iProb), lnDir, specT(iProb), mag, sigDirY, RupLen, RupWidth )

                       lgInten = lgIntenscl + lnDir
                       sigmaTotal = sigmaY - sigDirY
C     Apply JWL Directivity model. For Median adjustment use (101) for Sigma adjustment use (102)
C     Note: currently based on if statement above this will only get applied for periods greater than 0.5 sec. 
                    elseif (dirflag(iProb) .ge. 100 ) then

                       if (dirflag(iProb) .le. 103) then
                          call DirJWL_V2 (specT(iProb), DistRup, Rx, Ry, Ruplen, Mag, ftype(iFlt,iFtype), 
     1                      RupWidth, Dipaverage(1), HWflag, medadj, sigadj )
                       elseif (dirflag(iProb) .gt. 103) then
                          call DirJWL_V3 (specT(iProb), DistRup, Rx, Ry, Ruplen, Mag, ftype(iFlt,iFtype), 
     1                      RupWidth, Dipaverage(1), HWflag, medadj, sigadj )
                       endif

C     Call to Version 1 of Directivity Model (Not currently used but placeholder for Dirflag kept).
                          if (dirflag(iProb) .eq. 101) then
                             lgInten = lgIntenscl + medadj 
                          elseif (dirflag(iProb) .eq. 102) then  
                             sigmaTotal = sqrt(sigmaY*sigmaY + sigadj*sigadj)
                             if (sigmaTotal .lt. 0.0) Sigmatotal = 0.0                             
                          elseif (dirflag(iProb) .ge. 103) then
                             lgInten = lgIntenscl + medadj 
                             sigmaTotal = sqrt(sigmaY*sigmaY + sigadj*sigadj)
                          endif
                    endif

                  else
                    lgInten = lgIntenscl
                  endif
                else
                  lgInten = lgIntenscl
                endif

                  nSyn_Case(iFlt) = 1
                  probSyn = 1.
c                 Loop over synchronous ruptures (aleatory)
                  if (runflag .eq. 4) nSyn_Case(iFlt)=1
                  do 520 iSyn=1,nSyn_Case(iFlt)
c                  Loop over test ground motion values                  
                   if (runflag .eq. 4) nInten(iProb)=1

                   do 510 jInten = 1, nInten(iProb)

c                   Compute Probability of exceeding test  
                    pRock1 = pxceed3 (lgInten, lgTestInten, sigmaTotal, iProb,jInten,sigTrunc(iProb))
                    pRock2 = 0.5 * pxceed3 (lgInten, lgTestInten, sigma1, iProb,jInten,sigTrunc(iProb))
     1                       + 0.5 * pxceed3 (lgInten, lgTestInten, sigma2, iProb,jInten,sigTrunc(iProb))

                    if ( iMixture .eq. 0 ) then
                      pRock = pxceed3 (lgInten, lgTestInten, sigmaTotal, iProb,jInten,sigTrunc(iProb))

                    else 
                      sigma1 = sigmaTotal*0.8
                      sigma2 = sigmaTotal*1.2
                      pRock = 0.5 * pxceed3 (lgInten, lgTestInten, sigma1, iProb,jInten,sigTrunc(iProb))
     1                        +0.5 * pxceed3 (lgInten, lgTestInten, sigma2, iProb,jInten,sigTrunc(iProb))    
                    endif

c                   Compute number of standard deviations (epsilon) to reach test level
                    epsilon1 = ( lgtestInten(iProb,jInten)-lgInten )/sigmaTotal
c                    write (*,'( i5,2f12.5)') iMixture, prock2/prock1, epsilon1
   
c                   Set bin for epsilon deaggregation
                    call SetBin_Eps ( nEpsBins, epsBins, epsilon1, iepsBin)
                        
c                   Loop over parameter variations (epistemic)
                    if (runflag .eq. 4) nParamVar(iFlt,iFltWidth) = 1
                    do 500 iParam=1,nParamVar(iFlt,iFltWidth)

c                    Set the weight for this set of parameters (epistemic)
                     wt = RateParamWt(iFlt,iParam,iFltWidth) 
     1                 * magRecurWt(iFlt,iParam,iFltWidth) 
     2                 * faultWidthWt(iFlt,iFltWidth)
     2                 * maxMagWt(iFlt,iParam,iFltWidth) 
     2                 * ftype_wt(iFlt,iFtype) 

c                    Set the weight for this set of parameters (epistemic)
                     wt2(iFlt,iParam,iFltWidth) = RateParamWt(iFlt,iParam,iFltWidth) 
     1                 * magRecurWt(iFlt,iParam,iFltWidth) 
     2                 * faultWidthWt(iFlt,iFltWidth)
     2                 * maxMagWt(iFlt,iParam,iFltWidth) 

c                    Set up weight array for later output.
                     wtout(iFlt,iParam,iFltWidth,iFtype) = wt
         
c                    Set probability of this earthquake (w/o gm)

                     p1 = pMag(iParam,iFltWidth)*pArea*pWidth*pLocX*pLocY(iLocY)*phypoX*phypoZ*probSyn
                    
c                    Sum up probability (w/o ground motion) as a check
                     if ( iAtten .eq. 1 .and. iProb .eq. 1
     1                   .and. jInten .eq. 1) then
                       p1_sum = p1_sum + wt*p1                       
                     endif
                           
c                    Add weight for aleatory rupture segmentation
                     wt = wt * al_segWt(iFlt)

c                    Compute Marginal Rate of Occurance
                     mHaz = rate(iParam,iFltWidth) * prock * p1 * probAct(iFlt)                     
                     wt = wt * gm_wt(iProb,jType,iAtten)*segwt1(iFLt)

c                    Add marginal rate of exceed to total 
                     Haz(jInten,iProb,iFlt) = Haz(jInten,iProb,iFlt) + mHaz*wt

                     HazBins(iMagBin,iDistBin,iEpsBin,iProb,jInten) = 
     1                      HazBins(iMagBin,iDistBin,iEpsBin,iProb,jInten) + dble(mHaz*wt)

c  Note: directivity deaggregation was removed from this version
c                     HazBinsX(iXcost,iProb,jInten) = HazBinsX(iXcost,iProb,jInten) + dble(mHaz*wt)
     
c                    Add to mean deagg 
                     m_bar(iProb,jInten) = m_bar(iProb,jInten) + mHaz*wt*magTotal
                     d_bar(iProb,jInten) = d_bar(iProb,jInten) + mHaz*wt*distRup
                     e_bar(iProb,jInten) = e_bar(iProb,jInten) + mHaz*wt*epsilon1
                     Xcost_bar(iProb,jInten) = Xcost_bar(iProb,jInten) + mHaz*wt*Xcost

c                    Set up branch hazard curves for later output for fractile analysis.
c                    Skip this step for distance case (i.e., runflag=4).
                     if (runflag .ne. 4) then
                        call Set_Br_Haz (nBr, Br_Index, Br_wt, Br_Haz, Br_wt1,     
     1                                iFtype, ftype_Wt, nSegModel, segModelWt1, iflt, ifltwidth, 
     2                                iParam, nNode, jInten, iProb, iSeg )                      
                     endif

c                    Save Marginal Hazard to temp array for fractile output
c                     tempHaz(iParam,iFltWidth,jInten,iProb,iAtten) = mHaz
c     1                        + tempHaz(iParam,iFltWidth,jInten,iProb,iAtten)
                     tempHaz(iParam,iFltWidth,jInten,iProb,iAtten,iFtype) = mHaz
     1                        + tempHaz(iParam,iFltWidth,jInten,iProb,iAtten,iFtype)

 500                continue
 510               continue
 520              continue
 530             continue
 540            continue
 550           continue
 560          continue
 561          continue
 600         continue
 650        continue
 651        continue
 700       continue
 750      continue


c     Write out prob. density functions with rate and wts applied.
       do iWidth=1,nWidth(iFlt)
          do iParam=1,nParamVar(iFlt,iWidth)
c             rout(iParam,iWidth) = (rate(iParam,iWidth)*pmag(iParam,iWidth))
             rout(iParam,iWidth) = (rate(iParam,iWidth)*pmag(iParam,iWidth))
          enddo

          write (17,'(2i4,f8.3,70e15.8)') iFlt, iwidth, mag, (wt2(iFlt,iParam,iWidth),
     1                                    iParam=1,nParamVar(iFlt,iWidth)),
     2                                    (rout(iParam,iWidth), iParam=1,
     3                                    nParamVar(iFlt,iWidth)) 

       enddo

 800     continue

          write (17,'(2i4,f8.3,2e15.8)') -9, -9, -9.999, 0.0, 0.0
        
 850     shortDist(iFlt) = minDist

c     Write out distance values for runflag=4 case.
         if (runflag .eq. 4) then 
            if (sourcetype(iFlt) .eq. 1 .or. sourcetype(iflt) .eq. 5) then
               write (19,'(i3,3x,a30,4f8.2,i5,6f8.2,4i5,3f8.2,2i5,14f8.2)') iSite, fname(iflt), FaultDist(iFlt,iFltWidth,1),
     1             FaultDist(iFlt,iFltWidth,2), FaultDist(iFlt,iFltWidth,3), Rx, HWFlag, s2site,
     2             faultlen, avewidth, astrike*180.0/3.14159, dipaverage(1), xstep(iflt), icellRupDip, nFltGrid(1), 
     3             icellRupstrike, nFltGrid(2),fltGrid_fLen(iCellRupDip,iCellRupstrike), hDD(iFlt), hAS(iFlt), 
     4             hDDcell, hAScell, eDist, hDist, slit, azp1p2, slit/faultlen,
     6             dlit, phiang, dlit/avewidth,
     5             fltGrid_x(hDDcell,hAScell), fltGrid_y(hDDcell,hAScell), fltGrid_z(hDDcell,hAScell), 
     5             fltGrid_x(iCellRupDip,iCellRupStrike), fltGrid_y(iCellRupDip,iCellRupStrike), 
     6             fltGrid_z(iCellRupDip,iCellRupStrike)
            else
               write (19,'(i3,3x,a30,4f8.2,i5,48x,20x,24x,10x,2f8.2)') iSite, fname(iflt), FaultDist(iFlt,iFltWidth,1),
     1             FaultDist(iFlt,iFltWidth,2), FaultDist(iFlt,iFltWidth,3), FaultDist(iFlt,iFltWidth,2), HWFlag, 
     2             FaultDist(iFlt,iFltWidth,2), FaultDist(iFlt,iFltWidth,1)
            endif
         endif

 860    continue

c       Write temp Haz array to file
        if (runflag. ne. 4) then
           call WriteTempHaz ( tempHaz, nParamVar, nWidth, nInten, nProb, 
     1        nAtten, iFlt, attenType(iFlt), nFtype )
        endif

c       Write p1_sum as a check
        if (runflag .ne. 4) then
           write (*,'( 2x,'' Site = '',i5,2x,'' iFlt = '',i5,'' p1sum ='',f10.5, i5)') iSite, iflt, p1_sum, nFLt
        endif
c        write (*,'( 2x,''mindist ='',f10.4)') mindist

c       Write out the tornado plot numbers for this fault
        if (runflag .ne. 4) then
           call Write_Br_Haz ( iFlt, nInten, nBr, nNode, Br_haz, Br_wt1, nProb )
        endif

 900   continue

c      close outfiles 1 and 2 for this site
       close (11)
       close (17)

c      Write out the mean Haz, skipping this step for distance runs (i.e., runflag=4)
       if (runflag .eq. 4) then
          read (13,*) 
          read (13,*) 
          write (19,*)
       else          
         call output_TotalHaz ( isite, sitex, sitey, testInten, nInten,
     1       nFlt, nProb, Haz, fName, jCalc, sigTrunc, csrflag,
     2       attenName, period1, probAct, nWidth, m_bar, d_bar, e_bar,
     3       HazBins, nMagBins, nDistBins, nEpsBins, magBins, distBins,
     4       epsBins, al_segWt, shortDist, nAttenType, attenType,
     5       segwt1, dirflag, tapflag,intflag, fsys, faultdist,
     6       mMagout, hwflagout, ftype, vs, nMaxmag2, mmagoutWt, specT)
    
c      Write out the deagrregated hazard
          call output_HazBins ( isite, sitex, sitey, testInten, nInten,
     1       nProb, HazBins, jCalc, sigTrunc, csrflag,
     1       nMagBins, nDistBins,
     2       nEpsBins, magBins, distBins, epsBins,
     3       attenName, period1, m_bar, d_bar, e_bar,
     4       nAttenType, attenType, Xcost_bar, nXcostBins, XcostBins,
     5       HazBinsX)
       endif

 1000 continue

      if (runflag .eq. 4) then
         close (19)
         close (29)
      endif

 2000 continue
      close (77)

      stop
      end

      
      
      
