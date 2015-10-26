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
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_RJB(MAXFLT_DD,MAXFLT_AS)
      real testsum(1000), sum1(1000,10), dipaverage(1), Rx, Ry, Ry0
      real*8 p1_sum, wt, p1
      real*8 BR_haz(MAX_INTEN, MAX_PROB,MAX_BRANCH,MAX_NODE)
      integer BR_index(MAX_FLT,20,MAX_WIDTH,MAXPARAM), nNode(MAX_NODE)
      integer segModelFlag(MAX_FLT,100), nSegModel(MAX_FLT)
      integer icellRupStrike, icellRupDip, runflag
      real BR_wt(MAX_FLT,20,MAX_WIDTH,MAXPARAM), br_wt1(MAX_BRANCH,MAX_NODE)
      real segModelWt1(MAX_FLT,100), distDensity2(MAX_GRID), lnDir, lgIo, lgIntenscl
      real sigDirY, sigtemp, lg1, sig1, lat1
      integer n1AS(MAXFLT_AS), n2AS(MAXFLT_AS)
      real phi, tau, medadj, sigadj, phiSSS
      character*80 filebmode
      integer bnum, bnumflag, coefcountRrup, coefcountRjb
      integer iMixture(4, MAX_PROB, MAX_ATTEN)
      real pLocY(MAXFLT_AS), sigmaTotal, sigma1, sigma2
      real*8 prock1, prock2
      real*8 sum0_Mo(MAXPARAM), sum1_Mo(MAXPARAM)
      real Pmag_all(MAXPARAM)
      integer rup1_flag
      
      real*8 tempHaz1(MAXPARAM,MAX_WIDTH,MAX_INTEN, MAX_PROB,MAX_FTYPE)
      real*8 tempHaz2(4, MAX_INTEN, MAX_PROB, MAX_ATTEN)

   
c     Write Program information to the screen.
      write (*,*) '*********************************'
      write (*,*) '*   Hazard Code: Version 45j    *'
      write (*,*) '*         August, 2015          *'
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
     1     mMagout, mMagoutWt, fltDirect, synchron, nsyn_Case, synjcalc,
     1     synmag, syndistRup, syndistJB, synDistSeismo, synHypo,
     2     synftype, synhwflag, synwt, RateType, iDepthModel, depthParam, 
     3     nMaxmag2, segWt1, faultFlag, nDD, nFtype, ftype_wt, 
     4     br_index, br_wt, segModelFlag, nSegModel, segModelWt1, runflag, 
     7     syn_dip, syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0 )
      
c     Loop Over Number of Sites
      read (13,*) nSite    
      do 1000 iSite = 1, nSite      
            
c      Read site coordinates and properties
       read (13,*) SiteX, SiteY, vs, depthvs10, depthvs15, D25, vrup, forearc
       
c      Read site-specific site amplification
       if ( soilAmpFlag .eq. 1 ) then
        call RdSoilAmpModel ( refPeriod, nRefPer, RefGM, nRefGM, RefGM_mag, nRefMag, amp )
       endif

c      Output1 file which will contain the individual hazard curves. 
       read (13,'( a80)') file1
       open (11,file=file1,status='unknown')
       
c      Open Output2 file for Probability of Magnitude Density for each parameter combination.
       read (13,'( a80)') file2
       open (17,file=file2,status='unknown')
       write (17,'(i15, 3x,''nFlt, nWidth'')') nFlt
       write (17,'( 20i5)') (nWidth(iFlt), iFlt=1,nFlt)
       
c      Open Output5 file which will contain the individual source hazard curves averaged over GMPEs.
       read (13,'( a80)') file1
       open (27,file=file1,status='unknown')
       
c      Open Output6 file which will contain the individual GMPE hazard curves over SSC models.
       read (13,'( a80)') file1
       open (28,file=file1,status='unknown')

c      Initialize Haz Arrays to zero
       call InitHaz ( Haz )
       call InitHaz ( magbar1 )
       call InitHazBins ( HazBins )
       call InitHazBinsX ( HazBinsX )

c      Initialize Mean Deagg values for this site
       call InitDeagg ( m_bar, d_bar, e_bar, Xcost_bar )

C      Initialize temp hazard array for GM sensitivity
       call Init_tempHaz2 ( tempHaz2 )
        
c      Sum Over Number of sources 
       do 900 iFlt = 1, nFlt
        write (*,'(2x,''Site = '',i4,'', '',''iFlt ='',4i5)') iSite, iFlt, nFlt, sourceType(iFlt), ibnum

c       Write fault width weights to output2 file
        write (17,'( 20f10.6)') (faultWidthWt(iFlt,iWidth),iWidth=1,nWidth(iFlt))

c       Initialize p1_sum (check of the integration over all source pdfs)
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
        call Init_tempHaz1 ( tempHaz1 )

c       Loop over alternative Fault Widths (epistemic)
c       (This changes the geometry of the source)
        do 860 iFltWidth=1,nWidth(iFlt)
        	
c        Set bottom of fault for standard faults (source type 1)
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
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               minDist, xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )   
         elseif ( sourceType(iFlt) .eq. 5 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               minDist, xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )   
         elseif ( sourceType(iFlt) .eq. 6 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               minDist, xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 ) 
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
         
         do 800 iMag = 1, nMag(iFlt)
          mag = minMag(iFlt) + (iMag-0.5) * magStep(iFlt)
          magTotal = mag

c         Set the magnitude bin for deagregating
          call SetBin ( nMagBins, magBins, mag, iMagBin )

c         Compute Probability of mag between mag-magStep/2 and mag+magStep/2 
c         using the magnitude pdf for each parameter variation
          call magProb ( mag, maxMag, minMag, magStep, beta, iFlt, 
     1           pMag, nParamVar, nWidth, MagRecur, 
     2           mpdf_param, ExpMeanMo, CharMeanMo, rup1_flag )


         do iParam=1,nparamVar(iFlt,iFltWidth)
           sum1(iParam,iFltWidth) = sum1(iParam,iFltWidth) + pmag(iParam,iFltWidth)
         enddo

c         Echo magnitude integration step over magnitude to the screen 
c         as a check of the programs progress.
          write (*,'( 2x,2I5,f10.3)') IFLT, ifltWidth,MAG
          
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
     1                          faultWidth(iFlt,iFltWidth), yStep(iFlt), distDensity2, grid_x,
     2                          grid_y, x0, y0, nLocY, pLocX, r_horiz)

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
c            Pass along fault grid locations for calculation of HW and Rx values within CalcDist subroutine.     
             call CalcDist (sourceType(iflt), pscorflag, hypoDepth, RupWidth, RupLen, 
     1             r_horiz, mindepth(iflt), nFltGrid, n1AS, iLocX, iLocY, n2AS,
     2             fltGrid_x, fltGrid_y, fltGrid_z, fltgrid_x1, fltgrid_y1, 
     3             fltgrid_z1, fltgrid_x2, fltgrid_y2, fltgrid_x3, fltgrid_y3,
     4             fltgrid_x4, fltgrid_y4, fltgrid_z4, fltGrid_Rrup, fltGrid_Rjb,
     5             distJB, distRup, ZTOR, distSeismo, distepi, disthypo, HWFlag,
     6             dipavg, n1, n2, Rx, Ry, Ry0, icellRupstrike, icellRupdip, dip, iFltWidth, iFlt) 

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
              do 561 iFtype=1,nFtype(iFlt)

c             Loop Over Number of Problems (e.g. spectral periods)
              do 560 iProb=1,nProb
               jType = attenType(iFlt)
               do 550 iAtten = 1,nAtten(iProb,jType)

c               Check for negative jcalc values which will set the corresponding sigma to 
C               either a fixed value or sigma from another model.
                sigflag = 0
                if (jcalc(iProb,jType,iAtten) .lt. 0) then
                   jcalc1 = abs(jcalc(iProb,jType,iAtten) )
                   scalc1 = scalc(iProb,jtype,iAtten) 
                   sigfix1 = sigfix(iProb,jType,iAtten)
c                   ssscalc1 = ssscalc(iProb,jType,iAtten)
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

c              Compute the median and sigma of the ground motions
               call meanInten ( distRup, distJB, distSeismo,
     1               HWFlag, mag, jcalc1, specT(iProb),  
     2               lgInten,sigmaY, ftype(iFlt,iFtype), attenName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipaverage(1),
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefrrup, cfcoefrjb, Ry0 )

c               Add epistemic uncertainty term (constant shift) to median
                lgInten = lgInten + gmScale(iProb,jType,iAtten)

C               Second Median Call for Complex/Splay Rupture Median
C               Ground Motions will be SRSS with main rupture median ground motions.
                if (synchron(iFlt) .eq. 0) then
                  nSyn_Case(iFlt) = 1
                  synwt(iFlt,1)= 1.0
                endif

c               Temp: set probability of syn rupture to unity
c               Later, this will be an input
c           write (*,'( 2i5)') iFlt, synchron(iFlt)
                probSyn = 1.

c               Loop over synchronous ruptures (aleatory)          
                do 549 isyn=1,nSyn_Case(iFlt)
                 if (synchron(iFlt) .gt. 0 .and. rup1_flag .eq. 1) then
                  call meanInten ( synDistRup(iFlt,isyn), syndistJB(iFlt,isyn), 
     1                syndistSeismo(iFlt,isyn),
     2                synhwflag(iFlt,isyn), synmag(iFlt,isyn), synjcalc(iFlt), specT(iProb),  
     3                lgIntenS, temp, synftype(iFlt,isyn), attenName, period1, 
     4                iAtten, iProb, jType, vs, synhypo(iflt,1), intflag, AR, syn_dip(iFlt,isyn),
     5                disthypo, depthvs10, depthvs15, D25, tau,
     6                syn_zTOR(iFlt,isyn), theta_site, syn_RupWidth(iFlt,isyn), 
     7                vs30_class, forearc, syn_Rx, phi,
     8                cfcoefrrup, cfcoefrjb, syn_Ry0(iFlt,isyn) )
c                  write (*,'( 2f10.4)') lgInten, lgIntenS
c           pause 'syn flag'

c                 Compute SRSS of median                
                  lgInten = 0.5* alog( exp(lgInten)**2 + exp(lgIntenS)**2 )
                endif

C               Second call got GPE for different sigma model 
                if (sigflag .eq. 1) then
                  call meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, scalc1, specT(iProb),  
     2               temp, sigmaY, ftype(iFlt,iFtype), sigmaName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipaverage(1),
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi, 
     6               cfcoefrrup, cfcoefrjb, Ry0 )

c               Check if a constant, user input sigma, is selected
                elseif (sigflag .eq. 2) then
                  sigmaY = sigfix1
                endif

C               Adjust sigma value if sigvaradd .ne. 0.0
                if (sigvaradd(iProb,jType,iAtten) .ne. 0.0) then
                  sigmaY = sqrt(sigmaY*sigmaY + sigvaradd(iProb,jType,iAtten) )
                endif

c               Check that sigma is not less than zero (0.0001)
                if (sigmaY .lt. 0.0001 ) sigmaY = 0.0001

c               Reset SigmaTotal variable
                sigmaTotal = sigmaY

C               Application of Directivity model. 
                if ( fltDirect(iFlt) .eq. 1 .and. dirflag(iProb) .ge. 1
     1              .and. dirflag(iProb) .lt. 100 .and. mag .gt. 5.6          
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
                    
c               Loop over hypocenter location along strike (aleatory)
                do 540 iHypoX=1,nHypoX,nHypoXstep

c                Loop over hypocenter location down dip (aleatory)
                 do 530 iHypoZ=1,nHypoZ,nHypoZstep

C                 Call to the rupture directivity Subroutine if applicable
                  if ( dirflag(iProb) .ge. 1 .and. fltDirect(iFlt) .eq. 1) then
c                    call Directivity ( dirFlag(iProb), specT, DistRup, 
c     1                    Rx, Ry, Ruplen, mag, ftype(iFlt,iFtype), 
c     2                    RupWidth, Dipaverage(1), HWflag, lgIntenscl,
c     3                    lgInten, sigmaTotal )
                  endif 
  
c                  Loop over test ground motion values                  
                   do 510 jInten = 1, nInten(iProb)

c                   Compute Probability of exceeding test  
                    if ( iMixture(jType,iProb,iAtten)  .eq. 0 ) then
                      pRock = pxceed3 (lgInten, lgTestInten, sigmaTotal, iProb,jInten,sigTrunc(iProb))

                    else
                      sigma1 = sigmaTotal*0.8
                      sigma2 = sigmaTotal*1.2
                      pRock = 0.5 * pxceed3 (lgInten, lgTestInten, sigma1, iProb,jInten,sigTrunc(iProb))
     1                        +0.5 * pxceed3 (lgInten, lgTestInten, sigma2, iProb,jInten,sigTrunc(iProb))    
                    endif

c                   Compute number of standard deviations (epsilon) to reach test level
                    if (sigmatotal .le. 0.0001) then
                      epsilon1 = 0.0
                    else
                      epsilon1 = ( lgtestInten(iProb,jInten)-lgInten )/sigmaTotal
                    endif
   
c                   Set bin for epsilon deaggregation
                    call SetBin_Eps ( nEpsBins, epsBins, epsilon1, iepsBin)
                        
c                   Loop over parameter variations (epistemic)
                    do 500 iParam=1,nParamVar(iFlt,iFltWidth)

c                    Set the weight for this set of parameters (epistemic)
                     wt = RateParamWt(iFlt,iParam,iFltWidth) 
     1                 * magRecurWt(iFlt,iParam,iFltWidth) 
     2                 * faultWidthWt(iFlt,iFltWidth)
     2                 * maxMagWt(iFlt,iParam,iFltWidth) 
     2                 * ftype_wt(iFlt,iFtype) 

c                    Set up weight array for later output.
                     wtout(iFlt,iParam,iFltWidth,iFtype) = wt
         
c                    Set probability of this earthquake (w/o gm) - (aleatory)
                     p1 = pMag(iParam,iFltWidth)*pArea*pWidth*pLocX*pLocY(iLocY)
     1                    *phypoX*phypoZ*probSyn*synwt(iFlt,isyn)
                    
c                    Sum up probability (w/o ground motion) as a check
                     if ( iAtten .eq. 1 .and. iProb .eq. 1 .and. jInten .eq. 1) then
                       p1_sum = p1_sum + wt*p1                       
                     endif
                           
c                    Add weight for aleatory rupture segmentation
                     wt = wt * al_segWt(iFlt)

c                    Compute Marginal Rate of Occurance
                     mHaz = rate(iParam,iFltWidth) * prock * p1 * probAct(iFlt)
c                     write (*,'( 4e12.4)') rate(iParam,iFltWidth) , prock , p1 , probAct(iFlt)
                     wt = wt *segwt1(iFLt)

c                    Add marginal rate of exceed to total
                     Haz(jInten,iProb,iFlt) = Haz(jInten,iProb,iFlt) + mHaz*wt* gm_wt(iProb,jType,iAtten)
c                     write (*,'( 3e12.4)') mHaz, wt,  gm_wt(iProb,jType,iAtten)
c                    pause
                     
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
                     call Set_Br_Haz (nBr, Br_Index, Br_wt, Br_Haz, Br_wt1,     
     1                                iFtype, ftype_Wt, nSegModel, segModelWt1, iflt, ifltwidth, 
     2                                iParam, nNode, jInten, iProb, iSeg )                      

c                    Save Marginal Hazard to temp array for fractile output
                     tempHaz(iParam,iFltWidth,jInten,iProb,iAtten,iFtype) = mHaz
     1                        + tempHaz(iParam,iFltWidth,jInten,iProb,iAtten,iFtype)

                     tempHaz1(iParam,iFltWidth,jInten,iProb,iFtype) = mHaz* gm_wt(iProb,jType,iAtten)
     1                        + tempHaz1(iParam,iFltWidth,jInten,iProb,iFtype)

                     tempHaz2(jType, jInten,iProb,iAtten) = mHaz*wt
     1                        + tempHaz2(jType, jInten,iProb,iAtten)



 500                continue
 510               continue
 530              continue
 540             continue
 549            continue
 550           continue
 560          continue
 561          continue
 600         continue
 650        continue
 651        continue
 700       continue
 750      continue

c         Compute the rate for this magnitude for each epistemic parameter variation
          do iParam=1,nParamVar(iFlt,iFltWidth)
            rout(iParam,iFltWidth) = (rate(iParam,iFltWidth)*pmag(iParam,iFltWidth))
c            write (*,'( i5,2e12.4)') iparam, (rate(iParam,iFltWidth),pmag(iParam,iFltWidth))
          enddo
          mag = minMag(iFlt) + (iMag-0.5) * magStep(iFlt)

c         Write out magnitude rates of occurrence (out2 file)
c         Write header for the first magnitude only
          if (iMag .eq. 1 ) then
            write (17,'( i5,2x,a80)') iFlt, fname(iflt)
            write (17,'( 2i5,''  nMag, nParamVar'')') nMag(iFlt), nParamVar(iFlt,iFltWidth)
c           Set the weight for this set of parameters (epistemic)
            do iParam=1,nParamVar(iFlt,iFltWidth)
              wt2(iParam) = RateParamWt(iFlt,iParam,iFltWidth) 
     1                      * magRecurWt(iFlt,iParam,iFltWidth) 
     2                      * maxMagWt(iFlt,iParam,iFltWidth)
            enddo
            write (17,'( 100f10.6)') (wt2(iParam),
     1            iParam=1,nParamVar(iFlt,iFltWidth)) 
          endif
          write (17,'(2i4,f8.3,500e12.4)') iFlt, iFltwidth, mag, 
     1      (rout(iParam,iFltWidth), iParam=1,nParamVar(iFlt,iFltWidth)) 

 800     continue
 
 850     shortDist(iFlt) = minDist


 860    continue

c       Write temp Haz array to file
        call WriteTempHaz ( tempHaz, nParamVar, nWidth, nInten, nProb, 
     1        nAtten, iFlt, attenType(iFlt), nFtype )
        call WriteTempHaz1 ( tempHaz1, nParamVar, nWidth, nInten, nProb, 
     1        nAtten, iFlt, attenType(iFlt), nFtype )

c       Write p1_sum as a check
        write (*,'( 2x,'' Site = '',i5,2x,'' iFlt = '',i5,'' p1sum ='',f10.5, i5)') iSite, iflt, p1_sum, nFLt

 900   continue

c      close outfiles 1 and 2 for this site
       close (11)
       close (17)

c      Write out the mean Haz
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
     
        call WriteTempHaz2 ( tempHaz2, nInten, nProb, nAtten, nattenType )

 1000 continue

 2000 continue
      close (77)

      stop
      end

      
      
      
