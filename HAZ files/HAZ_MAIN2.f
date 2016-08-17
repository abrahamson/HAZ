      program haz45_2

c     Probabilisitic Seismic Hazard Program (PSHA) 

      implicit none
      include 'pfrisk.h'
      include 'declare1.h' 
   
c     Write Program information to the screen.
      write (*,*) '*********************************'
      write (*,*) '*   Hazard Code: Version 45.2   *'
      write (*,*) '*           Aug, 2016           *'
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
     5               sCalc, sigfix, bnumflag, cfcoefRrup, cfcoefRjb, 
     6               coefcountRrup, coefcountRjb, iMixture, version )

c     read fault File
      call Rd_Fault_Data ( nFlt, fName, minMag, magStep, xStep,
     1     yStep, segModelWt, rateParam, rateParamWt, beta,
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, maxMag,  
     3     maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, ftype, 
     4     sourceType, nRupArea, coeff_area, sigArea, nRupWidth, 
     5     coeff_width, sigWidth, nParamVar, iCoor, minDepth,
     6     fIndex, probAct, nWidth, mpdf_param, al_segWt, attenType, 
     7     sampleStep, grid_a, grid_dlong, grid_dlat, grid_n,
     8     grid_long, grid_lat, grid_top, minlat, maxlat, minlong, 
     9     maxlong, scaleRate, fsys, mMagout, mMagoutWt, fltDirect, 
     1     synchron, nsyn_Case, synjcalc, synmag, syndistRup, 
     2     syndistJB, synDistSeismo, synHypo, synftype, synhwflag, 
     3     synwt, RateType, iDepthModel, depthParam, nMaxmag2, segWt1, 
     4     faultFlag, nDD, nFtype, ftype_wt, br_index, br_wt, 
     5     segModelFlag, nSegModel, segModelWt1, runflag, syn_dip, 
     6     syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, magS7, rateS7,  
     7     DistS7, DipS7, mechS7, ncountS7, version )             
     
c     Loop Over Number of Sites
      read (13,*,err=2100) nSite    
      do 1000 iSite = 1, nSite      
            
c      Read site coordinates and properties
       read (13,*,err=2101) SiteX, SiteY, vs, depthvs10, depthvs15, D25, vrup, forearc
       
c      Read site-specific site amplification
       if ( soilAmpFlag .eq. 1 ) then
        call RdSoilAmpModel ( refPeriod, nRefPer, RefGM, nRefGM, RefGM_mag, nRefMag, amp )
       endif

c      Open Output1 file which will contain the individual hazard curves. 
       read (13,'( a80)',err=2102) file1
       open (11,file=file1,status='unknown')
       write (11, *) '45.2 Haz45.2 Out1 file - individual hazard curves'
       
c      Open Output2 file for Probability of Magnitude Density for each parameter combination.
       read (13,'( a80)',err=2104) file2
       open (17,file=file2,status='unknown')
       write (17, *) '45.2 Haz45.2 Out2 file - magnitude recurrence curves for SSC fractile code'
       write (17,'(i15, 3x,''nFlt, nWidth'')') nFlt
       write (17,'( 20i5)') (nWidth(iFlt), iFlt=1,nFlt)

c      Open Output3 file
       read (13,'( a80)',err=2106) file1
       open (12,file=file1,status='unknown')
       write (12, *) '45.2 Haz45.2 Out3 file - mean hazard curves'
       
c      Open Output4 file
       read (13,'( a80)') file1
       open (14,file=file1,status='unknown')
       write (14, *) '45.2 Haz45.2 Out4 file - deaggregation output'
       
c      Open Output5 file which will contain the individual source hazard curves averaged over GMPEs.
       read (13,'( a80)',err=2105) file1
       open (27,file=file1,status='unknown')
       write (27, *) '45.2 Haz45.2 Out5 file - SSC tornado output file'
       
c      Open Output6 file which will contain the individual GMPE hazard curves over SSC models.
       read (13,'( a80)',err=2106) file1
       open (28,file=file1,status='unknown')
       write (28, *) '45.2 Haz45.2 Out6 file - GMC tornado output file'

c      Open Output7 file which will contain deaggregations for each source.
       read (13,'( a80)',err=2107) file1
       open (29,file=file1,status='unknown')
       write (29, *) '45.2 Haz45.2 Out7 file - deaggregation by source'
       
c      Initialize Haz Arrays to zero
       call InitHaz ( Haz )
       call InitHaz ( magbar1 )
       call InitHazBins ( HazBins )
       call InitHazBinsX ( HazBinsX )

c      Initialize Mean Deagg values for this site
       call InitDeagg ( m_bar, d_bar, e_bar, Xcost_bar )
       call InitDeagg2 ( m_bar_s, rrup_bar_s, rjb_bar_s, rx_bar_s, e_bar_s)

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
        call InitMinDis (iFlt, nWidth(iFlt), MinRrup_temp, MinRjb_temp, MinSeismo_temp, SourceDist)

c       Loop over alternative Fault Widths (epistemic)
c       (This changes the geometry of the source)
        do 860 iFltWidth=1,nWidth(iFlt)

C         Initialize temp hazard array for this source
          call Init_tempHaz ( tempHaz )
          call Init_tempHaz1 ( tempHaz1 )
        	
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
         if ( sourceType(iFlt) .eq. 1 .or. sourceType(iFlt) .eq. 5 .or. sourceType(iFlt) .eq. 6 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )   
         endif
         if ( sourceType(iFlt) .eq. 1 .or. sourceType(iFlt) .eq. 5 ) 
     1      write (18,'( 2x,''fault area (km^2) = '',e12.3)') faultArea

c        Set Sampling of Rupture Area and Rupture Width Distributions
         call initRup ( sigArea, nRupArea, sigMaxArea, areaStep, iFlt)
         call initRup ( sigWidth, nRupWidth, sigMaxWidth, widthStep, iFlt)

c        Compute horizontal distance density function for areal sources (polygons or gridded seismicity)
         if ( sourceType(iFlt) .eq. 2 ) then        
           call CalcDistDensity (nPts, xFlt, yFlt, distDensity,
     1         xStep(iFlt), nLocXAS, x0, y0, sampleStep(iFlt))         
         elseif ( sourceType(iFlt) .eq. 3 ) then
           call CalcDistDensity1 ( iFlt, grid_a, grid_x, grid_y, grid_dx,
     1             grid_dy, grid_n, distDensity, xStep(iFlt), nLocXAS,
     2             x0, y0, sampleStep(iFlt))
         elseif ( sourceType(iFlt) .eq. 4 ) then
           call CalcDistDensity2 ( iFlt, grid_a, grid_n, distDensity2 )
         endif  
          
c        Compute activity rate: N(Mmin)
           call Set_Rates ( sourceType(iFlt), nParamVar, MagRecur, rate, beta, minMag,
     1         maxMag, iFlt, iFltWidth, faultArea, 
     2         RateParam, mpdf_param, magStep, RateType, 
     1         charMeanMo, expMeanMo )

c        Intergrate Over Magnitude (from minMag to maxMag) (Aleatory)
         
c        Set nMag(iFlt) = ncount for Source Type 7 case
         if (sourceType(iFlt) .eq. 7) then
           write (*,*) 'Number of sources for SourceType7 = ', iFlt,nMag(iFlt)
         endif
         
         do 800 iMag = 1, nMag(iFlt)
           if (sourceType(iFLt) .ne. 7) then
             mag = minMag(iFlt) + (iMag-0.5) * magStep(iFlt)
           elseif (sourceType(iFlt) .eq. 7) then
             mag = magS7(iFlt,iMag)
           endif
           magTotal = mag

c         Set the magnitude bin for deagregating
          call SetBin ( nMagBins, magBins, mag, iMagBin )

c         Compute Probability of mag between mag-magStep/2 and mag+magStep/2 
c         using the magnitude pdf for each parameter variation
            call magProb (sourceType(iFlt), mag, maxMag, minMag, magStep, beta, 
     1                    iFlt, iFltWidth, pMag, nParamVar, nWidth, MagRecur, 
     2                    mpdf_param, ExpMeanMo, CharMeanMo, rup1_flag)      

c         Echo magnitude integration step over magnitude to the screen 
c         as a check of the programs progress.
          if (sourceType(iFlt) .ne. 7) then
            write (*,'( 2x,2I5,f10.3)') iflt, ifltWidth, mag
            write (18,'( 2x,2I5,f10.3)') iflt, ifltWidth, mag
          elseif (sourceType(iFlt) .eq. 7) then
            nn10000 = (iMag/10000)*10000
            if (nn10000 .eq. iMag) then
              write (*,'( 2x,2I10,f10.3)') iMag, ncountS7(iFlt), mag                          
            endif
          endif
          
c         Intergrate Over Rupture Area for this mag (aleatory)
          do 750 iArea = 1, nRupArea(iFlt)

c          Compute Rupture Area and Probability of Rupture Area
           call rupDimProb ( sourceType(iFlt), mag, coeff_area, sigArea, 
     1          areaStep, sigMaxArea, rupArea, pArea, iFlt, iArea )

c          Intergrate Over Rupture Width for this mag (aleatory)
           do 700 iWidth = 1, nRupWidth(iFlt)

c           Compute Rupture Width and Probability of Rupture Width
            call rupDimProb ( sourceType(iFlt), mag, coeff_width, sigWidth, 
     1           widthStep, sigMaxWidth, rupWidth, pWidth, iFlt, iWidth)

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
               if (sourceType(iFlt) .le. 2 .or. sourceType(iFlt) .eq. 7) then
                 iDepthFlag = 1
               endif
             endif

c            Integrate Over Rupture Location - Down Dip (aleatory)
             do 600 iLocY = 1, nLocY

c            Find the Closest Distances for this rupture
c            Pass along fault grid locations for calculation of HW and Rx values within CalcDist subroutine.     
             call CalcDist (sourceType(iFlt), pscorflag, nFltGrid, n1AS, iLocX, iLocY, n2AS,
     1             iFltWidth, iFlt, iMag, ystep(iFlt), grid_top, RupWidth, RupLen, r_horiz, mindepth(iFlt), 
     2             fltGrid_x, fltGrid_y, fltGrid_z, fltgrid_x1, fltgrid_y1, fltgrid_z1, 
     3             fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3, fltgrid_z3,
     4             fltgrid_x4, fltgrid_y4, fltgrid_z4, fltGrid_Rrup, fltGrid_Rjb, dip, dipS7,
     5             distS7, HWFlag, n1, n2, icellRupstrike, icellRupdip, hypoDepth, distJB, 
     6             distRup, ZTOR, distSeismo, distepi, disthypo, dipavgd, Rx, Ry, Ry0)
        
c            Set minimum distances for output files.
             call Set_MinDist (sourceType(iFlt), iFlt, iFltWidth, distRup, distJB, distSeismo, 
     1                         SourceDist, MinRrup_temp, MinRjb_temp, MinSeismo_temp)
              
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
C                  Check for either fixed sigma value (scalc1<0) or other sigma model
                   if (scalc1 .lt. 0) then
                     sigflag = 2
                   else
                     sigflag = 1
                   endif
                else
                   jcalc1 = jcalc(iProb,jType,iAtten) 
                endif 

c              Compute the median and sigma of the ground motions
               if (sourceType(iFlt) .ne. 7) then
                 call meanInten ( distRup, distJB, distSeismo,
     1               HWFlag, mag, jcalc1, specT(iProb),  
     2               lgInten, sigmaY, ftype(iFlt,iFtype), attenName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipavgd,
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefrrup, cfcoefrjb, Ry0 )
               elseif (sourceType(iFlt) .eq. 7) then
                 call meanInten ( distRup, distJB, distSeismo,
     1               HWFlag, mag, jcalc1, specT(iProb),  
     2               lgInten, sigmaY, mechS7(iFlt,iMag), attenName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipavgd,
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefrrup, cfcoefrjb, Ry0 )
              endif

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

c                 Compute SRSS of median                
                  lgInten = 0.5* alog( exp(lgInten)**2 + exp(lgIntenS)**2 )
                endif

C               Second call got GPE for different sigma model 
                if (sigflag .eq. 1) then
                  if (sourceType(iFlt) .ne. 7) then
                    call meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, scalc1, specT(iProb),  
     2               temp, sigmaY, ftype(iFlt,iFtype), sigmaName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipavgd,
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi, 
     6               cfcoefrrup, cfcoefrjb, Ry0 )
                  elseif (sourceType(iFlt) .eq. 7) then
                    call meanInten ( distRup, distJB, distSeismo,
     1               HWFlag, mag, jcalc1, specT(iProb),  
     2               lgInten, sigmaY, mechS7(iFlt,iMag), attenName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipavgd,
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefrrup, cfcoefrjb, Ry0 )
                  endif

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

c               Set values for use with directivity                
                lgInten0 = lgInten
                sigma0 = sigmaTotal
                
c               Set default for no randomization of hypocenters                
                nHypoX = 1
                pHypoX = 1.
                nHypoXStep = 1
                nHypoZ = 1
                pHypoZ = 1.
                nHypoZStep = 1

C               Application of Directivity model. 
                if ( fltDirect(iFlt) .eq. 1 .and. dirflag(iProb) .ge. 1
     1              .and. mag .gt. 5.6 .and. specT(iProb) .ge. 0.50 ) then
                 dirFlag1 = 1
                 if ( dirflag(iProb) .lt. 100 ) then     
                   nHypoX = 9
                   nHypoZ = 9
                   pHypoX = 1./ 9.
                   pHypoZ = 1./ 9.
                 endif
                else
                 dirFlag1 = 0
                endif
                    
c               Loop over hypocenter location along strike (aleatory)
                do 540 iHypoX=1,nHypoX,nHypoXstep
                 fs = float(iHypoX) / (nHypoX + 1.)

c                Loop over hypocenter location down dip (aleatory)
                 do 530 iHypoZ=1,nHypoZ,nHypoZstep
                  fd = float(iHypoZ) / (nHypoZ + 1.)

C                 Call to the rupture directivity Subroutine if applicable
                  if ( dirflag1 .eq. 1) then
                    call Directivity ( dirFlag(iProb), specT(iProb), DistRup, zTOR, 
     1                 x0, y0, z0, Rx, Ry, Ry0, mag, ftype(iFlt,iFtype), RupWidth, 
     2                 RupLen, dipavgd, HWflag, dirMed, dirSigma, fltgrid_x, 
     3                 fltgrid_y, fltgrid_z, n1, n2, fs, fd, dpp_flag, 
     4                 iLocX, iLocY)
     
                       write (44,'( 6f8.2 )') mag, RupLen, fs, fd, dirMed, dirSigma

c                   Add directivity to median and sigma
                    lgInten = lgInten0 + dirMed
                    if ( dirSigma .lt. 0. ) then
                      t1 = sigma0**2 - dirSigma**2
                    else
                      t1 = sigma0**2 + dirSigma**2
                    endif
                    if ( t1 .lt. 0. ) t1 = 0.01
                    sigmaTotal = sqrt( t1 )  
                  endif 
  
c                  Loop over test ground motion values                  
                   do 510 jInten = 1, nInten(iProb)

c                   Compute Probability of exceeding test  
                    if ( iMixture(iProb,jType,iAtten)  .eq. 0 ) then
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
                    
c                    Set weights and probabilities for SourceType 7
                     if (sourceType(iFlt) .eq. 7) then
                        p1 = 1.0
                        wt = 1.0
                     endif
                    
c                    Sum up probability (w/o ground motion) as a check
                     if ( iAtten .eq. 1 .and. iProb .eq. 1 .and. jInten .eq. 1) then
                       if (sourcetype(iFlt) .ne. 7) then
                         p1_sum = p1_sum + wt*p1       
                       elseif (sourcetype(iFlt) .eq. 7) then
                         p1_sum = p1_sum + wt*p1/ncounts7(iFlt)       
                       endif                       
                     endif
                           
c                    Add weight for aleatory rupture segmentation
                     wt = wt * al_segWt(iFlt)
                     
c                    Compute Marginal Rate of Occurance
                     if (sourcetype(iFlt) .ne. 7) then
                       mHaz = rate(iParam,iFltWidth) * prock * p1 * probAct(iFlt)                       
                       wt = wt *segwt1(iFLt)
                     elseif (sourcetype(iFlt) .eq. 7) then
                       mHaz = rateS7(iFlt,iMag) * prock * p1 * probAct(iFlt)
                     endif

c                    Add marginal rate of exceed to total
                     Haz(jInten,iProb,iFlt) = Haz(jInten,iProb,iFlt) + mHaz*wt* gm_wt(iProb,jType,iAtten)
                     
                     HazBins(iMagBin,iDistBin,iEpsBin,iProb,jInten) = 
     1                      HazBins(iMagBin,iDistBin,iEpsBin,iProb,jInten) + dble(mHaz*wt)

c  Note: directivity deaggregation was removed from this version
c                     HazBinsX(iXcost,iProb,jInten) = HazBinsX(iXcost,iProb,jInten) + dble(mHaz*wt)
     
c                    Add to mean deagg 
                     wt1 = wt * gm_wt(iProb,jType,iAtten)
                     m_bar(iProb,jInten) = m_bar(iProb,jInten) + mHaz*wt1*magTotal
                     d_bar(iProb,jInten) = d_bar(iProb,jInten) + mHaz*wt1*distRup
                     e_bar(iProb,jInten) = e_bar(iProb,jInten) + mHaz*wt1*epsilon1
                     Xcost_bar(iProb,jInten) = Xcost_bar(iProb,jInten) + mHaz*wt1*Xcost

c
c                    Add to source deagg 
                     wt1 = wt * gm_wt(iProb,jType,iAtten)
                     m_bar_s(iFlt,iProb,jInten) = m_bar_s(iFlt,iProb,jInten) + mHaz*wt1*magTotal
                     rrup_bar_s(iFlt,iProb,jInten) = rrup_bar_s(iFlt,iProb,jInten) + mHaz*wt1*distRup
                     rjb_bar_s(iFlt,iProb,jInten) = rjb_bar_s(iFlt,iProb,jInten) + mHaz*wt1*distjb
                     rx_bar_s(iFlt,iProb,jInten) = rx_bar_s(iFlt,iProb,jInten) + mHaz*wt1*Rx
                     e_bar_s(iFlt,iProb,jInten) = e_bar_s(iFlt,iProb,jInten) + mHaz*wt1*epsilon1
c                     write (*,*) m_bar_s(iFlt,iProb,jInten),rrup_bar_s(iFlt,iProb,jInten), rjb_bar_s(iFlt,iProb,jInten), 
c     1                           rx_bar_s(iFlt,iProb,jInten), e_bar_s(iFlt,iProb,jInten) 

c                    Save Marginal Hazard to temp array for fractile output
                     tempHaz(iParam,jInten,iProb,iAtten,iFtype) = mHaz
     1                        + tempHaz(iParam,jInten,iProb,iAtten,iFtype)

                     tempHaz1(iParam,jInten,iProb,iFtype) = mHaz* gm_wt(iProb,jType,iAtten)
     1                        + tempHaz1(iParam,jInten,iProb,iFtype)

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
          enddo
          mag = minMag(iFlt) + (iMag-0.5) * magStep(iFlt)

c         Write out magnitude rates of occurrence (out2 file)
c         Not currently implemented for SourceType 7 cases.
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
 
 850     MinRrup(iFlt) = MinRrup_temp

c        Write temp Haz array to file
         call WriteTempHaz ( tempHaz, nParamVar, nInten, nProb, 
     1        nAtten, iFlt, attenType(iFlt), nFtype, iFltWidth, nWidth )
     
         call WriteTempHaz1 ( tempHaz1, nParamVar, nInten, nProb, 
     1        nAtten, iFlt, attenType(iFlt), nFtype, iFltWidth, nWidth )

 860    continue

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
     4       epsBins, al_segWt, MinRrup, nAttenType, attenType,
     5       segwt1, dirflag, tapflag,intflag, fsys, SourceDist,
     6       mMagout, hwflagout, ftype, vs, nMaxmag2, mmagoutWt, specT)

c      Write out the deagrregated hazard
        call output_HazBins ( isite, sitex, sitey, testInten, nInten,
     1       nProb, HazBins, jCalc, sigTrunc, csrflag,
     1       nMagBins, nDistBins,
     2       nEpsBins, magBins, distBins, epsBins,
     3       attenName, period1, m_bar, d_bar, e_bar,
     4       nAttenType, attenType, Xcost_bar, nXcostBins, XcostBins,
     5       HazBinsX)
     

       call output_sourcedeagg ( isite, sitex, sitey, testInten, nInten, 
     1           nFlt, Haz, fName, m_bar_s, rrup_bar_s, rjb_bar_s, 
     2           rx_bar_s, e_bar_s, specT, nProb)
     

        call WriteTempHaz2 ( tempHaz2, nInten, nProb, nAtten, nattenType )

 1000 continue

 2000 continue
      close (77)

      write (*,'( 2x,''Normal termination'')')
      stop
      
 2100 write (*,'( 2x,''input file error: number of sites'')')
      stop 99
 2101 write (*,'( 2x,''input file error: site long lat line'')')
      stop 99
 2102 write (*,'( 2x,''input file error: output file name'')')
      stop 99
 2103 write (*,'( 2x,''input file error: output file name'')')
      stop 99
 2104 write (*,'( 2x,''input file error: output file name'')')
      stop 99
 2105 write (*,'( 2x,''input file error: output file name'')')
      stop 99
 2106 write (*,'( 2x,''input file error: output file name'')')
      stop 99
 2107 write (*,'( 2x,''input file error: output file name'')')
      stop 99 
      
      end
