      subroutine S29_Rd_Fault_Data ( nFlt, fName, minMag, magStep, hxStep,
     1     hyStep, segModelWt, rateParam, rateParamWt, beta, 
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, maxMag, 
     3     maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, ftype, 
     4     sourceType, nRupArea, coef_area, sigArea, nRupWidth, 
     5     coef_width, sigWidth, nParamVar, iCoor, minDepth, 
     6     fIndex, probAct, nWidth, mpdf_param, al_segWt, attenType, 
     7     sampleStep, grid_a, grid_dlong, grid_dlat, grid_n, 
     8     grid_long, grid_lat, grid_top, minlat, maxlat, minlong, 
     9     maxlong, scaleRate, fsys, mMagout, mMagoutWt, fltDirect, 
     1     synchron, nsyn_Case, synjcalc, synmag, syndistRup, 
     2     syndistJB, synDistSeismo, synHypo, synftype, synhwflag, 
     3     synwt, RateType, iDepthModel, depthParam, nMaxMag2, segwt1, 
     4     faultFlag, nDownDip, nFtype, ftype_wt, 
     5     segModelFlag, nSegModel0, segModelWt1, syn_dip, 
     6     syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, magS7, rateS7, 
     7     DistS7, DipS7, mechS7, ncountS7, version )

      implicit none
      include 'pfrisk.h'

      integer synHWFlag(MAX_FLT,MAX_SYN), nsyn_Case(MAX_FLT), nFlt, 
     1        synjcalc(MAX_FLT), fltDirect(MAX_FLT), synchron(MAX_FLT),      
     2        fIndex(3,MAX_FLT), nWidth(MAX_FLT), sourceType(MAX_FLT), 
     3        attenType(MAX_FLT), grid_n(MAX_FLT), nMaxMag2(MAX_FLT),
     4        nfp(MAX_FLT), nDownDip(MAX_FLT), nMag(MAX_FLT), iCoor 
      integer nRupArea(MAX_FLT), nRupWidth(MAX_FLT), ncountS7(MAX_FLT),
     1        nParamVar(MAX_FLT,MAX_WIDTH), iDepthModel(MAX_FLT), 
     2        nFtype(MAX_FLT), faultFlag(MAX_FLT,100,MAX_FLT), 
     3        fsys(MAX_FLT),
     4        segModelFlag(MAX_FLT,100), nSegModel0(MAX_FLT) 
      real synmag(MAX_FLT,MAX_SYN), syndistRup(MAX_FLT,MAX_SYN),
     1     syndistJB(MAX_FLT,MAX_SYN), syndistSeismo(MAX_FLT,MAX_SYN), 
     2     synftype(MAX_FLT,MAX_SYN), synhypo(MAX_FLT,MAX_SYN), 
     3     synwt(MAX_FLT,MAX_SYN), syn_dip(MAX_FLT,MAX_SYN), 
     4     syn_zTOR(MAX_FLT,MAX_SYN), syn_RupWidth(MAX_FLT,MAX_SYN)
      real syn_RX(MAX_FLT,MAX_SYN), syn_Ry0(MAX_FLT,MAX_SYN),
     1     al_segWt(MAX_FLT), probAct(MAX_FLT), sigWidth(MAX_FLT),
     2     minlat, maxlat, minlong, maxlong, grid_a(MAX_FLT,MAX_GRID),
     3     grid_lat(MAX_FLT,MAX_GRID), grid_long(MAX_FLT,MAX_GRID), 
     4     grid_top(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT), version 
      real grid_dlat(MAX_FLT), minMag(MAX_FLT), magStep(MAX_FLT), 
     1     hxStep(MAX_FLT), hyStep(MAX_FLT), minDepth(MAX_FLT), 
     2     segModelWt(MAX_FLT), sampleStep(MAX_FLT), scaleRate(MAX_FLT), 
     3     rateParam(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_width(2,MAX_FLT),
     4     rateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH), segWt1(MAX_FLT)    
      real beta(MAX_FLT,MAXPARAM,MAX_WIDTH), ftype_wt(MAX_FLT,MAX_N1),
     1     magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH), magS7(MAX_FLT,MAX_S7), 
     2     magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH), mechS7(MAX_FLT,MAX_S7), 
     3     faultWidth(MAX_FLT,MAX_WIDTH), segModelWt1(MAX_FLT,100), 
     4     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6), DipS7(MAX_FLT,MAX_S7)
      real maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), rateS7(MAX_FLT,MAX_S7), 
     1     maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_area(2,MAX_FLT), 
     2     fLong(MAX_FLT,MAX_DD,MAX_SEG), fZ(MAX_FLT,MAX_DD,MAX_SEG), 
     3     fLat(MAX_FLT,MAX_DD,MAX_SEG), dip(MAX_FLT,MAX_WIDTH,MAX_SEG), 
     4     faultWidthWt(MAX_FLT,MAX_WIDTH), ftype(MAX_FLT,MAX_N1)
      real mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM), distS7(MAX_FLT,MAX_S7),  
     1     depthParam(MAX_FLT,5), mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM), 
     2     sigArea(MAX_FLT), RateType(MAX_FLT,MAXPARAM,MAX_WIDTH) 
      character*80 fName(MAX_FLT)
      
      if ( version .eq. 45.1 ) then
        call S29_Rd_Fault_Data_45_1 ( nFlt, fName, minMag, magStep, hxStep,
     1     hyStep, segModelWt, rateParam, rateParamWt, beta, 
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, maxMag, 
     3     maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, ftype, 
     4     sourceType, nRupArea, coef_area, sigArea, nRupWidth, 
     5     coef_width, sigWidth, nParamVar, iCoor, minDepth, 
     6     fIndex, probAct, nWidth, mpdf_param, al_segWt, attenType, 
     7     sampleStep, grid_a, grid_dlong, grid_dlat, grid_n, 
     8     grid_long, grid_lat, grid_top, minlat, maxlat, minlong, 
     9     maxlong, scaleRate, fsys, mMagout, mMagoutWt, fltDirect, 
     1     synchron, nsyn_Case, synjcalc, synmag, syndistRup, 
     2     syndistJB, synDistSeismo, synHypo, synftype, synhwflag, 
     3     synwt, RateType, iDepthModel, depthParam, nMaxMag2, segwt1, 
     4     faultFlag, nDownDip, nFtype, ftype_wt, 
     5     segModelFlag, nSegModel0, segModelWt1, syn_dip, 
     6     syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, magS7, rateS7, 
     7     DistS7, DipS7, mechS7, ncountS7 )
      elseif (version .eq. 45.2 .or. version .eq. 45.3 ) then
        call S29_Rd_Fault_Data_45_2 ( nFlt, fName, minMag, magStep, hxStep,
     1     hyStep, segModelWt, rateParam, rateParamWt, beta, 
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, maxMag, 
     3     maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, ftype, 
     4     sourceType, nRupArea, coef_area, sigArea, nRupWidth, 
     5     coef_width, sigWidth, nParamVar, iCoor, minDepth, 
     6     fIndex, probAct, nWidth, mpdf_param, al_segWt, attenType, 
     7     sampleStep, grid_a, grid_dlong, grid_dlat, grid_n, 
     8     grid_long, grid_lat, grid_top, minlat, maxlat, minlong, 
     9     maxlong, scaleRate, fsys, mMagout, mMagoutWt, fltDirect, 
     1     synchron, nsyn_Case, synjcalc, synmag, syndistRup, 
     2     syndistJB, synDistSeismo, synHypo, synftype, synhwflag, 
     3     synwt, RateType, iDepthModel, depthParam, nMaxMag2, segwt1, 
     4     faultFlag, nDownDip, nFtype, ftype_wt, 
     5     segModelFlag, nSegModel0, segModelWt1, syn_dip, 
     6     syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, magS7, rateS7, 
     7     DistS7, DipS7, mechS7, ncountS7 )
      else
        write (*,'( 2x,''Incompatible fault file, use Haz45.3, Haz45.2, or Haz45.1'')')
        stop 99
      endif
      
      return
      end

c  --------------------------------------------------------------------

      subroutine S29_Rd_Fault_Data_45_2 ( nFlt, fName, minMag, magStep, hxStep,
     1     hyStep, segModelWt, rateParam, rateParamWt, beta, 
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, maxMag, 
     3     maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, ftype, 
     4     sourceType, nRupArea, coef_area, sigArea, nRupWidth, 
     5     coef_width, sigWidth, nParamVar, iCoor, minDepth, 
     6     fIndex, probAct, nWidth, mpdf_param, al_segWt, attenType, 
     7     sampleStep, grid_a, grid_dlong, grid_dlat, grid_n, 
     8     grid_long, grid_lat, grid_top, minlat, maxlat, minlong, 
     9     maxlong, scaleRate, fsys, mMagout, mMagoutWt, fltDirect, 
     1     synchron, nsyn_Case, synjcalc, synmag, syndistRup, 
     2     syndistJB, synDistSeismo, synHypo, synftype, synhwflag, 
     3     synwt, RateType, iDepthModel, depthParam, nMaxMag2, segwt1, 
     4     faultFlag, nDownDip, nFtype, ftype_wt, 
     5     segModelFlag, nSegModel0, segModelWt1, syn_dip, 
     6     syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, magS7, rateS7, 
     7     DistS7, DipS7, mechS7, ncountS7 )

      implicit none
      include 'pfrisk.h'

      integer synHWFlag(MAX_FLT,MAX_SYN), nsyn_Case(MAX_FLT), 
     1        synjcalc(MAX_FLT), fltDirect(MAX_FLT), synchron(MAX_FLT),      
     2        fIndex(3,MAX_FLT), nWidth(MAX_FLT), sourceType(MAX_FLT), 
     3        attenType(MAX_FLT), grid_n(MAX_FLT), nMaxMag2(MAX_FLT),
     4        nfp(MAX_FLT), nDownDip(MAX_FLT), HWsource7, iRate, nb1 
      integer nMag(MAX_FLT), nRupArea(MAX_FLT), nRupWidth(MAX_FLT),
     1        nRate, n_bValue, nRefMag(MAX_N1), ifsystem, isyn, 
     2        nParamVar(MAX_FLT,MAX_WIDTH), nfsystem,  
     3        fsys(MAX_FLT), iDepthModel(MAX_FLT),
     4        nFtype(MAX_FLT), faultFlag(MAX_FLT,100,MAX_FLT) 
      integer nFtype1(MAX_FLT),
     2        iDip, iWidth, nThick1, nSR, nMoRate, nRecInt, ii, ipt, 
     3        nFlt, iCoor, iFlt0, k, nFlt2, i, iflt2, igrid, n_Dip, 
     4        nActRate, iRecur, iThickDip, iThick1, nRefMag0, iFM 
      integer iflt, nSegModel, nMagRecur, nFtypeModels, 
     1        nFM, iRefMag, i_bValue, segModelFlag(MAX_FLT,100), 
     2        nSegModel0(MAX_FLT), ncountS7(MAX_FLT)
      real synmag(MAX_FLT,MAX_SYN), syndistRup(MAX_FLT,MAX_SYN),
     1     syndistJB(MAX_FLT,MAX_SYN), syndistSeismo(MAX_FLT,MAX_SYN), 
     2     synftype(MAX_FLT,MAX_SYN), synhypo(MAX_FLT,MAX_SYN), 
     3     synwt(MAX_FLT,MAX_SYN), syn_dip(MAX_FLT,MAX_SYN), 
     4     syn_zTOR(MAX_FLT,MAX_SYN), syn_RupWidth(MAX_FLT,MAX_SYN)
      real syn_RX(MAX_FLT,MAX_SYN), syn_Ry0(MAX_FLT,MAX_SYN),
     1     rateParam1(MAX_N1), RateWt1(MAX_N1), al_segWt(MAX_FLT), 
     2     dipWt1(MAX_N1), deltaDip1(MAX_N1), bValue1(MAX_N1), 
     3     bValueWt1(MAX_N1), bValue2(MAX_N1), faultThick1(MAX_N1), 
     4     faultThickWt1(MAX_N1), refMag1(MAX_N1,MAX_N1)
      real refMagWt1(MAX_N1,MAX_N1), refMag0(MAX_N1), refMagWt0(MAX_N1),
     1     magRecurWt1(MAX_N1), magRecur1(MAX_N1), probAct(MAX_FLT),
     2     minlat, maxlat, minlong, maxlong, grid_a(MAX_FLT,MAX_GRID),
     3     grid_lat(MAX_FLT,MAX_GRID), grid_long(MAX_FLT,MAX_GRID), 
     4     grid_top(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT) 
      real grid_dlat(MAX_FLT), minMag(MAX_FLT), magStep(MAX_FLT), 
     1     hxStep(MAX_FLT), hyStep(MAX_FLT), minDepth(MAX_FLT), 
     2     segModelWt(MAX_FLT), sampleStep(MAX_FLT), rP1(MAXPARAM), 
     3     rateParam(MAX_FLT,MAXPARAM,MAX_WIDTH), rP2(MAXPARAM),
     4     rateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH), rp3(MAXPARAM)
      real beta(MAX_FLT,MAXPARAM,MAX_WIDTH), rp4(MAXPARAM), mtest,
     1     magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH), rp5(MAXPARAM),
     2     magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH), sigWidth(MAX_FLT),
     3     faultWidth(MAX_FLT,MAX_WIDTH), ftype_wt(MAX_FLT,MAX_N1),
     4     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6), rp6(MAXPARAM)
      real maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_width(2,MAX_FLT),
     1     maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_area(2,MAX_FLT), 
     2     sigArea(MAX_FLT), fLong(MAX_FLT,MAX_DD,MAX_SEG), 
     3     fLat(MAX_FLT,MAX_DD,MAX_SEG), fZ(MAX_FLT,MAX_DD,MAX_SEG), 
     4     ftype(MAX_FLT,MAX_N1), faultWidthWt(MAX_FLT,MAX_WIDTH)
      real ftype1(MAX_FLT,MAX_N1), ftype_wt1(MAX_FLT,MAX_N1), 
     1     ftmodelwt(MAX_N1), wt_srBranch, wt_ActrateBranch,
     2     mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM), dip1, top, 
     3     scaleRate(MAX_FLT), segWt(MAX_FLT,MAX_FLT)    
      real mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM), wt_recIntBranch,
     1     sr(MAXPARAM), wt_sr(MAXPARAM), actRate(MAXPARAM), 
     2     actRateWt1(MAXPARAM), MoRate(MAXPARAM), wt_MoRate(MAXPARAM),
     3     MoRateDepth(MAXPARAM), MoRDepth(MAXPARAM), segWt1(MAX_FLT),
     4     rec_Int(MAXPARAM), wt_recInt(MAXPARAM), dip2, testMaxMag    
      real dip(MAX_FLT,MAX_WIDTH, MAX_SEG), segModelWt1(MAX_FLT,100), 
     1     wt_MoRateBranch, sum, ProbAct0, depthParam(MAX_FLT,5), 
     2     rateType1(MAXPARAM), RateType(MAX_FLT,MAXPARAM,MAX_WIDTH), 
     3     magS7(MAX_FLT,MAX_S7), rateS7(MAX_FLT,MAX_S7), 
     4     distS7(MAX_FLT,MAX_S7), DipS7(MAX_FLT,MAX_S7) 
      real mechS7(MAX_FLT,MAX_S7)
      character*80 fName(MAX_FLT), fName1

c     Input Fault Parameters
      read (10,*,err=3001) iCoor
      read (10,*,err=3002) NFLT

      iflt = 0
      ifsystem = 1
      
      DO iFlt0=1,NFLT
       read (10,'( a80)',err=3003) fName1
       read (10,*,err=3004) probAct0

c      Read number of segmentation models for this fault system
       read (10,*,err=3005) nSegModel
       read (10,*,err=3006) (segWt(iFlt0,k),k=1,nSegModel)

c      Read total number of fault segments
       read (10,*,err=3007) nFlt2
      
       do i=1,nSegModel
         read (10,*,err=3008) (faultFlag(iFlt0,i,k),k=1,nFlt2)
       enddo

       do iflt2=1,nflt2
        iFlt = iFlt + 1
        call S21_CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c       Set fault indexes
        fIndex(1,iflt) = iflt0
        fIndex(2,iflt) = iflt2

c       Find total weight for this segment
        sum = 0.
        do i=1,nSegModel
          sum = sum + segWt(iFlt0,i) * faultFlag(iFlt0,i,iFlt2)
        enddo
        segWt1(iFlt) = sum

c       Set segment flags for this fault system
        sum = 0.
        do i=1,nSegModel
          segModelFlag(iFlt,i) = faultFlag(iFlt0,i,iFlt2)
          segModelWt1(iFlt,i) = segWt(iFlt0,i)
        enddo
        nSegModel0(iFlt) = nSegModel

c       Read name of this segment
        read(10,'( a80)',err=3009) fName(iFlt)
        write (*,'( i5,2x,a80)') iFlt, fName(iFlt)

c       Read type of source (planar, areal, grid1, grid2, irregular)
        read(10,*,err=3010) sourceType(iFlt),attenType(iFlt),sampleStep(iFlt),
     1                 fltDirect(iFlt), synchron(iFlT)

c       Now read in the synchronous rupture case parameters if needed.
        if (synchron(iFlt) .gt. 0) then
           read (10,*,err=3011) nsyn_Case(iFlt), synjcalc(iFlt)

c          Now read in the magnitude, distance and weigths for each synchronous rupture case.
           do isyn=1,nsyn_Case(iFlt)
             read (10,*,err=3012) synmag(iFlt,isyn),syndistRup(iFlt,isyn),
     1                 syndistJB(iFlt,isyn),syndistSeismo(iFlt,isyn),
     2                 synHWFlag(iFlt,isyn),synhypo(iFlt,isyn),
     3                 synftype(iFlt,isyn), 
     4                 syn_dip(iFlt,isyn), syn_zTOR(iFlt,isyn), 
     5                 syn_RupWidth(iFlt,isyn), syn_RX(iFlt,isyn), syn_Ry0(iFlt,isyn),
     5                 synwt(iFlt,isyn)

           enddo
        endif

c       Read aleatory segmentation wts
        read (10,*,err=3013) al_segWt(iFlt)

c       Check for standard fault source or areal source
        if ( sourceType(iFlt) .eq. 1 .or. sourceType(iFLt) .eq. 2) then
          read (10,*,err=3014) dip1, top
          read(10,*,err=3015) nfp(iFlt)

          call S21_CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
          do ipt=1,nfp(iFlt)
            read (10,*,err=3016) fLong(iFlt,1,ipt), fLat(iFlt,1,ipt)
            fZ(iFlt,1,ipt) = top
            if (sourceType(iFlt) .eq. 2) then
               grid_top(iFlt,ipt) = top
            endif
          enddo
          nDownDip(iFlt) = 1
        endif

c        Check for grid source (w/o depth)
         if ( sourceType(iFlt) .eq. 3 ) then
           read (10,*,err=3014)  dip1, top

           call S30_RdGrid1 (iFlt, grid_a, grid_dlong, grid_dlat, grid_n,
     1             grid_long, grid_lat, minLat, minLong,  maxLat, maxLong, scaleRate(iFlt) )

           do igrid=1,grid_n(iFlt)
             grid_top(iFlt,igrid) = top
           enddo
           nDownDip(iFlt) = 1
         endif

c        Check for grid source (w/ depth)
         if ( sourceType(iFlt) .eq. 4 ) then
              read (10,*,err=3014)  dip1

              call S30_RdGrid2 (iFlt,grid_a,grid_dlong,grid_dlat,grid_n,
     1             grid_long, grid_lat, minLat, minLong, maxLat, maxLong, scaleRate(iFlt),
     2             grid_top)

              nDownDip(iFlt) = 1
         endif

c        Check for custom fault source
         if ( sourceType(iFlt) .eq. 5 .or. sourceType(iFlt) .eq. 6) then
           read(10,*,err=3017) nDownDip(iFLt), nfp(iFlt)  
   
           call S21_CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
           do ipt=1,nfp(iFlt)
              read (10,*,err=3018) (fLong(iFlt,k,ipt), fLat(iFlt,k,ipt), fZ(iflt,k,ipt), k=1,nDownDip(iFlt) ) 

          enddo
         endif

c        Check for Sourcetype 7 (UCERF)
         if ( sourceType(iFlt) .eq. 7 ) then
            read (10,*)  top, HWsource7
            call S30_RdSource7 (iFlt, mags7, rates7, dists7, dips7, mechs7, ncountS7 )
            write (*,*) 'NcountS7 =', ncountS7(iFlt)
         endif

c        Read dip Variation
         if ( sourceType(iFlt) .lt. 5 ) then
           read (10,*,err=3019) n_Dip
           call S21_CheckDim ( n_Dip, MAX_N1, 'MAX_N1    ' )
           read (10,*,err=3020) (deltaDip1(i),i=1,n_Dip)
           read (10,*,err=3021) (dipWt1(i),i=1,n_Dip)
           call S21_CheckWt ( dipWt1, n_Dip, fName(iFlt), 'Dip                 ' )
         else
           n_Dip = 1
           deltaDip1(1) = 0.
           dipWt1(1) = 1.
         endif

c        Read b-values (not for activity rate cases)
         read (10,*,err=3022) n_bValue

         call S21_CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
         if ( n_bValue .gt. 0 ) then
           read (10,*,err=3023) (bValue1(i),i=1,n_bValue)
           read (10,*,err=3024) (bValueWt1(i),i=1,n_bValue)
         endif

c        Read activity rate - b-value pairs
         read (10,*,err=3025) nActRate 

         if ( nActRate .ne. 0 ) then
           do ii=1,nActRate
             read (10,*,err=3026) bValue2(ii), actRate(ii), actRateWt1(ii)
C     Scale activity rate for gridded source (i.e., sourcetype 3 or 4) based on limited lat and long values.
             if (sourcetype(iFlt) .eq. 3 .or. sourcetype(iFlt) .eq. 4) then
                if (scalerate(iFlt) .ne. 1.0) then
                   actRate(ii) = actRate(ii)*Scalerate(iFlt)
                endif
             endif
           enddo

           call S21_CheckWt ( actRateWt1, nActRate, fName(iFlt), 'ActRate             ' )
         endif

c        Read weights for rate methods
         read (10,*,err=3027) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
         sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
         if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFlt)
              stop 99
         endif

c        Read slip-rates
         read (10,*,err=3028) nSR
         if ( nSR .gt. 0 ) then
           read (10,*,err=3029) (sr(k),k=1,nSR)
           read (10,*,err=3030) (wt_sr(k),k=1,nSR)
           call S21_CheckWt (wt_sr, nSR, fName(iFlt), 'Slip Rates          ')
         endif

c        Read recurrence intervals
         read (10,*,err=3031) nRecInt
         if ( nRecInt .gt. 0 ) then
           read (10,*,err=3032) (rec_Int(k),k=1,nRecInt)
           read (10,*,err=3033) (wt_recInt(k),k=1,nRecInt)
           call S21_CheckWt (wt_recInt, nRecInt, fName(iFlt), 'Recurrence Intervals')
         endif

c        Read moment-rates
         read (10,*,err=3034) nMoRate
         if ( nMoRate .gt. 0 ) then
           read (10,*,err=3035) (MoRate(k),k=1,nMoRate)
           read (10,*,err=3036) (MoRateDepth(k),k=1,nMoRate)
           read (10,*,err=3037) (wt_MoRate(k),k=1,nMoRate)
           call S21_CheckWt (wt_MoRate, nMoRate, fName(iFlt), 'Moment Rates        ')
         endif

c        Check that slip-rates are not used for areal sources
         if ( sourceType(iFlt) .ge. 2 .and. sourceType(iFlt) .le. 4 ) then               
           if (nSR .ne. 0) then
              write (*,'( 2x,''Error: slip-rates not allowed for areal sources'')')
              write (*,'( 2x,''Flt: '',a60)') fname(iFlt)
              stop 99
           endif
         endif

c        Load into single array called "rate_param"
c        Set MoRDepth=1.0 or the inverse for latter scaling of MoRates
         nRate = nSR + nActRate + nRecInt + nMoRate
         call S21_CheckDim ( nRate, MAX_N1, 'MAX_N1    ' )
         do k=1,nSR
            rateParam1(k) = sr(k)
            rateWt1(k) = wt_sr(k)*wt_srbranch              
            rateType1(k) = 1
            MoRDepth(k) = 1.0
         enddo
         do k=1,nActRate
            rateParam1(k+nSR) = actRate(k)
            rateWt1(k+nSR) = actRateWt1(k) * wt_actRateBranch             
            rateType1(k+nSR) = 2
            MoRDepth(k+nSR) = 1.0
         enddo
         do k=1,nRecInt
            rateParam1(k+nSR+nActRate) = rec_Int(k)
            rateWt1(k+nSR+nActRate) = wt_recInt(k) *wt_recIntBranch              
            rateType1(k+nSR+nActRate) = 3
            MoRDepth(k+nSR+nActRate) = 1.0
         enddo
         do k=1,nMoRate
            rateParam1(k+nSR+nActRate+nRecInt) = MoRate(k)
            rateWt1(k+nSR+nActRate+nRecInt) = wt_MoRate(k) *wt_MoRateBranch              
            rateType1(k+nSR+nActRate+nRecInt) = 4
            MoRDepth(k+nSR+nActRate+nRecInt) = 1.0/MoRateDepth(nMoRate)
         enddo

c        Read Mag recurrence weights (char, exp, etc.)
         read (10,*,err=3038) nMagRecur
         call S21_CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
         read (10,*,err=3039) (magRecur1(i),i=1,nMagRecur)
         read (10,*,err=3040) (magRecurWt1(i),i=1,nMagRecur)
         call S21_CheckWt ( magRecurWt1,nMagRecur, fName(iFlt), 'Mag Recur           ' )

c        Read in corresponding magnitude parameters. 
         do iRecur=1,nMagRecur
           if (magRecur1(iRecur) .eq. 4) then
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
             rp6(iRecur) = 0.0
c          Read in necessary parameters for bi-exponential distribution.
           elseif (magRecur1(iRecur) .eq. 5) then
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
             rp6(iRecur) = 0.0
           elseif (magRecur1(iRecur) .eq. 10) then
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur),
     1       rp6(iRecur)
           else              
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur)
             rp4(iRecur) = 0.0
             rp5(iRecur) = 0.0
             rp6(iRecur) = 0.0
           endif
         enddo

c        Read seismogenic thickness
         if ( sourceType(iFlt) .lt. 5) then
           read (10,*,err=3042) nThick1
           call S21_CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
           read (10,*,err=3043) (faultThick1(i),i=1,nThick1)
           read (10,*,err=3044) (faultThickWt1(i),i=1,nThick1)         
           call S21_CheckWt (faultThickWt1, nThick1, fName(iFlt), 'Seismo Thick        ')
         else
           nThick1 = 1
           faultTHick1(1) = -99.
           faultThickWt1(1) = 1.
         endif
         
c        Read depth pdf
         read (10,*,err=3045) iDepthModel(iFlt), (depthParam(iflt,k),k=1,3)     

c        Read reference mags for each fault thickness
         iThickDip = 1
         do iThick1=1,nThick1
           read (10,*,err=3047) nRefMag0
           read (10,*,err=3048) (refMag0(i),i=1,nRefMag0)
           read (10,*,err=3049) (refMagWt0(i),i=1,nRefMag0)

           if ( nRefMag0 .ne. 0 ) then
             call S21_CheckWt ( refMagWt0, nrefMag0, fName(iFlt), 'Ref Mag                 ' )
           endif
              
c          Copy these ref magnitudes for each addtional dip (no correlation allowed)
           do iDip=1,n_Dip
             nRefMag(iThickDip) = nRefMag0
             do i=1,nRefMag0
               refMag1(iThickDip,i) = refMag0(i)
               refMagWt1(iThickDip,i) = refMagWt0(i)
             enddo
             iThickDip = iThickDip + 1
           enddo
         enddo

c        Read min mag, step sizes, and rupture variability info
         read (10,*,err=3050) minMag(iFlt), magStep(iFlt), hxStep(iFlt), 
     1              hyStep(iFlt), nRupArea(iFlt), nRupWidth(iFlt), minDepth(iFlt)

C       Allow hxstep and hystep to be different for sourceType 2
        if (sourceType(iFlt).ne.2) then
C         Check that Hxstep = Hystep
          if (hxstep(iFlt) .ne. hystep(iFlt) ) then
            write (*,*) 'Hxstep not equal to Hystep for fault: '
            write (*,*) fname1
            write (*,*) 'Hxstep = ', hxstep(iFlt) 
            write (*,*) 'Hystep = ', hystep(iFlt) 
            write (*,*) 'These values must be equal or errors can occur!'
            write (*,*) 'Check input fault file.'
            stop 99
          endif
        endif

         read (10,*,err=3051) (coef_area(k,iFlt),k=1,2), sigArea(iFlt)
         read (10,*,err=3052) (coef_width(k,iFlt),k=1,2), sigWidth(iFlt)

c        Read ftype Models
         read (10,*,err=3053) nFtypeModels
         do iFM=1,nFtypeModels
           read (10,*,err=3054) ftmodelwt(iFM)
           read (10,*,err=3055) nFtype1(iFM)
           read (10,*,err=3056) (ftype1(iFM,k),k=1,nFtype1(iFM))
           read (10,*,err=3057) (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
           call S21_CheckWt1a (Ftype_Wt1, nFtype1(iFM), iFM, 
     1          MAX_N1, fName(iFlt), 'Fault Mech          ' )
         enddo

c        Load up Ftype Models and weights.
         nFM = 1
         do iFM=1,nFtypeModels
           do k=1,nFtype1(iFM)
             ftype(iFlt,nFM) = ftype1(iFM,k)
             ftype_wt(iFlt,nFM) = ftype_wt1(iFM,k)*ftmodelwt(iFM)
             nFm = nFm + 1
           enddo
         enddo
         nFm = nFm - 1
         nFtype(iFlt) = nFm
         call S21_CheckDim ( nFtype(iFlt), MAX_FTYPE, 'MAX_FTYPE ' )

c        Load up parameter variations into large single dimension arrays
         testMaxMag = 0.
         iWidth = 0
         do iThick1=1,nThick1
        
           do iDip=1,n_Dip
           iWidth = iWidth + 1
           call S21_CheckDim ( iWidth, MAX_WIDTH, 'MAX_WIDTH' )   
           
             dip2 = dip1 + deltaDip1(iDip)
             dip(iFlt,iWidth,1) = dip2
             faultWidth(iFlt,iWidth) = faultThick1(iThick1)
             faultWidthWt(iFlt,iWidth) = faultThickWt1(iThick1) * dipWt1(iDip)

           i = 0
           mtest = 0.0

           do iRecur=1,nMagRecur    
           
            do iRate=1,nRate
              if ( rateType1(iRate) .eq. 2 ) then
                nb1 = 1
              else
                nb1 = n_bvalue
              endif
              do i_bValue=1,nb1
               do iRefMag=1,nRefMag(iWidth)
                 i = i + 1
                 call S21_CheckDim ( i, MAXPARAM, 'MAXPARAM  ' )
                 magRecur(iFlt,i,iWidth) = magRecur1(iRecur)
                 magRecurWt(iFlt,i,iWidth) = magRecurWt1(iRecur)

c                Scale moment rate by reference thickness.
                 if (rateType1(iRate) .eq. 4) then
                    rateParam(iFlt,i,iWidth) = rateParam1(iRate)*MoRDepth(iRate)*faultWidth(iFlt,iWidth)
                 else
                    rateParam(iFlt,i,iWidth) = rateParam1(iRate)
                 endif

                 rateType(iFlt,i,iWidth) = rateType1(iRate)
                 if ( rateType1(iRate) .eq. 2 ) then
                   beta(iFlt,i,iWidth) = bValue2(iRate - nSR)*alog(10.)
                 else
                   beta(iFlt,i,iWidth) = bValue1(i_bValue)*alog(10.)
                 endif
                 if ( rateType1(iRate) .eq. 2 ) then
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate)
                 else
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate) * bValueWt1(i_bValue) 
                 endif
                 maxMagWt(iFlt,i,iWidth) = refMagWt1(iWidth,iRefMag)

c                Set max mag
                 if (magRecur1(iRecur) .eq. 0 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 1 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP1(iRecur)
                 elseif (magRecur1(iRecur) .eq. 3 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 6 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 elseif (magRecur1(iRecur) .eq. 10 ) then
                   if ( rp5(iRecur) .eq. 0. ) then
                     maxMag(iFlt,i,iWidth) = rp2(iRecur)
                   else
                     maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rp5(iRecur)
                   endif
                 endif
                 if ( maxMag(iFlt,i,iWidth) .gt. testMaxMag ) then
                   nMag(iFlt) = ceiling((maxMag(iFlt,i,iWidth) - minMag(iFLt) ) / magStep(iFlt))
                   if (sourceType(iFlt) .eq. 7) then
                     nMag(iFlt) = ncountS7(iFlt)
                   endif
                     if (nMag(iFlt) .eq. 0 ) then
                       nMag(iFlt) = 1
                     endif    
                   testMaxMag = maxMag(iFlt,i,iWidth)
                 endif               

c  temp fix - the maxmag, refmag problem for waacy
                 if (magRecur1(iRecur) .eq. 10 ) then
                     maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 endif
                          
                 mpdf_param(iFlt,i,iWidth,1) = rP1(iRecur)
                 mpdf_param(iFlt,i,iWidth,2) = rP2(iRecur)
                 mpdf_param(iFlt,i,iWidth,3) = rP3(iRecur)
                 mpdf_param(iFlt,i,iWidth,4) = rP4(iRecur)
                 mpdf_param(iFlt,i,iWidth,5) = rP5(iRecur)
                 mpdf_param(iFlt,i,iWidth,6) = rP6(iRecur)

c      Test for the maximum magnitude value for all realizations of this fault.
                 if (maxMag(iFlt,i,iWidth) .ge. mtest) then
                   mtest = maxMag(iFlt,i,iWidth)
                 endif

                 mmagout(iFlt,iWidth,iRefMag) = refMag1(iWidth,iRefMag)
                 mmagoutwt(iFlt,iWidth,iRefMag) = refMagWt1(iWidth,iRefMag)
     
c     End loop over iRefMag
                enddo
c     End loop over ib_value
               enddo
c     End loop over iRate
              enddo
c     End loop over iRecur
             enddo
             nParamVar(iFlt,iWidth) = i
c     End loop over iDip
            enddo
c     End loop over iThick1
           enddo
           probAct(iFlt) = probAct0
           nWidth(iFlt) = iWidth
           fsys(iFlt) = ifsystem
c     End loop over iFlt2 (# segments)
          enddo
          ifsystem = ifsystem + 1          
c     End loop over iFlt
      enddo
      nfsystem = ifsystem - 1
      nFlt = iflt
      
      close (10)
      
      return
 3001 write (*,'( 2x,''Flt file error:  iCorr'')')
      stop 99
 3002 write (*,'( 2x,''Flt file error:  nFlt'')')
      stop 99
 3003 write (*,'( 2x,''Flt file error:  flt sys name'')')
      stop 99
 3004 write (*,'( 2x,''Flt file error:  prob Act'')')
      stop 99
 3005 write (*,'( 2x,''Flt file error:  nSegModel'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3006 write (*,'( 2x,''Flt file error:  Seg wts'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3007 write (*,'( 2x,''Flt file error:  number of segments'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3008  write (*,'( 2x,''Flt file error: seg flags for seg model'', i5)') i
       write (*,'( 2x,''From fault: '',a80)') fName1
       stop 99
 3009 write (*,'( 2x,''Flt file error: seg name for seg model '', i5)') i
      write (*,'( 2x,''From fault: '',a80)') fName1
      stop 99
 3010 write (*,'( 2x,''Flt file error: source type line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3011 write (*,'( 2x,''Flt file error: nSynRup line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3012 write (*,'( 2x,''Flt file error: syn Rup param line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3013 write (*,'( 2x,''Flt file error: Aleatory seg wt '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3014 write (*,'( 2x,''Flt file error: dip and top '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3015 write (*,'( 2x,''Flt file error: number flt point '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3016 write (*,'( 2x,''Flt file error: long lat values '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3017 write (*,'( 2x,''Flt file error: nDowndip, nFP '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3018 write (*,'( 2x,''Flt file error: long lat values '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3019 write (*,'( 2x,''Flt file error: nDip '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3020 write (*,'( 2x,''Flt file error: delta Dip '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3021 write (*,'( 2x,''Flt file error: Dip wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3022 write (*,'( 2x,''Flt file error: n b-value '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3023 write (*,'( 2x,''Flt file error: delta b-value '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3024 write (*,'( 2x,''Flt file error: b-value wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3025 write (*,'( 2x,''Flt file error: nActRate'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3026 write (*,'( 2x,''Flt file error: b, activity, wt'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3027 write (*,'( 2x,''Flt file error: wts for activity rate approach'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3028 write (*,'( 2x,''Flt file error: number Slip rates '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3029 write (*,'( 2x,''Flt file error: slip rates '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3030 write (*,'( 2x,''Flt file error: slip-rate wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3031 write (*,'( 2x,''Flt file error: number Rec Intervals '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3032 write (*,'( 2x,''Flt file error: Rec Intervals '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3033 write (*,'( 2x,''Flt file error: Rec Interval wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3034 write (*,'( 2x,''Flt file error: Number Moment rate'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3035 write (*,'( 2x,''Flt file error: Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3036 write (*,'( 2x,''Flt file error: depths for Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3037 write (*,'( 2x,''Flt file error: wts for Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3038 write (*,'( 2x,''Flt file error: number mag pdf'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3039 write (*,'( 2x,''Flt file error: mag pdf flag'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3040 write (*,'( 2x,''Flt file error: mag pdf wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3041 write (*,'( 2x,''Flt file error: param for mag pdf'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3042 write (*,'( 2x,''Flt file error: number crustal thickness '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3043 write (*,'( 2x,''Flt file error: crustal thicknesses '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3044 write (*,'( 2x,''Flt file error: crustal thickness wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3045 write (*,'( 2x,''Flt file error: depth pdf and param'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3047 write (*,'( 2x,''Flt file error: depth pdf and param'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3048 write (*,'( 2x,''Flt file error: number Ref Mag'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3049 write (*,'( 2x,''Flt file error: Ref mags'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3050 write (*,'( 2x,''Flt file error: Ref mag wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3051 write (*,'( 2x,''Flt file error: coeff A(M) model'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3052 write (*,'( 2x,''Flt file error: coeff W(M) model'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3053 write (*,'( 2x,''Flt file error: number Ftype models'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3054 write (*,'( 2x,''Flt file error: Ftype mode wt'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3055 write (*,'( 2x,''Flt file error: number Ftype'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3056 write (*,'( 2x,''Flt file error: Ftypes'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3057 write (*,'( 2x,''Flt file error: Ftype Aleatory wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99

      end

c ----------------------------

      subroutine S29_Rd_Fault_Data_45_1 ( nFlt, fName, minMag, magStep, hxStep,
     1     hyStep, segModelWt, rateParam, rateParamWt, beta, 
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, maxMag, 
     3     maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, ftype, 
     4     sourceType, nRupArea, coef_area, sigArea, nRupWidth, 
     5     coef_width, sigWidth, nParamVar, iCoor, minDepth, 
     6     fIndex, probAct, nWidth, mpdf_param, al_segWt, attenType, 
     7     sampleStep, grid_a, grid_dlong, grid_dlat, grid_n, 
     8     grid_long, grid_lat, grid_top, minlat, maxlat, minlong, 
     9     maxlong, scaleRate, fsys, mMagout, mMagoutWt, fltDirect, 
     1     synchron, nsyn_Case, synjcalc, synmag, syndistRup, 
     2     syndistJB, synDistSeismo, synHypo, synftype, synhwflag, 
     3     synwt, RateType, iDepthModel, depthParam, nMaxMag2, segwt1, 
     4     faultFlag, nDownDip, nFtype, ftype_wt, 
     5     segModelFlag, nSegModel0, segModelWt1, syn_dip, 
     6     syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, magS7, rateS7, 
     7     DistS7, DipS7, mechS7, ncountS7 )

      implicit none
      include 'pfrisk.h'

      integer synHWFlag(MAX_FLT,MAX_SYN), nsyn_Case(MAX_FLT), 
     1        synjcalc(MAX_FLT), fltDirect(MAX_FLT), synchron(MAX_FLT),      
     2        fIndex(3,MAX_FLT), nWidth(MAX_FLT), sourceType(MAX_FLT), 
     3        attenType(MAX_FLT), grid_n(MAX_FLT), nMaxMag2(MAX_FLT),
     4        nfp(MAX_FLT), nDownDip(MAX_FLT), HWsource7, iRate, nb1 
      integer nMag(MAX_FLT), nRupArea(MAX_FLT), nRupWidth(MAX_FLT),
     1        nRate, n_bValue, nRefMag(MAX_N1), ifsystem, isyn, 
     2        nParamVar(MAX_FLT,MAX_WIDTH), nfsystem,  
     3        fsys(MAX_FLT), iDepthModel(MAX_FLT), iOverRideMag,
     4        nFtype(MAX_FLT), faultFlag(MAX_FLT,100,MAX_FLT) 
      integer nFtype1(MAX_FLT),
     2        iDip, iWidth, nThick1, nSR, nMoRate, nRecInt, ii, ipt, 
     3        nFlt, iCoor, iFlt0, k, nFlt2, i, iflt2, igrid, n_Dip, 
     4        nActRate, iRecur, iThickDip, iThick1, nRefMag0, iFM 
      integer iflt, nSegModel, nMagRecur, nFtypeModels, 
     1        nFM, iRefMag, i_bValue, segModelFlag(MAX_FLT,100), 
     2        nSegModel0(MAX_FLT), ncountS7(MAX_FLT)
      real synmag(MAX_FLT,MAX_SYN), syndistRup(MAX_FLT,MAX_SYN),
     1     syndistJB(MAX_FLT,MAX_SYN), syndistSeismo(MAX_FLT,MAX_SYN), 
     2     synftype(MAX_FLT,MAX_SYN), synhypo(MAX_FLT,MAX_SYN), 
     3     synwt(MAX_FLT,MAX_SYN), syn_dip(MAX_FLT,MAX_SYN), 
     4     syn_zTOR(MAX_FLT,MAX_SYN), syn_RupWidth(MAX_FLT,MAX_SYN)
      real syn_RX(MAX_FLT,MAX_SYN), syn_Ry0(MAX_FLT,MAX_SYN),
     1     rateParam1(MAX_N1), RateWt1(MAX_N1), al_segWt(MAX_FLT), 
     2     dipWt1(MAX_N1), deltaDip1(MAX_N1), bValue1(MAX_N1), 
     3     bValueWt1(MAX_N1), bValue2(MAX_N1), faultThick1(MAX_N1), 
     4     faultThickWt1(MAX_N1), refMag1(MAX_N1,MAX_N1)
      real refMagWt1(MAX_N1,MAX_N1), refMag0(MAX_N1), refMagWt0(MAX_N1),
     1     magRecurWt1(MAX_N1), magRecur1(MAX_N1), probAct(MAX_FLT),
     2     minlat, maxlat, minlong, maxlong, grid_a(MAX_FLT,MAX_GRID),
     3     grid_lat(MAX_FLT,MAX_GRID), grid_long(MAX_FLT,MAX_GRID), 
     4     grid_top(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT) 
      real grid_dlat(MAX_FLT), minMag(MAX_FLT), magStep(MAX_FLT), 
     1     hxStep(MAX_FLT), hyStep(MAX_FLT), minDepth(MAX_FLT), 
     2     segModelWt(MAX_FLT), sampleStep(MAX_FLT), rP1(MAXPARAM), 
     3     rateParam(MAX_FLT,MAXPARAM,MAX_WIDTH), rP2(MAXPARAM),
     4     rateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH), rp3(MAXPARAM)
      real beta(MAX_FLT,MAXPARAM,MAX_WIDTH), rp4(MAXPARAM), mtest,
     1     magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH), rp5(MAXPARAM),
     2     magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH), sigWidth(MAX_FLT),
     3     faultWidth(MAX_FLT,MAX_WIDTH), ftype_wt(MAX_FLT,MAX_N1),
     4     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,6), rp6(MAXPARAM)
      real maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_width(2,MAX_FLT),
     1     maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH), coef_area(2,MAX_FLT), 
     2     sigArea(MAX_FLT), fLong(MAX_FLT,MAX_DD,MAX_SEG), 
     3     fLat(MAX_FLT,MAX_DD,MAX_SEG), fZ(MAX_FLT,MAX_DD,MAX_SEG), 
     4     ftype(MAX_FLT,MAX_N1), faultWidthWt(MAX_FLT,MAX_WIDTH)
      real ftype1(MAX_FLT,MAX_N1), ftype_wt1(MAX_FLT,MAX_N1), 
     1     ftmodelwt(MAX_N1), wt_srBranch, wt_ActrateBranch,
     2     mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM), dip1, top, 
     3     scaleRate(MAX_FLT), segWt(MAX_FLT,MAX_FLT)    
      real mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM), wt_recIntBranch,
     1     sr(MAXPARAM), wt_sr(MAXPARAM), actRate(MAXPARAM), 
     2     actRateWt1(MAXPARAM), MoRate(MAXPARAM), wt_MoRate(MAXPARAM),
     3     MoRateDepth(MAXPARAM), MoRDepth(MAXPARAM), segWt1(MAX_FLT),
     4     rec_Int(MAXPARAM), wt_recInt(MAXPARAM), dip2, testMaxMag    
      real dip(MAX_FLT,MAX_WIDTH, MAX_SEG), segModelWt1(MAX_FLT,100), 
     1     wt_MoRateBranch, sum, ProbAct0, depthParam(MAX_FLT,5), 
     2     rateType1(MAXPARAM), RateType(MAX_FLT,MAXPARAM,MAX_WIDTH), 
     3     magS7(MAX_FLT,MAX_S7), rateS7(MAX_FLT,MAX_S7), 
     4     distS7(MAX_FLT,MAX_S7), DipS7(MAX_FLT,MAX_S7) 
      real mechS7(MAX_FLT,MAX_S7)
      character*80 fName(MAX_FLT), fName1

c     Input Fault Parameters
      read (10,*,err=3001) iCoor
      read (10,*,err=3002) NFLT

      iflt = 0
      ifsystem = 1
      
      DO iFlt0=1,NFLT
       read (10,'( a80)',err=3003) fName1
       read (10,*,err=3004) probAct0

c      Read number of segmentation models for this fault system
       read (10,*,err=3005) nSegModel
       read (10,*,err=3006) (segWt(iFlt0,k),k=1,nSegModel)

c      Read total number of fault segments
       read (10,*,err=3007) nFlt2
      
       do i=1,nSegModel
         read (10,*,err=3008) (faultFlag(iFlt0,i,k),k=1,nFlt2)
       enddo

       do iflt2=1,nflt2
        iFlt = iFlt + 1
        call S21_CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c       Set fault indexes
        fIndex(1,iflt) = iflt0
        fIndex(2,iflt) = iflt2

c       Find total weight for this segment
        sum = 0.
        do i=1,nSegModel
          sum = sum + segWt(iFlt0,i) * faultFlag(iFlt0,i,iFlt2)
        enddo
        segWt1(iFlt) = sum

c       Set segment flags for this fault system
        sum = 0.
        do i=1,nSegModel
          segModelFlag(iFlt,i) = faultFlag(iFlt0,i,iFlt2)
          segModelWt1(iFlt,i) = segWt(iFlt0,i)
        enddo
        nSegModel0(iFlt) = nSegModel

c       Read name of this segment
        read(10,'( a80)',err=3009) fName(iFlt)
        write (*,'( i5,2x,a80)') iFlt, fName(iFlt)

c       Read type of source (planar, areal, grid1, grid2, irregular)
        read(10,*,err=3010) sourceType(iFlt),attenType(iFlt),sampleStep(iFlt),
     1                 fltDirect(iFlt), synchron(iFlT)

c       Now read in the synchronous rupture case parameters if needed.
        if (synchron(iFlt) .gt. 0) then
           read (10,*,err=3011) nsyn_Case(iFlt), synjcalc(iFlt)

c          Now read in the magnitude, distance and weigths for each synchronous rupture case.
           do isyn=1,nsyn_Case(iFlt)
             read (10,*,err=3012) synmag(iFlt,isyn),syndistRup(iFlt,isyn),
     1                 syndistJB(iFlt,isyn),syndistSeismo(iFlt,isyn),
     2                 synHWFlag(iFlt,isyn),synhypo(iFlt,isyn),
     3                 synftype(iFlt,isyn), 
     4                 syn_dip(iFlt,isyn), syn_zTOR(iFlt,isyn), 
     5                 syn_RupWidth(iFlt,isyn), syn_RX(iFlt,isyn), syn_Ry0(iFlt,isyn),
     5                 synwt(iFlt,isyn)

           enddo
        endif

c       Read aleatory segmentation wts
        read (10,*,err=3013) al_segWt(iFlt)

c       Check for standard fault source or areal source
        if ( sourceType(iFlt) .eq. 1 .or. sourceType(iFLt) .eq. 2) then
          read (10,*,err=3014) dip1, top
          read(10,*,err=3015) nfp(iFlt)

          call S21_CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
          do ipt=1,nfp(iFlt)
            read (10,*,err=3016) fLong(iFlt,1,ipt), fLat(iFlt,1,ipt)
            fZ(iFlt,1,ipt) = top
            if (sourceType(iFlt) .eq. 2) then
               grid_top(iFlt,ipt) = top
            endif
          enddo
          nDownDip(iFlt) = 1
        endif

c        Check for grid source (w/o depth)
         if ( sourceType(iFlt) .eq. 3 ) then
           read (10,*,err=3014)  dip1, top

           call S30_RdGrid1 (iFlt, grid_a, grid_dlong, grid_dlat, grid_n,
     1             grid_long, grid_lat, minLat, minLong,  maxLat, maxLong, scaleRate(iFlt) )

           do igrid=1,grid_n(iFlt)
             grid_top(iFlt,igrid) = top
           enddo
           nDownDip(iFlt) = 1
         endif

c        Check for grid source (w/ depth)
         if ( sourceType(iFlt) .eq. 4 ) then
              read (10,*,err=3014)  dip1

              call S30_RdGrid2 (iFlt,grid_a,grid_dlong,grid_dlat,grid_n,
     1             grid_long, grid_lat, minLat, minLong, maxLat, maxLong, scaleRate(iFlt),
     2             grid_top)

              nDownDip(iFlt) = 1
         endif

c        Check for custom fault source
         if ( sourceType(iFlt) .eq. 5 .or. sourceType(iFlt) .eq. 6) then
           read(10,*,err=3017) nDownDip(iFLt), nfp(iFlt)  
   
           call S21_CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
           do ipt=1,nfp(iFlt)
              read (10,*,err=3018) (fLong(iFlt,k,ipt), fLat(iFlt,k,ipt), fZ(iflt,k,ipt), k=1,nDownDip(iFlt) ) 

          enddo
         endif

c        Check for Sourcetype 7 (UCERF)
         if ( sourceType(iFlt) .eq. 7 ) then
            read (10,*)  top, HWsource7
            call S30_RdSource7 (iFlt, mags7, rates7, dists7, dips7, mechs7, ncountS7 )
            write (*,*) 'NcountS7 =', ncountS7(iFlt)
         endif

c        Read dip Variation
         if ( sourceType(iFlt) .lt. 5 ) then
           read (10,*,err=3019) n_Dip
           call S21_CheckDim ( n_Dip, MAX_N1, 'MAX_N1    ' )
           read (10,*,err=3020) (deltaDip1(i),i=1,n_Dip)
           read (10,*,err=3021) (dipWt1(i),i=1,n_Dip)
           call S21_CheckWt ( dipWt1, n_Dip, fName(iFlt), 'Dip                 ' )
         else
           n_Dip = 1
           deltaDip1(1) = 0.
           dipWt1(1) = 1.
         endif

c        Read b-values (not for activity rate cases)
         read (10,*,err=3022) n_bValue

         call S21_CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
         if ( n_bValue .gt. 0 ) then
           read (10,*,err=3023) (bValue1(i),i=1,n_bValue)
           read (10,*,err=3024) (bValueWt1(i),i=1,n_bValue)
         endif

c        Read activity rate - b-value pairs
         read (10,*,err=3025) nActRate 

         if ( nActRate .ne. 0 ) then
           do ii=1,nActRate
             read (10,*,err=3026) bValue2(ii), actRate(ii), actRateWt1(ii)
C     Scale activity rate for gridded source (i.e., sourcetype 3 or 4) based on limited lat and long values.
             if (sourcetype(iFlt) .eq. 3 .or. sourcetype(iFlt) .eq. 4) then
                if (scalerate(iFlt) .ne. 1.0) then
                   actRate(ii) = actRate(ii)*Scalerate(iFlt)
                endif
             endif
           enddo

           call S21_CheckWt ( actRateWt1, nActRate, fName(iFlt), 'ActRate             ' )
         endif

c        Read weights for rate methods
         read (10,*,err=3027) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch
         sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
         if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFlt)
              stop 99
         endif

c        Read slip-rates
         read (10,*,err=3028) nSR
         if ( nSR .gt. 0 ) then
           read (10,*,err=3029) (sr(k),k=1,nSR)
           read (10,*,err=3030) (wt_sr(k),k=1,nSR)
           call S21_CheckWt (wt_sr, nSR, fName(iFlt), 'Slip Rates          ')
         endif

c        Read recurrence intervals
         read (10,*,err=3031) nRecInt
         if ( nRecInt .gt. 0 ) then
           read (10,*,err=3032) (rec_Int(k),k=1,nRecInt)
           read (10,*,err=3033) (wt_recInt(k),k=1,nRecInt)
           call S21_CheckWt (wt_recInt, nRecInt, fName(iFlt), 'Recurrence Intervals')
         endif

c        Read moment-rates
         read (10,*,err=3034) nMoRate
         if ( nMoRate .gt. 0 ) then
           read (10,*,err=3035) (MoRate(k),k=1,nMoRate)
           read (10,*,err=3036) (MoRateDepth(k),k=1,nMoRate)
           read (10,*,err=3037) (wt_MoRate(k),k=1,nMoRate)
           call S21_CheckWt (wt_MoRate, nMoRate, fName(iFlt), 'Moment Rates        ')
         endif

c        Check that slip-rates are not used for areal sources
         if ( sourceType(iFlt) .ge. 2 .and. sourceType(iFlt) .le. 4 ) then               
           if (nSR .ne. 0) then
              write (*,'( 2x,''Error: slip-rates not allowed for areal sources'')')
              write (*,'( 2x,''Flt: '',a60)') fname(iFlt)
              stop 99
           endif
         endif

c        Load into single array called "rate_param"
c        Set MoRDepth=1.0 or the inverse for latter scaling of MoRates
         nRate = nSR + nActRate + nRecInt + nMoRate
         call S21_CheckDim ( nRate, MAX_N1, 'MAX_N1    ' )
         do k=1,nSR
            rateParam1(k) = sr(k)
            rateWt1(k) = wt_sr(k)*wt_srbranch              
            rateType1(k) = 1
            MoRDepth(k) = 1.0
         enddo
         do k=1,nActRate
            rateParam1(k+nSR) = actRate(k)
            rateWt1(k+nSR) = actRateWt1(k) * wt_actRateBranch             
            rateType1(k+nSR) = 2
            MoRDepth(k+nSR) = 1.0
         enddo
         do k=1,nRecInt
            rateParam1(k+nSR+nActRate) = rec_Int(k)
            rateWt1(k+nSR+nActRate) = wt_recInt(k) *wt_recIntBranch              
            rateType1(k+nSR+nActRate) = 3
            MoRDepth(k+nSR+nActRate) = 1.0
         enddo
         do k=1,nMoRate
            rateParam1(k+nSR+nActRate+nRecInt) = MoRate(k)
            rateWt1(k+nSR+nActRate+nRecInt) = wt_MoRate(k) *wt_MoRateBranch              
            rateType1(k+nSR+nActRate+nRecInt) = 4
            MoRDepth(k+nSR+nActRate+nRecInt) = 1.0/MoRateDepth(nMoRate)
         enddo

c        Read Mag recurrence weights (char, exp, etc.)
         read (10,*,err=3038) nMagRecur
         call S21_CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
         read (10,*,err=3039) (magRecur1(i),i=1,nMagRecur)
         read (10,*,err=3040) (magRecurWt1(i),i=1,nMagRecur)
         call S21_CheckWt ( magRecurWt1,nMagRecur, fName(iFlt), 'Mag Recur           ' )

c        Read in corresponding magnitude parameters. 
         do iRecur=1,nMagRecur
           if (magRecur1(iRecur) .eq. 4) then
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
             rp6(iRecur) = 0.0
c          Read in necessary parameters for bi-exponential distribution.
           elseif (magRecur1(iRecur) .eq. 5) then
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
             rp6(iRecur) = 0.0
           elseif (magRecur1(iRecur) .eq. 10) then
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur),
     1       rp6(iRecur)
           else              
             read (10,*,err=3041) rP1(iRecur), rP2(iRecur), rP3(iRecur)
             rp4(iRecur) = 0.0
             rp5(iRecur) = 0.0
             rp6(iRecur) = 0.0
           endif
         enddo

c        Read seismogenic thickness
         if ( sourceType(iFlt) .lt. 5) then
           read (10,*,err=3042) nThick1
           call S21_CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
           read (10,*,err=3043) (faultThick1(i),i=1,nThick1)
           read (10,*,err=3044) (faultThickWt1(i),i=1,nThick1)         
           call S21_CheckWt (faultThickWt1, nThick1, fName(iFlt), 'Seismo Thick        ')
         else
           nThick1 = 1
           faultTHick1(1) = -99.
           faultThickWt1(1) = 1.
         endif
         
c        Read depth pdf
         read (10,*,err=3045) iDepthModel(iFlt), (depthParam(iflt,k),k=1,3)     

c        Read Mag method (scaling relations or set values)
         read (10,*,err=3046) iOverRideMag
         if ( iOverRideMag .ne. 1 ) then
           write (*,'( 2x,''iOverRideMag flag option not working'')') 
           stop 99
         endif

c        Read reference mags for each fault thickness
         iThickDip = 1
         do iThick1=1,nThick1
           read (10,*,err=3047) nRefMag0
           read (10,*,err=3048) (refMag0(i),i=1,nRefMag0)
           read (10,*,err=3049) (refMagWt0(i),i=1,nRefMag0)

           if ( nRefMag0 .ne. 0 ) then
             call S21_CheckWt ( refMagWt0, nrefMag0, fName(iFlt), 'Ref Mag                 ' )
           endif
              
c          Copy these ref magnitudes for each addtional dip (no correlation allowed)
           do iDip=1,n_Dip
             nRefMag(iThickDip) = nRefMag0
             do i=1,nRefMag0
               refMag1(iThickDip,i) = refMag0(i)
               refMagWt1(iThickDip,i) = refMagWt0(i)
             enddo
             iThickDip = iThickDip + 1
           enddo
         enddo

c        Read min mag, step sizes, and rupture variability info
         read (10,*,err=3050) minMag(iFlt), magStep(iFlt), hxStep(iFlt), 
     1              hyStep(iFlt), nRupArea(iFlt), nRupWidth(iFlt), minDepth(iFlt)

C       Allow hxstep and hystep to be different for sourceType 2
        if (sourceType(iFlt).ne.2) then
C         Check that Hxstep = Hystep
          if (hxstep(iFlt) .ne. hystep(iFlt) ) then
            write (*,*) 'Hxstep not equal to Hystep for fault: '
            write (*,*) fname1
            write (*,*) 'Hxstep = ', hxstep(iFlt) 
            write (*,*) 'Hystep = ', hystep(iFlt) 
            write (*,*) 'These values must be equal or errors can occur!'
            write (*,*) 'Check input fault file.'
            stop 99
          endif
        endif

         read (10,*,err=3051) (coef_area(k,iFlt),k=1,2), sigArea(iFlt)
         read (10,*,err=3052) (coef_width(k,iFlt),k=1,2), sigWidth(iFlt)

c        Read ftype Models
         read (10,*,err=3053) nFtypeModels
         do iFM=1,nFtypeModels
           read (10,*,err=3054) ftmodelwt(iFM)
           read (10,*,err=3055) nFtype1(iFM)
           read (10,*,err=3056) (ftype1(iFM,k),k=1,nFtype1(iFM))
           read (10,*,err=3057) (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
           call S21_CheckWt1a (Ftype_Wt1, nFtype1(iFM), iFM, 
     1          MAX_N1, fName(iFlt), 'Fault Mech          ' )
         enddo

c        Load up Ftype Models and weights.
         nFM = 1
         do iFM=1,nFtypeModels
           do k=1,nFtype1(iFM)
             ftype(iFlt,nFM) = ftype1(iFM,k)
             ftype_wt(iFlt,nFM) = ftype_wt1(iFM,k)*ftmodelwt(iFM)
             nFm = nFm + 1
           enddo
         enddo
         nFm = nFm - 1
         nFtype(iFlt) = nFm

c        Load up parameter variations into large single dimension arrays
         testMaxMag = 0.
         iWidth = 0
         do iThick1=1,nThick1
        
           do iDip=1,n_Dip
           iWidth = iWidth + 1
           call S21_CheckDim ( iWidth, MAX_WIDTH, 'MAX_WIDTH' )   
           
             dip2 = dip1 + deltaDip1(iDip)
             dip(iFlt,iWidth,1) = dip2
             faultWidth(iFlt,iWidth) = faultThick1(iThick1)
             faultWidthWt(iFlt,iWidth) = faultThickWt1(iThick1) * dipWt1(iDip)

           i = 0
           mtest = 0.0

           do iRecur=1,nMagRecur    
           
            do iRate=1,nRate
              if ( rateType1(iRate) .eq. 2 ) then
                nb1 = 1
              else
                nb1 = n_bvalue
              endif
              do i_bValue=1,nb1
               do iRefMag=1,nRefMag(iWidth)
                 i = i + 1
                 call S21_CheckDim ( i, MAXPARAM, 'MAXPARAM  ' )
                 magRecur(iFlt,i,iWidth) = magRecur1(iRecur)
                 magRecurWt(iFlt,i,iWidth) = magRecurWt1(iRecur)

c                Scale moment rate by reference thickness.
                 if (rateType1(iRate) .eq. 4) then
                    rateParam(iFlt,i,iWidth) = rateParam1(iRate)*MoRDepth(iRate)*faultWidth(iFlt,iWidth)
                 else
                    rateParam(iFlt,i,iWidth) = rateParam1(iRate)
                 endif

                 rateType(iFlt,i,iWidth) = rateType1(iRate)
                 if ( rateType1(iRate) .eq. 2 ) then
                   beta(iFlt,i,iWidth) = bValue2(iRate - nSR)*alog(10.)
                 else
                   beta(iFlt,i,iWidth) = bValue1(i_bValue)*alog(10.)
                 endif
                 if ( rateType1(iRate) .eq. 2 ) then
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate)
                 else
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate) * bValueWt1(i_bValue) 
                 endif
                 maxMagWt(iFlt,i,iWidth) = refMagWt1(iWidth,iRefMag)

c                Set max mag
                 if (magRecur1(iRecur) .eq. 0 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 1 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP1(iRecur)
                 elseif (magRecur1(iRecur) .eq. 3 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 6 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 elseif (magRecur1(iRecur) .eq. 10 ) then
                   if ( rp5(iRecur) .eq. 0. ) then
                     maxMag(iFlt,i,iWidth) = rp2(iRecur)
                   else
                     maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rp5(iRecur)
                   endif
                 endif
                 if ( maxMag(iFlt,i,iWidth) .gt. testMaxMag ) then
                   nMag(iFlt) = ceiling((maxMag(iFlt,i,iWidth) - minMag(iFLt) ) / magStep(iFlt))
                   if (sourceType(iFlt) .eq. 7) then
                     nMag(iFlt) = ncountS7(iFlt)
                   endif
                     if (nMag(iFlt) .eq. 0 ) then
                       nMag(iFlt) = 1
                     endif    
                   testMaxMag = maxMag(iFlt,i,iWidth)
                 endif               

c  temp fix - the maxmag, refmag problem for waacy
                 if (magRecur1(iRecur) .eq. 10 ) then
                     maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 endif
                          
                 mpdf_param(iFlt,i,iWidth,1) = rP1(iRecur)
                 mpdf_param(iFlt,i,iWidth,2) = rP2(iRecur)
                 mpdf_param(iFlt,i,iWidth,3) = rP3(iRecur)
                 mpdf_param(iFlt,i,iWidth,4) = rP4(iRecur)
                 mpdf_param(iFlt,i,iWidth,5) = rP5(iRecur)
                 mpdf_param(iFlt,i,iWidth,6) = rP6(iRecur)

c      Test for the maximum magnitude value for all realizations of this fault.
                 if (maxMag(iFlt,i,iWidth) .ge. mtest) then
                   mtest = maxMag(iFlt,i,iWidth)
                 endif

                 mmagout(iFlt,iWidth,iRefMag) = refMag1(iWidth,iRefMag)
                 mmagoutwt(iFlt,iWidth,iRefMag) = refMagWt1(iWidth,iRefMag)
     
c     End loop over iRefMag
                enddo
c     End loop over ib_value
               enddo
c     End loop over iRate
              enddo
c     End loop over iRecur
             enddo
             nParamVar(iFlt,iWidth) = i
c     End loop over iDip
            enddo
c     End loop over iThick1
           enddo
           probAct(iFlt) = probAct0
           nWidth(iFlt) = iWidth
           fsys(iFlt) = ifsystem
c     End loop over iFlt2 (# segments)
          enddo
          ifsystem = ifsystem + 1          
c     End loop over iFlt
      enddo
      nfsystem = ifsystem - 1
      nFlt = iflt
      
      close (10)
      
      return
 3001 write (*,'( 2x,''Flt file error:  iCorr'')')
      stop 99
 3002 write (*,'( 2x,''Flt file error:  nFlt'')')
      stop 99
 3003 write (*,'( 2x,''Flt file error:  flt sys name'')')
      stop 99
 3004 write (*,'( 2x,''Flt file error:  prob Act'')')
      stop 99
 3005 write (*,'( 2x,''Flt file error:  nSegModel'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3006 write (*,'( 2x,''Flt file error:  Seg wts'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3007 write (*,'( 2x,''Flt file error:  number of segments'')')
      write (*,'( 2x,''fault sys: '',a80)') fName1
      stop 99
 3008  write (*,'( 2x,''Flt file error: seg flags for seg model'', i5)') i
       write (*,'( 2x,''From fault: '',a80)') fName1
       stop 99
 3009 write (*,'( 2x,''Flt file error: seg name for seg model '', i5)') i
      write (*,'( 2x,''From fault: '',a80)') fName1
      stop 99
 3010 write (*,'( 2x,''Flt file error: source type line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3011 write (*,'( 2x,''Flt file error: nSynRup line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3012 write (*,'( 2x,''Flt file error: syn Rup param line '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3013 write (*,'( 2x,''Flt file error: Aleatory seg wt '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3014 write (*,'( 2x,''Flt file error: dip and top '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3015 write (*,'( 2x,''Flt file error: number flt point '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3016 write (*,'( 2x,''Flt file error: long lat values '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3017 write (*,'( 2x,''Flt file error: nDowndip, nFP '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3018 write (*,'( 2x,''Flt file error: long lat values '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3019 write (*,'( 2x,''Flt file error: nDip '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3020 write (*,'( 2x,''Flt file error: delta Dip '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3021 write (*,'( 2x,''Flt file error: Dip wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3022 write (*,'( 2x,''Flt file error: n b-value '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3023 write (*,'( 2x,''Flt file error: delta b-value '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3024 write (*,'( 2x,''Flt file error: b-value wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3025 write (*,'( 2x,''Flt file error: nActRate'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3026 write (*,'( 2x,''Flt file error: b, activity, wt'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3027 write (*,'( 2x,''Flt file error: wts for activity rate approach'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3028 write (*,'( 2x,''Flt file error: number Slip rates '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3029 write (*,'( 2x,''Flt file error: slip rates '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3030 write (*,'( 2x,''Flt file error: slip-rate wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3031 write (*,'( 2x,''Flt file error: number Rec Intervals '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3032 write (*,'( 2x,''Flt file error: Rec Intervals '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3033 write (*,'( 2x,''Flt file error: Rec Interval wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3034 write (*,'( 2x,''Flt file error: Number Moment rate'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3035 write (*,'( 2x,''Flt file error: Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3036 write (*,'( 2x,''Flt file error: depths for Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3037 write (*,'( 2x,''Flt file error: wts for Moment rates'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3038 write (*,'( 2x,''Flt file error: number mag pdf'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3039 write (*,'( 2x,''Flt file error: mag pdf flag'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3040 write (*,'( 2x,''Flt file error: mag pdf wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3041 write (*,'( 2x,''Flt file error: param for mag pdf'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3042 write (*,'( 2x,''Flt file error: number crustal thickness '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3043 write (*,'( 2x,''Flt file error: crustal thicknesses '', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3044 write (*,'( 2x,''Flt file error: crustal thickness wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3045 write (*,'( 2x,''Flt file error: depth pdf and param'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99 
 3046 write (*,'( 2x,''Flt file error: iOverRideMag'', 2i5)') iFlt0, iFlt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99       
 3047 write (*,'( 2x,''Flt file error: depth pdf and param'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3048 write (*,'( 2x,''Flt file error: number Ref Mag'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3049 write (*,'( 2x,''Flt file error: Ref mags'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3050 write (*,'( 2x,''Flt file error: Ref mag wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3051 write (*,'( 2x,''Flt file error: coeff A(M) model'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3052 write (*,'( 2x,''Flt file error: coeff W(M) model'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3053 write (*,'( 2x,''Flt file error: number Ftype models'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3054 write (*,'( 2x,''Flt file error: Ftype mode wt'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3055 write (*,'( 2x,''Flt file error: number Ftype'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3056 write (*,'( 2x,''Flt file error: Ftypes'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99
 3057 write (*,'( 2x,''Flt file error: Ftype Aleatory wts'', 2i5)') iFlt0, iFLt2
      write (*,'( 2x,''From fault segment: '',a80)') fName(iFlt)
      stop 99

      end
