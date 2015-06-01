c  --------------------------------------------------------------------

      subroutine Rd_Fault_Data ( nFlt, fName, minMag, magStep, hxStep,
     1     hyStep, segModelWt, rateParam, rateParamWt, beta, 
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, 
     3     maxMag,  maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, 
     4     ftype, sourceType, nRupArea, coef_area, sigArea, nRupWidth, 
     5     coef_width, sigWidth, nParamVar, iCoor, minDepth, 
     6     fIndex, probAct, nWidth, mpdf_param, 
     7     al_segWt, attenType, sampleStep,
     8     grid_a, grid_dlong, grid_dlat, grid_n, grid_long, grid_lat,
     9     grid_top, minlat, maxlat, minlong, maxlong, scaleRate,
     1     fsys, 
     2     mMagout, mMagoutWt, fltDirect, synchron, nsyn_Case, synatten,
     3     synmag, syndistRup, syndistJB, synDistSeismo, synHypo,
     4     synftype, synhwflag, synwt, RateType, iDepthModel, depthParam,
     5     nMaxMag2, segwt1, faultFlag, nDownDip, nFtype, ftype_wt, 
     6     br_index, br_wt, segModelFlag, nSegModel0, segModelWt1, runflag )

      include 'pfrisk.h'
      
      real synmag(MAX_FLT,MAX_SYN), syndistRup(MAX_FLT,MAX_SYN)
      real syndistJB(MAX_FLT,MAX_SYN)
      real syndistSeismo(MAX_FLT,MAX_SYN), synftype(MAX_FLT,MAX_SYN)
      real synhypo(MAX_FLT,MAX_SYN), synwt(MAX_FLT,MAX_SYN)
      integer synHWFlag(MAX_FLT,MAX_SYN)
      integer nsyn_Case(MAX_FLT), synatten(MAX_FLT)
      integer fltDirect(MAX_FLT), synchron(MAX_FLT)
      real rateParam1(MAX_N1), RateWt1(MAX_N1)
      real al_segWt(1), dipWt1(MAX_N1), deltaDip1(MAX_N1)
      real bValue1(MAX_N1), bValueWt1(MAX_N1),  bValue2(MAX_N1)
      real faultThick1(MAX_N1), faultThickWt1(MAX_N1)
      real refMag1(MAX_N1,MAX_N1),  refMagWt1(MAX_N1,MAX_N1)
      real refMag0(MAX_N1), refMagWt0(MAX_N1)
      real magRecurWt1(MAX_N1), magRecur1(MAX_N1), probAct(MAX_FLT)
      integer fIndex(3,1), nWidth(1), sourceType(1), attenType(1)
      integer grid_n(MAX_FLT) 
      integer nMaxMag2(MAX_FLT), nMagRecur2(MAX_FLT), n_bValue2(MAX_FLT), nRate2(MAX_FLT)
      real minlat, maxlat, minlong, maxlong
      real grid_a(MAX_FLT,MAX_GRID)
      real grid_lat(MAX_FLT,MAX_GRID), grid_long(MAX_FLT,MAX_GRID), grid_top(MAX_FLT,MAX_GRID)
      real grid_dlong(MAX_FLT), grid_dlat(MAX_FLT)
      real minMag(1), magStep(1), hxStep(1), hyStep(1), minDepth(1)
      real segModelWt(1), sampleStep(1)
      real rateParam(MAX_FLT,MAXPARAM,MAX_WIDTH), rateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real beta(MAX_FLT,MAXPARAM,MAX_WIDTH),  betaWt(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH), magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real faultWidth(MAX_FLT,MAX_WIDTH), faultWidthWt(MAX_FLT,MAX_WIDTH)
      real mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,5)
      real rP1(MAXPARAM), rP2(MAXPARAM), rp3(MAXPARAM), rp4(MAXPARAM), rp5(MAXPARAM)
      real maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real coef_area(2,1), sigArea(1), coef_width(2,1), sigWidth(1)
      real fLong(MAX_FLT,MAX_DD,MAX_SEG), fLat(MAX_FLT,MAX_DD,MAX_SEG),
     1     fZ(MAX_FLT,MAX_DD,MAX_SEG)
      real ftype(MAX_FLT,MAX_N1), ftype_wt(MAX_FLT,MAX_N1)
      real ftype1(MAX_FLT,MAX_N1), ftype_wt1(MAX_FLT,MAX_N1), ftmodelwt(MAX_N1)
      integer nfp(MAX_FLT), nDownDip(MAX_FLT) 
      integer nMag(MAX_FLT), nRupArea(MAX_FLT), nRupWidth(MAX_FLT) 
      integer nRate, n_bValue, nRefMag(MAX_N1), nParamVar(MAX_FLT,1)
      character*80 fName(MAX_FLT), fName1
      integer ifsystem, nfsystem, fsys(MAX_FLT), nRuplength, nMaglength
      real rupLength(MAX_N1), rupLength_wt(MAX_N1), magApproach_wt(3)
      real c_magLength(4,MAX_N1), magLength_wt(MAX_N1)
      real c_magArea(4,MAX_N1), magArea_wt(MAX_N1)
      real mtest, mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real wt_srBranch, wt_ActrateBranch, wt_recIntBranch
      real sr(MAXPARAM), wt_sr(MAXPARAM)
      real actRate(MAXPARAM), actRateWt1(MAXPARAM)
      real MoRate(MAXPARAM), wt_MoRate(MAXPARAM)
      real MoRateDepth(MAXPARAM), MoRDepth(MAXPARAM)
      real rec_Int(MAXPARAM), wt_recInt(MAXPARAM)
      real rateType1(MAXPARAM), RateType(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real c_magDisp(4,MAX_N1), magDisp_wt(MAX_N1), aspectMin(MAX_N1)
      real disp(MAX_N1), disp_wt(MAX_N1)
      real depthParam(MAX_FLT,5)
      integer iDepthModel(MAX_FLT)
      real scaleRate(MAX_FLT)
      integer nFtype(1), faultFlag(MAX_FLT,100,MAX_FLT), nFtype1(MAX_FLT)
      real segWt(MAX_FLT,MAX_FLT), segWt1(MAX_FLT), lat1
      real dip(MAX_FLT,MAX_WIDTH, MAX_SEG)
      integer temp_BR(MAXPARAM), BR_index(MAX_FLT,20,MAX_WIDTH,MAXPARAM)
      real temp_BR_wt(MAXPARAM), BR_wt(MAX_FLT,20,MAX_WIDTH,MAXPARAM)
      integer segModelFlag(MAX_FLT,100), nSegModel0(1), runflag
      real segModelWt1(MAX_FLT,100)
      

c     Input Fault Parameters
      read (10,*) iCoor
      read (10,*) NFLT

      if (runflag .eq. 3) then
         write (*,*) 'iCoor = ', iCoor
         write (17,*) 'iCoor = ', iCoor
         write (*,*) 'NFLT = ' , nFLT
         write (17,*) 'NFLT = ' , nFLT
         write (*,*) 
         write (17,*) 
      endif  
     
      iflt = 0
      ifsystem = 1
      
      DO iFlt0=1,NFLT
       read (10,'( a80)') fName1
       read (10,*) probAct0

      if (runflag .eq. 3) then
         write (*,*) 'fName = ', fname1
         write (17,*)
         write (17,*) '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
         write (17,'(a13,i6,a4,i6)') 'Fault Number ', iFlt0, ' of ', NFLT
         write (17,*) 'fName1 = ', fName1
         write (*,*) 'ProbAct0 = ', ProbAct0
         write (17,*) 'ProbAct0 = ', ProbAct0
      endif  

c      Read number of segmentation models for this fault system       
       read (10,*) nSegModel
       read (10,*) (segWt(iFlt0,k),k=1,nSegModel)

      if (runflag .eq. 3) then
         write (*,*) 'nSegModel = ', nSegModel
         write (17,*) 'nSegModel = ', nSegModel
         do k=1,nSegModel
            write (*,*) 'SegWt = ', k,segWt(iFlt0,k)
            write (17,*) 'SetWt = ', k,segWt(iFlt0,k)
         enddo
      endif  

c      Read total number of fault segments 
       read (10,*) nFlt2
      if (runflag .eq. 3) then
         write (*,*) 'nFlt2 = ', nFlt2
         write (17,*) 'nFlt2 = ', nFlt2
      endif  

       do i=1,nSegModel
         read (10,*) (faultFlag(iFlt0,i,k),k=1,nFlt2)

      if (runflag .eq. 3) then
         write (*,*) 'faultFlag = ', i,(faultFlag(iFlt0,i,k),k=1,nFlt2)
         write (17,*) 'faultFlag = ', i,(faultFlag(iFlt0,i,k),k=1,nFlt2)
      endif  

       enddo
                        
       do iflt2=1,nflt2
        iFlt = iFlt + 1
        if (runflag .ne. 3) then
           call CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )
        endif

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
        read(10,'( a80)') fName(iFlt)
        write (*,'( i5,2x,a80)') iFlt, fName(iFlt)
        if (runflag .eq. 3) then
           write (17,'( a13,i5,2x,a80)') 'iFlt, Fname: ', iFlt, fName(iFlt)
        endif

c       Read type of source (planar, areal, grid1, grid2, irregular)
        read(10,*) sourceType(iFlt),attenType(iFlt),sampleStep(iFlt),
     1                 fltDirect(iFlt), synchron(iFlT)

        if (runflag .eq. 3) then
           write (*,*) iFlt,sourceType(iFlt),attenType(iFlt),sampleStep(iFlt),
     1                 fltDirect(iFlt), synchron(iFlT)
           write (17,*) iFlt,sourceType(iFlt),attenType(iFlt),sampleStep(iFlt),
     1                 fltDirect(iFlt), synchron(iFlT)
        endif
          
c       Now read in the synchronous rupture case parameters if needed.
        if (synchron(iFlt) .gt. 0) then
           read (10,*) nsyn_Case(iFlt), synatten(iFlt)

        if (runflag .eq. 3) then
           write (*,*) 'NSyn_Case, snyatten = ', nsyn_Case(iFlt), synatten(iFlt)
           write (17,*) 'NSyn_Case, snyatten = ', nsyn_Case(iFlt), synatten(iFlt)
        endif

c          Now read in the magnitude, distance and weigths for each synchronous rupture case.
           do isyn=1,nsyn_Case(iFlt)
             read (10,*) synmag(iFlt,isyn),syndistRup(iFlt,isyn),
     1                 syndistJB(iFlt,isyn),syndistSeismo(iFlt,isyn),
     2                 synHWFlag(iFlt,isyn),synhypo(iFlt,isyn),
     3                 synftype(iFlt,isyn),synwt(iFlt,isyn)

        if (runflag .eq. 3) then
           write (*,*) 'Read Synmag information: ', synmag(iFlt,isyn),syndistRup(iFlt,isyn),
     1                 syndistJB(iFlt,isyn),syndistSeismo(iFlt,isyn),
     2                 synHWFlag(iFlt,isyn),synhypo(iFlt,isyn),
     3                 synftype(iFlt,isyn),synwt(iFlt,isyn)
           write (17,*) 'Read Synmag information: ', synmag(iFlt,isyn),syndistRup(iFlt,isyn),
     1                 syndistJB(iFlt,isyn),syndistSeismo(iFlt,isyn),
     2                 synHWFlag(iFlt,isyn),synhypo(iFlt,isyn),
     3                 synftype(iFlt,isyn),synwt(iFlt,isyn)
        endif

           enddo
        endif

c       Read aleatory segmentation wts
        read (10,*) al_segWt(iFlt)

        if (runflag .eq. 3) then
           write (*,*) 'Al_Segwt = ', al_segWt(iFlt)
           write (17,*) 'Al_Segwt = ', al_segWt(iFlt)
        endif

c       Check for standard fault source or areal source
        if ( sourceType(iFlt) .eq. 1 .or. sourceType(iFLt) .eq. 2) then
          read (10,*) dip1, top

        if (runflag .eq. 3) then
           write (*,*) 'Dip1, Top = ', dip1, top
           write (17,*) 'Dip1, Top = ', dip1, top
        endif

          read(10,*) nfp(iFlt)     

        if (runflag .eq. 3) then
           write (*,*) 'nFp = ', nfp(iFlt)
           write (17,*) 'nFp = ', nfp(iFlt)
        endif

          call CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
          do ipt=1,nfp(iFlt)
            read (10,*) fLong(iFlt,1,ipt), fLat(iFlt,1,ipt)
            fZ(iFlt,1,ipt) = top
            if (sourceType(iFlt) .eq. 2) then
               grid_top(iFlt,ipt) = top
            endif

        if (runflag .eq. 3) then
           write (*,*) 'ipt, Flong, fLat = ',ipt,fLong(iFlt,1,ipt), fLat(iFlt,1,ipt)
           write (17,*) 'ipt, Flong, fLat = ', ipt,fLong(iFlt,1,ipt), fLat(iFlt,1,ipt)
        endif

          enddo
          nDownDip(iFlt) = 1
         endif

c        Check for grid source (w/o depth)
         if ( sourceType(iFlt) .eq. 3 ) then
           read (10,*)  dip1, top

        if (runflag .eq. 3) then
           write (*,*) 'Dip1, top = ', dip1, top
           write (17,*) 'Dip1, top = ', dip1, top
        endif

           call RdGrid1 (iFlt, grid_a, grid_dlong, grid_dlat, grid_n,
     1             grid_long, grid_lat, minLat, minLong,  maxLat, maxLong, scaleRate(iFlt) )

        if (runflag .eq. 3) then
           write (*,*) 'Successful call to RdGrid1 Subroutine.'
           write (17,*) 'Successful call to RdGrid1 Subroutine.'
        endif

           do igrid=1,grid_n(iFlt)
             grid_top(iFlt,igrid) = top
           enddo
           nDownDip(iFlt) = 1
         endif

c        Check for grid source (w/ depth)
         if ( sourceType(iFlt) .eq. 4 ) then
              read (10,*)  dip1

        if (runflag .eq. 3) then
           write (*,*) 'Dip1 = ', dip1
           write (17,*) 'Dip1 = ', dip1
        endif

              call RdGrid2 (iFlt,grid_a,grid_dlong,grid_dlat,grid_n,
     1             grid_long, grid_lat, minLat, minLong, maxLat, maxLong, scaleRate(iFlt),
     2             grid_top)

        if (runflag .eq. 3) then
           write (*,*) 'Successful call to RdGrid2 Subroutine.'
           write (17,*) 'Successful call to RdGrid2 Subroutine.'
        endif

              nDownDip(iFlt) = 1
         endif

c        Check for custom fault source
         if ( sourceType(iFlt) .eq. 5 .or. sourceType(iFlt) .eq. 6) then
           read(10,*) nDownDip(iFLt), nfp(iFlt)  
   
        if (runflag .eq. 3) then
           write (*,*) 'nDownDip, nFp = ', nDownDip(iFLt), nfp(iFlt)
           write (17,*) 'nDownDip, nFp = ', nDownDip(iFLt), nfp(iFlt)
        endif

           call CheckDim ( nfp(iFlt), MAX_SEG, 'MAX_SEG   ' )
           do ipt=1,nfp(iFlt)
              read (10,*) (fLong(iFlt,k,ipt), fLat(iFlt,k,ipt), fZ(iflt,k,ipt), k=1,nDownDip(iFlt) ) 

         if (runflag .eq. 3) then
           write (*,*) 'Flong, Flat, Fz: ', (fLong(iFlt,k,ipt), fLat(iFlt,k,ipt), fZ(iflt,k,ipt), k=1,nDownDip(iFlt) ) 
           write (17,*) 'Flong, Flat, Fz: ', (fLong(iFlt,k,ipt), fLat(iFlt,k,ipt), fZ(iflt,k,ipt), k=1,nDownDip(iFlt) ) 
        endif

          enddo
         endif
               
c        Read dip Variation
         if ( sourceType(iFlt) .lt. 5 ) then
           read (10,*) n_Dip

         if (runflag .eq. 3) then
           write (*,*) 'n_Dip = ',n_Dip
           write (17,*) 'n_Dip = ', n_Dip 
        endif

           call CheckDim ( n_Dip, MAX_N1, 'MAX_N1    ' )
           read (10,*) (deltaDip1(i),i=1,n_Dip)
           read (10,*) (dipWt1(i),i=1,n_Dip)

        if (runflag .eq. 3) then
           do i=1,n_dip
              write (*,*) 'DeltaDip, DipWt = ', i, deltaDip1(i), dipWt1(i)
              write (17,*) 'DeltaDip, DipWt = ', i, deltaDip1(i), dipWt1(i)
           enddo
        endif

           call CheckWt ( dipWt1, n_Dip, fName(iFlt), 'Dip                 ' )
         else
           n_Dip = 1
           deltaDip1(1) = 0.
           dipWt1(1) = 1.
         endif

c        Read b-values (not for activity rate cases)
         read (10,*) n_bValue

        if (runflag .eq. 3) then
           write (*,*) 'nb_value = ', n_bValue
           write (17,*) 'nb_value = ', n_bValue
        endif

         call CheckDim ( n_bValue, MAX_N1, 'MAX_N1    ' )
         if ( n_bValue .gt. 0 ) then
           read (10,*) (bValue1(i),i=1,n_bValue)
           read (10,*) (bValueWt1(i),i=1,n_bValue)
         endif

        if (runflag .eq. 3) then
           do i=1,n_bValue
              write (*,*) 'bValue, bValueWt = ', i, bValue1(i), bValueWt1(i)
              write (17,*) 'bValue, bValueWt = ', i, bValue1(i), bValueWt1(i)
           enddo
        endif
                
c        Read activity rate - b-value pairs
         read (10,*) nActRate 

        if (runflag .eq. 3) then
           write (*,*) 'nActRate = ', nActRate
           write (17,*) 'nActRate = ', nActRate
        endif

         if ( nActRate .ne. 0 ) then
           do ii=1,nActRate
             read (10,*) bValue2(ii), actRate(ii), actRateWt1(ii)
C     Scale activity rate for gridded source (i.e., sourcetype 3 or 4) based on limited lat and long values.
             if (sourcetype(iFlt) .eq. 3 .or. sourcetype(iFlt) .eq. 4) then
                if (scalerate(iFlt) .ne. 1.0) then
                   actRate(ii) = actRate(ii)*Scalerate(iFlt)
                endif
             endif
           enddo

        if (runflag .eq. 3) then
           do ii=1,nActRAte
              write (*,*) 'bValue, ActRate, ActRateWt = ', ii,bValue2(ii), actRate(ii), actRateWt1(ii)
              write (17,*) 'bValue, ActRate, ActRateWt = ', ii,bValue2(ii), actRate(ii), actRateWt1(ii)
           enddo
        endif

           call CheckWt ( actRateWt1, nActRate, fName(iFlt), 'ActRate             ' )
         endif

c        Read weights for rate methods
         read (10,*) wt_srBranch, wt_ActRateBranch, wt_recIntBranch, wt_MoRateBranch

        if (runflag .eq. 3) then
           write (*,*) 'wtSR, wtActRate, wtRecInt, wtMoRate = ',  wt_srBranch, wt_ActRateBranch, 
     1                        wt_recIntBranch, wt_MoRateBranch
           write (17,*) 'wtSR, wtActRate, wtRecInt, wtMoRate = ',  wt_srBranch, wt_ActRateBranch, 
     1                        wt_recIntBranch, wt_MoRateBranch
        endif

         sum = wt_srBranch + wt_ActRateBranch + wt_recIntBranch + wt_MoRateBranch
         if ( sum .lt. 0.999 .or. sum .gt. 1.001 ) then
              write (*,'( 2x,''rate method weights do not sum to unity for fault, '',a30)') fName(iFlt)
              stop 99
         endif

c        Read slip-rates
         read (10,*) nSR
        if (runflag .eq. 3) then
           write (*,*) 'nSR = ',  nSR
           write (17,*) 'nSR = ',  nSR
        endif

         if ( nSR .gt. 0 ) then
           read (10,*) (sr(k),k=1,nSR)
           read (10,*) (wt_sr(k),k=1,nSR)

        if (runflag .eq. 3) then
           do k=1,nSr
              write (*,*) 'SR, wt = ',  sr(k),wt_sr(k)
              write (17,*) 'SR, wt = ', sr(k),wt_sr(k)
           enddo
        endif

           call CheckWt (wt_sr, nSR, fName(iFlt), 'Slip Rates          ')
         endif

c        Read recurrence intervals
         read (10,*) nRecInt

        if (runflag .eq. 3) then
           write (*,*) 'nRecInt = ',  nRecInt
           write (17,*) 'nRecInt = ',  nRecInt
        endif

         if ( nRecInt .gt. 0 ) then
           read (10,*) (rec_Int(k),k=1,nRecInt)
           read (10,*) (wt_recInt(k),k=1,nRecInt)

        if (runflag .eq. 3) then
           do k=1,nRecInt
              write (*,*) 'RecInt, wt = ',  rec_Int(k),wt_recInt(k)
              write (17,*) 'RecInt, wt = ', rec_Int(k),wt_recInt(k)
           enddo
        endif

           call CheckWt (wt_recInt, nRecInt, fName(iFlt), 'Recurrence Intervals')
         endif

c        Read moment-rates
         read (10,*) nMoRate

        if (runflag .eq. 3) then
           write (*,*) 'nMoRate = ',  nMoRate
           write (17,*) 'nMoRate = ',  nMoRate
        endif

          if ( nMoRate .gt. 0 ) then
           read (10,*) (MoRate(k),k=1,nMoRate)
           read (10,*) (MoRateDepth(k),k=1,nMoRate)
           read (10,*) (wt_MoRate(k),k=1,nMoRate)

        if (runflag .eq. 3) then
           do k=1,nMoRate
              write (*,*) 'NoRate, MoRateDepth, wt = ',  MoRate(k),MoRateDepth(k),wt_MoRate(k)
              write (17,*) 'MoRate, MoRateDepth, wt = ', MoRate(k),MoRateDepth(k),wt_MoRate(k)
           enddo
        endif

           call CheckWt (wt_MoRate, nMoRate, fName(iFlt), 'Moment Rates        ')
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
         call CheckDim ( nRate, MAX_N1, 'MAX_N1    ' )
         do k=1,nSR
            rateParam1(k) = sr(k)
            rateWt1(k) = wt_sr(k)*wt_srbranch              
            rateType1(k) = 1
            temp_BR(k) = k
            temp_BR_wt(k) = wt_sr(k)
            MoRDepth(k) = 1.0
         enddo
         do k=1,nActRate
            rateParam1(k+nSR) = actRate(k)
            rateWt1(k+nSR) = actRateWt1(k) * wt_actRateBranch             
            rateType1(k+nSR) = 2
            temp_BR(k+nSR) = k
            temp_BR_wt(k+nSR) = actRateWt1(k)
            MoRDepth(k+nSR) = 1.0
         enddo
         do k=1,nRecInt
            rateParam1(k+nSR+nActRate) = rec_Int(k)
            rateWt1(k+nSR+nActRate) = wt_recInt(k) *wt_recIntBranch              
            rateType1(k+nSR+nActRate) = 3
            temp_BR(k+nSR+nActRate) = k
            temp_BR_wt(k+nSR+nActRate) = wt_recInt(k)
            MoRDepth(k+nSR+nActRate) = 1.0
         enddo
         do k=1,nMoRate
            rateParam1(k+nSR+nActRate+nRecInt) = MoRate(k)
            rateWt1(k+nSR+nActRate+nRecInt) = wt_MoRate(k) *wt_MoRateBranch              
            rateType1(k+nSR+nActRate+nRecInt) = 4
            temp_BR(k+nSR+nActRate+nRecInt) = k
            temp_BR_wt(k+nSR+nActRate+nRecInt) =  wt_MoRate(k)
            MoRDepth(k+nSR+nActRate+nRecInt) = 1.0/MoRateDepth(nMoRate)
         enddo
                     
c        Read Mag recurrence weights (char, exp, etc.)
         read (10,*) nMagRecur

        if (runflag .eq. 3) then
              write (*,*) 'nMagRecur = ',  nMagRecur
              write (17,*) 'nMagRecur = ', nMagRecur
        endif

         call CheckDim ( nMagRecur, MAX_N1, 'MAX_N1    ' )
         read (10,*) (magRecur1(i),i=1,nMagRecur)
         read (10,*) (magRecurWt1(i),i=1,nMagRecur)

        if (runflag .eq. 3) then
           do k=1,nMagRecur
              write (*,*) 'MagRecur1, wt = ',  magRecur1(k),magRecurWt1(k)
              write (17,*) 'MagRecur1, wt = ', magRecur1(k),magRecurWt1(k)
           enddo
        endif

         call CheckWt ( magRecurWt1,nMagRecur, fName(iFlt), 'Mag Recur           ' )

c        Read in corresponding magnitude parameters. 
         do iRecur=1,nMagRecur
           if (magRecur1(iRecur) .eq. 4) then
             read (10,*) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
c          Read in necessary parameters for bi-exponential distribution.
           elseif (magRecur1(iRecur) .eq. 5) then
             read (10,*) rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
           else              
             read (10,*) rP1(iRecur), rP2(iRecur), rP3(iRecur)
             rp4(iRecur) = 0.0
             rp5(iRecur) = 0.0
           endif

        if (runflag .eq. 3) then
           write (*,*) 'rP Array = ', rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
           write (17,*) 'rP Array = ', rP1(iRecur), rP2(iRecur), rP3(iRecur), rp4(iRecur), rp5(iRecur)
        endif

         enddo

c        Read seismogenic thickness
         if ( sourceType(iFlt) .lt. 5) then
           read (10,*) nThick1

        if (runflag .eq. 3) then
           write (*,*) 'nThick1 = ', nThick1
           write (17,*) 'nThick1 = ', nThick1
        endif

           call CheckDim ( nThick1, MAX_WIDTH, 'MAX_WIDTH ' )
           read (10,*) (faultThick1(i),i=1,nThick1)
           read (10,*) (faultThickWt1(i),i=1,nThick1)

        if (runflag .eq. 3) then
           do i=1,nThick1
              write (*,*) 'FaultThick1, wt = ', faultThick1(i), faultThickWt1(i)
              write (17,*) 'FaultThick1, wt = ', faultThick1(i), faultThickWt1(i)
           enddo
        endif
           
           call CheckWt (faultThickWt1, nThick1, fName(iFlt), 'Seismo Thick        ')
         else
           nThick1 = 1
           faultTHick1(1) = -99.
           faultThickWt1(1) = 1.
         endif
         
c        Read depth pdf
         read (10,*) iDepthModel(iFlt), (depthParam(iflt,k),k=1,3)     
  
        if (runflag .eq. 3) then
              write (*,*) 'Depth Model = ', iDepthModel(iFlt), (depthParam(iflt,k),k=1,3)    
              write (17,*) 'Depth Model = ', iDepthModel(iFlt), (depthParam(iflt,k),k=1,3)    
        endif

c        Read Mag method (scaling relations or set values)
         read (10,*) iOverRideMag

        if (runflag .eq. 3) then
              write (*,*) 'OverRideMag = ', iOverRideMag   
              write (17,*) 'OverRideMag = ', iOverRideMag   
        endif

         if ( iOverRideMag .eq. 1 ) then
c         Read reference mags for each fault thickness
          iThickDip = 1
          do iThick1=1,nThick1
            read (10,*) nRefMag0
            read (10,*) (refMag0(i),i=1,nRefMag0)
            read (10,*) (refMagWt0(i),i=1,nRefMag0)

        if (runflag .eq. 3) then
              do i=1,nRefMag0
                 write (*,*) 'Max Mags = ', refMag0(i), refMagWt0(i)  
                 write (17,*) 'Max Mags = ', refMag0(i), refMagWt0(i)    
              enddo
        endif

            if ( nRefMag0 .ne. 0 ) then
              call CheckWt ( refMagWt0, nrefMag0, fName(iFlt), 'Ref Mag                 ' )
            endif
              
c           Copy these ref magnitudes for each addtional dip (no correlation allowed)
            do iDip=1,n_Dip
              nRefMag(iThickDip) = nRefMag0
              do i=1,nRefMag0
                refMag1(iThickDip,i) = refMag0(i)
                refMagWt1(iThickDip,i) = refMagWt0(i)
              enddo
              iThickDip = iThickDip + 1
            enddo
          enddo
         else

c   *** New read for Mean Char Eqk computed in Code 	
c       (i.e. OverridgeMag .ne. 1) 
              read (10,*) nRupLength
              read (10,*) (rupLength(k),k=1,nRupLength)
              read (10,*) (rupLength_wt(k),k=1,nRupLength)
              read (10,*) (magApproach_wt(k),k=1,3)

        if (runflag .eq. 3) then
              write (*,*) 'Rupture Lenght = ', nRupLength, (rupLength(k),k=1,nRupLength),
     1             (rupLength_wt(k),k=1,nRupLength),(magApproach_wt(k),k=1,3)
              write (17,*) 'Rupture Lenght = ', nRupLength, (rupLength(k),k=1,nRupLength),
     1             (rupLength_wt(k),k=1,nRupLength),(magApproach_wt(k),k=1,3) 
        endif

c             read mag-length models
              read (10,*) nMaglength

        if (runflag .eq. 3) then
              write (*,*) 'nMagLength = ', nMaglength
              write (17,*) 'nMagLength = ', nMaglength
        endif

              if ( nMagLength .ne. 0 ) then
                do k=1,nMagLength
                  read (10,*) (c_magLength(kk,k), kk=1,4)

        if (runflag .eq. 3) then
              write (*,*) 'c_magLength = ', (c_magLength(kk,k), kk=1,4)
              write (17,*) 'c_magLength = ', (c_magLength(kk,k), kk=1,4)
        endif

                enddo
                read (10,*) (magLength_wt(k),k=1,nMagLength)

        if (runflag .eq. 3) then
              write (*,*) 'magLength_wt = ', (magLength_wt(k),k=1,nMagLength)
              write (17,*) 'magLength_wt = ', (magLength_wt(k),k=1,nMagLength)
        endif

              endif

c             read mag-area models
              read (10,*) nMagArea

        if (runflag .eq. 3) then
              write (*,*) 'nMagArea = ', nMagArea
              write (17,*) 'nMagArea = ', nMagArea
        endif

              if ( nMagArea .ne. 0 ) then
                do k=1,nMagArea
                  read (10,*) (c_magArea(kk,k), kk=1,4), aspectMin(k)

        if (runflag .eq. 3) then
              write (*,*) 'c_magArea, Aspect = ', (c_magArea(kk,k), kk=1,4), aspectMin(k)
              write (17,*) 'c_magArea, Aspect = ', (c_magArea(kk,k), kk=1,4), aspectMin(k)
        endif

                enddo
                read (10,*) (magArea_wt(k),k=1,nMagArea)

        if (runflag .eq. 3) then
              write (*,*) 'c_magAreaWt = ', (magArea_wt(k),k=1,nMagArea)
              write (17,*) 'c_magAreaWt = ', (magArea_wt(k),k=1,nMagArea)
        endif

              endif
               
c             read mag-displacement models
              read (10,*) nMagDisp
              
        if (runflag .eq. 3) then
              write (*,*) 'nMagDisp = ', nMagDisp
              write (17,*) 'nMagDisp = ', nMagDisp
        endif

              if ( nMagDisp .ne. 0 ) then
                do k=1,nMagDisp
                  read (10,*) (c_magDisp(kk,k), kk=1,4)

        if (runflag .eq. 3) then
              write (*,*) 'c_magDisp = ', (c_magDisp(kk,k), kk=1,4)
              write (17,*) 'c_magDisp = ',  (c_magDisp(kk,k), kk=1,4)
        endif

                enddo
                read (10,*) (magDisp_wt(k),k=1,nMagDisp)
                read (10,*) nDisp
                read (10,*) (disp(k), k=1,nDisp)
                read (10,*) (Disp_wt(k),k=1,nDisp)

        if (runflag .eq. 3) then
              write (*,*) 'magDispWt = ', (magDisp_wt(k),k=1,nMagDisp)
              write (*,*) 'nDips = ', nDisp
              write (*,*) 'Disp = ', (disp(k), k=1,nDisp)
              write (*,*) 'DispWt = ', (disp_wt(k), k=1,nDisp)
              write (17,*) 'magDispWt = ', (magDisp_wt(k),k=1,nMagDisp)
              write (17,*) 'nDips = ', nDisp
              write (17,*) 'Disp = ', (disp(k), k=1,nDisp)
              write (17,*) 'DispWt = ', (disp_wt(k), k=1,nDisp)
        endif

               endif
         endif

         read (10,*) minMag(iFlt), magStep(iFlt), hxStep(iFlt), 
     1              hyStep(iFlt), nRupArea(iFlt), nRupWidth(iFlt), minDepth(iFlt)

        if (runflag .eq. 3) then
              write (*,*) 'Fault parameters = ', minMag(iFlt), magStep(iFlt), hxStep(iFlt), 
     1              hyStep(iFlt), nRupArea(iFlt), nRupWidth(iFlt), minDepth(iFlt)
              write (17,*) 'Fault parameters = ', minMag(iFlt), magStep(iFlt), hxStep(iFlt), 
     1              hyStep(iFlt), nRupArea(iFlt), nRupWidth(iFlt), minDepth(iFlt)
        endif

C     Allow hxstep and hystep to be different for sourceType 2
C     cdh 4/29/2015
        if (sourceType(iFlt).ne.2) then
C     Check that Hxstep = Hystep
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

         read (10,*) (coef_area(k,iFlt),k=1,2), sigArea(iFlt)
         read (10,*) (coef_width(k,iFlt),k=1,2), sigWidth(iFlt)

        if (runflag .eq. 3) then
              write (*,*) 'Area and Width Model = ', (coef_area(k,iFlt),k=1,2), sigArea(iFlt)
              write (17,*) 'Area and Width Model = ', (coef_area(k,iFlt),k=1,2), sigArea(iFlt)
        endif

c        Read ftype Models
         read (10,*) nFtypeModels

        if (runflag .eq. 3) then
              write (*,*) 'nFtypeModels = ', nFtypeModels
              write (17,*) 'nFtypeModels = ', nFtypeModels
        endif

         do iFM=1,nFtypeModels
            read (10,*) ftmodelwt(iFM)
            read (10,*) nFtype1(iFM)

        if (runflag .eq. 3) then
              write (*,*) 'ftmodelwt, nFtype1 = ', ftmodelwt(iFM),nFtype1(iFM)
              write (17,*) 'ftmodelwt, nFtype1 = ', ftmodelwt(iFM),nFtype1(iFM)
        endif

            read (10,*) (ftype1(iFM,k),k=1,nFtype1(iFM))
            read (10,*) (ftype_wt1(iFM,k), k=1,nFtype1(iFM))

        if (runflag .eq. 3) then
              write (*,*) 'ftype1, ftype_wt = ', (ftype1(iFM,k),k=1,nFtype1(iFM)),
     1             (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
              write (17,*) 'ftype1, ftype_wt = ', (ftype1(iFM,k),k=1,nFtype1(iFM)),
     1             (ftype_wt1(iFM,k), k=1,nFtype1(iFM))
        endif

            call CheckWt1a (Ftype_Wt1, nFtype1(iFM), iFM, 
     1                  MAX_N1, fName(iFlt), 'Fault Mech          ' )
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

c        Build up distribution of mean Char Mag for each dip,thickness
         if ( iOverRideMag .ne. 1 ) then
           write (*,'( 2x,''this option not working'')')
           stop 99
c          call rdMagModel ()
         endif
 
c        Load up parameter variations into large single dimension arrays        
         testMaxMag = 0.
         iWidth = 0
         do iDip=1,n_Dip
          do iThick1=1,nThick1
           iWidth = iWidth + 1
           call CheckDim ( iWidth, MAX_WIDTH, 'MAX_WIDTH' )
 
           if ( sourceType(iFlt) .eq. 1 ) then
             
             dip2 = dip1 + deltaDip1(iDip)
             dip(iFlt,iWidth,1) = dip2
             faultWidth(iFlt,iWidth) = faultThick1(iThick1)
             faultWidthWt(iFlt,iWidth) = faultThickWt1(iThick1) * dipWt1(iDip)
           else

c            For areal source only set the first point and use thickness not width
             faultWidth(iFlt,iWidth) = faultThick1(iThick1)
             faultWidthWt(iFlt,iWidth) = faultThickWt1(iThick1) * dipWt1(iDip)
                
           endif
                
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
               do iRefMag=1,nRefMag(iThick1)
                 i = i + 1
                 if (runflag .ne. 3) then
                    call CheckDim ( i, MAXPARAM, 'MAXPARAM  ' )
                 endif
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
c                refMag(iFlt,i,iWidth) = refMag1(iThick1,iRefMag)
                 if ( rateType1(iRate) .eq. 2 ) then
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate)
                 else
                   RateParamWt(iFlt,i,iWidth) = RateWt1(iRate) * bValueWt1(i_bValue) 
                 endif
                 maxMagWt(iFlt,i,iWidth) = refMagWt1(iThick1,iRefMag)

c                Set max mag
                 if (magRecur1(iRecur) .eq. 0 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 1 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP1(iRecur)
                 elseif (magRecur1(iRecur) .eq. 3 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag) + rP3(iRecur)
                 elseif (magRecur1(iRecur) .eq. 6 ) then
                   maxMag(iFlt,i,iWidth) = refMag1(iWidth,iRefMag)
                 endif
                 if ( maxMag(iFlt,i,iWidth) .gt. testMaxMag ) then
                   nMag(iFlt) = nint((maxMag(iFlt,i,iWidth) - minMag(iFLt) ) / magStep(iFlt))
                     if (nMag(iFlt) .eq. 0 ) then
                       nMag(iFlt) = 1
                     endif    
                   testMaxMag = maxMag(iFlt,i,iWidth)
                 endif               
                           
                 mpdf_param(iFlt,i,iWidth,1) = rP1(iRecur)
                 mpdf_param(iFlt,i,iWidth,2) = rP2(iRecur)
                 mpdf_param(iFlt,i,iWidth,3) = rP3(iRecur)
                 mpdf_param(iFlt,i,iWidth,4) = rP4(iRecur)
                 mpdf_param(iFlt,i,iWidth,5) = rP5(iRecur)

c      Test for the maximum magnitude value for all realizations of this fault.
                 if (maxMag(iFlt,i,iWidth) .ge. mtest) then
                   mtest = maxMag(iFlt,i,iWidth)
                 endif

                 mmagout(iFlt,iWidth,iRefMag) = refMag1(iWidth,iRefMag)
                 mmagoutwt(iFlt,iWidth,iRefMag) = refMagWt1(iWidth,iRefMag)

c                Set branches for tornado plots
                 BR_index(iFlt,1,iWidth,i) = iDip
                 if ( rateType1(iRate) .eq. 1. ) then
                   BR_index(iFlt,2,iWidth,i) = 1
                   BR_index(iFlt,3,iWidth,i) = i_Bvalue
                   BR_index(iFlt,4,iWidth,i) = temp_BR(iRate)
                   BR_index(iFlt,5,iWidth,i) = 0
                   BR_index(iFlt,6,iWidth,i) = 0
                   BR_index(iFlt,7,iWidth,i) = 0
                 elseif (rateType1(iRate) .eq. 2. ) then
                   BR_index(iFlt,2,iWidth,i) = 2
                   BR_index(iFlt,3,iWidth,i) = 0
                   BR_index(iFlt,4,iWidth,i) = 0
                   BR_index(iFlt,5,iWidth,i) = temp_BR(iRate)
                   BR_index(iFlt,6,iWidth,i) = 0
                   BR_index(iFlt,7,iWidth,i) = 0
                 elseif (rateType1(iRate) .eq. 3. ) then
                   BR_index(iFlt,2,iWidth,i) = 3
                   BR_index(iFlt,3,iWidth,i) = i_Bvalue
                   BR_index(iFlt,4,iWidth,i) = 0
                   BR_index(iFlt,5,iWidth,i) = 0
                   BR_index(iFlt,6,iWidth,i) = temp_BR(iRate)
                   BR_index(iFlt,7,iWidth,i) = 0
                 elseif (rateType1(iRate) .eq. 4. ) then
                   BR_index(iFlt,2,iWidth,i) = 4
                   BR_index(iFlt,3,iWidth,i) = i_Bvalue
                   BR_index(iFlt,4,iWidth,i) = 0
                   BR_index(iFlt,5,iWidth,i) = 0
                   BR_index(iFlt,6,iWidth,i) = 0
                   BR_index(iFlt,7,iWidth,i) = temp_BR(iRate)
                 endif
                 BR_index(iFlt,8,iWidth,i) = iRecur
                 BR_index(iFlt,9,iWidth,i) = iThick1
                 BR_index(iFlt,10,iWidth,i) = 1
                 BR_index(iFlt,11,iWidth,i) = iRefMag

c                Set branches wts for tornado plots
                 BR_wt(iFlt,1,iWidth,i) = dipWt1(iDip)
                 if ( rateType1(iRate) .eq. 1. ) then
                   BR_wt(iFlt,2,iWidth,i) = wt_srbranch 
                   BR_wt(iFlt,3,iWidth,i) = bValueWt1(i_bValue)
                   BR_wt(iFlt,4,iWidth,i) = temp_BR_wt(iRate)
                   BR_wt(iFlt,5,iWidth,i) = 0
                   BR_wt(iFlt,6,iWidth,i) = 0
                   BR_wt(iFlt,7,iWidth,i) = 0
                 elseif (rateType1(iRate) .eq. 2. ) then
                   BR_wt(iFlt,2,iWidth,i) = wt_actRateBranch 
                   BR_wt(iFlt,3,iWidth,i) = 0
                   BR_wt(iFlt,4,iWidth,i) = 0
                   BR_wt(iFlt,5,iWidth,i) = temp_BR_wt(iRate)
                   BR_wt(iFlt,6,iWidth,i) = 0
                   BR_wt(iFlt,7,iWidth,i) = 0
                 elseif (rateType1(iRate) .eq. 3. ) then
                   BR_wt(iFlt,2,iWidth,i) = wt_recIntBranch 
                   BR_wt(iFlt,3,iWidth,i) = bValueWt1(i_bValue)
                   BR_wt(iFlt,4,iWidth,i) = 0
                   BR_wt(iFlt,5,iWidth,i) = 0
                   BR_wt(iFlt,6,iWidth,i) = temp_BR_wt(iRate)
                   BR_wt(iFlt,7,iWidth,i) = 0
                 elseif (rateType1(iRate) .eq. 4. ) then
                   BR_wt(iFlt,2,iWidth,i) = wt_MoRateBranch
                   BR_wt(iFlt,3,iWidth,i) = bValueWt1(i_bValue)
                   BR_wt(iFlt,4,iWidth,i) = 0
                   BR_wt(iFlt,5,iWidth,i) = 0
                   BR_wt(iFlt,6,iWidth,i) = 0
                   BR_wt(iFlt,7,iWidth,i) = temp_BR_wt(iRate)
                 endif
                 BR_wt(iFlt,8,iWidth,i) = magRecurWt1(iRecur)
                 BR_wt(iFlt,9,iWidth,i) = faultThickWt1(iThick1)
                 BR_wt(iFlt,10,iWidth,i) = 1.
                 BR_wt(iFlt,11,iWidth,i) = refMagWt1(iWidth,iRefMag)
     
c     End loop over iRefMag
                enddo
c     End loop over ib_value
               enddo
c     End loop over iRate
              enddo
c     End loop over iRecur
             enddo
             nParamVar(iFlt,iWidth) = i
c     End loop over iThick1
            enddo
c     End loop over iDip
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
      end
