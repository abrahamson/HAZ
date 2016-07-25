c    Declarations for Main file

      real*8 Haz(MAX_INTEN,MAX_PROB,MAX_FLT), p1_sum, wt,
     1       magbar1(MAX_INTEN,MAX_PROB,MAX_FLT), prock,
     2       tempHaz(MAXPARAM,MAX_INTEN,MAX_PROB,MAX_ATTEN,MAX_FTYPE),
     3       tempHaz1(MAXPARAM,MAX_INTEN,MAX_PROB,MAX_FTYPE),
     4       tempHaz2(4,MAX_INTEN,MAX_PROB,MAX_ATTEN), p1
      real*8 mHaz, pxceed3, d_bar(MAX_PROB,MAX_INTEN),
     1       m_bar(MAX_PROB,MAX_INTEN), e_bar(MAX_PROB,MAX_INTEN),
     2       Xcost_bar(MAX_PROB,MAX_INTEN)
      real*8 m_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN), rrup_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN),
     1       rjb_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN), rx_bar_s(MAX_FLT,MAX_PROB,MAX_INTEN),
     2       e_bar_s(MAX_FLT, MAX_PROB,MAX_INTEN)

      integer faultFlag(MAX_FLT,100,MAX_FLT), nDD(MAX_FLT), nfltGrid(2),
     1        BR_index(MAX_FLT,20,MAX_WIDTH,MAXPARAM),
     2        segModelFlag(MAX_FLT,100), nSegModel(MAX_FLT), runflag,
     3        n1AS(MAXFLT_AS), n2AS(MAXFLT_AS), icellRupStrike, icellRupDip,
     4        bnum, bnumflag, coefcountRrup, coefcountRjb, rup1_flag
      integer dirFlag1, iMixture(MAX_PROB,4,MAX_ATTEN), ibnum, nProb,
     1        nRefPer, nRefGM, nRefMag, nLocXAS, nLocYST1, nLocX,
     2        iDepthFlag, iLocX, nLocY, iLocY, n1, n2, iFtype, iProb,
     3        jType, isyn, nHypoX, iHypoX, iHypoZ, nHypoXStep,
     4        nHypoZ, nHypoZStep, dpp_flag, grid_n(MAX_FLT), iAtten
      integer csrflag(MAX_PROB),hwflag, psCorFlag, dirflag(MAX_PROB),
     1        tapflag(MAX_PROB), intflag(4,MAX_PROB), isite, nFlt, jInten,
     2        nInten(MAX_PROB), nfp(MAX_FLT), nMag(MAX_FLT), jcalc1,
     3        jcalc(MAX_PROB,4,MAX_ATTEN), nAtten(MAX_PROB,4), nSite,
     4        iMag, iflt, nRupArea(MAX_FLT), nRupWidth(MAX_FLT), nWidth(MAX_FLT)
      integer scalc(MAX_PROB,4,MAX_ATTEN), scalc1,
     1        iWidth, iArea, hwflagout(MAX_FLT), sourceType(MAX_FLT),
     2        attenType(MAX_FLT), iParam, nParamVar(MAX_FLT,MAX_WIDTH),
     3        iMagBin, iDistBin, iEpsBin, nMagBins, nDistBins, nEpsBins,
     4        nPts, iCoor, nXcostBins, nAttenType, vs30_class, foreARc
      integer fIndex(3,MAX_FLT), iFltWidth, fsys(MAX_FLT), sigflag,
     1        fltDirect(MAX_FLT), synchron(MAX_FLT), synjcalc(MAX_FLT),
     2        nSyn_Case(MAX_FLT), synHWFlag(MAX_FLT,MAX_SYN),
     3        soilAmpFlag, iDepthModel(MAX_FLT), nfType(MAX_FLT),
     4        nMaxMag2(MAX_FLT), ncountS7(MAX_FLT), nn10000

      real segWt1(MAX_FLT), fltGrid_X(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_y(MAXFLT_DD,MAXFLT_AS), fltGrid_z(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_fLen(MAXFLT_DD,MAXFLT_AS), fltGrid_w(MAXFLT_DD,MAXFLT_AS),
     3     fltGrid_a(MAXFLT_DD,MAXFLT_AS), fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS),
     4     fltGrid_RJB(MAXFLT_DD,MAXFLT_AS), dipavgd
      real Rx, Ry, Ry0, BR_wt(MAX_FLT,20,MAX_WIDTH,MAXPARAM),
     1     lgInten0, pLocY(MAXFLT_AS),
     2     sigmaTotal, sigma1, sigma2, wt1, phi, tau, distDensity2(MAX_GRID),
     3     segModelWt1(MAX_FLT,100), distmax, grid_dx, grid_dy, faultArea,
     4     faultLen, pLocX, hypoDepth, ZTOR, t1, aveWidth, probSyn
      real dirMed, dirSigma, fs, pHypoX, pHypoZ, temp, sigma0, AR,
     1     siteX, siteY, testInten(MAX_PROB, MAX_INTEN), lgInten,
     2     fLong(MAX_FLT, MAX_DD, MAX_SEG), fLat(MAX_FLT, MAX_DD, MAX_SEG),
     3     fZ(MAX_FLT, MAX_DD, MAX_SEG), scaleRate(MAX_FLT),
     4     dip(MAX_FLT,MAX_WIDTH, MAX_SEG), lgTestInten(MAX_PROB, MAX_INTEN)
      real grid_a(MAX_FLT,MAX_GRID), grid_lat(MAX_FLT,MAX_GRID),
     1     grid_long(MAX_FLT,MAX_GRID), grid_dlong(MAX_FLT),
     2     grid_dlat(MAX_FLT), grid_top(MAX_FLT,MAX_GRID), maxlat,
     3     grid_x(MAX_GRID), grid_y(MAX_GRID), xFlt(MAX_DD,MAX_SEG),
     4     yFlt(MAX_DD,MAX_SEG), zFlt(MAX_DD,MAX_SEG), minDepth(MAX_FLT)
      real coeff_area(2,MAX_FLT), sigArea(MAX_FLT), sigmaY,
     1     coeff_width(2,MAX_FLT), sigWidth(MAX_FLT), sigMaxArea,
     2     faultWidth(MAX_FLT,MAX_WIDTH), magStep(MAX_FLT), ruplen,
     3     faultWidthWt(MAX_FLT,MAX_WIDTH), sigTrunc(MAX_PROB),
     4     maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), minMag(MAX_FLT)
      real sigMaxWidth, areaStep, widthStep, segModelWt(MAX_FLT),
     1     beta(MAX_FLT,MAXPARAM,MAX_WIDTH), rate(MAXPARAM,MAX_WIDTH),
     2     maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH), maxmag1, mag,
     3     RateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH), minDist,
     4     RateParam(MAX_FLT,MAXPARAM,MAX_WIDTH), maxInten, sigfix1
      real MagRecur(MAX_FLT,MAXPARAM,MAX_WIDTH), xStep(MAX_FLT),
     1     magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH), wt2(MAXPARAM),
     2     mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,5), epsilon1,
     3     yStep(MAX_FLT), pMag(MAXPARAM,MAX_WIDTH), xcost, fd,
     4     wtout(MAX_FLT,MAXPARAM,MAX_WIDTH,MAX_FTYPE), vs, vrup
      real rout(MAXPARAM,MAX_WIDTH), period1(4,MAX_PROB), D25,
     1     MinRrup(MAX_FLT), xCostBins(MAX_Xcost), r_horiz, MinSeismo_temp,
     2     HazBinsX(MAX_Xcost,MAX_PROB,MAX_INTEN), rupArea, minlat,
     3     HazBins(MAX_MAG,MAX_DIST,MAX_EPS,MAX_PROB,MAX_INTEN),
     4     magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS)
      real gmScale(MAX_PROB,4,MAX_ATTEN), gm_wt(MAX_PROB,4,MAX_ATTEN),
     1     sigvaradd(MAX_PROB,4,MAX_ATTEN), rupWidth, x0, y0, z0,
     2     probAct(MAX_FLT), pWidth, pArea, distRup, distJB,
     3     distSeismo, disthypo, distepi, distDensity(MAX_DIST1),
     4     al_segWt(MAX_FLT), sampleStep(MAX_FLT), specT(MAX_PROB)
      real SourceDist(MAX_FLT,MAX_WIDTH,3), RefGM_Mag(MAX_AMPMAG),
     1     sigfix(MAX_PROB,4,MAX_ATTEN), lgIntenS, magTotal,
     2     mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM), depthvs10,
     3     mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM), depthvs15,
     4     SynMag(MAX_FLT,MAX_SYN), SynDistRup(MAX_FLT,MAX_SYN)
      real synDistJB(MAX_FLT,MAX_SYN), synDistSeismo(MAX_FLT,MAX_SYN),
     1     synhypo(MAX_FLT,MAX_SYN), synFtype(MAX_FLT,MAX_SYN),
     2     SynWt(MAX_FLT,MAX_SYN), syn_dip(MAX_FLT,MAX_SYN),
     3     syn_zTOR(MAX_FLT,MAX_SYN), syn_RupWidth(MAX_FLT,MAX_SYN),
     4     syn_RX(MAX_FLT,MAX_SYN), syn_Ry0(MAX_FLT,MAX_SYN)
      real refPeriod(MAX_AMPPER), RefGM(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM),
     1     amp(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM), Theta_site, MinRjb_temp,
     2     RateType(MAX_FLT,MAXPARAM,MAX_WIDTH), minlong, maxlong,
     3     charMeanMo(MAXPARAM,MAX_WIDTH), expMeanMo(MAXPARAM,MAX_WIDTH),
     4     depthParam(MAX_FLT,5), ftype(MAX_FLT,MAX_N1), ftype_wt(MAX_FLT,MAX_N1)
      real fltGrid_x1(MAXFLT_DD,MAXFLT_AS), fltGrid_y1(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_z1(MAXFLT_DD,MAXFLT_AS), fltGrid_x2(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_y2(MAXFLT_DD,MAXFLT_AS), fltGrid_z2(MAXFLT_DD,MAXFLT_AS),
     3     fltGrid_x3(MAXFLT_DD,MAXFLT_AS), fltGrid_y3(MAXFLT_DD,MAXFLT_AS),
     4     fltGrid_z3(MAXFLT_DD,MAXFLT_AS), fltGrid_x4(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_y4(MAXFLT_DD,MAXFLT_AS), fltGrid_z4(MAXFLT_DD,MAXFLT_AS),
     1     cfcoefrrup(MAX_ATTEN,11), cfcoefRjb(MAX_ATTEN,11),
     2     magS7(MAX_FLT,MAX_S7), rateS7(MAX_FLT,MAX_S7),
     3     distS7(MAX_FLT,MAX_S7), DipS7(MAX_FLT,MAX_S7),
     4     mechS7(MAX_FLT,MAX_S7)

      character*80 fName(MAX_FLT), attenName(4,MAX_PROB), file1, file2,
     1             sigmaName(4,MAX_PROB), filebmode
