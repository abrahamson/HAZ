      real*8 Haz(MAX_INTEN, MAX_PROB, MAX_FLT)
      real*8 magbar1(MAX_INTEN, MAX_PROB, MAX_FLT)
      real*8 tempHaz(MAXPARAM,MAX_WIDTH,MAX_INTEN, MAX_PROB, MAX_ATTEN,MAX_FTYPE)
      real*8 mnHaz(MAX_INTEN, MAX_PROB, MAX_ATTEN)
      real*8 varHaz(MAX_INTEN, MAX_PROB, MAX_ATTEN)
      
      real siteX, siteY, 
     1     fLong(MAX_FLT, MAX_DD, MAX_SEG), fLat(MAX_FLT, MAX_DD, MAX_SEG), 
     2     fZ(MAX_FLT, MAX_DD, MAX_SEG), 
     3     dip(MAX_FLT,MAX_WIDTH, MAX_SEG), lgTestInten(MAX_PROB, MAX_INTEN), lgInten,
     4     testInten(MAX_PROB, MAX_INTEN) 
     
      real scaleRate(MAX_FLT)
      integer grid_n(MAX_FLT), n_bvalue, nMagrecur, npar
      integer nSliprate(MAX_N1), nMaxMag(MAX_N1)
      real grid_a(MAX_FLT,MAX_GRID), grid_lat(MAX_FLT,MAX_GRID), grid_long(MAX_FLT,MAX_GRID)
      real grid_dlong(MAX_FLT), grid_dlat(MAX_FLT), grid_top(MAX_FLT,MAX_GRID), phypo, phypod
      real grid_x(MAX_GRID), grid_y(MAX_GRID)
      real xFlt(MAX_DD,MAX_SEG), yFlt(MAX_DD,MAX_SEG), zFlt(MAX_DD,MAX_SEG)
      real xRup(4,MAX_SEG), yRup(4,MAX_SEG), zRup(4,MAX_SEG),  dipRup(MAX_SEG)
      real xRupCamp(4,MAX_SEG), yRupCamp(4,MAX_SEG), zRupCamp(4,MAX_SEG)
      real xRupHypo(4,MAX_SEG), yRupHypo(4,MAX_SEG), zRupHypo(4,MAX_SEG)
      real dip1(MAX_SEG), dip0(MAX_SEG), strike(MAX_SEG), minDepth(MAX_FLT)
      real cumLength(MAX_SEG), coeff_area(2,MAX_FLT), sigArea(MAX_FLT),
     1     coeff_width(2,MAX_FLT), sigWidth(MAX_FLT),
     2     sigmaY, faultWidth(MAX_FLT,MAX_WIDTH),  
     3     faultWidthWt(MAX_FLT,MAX_WIDTH),
     4     maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH), minMag(MAX_FLT), 
     4     magStep(MAX_FLT), flength(MAX_SEG),
     5     beta(MAX_FLT,MAXPARAM,MAX_WIDTH), rate(MAXPARAM,MAX_WIDTH)
      real sigMaxArea, sigMaxWidth, areaStep, widthStep
      real maxMagWt(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real RateParamWt(MAX_FLT,MAXPARAM,MAX_WIDTH),
     1     RateParam(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real segModelWt(MAX_FLT)
      real MagRecur(MAX_FLT,MAXPARAM,MAX_WIDTH),
     1     magRecurWt(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real mpdf_param(MAX_FLT,MAXPARAM,MAX_WIDTH,5)
      real xStep(MAX_FLT), yStep(MAX_FLT), 
     1     maxmag1, minDist, maxInten, mag, pMag(MAXPARAM,MAX_WIDTH),
     1     mindist1, wtout(MAX_FLT,MAXPARAM,MAX_WIDTH,MAX_FTYPE), rout(MAXPARAM,MAX_WIDTH),
     1     wt2(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real*8 probX, mHaz, pxceed3
      real pxceed, faultLength, period1(4,MAX_PROB),
     1     shortDist(MAX_FLT)
      real HazBins(MAX_MAG,MAX_DIST,MAX_EPS,MAX_PROB,MAX_INTEN)
      real HazBinsX(MAX_Xcost,MAX_PROB,MAX_INTEN)
      real magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS)
      real xCostBins(MAX_Xcost), xcost, fd
      real*8 m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN)
      real*8 Xcost_bar(MAX_PROB,MAX_INTEN)
      real*8 e_bar(MAX_PROB,MAX_INTEN)
      real epsilon1 
      real sigTrunc(MAX_PROB), nu_poisson 
      real gmScale(MAX_PROB,4,MAX_ATTEN), gm_wt(MAX_PROB,4,MAX_ATTEN), sigvaradd(MAX_PROB,4,MAX_ATTEN)
      real rupArea, rupWidth, hx, hy, x0, y0, z0, probAct(MAX_FLT)
      real pWidth, pArea, pHx, pHy, distRup, distJB, distSeismo, disthypo, distepi
      real csratio, distDensity(MAX_DIST1), r_horiz
      real al_segWt(MAX_FLT), sampleStep(MAX_FLT), specT(MAX_PROB)
      real faultdist(MAX_FLT,MAX_WIDTH,3)
      integer csrflag(MAX_PROB),hwflag, psCorFlag, dirflag(MAX_PROB)
      integer tapflag(MAX_PROB), intflag(4,MAX_PROB), fwflag
      integer nInten(MAX_PROB), nfp(MAX_FLT), dirindex(MAX_PROB),
     1        nMag(MAX_FLT), jcalc1, jcalc(MAX_PROB,4,MAX_ATTEN), nAtten(MAX_PROB,4), nSite, nFlt, isite,
     3        iMag, iflt, nRupArea(MAX_FLT), nRupWidth(MAX_FLT), iAtten, iInten, nWidth(MAX_FLT),
     4        scalc(MAX_PROB,4,MAX_ATTEN), scalc1, ssscalc(MAX_PROB,4,MAX_ATTEN), ssscalc1, ssscalctemp
      integer iWidth, iArea, iHx, iHy, nSegRup
      integer hwflagout(MAX_FLT)
      integer sourceType(MAX_FLT), attenType(MAX_FLT)
      integer iParam, nParamVar(MAX_FLT,MAX_WIDTH), 
     1        iMagBin, iDistBin, iEpsBin, iXcost,
     1        nMagBins, nDistBins, nEpsBins, nPts, nHx, nHy, iCoor
      integer nXcostBins, nAttenType
      real sigfix(MAX_PROB,4,MAX_ATTEN), sigfix1
      real vs, mMagout(MAX_FLT,MAX_WIDTH,MAXPARAM)
      real mMagoutWt(MAX_FLT,MAX_WIDTH,MAXPARAM)
      integer fIndex(3,MAX_FLT), iFltWidth
      integer fsys(MAX_FLT)
      character*80 fName(MAX_FLT), attenName(4,MAX_PROB), file1, file2, file4, sigmaName(4,MAX_PROB)
      integer jInten, iSoilBin, nepsrock
      real lgRockSa, lgSoilSa, epsrock1, epsrock2
      real*8 prock
      integer fltDirect(MAX_FLT), synchron(MAX_FLT), synatten(MAX_FLT)
      integer nSyn_Case(MAX_FLT)
      real SynMag(MAX_FLT,MAX_SYN), SynDistRup(MAX_FLT,MAX_SYN)
      real synDistJB(MAX_FLT,MAX_SYN), synDistSeismo(MAX_FLT,MAX_SYN)
      real synhypo(MAX_FLT,MAX_SYN), synFtype(MAX_FLT,MAX_SYN)
      real SynWt(MAX_FLT,MAX_SYN), sigma, lnY, distRupS, distJBS
      real distSeismoS, ftypeS, lgIntenS, sigmaYS, magS, magTotal
      integer synHWFlag(MAX_FLT,MAX_SYN), hwflagS, jtypeS, iSegClose
      real refPeriod(MAX_AMPPER), RefGM_Mag(MAX_AMPMAG)
      real RefGM(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM)
      real amp(MAX_AMPMAG,MAX_AMPPER,MAX_AMPGM), depthvs10, depthvs15
      real RateType(MAX_FLT,MAXPARAM,MAX_WIDTH), D25
      integer soilAmpFlag, sigflag, attenflag
      character*80 attentitle, attenoutfile 
      real attenrupdist(200), attenjbdist(200), attenseisdist(200)
      real Theta_site          
      real charMeanMo(MAXPARAM,MAX_WIDTH), expMeanMo(MAXPARAM,MAX_WIDTH)
      integer anper, vs30_class, foreARc
      real minaper, maxaper, aper(MAX_PER), vrup
      real*8 tempTotal
      real depthParam(MAX_FLT,5), ep1, ep2, temp1, hy1, hy2, ruplen
      integer iDepthModel(MAX_FLT), nfType(MAX_FLT), nMaxMag2(MAX_FLT), nFtype1(MAX_FLT)
      real ftype(MAX_FLT,MAX_N1), ftype_wt(MAX_FLT,MAX_N1)
      real ftype1(MAX_FLT,MAX_N1), ftype_wt1(MAX_FLT,MAX_N1)
      real minlat, maxlat, minlong, maxlong
      real cfcoefrrup(MAX_ATTEN,11), cfcoefRjb(MAX_ATTEN,11)
      real fltGrid_x1(MAXFLT_DD,MAXFLT_AS), fltGrid_y1(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_z1(MAXFLT_DD,MAXFLT_AS), fltGrid_x2(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_y2(MAXFLT_DD,MAXFLT_AS), fltGrid_z2(MAXFLT_DD,MAXFLT_AS),
     3     fltGrid_x3(MAXFLT_DD,MAXFLT_AS), fltGrid_y3(MAXFLT_DD,MAXFLT_AS),
     4     fltGrid_z3(MAXFLT_DD,MAXFLT_AS), fltGrid_x4(MAXFLT_DD,MAXFLT_AS),
     5     fltGrid_y4(MAXFLT_DD,MAXFLT_AS), fltGrid_z4(MAXFLT_DD,MAXFLT_AS)
