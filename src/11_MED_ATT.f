
      subroutine S11_meanInten (rupdist, jbdist, seismodist, hwflag,
     1  mag, jcalc, specT, lnY, siga, ftype, attenName, period1, iAtten, iProb,
     2  jType, vs, depth,intflag, AR, dipavgd, disthypo, depthvs10,
     3  depthvs15, D25, tau, depthTop, Theta_Site, RupWidth, vs30_class,
     4  foreArc, Rx, phi, cfcoefrrup, cfcoefrjb, Ry0 )

      implicit none
      include 'pfrisk.h'

      real mag, lnY, siga, ftype, period1(4,MAX_PROB), mag1, period2, lnSa, m
      integer jcalc, soilflag, softrock, hardRock, hwflag, imod
      real jbdist, rupdist, seismodist, baseDepth, depth, specT, factor, a, b, faddmag
      character*80 attenName(4,MAX_ATTEN), attenname1, attenname0
      character*10 number
      integer intflag(4,MAX_PROB), iflag, vs30_class, iBranch
      integer iflag01, iflag02, iflag04, iflag10, foreArc, regionflag, basinflag, msasflag
      integer coefcountRrup, coefcountRjb, iAtten, iProb, jType, region

      real sc, sd, se, D25, RupWidth, Sa, Ss, Sr, Q0, SCa, SCb, SCc, SCd, SCe
      real magc, rupdistc, ftypec, lnYc, depthTop, phiSS
      real disthypo, depthvs10, depthvs15, AR, dipavgd
      real Theta_Site, sigmaEps, facFena, sdscale, BA08lnY, deltaC1
      real lnY01, lnY02, lnY04, lnY10, lnY02p, lnY04p
      real sigma01, sigma02, sigma04, sigma10
      real specT01, specT02, specT04, specT10, period02, Rx, SF2, Ry0
      real phi, tau, lnYH, sourceclass, sigmac, fth, frv, vs, sigma
      real cfcoefrrup(MAX_Atten,11), cfcoefrjb(MAX_Atten,11), c1
      real depthTop1, s03vfs, s03sr, s03fr, soil, GB, GC, pga4nl
      real sigmaH, phiH, tauH, sclass, pgaref, sjb

C LNY IS EXPECTED INTENSITY FOR THIS MAGNITUDE AND CLOSEST DISTANCE
      lnY = 1.e30
      tau = 0.0
      phi = 0.0
      iflag = 0

c  *** Turkey adjusted NGA1 models
      if (jcalc .eq. 9501 ) then
          call S08_AS_NGA_2008TR ( mag, dipavgd, fType, RupWidth, rupDist, jbdist,
     1                     vs, hwflag, lnY, sigma,
     2                     specT, period2, depthtop, iflag, vs30_class,
     3                     depthvs10, Rx )
      endif

      if (jcalc .eq. 9502 ) then
         call S08_BA_NGA_2008TR ( mag, jbdist, specT, period2, lnY, sigma, iflag, vs, ftype, pga4nl )
      endif

      if (jcalc .eq. 9503 ) then
         call S08_CB_NGA_2008TR ( mag, rupDist, jbdist, Ftype, specT,
     1                     period2, lnY, sigma, iflag, vs,
     2                     depthtop, D25, dipavgd )
      endif

      if (jcalc .eq. 9504) then
          call S08_CY_NGA_2008TR ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, depthTop, Ftype, depthvs10, vs30_class,
     3                     hwflag, Rx )
      endif

      if (jcalc .eq. 9505 ) then
            call S08_I_NGA_2008TR ( mag, rupDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
      endif
      if (jcalc .eq. 9506) then
          if ( depthTop .gt. 20. ) then
             depthTop1 = 20.
          else
             depthTop1 = depthTop
          endif
          call S08_CY_NGA_2008TR ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, depthTop1, Ftype, depthvs10, vs30_class,
     3                     hwflag, Rx )
      endif

C     Kale, Ö., Akkar, S., Ansari, A., & Hamzehloo, H. (2015).
C     -A Ground‐Motion Predictive Model for Iran and Turkey for Horizontal PGA,
C     PGV, and 5% Damped Response Spectrum: Investigation of Possible Regional Effects
C     BSSA: Vol. 105, No.2A, pp. 963-980.
C     Model Number = 2150
      if ( jcalc .eq. 2150 ) then
C     For Turkey
        region = 0

         call S02_KAAH_2015 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pgaref,region )
         attenname1 = 'Kale et AL. (2015)_Hor_Turkey'
       endif

      if ( jcalc .eq. 3150 ) then
C     For Iran
        region = 1

         call S02_KAAH_2015 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pgaref,region )
         attenname1 = 'Kale et AL. (2015)_Hor_Iran'
      endif

c    Taiwan GMPEs


C     **** Taiwan TNGA attenuation models  ****
c-----------------------------------------------------------------------------
c     New Lin et al. (2011) Taiwan crustal model SOIL ************
C     Model Number = 316
      if ( jcalc .eq. 316 ) then
      	if (hwflag .eq. 1 ) then
         call S04_Lin_hw_soil
     1    ( mag, rupDist, specT, period2, lnY, sigma, iflag)
        else
         call S04_Lin_fw_soil
     1    ( mag, rupDist, specT, period2, lnY, sigma, iflag)
        endif
         attenname1 = 'Lin et al. (2011) , Crustal soil'
      endif

c-----------------------------------------------------------------------------
c     New Lin et al. (2011) Taiwan crustal model ROCK ************
C     Model Number = 317
      if ( jcalc .eq. 317 ) then
      	if (hwflag .eq. 1 ) then
         call S04_Lin_hw_rock
     1    ( mag, rupDist, specT, period2, lnY, sigma, iflag)
        else
         call S04_Lin_fw_rock
     1    ( mag, rupDist, specT, period2, lnY, sigma, iflag)
        endif
         attenname1 = 'Lin et al. (2011) , Crustal rock'
      endif

c-----------------------------------------------------------------------------
c     Lin 2009 Doctoral thesis
C     Model Number = 315
      if ( jcalc .eq. 315) then
         call S04_Lin2009
     1( mag, rupDist, specT, period2, lnY, sigma, vs, iflag, ftype)
         attenname1 = 'Lin 2009 Doctoral thesis, Crustal, VS30'
      endif
c-----------------------------------------------------------------------------
c     TG09221 2012 project
C     Model Number = 441
      if ( jcalc .eq. 441) then
         call S04_TG09221_2012
     1( mag, rupDist, specT, period2, lnY, sigma, vs, iflag, ftype)
         attenname1 = 'TG09221 Report 2012/06, Crustal, VS30'
      endif

c-----------------------------------------------------------------------------
c     NCREE 2011 project
C     Model Number = 451
      if ( jcalc .eq. 451) then
         call S04_NCREE_2011
     1    ( mag, rupDist, specT, period2, lnY, sigma)
         attenname1 = 'NCREE Report 2011/01, Vs30°Ÿ360m/sec'
      endif

C     **** End of Taiwan TNGA attenuation models  ****

C ******* PEER NGA Attenuation Models ****
c ******* Abrahamson and Silva Model *********
C     Abrahamson&Silva 2008 - horizontal, Estimated Vs30m
C     Model Number = 787
      if ( jcalc .eq. 787 ) then
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S07_AS_NGA_2008 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, phi, tau )
         attenname1 = 'A&S_NGA_2008-Hor,Estimated Vs30m'
      endif

C     Abrahamson&Silva 2008 - horizontal, Measured Vs30m
C     Model Number = 788
      if ( jcalc .eq. 788 ) then
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S07_AS_NGA_2008 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, phi, tau )
         attenname1 = 'A&S_NGA_2008-Hor,Measured Vs30m'
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2008 - Horizontal, estimated Vs30m
C     Model Number = 797
c
      if ( jcalc .eq. 797 ) then
         vs30_class = 0
         call S07_CY_NGA_2008 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, phi, tau )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Estimated Vs30m'
       endif

C     Chiou and Youngs 2008 - Horizontal, measured Vs30m
C     Model Number = 798
c
      if ( jcalc .eq. 798 ) then
         vs30_class = 1
         call S07_CY_NGA_2008 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, phi, tau )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Measured Vs30m'
       endif

c ******* Chiou and Youngs Model - Small Magnitude Models *********
C     Chiou and Youngs 2008 - Horizontal, estimated Vs30m
C     Southern California Small Magnitude Model (2010)
C     Model Number = 799
c
      if ( jcalc .eq. 799 ) then
         vs30_class = 0
         call S07_CY_NGA_2008SC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Estimated Vs30m,SCal SMM'
       endif

C     Chiou and Youngs 2008 - Horizontal, measured Vs30m
C     Southern California Small Magnitude Model (2010)
C     Model Number = 800
c
      if ( jcalc .eq. 800 ) then
         vs30_class = 1
         call S07_CY_NGA_2008SC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Measured Vs30m,SCal SMM'
       endif

C     Chiou and Youngs 2008 - Horizontal, estimated Vs30m
C     Central California Small Magnitude Model (2010)
C     Model Number = 801
c
      if ( jcalc .eq. 801 ) then
         vs30_class = 0
         call S07_CY_NGA_2008CC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Estimated Vs30m,CCal SMM'
       endif

C     Chiou and Youngs 2008 - Horizontal, measured Vs30m
C     Central California Small Magnitude Model (2010)
C     Model Number = 802
c
      if ( jcalc .eq. 802 ) then
         vs30_class = 1
         call S07_CY_NGA_2008CC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Measured Vs30m,CCal SMM'
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2008 - horizontal
C     Model Number = 836
      if ( jcalc .eq. 836 ) then
         call S07_CB_NGA_2008 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, phi, tau )
         attenname1 = 'Campbell&Bozorgnia_NGA_2008-Hor'
       endif

c ******* Idriss Model *********
C     Idriss 2008 - Horizontal, Vs=450-900 m/sec and Vs>900 m/sec
C     Model Number = 910
      if ( jcalc .eq. 910 ) then
         if (vs .ge. 900.0) then
            call S07_I_NGA_2008vsgt900 ( mag, rupDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGA_2008_Hor,Vs>900m/s'
         else
            call S07_I_NGA_2008 ( mag, rupDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGA_2008_Hor,Vs=450-900m/s'
         endif
       endif

c ******* Boore and Atkinson Model *********
C     Boore and Atkinson July 2008 - horizontal
C     Model Number = 922
      if ( jcalc .eq. 922 ) then
         call S07_BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pga4nl, phi, tau )
         attenname1 = 'Boore&Atkinson_NGA_2008_Hor'
       endif

c ******* Boore and Atkinson Model *********
C     Boore and Atkinson July 2008 - horizontal with Atkinson (2010) small magnitude adjustment
C     Model Number = 923
      if ( jcalc .eq. 923 ) then
         call S07_BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pga4nl, phi, tau )
C     Apply small magnitude adjustment is Mag<=5.75.
         if (mag .le. 5.75 ) then
            a = max(0.0, 3.888 - 0.674*mag)
            b = max(0.0, 2.933 - 0.510*mag)
            factor = a - b*log10(jbdist+10.0)
            factor = factor*alog(10.0)
         endif
         attenname1 = 'Boore&Atkinson_NGA_2008_Hor with small mag adj Atkinson (2010)'
         lnY = lnY + factor
       endif

C ******* End of PEER NGA Attenuation Models ****

C ******* PEER NGA-West2 Attenuation Models ****
c ******* Abrahamson, Silva, and Kamai Model *********
C     Note: GMPE is not programmed for Aftershock cases.
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Estimated Vs30m

C     Model Number = 2787
      if ( jcalc .eq. 2787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-EstVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Estimated Vs30m
C     Model Number = 2788
      if ( jcalc .eq. 2788 ) then
         regionflag = 1
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-EstVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Estimated Vs30m
C     Model Number = 2789
      if ( jcalc .eq. 2789 ) then
         regionflag = 2
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-EstVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Estimated Vs30m
C     Model Number = 2790
      if ( jcalc .eq. 2790 ) then
         regionflag = 3
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-EstVs'
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Measured Vs30m
C     Model Number = 2791
      if ( jcalc .eq. 2791 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-MesVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Measured Vs30m
C     Model Number = 2792
      if ( jcalc .eq. 2792 ) then
         regionflag = 1
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-MesVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Measured Vs30m
C     Model Number = 2793
      if ( jcalc .eq. 2793 ) then
         regionflag = 2
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-MesVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Measured Vs30m
C     Model Number = 2794
      if ( jcalc .eq. 2794 ) then
         regionflag = 3
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-MesVs'
      endif

C     Note: These calls are for Aftershocks cases but it not currently implemented.
C           For implementation new distance term CRjb will need to be computed and passed.
c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Aftershock, Estimated Vs30m
C     Model Number = 3787
      if ( jcalc .eq. 3787 ) then
         regionflag = 0
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-AS-EstVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Aftershock, Estimated Vs30m
C     Model Number = 3788
      if ( jcalc .eq. 3788 ) then
         regionflag = 1
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-AS-EstVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Aftershock, Estimated Vs30m
C     Model Number = 3789
      if ( jcalc .eq. 3789 ) then
         regionflag = 2
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-AS-EstVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Aftershock, Estimated Vs30m
C     Model Number = 3790
      if ( jcalc .eq. 3790 ) then
         regionflag = 3
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-AS-EstVs'
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Aftershock, Measured Vs30m
C     Model Number = 3791
      if ( jcalc .eq. 3791 ) then
         regionflag = 0
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-AS-MesVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Aftershock, Measured Vs30m
C     Model Number = 3792
      if ( jcalc .eq. 3792 ) then
         regionflag = 1
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-AS-MesVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Aftershock, Measured Vs30m
C     Model Number = 3793
      if ( jcalc .eq. 3793 ) then
         regionflag = 2
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-AS-MesVs'
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Aftershock, Mesaured Vs30m
C     Model Number = 3794
      if ( jcalc .eq. 3794 ) then
         regionflag = 3
         msasflag = 1
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-AS-MesVs'
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, estimated Vs30
C     Model Number = 2797
      if ( jcalc .eq. 2797 ) then
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Estimated Vs30m'
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, measured Vs30
C     Model Number = 2798
      if ( jcalc .eq. 2798 ) then
         vs30_class = 1
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Measured Vs30m'
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment, estimated Vs30
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 2799
      if ( jcalc .eq. 2799 ) then
         vs30_class = 0
         regionflag = 1
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Estimated Vs30m'
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment, measured Vs30
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 2800
      if ( jcalc .eq. 2800 ) then
         vs30_class = 1
         regionflag = 1
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Measured Vs30m'
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment, estimated Vs30
C      Note: Only valid for the Following magnitude: M=7.9
C     Model Number = 2801
      if ( jcalc .eq. 2801 ) then
         vs30_class = 0
         regionflag = 2
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchuan,Estimated Vs30m'
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment, measured Vs30
C      Note: Only valid for the Following magnitude: M=7.9
C     Model Number = 2802
      if ( jcalc .eq. 2802 ) then
         vs30_class = 1
         regionflag = 2
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchaun,Measured Vs30m'
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, California
C     Model Number = 2836
      if ( jcalc .eq. 2836 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Calif'
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Japan
C     Model Number = 2837
      if ( jcalc .eq. 2837 ) then
         regionflag = 1
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Japan'
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, China
C     Model Number = 2838
      if ( jcalc .eq. 2838 ) then
         regionflag = 2
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,China'
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Italy
C     Model Number = 2839
      if ( jcalc .eq. 2839 ) then
         regionflag = 3
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Italy'
       endif


c ******* Idriss Model *********
C     Idriss 2013 - Horizontal
C     Model Number = 2910
      if ( jcalc .eq. 2910 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 2922
      if ( jcalc .eq. 2922 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, No Basin'
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, No Basin Adjustment
C     Model Number = 2923
      if ( jcalc .eq. 2923 ) then
         regionflag = 1
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, No Basin'
       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, No Basin Adjustment
C     Model Number = 2924
      if ( jcalc .eq. 2924 ) then
         regionflag = 2
         Basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, No Basin'
       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, Basin Adjustments
C     Model Number = 2925
      if ( jcalc .eq. 2925 ) then
         regionflag = 0
         basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, Basin'
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, Basin Adjustment
C     Model Number = 2926
      if ( jcalc .eq. 2926 ) then
         regionflag = 1
         basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, Basin'
       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, Basin Adjustment
C     Model Number = 2927
      if ( jcalc .eq. 2927 ) then
         regionflag = 2
         Basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, Basin'
       endif

C ******* Implementation of Al-Atik and Youngs 2014 NGA West2 Epistemic Model.
C ******* Low branch models. (jcalc values in the 5000s).

c ******* Abrahamson, Silva, and Kamai Model *********
C     Note: GMPE is not programmed for Aftershock cases.
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Estimated Vs30m

C     Model Number = 5787
      if ( jcalc .eq. 5787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-EstVs, LowEps'

C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Estimated Vs30m
C     Model Number = 5788
      if ( jcalc .eq. 5788 ) then
         regionflag = 1
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-EstVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Estimated Vs30m
C     Model Number = 5789
      if ( jcalc .eq. 5789 ) then
         regionflag = 2
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-EstVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Estimated Vs30m
C     Model Number = 5790
      if ( jcalc .eq. 5790 ) then
         regionflag = 3
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-EstVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Measured Vs30m
C     Model Number = 5791
      if ( jcalc .eq. 5791 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-MesVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Measured Vs30m
C     Model Number = 5792
      if ( jcalc .eq. 5792 ) then
         regionflag = 1
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-MesVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Measured Vs30m
C     Model Number = 5793
      if ( jcalc .eq. 5793 ) then
         regionflag = 2
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-MesVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Measured Vs30m
C     Model Number = 5794
      if ( jcalc .eq. 5794 ) then
         regionflag = 3
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-MesVs, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal
C     Model Number = 5797
      if ( jcalc .eq. 5797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Estimated Vs30m, LowEps'
C     Compute the low bramch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal
C     Model Number = 5798
      if ( jcalc .eq. 5798 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Measured Vs30m, LowEps'
C     Compute the low bramch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 5799
      if ( jcalc .eq. 5799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 1
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Estimated Vs30m, LowEps'
C     Compute the low bramch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 5800
      if ( jcalc .eq. 5800 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 1
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Measured Vs30m, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment
C      Note: Only valid for the Following magnitude : 6<M<6.9
C     Model Number = 5801
      if ( jcalc .eq. 5801 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 2
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchuan,Estimated Vs30m, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment
C      Note: Only valid for the Following magnitude: M=7.9
C     Model Number = 5802
      if ( jcalc .eq. 5802 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 2
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchaun,Measured Vs30m, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, California
C     Model Number = 5836
      if ( jcalc .eq. 5836 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Calif, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Japan
C     Model Number = 5837
      if ( jcalc .eq. 5837 ) then
         regionflag = 1
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Japan, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, China
C     Model Number = 5838
      if ( jcalc .eq. 5838 ) then
         regionflag = 2
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,China, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Italy
C     Model Number = 5839
      if ( jcalc .eq. 5839 ) then
         regionflag = 3
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Italy, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif


c ******* Idriss Model *********
C     Idriss 2013 - Horizontal
C     Model Number = 5910
      if ( jcalc .eq. 5910 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor, LowEps'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor, LowEps'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 5922
      if ( jcalc .eq. 5922 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, No Basin, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, No Basin Adjustment
C     Model Number = 5923
      if ( jcalc .eq. 5923 ) then
         regionflag = 1
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, No Basin, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, No Basin Adjustment
C     Model Number = 5924
      if ( jcalc .eq. 5924 ) then
         regionflag = 2
         Basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, No Basin, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, Basin Adjustments
C     Model Number = 5925
      if ( jcalc .eq. 5925 ) then
         regionflag = 0
         basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, Basin, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, Basin Adjustment
C     Model Number = 5926
      if ( jcalc .eq. 5926 ) then
         regionflag = 1
         basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, Basin, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, Basin Adjustment
C     Model Number = 5927
      if ( jcalc .eq. 5927 ) then
         regionflag = 2
         Basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, Basin, LowEps'
C     Compute the low branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*0.083
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY - 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY - 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif


C ******* High branch models. (jcalc values in the 7000s).

c ******* Abrahamson, Silva, and Kamai Model *********
C     Note: GMPE is not programmed for Aftershock cases.
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Estimated Vs30m

C     Model Number = 7787
      if ( jcalc .eq. 7787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-EstVs, HighEps'

C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Estimated Vs30m
C     Model Number = 7788
      if ( jcalc .eq. 7788 ) then
         regionflag = 1
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-EstVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Estimated Vs30m
C     Model Number = 7789
      if ( jcalc .eq. 7789 ) then
         regionflag = 2
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-EstVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Estimated Vs30m
C     Model Number = 7790
      if ( jcalc .eq. 7790 ) then
         regionflag = 3
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-EstVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Measured Vs30m
C     Model Number = 7791
      if ( jcalc .eq. 7791 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-MesVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Measured Vs30m
C     Model Number = 7792
      if ( jcalc .eq. 7792 ) then
         regionflag = 1
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-MesVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Measured Vs30m
C     Model Number = 7793
      if ( jcalc .eq. 7793 ) then
         regionflag = 2
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-MesVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Measured Vs30m
C     Model Number = 7794
      if ( jcalc .eq. 7794 ) then
         regionflag = 3
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau )
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-MesVs, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal
C     Model Number = 7797
      if ( jcalc .eq. 7797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Estimated Vs30m, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal
C     Model Number = 7798
      if ( jcalc .eq. 7798 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Measured Vs30m, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 7799
      if ( jcalc .eq. 7799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 1
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Estimated Vs30m, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 7800
      if ( jcalc .eq. 7800 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 1
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Measured Vs30m, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment
C      Note: Only valid for the Following magnitude : 6<M<6.9
C     Model Number = 7801
      if ( jcalc .eq. 7801 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 2
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchuan,Estimated Vs30m, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment
C      Note: Only valid for the Following magnitude: M=7.9
C     Model Number = 7802
      if ( jcalc .eq. 7802 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 2
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchaun,Measured Vs30m, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, California
C     Model Number = 7836
      if ( jcalc .eq. 7836 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Calif, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Japan
C     Model Number = 7837
      if ( jcalc .eq. 7837 ) then
         regionflag = 1
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Japan, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, China
C     Model Number = 7838
      if ( jcalc .eq. 7838 ) then
         regionflag = 2
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,China, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Italy
C     Model Number = 7839
      if ( jcalc .eq. 7839 ) then
         regionflag = 3
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Italy, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif


c ******* Idriss Model *********
C     Idriss 2013 - Horizontal
C     Model Number = 7910
      if ( jcalc .eq. 7910 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor, HighEps'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor, HighEps'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 7922
      if ( jcalc .eq. 7922 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, No Basin, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, No Basin Adjustment
C     Model Number = 7923
      if ( jcalc .eq. 7923 ) then
         regionflag = 1
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, No Basin, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, No Basin Adjustment
C     Model Number = 7924
      if ( jcalc .eq. 7924 ) then
         regionflag = 2
         Basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, No Basin, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, Basin Adjustments
C     Model Number = 7925
      if ( jcalc .eq. 7925 ) then
         regionflag = 0
         basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, Basin, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, Basin Adjustment
C     Model Number = 7926
      if ( jcalc .eq. 7926 ) then
         regionflag = 1
         basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, Basin, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, Basin Adjustment
C     Model Number = 7927
      if ( jcalc .eq. 7927 ) then
         regionflag = 2
         Basinflag = 1
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, Basin, HighEps'
C     Compute the high branch epistemic adjustment for the median motions
         if (ftype .ge. 0.0) then
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*0.083
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT))
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT))
               endif
            endif
         else
            if (specT .le. 1.0) then
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.038)
               endif
            else
               if (mag .lt. 7.0) then
                  lnY = lnY + 1.645*(0.083 + 0.0171*alog(specT) + 0.038)
               else
                  lnY = lnY + 1.645*(0.056*(mag-7.0) + 0.083 + 0.0171*alog(specT) + 0.038)
               endif
            endif
         endif

       endif

C **** NGA West2 Vertical Models *********
c ******* Stewart, Seyhan, Boore, and Atkinson Model *********
C     Stewart, Seyhan, Boore, and Atkinson 2013 - vertical
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 4922
      if ( jcalc .eq. 4922 ) then
         regionflag = 0
         basinflag = 0
         call S10_SSBA_NGAWest2_2013_Vert ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, phi, tau )
         attenname1 = 'SSBA_NGAWest2_2013_Ver, DC3Global, No Basin'
       endif

C     Stewart, Seyhan, Boore, and Atkinson 2013 - V/H Ratio
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 6922
      if ( jcalc .eq. 6922 ) then
         regionflag = 0
         basinflag = 0
         call S10_SSBA_NGAWest2_2013_Vert ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, phi, tau )

         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnYH, sigmaH, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, phiH, tauH )
         attenname1 = 'SSBA_NGAWest2_2013_V/H, DC3Global, No Basin'
C     Divide Vertical by horizontal to get V/H ratio
         lnY = lnY - lnYH + 6.89
       endif

c ******* Gulerce, Abrahamson, Silva, and Kamai Model *********
C     Gulerce, Abrahamson, Silva, and Kamai 2013 - Vertical, Global, Mainshock, Estimated Vs30m

C     Model Number = 4787
      if ( jcalc .eq. 4787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S10_GKAS_NGAWest2_2013_Vert ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'GKAS_NGAWest2_2013-Ver-Glob-MS-EstVs'
      endif

C     Gulerce, Abrahamson, Silva, and Kamai 2013 - Vertical, Global, Mainshock, Measured Vs30m

C     Model Number = 4788
      if ( jcalc .eq. 4788 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S10_GKAS_NGAWest2_2013_Vert ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'GKAS_NGAWest2_2013-Ver-Glob-MS-MsrVs'
      endif

c ******* Gulerce, Abrahamson, Silva, and Kamai Model *********
C     Gulerce, Abrahamson, Silva, and Kamai 2013 - V/H Ratio, Global, Mainshock, Estimated Vs30m

C     Model Number = 6787
      if ( jcalc .eq. 6787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S10_GKAS_NGAWest2_2013_Vert ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnYH, sigmaH, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phiH, tauH)

C     Divide Vertical by Horizontal to get V/H Ratio
         lnY = lnY - lnYH + 6.89
         attenname1 = 'GKAS_NGAWest2_2013-V/H-Glob-MS-EstVs'
      endif

c ******* Gulerce, Abrahamson, Silva, and Kamai Model *********
C     Gulerce, Abrahamson, Silva, and Kamai 2013 - V/H Ratio, Global, Mainshock, Measured Vs30m

C     Model Number = 6787
      if ( jcalc .eq. 6787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 1
         call S10_GKAS_NGAWest2_2013_Vert ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)

         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnYH, sigmaH, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phiH, tauH)

C     Divide Vertical by Horizontal to get V/H Ratio
         lnY = lnY - lnYH + 6.89
         attenname1 = 'GKAS_NGAWest2_2013-V/H-Glob-MS-MsrVs'
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Vertical, estimated Vs30
C     Model Number = 4797
      if ( jcalc .eq. 4797 ) then
         vs30_class = 0
         regionflag = 0
         call S10_CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Ver,Estimated Vs30m'
       endif

C     Chiou and Youngs 2013 - Vertical, measured Vs30
C     Model Number = 4798
      if ( jcalc .eq. 4798 ) then
         vs30_class = 0
         regionflag = 0
         call S10_CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Ver,Measured Vs30m'
       endif

C     V/H Model, estimated Vs30
C     Model Number = 6797
      if ( jcalc .eq. 6797 ) then
         vs30_class = 0
         regionflag = 0
         call S10_CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
C     Call the horizontal model
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnYH, sigmaH, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phiH, tauH )

         attenname1 = 'CY_NGAWest2_2013-V/H,Estimated Vs30m'

C     Divide vertical by horizontal to get V/H Ratio.
         lnY = lnY - lnYH + 6.89
       endif

C     V/H Model, measured Vs30
C     Model Number = 6798
      if ( jcalc .eq. 6798 ) then
         vs30_class = 1
         regionflag = 0
         call S10_CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
C     Call the horizontal model
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnYH, sigmaH, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phiH, tauH )

         attenname1 = 'CY_NGAWest2_2013-V/H,Measured Vs30m'

C     Divide vertical by horizontal to get V/H Ratio.
         lnY = lnY - lnYH + 6.89
       endif

c ******* Campbell and Bozorgnia Model *********
C     Bozorgnia and Campbell 2013 - Vertical, California
C     Model Number = 4836
      if ( jcalc .eq. 4836 ) then
         regionflag = 0
         call S10_BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )
         attenname1 = 'BC_NGAWest2_2013-Ver,Calif'
       endif

C     Bozorgnia and Campbell 2013 - Vertical, Japan
C     Model Number = 4837
      if ( jcalc .eq. 4837 ) then
         regionflag = 1
         call S10_BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )
         attenname1 = 'BC_NGAWest2_2013-Ver,Japan'
       endif

C     Bozorgnia and Campbell 2013 - Implied Vertical/Horizontal Ratio, California
C     Model Number = 6836
      if ( jcalc .eq. 6836 ) then
         regionflag = 0
         call S10_BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )

         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnYH, sigmaH, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag, phiH, tauH )

C     Divide vertical by horizontal to get V/H Ratio.
         lnY = lnY - lnYH + 6.89
         attenname1 = 'BC_NGAWest2_2013-V/H,Calif'
       endif
C     Bozorgnia and Campbell 2013 - Implied Vertical/Horizontal Ratio, Japan
C     Model Number = 6837
      if ( jcalc .eq. 6837 ) then
         regionflag = 1
         call S10_BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )

         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnYH, sigmaH, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag, phiH, tauH )

C     Divide vertical by horizontal to get V/H Ratio.
         lnY = lnY - lnYH + 6.89
         attenname1 = 'BC_NGAWest2_2013-V/H,Japan'
       endif

C ******* End of vertical NGA West2 GMPEs.*******



C ******* End of PEER NGA-West2 Attenuation Models ****

C ***************************************************

c     Campbell and Bozorgnia (2003), Horizontal
C              StrikeSlip and Reverse Events, Firm Soil
C     Notes: This site condition roughly corresponds to NEHRP D.
C            Sources with Ftype=1 are classified as Reverse events.
C            HW sites are determined for those faults with dip angles
C            less than 80 degrees although the Campbell and Bozorgnia
C            relationship has a cut off dip angle of 70 degrees. The
C            80 degree angle is taken from the Abrahamson and Silva
C            relationship. Valid distance range is 0 - 60 km, with an
C            acceptable extrapolation to 100 km. For distances greater
C            than 100 km caution should be used.
C            This relationship allows for the separation of groun motions
C            from Reverse and Thrust events. Reverse events are those faults
C            with dip angles greater than or equal to 45 degrees and Thrust
C            events are those faults with dip angles less than 45 degrees.
C            The authors do note that this differnce in ground motions may
C            be statistically insignificant. Combining the Reverse and
C            Thrust events is modeled below by taking Frv=Fth=0.5.

C     Campbell and Bozorgnia (2003), Horizontal, Firm Soil, SS and Reverse
C     Model Number = 070
      if ( jcalc .eq. 70 ) then
c     Set site condition = Firm Soil (Vs~298m/s)
         s03vfs = 0.0
         s03sr  = 0.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_H ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Firm Soil'
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Very Firm Soil, SS and Reverse
C     Model Number = 071
      if ( jcalc .eq. 71 ) then
c     Set site condition = Very Firm Soil (Vs~368m/s)
         s03vfs = 1.0
         s03sr  = 0.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_H ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Very Firm Soil'
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Soft Rock, SS and Reverse
C     Model Number = 072
      if ( jcalc .eq. 72 ) then
c     Set site condition = Soft Rock (Vs~421)
         s03vfs = 0.0
         s03sr  = 1.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_H ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Soft Rock'
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Firm Rock, SS and Reverse
C     Model Number = 073
      if ( jcalc .eq. 73 ) then
c     Set site condition = Firm Rock (Vs~830m/s)
         s03vfs = 0.0
         s03sr  = 0.0
         s03fr  = 1.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_H ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Firm Rock'
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Generic Rock, SS and Reverse
C     Model Number = 074
      if ( jcalc .eq. 74 ) then
c     Set site condition = Generic Rock (Vs~620m/s)
         s03vfs = 0.0
         s03sr  = 0.5
         s03fr  = 0.5
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_H ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,GenericRock'
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Generic Soil, SS and Reverse
C     Model Number = 075
      if ( jcalc .eq. 75 ) then
c     Set site condition = Generic Soil (Vs~310m/s)
         s03vfs = 0.25
         s03sr  = 0.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_H ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,GenericSoil'
      endif

C     Campbell and Bozorgnia (2003), Vertical, Firm Soil, SS and Reverse
C     Model Number = 076
      if ( jcalc .eq. 76 ) then
c     Set site condition = Firm Soil (Vs~298m/s)
         s03vfs = 0.0
         s03sr  = 0.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_V ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Firm Soil'
      endif

C     Campbell and Bozorgnia (2003), Vertical, Very Firm Soil, SS and Reverse
C     Model Number = 077
      if ( jcalc .eq. 77 ) then
c     Set site condition = Very Firm Soil (Vs~368m/s)
         s03vfs = 1.0
         s03sr  = 0.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_V ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Very Firm Soil'
      endif

C     Campbell and Bozorgnia (2003), Vertical, Soft Rock, SS and Reverse
C     Model Number = 078
      if ( jcalc .eq. 78 ) then
c     Set site condition = Soft Rock (Vs~421)
         s03vfs = 0.0
         s03sr  = 1.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_V ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Soft Rock'
      endif

C     Campbell and Bozorgnia (2003), Vertical, Firm Rock, SS and Reverse
C     Model Number = 079
      if ( jcalc .eq. 79 ) then
c     Set site condition = Firm Rock (Vs~830m/s)
         s03vfs = 0.0
         s03sr  = 0.0
         s03fr  = 1.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_V ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Firm Rock'
      endif

C     Campbell and Bozorgnia (2003), Vertical, Generic Rock, SS and Reverse
C     Model Number = 080
      if ( jcalc .eq. 80 ) then
c     Set site condition = Generic Rock (Vs~620m/s)
         s03vfs = 0.0
         s03sr  = 0.5
         s03fr  = 0.5
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_V ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,GenericRock'
      endif

C     Campbell and Bozorgnia (2003), Vertical, Generic Soil, SS and Reverse
C     Model Number = 081
      if ( jcalc .eq. 81 ) then
c     Set site condition = Generic Soil (Vs~310m/s)
         s03vfs = 0.25
         s03sr  = 0.0
         s03fr  = 0.0
C     Set mechanism term for Reverse/Thrust events.
         if (ftype .eq. 1) then
             frv = 0.5
             fth = 0.5
         else
             frv = 0.0
             fth = 0.0
         endif

         call S02_Camp03_V ( mag, seismoDist, jbDist, lnY, sigma,
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag,
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,GenericSoil'
      endif


c ******* Abrahamson and Silva Models *********
C     Abrahamson&Silva 1997 (Rock) - horizontal
C     Model Number = 001
      if ( jcalc .eq. 1 ) then
         soil = 0.
         call S02_AS_97_H ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock'
       endif

C     Abrahamson&Silva 1997 (Rock) - vertical
C     Model Number = 002
      if ( jcalc .eq. 2 ) then
         soil = 0.
         call S02_AS_97_V ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,rock'
       endif

C     Abrahamson&Silva 1997 (Soil) - horizontal
C     Model Number = 003
      if ( jcalc .eq. 3 ) then
         soil = 1.
         call S02_AS_97_H ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,soil'
       endif

C     Abrahamson&Silva 1995 (Soil) - vertical
C     Model Number = 004
      if ( jcalc .eq. 4 ) then
         soil = 1.
         call S02_AS_97_V ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,soil'
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors.
C     Model Number = 005
      if ( jcalc .eq. 5 ) then
         soil = 0.
         call S02_AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock,NF'
       endif

C     Abrahamson&Silva 1997 (Rock) - vertical
C     Model Number = 006
      if ( jcalc .eq. 6 ) then
         soil = 0.
         call S02_AS_97_V_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,rock,NF'
       endif

C     Abrahamson&Silva 1997 (Soil) - horizontal with Normal faulting factors.
C     Model Number = 007
      if ( jcalc .eq. 7 ) then
         soil = 1.
         call S02_AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,soil,NF'
       endif

C     Abrahamson&Silva 1997 (Soil) - vertical with Normal faulting factors
C     Model Number = 008
      if ( jcalc .eq. 8 ) then
         soil = 1.
         call S02_AS_97_V_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,soil,NF'
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors scaled
C           by 1.0/1.67 factor.
C     Model Number = 2005
      if ( jcalc .eq. 2005 ) then
         soil = 0.
         call S02_AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock,NF(1/1.67)'
C      Scale the computed ground motion by 1.0/1.67
         lnY = lnY + alog(1.0/1.67)
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors scaled
C           by 1.67 factor.
C     Model Number = 3005
      if ( jcalc .eq. 3005 ) then
         soil = 0.
         call S02_AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock,NF(1.67)'
C      Scale the computed ground motion by 1.67
         lnY = lnY + alog(1.67)
       endif


c ********* Boore, Joyner and Fumal Models ******************
c     BJF94, Horizontal, Class A
C     Model Number = 010
      if ( jcalc .eq. 10 ) then
         GB = 0
         GC = 0
         call S02_bjf94 ( mag, jbDist, ftype, lnY, sigma, GB, GC, specT,
     1            attenName1, period2,iflag )
c     BJF94, Horizontal, Class B
C     Model Number = 011
      elseif ( jcalc .eq. 11 ) then
         GB = 1
         GC = 0
         call S02_bjf94 ( mag, jbDist, ftype, lnY, sigma, GB, GC, specT,
     1            attenName1, period2,iflag )
c     BJF94, Horizontal, Class c
C     Model Number = 012
      elseif ( jcalc .eq. 12 ) then
         GB = 0
         GC = 1
         call S02_bjf94 ( mag, jbDist, ftype, lnY, sigma, GB, GC, specT,
     1            attenName1, period2,iflag )
      endif

c     BJF97, Horizontal, Vs top 30 meters
C     Model Number = 013
      if ( jcalc .eq. 13 ) then
         call S02_bjf97 ( mag, jbDist, ftype, lnY, sigma, specT,
     1            attenName1, period2, vs,iflag )
      endif

c     BJF97, Horizontal, Vs top 30 meters scaled by factor 1.0/1.67
C     Model Number = 2013
      if ( jcalc .eq. 2013 ) then
         call S02_bjf97 ( mag, jbDist, ftype, lnY, sigma, specT,
     1            attenName1, period2, vs,iflag )
c     Scale ground motion by factor 1.0/1.67
         lnY = lnY + alog(1.0/1.67)
         attenName1='Boore, Joyner, Fumal (1997) (1.0/1.67)'
      endif

c     BJF97, Horizontal, Vs top 30 meters scaled by factor 1.67
C     Model Number = 3013
      if ( jcalc .eq. 3013 ) then
         call S02_bjf97 ( mag, jbDist, ftype, lnY, sigma, specT,
     1            attenName1, period2, vs,iflag )
c     Scale ground motion by factor 1.67
         lnY = lnY + alog(1.67)
         attenName1='Boore, Joyner, Fumal (1997) (1.67)'
      endif

c ******** Campbell Models ******
c     Campbell (1990), Horizontal, Rock
C     Model Number = 020
      if ( jcalc .eq. 20 ) then
         baseDepth = 2.0
         call S02_Camp90 ( mag, seismoDist, ftype, lnY, sigma, baseDepth,
     1            specT, attenName1, period2,iflag )
      endif

c     Campbell (1990) - vertical, Rock
C     Model Number = 021
      if ( jcalc .eq. 21 ) then
         baseDepth = 2.0
         call S02_Camp90v ( mag, seismoDist, ftype, lnY, sigma, baseDepth,
     1            specT, attenName1, period2,iflag )
      endif

c     Campbell (1990/1994), Horizontal, Rock
C     Model Number = 022
      if ( jcalc .eq. 22 ) then
         baseDepth = 2.0
         soilflag = 0
         call S02_Camp90_94 ( mag, seismoDist, ftype, lnY, sigma, baseDepth,
     1             specT, attenName1, period2,iflag )
      endif

c     Campbell (1993-1994), Horizontal Soil
C     Model Number = 023
      if ( jcalc .eq. 23 ) then
         baseDepth = 2.0
         soilFlag = 0
         softRock = 0
         hardRock = 0
         call S02_Campbell_94 ( mag, seismoDist, ftype, lnY, sigma, specT,
     1                soilFlag, softRock, hardRock, baseDepth,
     1                attenName1, period2,iflag )
c     Campbell (1993-1994), Horizontal Soft Rock
C     Model Number = 024
      elseif ( jcalc .eq. 24 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 1
         hardRock = 0
         call S02_Campbell_94 ( mag, seismoDist, ftype, lnY, sigma, specT,
     1                soilFlag, softRock, hardRock, baseDepth,
     1                attenName1, period2,iflag )
c     Campbell (1993-1994), Horizontal Hard Rock
C     Model Number = 025
      elseif ( jcalc .eq. 25 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 0
         hardRock = 1
         call S02_Campbell_94 ( mag, seismoDist, ftype, lnY, sigma, specT,
     1                soilFlag, softRock, hardRock, baseDepth,
     1                attenName1, period2,iflag )
      endif

c     Campbell (1997), Horizontal, Soil
C     Model Number = 026
      if ( jcalc .eq. 26 ) then
         baseDepth = 4.0
         soilFlag = 0
         softRock = 0
         hardRock = 0
         call S02_Camp97_H ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
c     Campbell (1997), Horizontal, Soft Rock
C     Model Number = 027
      elseif ( jcalc .eq. 27 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 1
         hardRock = 0
         call S02_Camp97_H ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
c     Campbell (1997), Horizontal, Hard Rock
C     Model Number = 028
      elseif ( jcalc .eq. 28 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 0
         hardRock = 1
         call S02_Camp97_H ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
      endif

c     Campbell (1997)  vertical, Soil
C     Model Number = 029
      if ( jcalc .eq. 29 ) then
         baseDepth = 4.0
         soilFlag = 0
         softRock = 0
         hardRock = 0
         call S02_Camp97_Z ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
c     Campbell (1997)  vertical, Soft Rock
C     Model Number = 030
      elseif ( jcalc .eq. 30 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 1
         hardRock = 0
         call S02_Camp97_Z ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
c     Campbell (1997)  vertical, Hard Rock
C     Model Number = 031
      elseif ( jcalc .eq. 31 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 0
         hardRock = 1
         call S02_Camp97_Z ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
      endif

c ******** Idriss Models *******
c     Idriss (1991), Horizontal, Rock
C     Model Number = 040
      if ( jcalc .eq. 40 ) then
         call S02_Idriss91_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
      endif

C     Idriss (1991), Horizontal, Soft-soil, PGA
C     Model Number = 041
c     soft-soil
      if ( jcalc .eq. 41 .and. specT .eq. 0.0 ) then
         call S02_Idriss91_soft (mag, rupDist, ftype, lnY, sigma,
     1            attenName1, period2 )
         iflag = 0
      elseif (jcalc .eq. 41 .and. specT .ne. 0.0) then
        write (*,*) 'Idriss (1991), Horizontal, Soft Soil'
        write (*,*) 'not defined for spectral acclereation!!!'
        write (*,*) 'Check input file.'
        stop 99
      endif

c     Idriss 1997 Horizontal, soft-soil, PGA
C     Model Number = 042
      if ( jcalc .eq. 42 .and. specT .eq. 0.0 ) then
         call S02_Idriss97_soft (mag, rupDist, ftype, lnY, sigma,
     1            attenName1, period2 )
         iflag = 0
      elseif (jcalc .eq. 42 .and. specT .ne. 0.0) then
         write (*,*) 'Idriss (1997), Horizontal, Soft Soil'
         write (*,*) 'not defined for spectral acclereation!!!'
         write (*,*) 'Check input file.'
         stop 99
      endif

c     Idriss (1991:1995), Horizontal, Rock
C     Model Number = 043
      if ( jcalc .eq. 43 ) then
         call S02_Idriss91_95_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
      endif

c     Idriss (1991:1995), Horizontal, Rock scaled by factor 1.0/1.67
C     Model Number = 2043
      if ( jcalc .eq. 2043 ) then
         call S02_Idriss91_95_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
c     Scaled ground motion by factor 1.0/1.67
         lnY = lnY + alog (1.0/1.67)
         attenName1 = 'Idriss (1991;1995), Rock (1.0/1.67)'
      endif

c     Idriss (1991:1995), Horizontal, Rock scaled by factor 1.67
C     Model Number = 3043
      if ( jcalc .eq. 3043 ) then
         call S02_Idriss91_95_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
c     Scaled ground motion by factor 1.67
         lnY = lnY + alog (1.67)
         attenName1 = 'Idriss (1991;1995), Rock (1.67)'
      endif

c  ******* Sadigh/Geomatrix Models *******
c     Geomatrix 93 (rock) vertical
C     Model Number = 050
      if ( jcalc .eq. 50 ) then
         call S02_Geomatrix93_V_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
      endif

c     Sadigh et al. 97 (rock) Horizontal
C     Model Number = 051
      if ( jcalc .eq. 51 ) then
         call S02_Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
      endif

c     Sadigh et al. 97 (rock) Horizontal scaled by factor 1.0/1.67
C     Model Number = 2051
      if ( jcalc .eq. 2051 ) then
         call S02_Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
c     Scale ground motion by factor 1.0/1.67
         lnY = lnY + alog(1.0/1.67)
         attenName1 = 'Sadigh et al. (1997), Horizontal, rock (1.0/1.67)'
      endif

c     Sadigh et al. 97 (rock) Horizontal scaled by factor 1.67
C     Model Number = 3051
      if ( jcalc .eq. 3051 ) then
         call S02_Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
c     Scale ground motion by factor 1.67
         lnY = lnY + alog(1.67)
         attenName1 = 'Sadigh et al. (1997), Horizontal, rock (1.67)'
      endif

c     Sadigh et al. 97 (soil) horizontal
C     Model Number = 052
      if ( jcalc .eq. 52 ) then
         call S02_Sadigh97_H_soil ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
      endif

c     Sadigh et al. 97 (rock) Horizontal - Sigma = 0.0
C     Model Number = 053
      if ( jcalc .eq. 53 ) then
         call S02_Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
c     Now set sigma = 0.0
         sigma = 1.0e-10
         attenName1 = 'Sadigh et al. (1997),Hor.,rock,Sigma=0.0'
      endif

c ******** Spudich et al. (1997) Models *******
C     Spudich et al. (1997), Horizontal, Rock, Extensional Regimes
C     Model Number = 060
      if (jcalc .eq. 60) then
         call S02_Spudich96 ( mag, JBdist, lnY, sigma, 0, specT,
     1                   attenName1, period2,iflag )
         attenname1 = 'Spudich et al. (1997), Horizontal, Rock'
      endif

C     Spudich et al. (1997), Horizontal, Soil, Extensional Regimes
C     Model Number = 061
      if (jcalc .eq. 61) then
         call S02_Spudich96 ( mag, JBdist, lnY, sigma, 1, specT,
     1                   attenName1, period2,iflag )
         attenname1 = 'Spudich et al. (1997), Horizontal, Soil'
      endif

c ******** Youngs Models *******
c     Youngs et al (1993) Horizontal, subduction, Rock
C     Model Number = 200
      if ( jcalc .eq. 200 ) then
         call S02_youngs93 ( mag, rupDist, lnY, sigma, attenName1,
     1       period2, specT, ftype,iflag )
         attenName1 ='Youngs (1993) Rock, Subduction'
      endif

c     Youngs et al (1997) Horizontal, subduction, rock
C     Model Number = 201
      if ( jcalc .eq. 201 ) then
         call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1       period2, specT, ftype, depth,iflag )
         attenName1 ='Youngs (1997) Rock, Subduction'
      endif

c     Youngs et al (1997) Horizontal, subduction, soil
C     Model Number = 202
      if ( jcalc .eq. 202 ) then
         call S02_youngs97_soil ( mag, rupDist, lnY, sigma, attenName1,
     1       period2, specT, ftype, depth,iflag )
         attenName1 ='Youngs (1997) Soil, Subduction'
      endif

c ***** Synchronous Rupture Ground Motion Models for HBIP *****
c       Model consists of the SRSS from the Subduction GM
c       and the Crustal GM.

c     This first suite of synchronous ruptures is for seismic source
c          Model A (Geomatrix). Fault parameters for crustal event are:
c          M =7.4, D=1.0 km, Ftype=1 (reverse).

C     Now call the synchronous rupture cases for rock ground motions
c     (i.e., no HBIP amp factors applied).
c     Set the crustal earthquake parameters
      if (jcalc .eq. 7201001 .or. jcalc .eq. 7201043 .or. jcalc .eq.
     1      7201051) then
         rupdistc = 1.0
         ftypec = 1
         magc = 7.4
C     Now call the specific Subduction and Crustal Attenuation Model.

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Abrahamson&Silva (1997) Horizontal Crustal
C     Model Number = 7201001
         if (jcalc .eq. 7201001) then
            call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            soil = 0.
            hwflag = 1
            call S02_AS_97_H ( magc, rupDistc, ftypec, soil, hwflag, specT,
     1                     period2, lnYc, sigmac,iflag )
            attenname1 = 'Syn. Rupture: Youngs+AS (M7.4)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Idriss (1991:1995) Horizontal Crustal
C     Model Number = 7201043
         elseif ( jcalc .eq. 7201043 ) then
            call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call S02_Idriss91_95_rock (magc, rupDistc, ftypec, lnYc, sigmac, specT,
     1            attenName1, period2,iflag )
            attenname1 = 'Synchronous Rupture: Youngs+Idriss (M7.4)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Sadigh et al. (1997) Horizontal Crustal
C     Model Number = 7201051
         elseif ( jcalc .eq. 7201051 ) then
            call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call S02_Geomatrix93_H_rock ( magc, rupDistc, ftypec, lnYc, sigmac,
     1             specT, attenName1, period2,iflag )
            attenname1 = 'Syn. Rupture: Youngs+Sadigh (M7.4)-Rock'
         endif

c     First Combine the Sigma values.
         sigma=sqrt(((sigma**2)*(exp(lnY-6.89))**4 +
     1          (sigmac**2)*(exp(lnYc-6.89))**4)
     2          / ( ((exp(lnY-6.89))**2 + (exp(lnYc-6.89))**2))**2 )

c     Now combine the Subduction GM and Crustal GM values.
         lnY = alog ( sqrt( (exp(lnY))**2 + (exp(lnYc))**2 ) )
      endif

c     The second suite of synchronous ruptures is for seismic source
c          Model B (Carver). Fault parameters for crustal event are:
c          M =7.7, D=1.0 km, Ftype=1 (reverse)

C     Now call the synchronous rupture cases for rock ground motions
c     (i.e., no HBIP amp factors applied).
c     Set the crustal earthquake parameters
      if (jcalc .eq. 8201001 .or. jcalc .eq. 8201043 .or. jcalc .eq.
     1      8201051) then
         rupdistc = 1.0
         ftypec = 1
         magc = 7.7
C     Now call the specific Subduction and Crustal Attenuation Model.

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Abrahamson&Silva (1997) Horizontal Crustal
C     Model Number = 8201001
         if (jcalc .eq. 8201001) then
            call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            soil = 0.
            hwflag = 1
            call S02_AS_97_H ( magc, rupDistc, ftypec, soil, hwflag, specT,
     1                     period2, lnYc, sigmac,iflag )
            attenname1 = 'Syn. Rupture: Youngs+AS (M7.7)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Idriss (1991:1995) Horizontal Crustal
C     Model Number = 8201043
         elseif ( jcalc .eq. 8201043 ) then
            call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call S02_Idriss91_95_rock (magc, rupDistc, ftypec, lnYc, sigmac, specT,
     1            attenName1, period2,iflag )
            attenname1 = 'Synchronous Rupture: Youngs+Idriss (M7.7)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Sadigh et al. (1997) Horizontal Crustal
C     Model Number = 8201051
         elseif ( jcalc .eq. 8201051 ) then
            call S02_youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call S02_Geomatrix93_H_rock ( magc, rupDistc, ftypec, lnYc, sigmac,
     1             specT, attenName1, period2,iflag )
            attenname1 = 'Syn. Rupture: Youngs+Sadigh (M7.7)-Rock'
         endif

c     First Combine the Sigma values.
         sigma=sqrt(((sigma**2)*(exp(lnY-6.89))**4 +
     1          (sigmac**2)*(exp(lnYc-6.89))**4)
     2          / ( ((exp(lnY-6.89))**2 + (exp(lnYc-6.89))**2))**2 )

c     Now combine the Subduction GM and Crustal GM values.
         lnY = alog ( sqrt( (exp(lnY))**2 + (exp(lnYc))**2 ) )
      endif

c ***** Atkinson and Boore Subduction Models *****
c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 210
      if (jcalc .eq. 210) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-B'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 211
      if (jcalc .eq. 211) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-C'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 212
      if (jcalc .eq. 212) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-D'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 213
      if (jcalc .eq. 213) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-E'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 220
      if (jcalc .eq. 220) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-B'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 221
      if (jcalc .eq. 221) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-C'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 222
      if (jcalc .eq. 222) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-D'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 223
      if (jcalc .eq. 223) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-E'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 230
      if (jcalc .eq. 230) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-B'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 231
      if (jcalc .eq. 231) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-C'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 232
      if (jcalc .eq. 232) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-D'
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 233
      if (jcalc .eq. 233) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-E'
      endif

C Atkinson and Boore (2008) Subduction Erratum Corrected for Interface Events ****
c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-B, Subduction
C     Model Number = 310
      if (jcalc .eq. 310) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-B'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the
c         interpolation for periods between these end points.
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call S02_AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call S02_AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call S02_AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call S02_AB03 ( mag, rupdist, lnY10, sigma10, specT10,
     1              attenName0, period02,iflag10, ftype, depth, Sc, Sd, Se)
               period2 = specT

C     SpecT falls between 0.1 and 0.2sec
               if (specT .gt. 0.1 .and. specT .le. 0.2) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Interpolate to given SpecT value with new 0.2sec ground motion value.
                 lnY = lnY01 + (lnY01-lnY02p)*(alog(specT)-alog(0.1))/(alog(0.1)-alog(0.2))
C     SpecT falls between 0.2 and 0.4sec
               elseif (specT .gt. 0.2 .and. specT .le. 0.4) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 0.2 and 0.4sec ground motion values.
                 lnY = lnY02p + (lnY02p-lnY04p)*(alog(specT)-alog(0.2))/(alog(0.2)-alog(0.4))
C     SpecT falls between 0.4 and 1.0sec
               elseif (specT .gt. 0.4 .and. specT .lt. 1.0) then
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 1.0sec ground motion value.
                 lnY = lnY04p + (lnY04p-lnY10)*(alog(specT)-alog(0.4))/(alog(0.4)-alog(1.0))
               endif
            endif
	   endif
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-C, Subduction
C     Model Number = 311
      if (jcalc .eq. 311) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-C'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the
c         interpolation for periods between these end points.
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call S02_AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call S02_AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call S02_AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call S02_AB03 ( mag, rupdist, lnY10, sigma10, specT10,
     1              attenName0, period02,iflag10, ftype, depth, Sc, Sd, Se)
               period2 = specT

C     SpecT falls between 0.1 and 0.2sec
               if (specT .gt. 0.1 .and. specT .le. 0.2) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Interpolate to given SpecT value with new 0.2sec ground motion value.
                 lnY = lnY01 + (lnY01-lnY02p)*(alog(specT)-alog(0.1))/(alog(0.1)-alog(0.2))
C     SpecT falls between 0.2 and 0.4sec
               elseif (specT .gt. 0.2 .and. specT .le. 0.4) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 0.2 and 0.4sec ground motion values.
                 lnY = lnY02p + (lnY02p-lnY04p)*(alog(specT)-alog(0.2))/(alog(0.2)-alog(0.4))
C     SpecT falls between 0.4 and 1.0sec
               elseif (specT .gt. 0.4 .and. specT .lt. 1.0) then
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 1.0sec ground motion value.
                 lnY = lnY04p + (lnY04p-lnY10)*(alog(specT)-alog(0.4))/(alog(0.4)-alog(1.0))
               endif
            endif
	   endif
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-D, Subduction
C     Model Number = 312
      if (jcalc .eq. 312) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-D'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the
c         interpolation for periods between these end points.
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call S02_AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call S02_AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call S02_AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call S02_AB03 ( mag, rupdist, lnY10, sigma10, specT10,
     1              attenName0, period02,iflag10, ftype, depth, Sc, Sd, Se)
               period2 = specT
C     SpecT falls between 0.1 and 0.2sec
               if (specT .gt. 0.1 .and. specT .le. 0.2) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Interpolate to given SpecT value with new 0.2sec ground motion value.
                 lnY = lnY01 + (lnY01-lnY02p)*(alog(specT)-alog(0.1))/(alog(0.1)-alog(0.2))
C     SpecT falls between 0.2 and 0.4sec
               elseif (specT .gt. 0.2 .and. specT .le. 0.4) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 0.2 and 0.4sec ground motion values.
                 lnY = lnY02p + (lnY02p-lnY04p)*(alog(specT)-alog(0.2))/(alog(0.2)-alog(0.4))
C     SpecT falls between 0.4 and 1.0sec
               elseif (specT .gt. 0.4 .and. specT .lt. 1.0) then
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 1.0sec ground motion value.
                 lnY = lnY04p + (lnY04p-lnY10)*(alog(specT)-alog(0.4))/(alog(0.4)-alog(1.0))
               endif
            endif
	   endif
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-E, Subduction
C     Model Number = 313
      if (jcalc .eq. 313) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call S02_AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-E'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the
c         interpolation for periods between these end points.
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call S02_AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call S02_AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call S02_AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call S02_AB03 ( mag, rupdist, lnY10, sigma10, specT10,
     1              attenName0, period02,iflag10, ftype, depth, Sc, Sd, Se)
               period2 = specT
C     SpecT falls between 0.1 and 0.2sec
               if (specT .gt. 0.1 .and. specT .le. 0.2) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Interpolate to given SpecT value with new 0.2sec ground motion value.
                 lnY = lnY01 + (lnY01-lnY02p)*(alog(specT)-alog(0.1))/(alog(0.1)-alog(0.2))
C     SpecT falls between 0.2 and 0.4sec
               elseif (specT .gt. 0.2 .and. specT .le. 0.4) then
c     Correct the 0.2sec ground motion value.
                 lnY02p = 0.667*lnY04 + 0.333*lnY02
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 0.2 and 0.4sec ground motion values.
                 lnY = lnY02p + (lnY02p-lnY04p)*(alog(specT)-alog(0.2))/(alog(0.2)-alog(0.4))
C     SpecT falls between 0.4 and 1.0sec
               elseif (specT .gt. 0.4 .and. specT .lt. 1.0) then
c     Correct the 0.4sec ground motion value.
                 lnY04p = 0.333*lnY04 + 0.667*lnY02
c     Interpolate to given SpecT value with new 1.0sec ground motion value.
                 lnY = lnY04p + (lnY04p-lnY10)*(alog(specT)-alog(0.4))/(alog(0.4)-alog(1.0))
               endif
            endif
	   endif
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-B, Subduction
C     Model Number = 320
      if (jcalc .eq. 320) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-B'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-C, Subduction
C     Model Number = 321
      if (jcalc .eq. 321) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-C'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-D, Subduction
C     Model Number = 322
      if (jcalc .eq. 322) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-D'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-E, Subduction
C     Model Number = 323
      if (jcalc .eq. 323) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call S02_AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-E'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-B, Subduction
C     Model Number = 330
      if (jcalc .eq. 330) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-B'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-C, Subduction
C     Model Number = 331
      if (jcalc .eq. 331) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-C'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-D, Subduction
C     Model Number = 332
      if (jcalc .eq. 332) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-D'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-E, Subduction
C     Model Number = 333
      if (jcalc .eq. 333) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call S02_AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-E'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

      endif

c ***** Gregor et al. Cascadia Subduction Model *****
c     Gregor et al. (2002) - Horizontal, Rock, Cascadia Subduction
C     Model Number = 240
C     Note: This model is only valid for the magnitude range of the
C           the data used in the regression: 8.0 - 9.0.  The program
C           will check to make sure that the magnitude range being
C           requested falls within this range.
      if (jcalc .eq. 240) then
C        Check the magnitude range.
c         if (mag .lt. 8.0 .or. mag .gt. 9.0) then
c            write (*,*) 'User has requested the Gregor et al. (2002)'
c            write (*,*) 'Cascadia Subduction attenuation relationship,'
c            write (*,*) 'however, the magnitude of = ', mag
c            write (*,*) 'falls outside of the applicable magnitude'
c            write (*,*) 'range for this attenuation relationship.'
c            write (*,*) 'Check your input fault parameter file.'
c            stop 99
c         endif
         call S02_Gregor02CasR ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag)
         attenname1 = 'Gregor et al.(2002),Rock,Cascadia Subduction'
      endif

c     Gregor et al. (2002) - Horizontal, Soil, Cascadia Subduction
C     Model Number = 241
C     Note: This model is only valid for the magnitude range of the
C           the data used in the regression: 8.0 - 9.0.  The program
C           will check to make sure that the magnitude range being
C           requested falls within this range.
      if (jcalc .eq. 241) then
C        Check the magnitude range.
c         if (mag .lt. 8.0 .or. mag .gt. 9.0) then
c            write (*,*) 'User has requested the Gregor et al. (2002)'
c            write (*,*) 'Cascadia Subduction attenuation relationship,'
c            write (*,*) 'however, the magnitude of = ', mag
c            write (*,*) 'falls outside of the applicable magnitude'
c            write (*,*) 'range for this attenuation relationship.'
c            write (*,*) 'Check your input fault parameter file.'
c            stop 99
c         endif
         call S02_Gregor02CasS ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag)
         attenname1 = 'Gregor et al.(2002),Soil,Cascadia Subduction'
      endif

c ***** Gregor et al. Cascadia Subduction Model *****
c     Gregor et al. (2006) - Horizontal, Cascadia Subduction
C     Model Number = 242

      if (jcalc .eq. 242) then
         call S02_Gregor06Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, vs, period2, iflag)
         attenname1 = 'Gregor et al.(2006), Cascadia Subduction'
      endif

c ***** Zhao et al. (2006) Crustal/Subduction Models *****
c     Zhao et al. (2006) - Horizontal, Hard Rock, Subduction
C     Model Number = 250
      if (jcalc .eq. 250) then
c        Set Site Class = Hard Rock and Source Type
         sclass = 0
         if (ftype .eq. 1) then
            sourceclass = 2
         elseif (ftype .eq. 0) then
            sourceclass = 1
         endif

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Sub., Hard Rock'
      endif

c     Zhao et al. (2006) - Horizontal, Rock, Subduction
C     Model Number = 251
      if (jcalc .eq. 251) then

c        Set Site Class = Rock and Source Type
         sclass = 1.0
         if (ftype .eq. 1) then
            sourceclass = 2.0
         elseif (ftype .eq. 0) then
            sourceclass = 1.0
         endif
         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Sub., Rock (SC I)'
      endif

c     Zhao et al. (2006) - Horizontal, Hard Soil, Subduction
C     Model Number = 252
      if (jcalc .eq. 252) then
c        Set Site Class = Hard Soil and Source Type
         sclass = 2
         if (ftype .eq. 1) then
            sourceclass = 2
         elseif (ftype .eq. 0) then
            sourceclass = 1
         endif

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Sub., Hard Soil (SC II)'
      endif

c     Zhao et al. (2006) - Horizontal, Medium Soil, Subduction
C     Model Number = 253
      if (jcalc .eq. 253) then
c        Set Site Class = Medium Soil and Source Type
         sclass = 3
         if (ftype .eq. 1) then
            sourceclass = 2
         elseif (ftype .eq. 0) then
            sourceclass = 1
         endif

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Sub., Medium Soil (SC III)'
      endif

c     Zhao et al. (2006) - Horizontal, Soft Soil, Subduction
C     Model Number = 254
      if (jcalc .eq. 254) then
c        Set Site Class = Soft Soil and Source Type
         sclass = 4
         if (ftype .eq. 1) then
            sourceclass = 2
         elseif (ftype .eq. 0) then
            sourceclass = 1
         endif

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Sub., Soft Soil (SC IV)'
      endif

c     Zhao et al. (2006) - Horizontal, Hard Rock, Crustal
C     Model Number = 255
      if (jcalc .eq. 255) then
c        Set Site Class = Hard Rock and Source Type
         sclass = 0
         sourceclass = 0

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Crust, Hard Rock'
      endif

c     Zhao et al. (2006) - Horizontal, Rock, Crustal
C     Model Number = 256
      if (jcalc .eq. 256) then
c        Set Site Class = Rock and Source Type
         sclass = 1.0
         sourceclass = 0

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Crust, Rock (SC I)'
      endif

c     Zhao-Lu (2011) - Horizontal, Rock, Crustal
C     Modified (2006) model with magnitude capped at 7.1.
C     Model Number = 2256
      if (jcalc .eq. 2256) then
c        Set Site Class = Rock and Source Type
         sclass = 1.0
         sourceclass = 0
         if (mag .gt. 7.1) then
            m = 7.1
         else
            m = mag
         endif
         call S02_Zhaoetal2006 ( m, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
         attenname1 = 'ZhaoLu(2011)-Crust, Rock (SC I)'
      endif



c     Zhao et al. (2006) - Horizontal, Hard Soil, Crustal
C     Model Number = 257
      if (jcalc .eq. 257) then
c        Set Site Class = Hard Soil and Source Type
         sclass = 2
         sourceclass = 0

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Crust, Hard Soil (SC II)'
      endif

c     Zhao et al. (2006) - Horizontal, Medium Soil, Crustal
C     Model Number = 258
      if (jcalc .eq. 258) then
c        Set Site Class = Medium Soil and Source Type
         sclass = 3
         sourceclass = 0

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Crust, Medium Soil (SC III)'
      endif

c     Zhao et al. (2006) - Horizontal, Soft Soil, Crustal
C     Model Number = 259
      if (jcalc .eq. 259) then
c        Set Site Class = Soft Soil and Source Type
         sclass = 4
         sourceclass = 0

         call S02_Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )

         attenname1 = 'Zhao etal(2006)-Crust, Soft Soil (SC IV)'
      endif

c ***** Kanno Subduction Models *****
c     Kanno (2006) - Horizontal, Subduction
C     Model Number = 260
      if ( jcalc .eq. 260 ) then

         call S02_kanno2006 ( mag, rupdist, specT,
     1                    period2, lnY, sigma, iflag, vs, depth )

         attenname1 = 'Kanno Subduction (2006)'
      endif

c ***** Garcia et al. (2005) Subduction Model-Inslab *****
c     Garcia et al. (2005) - Horizontal, Subduction-Inslab
C     Model Number = 270
      if ( jcalc .eq. 270 ) then
         if (ftype .ne. 1) then
            write (*,*)
            write (*,*) 'Garcia et al. (2005) Subduction atttenuation model'
            write (*,*) 'is defined only for inslab events'
            write (*,*) 'Mechanism Ftype = 1'
            write (*,*) 'Please check the input file.'
         else
            call S02_GarciaH05 ( mag, rupdist, specT,
     1                    period2, lnY, sigma, iflag, depth )
            attenname1 = 'Garcia et al. (2005), Hor-Inslab, Rock'
         endif
       endif

C     Model Number = 271
c     Garcia et al. (2005) - Vertical, Subduction-Inslab
      if ( jcalc .eq. 271 ) then
         if (ftype .ne. 1) then
            write (*,*)
            write (*,*) 'Garcia et al. (2005) Subduction atttenuation model'
            write (*,*) 'is defined only for inslab events'
            write (*,*) 'Mechanism Ftype = 1'
            write (*,*) 'Please check the input file.'
         else
            call S02_GarciaV05 ( mag, rupdist, specT,
     1                       period2, lnY, sigma, iflag, depth )
            attenname1 = 'Garcia et al. (2005), Ver-Inslab, Rock'
         endif
       endif

c ***** Lin and Lee (2008) Subduction Model *****
c      Lin and Lee (2008) - Horizontal, Subduction, Rock
C     Model Number = 280
      if ( jcalc .eq. 280 ) then
         call S02_LinLee08rock ( mag, rupdist, specT,
     1                    period2, lnY, sigma, iflag, depth, ftype )
         attenname1 = 'Lin and Lee (2008), Subduction, Rock'
      endif

c      Lin and Lee (2008) - Horizontal, Subduction, Soil
C     Model Number = 281
      if ( jcalc .eq. 281 ) then
         call S02_LinLee08soil ( mag, rupdist, specT,
     1                    period2, lnY, sigma, iflag, depth, ftype )
         attenname1 = 'Lin and Lee (2008), Subduction, Soil'
      endif

C  **** BCSubduction Model *******
C     Base Case Model, DeltaC1 = 0.0
C     Model Number = 350
      if ( jcalc .eq. 350 ) then
         deltaC1 = 0.0
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, DeltaC1=0'
      endif

C     Base Case Model, DeltaC1 = -0.5
C     Model Number = 351
      if ( jcalc .eq. 351 ) then
         deltaC1 = -0.5
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, DeltaC1=-0.5'
      endif

C     Base Case Model, DeltaC1 = 0.5
C     Model Number = 352
      if ( jcalc .eq. 352 ) then
         deltaC1 = 0.5
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, DeltaC1=0.5'
      endif

C     Base Case Model, DeltaC1 = 0.0, Single Station Sigma=0.60
C     Model Number = 353
      if ( jcalc .eq. 353 ) then
         deltaC1 = 0.0
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=0'
      endif

C     Base Case Model, DeltaC1 = -0.5, Single Station Sigma=0.60
C     Model Number = 354
      if ( jcalc .eq. 354 ) then
         deltaC1 = -0.5
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=-0.5'
      endif

C     Base Case Model, DeltaC1 = 0.5, Single Station Sigma=0.58
C     Model Number = 355
      if ( jcalc .eq. 355 ) then
         deltaC1 = 0.5
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=0.5'
      endif

C     Base Case Model, DeltaC1 = -0.2
C     Model Number = 356
      if ( jcalc .eq. 356 ) then
         deltaC1 = -0.2
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, DeltaC1=-0.2'
      endif

C     Base Case Model, DeltaC1 = 0.2
C     Model Number = 357
      if ( jcalc .eq. 357 ) then
         deltaC1 = 0.2
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, DeltaC1=0.2'
      endif

C     Base Case Model, DeltaC1 = -0.2, Single Station Sigma=0.60
C     Model Number = 358
      if ( jcalc .eq. 358 ) then
         deltaC1 = -0.2
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=-0.2'
      endif

C     Base Case Model, DeltaC1 = 0.2, Single Station Sigma=0.60
C     Model Number = 359
      if ( jcalc .eq. 359 ) then
         deltaC1 = 0.2
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=0.2'
      endif

C     Base Case Model, Variable DeltaC1 Adjustment - Central Values, Reg. Sigma
C     Model Number = 360
      if ( jcalc .eq. 360 ) then
c        Determine the DeltaC1 value based on recommended adjusted model and spectral period.
C        Period dependent model for interface events and constant for intraslab
         if (ftype .eq. 0) then
         if (specT .le. 0.3) then
            deltaC1 = 0.2
         elseif (specT .gt. 0.3 .and. specT .le. 0.5) then
            deltaC1 = 0.2 + (0.1-0.2)*(alog(specT)-alog(0.3)) / (alog(0.5)-alog(0.3))
         elseif (specT .gt. 0.5 .and. specT .le. 1.0) then
            deltaC1 = 0.1 + (0.0-0.1)*(alog(specT)-alog(0.5)) / (alog(1.0)-alog(0.5))
         elseif (specT .gt.1.0 .and. specT .le. 2.0) then
            deltaC1 = 0.0 + (-0.1-0.0)*(alog(specT)-alog(1.0)) / (alog(2.0)-alog(1.0))
         elseif (specT .gt. 2.0 .and. specT .le. 3.0) then
            deltaC1 = -0.1 + (-0.2+0.1)*(alog(specT)-alog(2.0)) / (alog(3.0)-alog(2.0))
         else
            deltaC1 = -0.2
         endif
         elseif (ftype .eq. 1) then
            deltaC1 = -0.3
         endif
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, Var. Central DeltaC1'
      endif

C     Base Case Model, Variable DeltaC1 Adjustment - Lower Values, Reg. Sigma
C     Model Number = 361
      if ( jcalc .eq. 361 ) then
c        Determine the DeltaC1 value based on recommended adjusted model and spectral period.
C        Period dependent model for interface events and constant for intraslab
         if (ftype .eq. 0) then
         if (specT .le. 0.3) then
            deltaC1 = 0.0
         elseif (specT .gt. 0.3 .and. specT .le. 0.5) then
            deltaC1 = 0.0 + (-0.1-0.0)*(alog(specT)-alog(0.3)) / (alog(0.5)-alog(0.3))
         elseif (specT .gt. 0.5 .and. specT .le. 1.0) then
            deltaC1 = -0.1 + (-0.2+0.1)*(alog(specT)-alog(0.5)) / (alog(1.0)-alog(0.5))
         elseif (specT .gt.1.0 .and. specT .le. 2.0) then
            deltaC1 = -0.2 + (-0.3+0.2)*(alog(specT)-alog(1.0)) / (alog(2.0)-alog(1.0))
         elseif (specT .gt. 2.0 .and. specT .le. 3.0) then
            deltaC1 = -0.3 + (-0.4+0.3)*(alog(specT)-alog(2.0)) / (alog(3.0)-alog(2.0))
         else
            deltaC1 = -0.4
         endif
         elseif (ftype .eq. 1) then
            deltaC1 = -0.5
         endif
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, Var. Lower DeltaC1'
      endif

C     Base Case Model, Variable DeltaC1 Adjustment - Upper Values, Reg. Sigma
C     Model Number = 362
      if ( jcalc .eq. 362 ) then
c        Determine the DeltaC1 value based on recommended adjusted model and spectral period.
C        Period dependent model for interface events and constant for intraslab
         if (ftype .eq. 0) then
         if (specT .le. 0.3) then
            deltaC1 = 0.4
         elseif (specT .gt. 0.3 .and. specT .le. 0.5) then
            deltaC1 = 0.4 + (0.3-0.4)*(alog(specT)-alog(0.3)) / (alog(0.5)-alog(0.3))
         elseif (specT .gt. 0.5 .and. specT .le. 1.0) then
            deltaC1 = 0.3 + (0.2-0.3)*(alog(specT)-alog(0.5)) / (alog(1.0)-alog(0.5))
         elseif (specT .gt.1.0 .and. specT .le. 2.0) then
            deltaC1 = 0.2 + (0.1-0.2)*(alog(specT)-alog(1.0)) / (alog(2.0)-alog(1.0))
         elseif (specT .gt. 2.0 .and. specT .le. 3.0) then
            deltaC1 = 0.1 + (0.0-0.1)*(alog(specT)-alog(2.0)) / (alog(3.0)-alog(2.0))
         else
            deltaC1 = 0.0
         endif
         elseif (ftype .eq. 1) then
            deltaC1 = -0.1
         endif
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         attenname1 = 'BCHydroSub_V3, Var. Upper DeltaC1'
      endif

C     Base Case Model, Variable DeltaC1 Adjustment - Central Values, Single Station Sigma
C     Model Number = 363
      if ( jcalc .eq. 363 ) then
c        Determine the DeltaC1 value based on recommended adjusted model and spectral period.
C        Period dependent model for interface events and constant for intraslab
         if (ftype .eq. 0) then
         if (specT .le. 0.3) then
            deltaC1 = 0.2
         elseif (specT .gt. 0.3 .and. specT .le. 0.5) then
            deltaC1 = 0.2 + (0.1-0.2)*(alog(specT)-alog(0.3)) / (alog(0.5)-alog(0.3))
         elseif (specT .gt. 0.5 .and. specT .le. 1.0) then
            deltaC1 = 0.1 + (0.0-0.1)*(alog(specT)-alog(0.5)) / (alog(1.0)-alog(0.5))
         elseif (specT .gt.1.0 .and. specT .le. 2.0) then
            deltaC1 = 0.0 + (-0.1-0.0)*(alog(specT)-alog(1.0)) / (alog(2.0)-alog(1.0))
         elseif (specT .gt. 2.0 .and. specT .le. 3.0) then
            deltaC1 = -0.1 + (-0.2+0.1)*(alog(specT)-alog(2.0)) / (alog(3.0)-alog(2.0))
         else
            deltaC1 = -0.2
         endif
         elseif (ftype .eq. 1) then
            deltaC1 = -0.3
         endif
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
C     Changed SSS to 0.60 (January 17, 2014)
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3, Var. Central DeltaC1, SSS'
      endif

C     Base Case Model, Variable DeltaC1 Adjustment - Lower Values, Single Station Sigma
C     Model Number = 364
      if ( jcalc .eq. 364 ) then
c        Determine the DeltaC1 value based on recommended adjusted model and spectral period.
C        Period dependent model for interface events and constant for intraslab
         if (ftype .eq. 0) then
         if (specT .le. 0.3) then
            deltaC1 = 0.0
         elseif (specT .gt. 0.3 .and. specT .le. 0.5) then
            deltaC1 = 0.0 + (-0.1-0.0)*(alog(specT)-alog(0.3)) / (alog(0.5)-alog(0.3))
         elseif (specT .gt. 0.5 .and. specT .le. 1.0) then
            deltaC1 = -0.1 + (-0.2+0.1)*(alog(specT)-alog(0.5)) / (alog(1.0)-alog(0.5))
         elseif (specT .gt.1.0 .and. specT .le. 2.0) then
            deltaC1 = -0.2 + (-0.3+0.2)*(alog(specT)-alog(1.0)) / (alog(2.0)-alog(1.0))
         elseif (specT .gt. 2.0 .and. specT .le. 3.0) then
            deltaC1 = -0.3 + (-0.4+0.3)*(alog(specT)-alog(2.0)) / (alog(3.0)-alog(2.0))
         else
            deltaC1 = -0.4
         endif
         elseif (ftype .eq. 1) then
            deltaC1 = -0.5
         endif
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3, Var. Lower DeltaC1, SSS'
      endif

C     Base Case Model, Variable DeltaC1 Adjustment - Upper Values, Single Station Sigma
C     Model Number = 365
      if ( jcalc .eq. 365 ) then
c        Determine the DeltaC1 value based on recommended adjusted model and spectral period.
C        Period dependent model for interface events and constant for intraslab
         if (ftype .eq. 0) then
         if (specT .le. 0.3) then
            deltaC1 = 0.4
         elseif (specT .gt. 0.3 .and. specT .le. 0.5) then
            deltaC1 = 0.4 + (0.3-0.4)*(alog(specT)-alog(0.3)) / (alog(0.5)-alog(0.3))
         elseif (specT .gt. 0.5 .and. specT .le. 1.0) then
            deltaC1 = 0.3 + (0.2-0.3)*(alog(specT)-alog(0.5)) / (alog(1.0)-alog(0.5))
         elseif (specT .gt.1.0 .and. specT .le. 2.0) then
            deltaC1 = 0.2 + (0.1-0.2)*(alog(specT)-alog(1.0)) / (alog(2.0)-alog(1.0))
         elseif (specT .gt. 2.0 .and. specT .le. 3.0) then
            deltaC1 = 0.1 + (0.0-0.1)*(alog(specT)-alog(2.0)) / (alog(3.0)-alog(2.0))
         else
            deltaC1 = 0.0
         endif
         elseif (ftype .eq. 1) then
            deltaC1 = -0.1
         endif
         call S05_BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY,
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 )
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3, Var. Upper DeltaC1, SSS'
      endif

c ***** Atkinson and Macias (2009) Cascadia Model, NEHRP B/C *****
C     Model Number = 370
      if ( jcalc .eq. 370 ) then
         call S02_AM09_Cas ( mag, rupDist, lnY,
     1            sigma, specT, period2, iflag )
         attenname1 = 'Atkinson&Macias, Cascadia, NEHRP B/C'
      endif

C ******  CEUS Models *********
C
c ***** Atkinson and Boore Models *****
c     Atkinson and Boore (1994) - Horizontal, EUS Hard Rock
C     Model Number = 100
      if ( jcalc .eq. 100 ) then
         call S02_AB95 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag )
         attenname1 = 'Atkinson&Boore 1994, EUS, Rock'
      endif

c     Atkinson and Boore (1994), Horizontal, EUS Hard Rock, magnitude Nuttli
C     Model Number = 101
      if ( jcalc .eq. 101 ) then
C         Convert Nuttli magnitude to Moment magnitude.
          if (mag .le. 5.5) then
            mag1 = -0.39 + 0.98*mag
          else
            mag1 = 2.715 - 0.277*mag + 0.127*mag*mag
          endif
          call S02_AB95Mn ( mag1, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag )
          attenname1 = 'Atkinson&Boore 1994, EUS, Rock, Nuttli Mag'
      endif

c     Atkinson and Boore (2006) - Horizontal, EUS Hard Rock
C     Model Number = 102
      if ( jcalc .eq. 102 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
         attenname1 = 'Atkinson&Boore 2006, EUS, Hard Rock'
      endif

c     Atkinson and Boore (2006) - Horizontal, EUS Vs=760m/sec
C     Model Number = 103
      if ( jcalc .eq. 103 ) then
         call S06_AB06vs760 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
         attenname1 = 'Atkinson&Boore 2006, EUS, Vs=760m/sec'
      endif

c     Atkinson (2008, weighted C0) - Horizontal, CEUS-NGA nga Vs=760m/sec
C     Model Number = 104
      if ( jcalc .eq. 104 ) then
c     First call BA08 to get ground motion values which will be adjusted.
         vs = 760.0
         call S07_BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, BA08lnY, sigma, iflag, vs, ftype, pga4nl,
     1                    phi, tau )
C     Now call S06_A08vs760 to adjust BA08 NGA value to CEUS Vs=760 value.
         call S06_A08vs760 ( mag, jbdist, specT, BA08lnY,
     1                    period2, lnY, sigma, iflag )
         attenname1 = 'Atkinson 2008 wt C0, EUS-NGA BA08, Vs=760m/sec'
      endif

c     Atkinson and Boore (2006) with Atikinson (2010) stress drop adjustment - Horizontal, EUS Hard Rock
C     Model Number = 105
      if ( jcalc .eq. 105 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 10**(3.45 - 0.2*mag)
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006/Atkinson2010, EUS, Hard Rock'
         LnY = lnY + sdscale*SF2*alog(10.0)
      endif

c     Atkinson and Boore (2006) with Atikinson (2010) stress drop adjustment - Horizontal, EUS Vs760m/sec
C     Model Number = 106
      if ( jcalc .eq. 106 ) then
         call S06_AB06vs760 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 10**(3.45 - 0.2*mag)
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006/Atkinson2010, EUS, Vs760m/s'
         LnY = lnY + sdscale*SF2*alog(10.0)
      endif

c     Atkinson (2010) - Horizontal, CEUS-NGA NGA Based
C     Model Number = 107
      if ( jcalc .eq. 107 ) then
c     First call BA08 to get ground motion values which will be adjusted.
         call S07_BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, BA08lnY, sigma, iflag, vs, ftype, pga4nl,
     1                    phi, tau )
C     Apply small magnitude adjustment is Mag<=5.75.
         if (mag .le. 5.75 ) then
            a = max(0.0, 3.888 - 0.674*mag)
            b = max(0.0, 2.933 - 0.510*mag)
            factor = a - b*log10(jbdist+10.0)
            factor = factor*alog(10.0)
         endif
C     Now apply Fena adjustment as given in Atkinson (2010).
         call S06_Fena ( jbdist, specT, facFena )
         attenname1 = 'Atkinson 2010, EUS-NGA BA08 based'
         lnY =  BA08lnY + factor + FacFena
      endif

c     Atkinson (2008, average C0) - Horizontal, CEUS-NGA nga Vs=760m/sec
C     Model Number = 108
      if ( jcalc .eq. 108 ) then
c     First call BA08 to get ground motion values which will be adjusted.
         vs = 760.0
         call S07_BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, BA08lnY, sigma, iflag, vs, ftype, pga4nl,
     1                    phi, tau )
C     Now Call A08vs760 to adjust BA08 NGA value to CEUS Vs=760 value.
         call S06_A08vs760C0 ( mag, jbdist, specT, BA08lnY,
     1                    period2, lnY, sigma, iflag )
         attenname1 = 'Atkinson 2008 avg C0, EUS-NGA BA08, Vs=760m/sec'
      endif

c     Atkinson and Boore (2006) with 2x stress drop adjustment (280bars) - Horizontal, EUS Hard Rock
C     Model Number = 130
      if ( jcalc .eq. 130 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 280.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (2x) StressDrop, EUS, Hard Rock'
         LnY = lnY + sdscale*SF2*alog(10.0)
      endif

c     Atkinson and Boore (2006) with 0.5x stress drop adjustment (70bars) - Horizontal, EUS Hard Rock
C     Model Number = 131
      if ( jcalc .eq. 131 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 70.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (0.5x) StressDrop, EUS, Hard Rock'
         LnY = lnY + sdscale*SF2*alog(10.0)
      endif

c ******* Toro et al. Models *******
C Toro et al. (1997) MidCon., Horizontal, Rock
C     Model Number = 110
      if (jcalc .eq. 110) then
         call S02_TAS96 ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Toro et al. (1997), Horizontal, MidCon.'
      endif

C Toro et al. (1997) MidCon., Horizontal, Rock, MLg magnitude
C     Model Number = 111
      if (jcalc .eq. 111) then
         call S02_TAS96MLg ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag)
         attenname1 = 'Toro et al. (1997), Horizontal, MidCon., MLg'
      endif

C Toro et al. (1997) Gulf, Horizontal, Rock
C     Model Number = 112
      if (jcalc .eq. 112) then
         call S02_TAS96Gulf ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Toro et al. (1997), Horizontal, Gulf'
      endif

C Toro et al. (1997) Gulf, Horizontal, Rock, MLg magnitude
C     Model Number = 113
      if (jcalc .eq. 113) then
         call S02_TAS96GulfMLg ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Toro et al. (1997), Horizontal, Gulf, MLg'
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal, Hard Rock
C     Model Number = 120
      if (jcalc .eq. 120) then
         call S06_CHY03 ( mag, rupdist, lnY, sigma, specT,
     1                  period2,iflag )
         attenname1 = 'Campbell (2003), Hor., CEUS-Hybrid, Hard Rock'
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal minus SigmaEps, Hard Rock
C     Model Number = 121
      if (jcalc .eq. 121) then
         call S06_CHY03 ( mag, rupdist, lnY, sigma, specT,
     1                  period2,iflag )
         call S06_CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )
         attenname1 = 'Campbell (2003), Hor - SigmaEps, CEUS-Hybrid, Hard Rock'
         lnY = lnY - sigmaEps
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal plus SigmaEps, Hard Rock
C     Model Number = 122
      if (jcalc .eq. 122) then
         call S06_CHY03 ( mag, rupdist, lnY, sigma, specT,
     1                  period2,iflag )
         call S06_CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )
         attenname1 = 'Campbell (2003), Hor + SigmaEps, CEUS-Hybrid, Hard Rock'
         lnY = lnY + sigmaEps
      endif

c ******* PE&A CEUS Models *******
C Silva et al. (2002) 2 Corner, Rock
C     Model Number = 401
      if (jcalc .eq. 401) then
         call S06_PEA2C ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 2-Corner, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 2 Corner-Saturation, Rock
C     Model Number = 402
      if (jcalc .eq. 402) then
         call S06_PEA2CS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 2-Corner-Sat, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-High, Rock
C     Model Number = 403
      if (jcalc .eq. 403) then
         call S06_PEA1CVH ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-High, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Medium, Rock
C     Model Number = 404
      if (jcalc .eq. 404) then
         call S06_PEA1CVM ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Med, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Low, Rock
C     Model Number = 405
      if (jcalc .eq. 405) then
         call S06_PEA1CVL ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Low, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High, Rock
C     Model Number = 406
      if (jcalc .eq. 406) then
         call S06_PEA1CCH ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Medium, Rock
C     Model Number = 407
      if (jcalc .eq. 407) then
         call S06_PEA1CCM ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low, Rock
C     Model Number = 408
      if (jcalc .eq. 408) then
         call S06_PEA1CCL ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High-Sat, Rock
C     Model Number = 409
      if (jcalc .eq. 409) then
         call S06_PEA1CCHS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High-Sat, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Med-Sat, Rock
C     Model Number = 410
      if (jcalc .eq. 410) then
         call S06_PEA1CCMS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med-Sat, Horizontal'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low-Sat, Rock
C     Model Number = 411
      if (jcalc .eq. 411) then
         call S06_PEA1CCLS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low-Sat, Horizontal'
      endif

C New PE&A CEUS Gulf Models
C Silva et al. (2002) 2 Corner, Rock, Gulf
C     Model Number = 501
      if (jcalc .eq. 501) then
         call S06_PEAG2C ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 2-Corner, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 2 Corner-Saturation, Rock, Gulf
C     Model Number = 502
      if (jcalc .eq. 502) then
         call S06_PEAG2CS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 2-Corner-Sat, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-High, Rock, Gulf
C     Model Number = 503
      if (jcalc .eq. 503) then
         call S06_PEAG1CVH ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-High, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Medium, Rock, Gulf
C     Model Number = 504
      if (jcalc .eq. 504) then
         call S06_PEAG1CVM ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Med, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Low, Rock, Gulf
C     Model Number = 505
      if (jcalc .eq. 505) then
         call S06_PEAG1CVL ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Low, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High, Rock, Gulf
C     Model Number = 506
      if (jcalc .eq. 506) then
         call S06_PEAG1CCH ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Medium, Rock, Gulf
C     Model Number = 507
      if (jcalc .eq. 507) then
         call S06_PEAG1CCM ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low, Rock, Gulf
C     Model Number = 508
      if (jcalc .eq. 508) then
         call S06_PEAG1CCL ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High-Sat, Rock, Gulf
C     Model Number = 509
      if (jcalc .eq. 509) then
         call S06_PEAG1CCHS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High-Sat, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Med-Sat, Rock, Gulf
C     Model Number = 510
      if (jcalc .eq. 510) then
         call S06_PEAG1CCMS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med-Sat, Hor. Gulf'
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low-Sat, Rock, Gulf
C     Model Number = 511
      if (jcalc .eq. 511) then
         call S06_PEAG1CCLS ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low-Sat, Hor. Gulf'
      endif

C *****  Misc Models ******
c     McVerry et al (1993) new zealand
C     Model Number = 300
      if ( jcalc .eq. 300 .and. specT .eq. 0.0 ) then
        call S02_mcverry93 ( mag, rupDist, lnY, sigma, attenName1,
     1       ftype )
          period2 = 0.0
          iflag = 0
      elseif ( jcalc .eq. 300 .and. specT .ne. 0.0) then
          write (*,*) 'McVerry (1993), Horizontal, Rock'
          write (*,*) 'not defined for spectral acclereation!!!'
          write (*,*) 'Check input file.'
          stop 99
      endif

c     fukushima (1990) rock
C     Model Number = 301
      if ( jcalc .eq. 301 .and. specT .eq. 0.0 ) then
          call S02_fukushima90 ( mag, rupDist, lnY, sigma, attenName1)
          period2 = 0.0
          iflag = 0
      elseif (jcalc .eq. 301 .and. specT .ne. 0.0) then
          write (*,*) 'Fukushima (1990), Horizontal, Rock'
          write (*,*) 'not defined for spectral acclereation!!!'
          write (*,*) 'Check input file.'
          stop 99
      endif

c     Loh high speed rail (New Joyner-Boore form)
C     Model Number = 302
      if ( jcalc .eq. 302 .and. specT .eq. 0.0 ) then
         call S02_HighSpeedRail ( mag, rupDist, lnY, sigma,
     1                      attenName1, period2 )
         iflag = 0
      elseif (jcalc .eq. 302 .and. specT .ne. 0.0 ) then
         write (*,*) 'Loh High Speed Rail, Horizontal, Rock'
         write (*,*) 'not defined for spectral acclereation!!!'
         write (*,*) 'Check input file.'
         stop 99
      endif

c     New Loh (1996) model (unpublished)
C     Model Number = 303
      if ( jcalc .eq. 303 .and. specT .eq. 0.0) then
         call S02_Loh96 ( mag, rupDist, lnY, sigma, attenName1, period2)
         iflag = 0
      elseif (jcalc .eq. 303 .and. specT .ne. 0.0 ) then
         write (*,*) 'Loh (1996), Horizontal, Rock'
         write (*,*) 'not defined for spectral acclereation!!!'
         write (*,*) 'Check input file.'
         stop 99
      endif

c ******* Ambraseys et al 2005 Model *********
C     Model Number = 601
      if ( jcalc .eq. 601 ) then
         call S03_Ambraseys_2005 ( mag, jbDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
         attenname1 = 'Ambraseys_et_al_2005_Hor'
      endif


C     **** BCHydro SCR attenuation models adjusted for Vs=760m/sec ****

c     Atkinson and Boore (2006) - Horizontal, CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1020
      if ( jcalc .eq. 1020 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Atkinson&Boore 2006, EUS, BCH Amps for Vs760m/s'
         lnY = lnSa
      endif

c     Atkinson and Boore (2006) with Atikinson (2010) stress drop adjustment - Horizontal, EUS BCH Amps for Vs760m/s
C     Model Number = 1050
      if ( jcalc .eq. 1050 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 10**(3.45 - 0.2*mag)
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006/Atkinson2010, EUS, BCH Amps for Vs760m/s'
         LnY = lnY + sdscale*SF2*alog(10.0)
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c     Atkinson and Boore (2006) with 2x stress drop adjustment (280bars) - Horizontal, EUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1300
      if ( jcalc .eq. 1300 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 280.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (2x) StressDrop, EUS, BCH Amps for Vs760'
         LnY = lnY + sdscale*SF2*alog(10.0)
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c     Atkinson and Boore (2006) with 0.5x stress drop adjustment (70bars) - Horizontal, EUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1310
      if ( jcalc .eq. 1310 ) then
         call S06_AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.
         call S06_AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 70.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (0.5x) StressDrop, EUS, BCH Amps for Vs760'
         LnY = lnY + sdscale*SF2*alog(10.0)
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal, BCH Amps for Vs760m/s
C     Model Number = 1200
      if (jcalc .eq. 1200) then
         call S06_CHY03 ( mag, rupdist, lnY, sigma, specT,
     1                  period2,iflag )
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Campbell (2003), Hor., CEUS-Hybrid, BCH Amps for Vs760m/s'
         lnY = lnSa
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal minus SigmaEps, BCH Amps for Vs760m/s
C     Model Number = 1210s
      if (jcalc .eq. 1210) then
         call S06_CHY03 ( mag, rupdist, lnY, sigma, specT,
     1                  period2,iflag )
         call S06_CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )
         attenname1 = 'Campbell (2003), Hor - SigmaEps, CEUS-Hybrid, BCH Amps for Vs760m/s'
         lnY = lnY - sigmaEps
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal plus SigmaEps, BCH Amps for Vs760m/s
C     Model Number = 1220
      if (jcalc .eq. 1220) then
         call S06_CHY03 ( mag, rupdist, lnY, sigma, specT,
     1                  period2,iflag )
         call S06_CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )
         attenname1 = 'Campbell (2003), Hor + SigmaEps, CEUS-Hybrid, BCH Amps for Vs760m/s'
         lnY = lnY + sigmaEps
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-High, BCH Amps for Vs760m/s
C     Model Number = 4030
      if (jcalc .eq. 4030) then
         call S06_PEA1CVH ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-High, Hor,BCH Amps for Vs760m/s'
         lnY = lnSA
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Medium, BCH Amps for Vs760m/s
C     Model Number = 4040
      if (jcalc .eq. 4040) then
         call S06_PEA1CVM ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Med, Hor,BCH Amps for Vs760m/s'
         lnY = lnSa
      endif

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Low, BCH Amps for Vs760m/s
C     Model Number = 4050
      if (jcalc .eq. 4050) then
         call S06_PEA1CVL ( mag, jbdist, lnY, sigma, specT,
     1                  attenName1, period2,iflag )
         call S05_BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Low, Hor,BCH Amps for Vs760m/s'
         lnY = lnSa
      endif

C     **** End of BCHydro SCR attenuation models adjusted for Vs=760m/sec ****

C     Akkar and Cayan (2010) - A Local Ground-Motion Predictive Model for Turkey, and
C           its Comparison with Other Regional and Global Ground-Motion Models
C     BSSA: Vol. 100, No.6, pp. 2978-2995.
C     Model Number = 150
      if ( jcalc .eq. 150 ) then
         call S02_AC_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pga4nl )
         attenname1 = 'Akkar&Cagan_2010_Hor'
       endif

C     Akkar and Bommer (2010) - Empirical Equations for the Prediction of PGA, PGV,
C             and Spectral Accelerations in Europe, the Mediterranean Region and
C             the Middle East
C     SRL: March/April, Vol.81, pp. 195-206
C     Rock Site conditions Vs>750m/s
C     Model Number = 151
      if ( jcalc .eq. 151 ) then
         Ss = 0.0
         Sa = 0.0
         call S02_AB_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Ss, Sa )
         attenname1 = 'Akkar&Bommer_2010_Hor, Rock'
       endif

C     Stiff Soils Site conditions 360<Vs<750m/s
C     Model Number = 152
      if ( jcalc .eq. 152 ) then
         Ss = 0.0
         Sa = 1.0
         call S02_AB_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Ss, Sa )
         attenname1 = 'Akkar&Bommer_2010_Hor, Stiff Soil'
       endif
C     Soft Soil Site conditions Vs<360m/s
C     Model Number = 153
      if ( jcalc .eq. 153 ) then
         Ss = 1.0
         Sa = 0.0
         call S02_AB_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Ss, Sa )
         attenname1 = 'Akkar&Bommer_2010_Hor, Soft Soil'
       endif

C     Akkar, Sandikkaya, and Bommer (2013) - Empirical ground-motion models for
C             point- and extended-source crustal earthquake scenarios in
C             Europe and the Middle East
C     Bull Earthquake Engineering: May 31, 2013
C         Applicable Range:
C            Mw = 4 - 8
C            Distance < 200km
C            Vs = 150 - 1200 (Note site response is constant for Vs>1000)
C                Reference Vs = 750m/s
C            Horizontal is Geomean
C            Defined for PGA, 0.01 - 4.0 sec and PGV
C     Model Number = 154
      if ( jcalc .eq. 154 ) then
         call S02_ASB_2013 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Vs, phi, tau )
         attenname1 = 'Akkar,Sandikkaya&Bommer_2013_Hor'
       endif





c ******* Bradley (2010) *********
C     Bradley (2010) - NZ-Specific Psuedo-Spectral Acceleration Ground
C         Motion Prediction Equation Based on Foreign Models
C     University of Canterbury Research Report 2010-03
C     Bradley 2010 - Horizontal, estimated Vs30
C     Model Number = 160
      if ( jcalc .eq. 160 ) then
         vs30_class = 0
         call S02_Bradley_2010 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Bradley-2010-Hor, Estimated Vs30'
       endif

C     Bradley 2010 - Horizontal, measured Vs30
C     Model Number = 161
      if ( jcalc .eq. 161 ) then
         vs30_class = 1
         call S02_Bradley_2010 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Bradley-2010-Hor, Measured Vs30'
       endif

c ******* McVerry et al. (2006) Crustal Events, Horizontal *********
C     McVerry, G. H., J. X. Zhao, N.A. Abrahamson and P.G. Somerville (2006)
C         New Zealand Acceleration Response Spectrum Attenuation
C         Reations for Crustal and Subduction Zone Earthquakes,
C         Bulletin New Zealand Society for Earthquake Engineering, Vol. 39,
C         No. 1, pp. 1 - 58.
C     Model Number = 140, Site Class A/B
      if ( jcalc .eq. 140 ) then
         Sc = 0.0
         Sd = 0.0
         call S02_McVerry_Crustal_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype,
     3                     hwflag, Sc, Sd)
         attenname1 = 'McVerryetal-Crustal-2006, Hor, Site A/B'
       endif
C     Model Number = 141, Site Class C
      if ( jcalc .eq. 141 ) then
         Sc = 1.0
         Sd = 0.0
         call S02_McVerry_Crustal_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype,
     3                     hwflag, Sc, Sd)
         attenname1 = 'McVerryetal-Crustal-2006, Hor, Site C'
       endif
C     Model Number = 142, Site Class D
      if ( jcalc .eq. 142 ) then
         Sc = 0.0
         Sd = 1.0
         call S02_McVerry_Crustal_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype,
     3                     hwflag, Sc, Sd)
         attenname1 = 'McVerryetal-Crustal-2006, Hor, Site D'
       endif

c ******* McVerry et al. (2006) Subduction Events, Horizontal *********
C     McVerry, G. H., J. X. Zhao, N.A. Abrahamson and P.G. Somerville (2006)
C        "New Zealand Acceleration Response Spectrum Attenuation
C         Reations for Crustal and Subduction Zone Earthquakes",
C         Bulletin New Zealand Society for Earthquake Engineering, Vol. 39,
C         No. 1, pp. 1 - 58.
C     Model Number = 143, Site Class A/B
      if ( jcalc .eq. 143 ) then
         Sc = 0.0
         Sd = 0.0
         call S02_McVerry_Subduction_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype,
     3                     depthtop, dipavgd, rupwidth, depth, Sc, Sd )

         attenname1 = 'McVerryetal-Subduction-2006, Hor, Site A/B'
       endif
C     Model Number = 144, Site Class C
      if ( jcalc .eq. 144 ) then
         Sc = 1.0
         Sd = 0.0
         call S02_McVerry_Subduction_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype,
     3                     depthtop, dipavgd, rupwidth, depth, Sc, Sd )
         attenname1 = 'McVerryetal-Subduction-2006, Hor, Site C'
       endif
C     Model Number = 145, Site Class D
      if ( jcalc .eq. 145 ) then
         Sc = 0.0
         Sd = 1.0
         call S02_McVerry_Subduction_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype,
     3                     depthtop, dipavgd, rupwidth, depth, Sc, Sd )
         attenname1 = 'McVerryetal-Subduction-2006, Hor, Site D'
       endif

c ******* Bindi et al. (2009) Crustal Events, Horizontal *********
C     Bindi, D., L. Luzi, M. Massa, and F. Pacor (2009). Horizontal and Vertical
C          ground motion prediction equaitons derived from the Italian
C          Accelerometric Archive (ITACA), Bull. Earthquake Eng.,
C          DOI 10.1007/s10518-009-9130-9
C      Note: Only the horizontal model is currently coded and the coefficients
C            independent of fault mechanism as preferred by the authors. Coefficients
C            are for the model as a function of JBDist.
C     Model Number = 95, Rock, Horizontal
      if ( jcalc .eq. 95 ) then
         Sr = 1.0
         Ss = 0.0
         Sd = 0.0
         call S02_Bindi_Hor_2009 ( mag, jbdist, specT,
     1                     period2, lnY, sigma, iflag, Sr, Ss, Sd )
         attenname1 = 'Bindietal-Hor-2009, Rock'
       endif
C     Model Number = 96, Shallow Alluvium
      if ( jcalc .eq. 96 ) then
         Sr = 0.0
         Ss = 1.0
         Sd = 0.0
         call S02_Bindi_Hor_2009 ( mag, jbdist, specT,
     1                     period2, lnY, sigma, iflag, Sr, Ss, Sd )
         attenname1 = 'Bindietal-Hor-2009, Shallow Alluvium'
       endif
C     Model Number = 97, Deep Alluvium
      if ( jcalc .eq. 97 ) then
         Sr = 0.0
         Ss = 0.0
         Sd = 1.0
         call S02_Bindi_Hor_2009 ( mag, jbdist, specT,
     1                     period2, lnY, sigma, iflag, Sr, Ss, Sd )
         attenname1 = 'Bindietal-Hor-2009, Deep Alluvium'
       endif

c ******* Bindi et al. (2011) Crustal Events, Horizontal *********
C     Bindi, D., F. Pacor, L. Luzi, R. Puglia, M. Massa, G. Ameri and R. Paolucci (2011).
C          Ground motion prediction equations derived from the Italian
C          strong motion database.
C          Bull. Earthquake Eng, Vol. 9, pp. 1899-1920. uake Eng.,
C          DOI 10.1007/s10518-011-9313-z
C       Dataset is for 4.1 < M < 6.9 and distance less than 200km.
C       Model is considered an update of the Bindi (2009) model
C          change in site conditions classes from 2009 to 2011 with new model
C          being consistent with Eurocode 8 Site Classes:
C             Class A: Vs>800m/s
C             Class B: Vs>360-800m/s
C             Class C: Vs>180-360m/s
C             Class D: Vs<180-800m/s
C             Class E: 5-20m of C or D type alluvium underlain by stiffer material with Vs>800m/s
C
C     Model Number = 195, Horizontal, Class A (Vs>800m/s)
      if ( jcalc .eq. 195 ) then
         SCa = 1.0
         SCb = 0.0
         SCc = 0.0
         SCd = 0.0
         SCe = 0.0
         call S02_Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class A(Vs>800m/s)'
       endif
C     Model Number = 196, Class B (Vs=360-800m/s)
      if ( jcalc .eq. 196 ) then
         SCa = 0.0
         SCb = 1.0
         SCc = 0.0
         SCd = 0.0
         SCe = 0.0
         call S02_Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class B(360<Vs<800m/s)'
       endif

C     Model Number = 197, Class C (Vs=180-360m/s)
      if ( jcalc .eq. 197) then
         SCa = 0.0
         SCb = 0.0
         SCc = 1.0
         SCd = 0.0
         SCe = 0.0
         call S02_Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class C(180<Vs<360m/s)'
       endif

C     Model Number = 198, Class D (Vs<180)
      if ( jcalc .eq. 198) then
         SCa = 0.0
         SCb = 0.0
         SCc = 0.0
         SCd = 1.0
         SCe = 0.0
         call S02_Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class D(Vs<180)'
       endif

C     Model Number = 199, Class E (Vs=180-360m/s)
      if ( jcalc .eq. 199) then
         SCa = 0.0
         SCb = 0.0
         SCc = 0.0
         SCd = 0.0
         SCe = 1.0
         call S02_Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class E(5-20m Class D/E over ClassA)'
       endif

c ******* Bindi et al. (2013) Crustal Events, Horizontal *********
C     Bindi, D., M. Massa, L. Luzi, G. Ameri, R. Puglia, and P. Augliera (2013).
C           Pan-European Ground-Motion Prediction Equations for the Average
C           Horizontal Component of PGA, PGV, and 5%-Damped PSA at Spectral
C           Periods up to 3.0 s using the RESORCE dataset
C     Several versions of the model are developed in terms of Vs or Eurocode
C           Site Classes and distance metrics of RJB and Rhypo. The version
C           code in the code is for the Vs and RJB parameters.
C     Recommended application range:
C           Magnitude = 4 - 7.6
C           Distance < 300 km
C           Hypocenteral Depths < 35km
C
C     Model Number = 295, Horizontal, Rjb, Vs
      if ( jcalc .eq. 295 ) then
         call S02_Bindi_Hor_2013 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, vs, phi, tau )
         attenname1 = 'Bindietal-Hor-2013, Rjb, Vs'
       endif

c ******* Graizer and Kalkan (Nov. 2012) *********
C     This is an unpublished update to the Jan/Feb. 2011 SRL version of the GMPE.
C          Fortran code for this model was obtained from V. Graizer and it is
C          currently being submitted for publication.
C     Note that currently a fixed Q values typical for California is set.
C     Currently depthvs15 is set equal to depthvs10 value
C     Model Number = 90
      if ( jcalc .eq. 90 ) then
         Q0 = 150.0
         call S02_GK_Nov2012 ( mag, RupDist, specT, ftype,
     1                     period2, lnY, sigma, iflag, Vs, Q0, depthvs15 )
         attenname1 = 'Graizer&Kalkan, Nov.2012'
       endif
      if ( jcalc .eq. 91 ) then
         Q0 = 75.0
         call S02_GK_Nov2012 ( mag, RupDist, specT, ftype,
     1                     period2, lnY, sigma, iflag, Vs, Q0, depthvs15 )
         attenname1 = 'Graizer&Kalkan, Nov.2012, Q0=75'
       endif
      if ( jcalc .eq. 92 ) then
         Q0 = 300.0
         call S02_GK_Nov2012 ( mag, RupDist, specT, ftype,
     1                     period2, lnY, sigma, iflag, Vs, Q0, depthvs15 )
         attenname1 = 'Graizer&Kalkan, Nov.2012, Q0=300'
       endif


c ******* DCPP Common Function Form Models *********
C     Preliminary Models for DCPP for Workshop 3, March 2014
C     Model Numbers 8001, ASK form
      if ( jcalc .eq. 8001 ) then
         call S02_DCPP_CommonASK ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model, ASK'
       endif
C     Model Numbers 8002, BSSA form
      if ( jcalc .eq. 8002 ) then
         call S02_DCPP_CommonBSSA ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model, BSSA'
       endif
C     Model Numbers 8003, Common Model 001
      if ( jcalc .eq. 8003 ) then
         call S02_DCPP_Common001 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model001'
       endif
C     Model Numbers 8004, Common Model 002
      if ( jcalc .eq. 8004 ) then
         call S02_DCPP_Common002 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model002'
       endif
C     Model Numbers 8005, Common Model 003
      if ( jcalc .eq. 8005 ) then
         call S02_DCPP_Common003 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model003'
       endif
C     Model Numbers 8006, Common Model 004
      if ( jcalc .eq. 8006 ) then
         call S02_DCPP_Common004 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model004'
       endif
C     Model Numbers 8007, Common Model 005
      if ( jcalc .eq. 8007 ) then
         call S02_DCPP_Common005 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model005'
       endif
C     Model Numbers 8008, Common Model 006
      if ( jcalc .eq. 8008 ) then
         call S02_DCPP_Common006 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model006'
       endif
C     Model Numbers 8009, Common Model 007
      if ( jcalc .eq. 8009 ) then
         call S02_DCPP_Common007 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model007'
       endif
C     Model Numbers 8010, Common Model 008
      if ( jcalc .eq. 8010 ) then
         call S02_DCPP_Common008 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model008'
       endif
C     Model Numbers 8011, Common Model 009
      if ( jcalc .eq. 8011 ) then
         call S02_DCPP_Common009 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model009'
       endif

c ******* PVNGS Common Function Form Models *********
C     Preliminary Models for PVNGS for Workshop 3, March 2014

C     Model Numbers 9001 Common Model ASK
      if ( jcalc .eq. 9001 ) then
         call S02_PVNGS_CommonASK ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model ASK'
      endif
C     Model Numbers 9002 Common Model ASK
      if ( jcalc .eq. 9002 ) then
         call S02_PVNGS_CommonBindi ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model Bindi'
      endif
C     Model Numbers 9003 Common Model ASK
      if ( jcalc .eq. 9003 ) then
         call S02_PVNGS_Common001 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model001'
      endif
C     Model Numbers 9004 Common Model ASK
      if ( jcalc .eq. 9004 ) then
         call S02_PVNGS_Common002 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model002'
      endif
C     Model Numbers 9005 Common Model ASK
      if ( jcalc .eq. 9005 ) then
         call S02_PVNGS_Common003 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model003'
      endif
C     Model Numbers 9006 Common Model ASK
      if ( jcalc .eq. 9006 ) then
         call S02_PVNGS_Common004 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model004'
      endif
C     Model Numbers 9007 Common Model ASK
      if ( jcalc .eq. 9007 ) then
         call S02_PVNGS_Common005 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model005'
      endif
C     Model Numbers 9008 Common Model ASK
      if ( jcalc .eq. 9008 ) then
         call S02_PVNGS_Common006 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model006'
      endif
C     Model Numbers 9009 Common Model ASK
      if ( jcalc .eq. 9009 ) then
         call S02_PVNGS_Common007 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model007'
      endif
C     Model Numbers 9010 Common Model ASK
      if ( jcalc .eq. 9010 ) then
         call S02_PVNGS_Common008 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model008'
      endif
C     Model Numbers 9011 Common Model ASK
      if ( jcalc .eq. 9011 ) then
         call S02_PVNGS_Common009 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model009'
      endif




C     SWUS Common Functional Form as a function of Rrup
C         10,000 < jcalc < 11,000
      if ( jcalc .gt. 10000 .and. jcalc .lt. 11000 ) then
         coefcountrrup = jcalc - 10000
         if (coefcountrrup .lt. 0) then
            write (*,*) 'Incorrect jcalc for SWUS Common Functional Model Rrup!!!'
            write (*,*) 'Check input file.'
            Stop 99
         endif
         call S02_SWUS_CFRrup ( mag, RupDist, jbDist, depthtop, ftype, dipavgd, RupWidth, Rx, HWFlag,
     1           specT, lnY, sigma, iflag, cfcoefrrup, coefcountrrup, phi, tau )
         attenname1 = 'SWUS Common Function Model-Rrup'
      endif


C     SWUS Common Functional Form as a function of Rjb
C         11,000 < jcalc < 12,000
      if ( jcalc .gt. 11000 .and. jcalc .lt. 12000 ) then
         coefcountrjb = jcalc - 11000
         if (coefcountrjb .lt. 0) then
            write (*,*) 'Incorrect jcalc for SWUS Common Functional Model Rjb!!!'
            write (*,*) 'Check input file.'
            Stop 99
         endif
         call S02_SWUS_CFRjb ( mag, RupDist, jbDist, depthtop, ftype, dipavgd, RupWidth, Rx, HWFlag,
     1           specT, lnY, sigma, iflag, cfcoefrjb, coefcountrjb, phi, tau )
         attenname1 = 'SWUS Common Function Model-Rjb'
      endif

C     SWUS Common Functional Form as a function of Rrup - for DCPP
C         12,000 < jcalc < 13,000
      if ( jcalc .gt. 12000 .and. jcalc .lt. 13000 ) then
         coefcountrrup = jcalc - 12000
         if (coefcountrrup .lt. 0) then
            write (*,*) 'Incorrect jcalc for SWUS Common Functional Model Rrup!!!'
            write (*,*) 'Check input file.'
            Stop 99
         endif
         call S02_SWUS_CFRrup_DCPP ( mag, RupDist, jbDist, depthtop, ftype, dipavgd, RupWidth, Rx, HWFlag,
     1           specT, lnY, sigma, iflag, cfcoefrrup, coefcountrrup, phi, tau )
         attenname1 = 'SWUS DCPP Common Function Model-Rrup'
      endif


C ******* PEER NGA-West2 Attenuation Models with Additional Mag Scaling Uncertainty ****
C ******* For use with Zones 1, 2, and 3  with path term accounted for *****
C         (Jcalcs start with 8000s)

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     Central Mag Uncertainty
C     Model Number = 8787
      if ( jcalc .eq. 8787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-Cent-EstVs'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     Low Mag Uncertainty
C     Model Number = 8788
      if ( jcalc .eq. 8788 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-Low-EstVs'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     High Mag Uncertainty
C     Model Number = 8789
      if ( jcalc .eq. 8789 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-High-EstVs'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
      endif


c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 8922
      if ( jcalc .eq. 8922 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-Cent, No Basin'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 8923
      if ( jcalc .eq. 8923 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-Low, No Basin'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     High Mag Uncertainty
C     Model Number = 8924
      if ( jcalc .eq. 8924 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-High, No Basin'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 8836
      if ( jcalc .eq. 8836 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-Cent'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 8837
      if ( jcalc .eq. 8837 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-Low'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     High Mag Uncertainty
C     Model Number = 8838
      if ( jcalc .eq. 8838 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-High'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif



c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 8797
      if ( jcalc .eq. 8797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-Cent,Est Vs30m'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 8798
      if ( jcalc .eq. 8798 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-Low,Est Vs30m'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 8799
      if ( jcalc .eq. 8799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-High,Est Vs30m'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif


c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 8910
      if ( jcalc .eq. 8910 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Cent'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Cent'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 8911
      if ( jcalc .eq. 8911 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Low'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Low'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 8912
      if ( jcalc .eq. 8912 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-High'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-High'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         if (mag .gt. 7.0) then
            if (specT .lt. 1.0) then
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0))**2.0 - 0.083**2.0)
            else
               faddmag = c1*sqrt( (0.083+0.056*(mag-7.0)+0.0171*alog(specT))**2.0 - 0.083**2.0)
            endif
            lnY = lnY + faddmag
         endif
       endif

C ******* PEER NGA-West2 Attenuation Models with Additional Mag Scaling Uncertainty ****
C ******* For use with Zones 1, 2, and 3  with no path term accounted for *****
C         (Jcalcs start with 9000s)

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     Central Mag Uncertainty
C     Model Number = 9787
      if ( jcalc .eq. 9787 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-Cent-EstVs'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     Low Mag Uncertainty
C     Model Number = 9788
      if ( jcalc .eq. 9788 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-Low-EstVs'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     High Mag Uncertainty
C     Model Number = 9789
      if ( jcalc .eq. 9789 ) then
         regionflag = 0
         msasflag = 0
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop).
         vs30_Class = 0
         call S09_ASK_NGAWest2_2013 ( mag, dipavgd, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag,
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-High-EstVs'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
      endif


c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 9922
      if ( jcalc .eq. 9922 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-Cent, No Basin'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 9923
      if ( jcalc .eq. 9923 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-Low, No Basin'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     High Mag Uncertainty
C     Model Number = 9924
      if ( jcalc .eq. 9924 ) then
         regionflag = 0
         basinflag = 0
         call S09_BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-High, No Basin'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 9836
      if ( jcalc .eq. 9836 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-Cent'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 9837
      if ( jcalc .eq. 9837 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-Low'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     High Mag Uncertainty
C     Model Number = 9838
      if ( jcalc .eq. 9838 ) then
         regionflag = 0
         call S09_CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs,
     2                    depthTop, D25, dipavgd, depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-High'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif



c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 9797
      if ( jcalc .eq. 9797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-Cent,Est Vs30m'

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 9798
      if ( jcalc .eq. 9798 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-Low,Est Vs30m'

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 9799
      if ( jcalc .eq. 9799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call S09_CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag,
     2                     vs, dipavgd, Depthtop, Ftype,
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-High,Est Vs30m'

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif


c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 9910
      if ( jcalc .eq. 9910 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Cent'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Cent'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif

C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 9911
      if ( jcalc .eq. 9911 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Low'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Low'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif

C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 9912
      if ( jcalc .eq. 9912 ) then
         if (vs .ge. 450.0) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-High'
         elseif (vs .gt. 1200) then
            call S09_I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-High'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif

C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag
       endif

C     Total Sigma Model for SWUS
C     Note that these sigma models should only be called with given scalc number.
C     The use of these values as jcalcs will be incorrect as the median ground motions is kept
C     at the large value of 1.0e10.

C     SWUS Total Sigma: DCPP Central
C     Model Number = 13001
      if ( jcalc .eq. 13001 ) then
         attenname1 = 'SWUS Total Sigma DCPP-Central'
         call S32_SWUS_Sigma_DCPP_Cen ( mag, specT, sigma, iflag )
C     Kepp median ground motions large since this is only for sigma model
         lnY = 1.0e10
       endif

C     SWUS Total Sigma: DCPP Low
C     Model Number = 13002
      if ( jcalc .eq. 13002 ) then
         attenname1 = 'SWUS Total Sigma DCPP-Low'
         call S32_SWUS_Sigma_DCPP_Low ( mag, specT, sigma, iflag )
C     Kepp median ground motions large since this is only for sigma model
         lnY = 1.0e10
       endif

C     SWUS Total Sigma: DCPP High
C     Model Number = 13003
      if ( jcalc .eq. 13003 ) then
         attenname1 = 'SWUS Total Sigma DCPP-High'
         call S32_SWUS_Sigma_DCPP_High ( mag, specT, sigma, iflag )
C     Kepp median ground motions large since this is only for sigma model
         lnY = 1.0e10
       endif

C     SWUS PHISS Sigma: PhiSS_CA1 - Low
C     Model Number = 13004
      if ( jcalc .eq. 13004 ) then
         attenname1 = 'SWUS phiSS_CA1 Low'
         iBranch = 1
         call S32_SWUS_PHISS_CA1 ( mag, specT, phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_CA1 - central
C     Model Number = 13005
      if ( jcalc .eq. 13005 ) then
         attenname1 = 'SWUS phiSS_CA1 Central'
         iBranch = 2
         call S32_SWUS_PHISS_CA1 ( mag, specT, phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_CA1 - high
C     Model Number = 13006
      if ( jcalc .eq. 13006 ) then
         attenname1 = 'SWUS phiSS_CA1 HIgh'
         iBranch = 3
         call S32_SWUS_PHISS_CA1 ( mag, specT, phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_CA2 - Low
C     Model Number = 13007
      if ( jcalc .eq. 13007 ) then
         attenname1 = 'SWUS phiSS_CA2 Low'
         iBranch = 1
         call S32_SWUS_PHISS_CA2 ( mag, specT, phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_CA2 - central
C     Model Number = 13008
      if ( jcalc .eq. 13008 ) then
         attenname1 = 'SWUS phiSS_CA2 Central'
         iBranch = 2
         call S32_SWUS_PHISS_CA2 ( mag, specT, phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_CA2 - high
C     Model Number = 13009
      if ( jcalc .eq. 13009 ) then
         attenname1 = 'SWUS phiSS_CA2 HIgh'
         iBranch = 3
         call S32_SWUS_PHISS_CA2 ( mag, specT, phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_Global_R50 - Low
C     Model Number = 13010
      if ( jcalc .eq. 13010 ) then
         attenname1 = 'SWUS phiSS_Global_R50 Low'
         iBranch = 1
         call S32_SWUS_PHISS_Global_R50 ( phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_Global_R50 - central
C     Model Number = 13011
      if ( jcalc .eq. 13011 ) then
         attenname1 = 'SWUS phiSS_Global_R50 Central'
         iBranch = 2
         call S32_SWUS_PHISS_Global_R50 ( phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif

C     SWUS PHISS Sigma: PhiSS_Global_R50 - high
C     Model Number = 13012
      if ( jcalc .eq. 13012 ) then
         attenname1 = 'SWUS phiSS_Global_R50 HIgh'
         iBranch = 3
         call S32_SWUS_PHISS_Global_R50 ( phiSS, iflag, iBranch )

C        set dummy value for median (this is used only for sigma)
         lnY = 1.0e10
         write (*,'( 2x,''phiSS only, need to add tau'')')
         stop 99
       endif




C********************************************************************
C ***** EPRI Updated (2013) GMPE Models *****
C    Nomeclature is as follows:
C      Number 1-4 = "2013" EPRI Update 2013 Mid Continent Models
C      Number 5-6 = Median model cases
C                     "01" = Cluster01-Low,  Functional Model 1&3
C                     "02" = Cluster01-Med,  Functional Model 1&3
C                     "03" = Cluster01-High, Functional Model 1&3
C                     "04" = Cluster02-Low,  Functional Model 2
C                     "05" = Cluster02-Med,  Functional Model 2
C                     "06" = Cluster02-High, Functional Model 2
C                     "07" = Cluster03-Low,  Functional Model 1&3
C                     "08" = Cluster03-Med,  Functional Model 1&3
C                     "09" = Cluster03-High, Functional Model 1&3
C                     "10" = Cluster04-Low (Rift),  Functional Model 4
C                     "11" = Cluster04-Med (Rift),  Functional Model 4
C                     "12" = Cluster04-High (Rift), Functional Model 4
C                     "13" = Cluster04-Low (NonRift),  Functional Model 4
C                     "14" = Cluster04-Med (NonRift),  Functional Model 4
C                     "15" = Cluster04-High (NonRift), Functional Model 4
C********************************************************************
C

C *** Cluster 01-Low, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock ***
C     Model Number = 201301
      if (jcalc .eq. 201301) then
         call S06_EPRI13C1Low ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster01-Low,MidC, Hor, HardRock'
      endif

C *** Cluster 01-Med, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock ***
C     Model Number = 201302
      if (jcalc .eq. 201302) then
         call S06_EPRI13C1Med ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster01-Med,MidC, Hor, HardRock'
      endif

C *** Cluster 01-High, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock ***
C     Model Number = 201303
      if (jcalc .eq. 201303) then
         call S06_EPRI13C1High ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster01-High,MidC, Hor, HardRock'
      endif

C *** Cluster 02-Low, Mid-Continent: Functional Model 2, Horizontal, CEUS Hard Rock ***
C     Model Number = 201304
      if (jcalc .eq. 201304) then
         call S06_EPRI13C2Low ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster02-Low,MidC, Hor, HardRock'
      endif

C *** Cluster 02-Med, Mid-Continent: Functional Model 2, Horizontal, CEUS Hard Rock ***
C     Model Number = 201305
      if (jcalc .eq. 201305) then
         call S06_EPRI13C2Med ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster02-Med,MidC, Hor, HardRock'
      endif

C *** Cluster 02-High, Mid-Continent: Functional Model 2, Horizontal, CEUS Hard Rock ***
C     Model Number = 201306
      if (jcalc .eq. 201306) then
         call S06_EPRI13C2High ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster02-High,MidC, Hor, HardRock'
      endif

C *** Cluster 03-Low, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock ***
C     Model Number = 201307
      if (jcalc .eq. 201307) then
         call S06_EPRI13C3Low ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster03-Low,MidC, Hor, HardRock'
      endif

C *** Cluster 03-Med, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock ***
C     Model Number = 201308
      if (jcalc .eq. 201308) then
         call S06_EPRI13C3Med ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster03-Med,MidC, Hor, HardRock'
      endif

C *** Cluster 03-High, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock ***
C     Model Number = 201309
      if (jcalc .eq. 201309) then
         call S06_EPRI13C3High ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster03-High,MidC, Hor, HardRock'
      endif

C *** Cluster 04-Low (Rift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock ***
C     Model Number = 201310
      if (jcalc .eq. 201310) then
         call S06_EPRI13C4RLow ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Low-Rift,MidC, Hor, HardRock'
      endif

C *** Cluster 04-Med (Rift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock ***
C     Model Number = 201311
      if (jcalc .eq. 201311) then
         call S06_EPRI13C4RMed ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Med-Rift,MidC, Hor, HardRock'
      endif

C *** Cluster 04-High (Rift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock ***
C     Model Number = 201312
      if (jcalc .eq. 201312) then
         call S06_EPRI13C4RHigh ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-High-Rift,MidC, Hor, HardRock'
      endif


C *** Cluster 04-Low (NonRift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock ***
C     Model Number = 201313
      if (jcalc .eq. 201313) then
         call S06_EPRI13C4NRLow ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Low-NonRift,MidC, Hor, HardRock'
      endif

C *** Cluster 04-Med (NonRift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock ***
C     Model Number = 201314
      if (jcalc .eq. 201314) then
         call S06_EPRI13C4NRMed ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Med-NonRift,MidC, Hor, HardRock'
      endif

C *** Cluster 04-High (NonRift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock ***
C     Model Number = 201315
      if (jcalc .eq. 201315) then
         call S06_EPRI13C4NRHigh ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-High-NonRift,MidC, Hor, HardRock'
      endif

C********************************************************************
C ***** EPRI Updated (2013) GMPE Models *****
C    Nomeclature is as follows:
C      Number 1-4 = "2013" EPRI Update 2013 Mid Continent Models
C      Number 5-6 = Median model cases with distance dependent sigma model
C                     "21" = Cluster01-Low,  Functional Model 1&3
C                     "22" = Cluster01-Med,  Functional Model 1&3
C                     "23" = Cluster01-High, Functional Model 1&3
C                     "24" = Cluster02-Low,  Functional Model 2
C                     "25" = Cluster02-Med,  Functional Model 2
C                     "26" = Cluster02-High, Functional Model 2
C                     "27" = Cluster03-Low,  Functional Model 1&3
C                     "28" = Cluster03-Med,  Functional Model 1&3
C                     "29" = Cluster03-High, Functional Model 1&3
C                     "30" = Cluster04-Low (Rift),  Functional Model 4
C                     "31" = Cluster04-Med (Rift),  Functional Model 4
C                     "32" = Cluster04-High (Rift), Functional Model 4
C                     "33" = Cluster04-Low (NonRift),  Functional Model 4
C                     "34" = Cluster04-Med (NonRift),  Functional Model 4
C                     "35" = Cluster04-High (NonRift), Functional Model 4
C********************************************************************
C

C *** Cluster 01-Low, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201321
      if (jcalc .eq. 201321) then
         call S06_EPRI13C1Low ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster01-Low,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 01-Med, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201322
      if (jcalc .eq. 201322) then
         call S06_EPRI13C1Med ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster01-Med,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 01-High, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201323
      if (jcalc .eq. 201323) then
         call S06_EPRI13C1High ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster01-High,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 02-Low, Mid-Continent: Functional Model 2, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201324
      if (jcalc .eq. 201324) then
         call S06_EPRI13C2Low ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster02-Low,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 02-Med, Mid-Continent: Functional Model 2, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201325
      if (jcalc .eq. 201325) then
         call S06_EPRI13C2Med ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster02-Med,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 02-High, Mid-Continent: Functional Model 2, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201326
      if (jcalc .eq. 201326) then
         call S06_EPRI13C2High ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster02-High,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 03-Low, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201327
      if (jcalc .eq. 201327) then
         call S06_EPRI13C3Low ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster03-Low,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 03-Med, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201328
      if (jcalc .eq. 201328) then
         call S06_EPRI13C3Med ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster03-Med,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 03-High, Mid-Continent: Functional Model 1&3, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201329
      if (jcalc .eq. 201329) then
         call S06_EPRI13C3High ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster03-High,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 04-Low (Rift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201330
      if (jcalc .eq. 201330) then
         call S06_EPRI13C4RLow ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Low-Rift,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 04-Med (Rift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201331
      if (jcalc .eq. 201331) then
         call S06_EPRI13C4RMed ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Med-Rift,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 04-High (Rift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201332
      if (jcalc .eq. 201332) then
         call S06_EPRI13C4RHigh ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-High-Rift,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif


C *** Cluster 04-Low (NonRift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201333
      if (jcalc .eq. 201333) then
         call S06_EPRI13C4NRLow ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Low-NonRift,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 04-Med (NonRift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201334
      if (jcalc .eq. 201334) then
         call S06_EPRI13C4NRMed ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-Med-NonRift,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

C *** Cluster 04-High (NonRift), Mid-Continent: Functional Model 4, Horizontal, CEUS Hard Rock, Rjb Sigma model ***
C     Model Number = 201335
      if (jcalc .eq. 201335) then
         call S06_EPRI13C4NRHigh ( mag, jbdist, lnY, specT,
     1                  attenName1, period2, iflag, sigma )
         attenname1 = 'EPRI(2013),Cluster04-High-NonRift,MidC, Hor, HardRock, Rjb Sigma model'
C     Adjust the Rjb distance sigma model
         if (jbdist .le. 10.0) then
            sigma = sqrt (sigma*sigma + 0.16*0.16)
         elseif (jbdist .lt. 20.0) then
            sjb = 0.16*(1.0 - alog(jbdist/10.0)/alog(20.0/10.0))
            sigma = sqrt(sigma*sigma + sjb*sjb)
         endif
      endif

c ******* NGA-East 2018 Models *******
c     Goulet et al., 2018 (PEER Report 2018/08) - Appendix H (17 median models)
c     Model Numbers 7001 - 7017
      if ( jcalc .ge. 7001 .and. jcalc .le. 7017 ) then
        imod = jcalc - 7000
        write(number, '(i0)') imod
        call S34_NGAEast_Med ( mag, rupDist, specT, imod, period2, lnY, iflag )
        attenname1 = 'NGAEast_2018_MedianModel_No'//trim(adjustl(number))
      endif

c     Goulet et al., 2018 (PEER Report 2018/08) - NGA-East Composite Ergodic Sigma, CENA, Low
c     Model Number = 7101
      if ( jcalc .eq. 7101 ) then
        call S32_NGAEast_CompErgSig_Low ( mag, specT, sigma, iflag )
        attenname1 = 'NGAEast_Composite_Ergodic_Sigma_CENA_Low'
c       keep median ground motions large since this is only for sigma model
        lnY = 1.0e10
      endif

c     Goulet et al., 2018 (PEER Report 2018/08) - NGA-East Composite Ergodic Sigma, CENA, Central
c     Model Number = 7102
      if ( jcalc .eq. 7102 ) then
        call S32_NGAEast_CompErgSig_Cen ( mag, specT, sigma, iflag )
        attenname1 = 'NGAEast_Composite_Ergodic_Sigma_CENA_Central'
c       keep median ground motions large since this is only for sigma model
        lnY = 1.0e10
      endif

c     Goulet et al., 2018 (PEER Report 2018/08) - NGA-East Composite Ergodic Sigma, CENA, High
c     Model Number = 7103
      if ( jcalc .eq. 7103 ) then
        call S32_NGAEast_CompErgSig_High ( mag, specT, sigma, iflag )
        attenname1 = 'NGAEast_Composite_Ergodic_Sigma_CENA_High'
c       keep median ground motions large since this is only for sigma model
        lnY = 1.0e10
      endif

c     Goulet et al., 2018 (PEER Report 2018/08) - NGA-East Composite Single-Station Sigma, CENA, Low
c     Model Number = 7104
      if ( jcalc .eq. 7104 ) then
        call S32_NGAEast_CompSSSig_Low ( mag, specT, sigma, iflag )
        attenname1 = 'NGAEast_Composite_SingleStation_Sigma_CENA_Low'
c       keep median ground motions large since this is only for sigma model
        lnY = 1.0e10
      endif

c     Goulet et al., 2018 (PEER Report 2018/08) - NGA-East Composite Single-Station Sigma, CENA, Central
c     Model Number = 7105
      if ( jcalc .eq. 7105 ) then
        call S32_NGAEast_CompSSSig_Cen ( mag, specT, sigma, iflag )
        attenname1 = 'NGAEast_Composite_SingleStation_Sigma_CENA_Central'
c       keep median ground motions large since this is only for sigma model
        lnY = 1.0e10
      endif

c     Goulet et al., 2018 (PEER Report 2018/08) - NGA-East Composite Single-Station Sigma, CENA, High
c     Model Number = 7106
      if ( jcalc .eq. 7106 ) then
        call S32_NGAEast_CompSSSig_High ( mag, specT, sigma, iflag )
        attenname1 = 'NGAEast_Composite_SingleStation_Sigma_CENA_High'
c       keep median ground motions large since this is only for sigma model
        lnY = 1.0e10
      endif


c     Check for valid jcalc
      if ( lnY .gt. 1.0e10 ) then
         write (*,'( 2x,''invalid jcalc:'',i7,3f10.4,e12.4,f12.4)') jcalc, mag,
     1      rupDist, ftype, lnY, sigma
         stop 99
      endif

      attenName(jType,iAtten) = attenname1
      period1(jType,iProb) = period2
      intflag(jType,iProb) = iflag
      siga = sigma

c     Convert to g
      lnY = lnY - 6.89

      return
      end
