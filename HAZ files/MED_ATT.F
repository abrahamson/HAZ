c meanpga.f

      subroutine meanInten (rupdist, jbdist, seismodist, hwflag,
     1  mag, jcalc, specT, lnY, siga, ftype, attenName, period1, iProb,
     2  jType, vs, depth,intflag, AR, dip1, disthypo, depthvs10, 
     3  depthvs15, D25, tau, depthTop, Theta_Site, RupWidth, vs30_class,
     4  foreArc, Rx, phi, cfcoefrrup, cfcoefrjb, Ry0 )        

      include 'pfrisk.h'
      
      real mag, lnY, siga, ftype, period1(4,1), mag1, period2, lnSa, m
      integer jcalc,soilflag,softrock,hardRock,hwflag
      real jbdist, rupdist, seismodist, baseDepth, depth, specT, factor, a, b, faddmag
      character*80 attenName(4,1), attenname1, attenname0
      integer intflag(4,1), iflag, vs30_class
      integer iflag01, iflag02, iflag04, iflag10, foreArc, regionflag, basinflag, msasflag
      integer coefcountRrup, coefcountRjb

      real sc, sd, se, D25, RupWidth, Ss, Sr, Q0, SCa, SCb, SCc, SCd, SCe
      real magc, rupdistc, ftypec, lnYc, mcb, depthTop
      real disthypo, depthvs10, depthvs15, AR, dip1(1), dip
      real F_RV, F_NM, Theta_Site, Z1
      real lnY01, lnY02, lnY04, lnY10, lnY02p, lnY04p
      real sigma01, sigma02, sigma04, sigma10
      real specT01, specT02, specT04, specT10, period02, Rx, SF2, Ry0
      real phi, tau, lnYH
      real cfcoefrrup(MAX_Atten,11), cfcoefrjb(MAX_Atten,11)
      
C LNY IS EXPECTED INTENSITY FOR THIS MAGNITUDE AND CLOSEST DISTANCE
      lnY = 1.e30
      tau = 0.0
      phi = 0.0
      iflag = 0
      efac = 0.0
       
c  *** Turkey adjusted NGA1 models
      if (jcalc .eq. 9501 ) then
          call AS_NGA_2008TR ( mag, dip, fType, RuoWidth, rupDist, jbdist,
     1                     vs, hwflag, lnY, sigma, 
     2                     specT, period2, depthtop, iflag, vs30_class,
     3                     depthvs10, Rx )     
      endif

      if (jcalc .eq. 9502 ) then
         call BA_NGA_2008TR ( mag, jbdist, specT, period2, lnY, sigma, iflag, vs, ftype, pga4nl ) 
      endif

      if (jcalc .eq. 9503 ) then
         call CB_NGA_2008TR ( mag, rupDist, jbdist, Ftype, specT, 
     1                     period2, lnY, sigma, iflag, vs,
     2                     depthtop, D25, Dip ) 
      endif

      if (jcalc .eq. 9504) then
          call CY_NGA_2008TR ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip, depthTop, Ftype, depthvs10, vs30_class,
     3                     hwflag, Rx )
      endif

      if (jcalc .eq. 9505 ) then
            call I_NGA_2008TR ( mag, rupDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
      endif
      if (jcalc .eq. 9506) then
          if ( depthTop .gt. 20. ) then
             depthTop1 = 20.
          else
             depthTop1 = depthTop
          endif
          call CY_NGA_2008TR ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip, depthTop1, Ftype, depthvs10, vs30_class,
     3                     hwflag, Rx )
      endif

     
C ******* PEER NGA Attenuation Models ****
c ******* Abrahamson and Silva Model *********
C     Abrahamson&Silva 2008 - horizontal, Estimated Vs30m
C     Model Number = 787
      if ( jcalc .eq. 787 ) then
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call AS_NGA_2008 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, phi, tau ) 
         attenname1 = 'A&S_NGA_2008-Hor,Estimated Vs30m'
         lnY = lnY + efac
      endif

C     Abrahamson&Silva 2008 - horizontal, Measured Vs30m
C     Model Number = 788
      if ( jcalc .eq. 788 ) then
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call AS_NGA_2008 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, phi, tau ) 
         attenname1 = 'A&S_NGA_2008-Hor,Measured Vs30m'
         lnY = lnY + efac
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2008 - Horizontal, estimated Vs30m
C     Model Number = 797
c
      if ( jcalc .eq. 797 ) then
         vs30_class = 0
         call CY_NGA_2008 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, phi, tau )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Estimated Vs30m'
         lnY = lnY + efac
       endif

C     Chiou and Youngs 2008 - Horizontal, measured Vs30m
C     Model Number = 798
c
      if ( jcalc .eq. 798 ) then
         vs30_class = 1
         call CY_NGA_2008 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, phi, tau )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Measured Vs30m'
         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model - Small Magnitude Models *********
C     Chiou and Youngs 2008 - Horizontal, estimated Vs30m
C     Southern California Small Magnitude Model (2010)
C     Model Number = 799
c
      if ( jcalc .eq. 799 ) then
         vs30_class = 0
         call CY_NGA_2008SC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Estimated Vs30m,SCal SMM'
         lnY = lnY + efac
       endif

C     Chiou and Youngs 2008 - Horizontal, measured Vs30m
C     Southern California Small Magnitude Model (2010)
C     Model Number = 800
c
      if ( jcalc .eq. 800 ) then
         vs30_class = 1
         call CY_NGA_2008SC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Measured Vs30m,SCal SMM'
         lnY = lnY + efac
       endif

C     Chiou and Youngs 2008 - Horizontal, estimated Vs30m
C     Central California Small Magnitude Model (2010)
C     Model Number = 801
c
      if ( jcalc .eq. 801 ) then
         vs30_class = 0
         call CY_NGA_2008CC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Estimated Vs30m,CCal SMM'
         lnY = lnY + efac
       endif

C     Chiou and Youngs 2008 - Horizontal, measured Vs30m
C     Central California Small Magnitude Model (2010)
C     Model Number = 802
c
      if ( jcalc .eq. 802 ) then
         vs30_class = 1
         call CY_NGA_2008CC ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx )
         attenname1 = 'Chiou&Youngs_NGA_2008-Hor,Measured Vs30m,CCal SMM'
         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2008 - horizontal
C     Model Number = 836
      if ( jcalc .eq. 836 ) then
         call CB_NGA_2008 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), phi, tau )
         attenname1 = 'Campbell&Bozorgnia_NGA_2008-Hor'
         lnY = lnY + efac
       endif

c ******* Idriss Model *********
C     Idriss 2008 - Horizontal, Vs=450-900 m/sec and Vs>900 m/sec
C     Model Number = 910
      if ( jcalc .eq. 910 ) then
         if (vs .ge. 900.0) then
            call I_NGA_2008vsgt900 ( mag, rupDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGA_2008_Hor,Vs>900m/s'
         else
            call I_NGA_2008 ( mag, rupDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGA_2008_Hor,Vs=450-900m/s'
         endif
         lnY = lnY + efac

       endif

c ******* Boore and Atkinson Model *********
C     Boore and Atkinson July 2008 - horizontal
C     Model Number = 922
      if ( jcalc .eq. 922 ) then
         call BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pga4nl, phi, tau )
         attenname1 = 'Boore&Atkinson_NGA_2008_Hor'
         lnY = lnY + efac
       endif

c ******* Boore and Atkinson Model *********
C     Boore and Atkinson July 2008 - horizontal with Atkinson (2010) small magnitude adjustment
C     Model Number = 923
      if ( jcalc .eq. 923 ) then
         call BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pga4nl, phi, tau )
C     Apply small magnitude adjustment is Mag<=5.75.
         if (mag .le. 5.75 ) then
            a = max(0.0, 3.888 - 0.674*mag)
            b = max(0.0, 2.933 - 0.510*mag)
            factor = a - b*log10(jbdist+10.0)
            factor = factor*alog(10.0)
         endif
         attenname1 = 'Boore&Atkinson_NGA_2008_Hor with small mag adj Atkinson (2010)'
         lnY = lnY + factor + efac
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
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-EstVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Estimated Vs30m
C     Model Number = 2788
      if ( jcalc .eq. 2788 ) then
         regionflag = 1
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-EstVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Estimated Vs30m
C     Model Number = 2789
      if ( jcalc .eq. 2789 ) then
         regionflag = 2
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-EstVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Estimated Vs30m
C     Model Number = 2790
      if ( jcalc .eq. 2790 ) then
         regionflag = 3
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-EstVs'
         lnY = lnY + efac
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Mainshock, Measured Vs30m
C     Model Number = 2791
      if ( jcalc .eq. 2791 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-MS-MesVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Mainshock, Measured Vs30m
C     Model Number = 2792
      if ( jcalc .eq. 2792 ) then
         regionflag = 1
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-MS-MesVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Mainshock, Measured Vs30m
C     Model Number = 2793
      if ( jcalc .eq. 2793 ) then
         regionflag = 2
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-MS-MesVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Mainshock, Measured Vs30m
C     Model Number = 2794
      if ( jcalc .eq. 2794 ) then
         regionflag = 3
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-MS-MesVs'
         lnY = lnY + efac
      endif

C     Note: These calls are for Aftershocks cases but it not currently implemented. 
C           For implementation new distance term CRjb will need to be computed and passed. 
c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Aftershock, Estimated Vs30m
C     Model Number = 3787
      if ( jcalc .eq. 3787 ) then
         regionflag = 0
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-AS-EstVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Aftershock, Estimated Vs30m
C     Model Number = 3788
      if ( jcalc .eq. 3788 ) then
         regionflag = 1
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-AS-EstVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Aftershock, Estimated Vs30m
C     Model Number = 3789
      if ( jcalc .eq. 3789 ) then
         regionflag = 2
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-AS-EstVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Aftershock, Estimated Vs30m
C     Model Number = 3790
      if ( jcalc .eq. 3790 ) then
         regionflag = 3
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-AS-EstVs'
         lnY = lnY + efac
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Global, Aftershock, Measured Vs30m
C     Model Number = 3791
      if ( jcalc .eq. 3791 ) then
         regionflag = 0
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Glob-AS-MesVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Taiwan, Aftershock, Measured Vs30m
C     Model Number = 3792
      if ( jcalc .eq. 3792 ) then
         regionflag = 1
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Taiw-AS-MesVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, China, Aftershock, Measured Vs30m
C     Model Number = 3793
      if ( jcalc .eq. 3793 ) then
         regionflag = 2
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-China-AS-MesVs'
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2013 - Horizontal, Japan, Aftershock, Mesaured Vs30m
C     Model Number = 3794
      if ( jcalc .eq. 3794 ) then
         regionflag = 3
         msasflag = 1
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau ) 
         attenname1 = 'ASK_NGAWest2_2013-Hor-Japan-AS-MesVs'
         lnY = lnY + efac
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal
C     Model Number = 2797
      if ( jcalc .eq. 2797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Estimated Vs30m'
         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal
C     Model Number = 2798
      if ( jcalc .eq. 2798 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor,Measured Vs30m'
         lnY = lnY + efac         
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 2799
      if ( jcalc .eq. 2799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 1
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Estimated Vs30m'
         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Japan and Italy Adjustment
C      Note: Only valid for the Following magnitude Range: 6<M<6.9
C     Model Number = 2800
      if ( jcalc .eq. 2800 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 1
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Jap/Ity,Measured Vs30m'
         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment
C      Note: Only valid for the Following magnitude : 6<M<6.9
C     Model Number = 2801
      if ( jcalc .eq. 2801 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 2
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchuan,Estimated Vs30m'
         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Horizontal, Wenchaun Adjustment
C      Note: Only valid for the Following magnitude: M=7.9
C     Model Number = 2802
      if ( jcalc .eq. 2802 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 2
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Hor-Wenchaun,Measured Vs30m'
         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, California
C     Model Number = 2836
      if ( jcalc .eq. 2836 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Calif'
         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Japan
C     Model Number = 2837
      if ( jcalc .eq. 2837 ) then
         regionflag = 1
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Japan'
         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, China
C     Model Number = 2838
      if ( jcalc .eq. 2838 ) then
         regionflag = 2
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,China'
         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2013 - horizontal, Italy
C     Model Number = 2839
      if ( jcalc .eq. 2839 ) then
         regionflag = 3
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2013-Hor,Italy'
         lnY = lnY + efac
       endif


c ******* Idriss Model *********
C     Idriss 2013 - Horizontal
C     Model Number = 2910
      if ( jcalc .eq. 2910 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2013_Hor'
         else
            write (*,*) 'Idriss NGA West 2 GMPE not defined'
            write (*,*) 'for Vs<450m/s.'
            stop 99
         endif
         lnY = lnY + efac
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 2922
      if ( jcalc .eq. 2922 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, No Basin'
         lnY = lnY + efac
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, No Basin Adjustment
C     Model Number = 2923
      if ( jcalc .eq. 2923 ) then
         regionflag = 1
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, No Basin'
         lnY = lnY + efac
       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, No Basin Adjustment
C     Model Number = 2924
      if ( jcalc .eq. 2924 ) then
         regionflag = 2
         Basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag,
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, No Basin'
         lnY = lnY + efac
       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Global Adjustments, Basin Adjustments
C     Model Number = 2925
      if ( jcalc .eq. 2925 ) then
         regionflag = 0
         basinflag = 1
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, 
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3Global, Basin'
         lnY = lnY + efac
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model China-Turkey Adjustments, Basin Adjustment
C     Model Number = 2926
      if ( jcalc .eq. 2926 ) then
         regionflag = 1
         basinflag = 1
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, 
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ChinaTurkey, Basin'
         lnY = lnY + efac
       endif
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2013 - horizontal
C          DeltaC3 Model Italy-Japan Adjustments, Basin Adjustment
C     Model Number = 2927
      if ( jcalc .eq. 2927 ) then
         regionflag = 2
         Basinflag = 1
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, Basinflag, 
     1               phi, tau  )
         attenname1 = 'BSSA_NGAWest2_2013_Hor, DC3ItalyJapan, Basin'
         lnY = lnY + efac
       endif



C **** NGA West2 Vertical Models *********
c ******* Stewart, Seyhan, Boore, and Atkinson Model *********
C     Stewart, Seyhan, Boore, and Atkinson 2013 - vertical
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 4922
      if ( jcalc .eq. 4922 ) then
         regionflag = 0
         basinflag = 0
         call SSBA_NGAWest2_2013_Vert ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, phi, tau )
         attenname1 = 'SSBA_NGAWest2_2013_Ver, DC3Global, No Basin'
         lnY = lnY + efac
       endif

C     Stewart, Seyhan, Boore, and Atkinson 2013 - V/H Ratio
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Model Number = 6922
      if ( jcalc .eq. 6922 ) then
         regionflag = 0
         basinflag = 0
         call SSBA_NGAWest2_2013_Vert ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, phi, tau )

         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnYH, sigmaH, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag, phiH, tauH )
         attenname1 = 'SSBA_NGAWest2_2013_V/H, DC3Global, No Basin'
C     Divide Vertical by horizontal to get V/H ratio
         lnY = lnY - lnYH + 6.89

         lnY = lnY + efac
       endif

c ******* Gulerce, Abrahamson, Silva, and Kamai Model *********
C     Gulerce, Abrahamson, Silva, and Kamai 2013 - Vertical, Global, Mainshock, Estimated Vs30m

C     Model Number = 4787
      if ( jcalc .eq. 4787 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call GKAS_NGAWest2_2013_Vert ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
         attenname1 = 'GKAS_NGAWest2_2013-Ver-Glob-MS-EstVs'
         lnY = lnY + efac
      endif

C     Gulerce, Abrahamson, Silva, and Kamai 2013 - Vertical, Global, Mainshock, Measured Vs30m

C     Model Number = 4788
      if ( jcalc .eq. 4788 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call GKAS_NGAWest2_2013_Vert ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
         attenname1 = 'GKAS_NGAWest2_2013-Ver-Glob-MS-MsrVs'
         lnY = lnY + efac
      endif

c ******* Gulerce, Abrahamson, Silva, and Kamai Model *********
C     Gulerce, Abrahamson, Silva, and Kamai 2013 - V/H Ratio, Global, Mainshock, Estimated Vs30m

C     Model Number = 6787
      if ( jcalc .eq. 6787 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call GKAS_NGAWest2_2013_Vert ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau)      
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnYH, sigmaH, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phiH, tauH) 

C     Divide Vertical by Horizontal to get V/H Ratio
         lnY = lnY - lnYH + 6.89

         attenname1 = 'GKAS_NGAWest2_2013-V/H-Glob-MS-EstVs'
         lnY = lnY + efac
      endif

c ******* Gulerce, Abrahamson, Silva, and Kamai Model *********
C     Gulerce, Abrahamson, Silva, and Kamai 2013 - V/H Ratio, Global, Mainshock, Measured Vs30m

C     Model Number = 6787
      if ( jcalc .eq. 6787 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 1
         call GKAS_NGAWest2_2013_Vert ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
     
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnYH, sigmaH, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phiH, tauH) 

C     Divide Vertical by Horizontal to get V/H Ratio
         lnY = lnY - lnYH + 6.89

         attenname1 = 'GKAS_NGAWest2_2013-V/H-Glob-MS-MsrVs'
         lnY = lnY + efac
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2013 - Vertical 
C     Model Number = 4797
      if ( jcalc .eq. 4797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Ver,Estimated Vs30m'
         lnY = lnY + efac
       endif

C     Chiou and Youngs 2013 - Vertical 
C     Model Number = 4798
      if ( jcalc .eq. 4798 ) then
c     Current model set for measured Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
         attenname1 = 'CY_NGAWest2_2013-Ver,Measured Vs30m'
         lnY = lnY + efac
       endif

C     V/H Model
C     Model Number = 6797
      if ( jcalc .eq. 6797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
C     Call the horizontal model
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnYH, sigmaH, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phiH, tauH )

         attenname1 = 'CY_NGAWest2_2013-V/H,Estimated Vs30m'

C     Divide vertical by horizontal to get V/H Ratio. 
         lnY = lnY - lnYH + 6.89    

         lnY = lnY  + efac
       endif

C     V/H Model
C     Model Number = 6798
      if ( jcalc .eq. 6798 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 1
         regionflag = 0
         call CY_NGAWest2_2013_V ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phi, tau )
C     Call the horizontal model
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnYH, sigmaH, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag, phiH, tauH )

         attenname1 = 'CY_NGAWest2_2013-V/H,Measured Vs30m'

C     Divide vertical by horizontal to get V/H Ratio. 
         lnY = lnY - lnYH + 6.89    

         lnY = lnY  + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Bozorgnia and Campbell 2013 - Vertical, California
C     Model Number = 4836
      if ( jcalc .eq. 4836 ) then
         regionflag = 0
         call BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )
         attenname1 = 'BC_NGAWest2_2013-Ver,Calif'
         lnY = lnY + efac
       endif

C     Bozorgnia and Campbell 2013 - Vertical, Japan
C     Model Number = 4837
      if ( jcalc .eq. 4837 ) then
         regionflag = 1
         call BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )
         attenname1 = 'BC_NGAWest2_2013-Ver,Japan'
         lnY = lnY + efac
       endif

C     Bozorgnia and Campbell 2013 - Implied Vertical/Horizontal Ratio, California
C     Model Number = 6836
      if ( jcalc .eq. 6836 ) then
         regionflag = 0
         call BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )
     
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnYH, sigmaH, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag, phiH, tauH )

C     Divide vertical by horizontal to get V/H Ratio. 
         lnY = lnY - lnYH + 6.89    

         attenname1 = 'BC_NGAWest2_2013-V/H,Calif'
         lnY = lnY + efac
       endif
C     Bozorgnia and Campbell 2013 - Implied Vertical/Horizontal Ratio, California
C     Model Number = 6837
      if ( jcalc .eq. 6837 ) then
         regionflag = 1
         call BC_NGAWest2_2013_Vert ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag, phi, tau )
     
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnYH, sigmaH, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag, phiH, tauH )

C     Divide vertical by horizontal to get V/H Ratio. 
         lnY = lnY - lnYH + 6.89    

         attenname1 = 'BC_NGAWest2_2013-V/H,Japan'
         lnY = lnY + efac
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

         call Camp03_H ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Firm Soil'
         lnY = lnY + efac
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

         call Camp03_H ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Very Firm Soil'
         lnY = lnY + efac
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

         call Camp03_H ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Soft Rock'
         lnY = lnY + efac
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

         call Camp03_H ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,Firm Rock'
         lnY = lnY + efac
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

         call Camp03_H ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,GenericRock'
         lnY = lnY + efac
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

         call Camp03_H ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Hor,GenericSoil'
         lnY = lnY + efac
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

         call Camp03_V ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Firm Soil'
         lnY = lnY + efac
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

         call Camp03_V ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Very Firm Soil'
         lnY = lnY + efac
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

         call Camp03_V ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Soft Rock'
         lnY = lnY + efac
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

         call Camp03_V ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,Firm Rock'
         lnY = lnY + efac
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

         call Camp03_V ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,GenericRock'
         lnY = lnY + efac
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

         call Camp03_V ( mag, seismoDist, jbDist, lnY, sigma, 
     1        specT, period2, s03vfs, s03sr, s03fr, frv, fth, hwflag, 
     2        iflag)
         attenname1 = 'Campbell&Bozorgnia(2003)-Ver,GenericSoil'
         lnY = lnY + efac
      endif


c ******* Abrahamson and Silva Models *********
C     Abrahamson&Silva 1997 (Rock) - horizontal
C     Model Number = 001
      if ( jcalc .eq. 1 ) then
         soil = 0.
         call AS_97_H ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Rock) - vertical
C     Model Number = 002
      if ( jcalc .eq. 2 ) then
         soil = 0.
         call AS_97_V ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,rock'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Soil) - horizontal
C     Model Number = 003
      if ( jcalc .eq. 3 ) then
         soil = 1.
         call AS_97_H ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,soil'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1995 (Soil) - vertical
C     Model Number = 004
      if ( jcalc .eq. 4 ) then
         soil = 1.
         call AS_97_V ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,soil'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors.
C     Model Number = 005
      if ( jcalc .eq. 5 ) then
         soil = 0.
         call AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock,NF'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Rock) - vertical
C     Model Number = 006
      if ( jcalc .eq. 6 ) then
         soil = 0.
         call AS_97_V_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,rock,NF'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Soil) - horizontal with Normal faulting factors.
C     Model Number = 007
      if ( jcalc .eq. 7 ) then
         soil = 1.
         call AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,soil,NF'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Soil) - vertical with Normal faulting factors
C     Model Number = 008
      if ( jcalc .eq. 8 ) then
         soil = 1.
         call AS_97_V_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Ver,soil,NF'
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors scaled
C           by 1.0/1.67 factor.
C     Model Number = 2005
      if ( jcalc .eq. 2005 ) then
         soil = 0.
         call AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock,NF(1/1.67)'
C      Scale the computed ground motion by 1.0/1.67
         lnY = lnY + alog(1.0/1.67)
         lnY = lnY + efac
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors scaled
C           by 1.67 factor.
C     Model Number = 3005
      if ( jcalc .eq. 3005 ) then
         soil = 0.
         call AS_97_H_NF ( mag, rupDist, ftype, soil, hwflag, specT,
     1                     period2, lnY, sigma,iflag )
         attenname1 = 'Abrahamson/Silva(97)-Hor,rock,NF(1.67)'
C      Scale the computed ground motion by 1.67
         lnY = lnY + alog(1.67)
         lnY = lnY + efac
       endif


c ********* Boore, Joyner and Fumal Models ******************
c     BJF94, Horizontal, Class A
C     Model Number = 010
      if ( jcalc .eq. 10 ) then
         GB = 0
         GC = 0
         call bjf94 ( mag, jbDist, ftype, lnY, sigma, GB, GC, specT,
     1            attenName1, period2,iflag )
         lnY = lnY + efac
c     BJF94, Horizontal, Class B
C     Model Number = 011
      elseif ( jcalc .eq. 11 ) then
         GB = 1
         GC = 0
         call bjf94 ( mag, jbDist, ftype, lnY, sigma, GB, GC, specT,
     1            attenName1, period2,iflag )
         lnY = lnY + efac
c     BJF94, Horizontal, Class c
C     Model Number = 012
      elseif ( jcalc .eq. 12 ) then
         GB = 0
         GC = 1
         call bjf94 ( mag, jbDist, ftype, lnY, sigma, GB, GC, specT,
     1            attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     BJF97, Horizontal, Vs top 30 meters
C     Model Number = 013
      if ( jcalc .eq. 13 ) then
         call bjf97 ( mag, jbDist, ftype, lnY, sigma, specT,
     1            attenName1, period2, vs,iflag )
         lnY = lnY + efac
      endif

c     BJF97, Horizontal, Vs top 30 meters scaled by factor 1.0/1.67
C     Model Number = 2013
      if ( jcalc .eq. 2013 ) then
         call bjf97 ( mag, jbDist, ftype, lnY, sigma, specT,
     1            attenName1, period2, vs,iflag )
c     Scale ground motion by factor 1.0/1.67
         lnY = lnY + alog(1.0/1.67)
         attenName1='Boore, Joyner, Fumal (1997) (1.0/1.67)'                               
         lnY = lnY + efac
      endif

c     BJF97, Horizontal, Vs top 30 meters scaled by factor 1.67
C     Model Number = 3013
      if ( jcalc .eq. 3013 ) then
         call bjf97 ( mag, jbDist, ftype, lnY, sigma, specT,
     1            attenName1, period2, vs,iflag )
c     Scale ground motion by factor 1.67
         lnY = lnY + alog(1.67)
         attenName1='Boore, Joyner, Fumal (1997) (1.67)'                               
         lnY = lnY + efac
      endif

c ******** Campbell Models ******
c     Campbell (1990), Horizontal, Rock
C     Model Number = 020
      if ( jcalc .eq. 20 ) then
         baseDepth = 2.0
         call Camp90 ( mag, seismoDist, ftype, lnY, sigma, baseDepth,
     1            specT, attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Campbell (1990) - vertical, Rock
C     Model Number = 021
      if ( jcalc .eq. 21 ) then
         baseDepth = 2.0
         call Camp90v ( mag, seismoDist, ftype, lnY, sigma, baseDepth,
     1            specT, attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Campbell (1990/1994), Horizontal, Rock
C     Model Number = 022
      if ( jcalc .eq. 22 ) then
         baseDepth = 2.0
         soilflag = 0
         call Camp90_94 ( mag, seismoDist, ftype, lnY, sigma, baseDepth,
     1             specT, attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Campbell (1993-1994), Horizontal Soil
C     Model Number = 023
      if ( jcalc .eq. 23 ) then
         baseDepth = 2.0
         soilFlag = 0
         softRock = 0
         hardRock = 0
         call Campbell_94 ( mag, seismoDist, ftype, lnY, sigma, specT,
     1                soilFlag, softRock, hardRock, baseDepth,
     1                attenName1, period2,iflag )
         lnY = lnY + efac
c     Campbell (1993-1994), Horizontal Soft Rock
C     Model Number = 024
      elseif ( jcalc .eq. 24 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 1
         hardRock = 0
         call Campbell_94 ( mag, seismoDist, ftype, lnY, sigma, specT,
     1                soilFlag, softRock, hardRock, baseDepth,
     1                attenName1, period2,iflag )
         lnY = lnY + efac
c     Campbell (1993-1994), Horizontal Hard Rock
C     Model Number = 025
      elseif ( jcalc .eq. 25 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 0
         hardRock = 1
         call Campbell_94 ( mag, seismoDist, ftype, lnY, sigma, specT,
     1                soilFlag, softRock, hardRock, baseDepth,
     1                attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Campbell (1997), Horizontal, Soil
C     Model Number = 026
      if ( jcalc .eq. 26 ) then
         baseDepth = 4.0
         soilFlag = 0
         softRock = 0
         hardRock = 0
         call Camp97_H ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
         lnY = lnY + efac
c     Campbell (1997), Horizontal, Soft Rock
C     Model Number = 027
      elseif ( jcalc .eq. 27 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 1
         hardRock = 0
         call Camp97_H ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
         lnY = lnY + efac
c     Campbell (1997), Horizontal, Hard Rock
C     Model Number = 028
      elseif ( jcalc .eq. 28 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 0
         hardRock = 1
         call Camp97_H ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,                  
     1                attenName1, period2, softrock, hardRock,iflag )
         lnY = lnY + efac
      endif

c     Campbell (1997)  vertical, Soil
C     Model Number = 029
      if ( jcalc .eq. 29 ) then
         baseDepth = 4.0
         soilFlag = 0
         softRock = 0
         hardRock = 0
         call Camp97_Z ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,                  
     1                attenName1, period2, softrock, hardRock,iflag )
         lnY = lnY + efac
c     Campbell (1997)  vertical, Soft Rock
C     Model Number = 030
      elseif ( jcalc .eq. 30 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 1
         hardRock = 0
         call Camp97_Z ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,
     1                attenName1, period2, softrock, hardRock,iflag )
         lnY = lnY + efac
c     Campbell (1997)  vertical, Hard Rock
C     Model Number = 031
      elseif ( jcalc .eq. 31 ) then
         baseDepth = 2.0
         soilFlag = 1
         softRock = 0
         hardRock = 1
         call Camp97_Z ( mag, seismoDist, ftype, lnY, sigma,
     1                baseDepth, specT,                  
     1                attenName1, period2, softrock, hardRock,iflag )
         lnY = lnY + efac
      endif

c ******** Idriss Models *******
c     Idriss (1991), Horizontal, Rock
C     Model Number = 040
      if ( jcalc .eq. 40 ) then
         call Idriss91_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
         lnY = lnY + efac
      endif

C     Idriss (1991), Horizontal, Soft-soil, PGA
C     Model Number = 041
c     soft-soil
      if ( jcalc .eq. 41 .and. specT .eq. 0.0 ) then
         call Idriss91_soft (mag, rupDist, ftype, lnY, sigma,
     1            attenName1, period2 )
         iflag = 0
         lnY = lnY + efac
      elseif (jcalc .eq. 41 .and. specT .ne. 0.0) then
        write (*,*) 'Idriss (1991), Horizontal, Soft Soil'
        write (*,*) 'not defined for spectral acclereation!!!'
        write (*,*) 'Check input file.'
        stop 99
      endif

c     Idriss 1997 Horizontal, soft-soil, PGA
C     Model Number = 042
      if ( jcalc .eq. 42 .and. specT .eq. 0.0 ) then
         call Idriss97_soft (mag, rupDist, ftype, lnY, sigma,
     1            attenName1, period2 )
         iflag = 0
         lnY = lnY + efac
      elseif (jcalc .eq. 42 .and. specT .ne. 0.0) then
         write (*,*) 'Idriss (1997), Horizontal, Soft Soil'
         write (*,*) 'not defined for spectral acclereation!!!'
         write (*,*) 'Check input file.'
         stop 99
      endif

c     Idriss (1991:1995), Horizontal, Rock
C     Model Number = 043
      if ( jcalc .eq. 43 ) then
         call Idriss91_95_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Idriss (1991:1995), Horizontal, Rock scaled by factor 1.0/1.67
C     Model Number = 2043
      if ( jcalc .eq. 2043 ) then
         call Idriss91_95_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
c     Scaled ground motion by factor 1.0/1.67
         lnY = lnY + alog (1.0/1.67)
         attenName1 = 'Idriss (1991;1995), Rock (1.0/1.67)'                                              
         lnY = lnY + efac
      endif

c     Idriss (1991:1995), Horizontal, Rock scaled by factor 1.67
C     Model Number = 3043
      if ( jcalc .eq. 3043 ) then
         call Idriss91_95_rock (mag, rupDist, ftype, lnY, sigma, specT,
     1            attenName1, period2,iflag )
c     Scaled ground motion by factor 1.67
         lnY = lnY + alog (1.67)
         attenName1 = 'Idriss (1991;1995), Rock (1.67)'                                              
         lnY = lnY + efac
      endif

c  ******* Sadigh/Geomatrix Models *******
c     Geomatrix 93 (rock) vertical
C     Model Number = 050
      if ( jcalc .eq. 50 ) then
         call Geomatrix93_V_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Sadigh et al. 97 (rock) Horizontal
C     Model Number = 051
      if ( jcalc .eq. 51 ) then
         call Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
         lnY = lnY + efac         
      endif

c     Sadigh et al. 97 (rock) Horizontal scaled by factor 1.0/1.67
C     Model Number = 2051
      if ( jcalc .eq. 2051 ) then
         call Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
c     Scale ground motion by factor 1.0/1.67
         lnY = lnY + alog(1.0/1.67)
         attenName1 = 'Sadigh et al. (1997), Horizontal, rock (1.0/1.67)'                              
         lnY = lnY + efac
      endif

c     Sadigh et al. 97 (rock) Horizontal scaled by factor 1.67
C     Model Number = 3051
      if ( jcalc .eq. 3051 ) then
         call Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
c     Scale ground motion by factor 1.67
         lnY = lnY + alog(1.67)
         attenName1 = 'Sadigh et al. (1997), Horizontal, rock (1.67)'                              
         lnY = lnY + efac
      endif

c     Sadigh et al. 97 (soil) horizontal
C     Model Number = 052
      if ( jcalc .eq. 52 ) then
         call Sadigh97_H_soil ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
         lnY = lnY + efac
      endif

c     Sadigh et al. 97 (rock) Horizontal - Sigma = 0.0
C     Model Number = 053
      if ( jcalc .eq. 53 ) then
         call Geomatrix93_H_rock ( mag, rupDist, ftype, lnY, sigma,
     1             specT, attenName1, period2,iflag )
c     Now set sigma = 0.0
         sigma = 1.0e-10
         attenName1 = 'Sadigh et al. (1997),Hor.,rock,Sigma=0.0'                              
         lnY = lnY + efac
      endif

c ******** Spudich et al. (1997) Models *******
C     Spudich et al. (1997), Horizontal, Rock, Extensional Regimes
C     Model Number = 060
      if (jcalc .eq. 60) then
         call Spudich96 ( mag, JBdist, lnY, sigma, 0, specT,           
     1                   attenName1, period2,iflag )  
         attenname1 = 'Spudich et al. (1997), Horizontal, Rock'
         lnY = lnY + efac
      endif

C     Spudich et al. (1997), Horizontal, Soil, Extensional Regimes
C     Model Number = 061
      if (jcalc .eq. 61) then
         call Spudich96 ( mag, JBdist, lnY, sigma, 1, specT,           
     1                   attenName1, period2,iflag )  
         attenname1 = 'Spudich et al. (1997), Horizontal, Soil'
         lnY = lnY + efac
      endif

c ******** Youngs Models *******
c     Youngs et al (1993) Horizontal, subduction, Rock
C     Model Number = 200
      if ( jcalc .eq. 200 ) then
         call youngs93 ( mag, rupDist, lnY, sigma, attenName1,
     1       period2, specT, ftype,iflag )
         attenName1 ='Youngs (1993) Rock, Subduction'                                     
         lnY = lnY + efac
      endif

c     Youngs et al (1997) Horizontal, subduction, rock
C     Model Number = 201
      if ( jcalc .eq. 201 ) then
         call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1       period2, specT, ftype, depth,iflag )
         attenName1 ='Youngs (1997) Rock, Subduction'                                     
         lnY = lnY + efac
      endif

c     Youngs et al (1997) Horizontal, subduction, soil
C     Model Number = 202
      if ( jcalc .eq. 202 ) then
         call youngs97_soil ( mag, rupDist, lnY, sigma, attenName1,
     1       period2, specT, ftype, depth,iflag )
         attenName1 ='Youngs (1997) Soil, Subduction'                                     
         lnY = lnY + efac
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
            call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            soil = 0.
            hwflag = 1 
            call AS_97_H ( magc, rupDistc, ftypec, soil, hwflag, specT,
     1                     period2, lnYc, sigmac,iflag )
            attenname1 = 'Syn. Rupture: Youngs+AS (M7.4)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Idriss (1991:1995) Horizontal Crustal
C     Model Number = 7201043
         elseif ( jcalc .eq. 7201043 ) then
            call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call Idriss91_95_rock (magc, rupDistc, ftypec, lnYc, sigmac, specT,
     1            attenName1, period2,iflag )
            attenname1 = 'Synchronous Rupture: Youngs+Idriss (M7.4)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Sadigh et al. (1997) Horizontal Crustal
C     Model Number = 7201051
         elseif ( jcalc .eq. 7201051 ) then
            call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call Geomatrix93_H_rock ( magc, rupDistc, ftypec, lnYc, sigmac,
     1             specT, attenName1, period2,iflag )
            attenname1 = 'Syn. Rupture: Youngs+Sadigh (M7.4)-Rock'
         endif

c     First Combine the Sigma values.
         sigma=sqrt(((sigma**2)*(exp(lnY-6.89))**4 + 
     1          (sigmac**2)*(exp(lnYc-6.89))**4)
     2          / ( ((exp(lnY-6.89))**2 + (exp(lnYc-6.89))**2))**2 )

c     Now combine the Subduction GM and Crustal GM values.
         lnY = alog ( sqrt( (exp(lnY))**2 + (exp(lnYc))**2 ) )
         lnY = lnY + efac
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
            call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            soil = 0.
            hwflag = 1 
            call AS_97_H ( magc, rupDistc, ftypec, soil, hwflag, specT,
     1                     period2, lnYc, sigmac,iflag )
            attenname1 = 'Syn. Rupture: Youngs+AS (M7.7)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Idriss (1991:1995) Horizontal Crustal
C     Model Number = 8201043
         elseif ( jcalc .eq. 8201043 ) then
            call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call Idriss91_95_rock (magc, rupDistc, ftypec, lnYc, sigmac, specT,
     1            attenName1, period2,iflag )
            attenname1 = 'Synchronous Rupture: Youngs+Idriss (M7.7)-Rock'
c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Sadigh et al. (1997) Horizontal Crustal
C     Model Number = 8201051
         elseif ( jcalc .eq. 8201051 ) then
            call youngs97_rock ( mag, rupDist, lnY, sigma, attenName1,
     1          period2, specT, ftype, depth,iflag )
            call Geomatrix93_H_rock ( magc, rupDistc, ftypec, lnYc, sigmac,
     1             specT, attenName1, period2,iflag )
            attenname1 = 'Syn. Rupture: Youngs+Sadigh (M7.7)-Rock'
         endif

c     First Combine the Sigma values.
         sigma=sqrt(((sigma**2)*(exp(lnY-6.89))**4 + 
     1          (sigmac**2)*(exp(lnYc-6.89))**4)
     2          / ( ((exp(lnY-6.89))**2 + (exp(lnYc-6.89))**2))**2 )

c     Now combine the Subduction GM and Crustal GM values.
         lnY = alog ( sqrt( (exp(lnY))**2 + (exp(lnYc))**2 ) )
         lnY = lnY + efac
      endif
           
c ***** Atkinson and Boore Subduction Models *****
c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 210
      if (jcalc .eq. 210) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-B'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 211
      if (jcalc .eq. 211) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-C'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 212
      if (jcalc .eq. 212) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-D'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 213
      if (jcalc .eq. 213) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Subduction, NEHRP-E'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 220
      if (jcalc .eq. 220) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-B'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 221
      if (jcalc .eq. 221) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-C'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 222
      if (jcalc .eq. 222) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-D'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 223
      if (jcalc .eq. 223) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Cascadia, NEHRP-E'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 230
      if (jcalc .eq. 230) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-B'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 231
      if (jcalc .eq. 231) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-C'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 232
      if (jcalc .eq. 232) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-D'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 233
      if (jcalc .eq. 233) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003, Sub-Japan, NEHRP-E'
         lnY = lnY + efac
      endif

C Atkinson and Boore (2008) Subduction Erratum Corrected for Interface Events ****
c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-B, Subduction
C     Model Number = 310
      if (jcalc .eq. 310) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-B'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the 
c         interpolation for periods between these end points. 
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call AB03 ( mag, rupdist, lnY10, sigma10, specT10,
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
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-C, Subduction
C     Model Number = 311
      if (jcalc .eq. 311) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-C'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the 
c         interpolation for periods between these end points. 
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call AB03 ( mag, rupdist, lnY10, sigma10, specT10,
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
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-D, Subduction
C     Model Number = 312
      if (jcalc .eq. 312) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-D'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the 
c         interpolation for periods between these end points. 
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call AB03 ( mag, rupdist, lnY10, sigma10, specT10,
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
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-E, Subduction
C     Model Number = 313
      if (jcalc .eq. 313) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call AB03 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Subduction, NEHRP-E'

c     Erratum correction only for interface events.
C     Need to check ground motions between 0.1-1.0sec because of the 
c         interpolation for periods between these end points. 
         if (ftype .eq. 0) then
            if (specT .gt. 0.1 .and. specT .lt. 1.0) then
c     Call the attenuation model with needed spectral periods
               specT01 = 0.1
               call AB03 ( mag, rupdist, lnY01, sigma01, specT01,
     1              attenName0, period02,iflag01, ftype, depth, Sc, Sd, Se)
               specT02 = 0.2
               call AB03 ( mag, rupdist, lnY02, sigma02, specT02,
     1              attenName0, period02,iflag02, ftype, depth, Sc, Sd, Se)
               specT04 = 0.4
               call AB03 ( mag, rupdist, lnY04, sigma04, specT04,
     1              attenName0, period02,iflag04, ftype, depth, Sc, Sd, Se)
               specT10 = 1.0
               call AB03 ( mag, rupdist, lnY10, sigma10, specT10,
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
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-B, Subduction
C     Model Number = 320
      if (jcalc .eq. 320) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-B'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-C, Subduction
C     Model Number = 321
      if (jcalc .eq. 321) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-C'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-D, Subduction
C     Model Number = 322
      if (jcalc .eq. 322) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-D'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-E, Subduction
C     Model Number = 323
      if (jcalc .eq. 323) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call AB03Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Cascadia, NEHRP-E'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-B, Subduction
C     Model Number = 330
      if (jcalc .eq. 330) then
C        Set NEHRP-B site class by setting Sc=Sd=Se=0
         Sc = 0
         Sd = 0
         Se = 0
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-B'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-C, Subduction
C     Model Number = 331
      if (jcalc .eq. 331) then
C        Set NEHRP-C site class by setting Sc=1, Sd=Se=0
         Sc = 1
         Sd = 0
         Se = 0
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-C'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-D, Subduction
C     Model Number = 332
      if (jcalc .eq. 332) then
C        Set NEHRP-D site class by setting Sc=Se=0, Sd=1
         Sc = 0
         Sd = 1
         Se = 0
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-D'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
      endif

c     Atkinson and Boore (2003/08) - Horizontal, NEHRP-E, Subduction
C     Model Number = 333
      if (jcalc .eq. 333) then
C        Set NEHRP-E site class by setting Sc=Sd=0, Se=1
         Sc = 0
         Sd = 0
         Se = 1
         call AB03Jap ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag, ftype, depth, Sc, Sd, Se)
         attenname1 = 'Atkinson&Boore 2003/08, Sub-Japan, NEHRP-E'

c     Erratum correction only for world-wide model case interface events.
c     Note that this model is the same as the 2003 model.

         lnY = lnY + efac
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
         call Gregor02CasR ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag)
         attenname1 = 'Gregor et al.(2002),Rock,Cascadia Subduction'
         lnY = lnY + efac
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
         call Gregor02CasS ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag)
         attenname1 = 'Gregor et al.(2002),Soil,Cascadia Subduction'
         lnY = lnY + efac
      endif

c ***** Gregor et al. Cascadia Subduction Model *****
c     Gregor et al. (2006) - Horizontal, Cascadia Subduction
C     Model Number = 242

      if (jcalc .eq. 242) then
         call Gregor06Cas ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, vs, period2, iflag)
         attenname1 = 'Gregor et al.(2006), Cascadia Subduction'
         lnY = lnY + efac
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

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Sub., Hard Rock'
         lnY = lnY + efac
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
         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Sub., Rock (SC I)'
         lnY = lnY + efac
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

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Sub., Hard Soil (SC II)'
         lnY = lnY + efac
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

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Sub., Medium Soil (SC III)'
         lnY = lnY + efac
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

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Sub., Soft Soil (SC IV)'
         lnY = lnY + efac
      endif

c     Zhao et al. (2006) - Horizontal, Hard Rock, Crustal
C     Model Number = 255
      if (jcalc .eq. 255) then
c        Set Site Class = Hard Rock and Source Type
         sclass = 0
         sourceclass = 0

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Crust, Hard Rock'
         lnY = lnY + efac
      endif

c     Zhao et al. (2006) - Horizontal, Rock, Crustal
C     Model Number = 256
      if (jcalc .eq. 256) then
c        Set Site Class = Rock and Source Type
         sclass = 1.0
         sourceclass = 0

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Crust, Rock (SC I)'
         lnY = lnY + efac
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
         call Zhaoetal2006 ( m, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
         attenname1 = 'ZhaoLu(2011)-Crust, Rock (SC I)'
         lnY = lnY + efac
      endif



c     Zhao et al. (2006) - Horizontal, Hard Soil, Crustal
C     Model Number = 257
      if (jcalc .eq. 257) then
c        Set Site Class = Hard Soil and Source Type
         sclass = 2
         sourceclass = 0

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Crust, Hard Soil (SC II)'
         lnY = lnY + efac
      endif

c     Zhao et al. (2006) - Horizontal, Medium Soil, Crustal
C     Model Number = 258
      if (jcalc .eq. 258) then
c        Set Site Class = Medium Soil and Source Type
         sclass = 3
         sourceclass = 0

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Crust, Medium Soil (SC III)'
         lnY = lnY + efac
      endif

c     Zhao et al. (2006) - Horizontal, Soft Soil, Crustal
C     Model Number = 259
      if (jcalc .eq. 259) then
c        Set Site Class = Soft Soil and Source Type
         sclass = 4
         sourceclass = 0

         call Zhaoetal2006 ( mag, rupdist, ftype, lnY, sigma, sclass, specT,
     1            attenName1, period2,iflag, sourceclass, depth, phi, tau )
     
         attenname1 = 'Zhao etal(2006)-Crust, Soft Soil (SC IV)'
         lnY = lnY + efac
      endif

c ***** Kanno Subduction Models *****
c     Kanno (2006) - Horizontal, Subduction
C     Model Number = 260
      if ( jcalc .eq. 260 ) then

         call kanno2006 ( mag, rupdist, specT, 
     1                    period2, lnY, sigma, iflag, vs, depth ) 

         attenname1 = 'Kanno Subduction (2006)'
         lnY = lnY + efac
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
            call GarciaH05 ( mag, rupdist, specT, 
     1                    period2, lnY, sigma, iflag, depth ) 
            attenname1 = 'Garcia et al. (2005), Hor-Inslab, Rock'
         endif
         lnY = lnY + efac
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
            call GarciaV05 ( mag, rupdist, specT, 
     1                       period2, lnY, sigma, iflag, depth ) 
            attenname1 = 'Garcia et al. (2005), Ver-Inslab, Rock'
         endif
         lnY = lnY + efac
       endif

c ***** Lin and Lee (2008) Subduction Model *****
c      Lin and Lee (2008) - Horizontal, Subduction, Rock
C     Model Number = 280
      if ( jcalc .eq. 280 ) then
         call LinLee08rock ( mag, rupdist, specT, 
     1                    period2, lnY, sigma, iflag, depth, ftype ) 
         attenname1 = 'Lin and Lee (2008), Subduction, Rock'
         lnY = lnY + efac
      endif

c      Lin and Lee (2008) - Horizontal, Subduction, Soil
C     Model Number = 281
      if ( jcalc .eq. 281 ) then
         call LinLee08soil ( mag, rupdist, specT, 
     1                    period2, lnY, sigma, iflag, depth, ftype ) 
         attenname1 = 'Lin and Lee (2008), Subduction, Soil'
         lnY = lnY + efac
      endif

C  **** BCSubduction Model *******
C     Base Case Model, DeltaC1 = 0.0
C     Model Number = 350
      if ( jcalc .eq. 350 ) then
         deltaC1 = 0.0
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         attenname1 = 'BCHydroSub_V3, DeltaC1=0'      
      endif

C     Base Case Model, DeltaC1 = -0.5
C     Model Number = 351
      if ( jcalc .eq. 351 ) then
         deltaC1 = -0.5
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         attenname1 = 'BCHydroSub_V3, DeltaC1=-0.5'
      endif

C     Base Case Model, DeltaC1 = 0.5
C     Model Number = 352
      if ( jcalc .eq. 352 ) then
         deltaC1 = 0.5
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         attenname1 = 'BCHydroSub_V3, DeltaC1=0.5'
      endif

C     Base Case Model, DeltaC1 = 0.0, Single Station Sigma=0.60
C     Model Number = 353
      if ( jcalc .eq. 353 ) then
         deltaC1 = 0.0
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=0'
      endif

C     Base Case Model, DeltaC1 = -0.5, Single Station Sigma=0.60
C     Model Number = 354
      if ( jcalc .eq. 354 ) then
         deltaC1 = -0.5
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=-0.5'
      endif

C     Base Case Model, DeltaC1 = 0.5, Single Station Sigma=0.58
C     Model Number = 355
      if ( jcalc .eq. 355 ) then
         deltaC1 = 0.5
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=0.5'
      endif

C     Base Case Model, DeltaC1 = -0.2
C     Model Number = 356
      if ( jcalc .eq. 356 ) then
         deltaC1 = -0.2
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         attenname1 = 'BCHydroSub_V3, DeltaC1=-0.2'
      endif

C     Base Case Model, DeltaC1 = 0.2
C     Model Number = 357
      if ( jcalc .eq. 357 ) then
         deltaC1 = 0.2
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         attenname1 = 'BCHydroSub_V3, DeltaC1=0.2'
      endif

C     Base Case Model, DeltaC1 = -0.2, Single Station Sigma=0.60
C     Model Number = 358
      if ( jcalc .eq. 358 ) then
         deltaC1 = -0.2
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3-SSS, DeltaC1=-0.2'
      endif

C     Base Case Model, DeltaC1 = 0.2, Single Station Sigma=0.60
C     Model Number = 359
      if ( jcalc .eq. 359 ) then
         deltaC1 = 0.2
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
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
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
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
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
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
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
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
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
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
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
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
         call BCHydroSub_V3 ( mag, ftype, rupDist, vs, lnY, 
     1            sigma, specT, period2, iflag, foreArc, depth, disthypo, deltaC1 ) 
         sigma = 0.60
         attenname1 = 'BCHydroSub_V3, Var. Upper DeltaC1, SSS'      
      endif

c ***** Atkinson and Macias (2009) Cascadia Model, NEHRP B/C *****
C     Model Number = 370
      if ( jcalc .eq. 370 ) then
         call AM09_Cas ( mag, rupDist, lnY, 
     1            sigma, specT, period2, iflag ) 
         attenname1 = 'Atkinson&Macias, Cascadia, NEHRP B/C'      
      endif
      
C ******  CEUS Models *********
C
c ***** Atkinson and Boore Models *****
c     Atkinson and Boore (1994) - Horizontal, EUS Hard Rock
C     Model Number = 100
      if ( jcalc .eq. 100 ) then
         call AB95 ( mag, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag )
         attenname1 = 'Atkinson&Boore 1994, EUS, Rock'
         lnY = lnY + efac
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
          call AB95Mn ( mag1, rupdist, lnY, sigma, specT,
     1            attenName1, period2,iflag )
          attenname1 = 'Atkinson&Boore 1994, EUS, Rock, Nuttli Mag'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2006) - Horizontal, EUS Hard Rock
C     Model Number = 102
      if ( jcalc .eq. 102 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
         attenname1 = 'Atkinson&Boore 2006, EUS, Hard Rock'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2006) - Horizontal, EUS Vs=760m/sec
C     Model Number = 103
      if ( jcalc .eq. 103 ) then
         call AB06vs760 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
         attenname1 = 'Atkinson&Boore 2006, EUS, Vs=760m/sec'
         lnY = lnY + efac
      endif

c     Atkinson (2008, weighted C0) - Horizontal, CEUS-NGA nga Vs=760m/sec
C     Model Number = 104 
      if ( jcalc .eq. 104 ) then
c     First call BA08 to get ground motion values which will be adjusted. 
         vs = 760.0
         call BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, BA08lnY, sigma, iflag, vs, ftype, pga4nl,
     1                    phi, tau )
C     Now Call A08vs760 to adjust BA08 NGA value to CEUS Vs=760 value.
         call A08vs760 ( mag, jbdist, specT, BA08lnY,
     1                    period2, lnY, sigma, iflag )
         attenname1 = 'Atkinson 2008 wt C0, EUS-NGA BA08, Vs=760m/sec'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2006) with Atikinson (2010) stress drop adjustment - Horizontal, EUS Hard Rock
C     Model Number = 105
      if ( jcalc .eq. 105 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 10**(3.45 - 0.2*mag)
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006/Atkinson2010, EUS, Hard Rock'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
      endif

c     Atkinson and Boore (2006) with Atikinson (2010) stress drop adjustment - Horizontal, EUS Vs760m/sec
C     Model Number = 106
      if ( jcalc .eq. 106 ) then
         call AB06vs760 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 10**(3.45 - 0.2*mag)
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006/Atkinson2010, EUS, Vs760m/s'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
      endif

c     Atkinson (2010) - Horizontal, CEUS-NGA NGA Based
C     Model Number = 107 
      if ( jcalc .eq. 107 ) then
c     First call BA08 to get ground motion values which will be adjusted. 
         call BA_NGA_2008 ( mag, jbdist, specT,
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
         call Fena ( jbdist, specT, facFena )
         attenname1 = 'Atkinson 2010, EUS-NGA BA08 based'
         lnY =  BA08lnY + factor + FacFena + efac
      endif

c     Atkinson (2008, average C0) - Horizontal, CEUS-NGA nga Vs=760m/sec
C     Model Number = 108
      if ( jcalc .eq. 108 ) then
c     First call BA08 to get ground motion values which will be adjusted. 
         vs = 760.0
         call BA_NGA_2008 ( mag, jbdist, specT,
     1                    period2, BA08lnY, sigma, iflag, vs, ftype, pga4nl,
     1                    phi, tau )
C     Now Call A08vs760 to adjust BA08 NGA value to CEUS Vs=760 value.
         call A08vs760C0 ( mag, jbdist, specT, BA08lnY,
     1                    period2, lnY, sigma, iflag )
         attenname1 = 'Atkinson 2008 avg C0, EUS-NGA BA08, Vs=760m/sec'
         lnY = lnY + efac
      endif

c     Atkinson and Boore (2006) with 2x stress drop adjustment (280bars) - Horizontal, EUS Hard Rock
C     Model Number = 130
      if ( jcalc .eq. 130 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 280.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (2x) StressDrop, EUS, Hard Rock'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
      endif

c     Atkinson and Boore (2006) with 0.5x stress drop adjustment (70bars) - Horizontal, EUS Hard Rock
C     Model Number = 131
      if ( jcalc .eq. 131 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 70.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (0.5x) StressDrop, EUS, Hard Rock'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
      endif

c ******* Toro et al. Models *******
C Toro et al. (1997) MidCon., Horizontal, Rock
C     Model Number = 110
      if (jcalc .eq. 110) then
         call TAS96 ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Toro et al. (1997), Horizontal, MidCon.'
         lnY = lnY + efac
      endif 

C Toro et al. (1997) MidCon., Horizontal, Rock, MLg magnitude
C     Model Number = 111
      if (jcalc .eq. 111) then
         call TAS96MLg ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag)                                    
         attenname1 = 'Toro et al. (1997), Horizontal, MidCon., MLg'
         lnY = lnY + efac
      endif 

C Toro et al. (1997) Gulf, Horizontal, Rock
C     Model Number = 112
      if (jcalc .eq. 112) then
         call TAS96Gulf ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Toro et al. (1997), Horizontal, Gulf'
         lnY = lnY + efac
      endif 

C Toro et al. (1997) Gulf, Horizontal, Rock, MLg magnitude
C     Model Number = 113
      if (jcalc .eq. 113) then
         call TAS96GulfMLg ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Toro et al. (1997), Horizontal, Gulf, MLg'
         lnY = lnY + efac
      endif 

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal, Hard Rock
C     Model Number = 120
      if (jcalc .eq. 120) then
         call CHY03 ( mag, rupdist, lnY, sigma, specT,                             
     1                  period2,iflag )                                    
         attenname1 = 'Campbell (2003), Hor., CEUS-Hybrid, Hard Rock'
         lnY = lnY + efac
      endif 

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal minus SigmaEps, Hard Rock
C     Model Number = 121
      if (jcalc .eq. 121) then
         call CHY03 ( mag, rupdist, lnY, sigma, specT,                             
     1                  period2,iflag )     
         call CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )                                                
         attenname1 = 'Campbell (2003), Hor - SigmaEps, CEUS-Hybrid, Hard Rock'
         lnY = lnY - sigmaEps + efac
      endif 

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal plus SigmaEps, Hard Rock
C     Model Number = 122
      if (jcalc .eq. 122) then
         call CHY03 ( mag, rupdist, lnY, sigma, specT,                             
     1                  period2,iflag )     
         call CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )                                                
         attenname1 = 'Campbell (2003), Hor + SigmaEps, CEUS-Hybrid, Hard Rock'
         lnY = lnY + sigmaEps + efac
      endif 

c ******* PE&A CEUS Models *******
C Silva et al. (2002) 2 Corner, Rock 
C     Model Number = 401
      if (jcalc .eq. 401) then
         call PEA2C ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 2-Corner, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 2 Corner-Saturation, Rock 
C     Model Number = 402
      if (jcalc .eq. 402) then
         call PEA2CS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 2-Corner-Sat, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-High, Rock 
C     Model Number = 403
      if (jcalc .eq. 403) then
         call PEA1CVH ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-High, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Medium, Rock 
C     Model Number = 404
      if (jcalc .eq. 404) then
         call PEA1CVM ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Med, Horizontal'
         lnY = lnY + efac
      endif 
     
C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Low, Rock 
C     Model Number = 405
      if (jcalc .eq. 405) then
         call PEA1CVL ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Low, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High, Rock 
C     Model Number = 406
      if (jcalc .eq. 406) then
         call PEA1CCH ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Medium, Rock 
C     Model Number = 407
      if (jcalc .eq. 407) then
         call PEA1CCM ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low, Rock 
C     Model Number = 408
      if (jcalc .eq. 408) then
         call PEA1CCL ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High-Sat, Rock 
C     Model Number = 409
      if (jcalc .eq. 409) then
         call PEA1CCHS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High-Sat, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Med-Sat, Rock 
C     Model Number = 410
      if (jcalc .eq. 410) then
         call PEA1CCMS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med-Sat, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low-Sat, Rock 
C     Model Number = 411
      if (jcalc .eq. 411) then
         call PEA1CCLS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low-Sat, Horizontal'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Gulf Models
C Silva et al. (2002) 2 Corner, Rock, Gulf 
C     Model Number = 501
      if (jcalc .eq. 501) then
         call PEAG2C ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 2-Corner, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 2 Corner-Saturation, Rock, Gulf 
C     Model Number = 502
      if (jcalc .eq. 502) then
         call PEAG2CS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 2-Corner-Sat, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-High, Rock, Gulf 
C     Model Number = 503
      if (jcalc .eq. 503) then
         call PEAG1CVH ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-High, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Medium, Rock, Gulf
C     Model Number = 504
      if (jcalc .eq. 504) then
         call PEAG1CVM ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Med, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Low, Rock, Gulf 
C     Model Number = 505
      if (jcalc .eq. 505) then
         call PEAG1CVL ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Low, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High, Rock, Gulf 
C     Model Number = 506
      if (jcalc .eq. 506) then
         call PEAG1CCH ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Medium, Rock, Gulf 
C     Model Number = 507
      if (jcalc .eq. 507) then
         call PEAG1CCM ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low, Rock, Gulf 
C     Model Number = 508
      if (jcalc .eq. 508) then
         call PEAG1CCL ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-High-Sat, Rock, Gulf 
C     Model Number = 509
      if (jcalc .eq. 509) then
         call PEAG1CCHS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-High-Sat, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Med-Sat, Rock, Gulf 
C     Model Number = 510
      if (jcalc .eq. 510) then
         call PEAG1CCMS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Med-Sat, Hor. Gulf'
         lnY = lnY + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Constant-Low-Sat, Rock, Gulf 
C     Model Number = 511
      if (jcalc .eq. 511) then
         call PEAG1CCLS ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         attenname1 = 'Silva et al. (2002), 1-Corner-Const-Low-Sat, Hor. Gulf'
         lnY = lnY + efac
      endif 

C *****  Misc Models ******
c     McVerry et al (1993) new zealand
C     Model Number = 300
      if ( jcalc .eq. 300 .and. specT .eq. 0.0 ) then
        call mcverry93 ( mag, rupDist, lnY, sigma, attenName1,
     1       ftype )
          period2 = 0.0
          iflag = 0
         lnY = lnY + efac
      elseif ( jcalc .eq. 300 .and. specT .ne. 0.0) then
          write (*,*) 'McVerry (1993), Horizontal, Rock'
          write (*,*) 'not defined for spectral acclereation!!!'
          write (*,*) 'Check input file.'
          stop 99
      endif

c     fukushima (1990) rock
C     Model Number = 301
      if ( jcalc .eq. 301 .and. specT .eq. 0.0 ) then
          call fukushima90 ( mag, rupDist, lnY, sigma, attenName1)
          period2 = 0.0
          iflag = 0
         lnY = lnY + efac
      elseif (jcalc .eq. 301 .and. specT .ne. 0.0) then
          write (*,*) 'Fukushima (1990), Horizontal, Rock'
          write (*,*) 'not defined for spectral acclereation!!!'
          write (*,*) 'Check input file.'
          stop 99
      endif

c     Loh high speed rail (New Joyner-Boore form)
C     Model Number = 302      
      if ( jcalc .eq. 302 .and. specT .eq. 0.0 ) then
         call HighSpeedRail ( mag, rupDist, lnY, sigma, 
     1                      attenName1, period2 )
         iflag = 0
         lnY = lnY + efac
      elseif (jcalc .eq. 302 .and. specT .ne. 0.0 ) then
         write (*,*) 'Loh High Speed Rail, Horizontal, Rock'
         write (*,*) 'not defined for spectral acclereation!!!'
         write (*,*) 'Check input file.'
         stop 99
      endif
          
c     New Loh (1996) model (unpublished)
C     Model Number = 303
      if ( jcalc .eq. 303 .and. specT .eq. 0.0) then
         call Loh96 ( mag, rupDist, lnY, sigma, attenName1, period2)
         iflag = 0
         lnY = lnY + efac
      elseif (jcalc .eq. 303 .and. specT .ne. 0.0 ) then
         write (*,*) 'Loh (1996), Horizontal, Rock'
         write (*,*) 'not defined for spectral acclereation!!!'
         write (*,*) 'Check input file.'
         stop 99
      endif

c ******* Ambraseys et al 2005 Model *********
C     Model Number = 601
      if ( jcalc .eq. 601 ) then
         call Ambraseys_2005 ( mag, jbDist, ftype, specT,
     1                     period2, lnY, sigma, iflag )
         attenname1 = 'Ambraseys_et_al_2005_Hor'
         lnY = lnY + efac
      endif


C     **** BCHydro SCR attenuation models adjusted for Vs=760m/sec ****

c     Atkinson and Boore (2006) - Horizontal, CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1020
      if ( jcalc .eq. 1020 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Atkinson&Boore 2006, EUS, BCH Amps for Vs760m/s'
         lnY = lnSa + efac
      endif

c     Atkinson and Boore (2006) with Atikinson (2010) stress drop adjustment - Horizontal, EUS BCH Amps for Vs760m/s
C     Model Number = 1050
      if ( jcalc .eq. 1050 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 10**(3.45 - 0.2*mag)
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006/Atkinson2010, EUS, BCH Amps for Vs760m/s'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c     Atkinson and Boore (2006) with 2x stress drop adjustment (280bars) - Horizontal, EUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1300
      if ( jcalc .eq. 1300 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 280.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (2x) StressDrop, EUS, BCH Amps for Vs760'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c     Atkinson and Boore (2006) with 0.5x stress drop adjustment (70bars) - Horizontal, EUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1310
      if ( jcalc .eq. 1310 ) then
         call AB06 ( mag, rupdist, lnY, sigma, specT,
     1            period2,iflag )
C     Now apply compute the SF2 factor from the AB06 model.  
         call AB06SF2 (mag, specT, period2, SF2)
c     Now compute the scale factor for different magnitude dependent stress drops.
         sd = 70.0
         sdscale = log10(sd/140.0)/log10(2.0)
         attenname1 = 'Atkinson&Boore 2006 (0.5x) StressDrop, EUS, BCH Amps for Vs760'
         LnY = lnY + sdscale*SF2*alog(10.0) + efac
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal, BCH Amps for Vs760m/s
C     Model Number = 1200
      if (jcalc .eq. 1200) then
         call CHY03 ( mag, rupdist, lnY, sigma, specT,                             
     1                  period2,iflag )                                    
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Campbell (2003), Hor., CEUS-Hybrid, BCH Amps for Vs760m/s'
         lnY = lnSa + efac
      endif 

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal minus SigmaEps, BCH Amps for Vs760m/s
C     Model Number = 1210s
      if (jcalc .eq. 1210) then
         call CHY03 ( mag, rupdist, lnY, sigma, specT,                             
     1                  period2,iflag )     
         call CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )                                                
         attenname1 = 'Campbell (2003), Hor - SigmaEps, CEUS-Hybrid, BCH Amps for Vs760m/s'
         lnY = lnY - sigmaEps + efac
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif 

c ******* Campbell Hybrid CEUS Models *******
C Campbell (2003) Hybrid CEUS Horizontal plus SigmaEps, BCH Amps for Vs760m/s
C     Model Number = 1220
      if (jcalc .eq. 1220) then
         call CHY03 ( mag, rupdist, lnY, sigma, specT,                             
     1                  period2,iflag )     
         call CHY03Eps ( mag, rupdist, sigmaeps, specT, period2, iflag )                                                
         attenname1 = 'Campbell (2003), Hor + SigmaEps, CEUS-Hybrid, BCH Amps for Vs760m/s'
         lnY = lnY + sigmaEps + efac
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         LnY = LnSa
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-High, BCH Amps for Vs760m/s 
C     Model Number = 4030
      if (jcalc .eq. 4030) then
         call PEA1CVH ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-High, Hor,BCH Amps for Vs760m/s'
         lnY = lnSA + efac
      endif 

C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Medium, BCH Amps for Vs760m/s 
C     Model Number = 4040
      if (jcalc .eq. 4040) then
         call PEA1CVM ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Med, Hor,BCH Amps for Vs760m/s'
         lnY = lnSa + efac
      endif 
     
C New PE&A CEUS Models
C Silva et al. (2002) 1 Corner Variable-Low, BCH Amps for Vs760m/s 
C     Model Number = 4050
      if (jcalc .eq. 4050) then
         call PEA1CVL ( mag, jbdist, lnY, sigma, specT,                             
     1                  attenName1, period2,iflag )                                    
         call BCHHR2Vs760 ( lnY, specT, lnSa )
         attenname1 = 'Silva et al. (2002), 1-Corner-Var-Low, Hor,BCH Amps for Vs760m/s'
         lnY = lnSa + efac
      endif 

C     **** End of BCHydro SCR attenuation models adjusted for Vs=760m/sec ****

C     Akkar and Cayan (2010) - "A Local Ground-Motion Predictive Model for Turkey, and
C           its Comparison with Other Regional and Global Ground-Motion Models"
C     BSSA: Vol. 100, No.6, pp. 2978-2995.
C     Model Number = 150
      if ( jcalc .eq. 150 ) then
         call AC_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, vs, ftype, pga4nl )
         attenname1 = 'Akkar&Cagan_2010_Hor'
         lnY = lnY + efac
       endif

C     Akkar and Bommer (2010) - "Empirical Equations for the Prediction of PGA, PGV,
C             and Spectral Accelerations in Europe, the Mediterranean Region and 
C             the Middle East"
C     SRL: March/April, Vol.81, pp. 195-206
C     Rock Site conditions Vs>750m/s
C     Model Number = 151
      if ( jcalc .eq. 151 ) then 
         Ss = 0.0
         Sa = 0.0   
         call AB_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Ss, Sa )
         attenname1 = 'Akkar&Bommer_2010_Hor, Rock'
         lnY = lnY + efac
       endif

C     Stiff Soils Site conditions 360<Vs<750m/s
C     Model Number = 152
      if ( jcalc .eq. 152 ) then 
         Ss = 0.0
         Sa = 1.0   
         call AB_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Ss, Sa )
         attenname1 = 'Akkar&Bommer_2010_Hor, Stiff Soil'
         lnY = lnY + efac
       endif
C     Soft Soil Site conditions Vs<360m/s
C     Model Number = 153
      if ( jcalc .eq. 153 ) then 
         Ss = 1.0
         Sa = 0.0   
         call AB_2010 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Ss, Sa )
         attenname1 = 'Akkar&Bommer_2010_Hor, Soft Soil'
         lnY = lnY + efac
       endif

C     Akkar, Sandikkaya, and Bommer (2013) - "Empirical ground-motion models for
C             point- and extended-source crustal earthquake scenarios in 
C             Europe and the Middle East"
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
         call ASB_2013 ( mag, jbdist, specT,
     1                    period2, lnY, sigma, iflag, ftype, Vs, phi, tau )
         attenname1 = 'Akkar,Sandikkaya&Bommer_2013_Hor'
         lnY = lnY + efac
       endif





c ******* Bradley (2010) *********
C     Bradley (2010) - "NZ-Specific Psuedo-Spectral Acceleration Ground
C         Motion Prediction Equation Based on Foreign Models"
C     University of Canterbury Research Report 2010-03
C     Model Number = 160
      if ( jcalc .eq. 160 ) then
         call Bradley_2010 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, hwflag, Rx )
         attenname1 = 'Bradley-2010, Hor'
         lnY = lnY + efacs
       endif

c ******* McVerry et al. (2006) Crustal Events, Horizontal *********
C     McVerry, G. H., J. X. Zhao, N.A. Abrahamson and P.G. Somerville (2006) 
C        "New Zealand Acceleration Response Spectrum Attenuation 
C         Reations for Crustal and Subduction Zone Earthquakes",
C         Bulletin New Zealand Society for Earthquake Engineering, Vol. 39,
C         No. 1, pp. 1 - 58. 
C     Model Number = 140, Site Class A/B
      if ( jcalc .eq. 140 ) then
         Sc = 0.0
         Sd = 0.0
         call McVerry_Crustal_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype, 
     3                     hwflag, Sc, Sd)
         attenname1 = 'McVerryetal-Crustal-2006, Hor, Site A/B'
         lnY = lnY + efacs
       endif
C     Model Number = 141, Site Class C
      if ( jcalc .eq. 141 ) then
         Sc = 1.0
         Sd = 0.0
         call McVerry_Crustal_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype, 
     3                     hwflag, Sc, Sd)
         attenname1 = 'McVerryetal-Crustal-2006, Hor, Site C'
         lnY = lnY + efacs
       endif
C     Model Number = 142, Site Class D
      if ( jcalc .eq. 142 ) then
         Sc = 0.0
         Sd = 1.0
         call McVerry_Crustal_2006 ( mag, rupDist, specT,
     1                     period2, lnY, sigma, iflag, Ftype, 
     3                     hwflag, Sc, Sd)
         attenname1 = 'McVerryetal-Crustal-2006, Hor, Site D'
         lnY = lnY + efacs
       endif

c ******* Bindi et al. (2009) Crustal Events, Horizontal *********
C     Bindi, D., L. Luzi, M. Massa, and F. Pacor (2009). "Horizontal and Vertical 
C          ground motion prediction equaitons derived from the Italian
C          Accelerometric Archive (ITACA)", Bull. Earthquake Eng., 
C          DOI 10.1007/s10518-009-9130-9
C      Note: Only the horizontal model is currently coded and the coefficients
C            independent of fault mechanism as preferred by the authors. Coefficients
C            are for the model as a function of JBDist.
C     Model Number = 95, Rock, Horizontal
      if ( jcalc .eq. 95 ) then
         Sr = 1.0
         Ss = 0.0
         Sd = 0.0
         call Bindi_Hor_2009 ( mag, jbdist, specT,
     1                     period2, lnY, sigma, iflag, Sr, Ss, Sd )
         attenname1 = 'Bindietal-Hor-2009, Rock'
         lnY = lnY + efacs
       endif
C     Model Number = 96, Shallow Alluvium
      if ( jcalc .eq. 96 ) then
         Sr = 0.0
         Ss = 1.0
         Sd = 0.0
         call Bindi_Hor_2009 ( mag, jbdist, specT,
     1                     period2, lnY, sigma, iflag, Sr, Ss, Sd )
         attenname1 = 'Bindietal-Hor-2009, Shallow Alluvium'
         lnY = lnY + efacs
       endif
C     Model Number = 97, Deep Alluvium
      if ( jcalc .eq. 97 ) then
         Sr = 0.0
         Ss = 0.0
         Sd = 1.0
         call Bindi_Hor_2009 ( mag, jbdist, specT,
     1                     period2, lnY, sigma, iflag, Sr, Ss, Sd )
         attenname1 = 'Bindietal-Hor-2009, Deep Alluvium'
         lnY = lnY + efacs
       endif

c ******* Bindi et al. (2011) Crustal Events, Horizontal *********
C     Bindi, D., F. Pacor, L. Luzi, R. Puglia, M. Massa, G. Ameri and R. Paolucci (2011). 
C          "Ground motion prediction equations derived from the Italian 
C          strong motion database."
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
         call Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class A(Vs>800m/s)'
         lnY = lnY + efacs
       endif
C     Model Number = 196, Class B (Vs=360-800m/s)
      if ( jcalc .eq. 196 ) then
         SCa = 0.0
         SCb = 1.0
         SCc = 0.0
         SCd = 0.0
         SCe = 0.0
         call Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class B(360<Vs<800m/s)'
         lnY = lnY + efacs
       endif

C     Model Number = 197, Class C (Vs=180-360m/s)
      if ( jcalc .eq. 197) then
         SCa = 0.0
         SCb = 0.0
         SCc = 1.0
         SCd = 0.0
         SCe = 0.0
         call Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class C(180<Vs<360m/s)'
         lnY = lnY + efacs
       endif

C     Model Number = 198, Class D (Vs<180)
      if ( jcalc .eq. 198) then
         SCa = 0.0
         SCb = 0.0
         SCc = 0.0
         SCd = 1.0
         SCe = 0.0
         call Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class D(Vs<180)'
         lnY = lnY + efacs
       endif

C     Model Number = 199, Class E (Vs=180-360m/s)
      if ( jcalc .eq. 199) then
         SCa = 0.0
         SCb = 0.0
         SCc = 0.0
         SCd = 0.0
         SCe = 1.0
         call Bindi_Hor_2011 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, SCa, SCb, SCc, SCD, SCe, phi, tau )
         attenname1 = 'Bindietal-Hor-2011, Class E(5-20m Class D/E over ClassA)'
         lnY = lnY + efacs
       endif

c ******* Bindi et al. (2013) Crustal Events, Horizontal *********
C     Bindi, D., M. Massa, L. Luzi, G. Ameri, R. Puglia, and P. Augliera (2013). 
C          "Pan-European Ground-Motion Prediction Equations for the Average
C           Horizontal Component of PGA, PGV, and 5%-Damped PSA at Spectral
C           Periods up to 3.0 s using the RESORCE dataset"
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
         call Bindi_Hor_2013 ( mag, jbdist, ftype, specT,
     1                     period2, lnY, sigma, iflag, vs, phi, tau )
         attenname1 = 'Bindietal-Hor-2013, Rjb, Vs'
         lnY = lnY + efacs
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
         call GK_Nov2012 ( mag, RupDist, specT, ftype,
     1                     period2, lnY, sigma, iflag, Vs, Q0, depthvs15 )
         attenname1 = 'Graizer&Kalkan, Nov.2012'
         lnY = lnY + efacs
       endif
      if ( jcalc .eq. 91 ) then
         Q0 = 75.0
         call GK_Nov2012 ( mag, RupDist, specT, ftype,
     1                     period2, lnY, sigma, iflag, Vs, Q0, depthvs15 )
         attenname1 = 'Graizer&Kalkan, Nov.2012, Q0=75'
         lnY = lnY + efacs
       endif
      if ( jcalc .eq. 92 ) then
         Q0 = 300.0
         call GK_Nov2012 ( mag, RupDist, specT, ftype,
     1                     period2, lnY, sigma, iflag, Vs, Q0, depthvs15 )
         attenname1 = 'Graizer&Kalkan, Nov.2012, Q0=300'
         lnY = lnY + efacs
       endif


c ******* DCPP Common Function Form Models *********
C     Preliminary Models for DCPP for Workshop 3, March 2014
C     Model Numbers 8001, ASK form
      if ( jcalc .eq. 8001 ) then
         call DCPP_CommonASK ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model, ASK'
         lnY = lnY + efacs
       endif
C     Model Numbers 8002, BSSA form
      if ( jcalc .eq. 8002 ) then
         call DCPP_CommonBSSA ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model, BSSA'
         lnY = lnY + efacs
       endif
C     Model Numbers 8003, Common Model 001
      if ( jcalc .eq. 8003 ) then
         call DCPP_Common001 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model001'
         lnY = lnY + efacs
       endif
C     Model Numbers 8004, Common Model 002
      if ( jcalc .eq. 8004 ) then
         call DCPP_Common002 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model002'
         lnY = lnY + efacs
       endif
C     Model Numbers 8005, Common Model 003
      if ( jcalc .eq. 8005 ) then
         call DCPP_Common003 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model003'
         lnY = lnY + efacs
       endif
C     Model Numbers 8006, Common Model 004
      if ( jcalc .eq. 8006 ) then
         call DCPP_Common004 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model004'
         lnY = lnY + efacs
       endif
C     Model Numbers 8007, Common Model 005
      if ( jcalc .eq. 8007 ) then
         call DCPP_Common005 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model005'
         lnY = lnY + efacs
       endif
C     Model Numbers 8008, Common Model 006
      if ( jcalc .eq. 8008 ) then
         call DCPP_Common006 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model006'
         lnY = lnY + efacs
       endif
C     Model Numbers 8009, Common Model 007
      if ( jcalc .eq. 8009 ) then
         call DCPP_Common007 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model007'
         lnY = lnY + efacs
       endif
C     Model Numbers 8010, Common Model 008
      if ( jcalc .eq. 8010 ) then
         call DCPP_Common008 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model008'
         lnY = lnY + efacs
       endif
C     Model Numbers 8011, Common Model 009
      if ( jcalc .eq. 8011 ) then
         call DCPP_Common009 ( mag, RupDist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'DCPP Common Model009'
         lnY = lnY + efacs
       endif

c ******* PVNGS Common Function Form Models *********
C     Preliminary Models for PVNGS for Workshop 3, March 2014

C     Model Numbers 9001 Common Model ASK
      if ( jcalc .eq. 9001 ) then
         call PVNGS_CommonASK ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model ASK'
         lnY = lnY + efacs
      endif
C     Model Numbers 9002 Common Model ASK
      if ( jcalc .eq. 9002 ) then
         call PVNGS_CommonBindi ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model Bindi'
         lnY = lnY + efacs
      endif
C     Model Numbers 9003 Common Model ASK
      if ( jcalc .eq. 9003 ) then
         call PVNGS_Common001 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model001'
         lnY = lnY + efacs
      endif
C     Model Numbers 9004 Common Model ASK
      if ( jcalc .eq. 9004 ) then
         call PVNGS_Common002 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model002'
         lnY = lnY + efacs
      endif
C     Model Numbers 9005 Common Model ASK
      if ( jcalc .eq. 9005 ) then
         call PVNGS_Common003 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model003'
         lnY = lnY + efacs
      endif
C     Model Numbers 9006 Common Model ASK
      if ( jcalc .eq. 9006 ) then
         call PVNGS_Common004 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model004'
         lnY = lnY + efacs
      endif
C     Model Numbers 9007 Common Model ASK
      if ( jcalc .eq. 9007 ) then
         call PVNGS_Common005 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model005'
         lnY = lnY + efacs
      endif
C     Model Numbers 9008 Common Model ASK
      if ( jcalc .eq. 9008 ) then
         call PVNGS_Common006 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model006'
         lnY = lnY + efacs
      endif
C     Model Numbers 9009 Common Model ASK
      if ( jcalc .eq. 9009 ) then
         call PVNGS_Common007 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model007'
         lnY = lnY + efacs
      endif
C     Model Numbers 9010 Common Model ASK
      if ( jcalc .eq. 9010 ) then
         call PVNGS_Common008 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model008'
         lnY = lnY + efacs
      endif
C     Model Numbers 9011 Common Model ASK
      if ( jcalc .eq. 9011 ) then
         call PVNGS_Common009 ( mag, jbdist, depthtop, specT, lnY, sigma, iflag )
         attenname1 = 'PVNGS Common Model009'
         lnY = lnY + efacs
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
         call SWUS_CFRrup ( mag, RupDist, jbDist, depthtop, ftype, dip1(1), RupWidth, Rx, HWFlag, 
     1           specT, lnY, sigma, iflag, cfcoefrrup, coefcountrrup, phi, tau )
         attenname1 = 'SWUS Common Function Model-Rrup'
         lnY = lnY + efacs
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
         call SWUS_CFRjb ( mag, RupDist, jbDist, depthtop, ftype, dip1(1), RupWidth, Rx, HWFlag, 
     1           specT, lnY, sigma, iflag, cfcoefrjb, coefcountrjb, phi, tau )
         attenname1 = 'SWUS Common Function Model-Rjb'
         lnY = lnY + efacs
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
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
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
         
         lnY = lnY + efac
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     Low Mag Uncertainty
C     Model Number = 8788
      if ( jcalc .eq. 8788 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
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
         
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     High Mag Uncertainty
C     Model Number = 8789
      if ( jcalc .eq. 8789 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
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
         
         lnY = lnY + efac
      endif


c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 8922
      if ( jcalc .eq. 8922 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
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

         lnY = lnY + efac
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 8923
      if ( jcalc .eq. 8923 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
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

         lnY = lnY + efac
       endif
       
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     High Mag Uncertainty
C     Model Number = 8924
      if ( jcalc .eq. 8924 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
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

         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 8836
      if ( jcalc .eq. 8836 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
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

         lnY = lnY + efac
       endif
       
c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 8837
      if ( jcalc .eq. 8837 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
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

         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     High Mag Uncertainty
C     Model Number = 8838
      if ( jcalc .eq. 8838 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
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

         lnY = lnY + efac
       endif



c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 8797
      if ( jcalc .eq. 8797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
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

         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 8798
      if ( jcalc .eq. 8798 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
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

         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 8799
      if ( jcalc .eq. 8799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
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

         lnY = lnY + efac
       endif


c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 8910
      if ( jcalc .eq. 8910 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Cent'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
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

         lnY = lnY + efac
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 8911
      if ( jcalc .eq. 8911 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Low'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
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

         lnY = lnY + efac
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 8912
      if ( jcalc .eq. 8912 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-High'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
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

         lnY = lnY + efac
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
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-Cent-EstVs'
         
C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         
         
         lnY = lnY + efac
      endif


c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     Low Mag Uncertainty
C     Model Number = 9788
      if ( jcalc .eq. 9788 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-Low-EstVs'
         
C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         
         
         lnY = lnY + efac
      endif

c ******* Abrahamson, Silva, and Kamai Model *********
C     Abrahamson, Silva, and Kamai 2014 - Horizontal, Global, Mainshock, Estimated Vs30m
C     High Mag Uncertainty
C     Model Number = 9789
      if ( jcalc .eq. 9789 ) then
         regionflag = 0
         msasflag = 0
         eta = 0.
C     Dip angle of the fault is assumed to be equal to the dip angle
C     assigned for the first two points of the fault.
         dip = dip1(1)
C     Depth to the top of the fault is set equal to the first depth
C     defined for the fault in the data file (i.e., Depthtop). 
         vs30_Class = 0
         call ASK_NGAWest2_2013 ( mag, dip, ftype, Rupwidth, rupDist, jbdist,
     1            vs, hwflag, lnY, sigma, specT, period2, depthtop, iflag, 
     2            vs30_class, depthvs10, Rx, Ry0, regionflag, msasflag, phi, tau) 
         attenname1 = 'ASK_NGAWest2_2014-Hor-Zone123-High-EstVs'
         
C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         
         
         lnY = lnY + efac
      endif


c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 9922
      if ( jcalc .eq. 9922 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-Cent, No Basin'
         
C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif

c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     Central Mag Uncertainty
C     Model Number = 9923
      if ( jcalc .eq. 9923 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-Low, No Basin'
         
C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif
       
c ******* Boore and Atkinson Model *********
C     Boore, Stewart, Seyhan and Atkinson 2014 - horizontal
C          DeltaC3 Model Global Adjustments, No Basin Adjustments
C     High Mag Uncertainty
C     Model Number = 9924
      if ( jcalc .eq. 9924 ) then
         regionflag = 0
         basinflag = 0
         call BSSA_NGAWest2_2013 ( mag, jbdist, specT,
     1               period2, lnY, sigma, iflag, vs, ftype, pga4nl, depthvs10, regionflag, basinflag,
     1               phi, tau )
         attenname1 = 'BSSA_NGAWest2_2014_Hor, Zone123-High, No Basin'
         
C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 9836
      if ( jcalc .eq. 9836 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-Cent'
         
C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif
       
c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     Central Mag Uncertainty
C     Model Number = 9837
      if ( jcalc .eq. 9837 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-Low'
         
C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif

c ******* Campbell and Bozorgnia Model *********
C     Campbell and Bozorgnia 2014 - horizontal, Zone1,2,3
C     High Mag Uncertainty
C     Model Number = 9838
      if ( jcalc .eq. 9838 ) then
         regionflag = 0
         call CB_NGAWest2_2013 ( mag, rupdist, jbdist, ftype, specT,
     1                    period2, lnY, sigma, iflag, vs, 
     2                    depthTop, D25, dip1(1), depth, HWFlag, Rx, rupwidth, regionflag,
     1                    phi, tau )
         attenname1 = 'CB_NGAWest2_2014-Hor,Zone123-High'
         
C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif



c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 9797
      if ( jcalc .eq. 9797 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-Cent,Est Vs30m'
         
C     Apply Mag Uncertainty
         c1 = 0.0
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 9798
      if ( jcalc .eq. 9798 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-Low,Est Vs30m'
         
C     Apply Mag Uncertainty
         c1 = -1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 9799
      if ( jcalc .eq. 9799 ) then
c     Current model set for estimated Vs30 values (only impacts sigma)
         vs30_class = 0
         regionflag = 0
         call CY_NGAWest2_2013 ( mag, rupDist, jbdist, specT,
     1                     period2, lnY, sigma, iflag, 
     2                     vs, dip1(1), Depthtop, Ftype, 
     3                     depthvs10, vs30_class, hwflag, Rx, regionflag,
     1                     phi, tau )
         attenname1 = 'CY_NGAWest2_2014-Hor,Zone123-High,Est Vs30m'
         
C     Apply Mag Uncertainty
         c1 = 1.6
         faddmag = 0.0
         faddmag = c1*( (0.083+0.056*max(0.0,(mag-7.0)) ) + 0.0171*max(0.0, alog(specT)))
         lnY = lnY + faddmag         

         lnY = lnY + efac
       endif


c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Central Mag Uncertainty
C     Model Number = 9910
      if ( jcalc .eq. 9910 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Cent'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
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

         lnY = lnY + efac
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     Low Mag Uncertainty
C     Model Number = 9911
      if ( jcalc .eq. 9911 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-Low'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
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

         lnY = lnY + efac
       endif

c ******* Idriss Model *********
C     Idriss 2014 - Horizontal
C     High Mag Uncertainty
C     Model Number = 9912
      if ( jcalc .eq. 9912 ) then
         if (vs .ge. 450.0) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, vs, specT,
     1                     period2, lnY, sigma, iflag )
            attenname1 = 'Idriss_NGAWest2_2014_Hor, Zone123-High'
         elseif (vs .gt. 1200) then
            call I_NGAWest2_2013 ( mag, rupDist, ftype, 1200.0, specT,
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

         lnY = lnY + efac
       endif



       
       

c     Check for valid jcalc
      if ( lnY .gt. 1.0e10 ) then
         write (*,'( 2x,''invalid jcalc:'',i7,3f10.4,e12.4,f12.4)') jcalc, mag, 
     1      rupDist, ftype, lnY, sigma
         stop 99
      endif

      attenName(jType,iProb) = attenname1
      period1(jType,iProb) = period2
      intflag(jType,iProb) = iflag
      siga = sigma        

c     Convert to g
      lnY = lnY - 6.89

      return
      end
