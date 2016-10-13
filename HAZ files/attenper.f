C Spectral Attenuation Model Period Subroutine 
C     This subroutine will return the array of 
C     spectral periods fo the computation of 
C     deterministic spectra.

      subroutine attenper ( jcalc, anper, minaper, maxaper )

      implicit none      
      include 'pfrisk.h'
      
      integer jcalc, anper
      real minaper, maxaper, specT

c     Set number of attenuation spectral periods to -99 for 
c     checking of valid jcalc values. 
      anper = -99

C ******* PEER NGA West 2 GMPE Models *********

c ******* Idriss (NGA West2 2013) *********
C     Idriss (NGA West2 2013) - Horizontal
C     Model Number = 2910
      if ( jcalc .eq. 2910 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* BSSA (NGA West2 2013) *********
C     Boore, Stewart, Seyhan and Atkinson (NGA West2 2013) - Horizontal
C        DeltaC3 Global adjustment model, No Basin Adjustment
C     Model Number = 2922
      if ( jcalc .eq. 2922 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* BSSA (NGA West2 2013) *********
C     Boore, Stewart, Seyhan and Atkinson (NGA West2 2013) - Horizontal
C        DeltaC3 China-Turkey adjustment model, No Basin Adjustment
C     Model Number = 2923
      if ( jcalc .eq. 2923 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* BSSA (NGA West2 2013) *********
C     Boore, Stewart, Seyhan and Atkinson (NGA West2 2013) - Horizontal
C        DeltaC3 Itatly-Japan adjustment model, No Basin Adjustment
C     Model Number = 2924
      if ( jcalc .eq. 2924 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif
c ******* BSSA (NGA West2 2013) *********
C     Boore, Stewart, Seyhan and Atkinson (NGA West2 2013) - Horizontal
C        DeltaC3 Global adjustment model, Basin Adjustment
C     Model Number = 2925
      if ( jcalc .eq. 2925 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* BSSA (NGA West2 2013) *********
C     Boore, Stewart, Seyhan and Atkinson (NGA West2 2013) - Horizontal
C        DeltaC3 China-Turkey adjustment model, Basin Adjustment
C     Model Number = 2926
      if ( jcalc .eq. 2926 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* BSSA (NGA West2 2013) *********
C     Boore, Stewart, Seyhan and Atkinson (NGA West2 2013) - Horizontal
C        DeltaC3 Itatly-Japan adjustment model, Basin Adjustment
C     Model Number = 2927
      if ( jcalc .eq. 2927 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Abrahamson, Silva, and Kamai (NGA West2 2013) *********
C     Abrahamson, Silva, and Kamai (NGA West2 2013) - Horizontal
C        Same spectral periods for all variations of the GMPE model
C     Model Number = 2787 - 2794
      if ( jcalc .ge. 2787 .and. jcalc .le. 2794 ) then
         anper = 24
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 3787 - 3794
      if ( jcalc .ge. 3787 .and. jcalc .le. 3794 ) then
         anper = 24
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Campbell and Bozorgnia (NGA West2 2013) *********
C     Campbell and Bozorgnia (NGA West2 2013) - Horizontal, California
C     Model Number = 2836
      if ( jcalc .eq. 2836 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Campbell and Bozorgnia (NGA West2 2013) *********
C     Campbell and Bozorgnia (NGA West2 2013) - Horizontal, Japan
C     Model Number = 2837
      if ( jcalc .eq. 2837 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif
c ******* Campbell and Bozorgnia (NGA West2 2013) *********
C     Campbell and Bozorgnia (NGA West2 2013) - Horizontal, China
C     Model Number = 2838
      if ( jcalc .eq. 2838 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif
c ******* Campbell and Bozorgnia (NGA West2 2013) *********
C     Campbell and Bozorgnia (NGA West2 2013) - Horizontal, Italy
C     Model Number = 2839
      if ( jcalc .eq. 2839 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Chiou and Youngs (NGA West2 2013) *********
C     Chiou and Youngs (NGA West2 2013) - Horizontal
C     Model Number = 2797
      if ( jcalc .eq. 2797 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 2798
      if ( jcalc .eq. 2798 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif

C     Chiou and Youngs (NGA West2 2013) - Horizontal, Japan/Italy Adjustment
C     Model Number = 2799
      if ( jcalc .eq. 2799 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 2800
      if ( jcalc .eq. 2800 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif

C     Chiou and Youngs (NGA West2 2013) - Horizontal, Wenchaun Adjustment
C     Model Number = 2801
      if ( jcalc .eq. 2801 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 2802
      if ( jcalc .eq. 2802 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif


C *** NGA West2 Vertical GMPEs *****
c ******* Gulrce, Kamai, Abrahamson and Silva (NGA West2 2013) *********
C     Gulrce, Kamai, Abrahamson and Silva (NGA West2 2013) - Vertical
C        Same spectral periods for all variations of the GMPE model
C     Model Number = 4787
      if ( jcalc .ge. 4787 ) then
         anper = 24
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 4788
      if ( jcalc .ge. 4788 ) then
         anper = 24
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Gulrce, Kamai, Abrahamson and Silva (NGA West2 2013) *********
C     Gulrce, Kamai, Abrahamson and Silva (NGA West2 2013) - V/H Ratip
C        Same spectral periods for all variations of the GMPE model
C     Model Number = 6787
      if ( jcalc .ge. 6787 ) then
         anper = 24
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 6788
      if ( jcalc .ge. 6788 ) then
         anper = 24
         minaper = 0.01
         maxaper = 10.0
       endif


c ******* SSBA (NGA West2 2013) *********
C     Stewart, Seyhan, Boore, and Atkinson (NGA West2 2013) - Vertical
C        DeltaC3 Global adjustment model, No Basin Adjustment
C     Model Number = 4922
      if ( jcalc .eq. 4922 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* SSBA (NGA West2 2013) *********
C     Stewart, Seyhan, Boore, and Atkinson (NGA West2 2013) - V/H Ratio
C        DeltaC3 Global adjustment model, No Basin Adjustment
C     Model Number = 6922
      if ( jcalc .eq. 6922 ) then
         anper = 107
         minaper = 0.01
         maxaper = 10.0
       endif



c ******* Bozorgnia and Campbell (NGA West2 2013) *********
C     Bozorgnia and Campbell (NGA West2 2013) - Vertical, California
C     Model Number = 4836
      if ( jcalc .eq. 4836 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 4837
      if ( jcalc .eq. 4837 ) then
         anper = 19
         minaper = 0.01
         maxaper = 3.0
       endif

C     Bozorgnia and Campbell (NGA West2 2013) - V/H Ratio, California
C     Model Number = 6836
      if ( jcalc .eq. 6836 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 6837
      if ( jcalc .eq. 6837 ) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Chiou and Youngs (NGA West2 2013) *********
C     Chiou and Youngs (NGA West2 2013) - Vertical
C     Model Number = 4797
      if ( jcalc .eq. 4797 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 4798
      if ( jcalc .eq. 4798 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Chiou and Youngs (NGA West2 2013) *********
C     Chiou and Youngs (NGA West2 2013) - V/H Ratio
C     Model Number = 6797
      if ( jcalc .eq. 6797 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif
C     Model Number = 6798
      if ( jcalc .eq. 6798 ) then
         anper = 25
         minaper = 0.01
         maxaper = 10.0
       endif







C ******* PEER NGA 2008 Attenuation Models ****
c ******* Abrahamson and Silva Models *********
C     Abrahamson&Silva 072007 Model - horizontal, No Soil Depth, Estimated Vs30m
C     Model Number = 787
      if ( jcalc .eq. 787 ) then
         anper = 105
         minaper = 0.01
         maxaper = 10.0
      endif

C     Abrahamson&Silva 072007 Model - horizontal, No Soil Depth, Measured Vs30m
C     Model Number = 788
      if ( jcalc .eq. 788 ) then
         anper = 105
         minaper = 0.01
         maxaper = 10.0
      endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs November 2007 - Horizontal, estimated Vs30m
C     Model Number = 797
c
      if ( jcalc .eq. 797 ) then
          anper = 105
          minaper = 0.01
          maxaper = 10.0
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs November 2007 - Horizontal, measured Vs30m
C     Model Number = 798
c
      if ( jcalc .eq. 798 ) then
          anper = 105
          minaper = 0.01
          maxaper = 10.0
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs November 2007 - Horizontal, estimated Vs30m
C     Southern California Small Magnitude Model
C     Model Number = 799
c
      if ( jcalc .eq. 799 ) then
          anper = 105
          minaper = 0.01
          maxaper = 10.0
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs November 2007 - Horizontal, measured Vs30m
C     Southern California Small Magnitude Model
C     Model Number = 800
c
      if ( jcalc .eq. 800 ) then
          anper = 105
          minaper = 0.01
          maxaper = 10.0
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs November 2007 - Horizontal, estimated Vs30m
C     Central California Small Magnitude Model
C     Model Number = 801
c
      if ( jcalc .eq. 801 ) then
          anper = 105
          minaper = 0.01
          maxaper = 10.0
       endif

c ******* Chiou and Youngs Model *********
C     Chiou and Youngs November 2007 - Horizontal, measured Vs30m
C     Central California Small Magnitude Model
C     Model Number = 802
c
      if ( jcalc .eq. 802 ) then
          anper = 105
          minaper = 0.01
          maxaper = 10.0
       endif


C     Campbell and Bozorgnia May 2007 - horizontal
C     Model Number = 836
      if ( jcalc .eq. 836 ) then
         anper = 21
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Idriss (June 2010) *********
C     Idriss June 2010 - Horizontal, Vs=540-900 m/sec and Vs>900m/sec
C     Model Number = 910
      if ( jcalc .eq. 910 ) then
         anper = 32
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Boore and Atkinson (July 2007) *********
C     Boore and Atkinson July 2007 - horizontal
C     Model Number = 922
      if ( jcalc .eq. 922 ) then
         anper = 21
         minaper = 0.01
         maxaper = 10.0
       endif

c ******* Boore and Atkinson (July 2007/2010) *********
C     Boore and Atkinson July 2007 - horizontal with small magnitude adjustment (Atkinson, 2010)
C     Model Number = 923
      if ( jcalc .eq. 923 ) then
         anper = 21
         minaper = 0.01
         maxaper = 10.0
       endif

C ******* End of PEER NGA 2008 Attenuation Models ****

C     Campbell and Bozorgnia (2003), Horizontal, Firm Soil, SS and Reverse
C     Model Number = 070
      if ( jcalc .eq. 70 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Very Firm Soil, SS and Reverse
C     Model Number = 071
      if ( jcalc .eq. 71 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Soft Rock, SS and Reverse
C     Model Number = 072
      if ( jcalc .eq. 72 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Firm Rock, SS and Reverse
C     Model Number = 073
      if ( jcalc .eq. 73 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Generic Rock, SS and Reverse
C     Model Number = 074
      if ( jcalc .eq. 74 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Horizontal, Generic Soil, SS and Reverse
C     Model Number = 075
      if ( jcalc .eq. 75 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Vertical, Firm Soil, SS and Reverse
C     Model Number = 076
      if ( jcalc .eq. 76 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Vertical, Very Firm Soil, SS and Reverse
C     Model Number = 077
      if ( jcalc .eq. 77 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Vertical, Soft Rock, SS and Reverse
C     Model Number = 078
      if ( jcalc .eq. 78 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Vertical, Firm Rock, SS and Reverse
C     Model Number = 079
      if ( jcalc .eq. 79 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Vertical, Generic Rock, SS and Reverse
C     Model Number = 080
      if ( jcalc .eq. 80 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

C     Campbell and Bozorgnia (2003), Vertical, Generic Soil, SS and Reverse
C     Model Number = 081
      if ( jcalc .eq. 81 ) then
         anper = 15
         minaper = 0.05
         maxaper = 4.0
      endif

c ******* Abrahamson and Silva Models *********
C     Abrahamson&Silva 1997 (Rock) - horizontal
C     Model Number = 001
      if ( jcalc .eq. 1 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Rock) - vertical
C     Model Number = 002
      if ( jcalc .eq. 2 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Soil) - horizontal
C     Model Number = 003
      if ( jcalc .eq. 3 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1995 (Soil) - vertical
C     Model Number = 004
      if ( jcalc .eq. 4 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors.
C     Model Number = 005
      if ( jcalc .eq. 5 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Rock) - vertical
C     Model Number = 006
      if ( jcalc .eq. 6 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Soil) - horizontal with Normal faulting factors.
C     Model Number = 007
      if ( jcalc .eq. 7 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Soil) - vertical with Normal faulting factors
C     Model Number = 008
      if ( jcalc .eq. 8 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors scaled
C           by 1.0/1.67 factor.
C     Model Number = 2005
      if ( jcalc .eq. 2005 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

C     Abrahamson&Silva 1997 (Rock) - horizontal with Normal Faulting factors scaled
C           by 1.67 factor.
C     Model Number = 3005
      if ( jcalc .eq. 3005 ) then
         anper = 27
         minaper = 0.03
         maxaper = 5.0
       endif

c ********* Boore, Joyner and Fumal Models ******************
c     BJF94, Horizontal, Class A
C     Model Number = 010
      if ( jcalc .eq. 10 ) then
         anper = 11
         minaper = 0.1
         maxaper = 2.0

c     BJF94, Horizontal, Class B
C     Model Number = 011
      elseif ( jcalc .eq. 11 ) then
         anper = 11
         minaper = 0.1
         maxaper = 2.0

c     BJF94, Horizontal, Class c
C     Model Number = 012
      elseif ( jcalc .eq. 12 ) then
         anper = 11
         minaper = 0.1
         maxaper = 2.0
      endif

c     BJF97, Horizontal, Vs top 30 meters
C     Model Number = 013
      if ( jcalc .eq. 13 ) then
         anper = 12
         minaper = 0.1
         maxaper = 2.0
      endif

c     BJF97, Horizontal, Vs top 30 meters scaled by factor 1.0/1.67
C     Model Number = 2013
      if ( jcalc .eq. 2013 ) then
         anper = 12
         minaper = 0.1
         maxaper = 2.0
      endif

c     BJF97, Horizontal, Vs top 30 meters scaled by factor 1.67
C     Model Number = 3013
      if ( jcalc .eq. 3013 ) then
         anper = 12
         minaper = 0.1
         maxaper = 2.0
      endif

c ******** Campbell Models ******
c     Campbell (1990), Horizontal, Rock
C     Model Number = 020
      if ( jcalc .eq. 20 ) then
         anper = 16
         minaper = 0.04
         maxaper = 4.0
      endif

c     Campbell (1990) - vertical, Rock
C     Model Number = 021
      if ( jcalc .eq. 21 ) then
         anper = 16
         minaper = 0.04
         maxaper = 4.0
      endif

c     Campbell (1990/1994), Horizontal, Rock
C     Model Number = 022
      if ( jcalc .eq. 22 ) then
         anper = 16
         minaper = 0.04
         maxaper = 4.0
      endif

c     Campbell (1993-1994), Horizontal Soil
C     Model Number = 023
      if ( jcalc .eq. 23 ) then
         anper = 16
         minaper = 0.04
         maxaper = 4.0

c     Campbell (1993-1994), Horizontal Soft Rock
C     Model Number = 024
      elseif ( jcalc .eq. 24 ) then
         anper = 16
         minaper = 0.04
         maxaper = 4.0

c     Campbell (1993-1994), Horizontal Hard Rock
C     Model Number = 025
      elseif ( jcalc .eq. 25 ) then
         anper = 16
         minaper = 0.04
         maxaper = 4.0
      endif

c     Campbell (1997), Horizontal, Soil
C     Model Number = 026
      if ( jcalc .eq. 26 ) then
         anper = 14
         minaper = 0.05
         maxaper = 4.0

c     Campbell (1997), Horizontal, Soft Rock
C     Model Number = 027
      elseif ( jcalc .eq. 27 ) then
         anper = 14
         minaper = 0.05
         maxaper = 4.0

c     Campbell (1997), Horizontal, Hard Rock
C     Model Number = 028
      elseif ( jcalc .eq. 28 ) then
         anper = 14
         minaper = 0.05
         maxaper = 4.0
      endif

c     Campbell (1997)  vertical, Soil
C     Model Number = 029
      if ( jcalc .eq. 29 ) then
         anper = 14
         minaper = 0.05
         maxaper = 4.0

c     Campbell (1997)  vertical, Soft Rock
C     Model Number = 030
      elseif ( jcalc .eq. 30 ) then
         anper = 14
         minaper = 0.05
         maxaper = 4.0

c     Campbell (1997)  vertical, Hard Rock
C     Model Number = 031
      elseif ( jcalc .eq. 31 ) then
         anper = 14
         minaper = 0.05
         maxaper = 4.0
      endif

c ******** Idriss Models *******
c     Idriss (1991), Horizontal, Rock
C     Model Number = 040
      if ( jcalc .eq. 40 ) then
         anper = 24
         minaper = 0.03
         maxaper = 5.0
      endif

C     Idriss (1991), Horizontal, Soft-soil, PGA
C     Model Number = 041
      if ( jcalc .eq. 41 .and. specT .eq. 0.0 ) then
         anper = 1
         minaper = 0.0
         maxaper = 0.0
      endif

c     Idriss 1997 Horizontal, soft-soil, PGA
C     Model Number = 042
      if ( jcalc .eq. 42 .and. specT .eq. 0.0 ) then
         anper = 1
         minaper = 0.0
         maxaper = 0.0
      endif

c     Idriss (1991:1995), Horizontal, Rock
C     Model Number = 043
      if ( jcalc .eq. 43 ) then
         anper = 24
         minaper = 0.03
         maxaper = 5.0
      endif

c     Idriss (1991:1995), Horizontal, Rock scaled by factor 1.0/1.67
C     Model Number = 2043
      if ( jcalc .eq. 2043 ) then
         anper = 24
         minaper = 0.03
         maxaper = 5.0
      endif

c     Idriss (1991:1995), Horizontal, Rock scaled by factor 1.67
C     Model Number = 3043
      if ( jcalc .eq. 3043 ) then
         anper = 24
         minaper = 0.03
         maxaper = 5.0
      endif

c  ******* Sadigh/Geomatrix Models *******
c     Geomatrix 93 (rock) vertical
C     Model Number = 050
      if ( jcalc .eq. 50 ) then
         anper = 22
         minaper = 0.04
         maxaper = 3.0
      endif

c     Sadigh et al. 97 (rock) Horizontal
C     Model Number = 051
      if ( jcalc .eq. 51 ) then
         anper = 22
         minaper = 0.04
         maxaper = 7.5
      endif

c     Sadigh et al. 97 (rock) Horizontal scaled by factor 1.0/1.67
C     Model Number = 2051
      if ( jcalc .eq. 2051 ) then
         anper = 22
         minaper = 0.04
         maxaper = 7.5
      endif

c     Sadigh et al. 97 (rock) Horizontal scaled by factor 1.67
C     Model Number = 3051
      if ( jcalc .eq. 3051 ) then
         anper = 22
         minaper = 0.04
         maxaper = 7.5
      endif

c     Sadigh et al. 97 (soil) horizontal
C     Model Number = 052
      if ( jcalc .eq. 52 ) then
         anper = 13
         minaper = 0.075
         maxaper = 4.0
      endif

c     Sadigh et al. 97 (rock) Horizontal - Sigma = 0.0
C     Model Number = 053
      if ( jcalc .eq. 53 ) then
         anper = 22
         minaper = 0.04
         maxaper = 7.5
      endif

c ******** Spudich et al. (1997) Models *******
C     Spudich et al. (1997), Horizontal, Rock, Extensional Regimes
C     Model Number = 060
      if (jcalc .eq. 60) then
         anper = 11
         minaper = 0.1
         maxaper = 2.0 
      endif

C     Spudich et al. (1997), Horizontal, Rock, Extensional Regimes
C     Model Number = 061
      if (jcalc .eq. 61) then
         anper = 11
         minaper = 0.1
         maxaper = 2.0
      endif

c ******** Youngs Models *******
c     Youngs et al (1993) Horizontal, subduction, Rock
C     Model Number = 200
      if ( jcalc .eq. 200 ) then
         anper = 10
         minaper = 0.067
         maxaper = 3.0
      endif

c     Youngs et al (1997) Horizontal, subduction, rock
C     Model Number = 201
      if ( jcalc .eq. 201 ) then
         anper = 12
         minaper = 0.075
         maxaper = 3.0
      endif

c     Youngs et al (1997) Horizontal, subduction, soil
C     Model Number = 202
      if ( jcalc .eq. 202 ) then
         anper = 13
         minaper = 0.075
         maxaper = 4.0
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
C     Now call the specific Subduction and Crustal Attenuation Model.

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Abrahamson&Silva (1997) Horizontal Crustal
C     Model Number = 7201001
         if (jcalc .eq. 7201001) then
           anper = 12
           minaper = 0.075
           maxaper = 3.0

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Idriss (1991:1995) Horizontal Crustal
C     Model Number = 7201043
         elseif ( jcalc .eq. 7201043 ) then
           anper = 12
           minaper = 0.075
           maxaper = 3.0

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Sadigh et al. (1997) Horizontal Crustal
C     Model Number = 7201051
         elseif ( jcalc .eq. 7201051 ) then
           anper = 12
           minaper = 0.075
           maxaper = 3.0
         endif
      endif

c     The second suite of synchronous ruptures is for seismic source
c          Model B (Carver). Fault parameters for crustal event are:
c          M =7.7, D=1.0 km, Ftype=1 (reverse)

C     Now call the synchronous rupture cases for rock ground motions
c     (i.e., no HBIP amp factors applied). 
c     Set the crustal earthquake parameters
      if (jcalc .eq. 8201001 .or. jcalc .eq. 8201043 .or. jcalc .eq.
     1      8201051) then
C     Now call the specific Subduction and Crustal Attenuation Model.

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Abrahamson&Silva (1997) Horizontal Crustal
C     Model Number = 8201001
         if (jcalc .eq. 8201001) then
           anper = 12
           minaper = 0.075
           maxaper = 3.0

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Idriss (1991:1995) Horizontal Crustal
C     Model Number = 8201043
         elseif ( jcalc .eq. 8201043 ) then
           anper = 12
           minaper = 0.075
           maxaper = 3.0

c     Youngs et al (1997) Horizontal, subduction, Rock
c     with Sadigh et al. (1997) Horizontal Crustal
C     Model Number = 8201051
         elseif ( jcalc .eq. 8201051 ) then
           anper = 12
           minaper = 0.075
           maxaper = 3.0
         endif
      endif
           
c ***** Atkinson and Boore Subduction Models *****
c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 210
      if (jcalc .eq. 210) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 211
      if (jcalc .eq. 211) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 212
      if (jcalc .eq. 212) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 213
      if (jcalc .eq. 213) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 220
      if (jcalc .eq. 220) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 221
      if (jcalc .eq. 221) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 222
      if (jcalc .eq. 222) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 223
      if (jcalc .eq. 223) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 230
      if (jcalc .eq. 230) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 231
      if (jcalc .eq. 231) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 232
      if (jcalc .eq. 232) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 233
      if (jcalc .eq. 233) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c ***** Atkinson and Boore Subduction Models with Errata Correction *****
c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 310
      if (jcalc .eq. 310) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 311
      if (jcalc .eq. 311) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 312
      if (jcalc .eq. 312) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 313
      if (jcalc .eq. 313) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 320
      if (jcalc .eq. 320) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 321
      if (jcalc .eq. 321) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 322
      if (jcalc .eq. 322) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 323
      if (jcalc .eq. 323) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-B, Subduction
C     Model Number = 330
      if (jcalc .eq. 330) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-C, Subduction
C     Model Number = 331
      if (jcalc .eq. 331) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-D, Subduction
C     Model Number = 332
      if (jcalc .eq. 332) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c     Atkinson and Boore (2003) - Horizontal, NEHRP-E, Subduction
C     Model Number = 333
      if (jcalc .eq. 333) then
         anper = 8
         minaper = 0.04
         maxaper = 3.0
      endif

c ***** Gregor et al. Cascadia Subduction Model *****
c     Gregor et al. (2002) - Horizontal, Rock, Cascadia Subduction
C     Model Number = 240
      if (jcalc .eq. 240) then
         anper = 25
         minaper = 0.02
         maxaper = 5.0
      endif

c     Gregor et al. (2002) - Horizontal, Soil, Cascadia Subduction
C     Model Number = 241
      if (jcalc .eq. 241) then
         anper = 25
         minaper = 0.02
         maxaper = 5.0
      endif

c ***** Gregor et al. Cascadia Subduction Model *****
c     Gregor et al. (2006) - Horizontal, Cascadia Subduction
C     Model Number = 242
      if (jcalc .eq. 242) then
         anper = 25
         minaper = 0.02
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Subduction, Hard Rock
C     Model Number = 250
      if ( jcalc .eq. 250 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Subduction, Rock
C     Model Number = 251
      if ( jcalc .eq. 251 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Subduction, Hard Soil
C     Model Number = 252
      if ( jcalc .eq. 252 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Subduction, Medium Soil
C     Model Number = 253
      if ( jcalc .eq. 253 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Subduction, Soft Soil
C     Model Number = 254
      if ( jcalc .eq. 254 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Crustal, Hard Rock
C     Model Number = 255
      if ( jcalc .eq. 255 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Crustal, Rock
C     Model Number = 256
      if ( jcalc .eq. 256 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao-Lu (2011), Horizontal, Crustal, Rock
C     Modified (2006) model with Magnitude>=7.1 capped at M=7.1
C     Model Number = 2256
      if ( jcalc .eq. 2256 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Crustal, Hard Soil
C     Model Number = 257
      if ( jcalc .eq. 257 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Crustal, Medium Soil
C     Model Number = 258
      if ( jcalc .eq. 258 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Zhao et al. (2006), Horizontal, Crustal, Soft Soil
C     Model Number = 259
      if ( jcalc .eq. 259 ) then
         anper = 22
         minaper = 0.01
         maxaper = 5.0
      endif

c     Kanno et al. (2006), Horizontal, Subduction, Vs30m
C     Model Number = 260
      if ( jcalc .eq. 260 ) then
         anper = 38
         minaper = 0.01
         maxaper = 5.0
      endif

c     Garcia et al. (2005), Horizontal, Subduction(intraslab), Rock
C     Model Number = 270
      if ( jcalc .eq. 270 ) then
         anper = 16
         minaper = 0.04
         maxaper = 5.0
      endif

c     Garcia et al. (2005), Vertical, Subduction(intraslab), Rock
C     Model Number = 271
      if ( jcalc .eq. 271 ) then
         anper = 16
         minaper = 0.04
         maxaper = 5.0
      endif

c     Lin and Lee (2008), Horizontal, Subduction, Rock
C     Model Number = 280
      if ( jcalc .eq. 280 ) then
         anper = 28
         minaper = 0.01
         maxaper = 5.0
      endif

c     Lin and Lee (2008), Horizontal, Subduction, Soil
C     Model Number = 281
      if ( jcalc .eq. 281 ) then
         anper = 28
         minaper = 0.01
         maxaper = 5.0
      endif

C  **** BCSubduction Models *******
C     Model Number = 350
      if (jcalc .eq. 350) then
         anper = 22
         minaper = 0.05
         maxaper = 10.0
      endif
C     Model Numbers = 351-365
      if (jcalc .ge. 351 .and. jcalc .le. 365) then
         anper = 23
         minaper = 0.02
         maxaper = 10.0
      endif

C  **** Atkinson&Macias (2009), Cascadia, NEHRP B/C *******
C     Model Number = 370
      if (jcalc .eq. 370) then
         anper = 24
         minaper = 0.05
         maxaper = 10.0
      endif

c ***** Atkinson and Boore Models *****
c     Atkinson and Boore (1994) - Horizontal, EUS Hard Rock
C     Model Number = 100
      if ( jcalc .eq. 100 ) then
         anper = 12
         minaper = 0.05
         maxaper = 2.0
      endif

c     Atkinson and Boore (1994), Horizontal, EUS Hard Rock, magnitude Nuttli
C     Model Number = 101
      if ( jcalc .eq. 101 ) then
         anper = 12
         minaper = 0.05
         maxaper = 2.0
      endif

c     Atkinson and Boore (2006), Horizontal, CEUS Hard Rock
C     Model Number = 102
      if ( jcalc .eq. 102 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Atkinson and Boore (2006), Horizontal, CEUS, Vs=760m/sec
C     Model Number = 103
      if ( jcalc .eq. 103 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Atkinson (2008 wt C0), Horizontal, CEUS, Vs=760m/sec
C     Model Number = 104
      if ( jcalc .eq. 104 ) then
         anper = 7
         minaper = 0.10
         maxaper = 5.0
      endif
      
c     Atkinson and Boore (2006) with Atkinson (2010) magnitude stress drop adjustment, Horizontal, CEUS Hard Rock
C     Model Number = 105
      if ( jcalc .eq. 105 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif
      
c     Atkinson and Boore (2006) with Atkinson (2010) magnitude stress drop adjustment, Horizontal, CEUS, Vs=760m/sec
C     Model Number = 106
      if ( jcalc .eq. 106 ) then
         anper = 26
         minaper = 0.01
         maxaper = 5.0
      endif
      
c     Atkinson (2010), Horizontal, CEUS, based on NGA BA08
C     Model Number = 107
      if (jcalc .eq. 107) then
         anper = 21
         minaper = 0.01
         maxaper = 10.0
      endif

c     Atkinson (2008 avg C0), Horizontal, CEUS, Vs=760m/sec
C     Model Number = 108
      if ( jcalc .eq. 108 ) then
         anper = 7
         minaper = 0.10
         maxaper = 5.0
      endif
      
c     Atkinson and Boore (2006) with (2x) stressDrop (280bars), Horizontal, CEUS Hard Rock
C     Model Number = 130
      if ( jcalc .eq. 130 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Atkinson and Boore (2006) with (0.5x) stressDrop (70bars), Horizontal, CEUS Hard Rock
C     Model Number = 131
      if ( jcalc .eq. 131 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c ******* Toro et al. Models *******
C Toro et al. (1996) MidCon., Horizontal, Rock
C     Model Number = 110
      if (jcalc .eq. 110) then
         anper = 8
         minaper = 0.029
         maxaper = 2.0                                   
      endif 

C Toro et al. (1996) MidCon., Horizontal, Rock, MLg magnitude
C     Model Number = 111
      if (jcalc .eq. 111) then
         anper = 8
         minaper = 0.029
         maxaper = 2.0                                   
      endif 

C Toro et al. (1996) Gulf, Horizontal, Rock
C     Model Number = 112
      if (jcalc .eq. 112) then
         anper = 8
         minaper = 0.029
         maxaper = 2.0                                  
      endif 

C Toro et al. (1996) Gulf, Horizontal, Rock, MLg magnitude
C     Model Number = 113
      if (jcalc .eq. 113) then
         anper = 8
         minaper = 0.029
         maxaper = 2.0                                   
      endif 

c     Campbell Hybrid (2003), Horizontal, CEUS Hard Rock
C     Model Number = 120
      if ( jcalc .eq. 120 ) then
         anper = 18
         minaper = 0.01
         maxaper = 5.0
      endif

c     Campbell Hybrid (2003) - SigmaEps, Horizontal, CEUS Hard Rock
C     Model Number = 121
      if ( jcalc .eq. 121 ) then
         anper = 18
         minaper = 0.01
         maxaper = 5.0
      endif

c     Campbell Hybrid (2003) + SigmaEps, Horizontal, CEUS Hard Rock
C     Model Number = 122
      if ( jcalc .eq. 122 ) then
         anper = 18
         minaper = 0.01
         maxaper = 5.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 2Corner
C     Model Number = 401
      if (jcalc .eq. 401) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 2Corner with Saturation
C     Model Number = 402
      if (jcalc .eq. 402) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Variable-High
C     Model Number = 403
      if (jcalc .eq. 403) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Variable-Med
C     Model Number = 404
      if (jcalc .eq. 404) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Variable-Low
C     Model Number = 405
      if (jcalc .eq. 405) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Const-High
C     Model Number = 406
      if (jcalc .eq. 406) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Const-Med
C     Model Number = 407
      if (jcalc .eq. 407) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Const-Low
C     Model Number = 408
      if (jcalc .eq. 408) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Const-High with Sat.
C     Model Number = 409
      if (jcalc .eq. 409) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Const-Med with Sat.
C     Model Number = 410
      if (jcalc .eq. 410) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Const-Low with Sat.
C     Model Number = 411
      if (jcalc .eq. 411) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 2Corner
C     Model Number = 501
      if (jcalc .eq. 501) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 2Corner with Saturation
C     Model Number = 502
      if (jcalc .eq. 502) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Variable-High
C     Model Number = 503
      if (jcalc .eq. 503) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Variable-Med
C     Model Number = 504
      if (jcalc .eq. 504) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Variable-Low
C     Model Number = 505
      if (jcalc .eq. 505) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Const-High
C     Model Number = 506
      if (jcalc .eq. 506) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Const-Med
C     Model Number = 507
      if (jcalc .eq. 507) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Const-Low
C     Model Number = 508
      if (jcalc .eq. 508) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Const-High with Sat.
C     Model Number = 509
      if (jcalc .eq. 509) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Const-Med with Sat.
C     Model Number = 510
      if (jcalc .eq. 510) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Gulf, 1Corner Const-Low with Sat.
C     Model Number = 511
      if (jcalc .eq. 511) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c *****  Misc Models ******
c     McVerry et al (1993) new zealand
C     Model Number = 300
      if ( jcalc .eq. 300 ) then
         anper = 1
         minaper = 0.0
         maxaper = 0.0
      endif

c     fukushima (1990) rock
C     Model Number = 301
      if ( jcalc .eq. 301 ) then
         anper = 1
         minaper = 0.0
         maxaper = 0.0
      endif

c     Loh high speed rail (New Joyner-Boore form)
C     Model Number = 302      
      if ( jcalc .eq. 302 ) then
         anper = 1
         minaper = 0.0
         maxaper = 0.0
      endif
         
c     New Loh (1996) model (unpublished)
C     Model Number = 303
      if ( jcalc .eq. 303 ) then
         anper = 1
         minaper = 0.0
         maxaper = 0.0
      endif

c ******* Ambraseys et al 2005 Model *********
C     Model Number = 601
      if ( jcalc .eq. 601 ) then
          anper = 62
          minaper = 0.05
          maxaper = 2.50
      endif


C     **** BCHydro SCR attenuation models adjusted for Vs=760m/sec ****
c     Atkinson and Boore (2006), Horizontal, CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1020
      if ( jcalc .eq. 1020 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Atkinson and Boore (2006) with Atkinson (2010) magnitude stress drop adjustment, Horizontal with BC Hydro Amps for Vs=760m/s
C     Model Number = 1050
      if ( jcalc .eq. 1050 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Atkinson and Boore (2006), Horizontal (2x) stressDrop (140), CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1300
      if ( jcalc .eq. 1300 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Atkinson and Boore (2006), Horizontal (0.5x) stressDrop (70), CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1310
      if ( jcalc .eq. 1310 ) then
         anper = 25
         minaper = 0.025
         maxaper = 5.0
      endif

c     Campbell Hybrid (2003), Horizontal, CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1200
      if ( jcalc .eq. 1200 ) then
         anper = 18
         minaper = 0.01
         maxaper = 5.0
      endif

c     Campbell Hybrid (2003) minus Espistemic, Horizontal, CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1210
      if ( jcalc .eq. 1210 ) then
         anper = 18
         minaper = 0.01
         maxaper = 5.0
      endif

c     Campbell Hybrid (2003) plus Espistemic, Horizontal, CEUS with BC Hydro Amps for Vs=760m/s
C     Model Number = 1220
      if ( jcalc .eq. 1220 ) then
         anper = 18
         minaper = 0.01
         maxaper = 5.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Variable-High with BC Hydro Amps for Vs=760m/s
C     Model Number = 4030
      if (jcalc .eq. 4030) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Variable-Med with BC Hydro Amps for Vs=760m/s
C     Model Number = 4040
      if (jcalc .eq. 4040) then
         anper = 27
         minaper = 0.01
         maxaper = 10.0
      endif

c     Silva et al. (2002) - Horizontal, CEUS-Midcont, 1Corner Variable-Low with BC Hydro Amps for Vs=760m/s
C     Model Number = 4050
      if (jcalc .eq. 4050) then
         anper = 26
         minaper = 0.02
         maxaper = 10.0
      endif

C     **** End of BCHydro SCR attenuation models adjusted for Vs=760m/sec ****

c     Akkar and Cagan (2010) - Horizontal, Turkey Crustal
C     Model Number = 150
      if (jcalc .eq. 150) then
         anper = 17
         minaper = 0.01
         maxaper = 2.0
      endif

c     Akkar and Bommer (2010) - Horizontal, Rock (Vs>750m/s)
C     Model Number = 151
      if (jcalc .eq. 151) then
         anper = 63
         minaper = 0.01
         maxaper = 3.0
      endif

c     Akkar and Bommer (2010) - Horizontal, Stiff Soil (360<Vs<750m/s)
C     Model Number = 152
      if (jcalc .eq. 152) then
         anper = 63
         minaper = 0.01
         maxaper = 3.0
      endif

c     Akkar and Bommer (2010) - Horizontal, Soft Soil (Vs<360m/s)
C     Model Number = 153
      if (jcalc .eq. 153) then
         anper = 63
         minaper = 0.01
         maxaper = 3.0
      endif

c     Akkar, Sandikkaya, and Bommer (2013) - Horizontal
C     Model Number = 154
      if (jcalc .eq. 154) then
         anper = 20
         minaper = 0.01
         maxaper = 4.0
      endif

c     Bradley (2010) - Horizontal
C     Model Number = 160
      if (jcalc .eq. 160) then
         anper = 23
         minaper = 0.01
         maxaper = 10.0
      endif

c     McVerry et al. (2006) - Crustal, Horizontal
C     Model Number = 140, Site Class A/B
      if (jcalc .eq. 140) then
         anper = 14
         minaper = 0.01
         maxaper = 3.0
      endif
C     Model Number = 141, Site Class C
      if (jcalc .eq. 141) then
         anper = 14
         minaper = 0.01
         maxaper = 3.0
      endif
C     Model Number = 142, Site Class D
      if (jcalc .eq. 142) then
         anper = 14
         minaper = 0.01
         maxaper = 3.0
      endif

c     McVerry et al. (2006) - Subduction, Horizontal
C     Model Number = 143, Site Class A/B
      if (jcalc .eq. 143) then
         anper = 14
         minaper = 0.01
         maxaper = 3.0
      endif
C     Model Number = 144, Site Class C
      if (jcalc .eq. 144) then
         anper = 14
         minaper = 0.01
         maxaper = 3.0
      endif
C     Model Number = 145, Site Class D
      if (jcalc .eq. 145) then
         anper = 14
         minaper = 0.01
         maxaper = 3.0
      endif      

c     Bindi et al. (2009) - Crustal, Horizontal
C     Model Number = 95, Hor, Rock
      if (jcalc .eq. 95) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif
C     Model Number = 96, Hor, Shallow alluvioum
      if (jcalc .eq. 96) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif
C     Model Number = 97, Hor, Deep alluvium
      if (jcalc .eq. 97) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif

c     Bindi et al. (2011) - Crustal, Horizontal, Class A
C     Model Number = 195, Hor, Class A
      if (jcalc .eq. 195) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif
C     Model Number = 196, Hor, Class B
      if (jcalc .eq. 196) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif
C     Model Number = 197, Hor, Class C
      if (jcalc .eq. 197) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif
C     Model Number = 198, Hor, Class D
      if (jcalc .eq. 198) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif
C     Model Number = 199, Hor, Class E
      if (jcalc .eq. 199) then
         anper = 23
         minaper = 0.01
         maxaper = 2.0
      endif

c     Bindi et al. (2013) - Crustal, Horizontal, Rjb, Vs
C     Model Number = 295, Hor, Class A
      if (jcalc .eq. 295) then
         anper = 26
         minaper = 0.01
         maxaper = 3.0
      endif

c     Grazier and Kalkan (Nov. 2012) Currently Unpublished Update 
c         to Jan/Feb. 2011 SRL model
C     Model Number = 90, Hor
      if (jcalc .eq. 90) then
         anper = 36
         minaper = 0.01
         maxaper = 10.0
      endif

C     Model Number = 91, Hor
      if (jcalc .eq. 91) then
         anper = 36
         minaper = 0.01
         maxaper = 10.0
      endif

C     Model Number = 92, Hor
      if (jcalc .eq. 92) then
         anper = 36
         minaper = 0.01
         maxaper = 10.0
      endif


C     DCPP Common Model Forms
C     ASK model
C     Model Number = 8001, Hor
      if (jcalc .eq. 8001) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     BSSA model
C     Model Number = 8002, Hor
      if (jcalc .eq. 8002) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 001
C     Model Number = 8003, Hor
      if (jcalc .eq. 8003) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 002
C     Model Number = 8004, Hor
      if (jcalc .eq. 8004) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 003
C     Model Number = 8005, Hor
      if (jcalc .eq. 8005) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 004
C     Model Number = 8006, Hor
      if (jcalc .eq. 8006) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 005
C     Model Number = 8007, Hor
      if (jcalc .eq. 8007) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 006
C     Model Number = 8008, Hor
      if (jcalc .eq. 8008) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 007
C     Model Number = 8009, Hor
      if (jcalc .eq. 8009) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 008
C     Model Number = 8010, Hor
      if (jcalc .eq. 8010) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif

C     Common model 009
C     Model Number = 8011, Hor
      if (jcalc .eq. 8011) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif


C     PVNGS Common Model Forms
C     ASK model
C     Model Number = 9001, Hor
      if (jcalc .eq. 9001) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     Bindi model
C     Model Number = 9002, Hor
      if (jcalc .eq. 9002) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     001 model
C     Model Number = 9003, Hor
      if (jcalc .eq. 9003) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     002 model
C     Model Number = 9004, Hor
      if (jcalc .eq. 9004) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     003 model
C     Model Number = 9005, Hor
      if (jcalc .eq. 9005) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     004 model
C     Model Number = 9006, Hor
      if (jcalc .eq. 9006) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     005 model
C     Model Number = 9007, Hor
      if (jcalc .eq. 9007) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     006 model
C     Model Number = 9008, Hor
      if (jcalc .eq. 9008) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     007 model
C     Model Number = 9009, Hor
      if (jcalc .eq. 9009) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     008 model
C     Model Number = 9010, Hor
      if (jcalc .eq. 9010) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif
C     009 model
C     Model Number = 9011, Hor
      if (jcalc .eq. 9011) then
         anper = 5
         minaper = 0.01
         maxaper = 2.0
      endif





c     Check for valid jcalc
      if ( anper .lt. 0 ) then
         write (*,'( 2x,''invalid jcalc:'',i7)') jcalc
         stop 99
      endif

      return
      end
