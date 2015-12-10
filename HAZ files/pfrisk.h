c     Set array dimensions

      integer MAX_SITE, MAX_FLT, MAX_SEG, MAX_INTEN, MAX_PROB, MAX_EPS,
     1        MAX_GM, MAX_Xcost, MAXPARAM, MAX_MAG, MAX_DIST, MAX_N1,
     2        MAX_WIDTH, MAX_DIST1, MAX_GRID, MAX_SYN, MAX_AMPMAG,
     3        MAX_AMPPER, MAX_AMPGM, MAX_PER, MAXDETM_DIST, MAXDIPSEG,
     4        MAX_DD, MAXFLT_DD, MAXFLT_AS, MAX_BRANCH, MAX_NODE,
     5        MAX_ATTEN, MAX_FTYPE


      PARAMETER ( MAX_SITE=1, MAX_FLT=60, MAX_SEG=300,
     1            MAX_INTEN=18,MAX_PROB=3, MAX_EPS=10, MAX_GM=20,
     2            MAX_Xcost=10, MAXPARAM=300, MAX_MAG=30, 
     3            MAX_DIST=25,MAX_N1=220, MAX_WIDTH=15, 
     4            MAX_DIST1=10000, MAX_GRID=32000, MAX_SYN=5)
      Parameter ( MAX_AMPMAG=25, MAX_AMPPER=15, MAX_AMPGM=15)
      parameter (MAX_PER=501)
      parameter (MAXDETM_DIST=2000)
      parameter (MAXDIPSEG=5,MAX_DD=12, MAXFLT_DD=2000, MAXFLT_AS=2000)
      parameter (MAX_BRANCH=30, MAX_NODE=100)
      parameter (MAX_ATTEN=60)
      parameter (MAX_FTYPE=3)
      
