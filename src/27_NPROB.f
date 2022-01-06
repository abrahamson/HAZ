
      subroutine S27_NDTR( X, P, D)
C          Reference: Abramowitz and Stegan equation 26.2.17
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).
C          D IS DENSITY VALUE (OUTPUT).

      implicit none

      real X, AX, P, D, T

      if (X .lt. 0.) then
        AX = -X
      elseif (X .ge. 0.) then
        AX = X
      endif

      if (AX-6.0 .lt. 0.) then
        T = 1. / (1.0 + 0.2316419 * AX)
        D = 0.3989423 * EXP(-X*X / 2.0)
        P = 1.0 - D*T*( (((1.330274*T - 1.821256)*T + 1.781478) * T-
     1    0.3565638) * T + 0.3193815)
      elseif (AX-6.0 .ge. 0.) then
        P = 1.
        D = 0.
      endif

      if (X .ge. 0) then
        P = 1.0 - P
      endif

      return
      end

c ------------------------------------------------------------------
      subroutine S27_NDTR3( X, P)
C          Reference: Abramowitz and Stegan equation 7.1.26
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).

      implicit none

      real x
      real*8 p, x1, x2, p1, a1, a2, a3, a4, a5, t
      data p1, a1, a2, a3, a4, a5 / 0.3275911, 0.254829592,
     1     -0.284496736, 1.421413741, -1.453152027, 1.061405429 /

      if ( x .lt. 0. ) then
	    x1 = abs(x)
	  else
	    x1 = x
      endif

        x2 = x1/(sqrt(2.))
        t = 1/(1+(p1*x2))
        p = 1-0.5*((a1*t)+(a2*(t**2))+(a3*(t**3))+(a4*(t**4))+(a5*(t**
     1      5)))*(exp(-(x2**2)))

      if ( x .gt. 0. ) then
	    p = 1. - p
      endif

      RETURN
      END
c ------------------------------------------------------------------

      real function pxceed ( eti, ti, siga, jj, i, sigTrunc )

      implicit none
      include 'pfrisk.h'

      integer jj, i
      real eti, ti(MAX_PROB,1), siga, sigTrunc, W, G, D, g1

      W = (ti(JJ,I)-ETI)/SIGA
      call S27_NDTR(W,G,D)

c     TRUNCATE DISTRIBUTION AT SigTrunc
      if (w .gt. sigTrunc) then
        pxceed = 0.0
      else
        call S27_NDTR(sigTrunc,g1,D)
        pxceed = (g-g1)*(1.+g1)
      endif
      return
      end

c ------------------------------------------------------------------

      real*8 function pxceed3 ( eti, ti, siga, jj, i, sigTrunc )

      implicit none
      include 'pfrisk.h'

      integer jj, i
      real eti, ti(MAX_PROB,MAX_INTEN), siga, sigTrunc, w
      real*8 g, g1

c     truncates and renormalizes on both low and high end

      w = (ti(jj,i)-eti)/siga
      call S27_NDTR3(w,g)

      if (w .gt. sigTrunc) then
        pxceed3 = 0.0
      else if (w .lt. (sigTrunc*(-1.))) then
        pxceed3 = 1.0
      else
        call S27_NDTR3(sigTrunc,g1)
        pxceed3 = (g-g1)/(1.-(2.*g1))
      endif

      return
      end

c ------------------------------------------------------------------

      real*8 function pxceed4 ( eti, ti, siga, sigTrunc )

      implicit none
      include 'pfrisk.h'

      real eti, ti, siga, sigTrunc, w
      real*8 g, g1

c     truncates and renormalizes on both low and high end

      w = (ti-eti)/siga
      call S27_NDTR3(w,g)

      if (w .gt. sigTrunc) then
        pxceed4 = 0.0
      else if (w .lt. (sigTrunc*(-1.))) then
        pxceed4 = 1.0
      else
        call S27_NDTR3(sigTrunc,g1)
        pxceed4 = (g-g1)/(1.-(2.*g1))
      endif

      return
      end


c ------------------------------------------------------------------

      subroutine S27_rupArea ( sourceType, mag, coef_area, sigArea,
     1           areastep, sigMaxArea, rupArea, pArea, iFlt, iArea )

      implicit none
      include 'pfrisk.h'

      integer sourceType, iFlt, iArea
      real mag, coef_area(2,MAX_FLT), sigArea(*), rupArea, pArea, nSigma,
     1     nSigma_plus, nSigma_minus, F0, F1, F2, D, areastep, sigMaxArea,
     2     C2, Mo, rigidity

      if (sourceType .eq. 7) then
        pArea = 1.0
      else
        nSigma = -sigMaxArea + (iArea-0.5)*areastep
        nSigma_plus = (nSigma + areastep/2.)
        nSigma_minus = (nSigma - areastep/2.)

c       rupture coefficient is C2 (Leonard)
        if (int(coef_area(1,iflt)) .eq. -999) then
          C2 = coef_area(2,iflt)* 10.**(-5.)
          rigidity = 3.0e11
          Mo = 10.**(16.05 + 1.5*mag)
          rupArea = 10.**(log10(Mo/(C2*rigidity))*(2./3.) + nSigma*sigArea(iflt)) *
     1              (10.**(-10.))
c       rupture coefficients are a and b (Wells and Coppersmith)
        else
          rupArea = 10.0**(coef_area(1,iflt)+coef_area(2,iflt)*mag+nSigma*sigArea(iflt))
        endif

c       Compute probability that (log) rupture dimension is between
c       dim_log_minus and dim_log_plus
        call S27_NDTR ( sigMaxArea, F0, D )
        call S27_NDTR ( nSigma_minus, F1, D )
        call S27_NDTR ( nSigma_plus, F2, D )
        pArea = (F1-F2)/(1-2*f0)
      endif

      return
      end



c ------------------------------------------------------------------

      subroutine S27_rupWidth ( sourceType, mag, rupArea, coef_width,
     1           sigWidth, widthstep, sigMaxWidth, rupWidth, pWidth, iFlt, iWidth)

      implicit none
      include 'pfrisk.h'

      integer sourceType, iFlt, iWidth
      real mag, coef_width(2,MAX_FLT), sigWidth(*), rupWidth, pWidth, nSigma,
     1     nSigma_plus, nSigma_minus, F0, F1, F2, D, widthstep, sigMaxWidth,
     2     rupArea, Beta, C1, AR

      if (sourceType .eq. 7) then
        pWidth = 1.0
      else
        nSigma = -sigMaxWidth + (iWidth-0.5)*widthstep
        nSigma_plus = (nSigma + widthstep/2.)
        nSigma_minus = (nSigma - widthstep/2.)

c       rupture coefficient is C1 (Leonard)
        if (int(coef_width(1,iflt)) .eq. -999) then
          C1 = coef_width(2,iflt)
          Beta = (2./3.)
          rupWidth = ((C1/10.)*(((rupArea/(C1/10.))**(3./5.))**Beta)) *
     1               (10.0**(nSigma*sigWidth(iflt)))
c       rupture coefficient is constant aspect ratio
        elseif (int(coef_width(1,iflt)) .eq. -888) then
          AR = coef_width(2,iflt)
          rupWidth = (rupArea/AR)**(1./2.) * (10.0**(nSigma*sigWidth(iflt)))
c       rupture coefficients are a and b (Wells and Coppersmith)
        else
          rupWidth = 10.0**(coef_width(1,iflt)+coef_width(2,iflt)*mag+nSigma*sigWidth(iflt))
        endif

c       Compute probability that (log) rupture dimension is between
c       dim_log_minus and dim_log_plus
        call S27_NDTR ( sigMaxWidth, F0, D )
        call S27_NDTR ( nSigma_minus, F1, D )
        call S27_NDTR ( nSigma_plus, F2, D )
        pWidth = (F1-F2)/(1-2*f0)
      endif

      return
      end
