C.....Subroutine for calling Deterministic Ground Motion Calculations.....
C       attenflag = 1 Attenuation curves
C       attenflag = 2 Spectra

      subroutine S15_determ (attenflag)

      implicit none
      include 'pfrisk.h'

      integer cfmodel, icalc, ncalc, jcalc1, ndist, hwflag,
     1        iAtten, jtype, scalc1, sigflag, anper, vs30_class,
     2        attenflag, intflag, jj, forearc, coefcountRrup,
     3        coefcountRjb, idist, iProb, nper, iper, icf, dcpp_GMPE_Flag
      real Theta_site, cfcoefRRup(MAX_Atten,11), cfcoefRjb(MAX_Atten,11),
     1     coefrr(MAX_PER,11), coefrj(MAX_PER,11), mag, distrup, distjb,
     2     distseismo, vs, sigfix1, lgInten, sigmaY, depthvs10, depthvs15,
     3     D25, minaper, maxaper, aper(MAX_PER), ftype(MAX_FLT,MAX_FTYPE),
     4     fTop(MAX_FLT, MAX_SEG), dipd, AR, tau, phi, hypodepth1
      real attenrupdist(MAXDETM_DIST), attenjbdist(MAXDETM_DIST),
     1     attenseisdist(MAXDETM_DIST), attenhypodist(MAXDETM_DIST),
     2     attenRx(MAXDETM_DIST), attenRy0(MAXDETM_DIST), hypodepth,
     3     specT(4,MAX_PROB), gfac1(MAX_PER), gfac2(MAX_PER), svad(MAX_PER),
     4     period1(4,MAX_PROB), disthypo, Rx, Ry0, temp, rupWidth
      character*80 attentitle, attenoutfile, attenName, sigmaName

C     If requested compute ground motion attenuation from input data file.
c        (i.e., attenflag = 1)
      if (attenflag .eq. 1) then
c     read input parameters for attenuation models.
         read (13, '(a80)' ) attentitle
         read (13, '(a80)' ) attenoutfile
         open (67, file=attenoutfile)
         write (67, '(a80)') attentitle
         read (13,*) ncalc

         do icalc=1, ncalc
            read (13,*) jcalc1, gfac1(1), gfac2(1), svad(1)

C     Check for different sigma models either fixed or GMPE.
C     Options are as follows:
C       jcalc > 0 (use GMPE sigma model)
C       jcalc < 0 (read in scalc and sigfix values)
C           if scalc > 0 (use scalc sigma model) sigflag = 1
C           if scalc < 0 (use fixed sigma value) sigflag = 2
            sigflag = 0
            if (jcalc1 .lt. 0) then
               backspace (13)
               read (13,*) jcalc1, gfac1(1), gfac2(1), svad(1), scalc1, sigfix1
               jcalc1 = abs(jcalc1)
               if (scalc1 .lt. 0) then
                  sigflag = 2
               else
                  sigflag = 1
               endif
            endif

C     Check for Common Functional Form with Rrup Distance (10000<jcalc<11000) selected and if so read in coefficients.
            if (jcalc1 .gt. 10000 .and. jcalc1 .lt. 11000) then
                coefcountrrup = jcalc1 - 10000
                read (13,*) (cfcoefrrup(coefcountRrup,jj),jj=1,11)
            endif
C     Check for Common Functional Form with RJB Distance (11000<jcalc<12000) selected and if so read in coefficients.
            if (jcalc1 .gt. 11000 .and. jcalc1 .lt. 12000) then
                coefcountrjb = jcalc1 - 11000
                read (13,*) (cfcoefrjb(coefcountrjb,jj),jj=1,11)
            endif

            read (13,*) mag
            read (13,*) ndist

            do idist=1,ndist
               read (13,*) attenrupdist(idist), attenjbdist(idist),
     1                     attenseisdist(idist), attenhypodist(idist),
     2                     attenRx(idist), attenRy0(idist)
            enddo
            read (13,*) ftype(1,1), dipd, hwflag

c     Check for dipping fault type with Dip angle = 90.
            if (ftype(1,1) .ne. 0.0 .and .dipd .eq. 90.0) then
               write (*,*) 'Fault Type Not Strike-Slip and Dip of '
               write (*,*) '      fault is 90 degrees!'
               write (*,*) '  *** Check input parameters ***'
            endif
            read (13,*) vs, vs30_class, forearc
            read (13,*) hypodepth, RupWidth
            read (13,*) depthvs10, depthvs15, D25, ftop(1,1)
            read (13,*) specT(1,1)

c     Write input parameters to output file.
            write (67, *) ' ***** Input Parameters *****'
            write (67,'(2x,a15,i8)')    'Jcalc      = ', jcalc1
            write (67,'(2x,a15,f8.3)')  'Constant1  = ', gfac1(1)
            write (67,'(2x,a15,f8.3)')  'Constant2  = ', gfac2(1)
            write (67,'(2x,a15,f8.3)')  'SigVarAdd  = ', svad(1)
            write (67,'(2x,a15,f8.3)')  'Mag        = ', mag
            write (67,'(2x,a15,f8.2)')  'Ftype      = ', ftype(1,1)
            write (67,'(2x,a15,f8.3)')  'Dip        = ', dipd
            write (67,'(2x,a15,i8)')    'HWFlag     = ', hwflag
            write (67,'(2x,a15,f8.3)')  'Vs30m      = ', vs
            if (vs30_class .eq. 1) then
               write (67,'(5x,a20)')       'Vs30m -- > measured '
            elseif (vs30_class .eq. 0) then
               write (67,'(5x,a20)')       'Vs30m -- > estimated'
            endif
            write (67,'(2x,a15,f8.3)')  'Hypodepth  = ', hypodepth
            write (67,'(2x,a15,f8.3)')  'Width      = ', RupWidth
            write (67,'(2x,a15,f8.3)')  'DepthVs10  = ', depthvs10
            write (67,'(2x,a15,f8.3)')  'DepthVs15  = ', depthvs15
            write (67,'(2x,a15,f8.3)')  'DepthVs25  = ', D25
            write (67,'(2x,a15,f8.3)')  'Depthtop   = ', ftop(1,1)
            if (forearc .eq. 0) then
               write (67,'(4x,a16)')     'Forearc Site'
            elseif (forearc .eq. 1) then
               write (67,'(4x,a16)')     'Backarc Site    '
            endif
            write (67, '(1x,a38,a110,10x,a80)') ' Period(s)   Mag    Rrup      Rjb     ',
     1            ' SeisDist HypoDist    Rx         Ry0          SA(g)         Sigma   Const1  Const2  SigAdd      Phi       Tau ',
     2            ' GMModel                                               SigmaModel'

C     Perform loop over distances for attenuation models.
            do idist = 1, ndist
               distrup = attenrupdist(idist)
               distJB = attenjbdist(idist)
               distSeismo = attenseisdist(idist)
               disthypo = attenhypodist(idist)
               Rx = attenRx(idist)
               Ry0 = attenRy0(idist)
               iAtten = 1
               jType = 1
               iProb = 1

               call S11_meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, jcalc1, specT(1,1),
     2               lgInten,sigmaY, ftype(1,1), attenName, period1,
     3               iAtten, iProb, jType, vs, hypoDepth,intflag, AR, dipd,
     4               disthypo, depthvs10, depthvs15, D25, tau, ftop(1,1),
     5               theta_Site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefRrup, cfcoefRjb, Ry0)

C             Adjust median ground motion by constant factors
              lgInten = lgInten + gfac1(1) + gfac2(1)

C     Check for sigma values different than requested GMPE.
               if (sigflag .eq. 1) then
                 call S11_meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, scalc1, specT(1,1),
     2               temp, sigmaY, ftype(1,1), sigmaName, period1,
     3               iAtten, iProb, jType, vs, hypoDepth,intflag,AR,dipd,
     4               disthypo, depthvs10, depthvs15, D25, tau,  ftop(1,1),
     5               theta_Site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefRrup, cfcoefRjb, Ry0 )

                 elseif (sigflag .eq. 2) then
                     sigmaY = sigfix1
                     sigmaname = 'Fixed Constant Sigma'
                 else
                     sigmaname = attenname
                 endif

C               Adjust sigma value is sigvaradd .ne. 0.0
                if (svad(1) .lt. 0.0) then
c               Check for reduction of sigma is greater than 0.0.
                   if ((SigmaY+svad(1)) .le. 0.0 ) then
                      sigmaY = 0.0
                   else
                      sigmaY = sqrt(sigmaY*sigmaY + svad(1) )
                   endif
                elseif (svad(1) .gt. 0.0) then
                   sigmaY = sqrt(sigmaY*sigmaY + svad(1) )
                endif

                write (67, 676) specT(1,1), mag, distrup, distJB,
     1                         distseismo, disthypo, Rx, Ry0, exp(lgInten), sigmaY,
     2                         gfac1(1), gfac2(1), svad(1), phi, tau,
     4                         attenname, sigmaname

            enddo
            write (67,*)
         enddo
         close (67)
         write (*,*) 'End of Attenuation Ground Motion Modeling.'

      endif
 676  format (f10.3,f8.3,6f10.3,e12.4,6f8.4,4x,a80,2x,a80)

C     If requested compute ground motion spectra from input data file.
c        (i.e., attenflag = 2)
      if (attenflag .eq. 2) then
c       read input parameters for attenuation models.
        read (13, '(a80)' ) attentitle
        read (13, '(a80)' ) attenoutfile
        open (67, file=attenoutfile)
        write (67, '(a80)') attentitle
        read (13,*) ncalc
        do icalc=1, ncalc
          read (13,*) jcalc1
          write (*,'( i5)') jcalc1

C         Check for different sigma models either fixed or GMPE.
          sigflag = 0
          if (jcalc1 .lt. 0) then
            backspace (13)
            read (13,*) jcalc1, scalc1, sigfix1
            jcalc1 = abs(jcalc1)
            if (scalc1 .lt. 0) then
              sigflag = 2
            else
              sigflag = 1
            endif
          endif

          read (13,*) mag
          read (13,*) ndist
          do idist=1,ndist
            read (13,*) attenrupdist(idist), attenjbdist(idist),
     1                     attenseisdist(idist), attenhypodist(idist),
     2                     attenRx(idist), attenRy0(idist)
          enddo
          read (13,*) ftype(1,1), dipd, hwflag

c         Check for dipping fault type with Dip angle = 90.
          if (ftype(1,1) .ne. 0.0 .and .dipd .eq. 90.0) then
             write (*,*) 'Fault Type Not Strike-Slip and Dip of '
             write (*,*) '      fault is 90 degrees!'
             write (*,*) '  *** Check input parameters ***'
          endif
          read (13,*) vs, vs30_class, forearc
          read (13,*) hypodepth, RupWidth
          read (13,*) depthvs10, depthvs15, D25, ftop(1,1)

C         Get the number of spectral periods for each attenuation relationship model
          read (13,*) nper
          do iper=1,nper
            read (13,*) aper(iper), gfac1(iper), gfac2(iper), svad(iper)

C           Check for Common Functional Form with Rrup Distance (10000<jcalc<11000) selected and if so read in coefficients.
            if (jcalc1 .gt. 10000 .and. jcalc1 .lt. 11000) then
              read (13,*) (coefrr(iper,jj),jj=1,11)
            endif

C           Check for Common Functional Form with RJB Distance (11000<jcalc<12000) selected and if so read in coefficients.
            if (jcalc1 .gt. 11000 .and. jcalc1 .lt. 12000) then
              read (13,*) (coefrj(iper,jj),jj=1,11)
            endif
          enddo

          call S12_attenper (jcalc1, anper, minaper, maxaper )

c         Write input parameters to output file.
          write (67, *) ' ***** Input Parameters *****'
          write (67,'(2x,a15,i8)')    'Jcalc      = ', jcalc1
          write (67,'(2x,a15,f8.3)')  'Mag        = ', mag
          write (67,'(2x,a15,f8.2)')  'Ftype      = ', ftype(1,1)
          write (67,'(2x,a15,f8.3)')  'Dip        = ', dipd
          write (67,'(2x,a15,i8)')    'HWFlag     = ', hwflag
          write (67,'(2x,a15,f8.3)')  'Vs30m      = ', vs
          if (vs30_class .eq. 1) then
            write (67,'(5x,a20)')       'Vs30m -- > measured '
          elseif (vs30_class .eq. 0) then
            write (67,'(5x,a20)')       'Vs30m -- > estimated'
          endif
          write (67,'(2x,a15,f8.3)')  'Hypodepth  = ', hypodepth
          write (67,'(2x,a15,f8.3)')  'Width      = ', RupWidth
          write (67,'(2x,a15,f8.3)')  'DepthVs10  = ', depthvs10
          write (67,'(2x,a15,f8.3)')  'DepthVs15  = ', depthvs15
          write (67,'(2x,a15,f8.3)')  'DepthVs25  = ', D25
          write (67,'(2x,a15,f8.3)')  'Depthtop   = ', ftop(1,1)
          if (forearc .eq. 0) then
            write (67,'(4x,a16)')     'Forearc Site'
          elseif (forearc .eq. 1) then
            write (67,'(4x,a16)')     'Backarc Site    '
          endif
          write (67, '(1x,a38,a110,2x,a80)') ' Period(s)   Mag    Rrup      Rjb     ',
     1             ' SeisDist HypoDist    Rx         Ry0          SA(g)         Sigma   Const1  Const2  SigAdd      Phi       Tau ',
     2             ' Per_Notes                     GMModel                     Sigma Model'

C         Perform loop over distances for attenuation models.
          do idist = 1, ndist
            distrup = attenrupdist(idist)
            distJB = attenjbdist(idist)
            distSeismo = attenseisdist(idist)
            disthypo = attenhypodist(idist)
            Rx = attenRx(idist)
            Ry0 = attenRy0(idist)
            iAtten = 1
            jType = 1

            do iper = 1, nper
              iprob = 1
              jtype = 1
              specT(1,1) = aper(iper)

              dcpp_GMPE_flag = 0

C             Set up coefficients for SWUS Common Functional Form.
              if (jcalc1. gt. 10000 .and. jcalc1 .lt. 11000) then
                cfmodel = jcalc1 - 10000
                dcpp_GMPE_flag = 1
              elseif (jcalc1 .gt. 11000 .and. jcalc1 .lt. 12000) then
                cfmodel = jcalc1 - 11000
                dcpp_GMPE_flag = 1
              endif
              if ( dcpp_GMPE_Flag .eq. 1) then
                do icf=1,11
                  cfcoefRrup(cfmodel,icf) = coefrr(iper,icf)
                  cfcoefRjb(cfmodel,icf) = coefrj(iper,icf)
                enddo
              endif

c             Check for valid period range between minaper and maxaper for the
c             requested attenuation model and or PGA value. Only call 'meanInten'
c             for valid spectral periods and skip cases in which spectral period
c             is not defined (output will be set to null value of -9.99999)
              if (specT(1,1) .ge. minaper .and. specT(1,1) .le. maxaper
     1                 .or. specT(1,1) .eq. 0.0 ) then

                call S11_meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, jcalc1, specT(1,1),
     2               lgInten, sigmaY, ftype(1,1), attenName,
     2               period1, iAtten, iProb, jType, vs, hypoDepth, intflag, AR, dipd,
     4               disthypo, depthvs10, depthvs15, D25, tau, ftop(1,1),
     5               theta_Site, RupWidth, vs30_class, foreArc, Rx, phi,
     6               cfcoefRrup, cfcoefRjb, Ry0 )

              elseif (specT(1,1) .eq. 0.0) then

                call S11_meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, jcalc1, specT(1,1),
     2               lgInten,sigmaY, ftype(1,1), attenName,
     2               period1, iAtten, iProb, jType, vs, hypoDepth, intflag, AR, dipd,
     4               hypodepth1, depthvs10, depthvs15, D25, tau,
     3               ftop(1,1), theta_Site, RupWidth, vs30_class, forearc, Rx, phi,
     4               cfcoefRrup, cfcoefRjb, Ry0 )

              elseif (specT(1,1) .eq. -1.0) then

                call S11_meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, jcalc1, specT(1,1),
     2               lgInten,sigmaY, ftype(1,1), attenName,
     2               period1, iAtten, iProb, jType, vs, hypoDepth, intflag, AR, dipd,
     4               hypodepth1, depthvs10, depthvs15, D25, tau,
     3               ftop(1,1), theta_Site, RupWidth, vs30_class, forearc, Rx, phi,
     4               cfcoefRrup, cfcoefRjb, Ry0 )

              elseif (specT(1,1) .eq. -2.0) then

                call S11_meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, jcalc1, specT(1,1),
     2               lgInten,sigmaY, ftype(1,1), attenName,
     2               period1, iAtten, iProb, jType, vs, hypoDepth, intflag, AR, dipd,
     4               hypodepth1, depthvs10, depthvs15, D25, tau,
     3               ftop(1,1), theta_Site, RupWidth, vs30_class, forearc, Rx, phi,
     4               cfcoefRrup, cfcoefRjb, Ry0 )

               else
                 intflag = -1
               endif

C              Adjust the median ground motions by const1 and const2
               lginten = lginten + gfac1(iper) + gfac2(iper)

C              Now compute the sigma if requested is different than jcalc GMPE.
               if (intflag .ge. 0) then
                 if (sigflag .eq. 1) then
                   call S11_meanInten ( distRup, distJB, distSeismo,
     1                     hwflag, mag, scalc1, specT(1,1),
     2                     temp,sigmaY, ftype(1,1), sigmaName,
     2                     period1, iAtten, iProb, jType, vs, hypoDepth,intflag, AR, dipd,
     4                     disthypo, depthvs10, depthvs15, D25, tau, ftop(1,1),
     5                     theta_Site, RupWidth, vs30_class, foreArc, Rx, phi,
     6                     cfcoefRrup, cfcoefRjb, Ry0 )
                 elseif (sigflag .eq. 2) then
                   sigmaY = sigfix1
                   sigmaname = 'Fixed Constant Sigma'
                 else
                   sigmaname = attenname
                 endif
               endif

C              Adjust sigma value is sigvaradd .ne. 0.0
               if (svad(iPer) .lt. 0.0) then

c                Check for reduction of sigma is greater than 0.0.
                 if (SigmaY**2+(svad(iPer)) .le. 0.0) then
                   sigmaY = 0.01
                 else
                   sigmaY = sqrt(sigmaY*sigmaY + svad(iPer) )
                 endif
               elseif (svad(iper) .gt. 0.0) then
                 sigmaY = sqrt(sigmaY*sigmaY + svad(iper) )
               endif

C              Write out spectra with notes for each spectral period.
               if (intflag .eq. 0 ) then
                 write (67, 678) specT(1,1), mag, distrup, distJB,
     1                            distseismo, disthypo, Rx, Ry0, exp(lgInten), sigmaY,
     4                            gfac1(iper), gfac2(iper), svad(iper), phi, tau,
     4                            'Defined', attenname, sigmaname
               elseif (intflag .eq. 1 ) then
                   write (67, 678) specT(1,1), mag, distrup, distJB,
     1                            distseismo, disthypo, Rx, Ry0, exp(lgInten), sigmaY,
     4                            gfac1(iper), gfac2(iper), svad(iper), phi, tau,
     4                            'Interp', attenname, sigmaname
               elseif (intflag .eq. -1 ) then
                   write (67, 678) specT(1,1), mag, distrup, distJB,
     1                            distseismo, disthypo, Rx, Ry0, -9.99999, -9.9999,
     4                            gfac1(iper), gfac2(iper), svad(iper), phi, tau,
     4                            'Outside', attenname, sigmaname
               endif
             enddo
             write (67,*)
           enddo
         enddo
         close (67)
         write (*,*) 'End of Spectra Ground Motion Modeling.'

      endif
 678  format (f10.3,f8.3,6f10.3,e12.4,6f8.4,4x,a7,3x,a80,2x,a80)

      return
      end
