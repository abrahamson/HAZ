       subroutine S28_RupDims (sourceType, rupWidth, aveWidth, rupArea, faultLen,
     1                     faultWidth, dip, nLocYST1, ystep, rupLen)

       implicit none

c      declarations passed in
       integer sourceType
       real aveWidth, rupArea, faultLen, faultWidth, ystep, dip

c      declarations passed in and out
       real rupWidth

c      declarations passed out
       integer nLocYST1
       real rupLen

c      declarations only used within subroutine
       real availddwidth

c           Check that rupture width doesn't exceed fault width and rupture area doesn't
c           exceed fault area. Set rupture length
            if (sourceType .eq. 1 ) then
              if (rupWidth .gt. aveWidth ) then
                rupWidth = aveWidth
              endif
              rupLen = rupArea / rupWidth
              if (rupLen .gt. faultLen) then
                rupLen = faultLen
              endif
              nLocYST1 = (aveWidth - rupWidth)/yStep + 1
            elseif (sourceType .eq. 5 .or. sourceType .eq. 6) then
              rupLen = rupArea / rupWidth
              if (rupLen .gt. faultLen) then
                rupLen = faultLen
              endif
            elseif (sourceType .eq. 7) then
              rupWidth = 12.0

c           Otherwise source is area
            else
              availddwidth = faultwidth/sin(abs(dip)*3.1415926/180.)
              if (rupWidth .gt. availddwidth) then
                rupWidth = availddwidth
              endif
              rupLen = rupArea / rupWidth
            endif

        return
       end


c ----------------------------------------------------------------------
       subroutine S28_nLocXcells (sourceType, nLocXAS, grid_n, nfltgrid, fltgrid_w,
     1                        rupWidth, fltgrid_a, ruparea, nLocYST1, nLocX, n1AS, n2AS)

       implicit none
       include 'pfrisk.h'

c      declarations passed in
       integer sourceType, nLocXAS, grid_n, nfltgrid(2), nLocYST1
       real fltgrid_w(MAXFLT_DD,MAXFLT_AS), rupWidth, fltgrid_a(MAXFLT_DD,MAXFLT_AS),
     1      ruparea

c      declarations passed out
       integer nLocX, n1AS(MAXFLT_AS), n2AS(MAXFLT_AS)

c      declarations only this subroutine
       integer countnLocX, in2m, in2, iw, in1
       real colarea, colarew

c      Set nLocX for grid and area sources
       if (sourceType .eq. 2 .or. sourceType .eq. 3) then
         nLocX = nLocXAS
       else if (sourceType .eq. 4) then
         nLocX = grid_n
       else if (sourcetype .eq. 7) then
         nLocX = 1

c      Set nLocX, n1AS, and n2AS for fault sources
       else

c         Need to do loop along strike adding up area for all cells with width
c         less than RupWidth

          countnLocX = 0
          do 1234 in2m=1,nfltgrid(2)
             colarea = 0.0

C            Check for case in which total area is less than RupArea.
             if (in2m .gt. 1 .and. countnLocX .eq. 0) then
                goto 1233
             endif
             do 1235 in2=in2m,nfltgrid(2)
                colarew = 0.0
                if (in2 .eq. in2m) then
                   do iw=1,nfltgrid(1)
                      colarew = colarew + fltgrid_w(iw,in2)

C                     Find the number of downdip cells for this along strike location to
C                     get the correct rupture width.
                      if (colarew .gt. rupWidth) then
                         colarew = 0.0
                         goto 2341
                      endif
                   enddo
 2341              continue
                   if (iw .gt. nfltgrid(1) ) iw = nfltgrid(1)
                endif
                do 1236 in1=1,iw
                   colarea = colarea + fltgrid_a(in1,in2)

 1236           continue
                   if (colarea .ge. ruparea) then
                      countnLocX = countnLocX + 1
                      n1AS(in2m) = nfltgrid(1) - iw + 1
                      n2AS(in2m) = in2

                      goto 1234
                   elseif (in2.eq.nfltgrid(2) .and. colarea.lt.ruparea) then
                      goto 1233
                   endif

 1235        continue
 1234     continue

 1233    continue
            if (countnLocX .eq. 0) then
               nLocX = 1
               if (sourceType .eq. 1) then
                 n1AS(1) = nLocYST1
               else
                 n1AS(1) = 1
               endif
               n2AS(1) = nfltgrid(2)
            else
               nLocX = countnLocX
            endif

       endif

       return
       end

c ----------------------------------------------------------------------
       subroutine S28_nLocYcells (iLocX, n1AS, sourceType, nLocX, distDensity, xStep,
     1                        faultWidth, ystep, distDensity2, grid_x, grid_y, x0, y0,
     2                        nLocY, pLocX, r_horiz, VarXstepFlag, VarYstepFlag)

       implicit none
       include 'pfrisk.h'

c      declarations passed in
       integer iLocX, n1AS(MAXFLT_AS), sourceType, nLocX, VarXstepFlag, VarYstepFlag
       real distDensity(MAX_DIST1), xStep, faultWidth, ystep, distDensity2(MAX_GRID),
     1      grid_x(MAX_GRID), grid_y(MAX_GRID), x0, y0

c      declarations passed out
       integer nLocY
       real pLocX, r_horiz

              if (sourceType .eq. 1 .or. sourceType .eq. 5 .or. sourceType .eq. 6) then
                nLocY = n1AS(iLocX)
                pLocX = 1./nLocX
              elseif ( sourceType .eq. 2) then
                pLocX = distDensity(iLocX)
                if ( pLocX .ne. 0. ) then
c                 variable horizontal discretization
                  if (VarXstepFlag .eq. 1) then
                    if (iLocX .lt. 11) then
                      xStep = 1
                      r_horiz = xStep * (iLocX-0.5)
                    else if (iLocX .lt. 19) then
                      xStep = 2
                      r_horiz = 10. + xStep * (iLocX-10.5)
                    else if (iLocX .lt. 30) then
                      xStep = 3
                      r_horiz = 26. + xStep * (iLocX-18.5)
                    else if (iLocX .lt. 53) then
                      xStep = 4
                      r_horiz = 59. + xStep * (iLocX-29.5)
                    else
                      xStep = 5
                      r_horiz = 151. + xStep * (iLocX-52.5)
                    endif
                  else
                    r_horiz = xStep * (iLocX-0.5)
                  endif
c                 variable depth discretization
                  if (VarYstepFlag .eq. 1) then
                    if (r_horiz .lt. 10.) then
                      yStep = 1.0
                    else if (r_horiz .lt. 25.) then
                      yStep = 2.0
                    else if (r_horiz .lt. 60.) then
                      yStep = 3.0
                    else if (r_horiz .lt. 150.) then
                      yStep = 4.0
                    else
                      yStep = 5.0
                    endif
                  endif
                  nLocY = nint(faultWidth / yStep)
                  if (nLocY.eq.0) then
                    nLocY = 1
                  endif
c                 adjust yStep to fit evenly into faultWidth if necessary
                  yStep = faultWidth / nLocY
                endif
              elseif ( sourceType .eq. 3 ) then
                pLocX = distDensity(iLocX)
                if ( pLocX .ne. 0. ) then
c                 variable horizontal discretization
                  if (VarXstepFlag .eq. 1) then
                    if (iLocX .lt. 11) then
                      xStep = 1
                      r_horiz = xStep * (iLocX-0.5)
                    else if (iLocX .lt. 19) then
                      xStep = 2
                      r_horiz = 10. + xStep * (iLocX-10.5)
                    else if (iLocX .lt. 30) then
                      xStep = 3
                      r_horiz = 26. + xStep * (iLocX-18.5)
                    else if (iLocX .lt. 53) then
                      xStep = 4
                      r_horiz = 59. + xStep * (iLocX-29.5)
                    else
                      xStep = 5
                      r_horiz = 151. + xStep * (iLocX-52.5)
                    endif
                  else
                    r_horiz = xStep * (iLocX-0.5)
                  endif
c                 variable depth discretization
                  if (VarYstepFlag .eq. 1) then
                    if (r_horiz .lt. 10.) then
                      yStep = 1.0
                    else if (r_horiz .lt. 25.) then
                      yStep = 2.0
                    else if (r_horiz .lt. 60.) then
                      yStep = 3.0
                    else if (r_horiz .lt. 150.) then
                      yStep = 4.0
                    else
                      yStep = 5.0
                    endif
                  endif
                  nLocY = nint(faultWidth / yStep)
                  if (nLocY.eq.0) then
                    nLocY = 1
                  endif
c                 adjust yStep to fit evenly into faultWidth if necessary
                  yStep = faultWidth / nLocY
                endif
              elseif ( sourceType .eq. 4 ) then
                pLocX = distDensity2(iLocX)
                if ( pLocX .ne. 0. ) then
                  r_horiz = sqrt( (grid_x(iLocX)-x0)**2 + (grid_y(iLocX)-y0)**2 )
c                 variable depth discretization
                  if (VarYstepFlag .eq. 1) then
                    if (r_horiz .lt. 10.) then
                      yStep = 1.0
                    else if (r_horiz .lt. 25.) then
                      yStep = 2.0
                    else if (r_horiz .lt. 60.) then
                      yStep = 3.0
                    else if (r_horiz .lt. 150.) then
                      yStep = 4.0
                    else
                      yStep = 5.0
                    endif
                  endif
                  nLocY = nint(faultWidth / yStep)
                  if (nLocY.eq.0) then
                    nLocY = 1
                  endif
c                 adjust yStep to fit evenly into faultWidth if necessary
                  yStep = faultWidth / nLocY
                endif
              elseif (sourceType .eq. 7) then
                nLocY = 1
                pLocX = 1.0
              endif

       return
       end
