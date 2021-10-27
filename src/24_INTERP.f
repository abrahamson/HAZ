       subroutine S24_interp (x1,x2,y1,y2,x,y,iflag)

C This subroutine will perform the Log-linear interpolation
C of the given input values. This routine is used to interpolate
C the regression cofficients of the attenuation models for
C spectral periods other than those defined in the model.

       implicit none

       integer iflag
       real x1, x2, y1, y2, x, y

C Check to see if the interpolation period is at an end point.
C Return the 'iflag' for output purposes with
C             iflag = 0  No interpolation
C                   = 1  Interpolation need.

       if (x .eq. x1 .or. x .eq. x2) then
          iflag = 0
       else
          iflag = 1
       endif

C Set the PGA period to 100 Hz (i.e., 0.01 Sec).
       if (x1 .eq. 0.0) then
           x1 = 0.01
       endif

C Take the Log of the Period values.
       x1 = alog(x1)
       x2 = alog(x2)
       x  = alog(x)
C Perform the log-linear interpolation.
       y = y1 + (y2-y1)*((x-x1)/(x2-x1))

C Convert the Log Periods back to period.
       x1 = exp(x1)
       x2 = exp(x2)
       x  = exp(x)

       return
       end

c-------------------------------------------------------------------------------

       subroutine S24_interp1 (x1,x2,y1,y2,x,y,iflag)

c This subroutine is the same as S24_interp but was modified for
c use with frequencies (NGA-East) rather than periods - the only
c change being that the check for PGA (T=0) and setting
c x1 = 0.01 s is not needed and was removed.

       implicit none

       integer iflag
       real x1, x2, y1, y2, x, y

c Check to see if the interpolation frequency is at an end point.
c Return the 'iflag' for output purposes with
c             iflag = 0  No interpolation
c                   = 1  Interpolation need

       if ( x .eq. x1 .or. x .eq. x2) then
          iflag = 0
       else
          iflag = 1
       endif

c Take the Log of the Frequency values
       x1 = alog(x1)
       x2 = alog(x2)
       x  = alog(x)
c Perform the log-linear interpolation
       y = y1 + (y2-y1)*((x-x1)/(x2-x1))

c Convert the Log Frequencies back to frequency.
       x1 = exp(x1)
       x2 = exp(x2)
       x  = exp(x)

       return
       end
