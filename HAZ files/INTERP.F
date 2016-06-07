       subroutine interp (x1,x2,y1,y2,x,y,iflag)

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
