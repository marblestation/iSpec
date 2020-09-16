
      subroutine crosscorr (ishift,rxy,xcross,ycross,npx)
c******************************************************************************
c     This routine determines the correlation of array 'x' and array 'y',
c     comparing array elements 'i' in 'x' with elements 'i+ishift';
c     the number of points in each array is assumed to be equal.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      real*4 xcross(100000),ycross(100000)
      real*8 xsum,ysum,xysum,x2sum,y2sum,covar,varx,vary
      real*8 xx,yy


      xsum = 0.
      ysum = 0.
      xysum = 0.
      x2sum = 0.
      y2sum = 0.
      if (ishift .ge. 0) then
         minpt = 1
         maxpt = npx - ishift
      else
         minpt = iabs(ishift) + 1
         maxpt = npx
      endif


      do 10 i=minpt,maxpt
         xx = xcross(i)
         yy = ycross(i+ishift)
         xsum = xsum + xx
         ysum = ysum + yy
         xysum = xysum + xx*yy
         x2sum = x2sum + xx*xx
10       y2sum = y2sum + yy*yy


      xn = real(maxpt - minpt + 1)
      covar = xn*xysum - xsum*ysum
      varx = dsqrt(xn*x2sum - xsum*xsum)
      vary = dsqrt(xn*y2sum - ysum*ysum)
      rxy = sngl(covar/(varx*vary))


      return
      end





