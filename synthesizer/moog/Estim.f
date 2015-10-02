      subroutine estim (xleft,right,a,wl,npts,ileft,iright)
c******************************************************************************
c     This routine returns channel (pixel) values
c******************************************************************************
 
      real*4 wl(1)
      real*8 a(9)

 
      if (a(9).eq.1 .or. a(9).eq.2 .or. a(9).eq.3) go to 40
      if (a(4).eq.0.0 .and. a(3).eq.0.0) go to 60
  

c*****do the estimation for quadratics and cubics
      ileft = int(sngl((-a(2) + dsqrt(a(2)**2-4.*a(3)*(a(1)-xleft)))/
     . (2.*a(3)) - 0.5))
      iright = int(sngl((-a(2) + dsqrt(a(2)**2-4.*a(3)*(a(1)-right)))/
     . (2.*a(3)) + 0.5))

      if (a(4) .ne. 0.) then
         ipt = max0(ileft-3,1)
         jpt = min0(ileft+3,npts)
         ileft = ipt + int(6.*(xleft-wl(ipt))/(wl(jpt)-wl(ipt)))
         ipt = max0(iright-3,1)
         jpt = min0(iright+3,npts)
         iright = ipt + int(6.*(right-wl(ipt))/(wl(jpt)-wl(ipt)))
      endif
      ileft = max0(1,ileft)
      iright = min0(npts,iright)
      return


c*****do the estimation for linears
60    if (a(2) .eq. 0.) then
         ileft = int(xleft)
         iright = int(right)
      else
         ileft = max0(1,int(sngl((xleft-a(1))/a(2))))
         iright = min0(npts,int(sngl((right-a(1))/a(2))))
      endif
      return


c*****do the estimation for Chebyshevs, Splines, and Legendres
40    if (xleft .lt. wl(1)) then
         ileft = 1
         go to 55
      endif
      do i=2,npts
         if (xleft .le. wl(i)) then
            if (abs(xleft-wl(i)) .gt. abs(xleft-wl(i-1))) then
               ileft = i-1
            else
               ileft = i
            endif
            go to 55
         endif
      enddo        
55    if (right .gt. wl(npts)) then
         iright = npts
         return
      endif
      do i=npts-1,1,-1
         if (right .ge. wl(i)) then
            if (abs(right-wl(i)) .gt. abs(right-wl(i+1))) then
               iright = i+1
            else
               iright = i
            endif
            return
         endif
      enddo     
      return
               

      end




