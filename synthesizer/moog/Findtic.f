
      subroutine findtic (end1,end2,bigtic,smltic)
c******************************************************************************
c     this routine makes tic marks at nice round intervals.
c******************************************************************************

      implicit real*8(a-h,o-z)
      real*4 end1, end2, bigtic, smltic, span, spanlog


      span = end2 - end1
      if (span .lt. 0.) span = - span
      spanlog = alog10(span)
      size = ifix(spanlog)
      if (spanlog .lt. 0.) size = size - 1
      chop = spanlog - size
      if (chop .lt. 0.31) then
         bigtic = 10.**(size)/2
         smltic = bigtic/5
      elseif (chop .lt. 0.71) then
         bigtic = 10.**(size)
         smltic = bigtic/5
      else
         bigtic = 10.**(size)*2
         smltic = bigtic/4
      endif
      return


      end






