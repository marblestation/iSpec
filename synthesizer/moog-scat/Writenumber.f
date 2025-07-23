
      subroutine writenumber (xnum)
c******************************************************************************
c     This subroutine decides on the order-of-magnitude of a number and
c     writes it to a character string in a reasonable format
c******************************************************************************

      include 'Pstuff.com'
      real*4 xnum, lognum
      integer numdec


      lognum = alog10(abs(xnum))
      if     (lognum .ge. 6.) then
         write (errmess,1002) 
      elseif (lognum .ge. 0.) then 
         numdec = 5 - nint(lognum)
         write (errmess,1001) numdec
      else   
         write (errmess,1003)
      endif
      write (array,errmess) xnum
      return


c*****format statements
1001  format ('(f8.',i1,'$)')
1002  format ('(1pe8.1$)')
1003  format ('(f8.5$)')


      end



