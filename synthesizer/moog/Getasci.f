
      subroutine getasci (num,line)
c******************************************************************************
c     this routine asks for information, by printing out the
c     characters contained in 'array', and returns what the
c     user types out in array 'chinfo'.  Variable 'num' on input
c     is the number of characters in 'message', and on output is
c     the number of characters in 'chinfo'
c******************************************************************************

      include 'Pstuff.com'

      write (chinfo,1005)
1005  format (80(' '))

      istat = ivmove(line-1,1)
      istat = ivcleol()
      istat = ivmove(line-1,1)
      if (num .lt. 10) then
         write (errmess,1001) num
1001     format ('(a',i1,'$)')
      else
         write (errmess,1002) num
1002     format ('(a',i2,'$)')
      endif
      write (*,errmess) array
      num = 80 - num
      if (num .lt. 10) then
         write (errmess,1003) num
1003     format ('(a',i1,')')
      else
         write (errmess,1004) num
1004     format ('(a',i2,')')
      endif
      read (*,errmess) chinfo
      call getcount (num,chinfo)
       return

      end



