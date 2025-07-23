
      subroutine putasci (num,line)
c******************************************************************************
c     this routine prints out the characters contained in 'array'
c******************************************************************************

      include 'Pstuff.com'

      istat = ivmove(line-1,1)
      istat = ivcleol()
      istat = ivmove(line-1,1)
      write (errmess,1001) num
      write (*,errmess) array
      return


c*****format statements
1001  format ('(a',i2,'$)')


      end







