
      subroutine getnum (nchars,line,xnum,xnumsngl)
c******************************************************************************
c     this routine gets a double precision number 'xnum' and its
c     equivalent single precision number 'xnumsngl' from screen input;
c     if a non-number is typed in, the user will be queried again
c******************************************************************************

      implicit real*8 (a-h,o-z)
      real*8 xnum
      real*4 xnumsngl
      integer line


      xnum = -9999.
1     call getasci (nchars,line)
      if (nchars .lt. 0) return
      call number (nchars,line,xnum)
         if (xnum .eq. -9999.) go to 1
      xnumsngl = sngl(xnum)
      return

      end

