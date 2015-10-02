
c***********************************************************************
      subroutine blankstring (string)
c***********************************************************************
c     this routine simply makes an 80-character string all blanks

      implicit real*8 (a-h,o-z)
      character*80 string


      do i=1,80
         string(i:i) = ' '
      enddo
      return


      end

