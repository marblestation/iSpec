
      subroutine defcolor (icolor)
c******************************************************************************
c     This routine decides on whether to call for a black & white hardcopy
c     plot or a colorful screen plot
c******************************************************************************

      include 'Pstuff.com'
      character*7 colors(8)


c*****assign colors to character arrays
      colors(1) = 'white  '
      colors(2) = 'red    '
      colors(3) = 'cyan   '
      colors(4) = 'yellow '
      colors(5) = 'green  '
      colors(6) = 'magenta'
      colors(7) = 'blue   '
      colors(8) = 'black  '


      if (choice.eq.'h' .or. choice.eq.'f' .or.
     .    choice.eq.'g') then
         call sm_ctype (colors(8))
      else
         call sm_ctype (colors(icolor))
      endif

      
      return
      end




