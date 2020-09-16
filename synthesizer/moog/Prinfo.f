
      subroutine prinfo (line)
c******************************************************************************
c     this routine prints out information on the text terminal
c******************************************************************************
 
      include 'Pstuff.com'
      include 'Atmos.com'
      character errm1*80,array1*80


c*****do not print out information in real time if the code is in 
c     batch mode
      if (silent .eq. 'y') return

       
      if (line .eq. 1) then
         istat = ivcleof(4,1)
      endif

      if (line .eq. maxline-5) then
         errm1 = errmess
         array1 = array
10       array = 'WANT TO SEE MORE ([y]/n)? '
         nchars = 26
         call getasci (nchars,4+line)
         if (chinfo(1:1).eq.'y' .or. nchars.le.0) then
            istat = ivcleof(4,1)
            line = 1
            array = array1
            errmess = errm1
         elseif (chinfo(1:1) .eq. 'n') then
            errmess = 'stopinfo!'
            return
         else
            go to 10
         endif
      endif


      istat = ivwrite(4+line,1,array,79)
      errmess = 'everything OK'
      return


      end












