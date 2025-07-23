
      subroutine makeplot (lscreen)
c******************************************************************************
c     This subroutine does the plot-package-specific commands to begin
c     a plot, then calls the specific plot drawing routine, then ends the
c     plot.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'
      integer sm_device, lscreen, nchars


c  open the plot device: screen terminal
      if (plotroutine(1:4) .eq. 'term') then
         if (sm_device(smterm) .lt. 0) then
            write (array,1001) smterm
            istat = ivwrite(lscreen+1,1,array,79)
            write (nf1out,1007) array(1:79)
            stop
         endif
      endif


c  open the plot device: hardcopy sent to printer
      if (plotroutine(1:4) .eq. 'hard') then
         if     (plotroutine(6:9) .eq. 'land') then
            if     (sm_device('postland') .lt. 0) then
               write (array,1002)
               istat = ivwrite(lscreen+1,1,array,34)
               write (nf1out,1007) array(1:34)
               stop
            endif
         elseif (plotroutine(6:9) .eq. 'port') then
            if (sm_device('postport') .lt. 0) then
               write (array,1009)
               istat = ivwrite(lscreen+1,1,array,34)
               write (nf1out,1007) array(1:34)
               write (nf1out,1009)
               stop
            endif
         endif
      endif


c  open the plot device: postscript file
      if (plotroutine(1:4) .eq. 'file') then
         if (f5out .eq. 'optional_output_file') then
            array = 'Give the file name for the POSTSRIPT plot image: '
            nchars = 49
            call getasci (nchars,maxline)
            f5out = chinfo(1:nchars)
         else
            nchars = 80
            call getcount (nchars,f5out)
         endif
         if     (plotroutine(6:9) .eq. 'land') then
            if (nchars .lt. 10) then
               write (errmess,1003) nchars
            else
               write (errmess,1004) nchars
            endif
         elseif (plotroutine(6:9) .eq. 'port') then
            if (nchars .lt. 10) then
               write (errmess,1005) nchars
            else
               write (errmess,1006) nchars
            endif
         endif
         write (array,errmess) f5out(1:nchars)
         if (sm_device(array(1:nchars+13)) .lt. 0) then
            write (nf1out,1007) array(1:nchars+9)
            istat = ivwrite(lscreen+1,1,array,nchars+9)
            stop
         endif
      endif


c  issue standard beginning commands
      call sm_graphics
      call sm_erase


c  call the routine that makes the desired plot
      if     (plotroutine(11:14) .eq. 'cog ') then
         call cogplot
      elseif (plotroutine(11:14) .eq. 'abun') then
         call abunplot
      elseif (plotroutine(11:14) .eq. 'spec') then
         call specplot
      elseif (plotroutine(11:14) .eq. 'bin ') then
         call binplot
      elseif (plotroutine(11:14) .eq. 'flux') then
         call fluxplot
      endif


c  issue standard ending commands; exit normally
      if (plotroutine(1:4) .eq. 'file') then
         f5out = 'optional_output_file'
      endif
      call sm_gflush
      if (plotroutine(1:4).eq.'hard' .or. 
     .    plotroutine(1:4).eq.'file') call sm_hardcopy
      call sm_alpha
      return


c*****format statements
1001  format ('DEVICE OPENING ERROR FOR:',a54)
1002  format ('DEVICE OPENING ERROR FOR: postland')
1009  format ('DEVICE OPENING ERROR FOR: postport')
1007  format (a80)
1003  format ('(13hpostlandfile ,a',i1,'$)')
1004  format ('(13hpostlandfile ,a',i2,'$)')
1005  format ('(13hpostportfile ,a',i1,'$)')
1006  format ('(13hpostportfile ,a',i2,'$)')


      end


