
      subroutine gridplo
c******************************************************************************
c     This routine plots a spectrum that already has been synthesized
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'


c*****examine the parameter file
1      call params


c*****open the files for: raw spectrum depths (as input), smoothed
c     spectra (output), and (if desired) IRAF-style smoothed spectra (output)
      nf2out = 21
      lscreen = 4
      array = 'RAW SYNTHESES INPUT'
      nchars = 19
      call infile ('input  ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)
      if (plotopt .ne. 0) then
         nf3out = 22
         lscreen = lscreen + 2
         array = 'SMOOTHED SYNTHESES OUTPUT'
         nchars = 25
         call infile ('output ',nf3out,'formatted  ',0,nchars,
     .                f3out,lscreen)
         nf5out = 26
         lscreen = lscreen + 2
         array = 'POSTSCRIPT PLOT OUTPUT'
         nchars = 22
         call infile ('output ',nf5out,'formatted  ',0,nchars,
     .                f5out,lscreen)
      endif
      if (plotopt .gt. 1) then
         nf6out = 27
         lscreen = lscreen + 2
         array = 'SPECTRUM COMPARISON OUTPUT'
         nchars = 27
         call infile ('output ',nf6out,'formatted  ',0,nchars,
     .                f6out,lscreen)
      endif
      if (iraf .ne. 0) then
         nf4out = 23
         lscreen = lscreen + 2
         array = 'IRAF ("rtext") OUTPUT'
         nchars = 24
         call infile ('output ',nf4out,'formatted  ',0,nchars,
     .                f4out,lscreen)
      endif


c*****now plot the spectrum
      if (plotopt.eq.2 .and. specfileopt.gt.0) then
         nfobs = 33               
         lscreen = lscreen + 2
         array = 'THE OBSERVED SPECTRUM'
         nchars = 21
         if     (specfileopt.eq.1 .or. specfileopt.eq.3) then
            call infile ('input  ',nfobs,'unformatted',2880,nchars,
     .                   fobs,lscreen)
         else
            call infile ('input  ',nfobs,'formatted  ',0,nchars,
     .                   fobs,lscreen)
         endif
      endif
      if (plotopt .ne. 0) then
         line = 10
         ncall = 1
         call pltspec (line,ncall)
      endif


c*****finish
      if (control .ne. 'gridend') then
         call finish (1)
         go to 1
      else
         call finish (0)
      endif
      return


      end   



