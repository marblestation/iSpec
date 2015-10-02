
      subroutine gridsyn
c******************************************************************************
c     This program can synthesize multiple sections of spectra for multiple
c     input model atmospheres
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'


c*****examine the parameter file
1     call params
      linprintopt = linprintalt


c*****open the files for: standard output, raw spectrum depths, smoothed 
c     spectra, and (if desired) IRAF-style smoothed spectra
      nf1out = 20     
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21               
      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS OUTPUT'
      nchars = 20
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)
      if (plotopt .gt. 0) then
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


c*****open and read the model atmosphere file
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****open the line list file and the strong line list file
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
      if (dostrong .gt. 0) then
         nfslines = 32
         lscreen = lscreen + 2
         array = 'THE STRONG LINE LIST'
         nchars = 20
         call infile ('input  ',nfslines,'formatted  ',0,nchars,
     .                 fslines,lscreen)
      endif
      

c*****do the syntheses
      if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
         isorun = 1
         isynth = 1
         nlines = 0
         mode = 3
         call inlines (1)
         call eqlib
         call nearly (1)
         call synspec
      else
         do n=1,numatomsyn
            isynth = n
            isorun = n
            start = oldstart
            sstop = oldstop
            mode = 3
            call inlines (1)
            molopt = 2
            call eqlib
            call nearly (1)
            call synspec
            linprintopt = 0
         enddo
      endif
         

c*****now plot the spectrum
      if (plotopt.eq.2 .and. specfileopt.gt.0) then
         nfobs = 33               
         lscreen = lscreen + 2
         array = 'THE OBSERVED SPECTRUM'
         nchars = 21
         if (plotopt.eq.1 .or. specfileopt.eq.3) then
            call infile ('input  ',nfobs,'unformatted',2880,nchars,
     .                   fobs,lscreen)
         else
            call infile ('input  ',nfobs,'formatted  ',0,nchars,
     .                   fobs,lscreen)
         endif
      endif
      if (plotopt .ne. 0) then
         ncall = 1
         call pltspec (lscreen,ncall)
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






