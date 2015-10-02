
      subroutine cogsyn
c******************************************************************************
c     This program produces a curve-of-growth of a blended feature,
c     altering only a specified species' abundance.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Pstuff.com'


c*****examine the parameter file
      call params


c*****open the files for standard output and summary curves-of-growth
      nf1out = 20
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21
      lscreen = lscreen + 2
      array = 'SUMMARY C-O-G OUTPUT'
      nchars = 20
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)


c*****open and read the model atmosphere
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****open and read the line list file; get ready for the line calculations
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
      call inlines (1)
      call eqlib
      call nearly (1)


c*****do the syntheses
      isynth = 1
      isorun = 1
      ncurve = 0
      iatom =nint(cogatom)
      pec(iatom) = 1
      numpecatom = 1
      pecabund(iatom,1) = 0.
      mode = 3
10    ncurve = ncurve + 1
      nlines = 0
      call synspec
      call total
      if (ncurve .eq. 1) then
         wstart = 10.**(rwlow)*wave1(lim1)
         wstop = 10.**(rwhigh)*wave1(lim1)
      endif
      molopt = 1
      linprintopt = 0
      if (w(ncurve) .gt. wstart) then
         pecabund(iatom,1) = pecabund(iatom,1) - rwstep
         go to 10
      endif
      pecabund(iatom,1) = rwstep
20    ncurve = ncurve + 1
      nlines = 0
      call synspec
      call total
      molopt = 1
      linprintopt = 0
      if (w(ncurve) .lt. wstop) then
         pecabund(iatom,1) = pecabund(iatom,1) + rwstep
         go to 20
      endif
      call pltcog

         
c*****finish
      call finish (0)
      end






