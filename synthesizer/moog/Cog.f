
      subroutine cog
c******************************************************************************
c     This routine produces sets of curves-of-growth                    
c******************************************************************************
 
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
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
      nf5out = 26
      lscreen = lscreen + 2
      array = 'POSTSCRIPT PLOT OUTPUT'
      nchars = 22
      call infile ('output ',nf5out,'formatted  ',0,nchars,
     .             f5out,lscreen)


c*****open and read the model atmosphere
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
102   nfmodel = 30
      lscreen = lscreen + 2
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****open and read the line list file; get ready for the line calculations
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .             flines,lscreen)
      isynth = 1
101   call inlines (1)
      call eqlib
      call nearly (1)


c*****define the range of lines (the whole list, in this case)
      mode = 1
      call linlimit
      if (lim1line .lt. 0) then
         call finish (0)
         return
      endif

 
c*****do the curves of growth, making plots if desired
      do lim1=lim1line,lim2line
         lim2 = lim1
         call curve
         call pltcog
         if (choice .eq. 'm') then
            close (unit=nfmodel)
            close (unit=nflines)
            rewind nf1out
            rewind nf2out
            rewind nf5out
            array = 'THE NEW MODEL ATMOSPHERE'
            nchars = 24
            fmodel =  'no_filename_given'
            go to 102
         endif
         if (choice .eq. 'v') go to 101
         array = 'DO ANOTHER CURVE-OF-GROWTH ([y]/n)? '
         nchars = 36
         lscreen = 16
         call getasci (nchars,lscreen)
         choice = chinfo(1:1)
         if (choice .eq. 'n') then
            call finish (0)
            return
         endif
      enddo


c****end the computations
      call finish (0)
      return


      end                                                               


