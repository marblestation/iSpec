
      subroutine abfind
c******************************************************************************
c     This routine derives abundances for individual atomic lines       
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'

  
c*****read the parameter file
      call params


c*****open the files for standard output and summary abundances
      nf1out = 20
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21
      lscreen = lscreen + 2
      array = 'SUMMARY ABUNDANCE OUTPUT'
      nchars = 24
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
      array = 'THE LINE LIST'
      nchars = 13
      nflines = 31
      lscreen = lscreen + 2
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)


c*****get ready for the line calculations: generate a curve-of-growth lookup
c     table, read the linelist, etc.
100   call fakeline
      call inlines (1)
      call eqlib
      call nearly (1)

 
c*****set some parameters
      ewsynthopt = -1
      mode = 2
      cogatom = 0.
      lim1line = 0


c*****find the range of lines of a species
5     call linlimit
      if (lim1line .lt. 0) then
         call finish (0)
         return
      endif
      lim1obs = lim1line
      lim2obs = lim2line


c*****find out whether molecular equilibrium is involved in the species
      call molquery


c*****force each abundance of a species member to predict the 
c     line equivalent width; here is the code for ordinary species
      if (molflag .eq. 0) then
         abundin =  dlog10(xabund(iabatom)) + 12.0
         do lim1=lim1line,lim2line
            call lineabund (abundin)
         enddo
         call stats
         call lineinfo (3)
      else


c*****and here is the code for species involved in molecular equilibrium;
c     this procedure is iterated until input and output abundances are in 
c     agreement
         iternumber = 1
10       abundin =  dlog10(xabund(iabatom)) + 12.0
         do lim1=lim1line,lim2line
            call lineabund (abundin)
         enddo
         call stats
         call lineinfo (3)
         if (t(jtau5).lt.3800                 .or. 
     .       atom1(lim1line).gt.100.0         .or.
     .       int(atom1(lim1line)+0.0001).eq.6 .or.
     .       int(atom1(lim1line)+0.0001).eq.8) then
            if (iternumber .lt. 6) then
               if (dabs(average-abundin) .gt. 0.02) then
                  xabund(iabatom) = 10.**(average-12.)
                  iternumber = iternumber + 1
                  call eqlib
                  call nearly (2)
                  go to 10
               else
                  write (array,1001) iternumber
                  ikount=kount+14
                  nchars = 53
                  call putasci (nchars,ikount)
               endif
            else
               write (array,1003) molecule, dlog10(abundin),
     .                            dlog10(average)
               lscreen = lscreen + 2
               call prinfo (lscreen)
               stop
            endif
         endif
      endif


c*****here a plot may be made on the terminal (and paper) if there 
c     are enough lines; then the user will be prompted on some
c     options concerning what is seen on the plot
      if (plotopt .ne. 0) then
         call pltabun
         if     (choice.eq.'v') then
            rewind nf1out
            rewind nf2out
            write (nf2out,1002) linitle,moditle
            choice = ' '
            go to 100
         elseif (choice .eq. 'm') then
            close (unit=nfmodel)
            close (unit=nflines)
            rewind nf1out
            rewind nf2out
            rewind nf5out
            array = 'THE NEW MODEL ATMOSPHERE'
            nchars = 24
            fmodel =  'no_filename_given'
            lim1line = 0
            lim2line = 0
            lim1obs = 0
            lim2obs = 0
            go to 102
         endif
      endif


 
             
c*****quit, or go on to another species?
      if (silent .eq. 'y') then
         choice = 'y'
         nchars = 0
      else
         array = 'DO ANOTHER SPECIES ([y]/n)? '
         nchars = 28
         call getasci (nchars,maxline)
         choice = chinfo(1:1)
      endif
      if (choice.eq.'y' .or. nchars.le.0) then
         if (mode .eq. 2) then
            go to 5
         else
            call finish (0)
            return
         endif
      else
         call finish (0)
         return
      endif


c*****format statements
1001  format ('THIS REQUIRED', i2,' ITERATIONS WITH MOLECULAR ',
     .        'EQUILIBRIUM')
1002  format (a80)
1003  format ('FOR SPECIES ', f10.1,' NO CONVERGENCE: '/
     .        'LAST ITERATIONS YIELD', 2f10.3, '  I QUIT!')


      end








