
      subroutine synpop
c******************************************************************************
c     Synthetic spectrum routine for stellar populations; usually meant for
c     integrated light spectra.  
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Quants.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Pstuff.com'
      include 'Multimod.com'
      real*8 tempspec(5,10000), rspec(10000)
      character*80 holdline(5,30)
      character*80 line

  
c*****read the parameter file 
      call params


c*****open the model table input file and the summary table output file;
c     read the information from the table input file
      call tablepop (2)


c*****open the standard output file
60    nf1out = 20
      lscreen = 8
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)


c*****FIRST PASS:  For each model, compute a raw synthetic spectrum;
c*****the starting do/if loops are for isotopic things only
      xhyd = 10.0**xsolar(1)
      do mmod=1,modtot
         if (numiso .gt. 0) then
            if (nisos .gt. 0) then
               do k=1,numiso
                  do l=1,nisos
                     if (isotope(k) .eq. isospecial(l)) then
                        do m=1,numisosyn
                           isoabund(k,m) = fracspecial(mmod,l)
                        enddo
                     endif
                  enddo
               enddo
            endif
         endif


c*****read in the model atmospheres and their summary output files
         line = synpre
         num = 80
         call getcount (num,line)
         if (mmod .lt. 10) then
            write (line(num+1:num+1),1013) mmod
         else
            write (line(num+1:num+2),1014) mmod
         endif
         nf2out = 21
         fmodoutput(mmod) = line
         f2out = fmodoutput(mmod)
         lscreen = 12
         array = 'INDIVIDUAL MODEL RAW SYNTHESIS OUTPUT'
         nchars = 37
         call infile ('output ',nf2out,'formatted  ',0,nchars,
     .                f2out,lscreen)
         line = modpre
         num = 80
         call getcount (num,line)
         if (mmod .lt. 10) then
            write (line(num+1:num+1),1013) mmod
         else
            write (line(num+1:num+2),1014) mmod
         endif
         nfmodel = 30
         fmodinput(mmod) = line
         fmodel = fmodinput(mmod)
         lscreen = 14
         array = 'INDIVIDUAL MODEL ATMOSPHERE'
         nchars = 27
         call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .                fmodel,lscreen)
         call inmodel
         write (nf2out,1001) moditle


c*****change the "default" abundances to be the ones read from the table file
         do i=1,nabs
            iabund = nint(elspecial(i))
            xabu(iabund) = 10.**abspecial(mmod,i)/xhyd
         enddo


c*****open the line list file and the strong line list file
         nflines = 31
         lscreen = 16
         array = 'THE LINE LIST'
         nchars = 13
         call infile ('input  ',nflines,'formatted  ',0,nchars,
     .                flines,lscreen)
         if (dostrong .gt. 0) then
            nfslines = 32
            lscreen = 18
            array = 'THE STRONG LINE LIST'
            nchars = 20
            call infile ('input  ',nfslines,'formatted  ',0,nchars,
     .                   fslines,lscreen)
         endif


c*****do the syntheses
         ncall = 1
         if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
            isorun = 1
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
               call eqlib
               call nearly (1)
               call synspec
               linprintopt = 0
            enddo
         endif
         fluxmod(mmod,1) = flux
         call finish (0)
      enddo


c*****clear the mean spectrum array and the information array
      do i=1,5
         do j=1,10000
            tempspec(i,j) = 0.
         enddo
      enddo
      do i=1,5
         do j=1,30
            write (holdline(i,j),1008)
         enddo
      enddo


c*****open the file where the mean raw spectrum will be written
      nf9out = 23
      lscreen = 10
      array = 'MEAN RAW SYNTHESIS OUTPUT'
      nchars = 25
      call infile ('output ',nf9out,'formatted  ',0,nchars,
     .             f9out,lscreen)


c*****compute normalized weights
      weighttot = 0.
      do mmod=1,modtot
         weightmod(mmod,1) = fluxmod(mmod,1)*radius(mmod)**2*      
     .                             relcount(mmod)
         weighttot = weighttot + weightmod(mmod,1)
      enddo
      do mmod=1,modtot
         weightmod(mmod,1) = weightmod(mmod,1)/weighttot
      enddo
      

c*****read back the syntheses, compute the weighted average 
      do mmod=1,modtot
         newunit = 26
         open (newunit,file=fmodoutput(mmod))
         read (newunit,1001) line
         do j=1,numatomsyn
            lincount = 0
50          read (newunit,1001) line
            lincount = lincount + 1
            if (mmod .eq. 1) write (holdline(j,lincount),1001) line
            if (line(1:5) .eq. 'MODEL') then
               read (newunit,1001) line
               lincount = lincount + 1
               if (mmod .eq. 1) write (holdline(j,lincount),1001) line
               read (line,*) wavemod1, wavemod2, wavestep
               nw1 = nint(1000.*wavemod1)
               nw2 = nint(1000.*wavemod2)
               nws = nint(1000.*wavestep)
               nnn = (nw2-nw1)/nws + 1
               read (newunit,*) (rspec(n),n=1,nnn)
               do n=1,nnn
                  tempspec(j,n) = tempspec(j,n) + weightmod(mmod,1)*
     .                            rspec(n)
               enddo
            else
               go to 50
            endif
         enddo
         close (unit=newunit)
      enddo


c     write the average spectrum back to disk.
      line(1:7) = 'MODEL: '
      line(8:80) = popitle(1:73) 
      do j=1,numatomsyn
         holdline(j,lincount-1) = line
         write (nf9out,1001) (holdline(j,l),l=1,lincount)
         write (nf9out,1003) (tempspec(j,n),n=1,nnn)
      enddo
      close (unit=nf9out)


c*****now plot the spectrum
      if (plotopt.eq.2 .and. specfileopt.gt.0) then
         nfobs = 33
         lscreen = 12
         array = 'THE OBSERVED SPECTRUM'
         nchars = 21
         if (specfileopt.eq.1 .or. specfileopt.eq.3) then
            call infile ('input  ',nfobs,'unformatted',2880,nchars,
     .                   fobs,lscreen)
         else
            call infile ('input  ',nfobs,'formatted  ',0,nchars,
     .                   fobs,lscreen)
         endif
      endif
      if (plotopt .ne. 0) then
         nf2out = nf9out
         f2out = f9out
         nf2out = 21
         lscreen = lscreen + 2
         array = 'RAW SYNTHESIS OUTPUT'
         nchars = 20
         call infile ('output ',nf2out,'formatted  ',0,nchars,
     .                f2out,lscreen)
         nf3out = 22
         lscreen = lscreen + 2
         array = 'SMOOTHED SYNTHESES OUTPUT'
         nchars = 25
         call infile ('output ',nf3out,'formatted  ',0,nchars,
     .                f3out,lscreen)
         if (f5out .ne. 'optional_output_file') then
            nf5out = 26
            lscreen = lscreen + 2
            array = 'POSTSCRIPT PLOT OUTPUT'
            nchars = 22
            call infile ('output ',nf5out,'formatted  ',0,nchars,
     .                   f5out,lscreen)
         endif
         if (plotopt .eq. 3) then
            call smooth (-1,ncall)
            choice = 'q'
         else
            call pltspec (lscreen,ncall)
         endif


c*****if needed, loop back with abundance changes
         if (choice .eq. 'n') then
            call chabund
            if (choice .eq. 'q') go to 60
         endif
      endif


c*****format statements
1001  format (a80)
1002  format ('POPULATION SYNTHESIS FOR INTEGRATED-LIGHT SPECTRA'/a80)
1003  format (10f7.4)
1004  format (a30, 2x, a10, 3f8.0)
1005  format ('#models =', i3, 5x, 'total weight =', f7.2//
     .        ('model#', i5, 5x, 'relative weight', f8.2))
1006  format (i3, 1x, f8.0, 2f8.2, 10f6.2)
1008  format (80(' '))
1009  format ('SPECIAL ELEMENTS OR NAMES OF ISOTOPES'/(10i8))
1010  format ('NO DECLARED SPECIAL ELEMENTS')
1011  format ('SPECIAL ISOTOPE NAMES'/((5f10.5)))
1012  format ('NO DECLARED SPECIAL ISOTOPES')
1013  format (i1)
1014  format (i2)
1015  format (6f13.5)


      end








