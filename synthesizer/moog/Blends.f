
      subroutine blends
c******************************************************************************
c     This program does abundance derivations from blended spectral features;
c     only one element will have its abundance derived per run.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Pstuff.com'
      include 'Dummy.com'
 

c*****examine the parameter file
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


c*****open and read the model atmosphere
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****initialize some variables
      isynth = 1
      isorun = 1
      iatom = nint(cogatom)
      pec(iatom) = 1
      numpecatom = 1
      pecabund(iatom,1) = 0.


c*****open and read the line list file; get ready for the line calculations
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
100   call inlines (1)
      call eqlib
      call nearly (1)


c*****start the large loop that will go through each blended feature
      ewsynthopt = -1
      mode = 4
30    do lll=1,1000
         if (lim2line .eq. nlines) exit


c*****define the set of lines responsible for a blended feature
         call linlimit
         lim1 = lim1line
         lim2 = lim2line
         write (99,1007) iatom, wave1(lim1), wave1(lim2)


c*****make sure that the element whose abundance is to be fit has
c     a representative line of the blend
         ifind = 0
         do j=lim1,lim2
            if (atom1(j) .lt. 100.) then
               if (iatom .eq. int(atom1(j))) then
                  ifind = 1
                  exit
               endif
            else
               call sunder (atom1(j),ia,ib)
               if (iatom.eq.ia .or. iatom.eq.ib) then
                  ifind = 1
                  exit
               endif
            endif
         enddo
         if (ifind .eq. 0) then
            do j=lim1,lim2
               abundout(j) = 999.99
            enddo
            write (nf1out,1002)
            lim1line = lim2line + 1
            if (lim1line .le. nlines+nstrong) cycle
         endif
 

c*****do the syntheses, forcing each abundance to predict the
c     feature equivalent width; a limit of 30 iterations is imposed
c     arbitrarily
         ncurve = 1
         start = wave1(lim1) - delwave
         sstop = wave1(lim2) + delwave
         delta = wave1(lim2) - wave1(lim1) + delwave
         oldstart = start
         oldstop = sstop
         oldstep = step
         olddelta = delta
         rwlgobs = dlog10(width(lim1)/wave1(lim1))
         gf1(ncurve) = 1.
         do k=1,30
            call synspec
            call total
            error = (w(ncurve)-width(lim1))/width(lim1)
            ratio = width(lim1)/w(ncurve)
            ncurve = ncurve + 1


c*****here we go for another iteration
            if (dabs(error) .ge. 0.0075) then
               rwlcomp = dlog10(w(ncurve-1)/wave1(lim1))
               if     (rwlcomp.lt.-5.2 .and. rwlgobs.lt.-5.2) then 
                  ratio = ratio
               elseif (rwlcomp.ge.-5.2 .and. rwlgobs.ge.-5.2) then 
                  ratio = ratio**2.0
               else   
                  ratio = ratio**1.5
               endif
               gf1(ncurve) = gf1(ncurve-1)*ratio
               do j=lim1,lim2
                  if (atom1(j) .gt. 100.) then
                     call sunder (atom1(j),ia,ib)
                     if (ia.eq.iatom .or. ib.eq.iatom) then
                        do i=1,ntau                                
                           kapnu0(j,i) = kapnu0(j,i)*ratio            
                        enddo
                     endif
                  elseif (int(atom1(j)) .eq. iatom) then
                     do i=1,ntau                                
                        kapnu0(j,i) = kapnu0(j,i)*ratio            
                     enddo
                  endif
               enddo
               if (k .eq. 20) then
                  write (*,1008)
                  stop
               endif
            else
               exit
            endif
         enddo


c*****here we do the final calculation when the predicted and observed are close
         gf1(ncurve) = gf1(ncurve-1)*ratio
         do j=lim1,lim2
            if (atom1(j) .gt. 100.) then
               call sunder (atom1(j),ia,ib)
               if (ia.eq.iatom .or. ib.eq.iatom) then
                  do i=1,ntau                                
                     kapnu0(j,i) = kapnu0(j,i)*ratio            
                  enddo
               endif
            elseif (int(atom1(j)) .eq. iatom) then
               do i=1,ntau                                
                  kapnu0(j,i) = kapnu0(j,i)*ratio            
               enddo
            endif
         enddo
         call synspec
         call total
         widout(lim1) = w(ncurve)
         diff = dlog10(gf1(ncurve))
         abundout(lim1) = dlog10(xabund(iatom)) + 12.0 + diff
         if (ncurve .ne. 1) then
            write (nf1out,1001) ncurve
         endif


c*****here is where some auxiliary things like mean depth are computed
         if (lim1.eq.lim2 .and. linprintopt.ge.3) then
            wave = wave1(lim1)
            call taukap
            call cdcalc(2)
            first = 0.4343*cd(1)
            d(1) = rinteg(xref,cd,dummy1,ntau,first)
            do i=1,ntau
               dummy1(i) = xref(i)*cd(i)
            enddo
            first = 0.
            cdmean = rinteg(xref,dummy1,dummy2,ntau,first)/
     .      rinteg(xref,cd,dummy2,ntau,first)
            do i=1,ntau
               if (cdmean .lt. cd(i)) exit
            enddo
            write (nf1out,1005) lim1, cdmean, i, xref(i)
            do i=1,ntau
               if (taunu(i)+taulam(i) .ge. 1.) exit
            enddo
            write (nf1out,1006) lim1, i, dlog10(tauref(i)),
     .                          dlog10(taulam(i)), dlog10(taunu(i))
         endif


c*****assign the abundance to the strongest line of the blend that
c     contains cogatom; put 999.99's as the abundances of all but 
c     this line; go back for another blended feature
         if (lim2 .gt. lim1) then
            abunblend = abundout(lim1)
            widblend = widout(lim1)
            strongest = 0.
            linstrongest = 0
            do j=lim1,lim2
               abundout(j) = 999.99
            enddo
            do j=lim1,lim2
               if (atom1(j) .lt. 100.) then
                  if (dint(atom1(j)) .eq. cogatom) then
                     if (kapnu0(j,jtau5) .gt. strongest) then
                        strongest = kapnu0(j,jtau5)
                        linstrongest = j
                     endif
                  endif
               else
                  call sunder (atom1(j),ia,ib)
                  if (dble(ia).eq.cogatom .or. dble(ib).eq.cogatom) then
                     if (kapnu0(j,jtau5) .gt. strongest) then
                        strongest = kapnu0(j,jtau5)
                        linstrongest = j
                     endif
                  endif
               endif
            enddo
            abundout(linstrongest) = abunblend
            widout(linstrongest) = widblend
         endif
         lim1line = lim2line + 1
         if (lim1line .le. nlines+nstrong) cycle
      enddo


c*****do abundance statistics; print out a summary of the abundances 
      lim1obs = 1
      lim2obs = nlines+nstrong
      call stats
      rewind nf2out
      call lineinfo (3)


c*****here a plot may be made on the terminal (and paper) if there
c     are enough lines; then the user will be prompted on some
c     options concerning what is seen on the plot
      if (plotopt .ne. 0) then
         call pltabun
         if (choice.eq.'v' .or. choice.eq.'m') go to 100
      endif


c*****now the option will be to redo the molecular equilibrium and redo
c     the last species, at the user's option.
      if (neq .ne. 0) then
         do n=1,neq
            if (iatom .eq. iorder(n)) then
               write (array,1004)
               ikount = min0(nlines+11,maxline)
               nchars = 70
35             call getasci (nchars,ikount)
               choice = chinfo(1:1)
               if (choice.eq.'y' .or. nchars.le.0) then
                  xabund(iatom) = 10.**(average-12.)
                  call eqlib
                  call nearly (1)
                  lim1line = 0
                  lim2line = 0
                  rewind nf2out
                  go to 30
               elseif (choice .eq. 'n') then
                  call finish (0)
               else
                  go to 35
                endif
            endif
         enddo
      endif
 

c*****exit the program
      call finish (0)


c*****format statements
1001  format (' This fit required ',i2,' iterations'/)
1002  format ('WARNING: NO ELEMENT LINE TO VARY IN BLEND')
1004  format ('SPECIES USES MOLECULAR EQUILIBRIUM!',
     .        '  REDO WITH NEW ABUNDANCE ([y]/n)? ')
1005  format (/'LINE ', i5, ':',
     .        ' weighted mean line contribution function C_d =',
     .        f6.2/ '  which occurs near level ', i3,
     .        ' with log tauref = ', f6.2)
1006  format (/'LINE ', i5, ':',
     .        '  tau(total) is greater than 1 at level',i3/
     .        '  logs of tauref, taulam, taunu =', 3f6.2)
1007  format (/'VARYING THIS ELEMENT: ', i8/
     .        'USING THE LINE GROUP IN THE RANGE: ', 2f10.3)
1008  format (' MAX OF 30 ITERATIONS REACHED; I QUIT!')


      end



