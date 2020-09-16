
      subroutine abpop
c******************************************************************************
c     Equivalent width matching routine for stellar populations; usually 
c     meant for integrated light spectra.  
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'
      include 'Multimod.com'
      character*80 line

  
c*****read the parameter file 
      call params


c*****open the model table input file and the summary table output file;
c     read the information from the table input file
      call tablepop (1)


c*****open the standard output file
      nf1out = 20
      lscreen = 8
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)

c*****FIRST PASS:  For each model, for each line, create a curve-of-growth
c     finely divided in steps of 0.01 in "Ngf"
      do mmod=1,modtot


c*****read in the model atmospheres and their summary output files
         line = synpre
         num = 80
         call getcount (num,line)
         if (mmod .lt. 10) then
            write(line(num+1:num+1),1013) mmod
         else
            write(line(num+1:num+2),1014) mmod
         endif
         nf2out = 21
         fmodoutput(mmod) = line
         f2out = fmodoutput(mmod)
         lscreen = 10
         array = 'SUMMARY C-O-G OUTPUT'
         nchars = 20
         call infile ('output ',nf2out,'formatted  ',0,nchars,
     .                f2out,lscreen)

         line = modpre
         num = 80
         call getcount (num,line)
         if (mmod .lt. 10) then
            write(line(num+1:num+1),1013) mmod
         else
            write(line(num+1:num+2),1014) mmod
         endif
         nfmodel = 30
         fmodinput(mmod) = line
         fmodel = fmodinput(mmod)
         array = 'THE MODEL ATMOSPHERE'
         nchars = 20
         lscreen = 12
         call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .                fmodel,lscreen)
         call inmodel
         write (nf2out,1001) moditle


c*****open and read the line list; this will internally be broken down
c     into 1-line lists for computational efficiency
         array = 'THE LINE LIST'
         nchars = 13
         nflines = 31
         lscreen = 14
         call infile ('input  ',nflines,'formatted  ',0,nchars,
     .                flines,lscreen)
         call inlines (1)
         if (nlines .gt. 1000) then
            write (*,1005)
            stop
         endif
         call eqlib
         call nearly (1)


c*****do the curves of growth; store the results
         rwlow = -6.5
         rwhigh = -3.5
         rwstep = 0.15
         do lim1=1,nlines
            lim2 = lim1
            waveold = 0.
            call curve
            if (ncurve .gt. 50) then
               nmodcurve(mmod,lim1) = 50
            else
               nmodcurve(mmod,lim1) = ncurve
            endif
            fluxmod(mmod,lim1) = flux
            lastcurve = nmodcurve(mmod,lim1)
            iatom = nint(atom1(lim1))
            xxab = dlog10(xabund(iatom))
            weightmod(mmod,lim1) = fluxmod(mmod,lim1)*radius(mmod)**2*
     .                             relcount(mmod)
            do icurve=1,lastcurve
               gfmodtab(mmod,lim1,icurve) = gf1(icurve) + xxab
               rwmodtab(mmod,lim1,icurve) = w(icurve)
            enddo
         enddo
      enddo


c*****set some parameters
         ewsynthopt = -1
         mode = 2
         cogatom = 0.
         lim1line = 0


c*****now, inside the 
c*****define the range of lines for a species
 5       call linlimit
         if (lim1line .lt. 0) then
            call finish (0)
            return
         endif
         lim1obs = lim1line
         lim2obs = lim2line


c*****find out whether molecular equilibrium is involved in the species
         call molquery


c*****for just the first model atmosphere, interpolate in the curve-of-growth
c     at very small log(Ngf) steps for just the first line in the list,
c     in order to make a lookup table for later abundance iterations
      ncurve = nmodcurve(1,1)
      do i=2,ncurve-2
         k = 30*(i-2)
         l = -1
         do m=k+1,k+30
           l = l + 1
           pp       = 0.005/0.15*(l-1)
           gftab(m) = gfmodtab(1,1,2) + 0.005*(m-1)
           rwtab(m) = rwmodtab(1,1,i-1)*(-pp)*(pp-1.)*(pp-2.)/6. +
     .                rwmodtab(1,1,i)*(pp*pp-1.)*(pp-2.)/2. +
     .                rwmodtab(1,1,i+1)*(-pp)*(pp+1.)*(pp-2.)/2. +
     .                rwmodtab(1,1,i+2)*pp*(pp*pp-1.)/6.
         enddo
      enddo
      ntabtot = m - 1


c*****now consider each line at a time; use these curve-of-growth to 
c     compute a weighted <EW_calc> 
      do lim1=lim1line,lim2line
         write (nf7out,1003) atom1(lim1), wave1(lim1), e(lim1,1),
     .                       dlog10(gf(lim1))
         rwlgerror = 50.
         deltangf = 0.
         do k=1,30
            call ewweighted
            rwlgcal = dlog10(ewmod(lim1)/wave1(lim1))
            do i=2,ntabtot
               if (rwtab(i) .gt. rwlgcal) then
                  gflgcal = gftab(i-1) + (gftab(i)-gftab(i-1))*
     .                      (rwlgcal-rwtab(i-1))/(rwtab(i)-rwtab(i-1))
                  exit
               endif
            enddo
            rwlgobs = dlog10(width(lim1)/wave1(lim1))
            do i=2,ntabtot
               if (rwtab(i) .gt. rwlgobs) then
                  gflgobs = gftab(i-1) + (gftab(i)-gftab(i-1))*
     .                      (rwlgobs-rwtab(i-1))/(rwtab(i)-rwtab(i-1))
                  exit
               endif
            enddo
            rwlgerror = rwlgobs - rwlgcal
            diffngf = gflgobs - gflgcal
            deltangf = deltangf + diffngf
            if (dabs(rwlgerror) .lt. 0.01) exit
            if (k .eq. 30) then
               write (*,1004)
               stop
            endif
         enddo


c*****here we go for a final iteration on a line, or finish
         call ewweighted
      enddo

      
c*****do the species statistics and output the results
      call stats
      call lineinfo (3)


c*****here a plot may be made on the terminal (and paper) if there 
c     are enough lines; then the user will be prompted on some
c     options concerning what is seen on the plot

      if (plotopt .ne. 0) then
         call blankstring (moditle)
         moditle(1:70) = popitle(1:70)
         moditle(57:80) = 'EW-POPULATION '
         call pltabun
      endif


*****quit, or go on to another species?
      array = 'DO ANOTHER SPECIES ([y]/n)? '
      if (silent .eq. 'n') then 
         nchars = 28
         call getasci (nchars,maxline)
         choice = chinfo(1:1)
      else
         choice = 'y'
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
1001  format (a80)
1003  format (/'species=',f5.1, 3x, 'lambda=', f8.3, 3x, 'EP=', f6.3,
     .        3x, 'log(gf)=', f7.3)
1004  format ('NO ABUNDANCE CONVERGENCE IN 30 TRIES; I QUIT!')
1005  format ('MORE THAN 1000 LINES FED TO ABPOP; I QUIT')
1013  format (i1)
1014  format (i2)
      end
