      subroutine ewfind_scat
c******************************************************************************
c     This routine predicts equivalent widths for individual atomic lines
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      include 'Scat.com'
      real*8 taunu0(100)
      character*4 ion(3)
      data ion/' I  ', ' II ', ' III'/


c*****read the parameter file
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
      array = 'SUMMARY PREDICTED EW OUTPUT'
      nchars = 27
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
      isynth = 1
      call inlines (1)
      call eqlib
      call nearly (1)

 
c*****set some parameters and write header stuff to output
      ewsynthopt = +1
      mode = 1
      call linlimit
      lim1obs = lim1line
      lim2obs = lim2line
      istat = ivcleof (4,1)
      write (nf2out,1002) linitle,moditle


c*****run single line computations once to predict the EW for each line
      do lim1=lim1line,lim2line
         call molquery
         write (array,1001)
         lscreen = lscreen + 2
c        call prinfo (lscreen)
         write (nf2out,1001)
         ncurve = lim1
         lim2 = lim1
         gf1(ncurve) = gf(lim1)
         call oneline_scat (1)
         widout(lim1) = w(ncurve)
         iatom = atom1(lim1)


         if (iatom .lt. 100) then
            xab = dlog10(xabund(iatom)) + 12.
         else
            xab = dlog10(xabund(iabtom)) + 12.
         endif
         ich = idint(charge(lim1) + 0.1)
         if (iatom .lt. 100) then
            write (array,1003) wave1(lim1), e(lim1,1), 
     .                         dlog10(gf(lim1)), names(iatom), 
     .                         ion(ich), xab, 1000.*widout(lim1)
            lscreen = lscreen + 2
c           call prinfo (lscreen)
            write (nf2out,1003) wave1(lim1), e(lim1,1), 
     .                          dlog10(gf(lim1)), names(iatom), 
     .                          ion(ich) ,xab, 1000.*widout(lim1)
         else
            write (array,1004) wave1(lim1), e(lim1,1), 
     .                         dlog10(gf(lim1)), names(iaa), 
     .                         names(ibb), xab, 1000.*widout(lim1)
            lscreen = lscreen + 2
c           call prinfo (lscreen)
            write (nf2out,1004) wave1(lim1), e(lim1,1), 
     .                          dlog10(gf(lim1)), names(iaa), 
     .                          names(ibb), xab, 1000.*widout(lim1)
         endif


c*****(re)compute the line optical depth at line center and the C_d curve

         do i=1,ntau
            kapnu(i) = kapnu0(lim1,i)
            dummy1(i) = tauref(i)*kapnu(i)/(0.4343*kapref(i))
         enddo
         first = tauref(1)*kapnu(1)/kapref(1)
         dummy2(1) = rinteg(xref,dummy1,taunu0,ntau,0.)
         taunu0(1) = first
         do i=2,ntau
            taunu0(i) = taunu0(i-1) + taunu0(i)
         enddo
         do i=1,ntau
            taunu(i) = taunu0(i)
         enddo
c         call cdcalc (2)
         call cdcalc_SCAT (2)
         if (linprintopt .ge. 2) then
            write (nf2out,1010) 
            write (nf2out,1011) (i, rhox(i), xref(i), int(t(i)), 
     .                          pgas(i), rho(i), kaplam(i), 
     .                          taulam(i), taunu0(i), adepth, i=1,ntau)
         endif


c*****compute layer where continuum optical depth > 1
         do i=1,ntau
            if (taulam(i) .ge. 1.) then
               write (nf2out,1013) tauref(i), i
               go to 10
            endif
         enddo


c     compute layer where line center optical depth > 1
10          if (taunu0(ntau) .lt. 1.) then
               write (nf2out,1016)
               go to 20
            endif
            do i=1,ntau
            if (taunu0(i) .ge. 1.) then
               write (nf2out,1014) tauref(i), i
               go to 20
            endif
         enddo


c     compute layer where line center plus continuum optical depth > 1
20       do i=1,ntau
            if (taunu0(i)+taulam(i) .ge. 1.) then
               write (nf2out,1015) tauref(i), i
               go to 30
            endif
         enddo


*****compute mean line-center formation level (weight: contribution function)
c    SCATTERING TO DO: create a two-dimensional array that stores adepth as a function of line
c    frequency; then consider the code below
30       do i=1,ntau
c            dummy1(i) = xref(i)*dabs(cd(i))
            dummy1(i) = xref(i)*dabs(adepth)
         enddo
         xrefcdinteg = rinteg(xref,dummy1,dummy2,ntau,0.)
         do i=1,ntau
c            dummy1(i) = dabs(cd(i))
            dummy1(i) = dabs(adepth)                                     | RE-CONSIDER
         enddo
         cdinteg = rinteg(xref,dummy1,dummy2,ntau,0.)
         xrefmean = xrefcdinteg/cdinteg
         do i=1,ntau
            if (xrefmean .le. xref(i)) then
               write (nf2out,1017) 10**(xrefmean), i
               go to 40
            endif
         enddo
40       continue
      enddo


c*****end the abundance computations
      call finish (0)
      return


c*****format statements
1001  format (/'wavelength        EP     logGF     ident',
     .        '     Abund    EWcalc')
1002  format (a80)
1003  format (f10.2,f10.2,f10.3,'     ',a2,a3,f10.2,f10.1/)
1004  format (f10.2, f10.2, f10.3, 6x, a2, a2, f10.2, f10.1/)
1010     format (' i', 5x, 'rhox', 2x, 'xref', 5x, 'T', 5x, 'Pgas', 
     .           6x, 'rho', 3x, 'kaplam', 3x, 'taulam',
     .           3x, 'taunu0', 8x, 'Cd')
1011  format (i2, 1pd9.2, 0pf6.2, i6, 1p5d9.2, d10.2)
1013           format (5x, 'tau(ref) =', 1pe10.2,
     .                 ' (level=',i2, ') for tau(cont) ~ 1')
1014           format (5x, 'tau(ref) =', 1pe10.2,
     .                 ' (level=',i2, ') for tau(line) ~ 1')
1015           format (5x, 'tau(ref) =', 1pe10.2,
     .                 ' (level=',i2, ') for tau(cont+line) ~ 1')
1016           format (7x,'  NOTE: line center tau(line) < 1',
     .                 '  at deepest atmosphere layer')
1017           format (5x, 'C_d weighted mean formation tau(ref) =',
     .                 1pe10.2, ' (level=',i2, ')')

      end

