
      subroutine ewweighted
c******************************************************************************
c     This routine computes a weighted mean EW for one line from a set 
c     of models
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Multimod.com'


c*****start the computations for a line
      if (dabs(rwlgerror) .lt. 0.01) write (nf7out,1001)
      ewweighttot = 0.
      weighttot = 0
      modcount = 0
      iatom = int(atom1(lim1)+0.0001)
      xngf = dlog10(xabund(iatom)*gf(lim1)) + deltangf


c*****for each model, interpolate in the curve-of-growth to get an EW_calc
c     for the assumed abundance
      do mmod=1,modtot
         ncurvetot = nmodcurve(mmod,lim1)
         do icurve=3,ncurvetot-2
            if (gfmodtab(mmod,lim1,icurve) .gt. xngf) then
               ic = icurve - 1
               pp = (xngf-gfmodtab(mmod,lim1,ic))/0.15
               rw = rwmodtab(mmod,lim1,ic-1)*(-pp)*(pp-1.)*(pp-2.)/6. +
     .              rwmodtab(mmod,lim1,ic)*(pp*pp-1.)*(pp-2.)/2. +
     .              rwmodtab(mmod,lim1,ic+1)*(-pp)*(pp+1.)*(pp-2.)/2. +
     .              rwmodtab(mmod,lim1,ic+2)*pp*(pp*pp-1.)/6.
               ew = 10**rw*wave1(lim1)
               go to 10
            endif
         enddo

c*****add this EW to the total, weighting it by flux*radius^2*relcount
10       ewweight = ew*weightmod(mmod,lim1)
         weighttot = weighttot + weightmod(mmod,lim1)
         ewweighttot = ewweighttot + ewweight
         if (dabs(rwlgerror) .lt. 0.01) 
     .       write (nf7out,1005) fmodinput(mmod), 
     .       fmodoutput(mmod), radius(mmod), relcount(mmod), 
     .       fluxmod(mmod,lim1), weightmod(mmod,lim1), 
     .       1000.*ew, 1000.*ewweight
      enddo


c*****write out the mean EW
      ewmod(lim1) = ewweighttot/weighttot
      abundout(lim1) =dlog10(xabund(iatom))+deltangf+12.
      if (dabs(rwlgerror) .lt. 0.01) then
         write(nf7out,1006) atom1(lim1), wave1(lim1), e(lim1,1),
     .                       dlog10(gf(lim1)), 1000.*ewmod(lim1), 
     .                       abundout(lim1)
      endif
c      write (nf7out,1006) 1000.*ewmod(lim1), abundout(lim1)
      return


c*****format statements
1001  format ('MODFILE', 5x, 'COGOUT', 9x, 'RADIUS', 3x, '*COUNT',
     .        5x, 'FLUX', 3x, 'WEIGHT', 5x, 'EW', ' EWWEIGHT')
1005  format (2a12, 1pe9.2, e9.2, e9.2, e9.2, 0pf7.1, 1pe9.2)
c1006  format ('EWmean, Abundance =', f8.1, f8.2)
1006  format ('Final: ', f5.1, 3x, f8.3, 3x, f6.3, 3x,f7.3,1x,f8.1,
     .        1x,f8.2)


      end













