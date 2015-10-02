
      subroutine curve
c******************************************************************************
c     This routine produces a curve-of-growth for a single line          
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'


c*****set up the parameters 
      ewsynthopt = -1
      ncurve = 1
      dec = 10.**(rwstep)
      wstart = 10.**(rwlow)*wave1(lim1)
      wstop = 10.**(rwhigh)*wave1(lim1)
      gf1(ncurve) = gf(lim1)

 
c*****increment the log(gf) value backward by rwstep and redo the calculations 
c     until the end (rwlow) is reached
31    call oneline (2)
      if (w(ncurve) .gt. wstart) then
         gf1(ncurve) = gf1(ncurve)/dec
         do i=1,ntau
            kapnu0(lim1,i) = kapnu0(lim1,i)/dec
         enddo
         go to 31
      endif
      go to 61
 

c*****increment the log(gf) value forward by rwstep and redo the calculations 
c     until the end (rwhigh) is reached
60    gf1(ncurve) = gf1(ncurve-1)*dec
      do i=1,ntau                                                     
         kapnu0(lim1,i) = kapnu0(lim1,i)*dec
      enddo
61    call oneline (2)
      if (w(ncurve) .lt. wstop) then
         ncurve = ncurve + 1
         go to 60
      endif
 

c*****end the computations with a summary print
      if (nf2out .ne. 0 .and. lim1 .eq. 1)
     .   write (nf2out,1001) moditle
      do i=1,ncurve
         w(i) = dlog10(w(i)/wave1(lim1))
         gf1(i) = dlog10(gf1(i))
      enddo
      iatom = atom1(lim1)
      if(iatom .ge. 100) iatom = 1
      abund = dlog10(xabund(iatom)) + 12.
      write (nf1out,1002) wave1(lim1),atom1(lim1),e(lim1,1),
     .                    abund,ncurve
      write (nf1out,1003) (gf1(i),w(i),i=1,ncurve)
      if (nf2out .eq. 0) return
      write (nf2out,1002) wave1(lim1),atom1(lim1),e(lim1,1),
     .       abund,ncurve
      write (nf2out,1003) (gf1(i),w(i),i=1,ncurve)
      return


c*****format statements
1001  format (a80)
1002  format(/'wavelength =', f9.3,5x, 'species =', f6.1,5x, 'ep =',
     .       f7.3, 'abundance =', f10.3, 5x, 'n =',i2)
1003  format('  curve of growth in (loggf,logrw) pairs'/
     .       (5(f7.3,',',f7.3)))
1004  format (10f7.3)


      end



