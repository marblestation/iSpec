
      subroutine lineabund (abundin)
c******************************************************************************
c     This routine iteratively determines the abundance for a single line;
c     for ease and speed of computation the actual changes are done to 
c     gf-values, and then at the end the adjustment is made to input
c     abundance by the amount that the gf had to be changed
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'

  
      lim2 = lim1
      ncurve = 1


c*****find the c-o-g gf that matches the observed RW
      gf1(ncurve) = gf(lim1)
      rwlgobs = dlog10(width(lim1)/wave1(lim1))
      gfobs = gftab(ntabtot)
      do i=2,ntabtot
         if (rwtab(i) .gt. rwlgobs) then
            gfobs = gftab(i-1) + (gftab(i)-gftab(i-1))*
     .              (rwlgobs-rwtab(i-1))/(rwtab(i)-rwtab(i-1))
            go to 15
         endif
      enddo


c*****using the current gf, compute a line and find the c-o-g gf that matches
c*****the computed RW
15    call oneline (1)
      rwlgcal = dlog10(w(ncurve)/wave1(lim1))
      gfcal = gftab(ntabtot)
      do i=2,ntabtot
         if (rwtab(i) .gt. rwlgcal) then
            gfcal = gftab(i-1) + (gftab(i)-gftab(i-1))* 
     .              (rwlgcal-rwtab(i-1))/(rwtab(i)-rwtab(i-1))
            go to 20
         endif
      enddo


c*****are the observed and computed RWs too different?  if so, iterate 
c*****by adjusting the assumed gf, and recompute the line.
c*****for strong lines, the iterations are slowed down by using the
c*****square root of the proposed gf shift.
20    error = (w(ncurve)-width(lim1))/width(lim1)
      ratio = 10.**(gfobs-gfcal)
      ncurve = ncurve + 1
      if (dabs(error) .ge. 0.0015 .and. ncurve .lt. 20) then
         rwlcomp = dlog10(w(ncurve-1)/wave1(lim1))
         if (rwlcomp .gt. -4.7) then
            gf1(ncurve) = gf1(ncurve-1)*dsqrt(ratio)
            do i=1,ntau                                
               kapnu0(lim1,i) = kapnu0(lim1,i)*dsqrt(ratio)
            enddo
         else
            gf1(ncurve) = gf1(ncurve-1)*ratio
            do i=1,ntau                                
               kapnu0(lim1,i) = kapnu0(lim1,i)*ratio            
            enddo
         endif
         go to 15
      endif


c*****if the observed and computed RWs are close, do a final gf adjustment
c*****and finish with one more line recomputation
      if (ncurve .eq. 30) then
         write (nf1out,1001)
         write (nf2out,1001)
      endif
      gf1(ncurve) = gf1(ncurve-1)*ratio
      do i=1,ntau
         kapnu0(lim1,i) = kapnu0(lim1,i)*ratio            
      enddo
      call oneline (2)
      widout(lim1) = w(ncurve)
      wid1comp(lim1) = w(1)
      diff = dlog10(gf1(ncurve)/gf(lim1))
      abundout(lim1) = abundin + diff
      if (ncurve .ne. 1) then
         write (nf1out,1002) ncurve
      endif


      return


c*****format statements
1001  format ('OH NO! ANOTHER FAILED ITERATION!')
1002  format (' This fit required ',i2,' iterations including ',
     .        'a final small adjustment'/)


      end


