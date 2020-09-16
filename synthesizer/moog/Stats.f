
      subroutine stats
c******************************************************************************
c     This routine does abundance statistics for a single species
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'

c*****compute the average and standard deviation
      average = 0.
      kount = 0
      do l=lim1obs,lim2obs
         if (abundout(l) .ne. 999.99) then
            average = average + abundout(l)
            kount = kount + 1
         endif
      enddo
      if (kount .gt. 0) then
         average = average/kount
      endif
      deviate = 0.
      if (kount .gt. 1) then
         do l=lim1obs,lim2obs
            if (abundout(l) .ne. 999.99) then
               deviate = deviate + (abundout(l)-average)**2
            endif
         enddo
         deviate = dsqrt(deviate/(kount-1))
      endif

         
c*****correlate the abundances with excitation potential, equivalent width,
c     and wavelength
      if (kount .gt. 2) then
         epmin = 999.
         epmax = -999.
         rwmin = 999.
         rwmax = -999.
         wvmin = 999999.
         wvmax = -999999.
         x1 = 0.
         x2 = 0.
         x3 = 0.
         x4 = 0.
         x5 = 0.
         x6 = 0.
         y1 = 0.
         y2 = 0.
         xy = 0.
         yz = 0.
         za = 0.

         do l=lim1obs,lim2obs
            if (abundout(l) .ne. 999.99) then
c              rw = dlog10(wid1comp(l)/wave1(l))
               rw = dlog10(width(l)/wave1(l))
               x1 = x1 + e(l,1)
               x2 = x2 + e(l,1)**2
               x3 = x3 + rw
               x4 = x4 + rw**2
               x5 = x5 + wave1(l)
               x6 = x6 + wave1(l)**2
               y1 = y1 + abundout(l)
               y2 = y2 + abundout(l)**2
               xy = xy + e(l,1)*abundout(l)
               yz = yz + rw*abundout(l)
               za = za + wave1(l)*abundout(l)
               if (e(l,1) .lt. epmin) epmin = e(l,1)
               if (e(l,1) .gt. epmax) epmax = e(l,1)
               if (rw .lt. rwmin) rwmin = rw
               if (rw .gt. rwmax) rwmax = rw
               if (wave1(l) .lt. wvmin) wvmin = wave1(l)
               if (wave1(l) .gt. wvmax) wvmax = wave1(l)
            endif
         enddo

         xxm1 =    (kount*xy - x1*y1)/(kount*x2 - x1**2)
         xxb1 =    (y1*x2 - xy*x1)/(kount*x2 - x1**2)
         xxr1 =    (kount*xy - x1*y1)/dsqrt((kount*x2 - x1**2)*
     .             (kount*y2 - y1**2))
         deltaep = epmax - epmin

         xxm2 =    (kount*yz - x3*y1)/(kount*x4 - x3**2)
         xxb2 =    (y1*x4 - yz*x3)/(kount*x4 - x3**2)
         xxr2 =    (kount*yz - x3*y1)/dsqrt((kount*x4 - x3**2)*
     .             (kount*y2 - y1**2))
         deltarw = rwmax - rwmin

         xxm3 =    (kount*za - x5*y1)/(kount*x6 - x5**2)
         xxb3 =    (y1*x6 - za*x5)/(kount*x6 - x5**2)
         xxr3 =    (kount*za - x5*y1)/dsqrt((kount*x6 - x5**2)*
     .             (kount*y2 - y1**2))
         deltawv = wvmax - wvmin

      endif

      return
      end
  

