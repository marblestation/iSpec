
      subroutine vargauss (line)
c******************************************************************************
c     This subroutine prepares synthesis data for plotting by smoothing
c     the data with a Gaussian whose FWHM has been specified by the user
c     at various wavelengths along the spectrum
c******************************************************************************
 
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      real*8 wavefwhm(100),fwhm(100)
      character*80 fsmooth
      data istart /0/


c*****initialize parameters
      write (abitle,1001) 
      nsyn = 1
      sigma = 0.   
      jrotdel = 0
      fsmooth = 'no_filename_given'


c*****on entering, rewind the unsmoothed and smoothed spectrum output 
c     files, and get the synthesis range parameters from the 'dump' file
      rewind nf2out
      rewind nf3out
55    read (nf2out,1002) moditle
      if (moditle(1:15).eq.'        element' .or.  
     .    moditle(1:15).eq.' ALL abundances' .or.
     .    moditle(1:15).eq.'Isotope Ratio: ') go to 55
      read (nf2out,*) start, sstop, step
      kount = nint((sstop - start + (step/4.0) )/step) + 1
      rewind nf2out


c*****the first time through, read in the Gaussian FWHM array
      if (istart .eq. 0) then
         istart = 1
         nfsmooth = 35
         array = 'SMOOTHING FWHM DATA'
         nchars = 19
         call infile ('input  ',nfsmooth,'formatted  ',0,nchars,
     .                fsmooth,line)
         j = 1
39       read (nfsmooth,*,end=40) wavefwhm(j), fwhm(j)
         j = j + 1
         if (j .le. 100) then
            go to 39
         else
            istat = ivcleof (line,1)
            istat = ivmove (line-1,1)
            write (*,*) 'ERROR: TOO MANY POINTS IN THE FWHM FILE!'
            smtype = 'e'
            return
         endif
40       jtotfwhm = j - 1
      endif


c*****now read in the raw spectrum and flip to a depth scale
45    noff = 80*(nsyn-1)
      abitle(noff+1:noff+12) = '[M/H] 0.00  '
      nabunds = 0
41    read (nf2out,1002,end=2000) array
      if     (array(1:15).eq.' ALL abundances') then
         abitle(noff+6:noff+10) = array(56:60)
         go to 41
      elseif (array(1:15).eq.'        element') then
         nabunds = nabunds + 1
         if (nabunds .le. 7) then
            ioff = noff + 12 + 9*(nabunds-1)
            abitle(ioff+1:ioff+2) = array(17:18)
            abitle(ioff+3:ioff+7) = array(34:38)
            abitle(ioff+8:ioff+9) = '  '
         endif
         go to 41
      elseif (array(1:15).eq.'Isotope Ratio: ') then
         nabunds = nabunds + 1
         if (nabunds .le. 7) then
            ioff = noff + 12 + 9*(nabunds-1)
            abitle(ioff+1:ioff+4) = array(27:30)
            abitle(ioff+5:ioff+5) = ' '
            do k=33,44
               if (array(k:k) .ne. ' ') then
                  abitle(ioff+6:ioff+8) = array(k:k+2)
                  go to 60
               endif
            enddo
60          abitle(ioff+9:ioff+9) = ' '
         endif
         go to 41
      endif
      read (nf2out,1002,end=2000)
      read (nf2out,1003,end=2000) (y(i),i=1,kount)
      do i=1,kount
         y(i) = 1.0 - y(i)
      enddo


c*****go through the spectrum wavelengths step by step; the 
c     Gaussian smoothing will need to be different at each step
      i = 0
      oldhalf = 0.
25    i = i + 1
      if (i .gt. kount) go to 90
      synstep = start + (i-1)*step


c*****interpolate linearly in the FWHM array to get the appropriate value 
c     for the current wavelength step
      if     (synstep .le. wavefwhm(1)) then
         half = fwhm(1)
      elseif (synstep .ge. wavefwhm(jtotfwhm)) then
         half = fwhm(jtotfwhm)
      else
         do j=2,jtotfwhm
            if (synstep .le. wavefwhm(j)) then
               half = fwhm(j-1) + (synstep-wavefwhm(j-1))*
     .                (fwhm(j)-fwhm(j-1))/(wavefwhm(j)-wavefwhm(j-1))
               if (half .gt. 0.) then
                  go to 10
               else
                  go to 15
               endif
            endif
         enddo
      endif


c*****compute the Gaussian smoothing function, if needed
10    if (dabs(half-oldhalf)/half .lt. 0.03) go to 50
      oldhalf = half
      sigma = half/2
      aa = 0.6932/sigma**2
      power = 1.0
      do k=1,1000
         p(k) = dexp(-aa*(step*k)**2 )
         power = power + 2*p(k)
         if (p(k) .lt. 0.05) then
            jdel = k
            min = jdel + 1
            max = kount - jdel
            p0 = 1.
            go to 50 
         endif
      enddo
      istat = ivcleof (line,1)
      istat = ivmove (line-1,1)
      write (*,1023) sigma, (p(i),i=1,1000)
      smtype = 'e'
      return


c*****if no smoothing, just equate the smoothed to the unsmoothed point
15    z(i) = y(i)
      go to 25

  
c*****otherwise smooth the spectrum
50    if (i.lt.min .or. i.gt.max) then
         z(i) = y(i)
      else
         z(i) = p0*y(i)
         do k=1,jdel
            z(i) = z(i) + p(k)*(y(i-k) + y(i+k))
         enddo
         z(i) = z(i)/power
      endif
      go to 25
 

c*****copy the smoothed spectrum to the appropriate array
90    do i=1,kount
         chunk(i,nsyn) = z(i)
      enddo


c*****compute the wavelength array; must be done for each synthetic
c     spectrum because of the way the equivalences were set up
      do i=1,kount
         xsyn(i) = start + (i-1)*step
      end do
 

c*****dump the smoothed spectrum in a MONGO-style set of 
c     (wavelength,flux) point pairs
      write (nf3out,1005) kount, start, sstop, step
      if (xsyn(1) .le. 100.0) then
         write (nf3out,1009) (xsyn(i),z(i),i=1,kount)
      else 
         write (nf3out,1008) (xsyn(i),z(i),i=1,kount)
      endif
      nsyn = nsyn + 1
      go to 45


c*****exit the routine normally
2000  nsyn = nsyn - 1
      return


c*****format statements
1001  format (400(' '))
1002  format (a80)
1003  format (10f7.4)
1005  format ('the number of points per synthesis = ',i5/
     .       'start = ',f10.3,5x,'stop = ',f10.3,5x,'step = ',f10.3)
1008  format (f10.3,'  ',f10.5)
1009  format (f10.6,'  ',f10.5)
1023  format ('ERROR: GAUSSIAN PROFILE TOO BIG! (HALF WIDTH=',
     .       f6.2,')      THE PROFILE:'/(15f5.2))

      end



