
      subroutine binplotprep
c******************************************************************************
c     This routine concatenates two syntheses in order to make a smoothed
c     binary synthetic spectrum
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      include 'Multistar.com'
      real*8 normalization
      integer pshift


c*****open the files 
      nf7out = 21
      lscreen = 4
      array = 'RAW SYNTHESIS INPUT FOR PRIMARY'
      nchars = 31
      call infile ('input  ',nf7out,'formatted  ',0,nchars,
     .             f7out,lscreen)
      nf8out = 22
      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS INPUT FOR SECONDARY'
      nchars = 33
      call infile ('input  ',nf8out,'formatted  ',0,nchars,
     .             f8out,lscreen)
      nf9out = 23
      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS OUTPUT FOR COMBINED BINARY'
      nchars = 40
      call infile ('output ',nf9out,'formatted  ',0,nchars,
     .             f9out,lscreen)
      if (plotopt .ne. 0) then
         nf10out = 24
         lscreen = lscreen + 2
         array = 'SMOOTHED SYNTHESIS OUTPUT FOR COMBINED BINARY'
         nchars = 45
         call infile ('output ',nf10out,'formatted  ',0,nchars,
     .                f10out,lscreen)
         nf5out = 26
         lscreen = lscreen + 2
         array = 'POSTSCRIPT PLOT OUTPUT FOR COMBINED BINARY'
         nchars = 22
         call infile ('output ',nf5out,'formatted  ',0,nchars,
     .                f5out,lscreen)
      endif
      if (plotopt .gt. 1) then
         nf6out = 27
         lscreen = lscreen + 2
         array = 'SPECTRUM COMPARISON OUTPUT'
         nchars = 27
         call infile ('output ',nf6out,'formatted  ',0,nchars,
     .                f6out,lscreen)
      endif


c*****compute and output the parameters related to the spectrum combination:
c     the declared velocity difference, the luminosity ratio, and
c     for information, the computed fluxes of the two stars; get the user's
c     agreement
      deltawave = deltaradvel/2.9979d+5*(start+sstop)/2.
      pshift = nint(deltawave/step)
      write (array,1002) fluxprimary, fluxsecondary
      lscreen = lscreen + 2
      nchars = 55
      call putasci (nchars,lscreen)
      write (array,1003) lumratio
      lscreen = lscreen + 1
      nchars = 49
      call putasci (nchars,lscreen)
      write (array,1001) deltaradvel, deltawave, pshift
      lscreen = lscreen + 1
      nchars = 71
      call putasci (nchars,lscreen)
      write (array,*) 'ARE THESE VALUES OK (y/n)? '
      lscreen = lscreen + 2
      nchars = 27


c*****read back the header information from the individual star raw 
c     synthetic spectra
99    read (nf7out,1004,end=100) array
      read (nf8out,1004) chinfo
      if     (array(1:7).eq.'Isotopi') then
         write (nf9out,1004) array
         go to 99
      elseif (array(1:7).eq.'ALL abu') then
         write (nf9out,1010) array(1:59), chinfo(54:59)
         go to 99
      elseif (array(1:7).eq.'Changin') then
         write (nf9out,1011) array(1:37), chinfo(31:37)
         go to 99
      elseif (array(1:7).eq.'element') then
         write (nf9out,1012) array(1:30), chinfo(25:30)
         go to 99
      elseif (array(1:7).eq.'MODEL: ') then 
         write (nf9out,1013) array(8:43), chinfo(8:43)
         modbin(1) = array(1:80)
         modbin(2) = chinfo(1:80)
      endif
      read (nf7out,1004) array
      read (nf8out,1004) chinfo
      read (array,*) start, sstop, step
      kount = nint((sstop - start + (step/4.0) )/step) + 1
      write (nf9out,1004) array


c*****read back the raw synthetic spectra from the individual stars;
c     shift by the appropriate point number, add the spectra (normalized),
c     and dump to combined raw synthetic spectrum file
      read (nf7out,1006) (y(i),i=1,kount)
      read (nf8out,1006) (z(i),i=1,kount)
      do i=1,kount
         y(i) = 1.0 - y(i)
         z(i) = (1.0 - z(i))/lumratio
      enddo
      normalization = 1.0 + 1.0/lumratio
      if (pshift .gt. 0) then
         do i=1,pshift
            dev(i) = 1.0
         enddo
         do i=pshift+1,kount
            dev(i) = (y(i)+z(i-pshift))/normalization
         enddo
      elseif (pshift .eq. 0) then
         do i=1,kount
            dev(i) = (y(i)+z(i))/normalization
         enddo
      else
         do i=1, kount-pshift
            dev(i) = (y(i)+z(i-pshift))/normalization
         enddo
         do i=kount-pshift+1,kount
            dev(i) = 1.0
         enddo
      endif
      do i=1,kount
         y(i) = 1.0 - dev(i)
      enddo
      write (nf9out,1006) (y(i),i=1,kount)
      go to 99
100   return


c*****format statements
1001  format ('INPUT DELTA VELOCITY = ', f8.3,
     .        ';  DELTA WAVELENGTH = ', f8.3,
     .        ';  POINT SHIFT = ', i4)
1002  format ('PRIMARY FLUX = ', 1pe10.3, 
     .        ';  SECONDARY FLUX = ', 1pe10.3)
1003  format ('INPUT LUMINOSITY RATIO = ', 0pf8.3)
1004  format (a80)
1006  format (10f7.4)
1010  format (a59, ',', a6, ' dex')
1011  format (a37, ',', a6, ' dex')
1012  format (a30, ',', a6)
1013  format ('MODEL: ', 2a36)


      end                                                            






