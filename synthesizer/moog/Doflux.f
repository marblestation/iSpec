
      subroutine doflux
c******************************************************************************
c     This routine produces flux curves
c******************************************************************************
 
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'
      include 'Pstuff.com'


c*****examine the parameter file
      call params

c*****open the files for standard output and summary curves-of-growth
      nf1out = 20
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21
      lscreen = 6
      array = 'SUMMARY FLUX OUTPUT'
      nchars = 19
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)
      nf5out = 26
      lscreen = lscreen + 2
      array = 'POSTSCRIPT PLOT OUTPUT'
      nchars = 22
      call infile ('output ',nf5out,'formatted  ',0,nchars,
     .             f5out,lscreen)


c*****open and read the model atmosphere
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****compute the flux curve
      wave = start
1     call opacit (2,wave)
      if (modprintopt .ge. 2) 
     .    write(nf1out,1002) wave,(kaplam(i),i=1,ntau)
      call cdcalc (1)
      first = 0.4343*cd(1)
      flux = rinteg(xref,cd,dummy1,ntau,first)
      if (flux .le. 0.1) flux = 0.
      if (iunits .eq. 1) then
         write (nf1out,1003) 1.d-4*wave,flux
      else
         write (nf1out,1004) wave,flux
      endif
      waveinv = 1.0d4/wave
      if (flux .gt. 0.) then
         fluxlog = dlog10(flux)
      else
         fluxlog = -1.0
      endif
      write (nf2out,1001) wave, flux, waveinv, fluxlog
      wave = wave + step
      if (wave .le. sstop) go to 1
      call pltflux


c****end the computations
      call finish (0)
      return


c*****format statements
1001  format (1p2d12.4,0p2f10.4)
1002  format ('  kaplam from 1 to ntau at wavelength',f11.3/
     1        (6(1pd12.4)))
1003  format ('AT WAVELENGTH/FREQUENCY =',f11.7,
     .        '  CONTINUUM FLUX/INTENSITY =', 1pd12.5)
1004  format ('AT WAVELENGTH/FREQUENCY =',f11.3,
     .        '  CONTINUUM FLUX/INTENSITY =', 1pd12.5)


      end                                                               


