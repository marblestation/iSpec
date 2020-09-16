
c******************************************************************************
c     this common block is used for storning information that will be
c     used in combining spectra of different stars in various ways
c******************************************************************************

      real*8          binpecabund(2,95,5), binpec(2,95),
     .                fluxprimary, fluxsecondary,
     .                lumratio, deltaradvel
      integer         begin
      character*80    modbin(2)

      common/binstuff/
     .                binpecabund, binpec,
     .                fluxprimary, fluxsecondary,
     .                lumratio, deltaradvel,
     .                begin
      common/binchars/
     .                modbin

