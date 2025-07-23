c******************************************************************************
c     this common block carries the data for population syntheses,
c     using multiple model calculations
c******************************************************************************


      real*8            gfmodtab(100,1000,50), rwmodtab(100,1000,50)
      real*8            weightmod(100,1000), fluxmod(100,1000),
     .                  ewmod(1000)
      real*8            radius(100), relcount(100)
      real*8            elspecial(10), abspecial(100,10)
      real*8            isospecial(10), fracspecial(100,10)
      real*8            deltangf, rwlgerror
      integer           nmodcurve(100,1000), modtot, nabs, nisos

      character*80      fmodinput(100), fmodoutput(100),
     .                  modpre, synpre


      common/multidata/ gfmodtab, rwmodtab,
     .                  weightmod, fluxmod,
     .                  ewmod,
     .                  radius, relcount,
     .                  elspecial, abspecial,
     .                  isospecial, fracspecial,
     .                  deltangf, rwlgerror,
     .                  nmodcurve, modtot, nabs, nisos

      common/multichar/ fmodinput, fmodoutput,
     .                  modpre, synpre


