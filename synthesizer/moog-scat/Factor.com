
c******************************************************************************
c     this common block is used for multiple syntheses, and contains
c     abundance information that MOOG uses to figure out what will be 
c     altered for the different syntheses.
c******************************************************************************

      real*8         pecabund(95,5), newpecabund(95,5), abfactor(5)
      integer        pec(95), newpec(95), 
     .               numpecatom, newnumpecatom,
     .               numatomsyn, newnumatomsyn,   
     .               isynth, ninetynineflag
      real*8         isotope(20), newisotope(20),
     .               isoabund(20,5), newisoabund(20,5)
      integer        numiso, newnumiso,
     .               numisosyn, newnumisosyn, isorun

      common/newfact/pecabund, newpecabund, abfactor,
     .               pec, newpec,
     .               numpecatom, newnumpecatom, 
     .               numatomsyn, newnumatomsyn, 
     .               isynth, ninetynineflag
      common/iso/    isotope, newisotope,
     .               isoabund, newisoabund,
     .               numiso, newnumiso,
     .               numisosyn, newnumisosyn, isorun

