c******************************************************************************
c     this common block carries the data to and from the molecular
c     equilibrium calculations; it is also used to hold isotopic
c     abundance information
c******************************************************************************

c     amol = names of the species
c     smallmolist = the small set of default molecule names
c     largemolist = the large set of default molecule names
c     xmol     = number density of the species at all atmopshere layers
c     iorder   = names of elements involved in molecular equilibrium,
c                in the order that they are encountered in the species list
c     xatom    = working array at a particular layer:  the number densities 
c                of neutral atomic species
c     patom    = working array at a particular layer:  the partial pressures 
c                of neutral atomic species
c     pmol     = working array at a particular layer:  the partial pressures
c                of molecules
c     xamol    = number densities of neutral atomic species at all layers
c     const    = molecular constants loaded in from Bmolec.
c     h2ocoeff = polynomial coefficients for H_2O partition functions
c     co2coeff = polynomial coefficients for H_2O partition functions


      real*8       pmol(110), xmol(110,100), xamol(30,100),
     .             xatom(30), patom(30),
     .             amol(110), smallmollist(110), largemollist(110),
     .             datmol(7,110), const(6,110), 
     .             h2ocoeff(5), co2coeff(5), xnh2o(100), xnco2(100),
     .             uh2o(100), uco2(100)
      integer      neq, lev, nmol, natoms, iorder(30), molopt, molset


      common /mol/ pmol, xmol, xamol,
     .             xatom, patom, 
     .             amol, smallmollist, largemollist,
     .             datmol, const,
     .             h2ocoeff, co2coeff, xnh2o, xnco2,
     .             uh2o, uco2,
     .             neq, lev, nmol, natoms, iorder, molopt, molset

