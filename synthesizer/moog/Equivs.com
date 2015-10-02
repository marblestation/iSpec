
c******************************************************************************
c     this is not really a common block, but it is a way to make
c     (in a consistent manner for several routines) definitions
c     of some temporary variables that are equivalenced (assigned,
c     or in other words, use the same memory space as) other variables. 
c     This is done so that plotting of large numbers of points (<=50000)
c     for synthetic spectra may be accomplished without allocating
c     large amounts of additional memory for the code.
c******************************************************************************

      real*4      chunk(100000,5)
      equivalence (chunk(1,1),a(1,1)),  
     .            (chunk(1,2),a(1,21)),
     .            (chunk(1,3),a(1,41)), 
     .            (chunk(1,4),a(1,61)),
     .            (chunk(1,5),a(1,81))

      real*4      xsyn(100000), dev(100000)
      equivalence (xsyn,kapnu0(1,1)), (dev,kapnu0(1,21))

      real*4      y(100000), z(100000)
      equivalence (y,kapnu0(1,41)), (z,kapnu0(1,61))

