
c******************************************************************************
c     this common block holds the internally-stored arrays for continuous
c     opacity calculations
c******************************************************************************

      real*8       aH1(100), aHminus(100)
      real*8       aHeminus(100)
      real*8       sigel(100), sigH(100), sigH2(100), sigHe(100)
      real*8       aC1(100), aMg1(100), aMg2(100), aAl1(100), aSi1(100), 
     .             aSi2(100), aFe1(100)
      real*8       evhkt(100)
      real*8       freq, freqlg

      common/kappa/aH1, aHminus, 
     .             aHeminus, 
     .             sigel, sigH, sigH2, sigHe,
     .             aC1, aMg1, aMg2, aAl1, aSi1, 
     .             aSi2, aFe1,
     .             evhkt,
     .             freq, freqlg


