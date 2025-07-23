
c******************************************************************************
c     this simply is a "scratch" common block that holds arrays
c     that are temporary and will be "recycled" for different
c     internal calculation purposes.
c******************************************************************************

      real*8       dummy1(100), dummy2(100), dummy3(5000), dummy4(5000)

      common/dummy/dummy1, dummy2, dummy3, dummy4

