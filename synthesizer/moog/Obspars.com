
c******************************************************************************
c     this common block is simply to share information between the
c     observation input routines
c******************************************************************************

      real*8 disp(9), bzero, bscale 
      integer ibits, nblock, naxis, naxis1, byteswap
      
      common/obspars/disp, bzero, bscale, 
     .               ibits, nblock, naxis, naxis1, byteswap

