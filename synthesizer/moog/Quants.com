
c******************************************************************************
c     this common block holds the internally-stored atomic data
c******************************************************************************

      real*8        xsolar(95), xam(95), newpartdata(50,6), 
     .              xchi1(95), xchi2(95), xchi3(95) 
      integer       nudata(6,380), partflag(95,4), nu 

      common/quants/xsolar, xam, newpartdata,
     .              xchi1, xchi2, xchi3,
     .              nudata, partflag, nu       

