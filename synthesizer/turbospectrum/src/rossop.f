      FUNCTION ROSSOP(T,PE,k)
*
*-----------------------------------------------------------------------
*
* 'ROSSOP' CALCULATES THE ROSSELAND MEAN OPACIY AS DEFINED ON THE
* FIRST WAVELENGTH SCALE IN 'INITAB' INPUT.  *NORD*
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      include 'spectrum.inc'
      common/CA5/ ab(mkomp), fakt(mkomp),ppe(ndp),tt(ndp),xla(20),
     &            xla3(20),ro,sumabs,sumsca,viktr,iset,nlb
 
*
cc      DATA NEWT/2/
*
      newt=1
      CALL JON(T,PE,1,PG,RO,DUM,0,k)
      CALL ABSKO(NEWT,1,T,PE,1,1,RSP,DUM)
cc      NEWT=1
      ROSSOP=RSP+DUM
*
      RETURN
      END
