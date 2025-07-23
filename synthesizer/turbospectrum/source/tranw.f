      SUBROUTINE TRANW(NTAU,TAU,XMY,BPLAN,X,XL,DELI)
*
*-----------------------------------------------------------------------
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      include 'spectrum.inc'
*
      DIMENSION TAU(ndp),BPLAN(ndp),X(ndp),XL(ndp),ARG(ndp)
*
      P=0.0
      R=0.0
      DO 1 K=2,NTAU
        DIF=(TAU(K)-TAU(K-1))/XMY
        P=(X(K)+X(K-1))*0.5*DIF+P
        R=(XL(K)+XL(K-1))*0.5*DIF+R
        ARG(K)=BPLAN(K)*EXP(-P)*((X(K)+XL(K))*(R-0.5*R*R)-XL(K))
    1 CONTINUE
      Q=0.0
      DO 2 K=2,NTAU
        DIF=(TAU(K)-TAU(K-1))/XMY
        Y=DIF*(ARG(K)+ARG(K-1))*0.5
        Q=Q+Y
    2 CONTINUE
      DELI=Q
*
      RETURN
      END
