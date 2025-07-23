      FUNCTION BPL(T,X)
*
*-----------------------------------------------------------------------
*
* Planc function
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      COMMON/BPLC/ EX,X5
*
      DATA CP/1.191044E27/,C2/1.438769E8/
*
      X5=((X**2)**2)*(X/CP)
      EX=EXP(C2/(T*X))
      BPL=1./((EX-1.)*X5)
*
      RETURN
      END
