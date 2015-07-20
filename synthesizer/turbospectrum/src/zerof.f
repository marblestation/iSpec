      SUBROUTINE ZEROF(F,DX,DFDX)
*
*-----------------------------------------------------------------------
*
* FIND DX=-F/DFDX, TO MAKE F ZERO. IF DX=0 AT ENTRY, THEN DFDX IS A
* START APPROXIMATION, AND IT IS THE FIRST CALL. OTHERWISE USE OLD
* INFO.  780926/NORDLUND.
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      IF(DX.NE.0.) DFDX=(F-FOLD)/DXOLD
      DX=-F/DFDX
      FOLD=F
      DXOLD=DX
*
      RETURN
      END
