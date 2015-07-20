      SUBROUTINE INP3(X,Y,XINT,YINT)
*
*-----------------------------------------------------------------------
*
* NEWTONINTERPOLATION, TREPUNKTS
* OBS *** INGEN SPAERR MOT EXTRAPOLATION *****
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      DIMENSION X(3),Y(3),F1(2)
*
      DO 1 K=1,2
        F1(K)=(Y(K+1)-Y(K))/(X(K+1)-X(K))
    1 CONTINUE
*      PRINT 100,K,X(K+1),X(K),X(3),X(1)
*  100 FORMAT(1X,I2,4F12.5)
      F2=(F1(2)-F1(1))/(X(3)-X(1))
      YINT=Y(1)+(XINT-X(1))*F1(1)+(XINT-X(1))*(XINT-X(2))*F2
*
      RETURN
      END
