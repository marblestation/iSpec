      SUBROUTINE TINT(N,X,Y,XINT,YINT)
*
*-----------------------------------------------------------------------
*
* THIS ONE DIMENSIONAL INTERPOLATION ROUTINE WORKS WITH SUCCESIVE
* SPLITTING IN HALVES.
* IN PRINCIPLE IT INTERPOLATES AND EXTRAPOLATES WILLINGLY.
* A WARNING IS GIVEN AT EXTRAPOLATION. THE INTERPOLATIONS ARE DONE
* WITH A THREE POINT FORMULA, SUBR.  I N P 3 .
* ***** OBSERVE. TABELS SHOULD BE  G R O W I N G *****
* SPECIAL VERSION FOR THE CANARIAS
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      DIMENSION X(N),Y(N),ARG(3),FUNK(3)
*
      COMMON/UTPUT/IREAD,IWRIT
*
      IF(N.EQ.2) GOTO 21
*
      NOEV=N
      NED=1
    1 NP=(NED+NOEV)/2
      IF(XINT-X(NP)) 2,2,3
    2 NOEV=NP
      GOTO 4
    3 NED=NP
    4 IF(NOEV-NED-2) 5,5,1
    5 IF(NOEV-2) 7,7,6
    6 J=NOEV-3
      GOTO 8
    7 J=NOEV-2
    8 DO 9 K=1,3
        JP=K+J
        ARG(K)=X(JP)
    9 FUNK(K)=Y(JP)
      IF(ARG(1)-XINT) 11,11,12
   11 IF(ARG(3)-XINT) 12,13,13
   12 continue
cc      WRITE(IWRIT,200) XINT,ARG
  200 FORMAT(' WARNING, EXTRAPOL. IN TINT. XINT=',E17.8,3X,/,
     &       ' ARG=',3E17.8)
   13 CALL INP3(ARG,FUNK,XINT,YINT)
      RETURN
*
   21 YINT=(Y(2)-Y(1))/(X(2)-X(1))*(XINT-X(1))+Y(1)
      RETURN
*
      END
