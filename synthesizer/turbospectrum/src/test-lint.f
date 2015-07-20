      program test_lint


      doubleprecision x(10),xint
      real y(10),yint

      data x/1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
      data y/2.,3.,4.,5.,6.,5.,4.,3.,2.,1./

      xint=0.
      call lint(10,x,y,xint,yint)
      print*,xint,yint
      do i=1,10
      do j=1,3
      xint=x(i)+float(j)/3.
      call lint(10,x,y,xint,yint)
      print*,xint,yint
      enddo
      enddo

      end

      SUBROUTINE LINT(N,X,Y,XINT,YINT)
*
*-----------------------------------------------------------------------
*
* THIS ONE DIMENSIONAL INTERPOLATION ROUTINE WORKS WITH SUCCESIVE
* SPLITTING IN HALVES.
* IN PRINCIPLE IT INTERPOLATES AND EXTRAPOLATES WILLINGLY.
* A WARNING IS GIVEN AT EXTRAPOLATION. THE INTERPOLATIONS ARE DONE
* linearly
* ***** OBSERVE. TABLES SHOULD BE  G R O W I N G *****
* SPECIAL VERSION FOR THE CANARIAS
*
* BPz 30/08-1999 from TINT.
*-----------------------------------------------------------------------
*
      doubleprecision x,xint,arg
      DIMENSION X(N),Y(N),ARG(2),FUNK(2)
*
      COMMON/UTPUT/IREAD,IWRIT
*
      IF(N.EQ.2) then
        arg(1)=x(1)
        arg(2)=x(2)
        funk(1)=y(1)
        funk(2)=y(2)
        GOTO 21
      endif
*
      NOEV=N
      NED=1
    1 NP=(NED+NOEV)/2
      IF(XINT-X(NP)) 2,2,3
    2 NOEV=NP
      GOTO 4
    3 NED=NP
    4 IF(NOEV-NED-1) 5,5,1
    5 IF(NOEV-1) 7,7,6
    6 J=NOEV-2
      GOTO 8
    7 J=NOEV-1
    8 DO 9 K=1,2
        JP=K+J
        ARG(K)=X(JP)
    9 FUNK(K)=Y(JP)
      IF(ARG(1)-XINT) 11,11,12
   11 IF(ARG(2)-XINT) 12,21,21
   12 WRITE(IWRIT,200) XINT,ARG
  200 FORMAT(' WARNING, EXTRAPOL. IN LINT. XINT=',E17.8,3X,/,
     &       ' ARG=',2E17.8)
*
   21 YINT=(funk(2)-funk(1))/(arg(2)-arg(1))*(XINT-arg(1))+funk(1)
      RETURN
*
      END
