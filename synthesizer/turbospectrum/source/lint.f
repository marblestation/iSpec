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
*
* BPz 30/08-1999 from TINT.
*-----------------------------------------------------------------------
*
      integer countextrap
      doubleprecision x,xint,arg
      DIMENSION X(N),Y(N),ARG(2),FUNK(2)
*
      COMMON/UTPUT/IREAD,IWRIT
      data countextrap/0/
      save countextrap
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
      countextrap=countextrap+1
  200 FORMAT(' WARNING, EXTRAPOL. IN LINT. XINT=',E17.8,3X,/,
     &       ' ARG=',2E17.8)
*
   21 YINT=(funk(2)-funk(1))/(arg(2)-arg(1))*(XINT-arg(1))+funk(1)

      if (countextrap.gt.10) then
        print*,
     & 'I told you 10 times that the continous opacity is extrapolated'
        print*,
     & ' Now that''s enough. Obviously you are running bsyn with'
        print*,
     & ' a wavelength set different from that you used for babsma.'
        print*,
     & ' Rerun babsma with the wavelength set you are using for bsyn,'
        print*,' and try again.'
        stop 'in lint.f'
      endif

      RETURN
*
      END
