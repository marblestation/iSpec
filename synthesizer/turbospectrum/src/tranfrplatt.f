      SUBROUTINE TRANFRplatt(ISTEP)
*
*-----------------------------------------------------------------------
*
* FORMAL SOLVES THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION 'SOURCE'.
* 'ERROR' IS THE RESULTING ERROR IN THE DEFINITION OF THE CONTINUUM
* SCATTERING SOURCE FUNCTION. TRANSFR CALCULATES THE MATRIX ELEMENTS
* OF THE PROBLEM. FLUX AND INTENSITIES AT TAU=0 ARE RETURNED IN /CSURF/.
* 79.06.21 *NORD*
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      include 'spectrum.inc'
*
      COMMON/CANGLE/ MMU,XMU(nrays),XMU2(nrays),H(nrays)
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),xh(NDP),XK(NDP),
     &  fillup(4*ndp+3)
      COMMON/TAUC/ TAU(ndp),DTAULN(ndp),JTAU
      COMMON/SPACE2/ SOURCE(ndp),ERROR(ndp),DUM(3*ndp),P(ndp),
     &               SP1(ndp,nrays),SP2(ndp,nrays),SP3(ndp,nrays),
     &               AN(ndp),AD(ndp),BD(ndp),
     &               FACT(ndp),DSO(ndp),C(nrays),T(nrays),EX(nrays),
     &               tomatch(4*nrays*ndp-ndp+1)
      COMMON/CSURF/ HSURF,Y1(nrays)
*
      IST=ISTEP
      IF(ISTEP.GT.0) GOTO 400
*
* MU LOOP
*
      JTAU1=JTAU-1
      JTAU2=JTAU-2
      DO 110 I=1,MMU
*
* K=1
*
      DTAUB=.5*(X(1)+S(1)+X(2)+S(2))*(TAU(2)-TAU(1))/XMU(I)
      A=1./DTAUB
      B=A**2
      SP2(1,I)=1.+2.*A
      SP3(1,I)=-2.*B
*
* LET P BE THE EVEN PART OF THE INTENSITY, THEN
*
*         P(2)= P(1) + D*P'(1) + .5*D2*P''(1)
* OR      P(2)= P(1) + D*(P(1)-I(1,-MU)) + .5*D2*(P(1)-S(1)) .
* WHERE   I(1,-MU) = S(1)*(1.-EXP(-T))
*
* THE DIFFERENCE AS COMPARED TO THE USUAL SECOND ORDER BOUNDARY 
* CONDITION IS THE ADDITIONAL TERM   I(1,-MU)=S(1)*(1.-EXP(-T)).
* THUS THE COEFFICIENT FOR S(1) IN THE FIRST EQUATION 
* SHOULD BE CHANGED AS FOLLOWS:
*         S(1)=S(1)*(1.+C*(1.-EXP(-T))
* WHERE   C=2./D
* *NORD* 751009
*
      C(I)=2.*A
      T(I)=TAU(1)*(X(1)+S(1))/XMU(I)
      IF(T(I).GT.70.) GOTO 700
      EX(I)=T(I)*(1.-.5*T(I)*(1.-.3333*T(I)))
      IF(T(I).GT.0.1) EX(I)=1.-EXP(-T(I))
      GOTO 701
  700 EX(I)=1.
  701 CONTINUE
*
* K=2,JTAU-1
*
      DO 100 K=2,JTAU1
        DTAUA=DTAUB
        DTAUB=.5*(X(K)+S(K)+X(K+1)+S(K+1))*(TAU(K+1)-TAU(K))/XMU(I)
        DTAUC=.5*(DTAUA+DTAUB)
        AD(K)=.166667*DTAUA/DTAUC
        BD(K)=.166667*DTAUB/DTAUC
        SP1(K,I)=-1./(DTAUA*DTAUC)+AD(K)
        SP2(K,I)=1.
        SP3(K,I)=-1./(DTAUB*DTAUC)+BD(K)
  100 CONTINUE
*
* K=JTAU
*
      SP2(JTAU,I)=1.
*
* END OF MU LOOP
*
  110 CONTINUE
*
* ELIMINATE SUBDIAGONAL, SAVE FACTORS IN SP1
*
      DO 121 I=1,MMU
        DO 120 K=1,JTAU2
          SP1(K,I)=-SP1(K+1,I)/(SP2(K,I)-SP3(K,I))
          SP2(K+1,I)=SP2(K+1,I)+SP1(K,I)*SP2(K,I)
          SP2(K,I)=SP2(K,I)-SP3(K,I)
  120   CONTINUE
        SP2(JTAU-1,I)=SP2(JTAU-1,I)-SP3(JTAU-1,I)
  121 CONTINUE
*
      RETURN
*
  400 CONTINUE
*
* ZEROSET
*
      DO 130 K=1,JTAU
        AN(K)=(3.*XK(K)-XJ(K))/8.*S(K)/(X(K)+S(K))
        XK(K)=0.
        XJ(K)=0.
  130 CONTINUE
*
* MU LOOP
*
      XH(1)=0.
      HSURF=0.
      DO 170 I=1,MMU
*
* INITIATE APPROXIMATIVE SOURCE FUNCTION
*
        P(1)=SOURCE(1)+AN(1)*(3.*XMU2(I)-1.)
*
* NOTE THE ANISOTROPIC SCATTERING CORRECTION
*
        S0=P(1)
        P(1)=P(1)*(1.+C(I)*EX(I))
        DO 140 K=2,JTAU1
          P(K)=(1.-AD(K)-BD(K))*(SOURCE(K)+AN(K)*(3.*XMU2(I)-1.))+
     &         AD(K)*(SOURCE(K-1)+AN(K-1)*(3.*XMU2(I)-1.))+
     &         BD(K)*(SOURCE(K+1)+AN(K+1)*(3.*XMU2(I)-1.))
  140   CONTINUE
        P(JTAU)=SOURCE(JTAU)
*
* ACCUMULATE RIGHT HAND SIDE
*
        DO 150 K=1,JTAU2
          P(K+1)=P(K+1)+SP1(K,I)*P(K)
  150   CONTINUE
*
* BACKSUBSTITUTE
*
        DO 160 K=1,JTAU1
          P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K,I)*P(JTAU-K+1))/SP2(JTAU-K,I)
          XK(JTAU-K)=XK(JTAU-K)+H(I)*P(JTAU-K)*XMU2(I)
          XJ(JTAU-K)=XJ(JTAU-K)+H(I)*P(JTAU-K)
  160   CONTINUE
*
* END OF MU LOOP
*
        XK(JTAU)=XK(JTAU)+H(I)*P(JTAU)*XMU2(I)
        R1=P(1)-S0*EX(I)
        XH(1)=XH(1)+H(I)*XMU(I)*R1
        P0=P(1)*(1.-EX(I))+.5*S0*EX(I)**2
        HSURF=HSURF+H(I)*XMU(I)*P0
        Y1(I)=2.*P0
*
* HSURF AND Y1(6) ARE THE FLUX AND INTENSITIES AT THE SURFACE
*
170   CONTINUE
*
      XJ(JTAU)=P(JTAU)
*
* 'XJ' IS THE NEW MEAN INTENSITY
*
      DO 180 K=1,JTAU
        ERROR(K)=(X(K)*BPLAN(K)+S(K)*XJ(K))/(X(K)+S(K))-SOURCE(K)
  180 CONTINUE
*
* FLUX AND SECOND MOMENT
*
      DO 190 K=2,JTAU
        XH(K)=2.*(XK(K)-XK(K-1))/(X(K)+S(K)+X(K-1)+S(K-1))/
     &        (TAU(K)-TAU(K-1))
  190 CONTINUE
*
      RETURN
      END
