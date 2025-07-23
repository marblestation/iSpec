      SUBROUTINE TRANSCplatt(ISTEP)
*
*-----------------------------------------------------------------------
*
* SCATTR SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING
* IN THE EDDINGTON APPROXIMATION, I.E., USING ONLY ONE MU POINT.
* 'ERROR' IS THE INHOMOGENEOUS TERM OF THE EQUATION, AND 'P' CONTAINS
* THE ESTIMATED MEAN INTENSITY ON EXIT. TRANSC CALCULATES THE MATRIX
* ELEMENTS FOR SCATTR.
* 79.06.21 *NORD*
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
* 
      include 'spectrum.inc'
*
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),hflux(NDP),XK(NDP),
     &        fillup(4*ndp+3)
      COMMON/TAUC/ TAU(ndp),DTAULN(ndp),JTAU
      COMMON/SPACE3/ SOURCE(ndp),ERROR(ndp),sp1(ndp),sp2(ndp),sp3(ndp),
     &               P(ndp),
     &               tomatch(7*nrays*ndp+4*ndp+3*nrays+1)
*
      DATA XMU,XMU2/0.5773503,0.3333333/
*
      IF(ISTEP.GT.0) GOTO 400
*
* K=1
*
      DTAUB=.5*(X(1)+S(1)+X(2)+S(2))*(TAU(2)-TAU(1))/XMU
      A=1./DTAUB
      B=A**2
      SP2(1)=2.*A+X(1)/(S(1)+X(1))
      SP3(1)=-2.*B
      C=2.*A
      T=TAU(1)*(X(1)+S(1))/XMU
      IF(T.GT.70) GOTO 700
      EX=T*(1.-.5*T*(1.-.33333*T))
      IF(T.GT.0.1) EX=1.-EXP(-T)
      GOTO 701
  700 EX=1.
  701 CONTINUE
      SP2(1)=SP2(1)-C*S(1)/(X(1)+S(1))*EX
*
* K=2,JTAU-1
*
      JTAU1=JTAU-1
      DO 100 K=2,JTAU1
        DTAUA=DTAUB
        DTAUB=.5*(X(K)+S(K)+X(K+1)+S(K+1))*(TAU(K+1)-TAU(K))/XMU
        DTAUC=.5*(DTAUA+DTAUB)
        A=1./(DTAUA*DTAUC)
        B=1./(DTAUB*DTAUC)
        SP1(K)=-A
        SP2(K)=X(K)/(S(K)+X(K))
        SP3(K)=-B
  100 CONTINUE
*
* K=JTAU
*
      SP2(JTAU)=X(JTAU)/(X(JTAU)+S(JTAU))
*
* ELIMINATE SUBDIAGONAL
*
      JTAU2=JTAU-2
      DO 110 K=1,JTAU2
        SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
        SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
        SP2(K)=SP2(K)-SP3(K)
  110 CONTINUE
      SP2(JTAU-1)=SP2(JTAU-1)-SP3(JTAU-1)
      RETURN
*
  400 CONTINUE
*
* INITIATE INHOMOGENOUS TERMS
*
      DO 120 K=1,JTAU
        P(K)=ERROR(K)
  120 CONTINUE
      DSDT=0.
*
* PRELIM
*
      P(1)=P(1)*(1.+C*EX)-DSDT*C*(EX-T*(1.-EX))
*
* ACCUMULATE INHOMOGENOUS TERMS
*
      DO 130 K=1,JTAU2
        P(K+1)=P(K+1)+SP1(K)*P(K)
  130 CONTINUE
*
* BACKSUBSTITUTE
*
      P(JTAU)=P(JTAU)/SP2(JTAU)
      DO 140 K=1,JTAU1
        P(JTAU-K)=(P(JTAU-K)-SP3(JTAU-K)*P(JTAU-K+1))/SP2(JTAU-K)
  140 CONTINUE
*
      RETURN
      END
