      SUBROUTINE TRANEQplatt(idebug)
*
*-----------------------------------------------------------------------
*
* TRANEQ SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING.
* FEATURES:
*
* 1. CANNONS PERTURBATION TECHNIQUE IS USED ON THE ANGULAR QUADRATURE.
*    THE BASIC IDEA IN THIS TECHNIQUE IS TO REPLACE THE INVERSION OF
*    A COMPLICATED (MMU ORDER) OPERATOR WITH THE INVERSION OF A SIMPLE
*    OPERATOR (ONE POINT=EDDINGTON APPROXIMATION), PLUS ITERATION ON
*    THE ERROR.
* 2. AITKEN EXTRAPOLATION ACCELLERATES THE CONVERGENCE.
* 3. A TRICK DUE TO ROBERT STEIN (PRIV. COMM., 1979) IS USED TO
*    ELIMINATE THE NEED FOR DOUBLE PRECISION STORAGE OF THE MATRIX
*    ELEMENTS. THE IDEA IS TO STORE THE (SMALL) SUM OF THE THREE
*    MATRIX ELEMENTS ON A ROW, INSTEAD OF THE (LARGE) DIAGONAL ELE-
*    MENT.
* 4. THE SOLUTION IS A CUBIC SPLINE, RATHER THAN A PIECE-WISE
*    QUADRATIC FUNCTION. THIS IS ACCOMPLISHED WITH THE CORRECTION
*    TERMS AD AND BD IN SUBROUTINE TRANFR.
* 5. THE SCATTERING IS TREATED AS DIPOLE SCATTERING INSTEAD OF THE 
*    NORMALLY USED ISOTROPIC APPROXIMATION. 
*    THIS CAN BE DONE VERY SIMPLY IN THE ITERATING CANNON SCHEME.
* 6. A BOUNDARY CONDITION WHICH INCLUDES AN ESTIMATED INFALLING
*    RADIATION MAKES THE SOLUTION GOOD ALSO FOR VALUES OF X+S
*    LARGE COMPARED WITH 1./TAU(1). A LOGARITHMIC TAU-SCALE
*    SHOULD BE USED.
*
* THIS VERSION OF TRANEQ IS COMPATIBLE WITH PREVIOUS TRANEQ'S.
* 79.06.21 *NORD*
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      include 'spectrum.inc'
      logical debug
      COMMON/TAUC/ TAU(ndp),DTAULN(ndp),JTAU
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),xh(NDP),XK(NDP),
     & fillup(4*ndp+3)
      COMMON/SPACE2/ SOURCE(ndp),ERROR(ndp),DUM(3*ndp),P(ndp),
     &               SP1(ndp,nrays),SP2(ndp,nrays),SP3(ndp,nrays),
     &               AN(ndp),AD(ndp),BD(ndp),FACT(ndp),DSO(ndp),
     &               tomatch(4*nrays*ndp-ndp+3*nrays+1)
*
      PARAMETER (ITMAX=15)
      DIMENSION A(ITMAX)

      if (idebug.eq.0) then 
        debug=.false.
      else
        debug=.true.
      endif
*
* INITIATE
*
      DO 100 K=1,JTAU
        FACT(K)=1.
        DSO(K)=0.
        XJ(K)=0.
        XK(K)=0.
        ERROR(K)=BPLAN(K)*X(K)/(X(K)+S(K))
        SOURCE(K)=0.
  100 CONTINUE
*
* CALCULATE THE MATRIX ELEMENTS
*
      CALL TRANFRplatt(0)
      CALL TRANSCplatt(0)
*
      DO 110 IT=1,ITMAX
        A(IT)=0.
  110 CONTINUE
*
* ITERATION LOOP
*
      DO 140 IT=1,ITMAX
*
        ITM=IT
*
* SOLVE THE CONTINUUM SCATTERING PROBLEM IN THE EDDINGTON APPROXIMATION
*
        CALL TRANSCplatt(1)
        DO 120 K=1,JTAU
          XJ(K)=XJ(K)+P(K)
          XK(K)=XK(K)+.333333*P(K)
*
* AITKEN EXTRAPOLATION USED FOR CONVERGENCE ACCELLERATION
*
          DS=ERROR(K)+P(K)*S(K)/(X(K)+S(K))
          IF(DSO(K).NE.0.) 
     &       FACT(K)=AMIN1(1.25,AMAX1(.8,FACT(K)-DS/DSO(K)))
          DS=DS/FACT(K)
          IF(IT.GE.2) DSO(K)=DS
          SOURCE(K)=SOURCE(K)+DS
  120   CONTINUE
*
* SOLVE THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
*
        CALL TRANFRplatt(1)
*
* CHECK ERROR IN SOURCE FUNCTION
*
        DO 130 K=1,JTAU
          A(IT)=AMAX1(A(IT),ABS(ERROR(K)/SOURCE(K)))
  130   CONTINUE
        IF(A(IT).LT.0.0001) GOTO 141
*
  140 CONTINUE
*
* END OF ITERATION LOOP
*
      PRINT 50,(A(IT),IT=1,ITM)
   50 FORMAT(' Traneq-MAXFEL =',8F10.7)
*
  141 CONTINUE

      if (debug) then
        print*,'traneqplatt: '
        print*,' k   tau(ref)   tau(lambda)    J    B     S    x   s'
        do k=1,1
          taut=TAU(1)*(X(1)+S(1))
          print*,k,tau(k),taut,xj(k),bplan(k),source(k),x(k),s(k)
        enddo
        do k=2,jtau
          taut=taut+.5*(X(K-1)+S(K-1)+X(K)+S(K))*(TAU(K)-TAU(K-1))
          print*,k,tau(k),taut,xj(k),bplan(k),source(k),x(k),s(k)
        enddo
      endif
*
      RETURN
      END
