C
      SUBROUTINE TRRAYS
C
C
      INCLUDE 'spectrum.inc'
C
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),XH(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),JTAU
      COMMON /SPACE2/ERROR(NDP),FACT(NDP),DSO(NDP),
     &  P(NDP),DUM(NDP,3),
     &  SP1(NDP,NRAYS),SP2(NDP,NRAYS),SP3(NDP,NRAYS),AD(NDP,NRAYS),
     &  BD(NDP,NRAYS),EX(NRAYS),
     &  NIMPAC,KIMPAC(NRAYS),PIMPAC(NRAYS),MMU(NDP),
     &  TAUT(NDP),DTAUT(NDP),
     &  PFEAU(NRAYS,NDP),XMU(NRAYS,NDP)
      COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP)
      COMMON /CSTYR/MIHAL  /CTAUM/TAUM
*
*  Distribute rays, based on two depth indices, JTAU0 and JTAU1.
*
*  JTAU0 is the largest depth index for which TAU < SQRT((X+S)/X),
*  and represents the "surface", where radiation is released.
*  Above JTAU0, linearized perturbations in the radiation field are
*  represented by a single ray, which is chosen in TRANFR as the ray
*  which is closest to the Eddington angle 1./sqrt(3.) at JTAU0.
*
*  JTAU1 is the largest depth index for which TAU < TAUM*SQRT((X+S)/X),
*  and represents the border to the "core", where diffusion is a good
*  approximation.  The detailed formal solution of the radiative transfer
*  on a set of parallel rays is only performed above JTAU1.  TAUM should
*  be at least 50 to 100.
*
      TAUT(1)=TAU(1)*(X(1)+S(1))
      JTAU0=1
      JTAU1=1
      ICASE=1
      DO 100 K=2,JTAU
      TAUT(K)=TAUT(K-1)+0.5*(X(K)+S(K)+X(K-1)+S(K-1))*(TAU(K)-TAU(K-1))
      GO TO (101,102,103,104,105),ICASE
101   IF (TAUT(K).LT.1.0) GO TO 105
      ICASE=2
102   IF (TAUT(K).LT.SQRT((X(K)+S(K))/X(K))) GO TO 105
      ICASE=3
103   IF (TAUT(K).LT.TAUM) GO TO 105
      ICASE=4
104   IF (TAUT(K).LT.TAUM*SQRT((X(K)+S(K))/X(K))) GO TO 105
      ICASE=5
105   IF (ABS(XH(K))/BPLAN(K).GT.0.10) ICASE=MIN0(ICASE,4)
      GO TO (106,106,107,107,100),ICASE
106   JTAU0=K
107   JTAU1=K
100   CONTINUE
      JTAU1=MAX0(JTAU1,3)
*
*  Distribute rays in the core.  This is done in such a way that the
*  mu values at depth JTAU0 are equidistant.  The rationale for this
*  is to get a good representation of the the "core" part of the
*  radiation as a function of mu, in the optically thin regions.
*
      RRK=RR(JTAU0)
      NCORE1=NCORE-1
      DR=(RR(JTAU0)-SQRT(RR(JTAU0)**2-RR(JTAU1)**2))/NCORE1
      DO 110 I=1,NCORE1
      PIMPAC(I)=SQRT(RR(JTAU0)**2-RRK**2)
      RRK=RRK-DR
110   KIMPAC(I)=JTAU1
C
C RAYS IN ATMOSPHERE
      I=NCORE
      KI=JTAU1
120   KIMPAC(I)=KI
      PIMPAC(I)=RR(KI)
      KIP=KI
121   KI=KI-1
* changed from tau to taut on  13-jan-92.
* the average number of rays is decreased from 50 to 37.
      IF (TAUt(KIP)/TAUt(KI).LT.DIFLOG.AND.KI.GT.1) GO TO 121
      I=I+1
      IF (KI.GE.3) GO TO 120
      KIMPAC(I)=0
      PIMPAC(I)=RR(1)
      NIMPAC=I-1
      IF (NIMPAC.LT.NRAYS) GO TO 131
      PRINT 122,NIMPAC,NCORE,DIFLOG
122   FORMAT('0**** NMBR OF RAYS TOO LARGE =',I3,'  NCORE,DIFLOG =',
     & I3,F6.3)
      STOP
C
C FIND THE NUMBER OF MU-PNTS FOR EACH K, PLUS ONE EXTRA FOR MU=0.0
131   II=NIMPAC+1
      DO 130 K=1,JTAU1
      MMU(K)=II
      XMU(II,K)=0.0
      IF (K+1.EQ.KIMPAC(II-1)) II=II-1
130   CONTINUE
      RETURN
      END
