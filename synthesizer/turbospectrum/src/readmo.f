      SUBROUTINE READMO(XL,X,S)
*
*-----------------------------------------------------------------------
*
*  READS MODEL ATMOSPHERE AND CALCULATES ABSORPTION AND SCATTERING
*  COEFFICIENTS AT RELEVANT OPTICAL DEPTHS AND WAVELENTH XL FROM
*  MODEL ATM FILE
*  TAU = TAU(5000)
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      INCLUDE 'spectrum.inc'
*
* for the dimensions, refer to babsma
      doubleprecision xxll,xlp(20*numbset)
      DIMENSION XS(20*numbset),SS(20*numbset),
     &   X(ndp),S(ndp),TAULN(NDP)
      character*50 MCODE
      REAL MUM
      common/babcont/ xlp,nlq
*
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),NTAU
      logical hydrovelo
      real velocity
      common/velo/velocity(ndp),hydrovelo
      COMMON/TAUC/ TAU(NDP),DTAULN(NDP),JTAU
      COMMON/ROSSC/ ROSS(NDP), cross(ndp)
      COMMON/CWAVES/ XLS
      COMMON/MODID/ MCODE
      COMMON/CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP)
      COMMON/RHOC/RHO(NDP)
* extension for large number of wavelengths and lines (monster II)
      doubleprecision xlambda
      common/large/ xlambda(lpoint),maxlam,ABSO(NDP,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint)
*
ccc      READ(12,100) MCODE,NTAU,XLS
      READ(12,*) MCODE,NTAU,XLS
      READ(12,103) NLQ
      if (nlq.gt.20*numbset) then
         print*,'readmo: nlq= ',nlq,' larger than numbset*20= ',
     &    numbset*20
         stop
      endif
      READ(12,104) (XLP(I), I=1,NLQ)
      DO 11 K=1,NTAU
        if (hydrovelo) then
          READ(12,*,err=99) RR(K),TAU(K),T(K),PE(K),PG(K),
     &                          RO(K),XI(K),
     &                          ROSS(K),velocity(k)
        else
          READ(12,*) RR(K),TAU(K),T(K),PE(K),PG(K),RO(K),XI(K),
     &                 ROSS(K)
          velocity(k)=0.0
        endif
* ross is standard opacity (not rosseland)
        READ(12,102) (XS(J),SS(J),J=1,NLQ)
*
* INTERPOLATE TO WAVELENGTH XL
*
         xxll=xl
         CALL LINT(NLQ,XLP,XS,xxll,X(K))
         CALL LINT(NLQ,XLP,SS,xxll,S(K))
* interpolate to all wavelengths
        do j=1,maxlam
         xxll=xlambda(j)
         CALL LINT(NLQ,XLP,XS,xxll,abso(K,j))
         CALL LINT(NLQ,XLP,SS,xxll,absos(K,j))
        enddo
*        
        XI(K)=XI(K)*1.E5
        TAULN(K)=ALOG(TAU(K))
        if (tau(k)-1..lt.0.1) mmm=k
        MUM(K)=(1.38*RO(K)*T(K))/(1.67E-08*PG(K))
        RHO(K)=RO(K)
        IF(K.EQ.1) GOTO 11
        DTAULN(K)=TAULN(K)-TAULN(K-1)
  11  CONTINUE
      RADIUS=RR(mmm)
      REWIND 12
      JTAU=NTAU
*
 100  FORMAT(1X,A50,I5,F10.2)
 101  FORMAT(1X,8E11.4)
 111  FORMAT(1X,9E11.4)
 102  FORMAT(1X,6E11.4)
 103  FORMAT(1X,I5)
 104  FORMAT(1X,10F11.3)
*
      RETURN

99    print*,'ERROR in reading model. You seem to want to compute'
      print*,'      the optical depth in a model with velocity shifts.'
      print*,'      Please verify that your input model (which is an'
      print*,'      output from babsma) has been prepared with the'
      print*,'      HYDRODYN_DEPTH switch on for babsma!'
      stop

      END
