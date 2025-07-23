C
      SUBROUTINE TRANSC
C
C SCATTR SOLVES THE TRANSFER EQUATION INCLUDING CONTINUUM SCATTERING
C IN THE EDDINGTON APPROXIMATION, I.E., USING ONLY ONE RAY.
C 'ERROR' IS THE INHOMOGENEOUS TERM OF THE EQUATION, AND 'P' GIVES THE
C ESTIMATED MEAN INTENSITY CORRECTION (FJ*P). TRANSC CALCULATES THE MATRIX
C ELEMENTS FOR SCATTR.
C 79.06.21 *NORD*
C
C
      implicit none
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
C
      INCLUDE 'spectrum.inc'
C
      real(sp)  :: x(ndp),s(ndp),bplan(ndp),xj(ndp),xh(ndp),xk(ndp)
      real(sp)  :: fj(ndp),source(ndp),taus(ndp),dtaus(ndp)
      integer   :: jtau0,jtau1,iscat 
      COMMON /CTRAN/ X,S,BPLAN,XJ,XH,XK,FJ,SOURCE,TAUS,DTAUS,
     &              JTAU0,JTAU1,ISCAT

      real(sp)  :: tau(ndp),dtauln(ndp)
      integer   :: jtau
      COMMON /TAUC/ TAU,DTAULN,JTAU

! this SPACE2 common block is not the same as in traneq.f or trrays.f
! SP1, SP2 and SP3 variables are not located at the same place.

!      real(sp)  :: error(ndp),fact(ndp),dso(ndp),p(ndp),dum(ndp,3)
!      real(sp)  :: sp1(ndp,nrays),sp2(ndp,nrays),sp3(ndp,nrays)
      real(sp)  :: error(ndp),fact(ndp),dso(ndp),p(ndp),dum(ndp,nrays,3)
      real(sp)  :: sp1(ndp),sp2(ndp),sp3(ndp)
      real(sp)  :: ad(ndp,nrays),bd(ndp,nrays),ex(nrays),pimpac(nrays)
      integer   :: nimpac,kimpac(nrays),mmu(ndp)
      real(sp)  :: taut(ndp),dtaut(ndp),pfeau(nrays,ndp),xmu(nrays,ndp)
!      COMMON /SPACE2/ ERROR,FACT,DSO,P,DUM,SP1,SP2,SP3,AD,BD,EX,
      COMMON /SPACE2/ ERROR,FACT,DSO,P,SP1,SP2,SP3,DUM,AD,BD,EX,
     &  NIMPAC,KIMPAC,PIMPAC,MMU,TAUT,DTAUT,PFEAU,XMU

      integer   :: ncore
      real(sp)  :: diflog,radius,rr(ndp)
      COMMON /CSPHER/ NCORE,DIFLOG,RADIUS,RR

      integer   :: k,ntau,ntau1
      real(sp)  :: a,b,dtauc,dz,dzdr,t,z,zold
C
C MU LOOP
      NTAU=KIMPAC(ISCAT)
      NTAU1=NTAU-1
C
C CALCULATE TAUS ALONG THE RAY, SPIRAL AT EDDINGTON ANGLE AT DEPTH.
! BPz 2024-12-05: gfortran with -O2 may result in <0 values of the argument of sqrt, when rr=pimpac
      Z=SQRT(max(0.0,RR(JTAU0)**2-PIMPAC(ISCAT)**2))
      ZOLD=Z
      DO 101 K=2,JTAU
      IF (JTAU-K+2.GT.JTAU0) GO TO 103
      Z=SQRT(max(0.0,RR(JTAU-K+1)**2-PIMPAC(ISCAT)**2))
      DZ=Z-ZOLD
      DZDR=DZ/(RR(JTAU-K+1)-RR(JTAU-K+2))
      GO TO 104
103   DZDR=1.732
104   DTAUS(JTAU-K+2)=DZDR*0.5*(X(JTAU-K+1)+S(JTAU-K+1)
     & +X(JTAU-K+2)+S(JTAU-K+2))*(TAU(JTAU-K+2)-TAU(JTAU-K+1))
101   ZOLD=Z
C
      TAUS(1)=DZDR*(X(1)+S(1))*TAU(1)
C
C K=1
      A=1./DTAUS(2)
      B=A**2
      SP2(1)=1.-FJ(1)*S(1)/(X(1)+S(1))+2.*A
      SP3(1)=-2.*B
      T=TAUS(1)
      EX(ISCAT)=TAUS(1)*(1.-TAUS(1)/2.*(1.-TAUS(1)/3.*(1.-taus(1)/4.)))
      IF (TAUS(1).GT.0.01) EX(ISCAT)=1.-EXP(-TAUS(1))
      SP2(1)=SP2(1)-2.*A*EX(ISCAT)*FJ(1)*S(1)/(X(1)+S(1))
      SP2(1)=SP2(1)/(1.+2.*A*EX(ISCAT))
      SP3(1)=SP3(1)/(1.+2.*A*EX(ISCAT))
C
C K=2,NTAU-1
      DO 100 K=2,NTAU1
      DTAUC=0.5*(DTAUS(K)+DTAUS(K+1))
      SP1(K)=-1./(DTAUS(K)*DTAUC)
      SP2(K)=1.-FJ(K)*S(K)/(X(K)+S(K))
100   SP3(K)=-1./(DTAUS(K+1)*DTAUC)
C
C K=NTAU
      SP1(NTAU)=0.0
      SP2(NTAU)=X(NTAU)/(X(NTAU)+S(NTAU))
      SP3(NTAU)=0.0
C
C ELIMINATE SUBDIAGONAL
      DO K=1,NTAU1
        SP1(K)=-SP1(K+1)/(SP2(K)-SP3(K))
        SP2(K+1)=SP2(K+1)+SP1(K)*SP2(K)
        SP2(K)=SP2(K)-SP3(K)
      ENDDO
      CONTINUE
      RETURN
C--------------------------------------------------------------------
C
      ENTRY SCATTR
C
C MU LOOP
      NTAU=KIMPAC(ISCAT)
      NTAU1=NTAU-1
C
C ACCUMULATE RIGHT HAND SIDE
      P(1)=ERROR(1)
      DO 150 K=1,NTAU1
150   P(K+1)=ERROR(K+1)+SP1(K)*P(K)
C
C BACKSUBSTITUTE
      P(NTAU)=P(NTAU)/SP2(NTAU)
      DO 160 K=1,NTAU1
160   P(NTAU-K)=(P(NTAU-K)-SP3(NTAU-K)*P(NTAU-K+1))/SP2(NTAU-K)
C
      RETURN
      END
