C
      SUBROUTINE TRRAYS
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

      real(sp)  :: error(ndp),fact(ndp),dso(ndp),p(ndp),dum(ndp,3)
      real(sp)  :: sp1(ndp,nrays),sp2(ndp,nrays),sp3(ndp,nrays)
      real(sp)  :: ad(ndp,nrays),bd(ndp,nrays),ex(nrays),pimpac(nrays)
      integer   :: nimpac,kimpac(nrays),mmu(ndp)
      real(sp)  :: taut(ndp),dtaut(ndp),pfeau(nrays,ndp),xmu(nrays,ndp)
      COMMON /SPACE2/ ERROR,FACT,DSO,P,DUM,SP1,SP2,SP3,AD,BD,EX,
     &  NIMPAC,KIMPAC,PIMPAC,MMU,TAUT,DTAUT,PFEAU,XMU

      integer   :: ncore
      real(sp)  :: diflog,radius,rr(ndp)
      COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR

      integer   :: mihal
      COMMON /CSTYR/MIHAL

      integer   :: taum
      common  /CTAUM/TAUM

      real(sp)  :: dr,rrk
      integer   :: i,ii,k,ki,kip,ncore1
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

      do k=2,jtau
        TAUT(K)=TAUT(K-1)+0.5*(X(K)+S(K)+X(K-1)+S(K-1))*
     &          (TAU(K)-TAU(K-1))
        if (TAUT(K).LT.1.0) then
          jtau0=k
          jtau1=k
        else if (TAUT(K).LT.SQRT((X(K)+S(K))/X(K))) then
          jtau0=k
          jtau1=k
        else if (TAUT(K).LT.TAUM) then
          jtau1=k
        else if (TAUT(K).LT.TAUM*SQRT((X(K)+S(K))/X(K))) then
          jtau1=k
        else if (ABS(XH(K))/BPLAN(K).GT.0.10) then
          jtau1=k
        endif
      enddo

      JTAU1=MAX0(JTAU1,3)
*
*  Distribute rays in the core.  This is done in such a way that the
*  mu values at depth JTAU0 are equidistant.  The rationale for this
*  is to get a good representation of the the "core" part of the
*  radiation as a function of mu, in the optically thin regions.
*
!!      print*,jtau0,jtau1

      RRK=RR(JTAU0)
      NCORE1=NCORE-1
      DR=(RR(JTAU0)-SQRT(RR(JTAU0)**2-RR(JTAU1)**2))/NCORE1

!!      print*,ncore, ncore1, rr(jtau1),dr

      do I=1,NCORE1
! BPz 2024-12-05: gfortran with -O2 may result in <0 values for the argument of sqrt, when rr=rrk
        PIMPAC(I)=SQRT(max(0.0,RR(JTAU0)**2-RRK**2))
!!        print*,'trrays',i,jtau0,rr(jtau0),rrk,pimpac(i)
        RRK=RRK-DR
        KIMPAC(I)=JTAU1
!!        print*,i,pimpac(i),rrk,kimpac(i)
      enddo
C
C RAYS IN ATMOSPHERE
      I=NCORE
      KI=JTAU1

!!      print*,i,ki
      do while (KI.ge.3)
        KIMPAC(I)=KI
        PIMPAC(I)=RR(KI)
        KIP=KI
* changed from tau to taut on  13-jan-92.
* the average number of rays is decreased from 50 to 37.
        do while (TAUt(KIP)/TAUt(KI).LT.DIFLOG.AND.KI.GT.1)
          KI=KI-1
        enddo
        I=I+1
!!        print*,ki,i
      enddo

      KIMPAC(I)=0
      PIMPAC(I)=RR(1)
      NIMPAC=I-1
!!      print*,'out',nimpac

      if (NIMPAC.ge.NRAYS) then
        PRINT 122,NIMPAC,NCORE,DIFLOG
122     FORMAT('0**** NMBR OF RAYS TOO LARGE =',I3,'  NCORE,DIFLOG =',
     &         I3,F6.3)
        STOP
      endif
C
C FIND THE NUMBER OF MU-PNTS FOR EACH K, PLUS ONE EXTRA FOR MU=0.0
      II=NIMPAC+1

      do K=1,JTAU1
        MMU(K)=II
        XMU(II,K)=0.0
        IF (K+1.EQ.KIMPAC(II-1)) II=II-1
      enddo

      RETURN
      END
