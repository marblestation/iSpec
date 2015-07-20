C Date: Thu, 09 Nov 1995 16:16:23 +0100
C From: "Uffe G. Jorgensen, Niels Bohr Institute" <UFFEGJ@nbivax.nbi.dk>
C Subject: os-subroutine for h2o
C To: plez@nordita.dk
* modifiee le 14/2-1996 pour incorporation dans Bsyn (BPz)

c
      subroutine  h2opartf(mtemp,temp,u)
      integer mtemp
      real temp(mtemp),u(mtemp)
C
      real SUM(MTEMP),Q(MTEMP),QS(MTEMP),G(5000),O(3),OA(3),OB(3),
     &     AC(3),AA(3),AB(3)
C

      DATA IPART/1/
      DATA CKMS /2.99792E5/
C
C
C Compute the partition function for all relevant temperatures
C 1.Summarized partition function over computed vibrational levels
C   
      ktemp = mtemp
      OPEN (UNIT=81,FILE=
     & 'DATA/h2olevels.morbid',
cc     & '/b1/plez/BIGGRID/DATA/h2olevels.morbid',
     & STATUS='OLD')
      DO 224 IT = 1,KTEMP
      QS(IT) = 0.
224   Q(IT) = 0.
      IV1MAX = 0
      IV2MAX = 0
      IV3MAX = 0
      ISUMMAX = 0
      DO 210 I=1,10000
      READ(81,*,END=299) NRUN,IV1,IV2,IV3,ELEV,DIP
      IF (IV1 .GT. IV1MAX) IV1MAX = IV1
      IF (IV2 .GT. IV2MAX) IV2MAX = IV2
      IF (IV3 .GT. IV3MAX) IV3MAX = IV3
      ISUM = IV1 + IV2 + IV3
      IF (ISUM.GT. ISUMMAX) ISUMMAX = ISUM
      HCKT3 = 1.4388/3500.
      BOLTZ3=EXP(-ELEV*HCKT3)
      Q3500 = Q3500 + BOLTZ3
      DO 220 IT = 1,KTEMP
      HCKT = 1.4388/temp(IT)
      BOLTZ=EXP(-ELEV*HCKT)
220   Q(IT) = Q(IT) + BOLTZ
210   CONTINUE
299   CONTINUE
C     HARMONIC PARTITION FUNCTION:
      DO 226 IT = 1,KTEMP
      HCKT = 1.4388/temp(IT)
      CPARTV=(1.-EXP(-3657.0*HCKT))*(1.-EXP(-1594.7*HCKT))   
     & *(1.-EXP(-3755.7*HCKT))
      QHARM = 1./CPARTV
      PART=Q(IT)*temp(IT)**1.5
226   WRITE(6,112) temp(IT),QHARM,Q(IT) ,part
112   FORMAT(' T, HARMONIC, SUMMARIZED QVIB, partf: ',F8.0,
     &       2F9.4,x,f10.1)
C
      OPEN (UNIT=82,FILE=
     &  'DATA/h2omolk.dat',
cc     & '/b1/plez/BIGGRID/DATA/h2omolk.dat',
     & STATUS='OLD')
      READ(82,*) (O(I),OA(I),OB(I),I=1,3)
     & ,X11,X11A,X11B,X22,X22A,X22B,X33,X33A,X33B
     & ,X12,X12A,X12B,X13,X13A,X13B,X23,X23A,X23B
     & ,BEA,(AA(I),I=1,3)
     & ,BEB,(AB(I),I=1,3)
     & ,BEC,(AC(I),I=1,3),QE
      CLOSE(82)
      V1=0.5
      V2=0.5
      V3=0.5
      G0=O(1)*V1 + O(2)*V2 + O(3)*V3 + X11*V1**2 + X22*V2**2 + X33*V3**2
     &   + X12*V1*V2 + X13*V1*V3 + X23*V2*V3
C
      NV = 0
      DO 222 I1 = 0,IV1MAX
      DO 222 I2 = 0,IV2MAX
      DO 222 I3 = 0,IV3MAX
      ISUM = IV1 + IV2 + IV3
      IF (ISUM.GT. ISUMMAX) GO TO 222
      NV = NV + 1
      IF (NV.GT.5000) STOP ' *** DIMENSION PROBLEMS FOR G(V) *** '
      V1=I1+0.5
      V2=I2+0.5
      V3=I3+0.5
      G(NV)=O(1)*V1 + O(2)*V2 + O(3)*V3 + X11*V1**2 + X22*V2**2
     &    + X33*V3**2 + X12*V1*V2 + X13*V1*V3 + X23*V2*V3  - G0
      DO 228 IT = 1,KTEMP
      HCKT = 1.4388/temp(IT)
      BOLTZ=EXP(-G(NV)*HCKT)
228   QS(IT) = QS(IT) + BOLTZ
222   CONTINUE
C
      WRITE(6,*)' IV1-3MAX,ISUMMAX= ',IV1MAX,IV2MAX,IV3MAX,ISUMMAX
      WRITE(6,*) ' Qvib summed over ',nv,' anharmonic energy levels:'
      DO 229 IT = 1,KTEMP
      PART=Qs(IT)*temp(IT)**1.5
      u(it)=part
229   WRITE(6,113) temp(IT),QS(IT) ,part
113   FORMAT(' T, QVIB  partf: ',F8.0,F11.4,x,f10.1)
      close(81)
      return
C
       END
C
C
