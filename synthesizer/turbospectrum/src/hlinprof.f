C
C   Table of contents:
C
C     SUBROUTINE HLINPROF(WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,
C    *                   H1FRC,HE1FRC,DOP,PROF)
C     SUBROUTINE STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
C     SUBROUTINE HGRID(HFILE,WAVE,WAVE0,TEMP,NHYD,NLOW,NUP,HLIN)
C     SUBROUTINE HSELF_PDWIDTH(LINE,TEMP,H1FRC,WIDTH)
C     SUBROUTINE HE_PDWIDTH(LINE,TEMP,HE1FRC,HEWID)
C     SUBROUTINE RAD_WIDTH(LINE,WIDTH)
C     SUBROUTINE ALI_GRIEM(LINE,WAVE0,H1FRC,WIDTH)
C     SUBROUTINE VCS(PR,XNE,T,DEL,II,N,M)                                      
C
C
C  Hydrogen Line Routines
C  All routines by Paul Barklem unless otherwise stated
C
      SUBROUTINE HLINPROF(WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,
     *                   H1FRC,HE1FRC,DOP,PROF)
C
C  Returns H line Profile at specific point
C
C  Inputs: 
C          WAVE   - wavelength in AA
C          WAVE0  - line centre in AA
C          TEMP   - temperature in K
C          XXNE   - electron number density in cm^-3
C          NLOW   - n quantum number of lower level
C          NUP    - n quantum number of upper level
C          H1FRC  - neutral hydrogen number density in cm^-3
C          HE1FRC - neutral helium number density in cm^-3
C          DOP    - the doppler width delta lambda / lambda
C                   (only used if HLINOP is called instead!)
C  Outputs:
C          PROF   - area normalised (in AA) profile for WAVE
C  Common Inputs: 
C     (CURRENTLY SET HARD!!!! see below)
C          ISTARK - switch for different recipes
C                  0=none
C                  1=Stehle
C                  2=Vidal, Cooper & Smith
C          ISELF  - switch for self broadening
C                  0=none
C                  1=Barklem, Piskunov and O'Mara (grid)
C                  2=Barklem, Piskunov and O'Mara (p-d)
C                  3=Ali & Griem
C          ICONV  - switch for convolution method
C                  1=adding
C                  2=trapezoidal
C          IHE    - switch for broadening by He
C                  0=none
C                  1=estimated from Barklem et al (1997)
C          NB: ICONV=2 is incompatible with either 0 option 
C  
C  The trapezoidal convolution is only for testing.  Use at own risk!
C  It's left here in as it might be useful to someone.
C      
      IMPLICIT NONE

      INTEGER NKMAX
      PARAMETER (NKMAX = 10001)
      INTEGER NLOW,NUP,LINE
      INTEGER ISTARK,ISELF,ICONV
      INTEGER IHE
      INTEGER NLOWW,NHIGH,NKRN,IKRN,J,NMID,I1,I2
      REAL*4 TEMP,XXNE,PROF,H1FRC,HE1FRC
      REAL*4 HLIN,STARK
      REAL*4 HEWID,RAD,DOP
      REAL*4 HLINOP
      REAL*8 TD,XNED,H1FRCD
      REAL*8 WAVE,WAVE0,DEL,HHLIN
      REAL*8 KRNSTEP,WLKRN(NKMAX),KERN(NKMAX),KMAX,WV2
      REAL*8 DEL1,DDEL,XNEXP,F1,HLIND
      REAL*8 STARKD
      REAL*8 LOR,HWID,RES
      REAL*8 C,PI
      CHARACTER*80 HFILE,HVCSFILE,HSELFFILE
      LOGICAL FIRST
C
      SAVE FIRST 
C
      DATA FIRST/.TRUE./
      PARAMETER (C = 2.997925E+18, PI = 3.14159265)

!      INCLUDE 'DATA.FILES'
!      parameter (HFILE='DATA/stehle.dat.DEC')              ! little_endian file
!      parameter (HSELFFILE = 'DATA/bpo_self.grid.DEC')     ! little_endian file
      parameter (HFILE = 'DATA/stehle.dat.SUN')              ! big_endian file
      parameter (HSELFFILE = 'DATA/bpo_self.grid.SUN')     ! big_endian file
      parameter (HVCSFILE = 'DATA/vcsbalmer.dat')          ! ascii file
C
C  Set the switches here 
C
      COMMON /HSWITCH/ISTARK,ISELF,ICONV,IHE
      ISTARK = 1
      ISELF = 2
      ICONV = 1
      IHE = 1
C
C We check we have a line which we can handle
C  if not then pass to HLINOP
C
      IF ((NLOW .NE. 2) .OR. (NUP .GT. 5)) THEN
        IF (FIRST) THEN
          FIRST = .FALSE.
          PRINT*,' Using Hlinop code!'
        ENDIF
        PROF = HLINOP(WAVE,NLOW,NUP,WAVE0,TEMP,XXNE,H1FRC,
     ;                                          HE1FRC,DOP)
        PROF = PROF * C/WAVE/WAVE
      RETURN
      END IF 
C
C For historical reasons copy some variables to real*8
C
      XNED = XXNE
      TD = TEMP
      H1FRCD = H1FRC
C      
C Precomputed things here
C
      DEL = DABS(WAVE0 - WAVE)
      LINE = NUP - NLOW
C
C Adding for convolution of Stark and self
C
      IF (ICONV .EQ. 1) THEN
        HLIN = 0. 
        HHLIN = 0.
        IF (ISTARK .EQ. 1) THEN
          CALL STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
C
C If outside table range use VCS
C
          IF (HLIN.LT.0.) THEN
            CALL VCS(HLIND,XNED,TD,DEL,1,NLOW,NUP,HVCSFILE)
            HLIN = HLIND
          ENDIF
        END IF
        IF (ISTARK .EQ. 2) THEN
          CALL VCS(HLIND,XNED,TD,DEL,1,NLOW,NUP,HVCSFILE)
          HLIN = HLIND
        END IF
        IF ((DEL .GT. .2d0) .AND. (H1FRC .GT. 0.)) THEN
          IF (ISELF .EQ. 1) THEN
             CALL HGRID(HSELFFILE,
     ;            WAVE,WAVE0,TD,H1FRCD,NLOW,NUP,HHLIN)
C
C If outside table range use p-d approx
C
            IF (HHLIN.LT.0.) THEN
              CALL HSELF_PDWIDTH(LINE,TEMP,H1FRC,HWID)
              HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
              HHLIN = HWID / (HWID * HWID + DEL * DEL) / PI
            ENDIF
          ENDIF
          IF (ISELF .EQ. 2) THEN
            CALL HSELF_PDWIDTH(LINE,TEMP,H1FRC,HWID)
            HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
            HHLIN = HWID / (HWID * HWID + DEL * DEL) / PI
          END IF
          IF (ISELF .EQ. 3) THEN
            CALL ALI_GRIEM(LINE,WAVE0,H1FRC,RES)
            HHLIN = RES / (RES * RES + DEL * DEL) / PI
          END IF 
          HLIN = HLIN + HHLIN
        END IF
      END IF
C
C Convolution by simple trapezoidal method
C
C   Notes:  Should have Stark profile > Lorentzian
C   If ever this is broken you should reverse the two or 
C   increase NKRN points significantly.
C
      IF (ICONV.EQ.2) THEN
        IF (H1FRC .EQ. 0.) THEN
          PRINT*,'STOPPING : ZERO H1FRC, CONVOLVING WILL FAIL'
          STOP
        END IF
C
C Use a logarithmic spread of kernel points 
C
        KMAX = 100.     !we go out KMAX AA into wings
        NKRN = 1001     !use an odd number, testing may be 
C                       !needed to see what this number should 
C                       !for a given case
        NMID = (NKRN + 1) / 2
        NLOWW = NMID - 1
        NHIGH = NMID + 1
        XNEXP = DLOG(1.d-4)   !this is the low limit of detuning 
C                             !(the log of it in AA)
C                             !should be smaller than smallest half width
C                             !which are typically ~0.01AA for solar
        KRNSTEP = (LOG(KMAX) - XNEXP) / DBLE(NLOWW - 1)
        DO J = NKRN, NHIGH, -1
          WLKRN(J) = EXP(LOG(KMAX) - KRNSTEP * DBLE(NKRN - J))
        END DO
        WLKRN(NMID) = 0.d0
        DO J = 1, NLOWW, 1
          WLKRN(J) = - EXP(LOG(KMAX) - KRNSTEP * DBLE(J - 1))
        END DO
C
        DO IKRN = 1, NKRN
          WV2 = WAVE0 + WLKRN(IKRN)
            HHLIN = 0.
            IF (ISELF .EQ. 1) CALL HGRID(HSELFFILE,
     ;              WV2,WAVE0,TD,H1FRCD,NLOW,NUP,HHLIN)
            IF (ISELF .EQ. 2) THEN
              CALL HSELF_PDWIDTH(LINE,TEMP,H1FRC,HWID)
              HWID = HWID * WAVE0 * WAVE0 / C / 2. / PI
              HHLIN = HWID / (HWID * HWID + WLKRN(IKRN) * WLKRN(IKRN)) 
     ,                     / PI
            END IF
            IF (ISELF .EQ. 3) THEN
              CALL ALI_GRIEM(LINE,WAVE0,H1FRC,RES)
              HHLIN = RES / (RES * RES + WLKRN(IKRN) * WLKRN(IKRN)) / PI
            END IF 
          KERN(IKRN) = HHLIN
        END DO
C
C  Do Convolution -- Integration is simple Trapezoidal here
C
        HLIND = 0.d0
        DO IKRN = 1, NKRN
          DEL1 = DEL + WLKRN(IKRN)
          I1 = MAX0(   1, IKRN - 1)
          I2 = MIN0(NKRN, IKRN + 1)
		  DDEL = WLKRN(I2) - WLKRN(I1)
          IF (ISTARK .EQ. 1) THEN
            CALL STEHLE(HFILE,WAVE0+DEL1,WAVE0,TEMP,XXNE,NLOW,NUP,STARK)
          END IF
          IF (ISTARK .EQ. 2) THEN
            CALL VCS(STARKD,XNED,TD,DEL1,1,NLOW,NUP,HVCSFILE)
            STARK = STARKD
          END IF
          F1 = KERN(IKRN) * STARK
          HLIND = HLIND + DDEL * F1
        END DO
        HLIND= HLIND * .5d0
        HLIN = HLIND
      END IF
C
C The He and radiative broadening are always assumed small enough to 
C convolve by adding 
C
      HEWID = 0.
      IF (IHE .EQ. 1 ) CALL HE_PDWIDTH(LINE,TEMP,HE1FRC,HEWID)
      CALL RAD_WIDTH(LINE,RAD)
C
C  Compute the width in lambda units (AA)
C
      LOR = RAD + HEWID
      LOR = LOR * WAVE0 * WAVE0 / C / 2. / PI
C
      IF (DEL .GT. .2)
     ;  HLIN = HLIN + LOR / (LOR * LOR + DEL * DEL) / PI
C
      PROF = HLIN
      RETURN
C
      END

      SUBROUTINE STEHLE(HFILE,WAVE,WAVE0,TEMP,XXNE,NLOW,NUP,HLIN)
C
C  N Piskunov's code to return Stehle Stark profiles
C
      DOUBLE PRECISION WAVE,WAVE0
      PARAMETER (NWL=50,NTEMP=10,NNE=15,NLINE=6,SQRTPI=1.77245385)
      REAL T(NTEMP),XNE(NNE)
      REAL WING(NNE,NLINE),F0(NNE,NLINE),WIDTH(NTEMP,NNE,NLINE),
     *     ALPHA(NWL,NTEMP,NNE,NLINE),SPROF(NWL,NTEMP,NNE,NLINE)
      INTEGER NALPHA(NTEMP,NNE,NLINE),NWL1,NTEMP1,NNE1,NLINE1,
     *     N1(NLINE),M1(NLINE)
      INTEGER*2 LINDEX(100,100)
      LOGICAL FIRST
      CHARACTER HFILE*(*)
      SAVE T,XNE,WING,F0,WIDTH,ALPHA,SPROF,NALPHA,
     *     NWL1,NTEMP1,NNE1,NLINE1,FIRST,LINDEX
      DATA FIRST/.TRUE./,LINDEX/10000*0/
      EXP10(X)=EXP(2.30258509299405E0*X)
C
C  Read in the table (this is done only once)
C
      IF(FIRST) THEN
        FIRST=.FALSE.
C
C  Open file and process errors
C
        OPEN(1,file=HFILE,IOSTAT=IERR,FORM='UNFORMATTED',STATUS='OLD')
        IF(IERR.NE.0) THEN
          WRITE(*,*) 'ERROR: STEHLE did not find Hydrogen line file'
          WRITE(*,*) HFILE
          STOP
        ENDIF
C
C  Read table dimensions and process errors
C
        READ(1,ERR=9,END=10) NWL1,NTEMP1,NNE1,NLINE1
        IF(NWL1.GT.NWL) THEN
          WRITE(*,*) 'ERROR: NWL must be at least',NWL1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
        IF(NTEMP1.GT.NTEMP) THEN
          WRITE(*,*) 'ERROR: NTEMP must be at least',NTEMP1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
        IF(NNE1.GT.NNE) THEN
          WRITE(*,*) 'ERROR: NNE must be at least',NNE1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
        IF(NLINE1.GT.NLINE) THEN
          WRITE(*,*) 'ERROR: NLINE must be at least',NLINE1,
     *               ' to hold the Hydrogen line table in STEHLE'
          STOP
        ENDIF
C
C  Read line ID and construct line index
C
        READ(1,ERR=9,END=10) (N1(LINE),M1(LINE),LINE=1,NLINE1)
        DO 1 LINE=1,NLINE1
   1    LINDEX(N1(LINE),M1(LINE))=LINE
C
C  Read temperatures and electron densities
C
        READ(1,ERR=9,END=10) (T(I),I=1,NTEMP1)
        READ(1,ERR=9,END=10) (XNE(J),J=1,NNE1)
C
C  Read wavelength conversion factor, wing extrapolation parameter
C  and line widths
C
        READ(1,ERR=9,END=10) ((F0(J,L),J=1,NNE1),L=1,NLINE1)
        READ(1,ERR=9,END=10) ((WING(J,L),J=1,NNE1),L=1,NLINE1)
        READ(1,ERR=9,END=10) (((WIDTH(I,J,L),I=1,NTEMP1),J=1,NNE1),
     *                       L=1,NLINE1)
C
C  Read the number of detuning points and the two big arrays
C  for detuning (ALPHA) and line profiles (SPROF)
C
        READ(1,ERR=9,END=10) (((NALPHA(I,J,L),I=1,NTEMP1),J=1,NNE1), 
     *                       L=1,NLINE1)
        DO 2 L=1,NLINE1
        READ(1,ERR=9,END=10) (((ALPHA(K,I,J,L),K=1,NWL1),
     *                       I=1,NTEMP1),J=1,NNE1)
   2    READ(1,ERR=9,END=10) (((SPROF(K,I,J,L),K=1,NWL1),
     *                       I=1,NTEMP1),J=1,NNE1)
        CLOSE(1)
        IF(LINDEX(NLOW,NUP).GT.0) WRITE(*,*) 
     *    'Will use Stehle''s tables for Hydrogen line profiles'
      END IF
C
C  Interpolation section
C
      LINE=LINDEX(NLOW,NUP)
C
C  Check if line is in the table
C
      IF(LINE.EQ.0) THEN
        HLIN=-1.
        RETURN
      END IF
C
C  Temperature interpolation (no extrapolation)
C
      IF(TEMP.LT.T(1).OR.TEMP.GT.T(NTEMP1)) THEN
        HLIN=-1.
        RETURN
      END IF
      ITEMP=2
   3  IF(T(ITEMP).LT.TEMP.AND.ITEMP.LT.NTEMP1) THEN
        ITEMP=ITEMP+1
        GO TO 3
      END IF
      DT=LOG10(T(ITEMP))-LOG10(T(ITEMP-1))
      DT1=(LOG10(TEMP)-LOG10(T(ITEMP-1)))/DT
      DT2=(LOG10(T(ITEMP))-LOG10(TEMP))/DT
C
C   Electron density interpolation (no extrapolation)
C
      ANE=LOG10(XXNE)
      IF(ANE.LT.XNE(1).OR.ANE.GT.XNE(NNE1)) THEN
        HLIN=-1.
        RETURN
      END IF
      INE=2
   4  IF(XNE(INE).LT.ANE.AND.INE.LT.NNE1) THEN
        INE=INE+1
        GO TO 4
      END IF
      DNE=XNE(INE)-XNE(INE-1)
      DNE1=(ANE-XNE(INE-1))/DNE
      DNE2=(XNE(INE)-ANE)/DNE
C
C  Interpolate over alpha. For each IT and IN we consider 3 cases:
C  1. alpha is in the first interval - simple linear interpolation
C  2. alpha is outside the table     - extrapolate using wing parameter
C  3. alpha is in between            - log interpolation
C
      IT=ITEMP-1
      IN=INE-1
      ALP=ABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.ALPHA(2,IT,IN,LINE)) THEN
        DAL=ALPHA(2,IT,IN,LINE)-ALPHA(1,IT,IN,LINE)
        DAL1=(ALP-ALPHA(1,IT,IN,LINE))/DAL
        DAL2=(ALPHA(2,IT,IN,LINE)-ALP)/DAL
        S11=DAL2*SPROF(1,IT,IN,LINE)+DAL1*SPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.ALPHA(NALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
        S11=WING(IN,LINE)/(ALP*ALP*SQRT(ALP))*SQRT(WAVE/WAVE0)
      ELSE
        ALP=LOG10(ALP)
        K=2
   5    IF(LOG10(ALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.NALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 5
        END IF
        DAL=LOG10(ALPHA(K,IT,IN,LINE))-LOG10(ALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-LOG10(ALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(LOG10(ALPHA(K,IT,IN,LINE))-ALP)/DAL
        S11=DAL2*LOG10(SPROF(K-1,IT,IN,LINE))+
     +      DAL1*LOG10(SPROF(K,  IT,IN,LINE))
        S11=EXP10(S11)
      END IF
      S11=LOG10(S11/F0(IN,LINE))
C
      IT=ITEMP
      IN=INE-1
      ALP=ABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.ALPHA(2,IT,IN,LINE)) THEN
        DAL=ALPHA(2,IT,IN,LINE)-ALPHA(1,IT,IN,LINE)
        DAL1=(ALP-ALPHA(1,IT,IN,LINE))/DAL
        DAL2=(ALPHA(2,IT,IN,LINE)-ALP)/DAL
        S21=DAL2*SPROF(1,IT,IN,LINE)+DAL1*SPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.ALPHA(NALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
        S21=WING(IN,LINE)/(ALP*ALP*SQRT(ALP))*SQRT(WAVE/WAVE0)
      ELSE
        ALP=LOG10(ALP)
        K=2
   6    IF(LOG10(ALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.NALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 6
        END IF
        DAL=LOG10(ALPHA(K,IT,IN,LINE))-LOG10(ALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-LOG10(ALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(LOG10(ALPHA(K,IT,IN,LINE))-ALP)/DAL
        S21=DAL2*LOG10(SPROF(K-1,IT,IN,LINE))+
     +      DAL1*LOG10(SPROF(K,  IT,IN,LINE))
        S21=EXP10(S21)
      END IF
      S21=LOG10(S21/F0(IN,LINE))
C
      IT=ITEMP-1
      IN=INE
      ALP=ABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.ALPHA(2,IT,IN,LINE)) THEN
        DAL=ALPHA(2,IT,IN,LINE)-ALPHA(1,IT,IN,LINE)
        DAL1=(ALP-ALPHA(1,IT,IN,LINE))/DAL
        DAL2=(ALPHA(2,IT,IN,LINE)-ALP)/DAL
        S12=DAL2*SPROF(1,IT,IN,LINE)+DAL1*SPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.ALPHA(NALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
        S12=WING(IN,LINE)/(ALP*ALP*SQRT(ALP))*SQRT(WAVE/WAVE0)
      ELSE
        ALP=LOG10(ALP)
        K=2
   7    IF(LOG10(ALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.NALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 7
        END IF
        DAL=LOG10(ALPHA(K,IT,IN,LINE))-LOG10(ALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-LOG10(ALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(LOG10(ALPHA(K,IT,IN,LINE))-ALP)/DAL
        S12=DAL2*LOG10(SPROF(K-1,IT,IN,LINE))+
     +      DAL1*LOG10(SPROF(K,  IT,IN,LINE))
        S12=EXP10(S12)
      END IF
      S12=LOG10(S12/F0(IN,LINE))
C
      IT=ITEMP
      IN=INE
      ALP=ABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.ALPHA(2,IT,IN,LINE)) THEN
        DAL=ALPHA(2,IT,IN,LINE)-ALPHA(1,IT,IN,LINE)
        DAL1=(ALP-ALPHA(1,IT,IN,LINE))/DAL
        DAL2=(ALPHA(2,IT,IN,LINE)-ALP)/DAL
        S22=DAL2*SPROF(1,IT,IN,LINE)+DAL1*SPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.ALPHA(NALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
        S22=WING(IN,LINE)/(ALP*ALP*SQRT(ALP))*SQRT(WAVE/WAVE0)
      ELSE
        ALP=LOG10(ALP)
        K=2
   8    IF(LOG10(ALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.NALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 8
        END IF
        DAL=LOG10(ALPHA(K,IT,IN,LINE))-LOG10(ALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-LOG10(ALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(LOG10(ALPHA(K,IT,IN,LINE))-ALP)/DAL
        S22=DAL2*LOG10(SPROF(K-1,IT,IN,LINE))+
     +      DAL1*LOG10(SPROF(K,  IT,IN,LINE))
        S22=EXP10(S22)
      END IF
      S22=LOG10(S22/F0(IN,LINE))
C
C  Interpolate over T and Ne
C
      HLIN=EXP10(DT2*DNE2*S11+DT1*DNE2*S21+DT2*DNE1*S12+DT1*DNE1*S22)
      RETURN
C
C  I/O error processing
C
   9  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'ERROR reading binary file '//HFILE(1:I)//' in HTABLE'
      STOP
  10  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'EOF found while reading binary file '//HFILE(1:I)//
     *           ' in HTABLE'
      STOP
      END

      SUBROUTINE HGRID(HFILE,WAVE,WAVE0,TEMP,NHYD,NLOW,NUP,HLIN)
C
C  N Piskunov's code to return Stehle Stark profiles
C  ADAPTED BY PAUL BARKLEM FOR MY GRIDS OF HYDROGEN PROFILES
C  ALL DOUBLE-PRECISION NOW ALSO
C  MAIN THING IS PROFILES ARE NOT ASSUMED SYMMETRIC
C
C  THERE IS NO WING APPROX YET, JUST SET TO A VERY SMALL NUMBER
C
      IMPLICIT NONE
      REAL*8 WAVE,WAVE0
      INTEGER MNPROFS,MNTEMP,MNNH,MNLINE
      PARAMETER (MNTEMP=10,MNNH=30,MNLINE=10,MNPROFS=200)
      REAL*8 T(MNTEMP),NH(MNNH),NHYD,TEMP,HLIN
      REAL*8 F0(MNNH,MNLINE),
     *     MALPHA(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     PALPHA(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     FALPHA(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     MSPROF(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     PSPROF(MNPROFS,MNTEMP,MNNH,MNLINE),
     *     FSPROF(MNPROFS,MNTEMP,MNNH,MNLINE)
      INTEGER MNALPHA(MNTEMP,MNNH,MNLINE),
     *     PNALPHA(MNTEMP,MNNH,MNLINE),
     *     FNALPHA(MNTEMP,MNNH,MNLINE)
      REAL*8 ALP,DAL,DAL1,DAL2,S11,S12,S21,S22
      REAL*8 DT,DT1,DT2,DNH,DNH1,DNH2
      REAL*8 SQRTPI
      INTEGER NPROFS,NTEMP,NNH,NLINE
      INTEGER NU(MNLINE),NL(MNLINE),NLOW,NUP
      INTEGER IERR,I,K,LINE,IN,IT,ITEMP,INH
      INTEGER MIDDLE
      INTEGER*2 LINDEX(100,100)
      LOGICAL FIRST
      CHARACTER HFILE*(*)
      SAVE T,NH,F0,MALPHA,MSPROF,PALPHA,PSPROF,PNALPHA,MNALPHA,
     *     NPROFS,NTEMP,NNH,NLINE,FIRST,LINDEX
      DATA FIRST/.TRUE./,LINDEX/10000*0/
      PARAMETER (SQRTPI=1.77245385)
      REAL*8 DEXP10,X
      DEXP10(X)=EXP(2.30258509299405D0*X)
C
C  Read in the table (this is done only once)
C
      IF(FIRST) THEN
        FIRST=.FALSE.
C
C  Open file and process errors
C
        OPEN(1,file=HFILE,IOSTAT=IERR,FORM='UNFORMATTED',STATUS='OLD')
        IF(IERR.NE.0) THEN
          WRITE(*,*) 'ERROR: HGRID did not find Hydrogen line file'
          WRITE(*,*) HFILE
          STOP
        ENDIF
C
C  Read table dimensions and process errors
C
        READ(1,ERR=9,END=10) NLINE
C
C  Read line ID and construct line index, read grid parameters
C
        READ(1,ERR=9,END=10) (NL(LINE),NU(LINE),LINE=1,NLINE)
        DO 401 I=1,NLINE
 401       LINDEX(NL(I),NU(I))=I
        READ(1,ERR=9,END=10) NNH
        READ(1,ERR=9,END=10) (NH(I),I=1,NNH)
        READ(1,ERR=9,END=10) NTEMP
        READ(1,ERR=9,END=10) (T(I),I=1,NTEMP)
C
C  Read the number of detuning points and the two big arrays
C  for detuning (ALPHA) and line profiles (SPROF)
C  Right Now ALPHA=LAMBDA, SPROF=INTENSITY (normalised)
C  hence F0=1 everywhere 
C
        READ(1,ERR=9,END=10) NPROFS
        MIDDLE=(NPROFS-1)/2+1
        DO LINE=1,NLINE
           DO IN=1,NNH
              DO IT=1,NTEMP
                 DO I=1,NPROFS
          READ(1,ERR=9,END=10) FALPHA(I,IT,IN,LINE),FSPROF(I,IT,IN,LINE)
                 END DO
              END DO
           END DO
        END DO
        DO LINE=1,NLINE
           DO IN=1,NNH
              F0(IN,LINE)=1.d0
           END DO
        END DO
        DO LINE=1,NLINE
           DO IN=1,NNH
              DO ITEMP=1,NTEMP
                 FNALPHA(ITEMP,IN,LINE)=NPROFS
                 PNALPHA(ITEMP,IN,LINE)=MIDDLE
                 MNALPHA(ITEMP,IN,LINE)=MIDDLE
              END DO
           END DO
        END DO
C
C  SET UP SEPARATE ARRAYS OF EACH SIDE OF THE PROFILE
C
        DO LINE=1,NLINE
           DO IN=1,NNH
              DO IT=1,NTEMP
                 DO K=1,MIDDLE
                    PALPHA(K,IT,IN,LINE)=FALPHA(MIDDLE+K-1,IT,IN,LINE)
              MALPHA(K,IT,IN,LINE)=DABS(FALPHA(MIDDLE-K+1,IT,IN,LINE))
                    PSPROF(K,IT,IN,LINE)=FSPROF(MIDDLE+K-1,IT,IN,LINE)
                    MSPROF(K,IT,IN,LINE)=FSPROF(MIDDLE-K+1,IT,IN,LINE)
                 END DO
              END DO
           END DO
        END DO                    
C
C  IF OUR LINE EXIST THEN WRITE OUT FIRST TIME
C
        IF(LINDEX(NLOW,NUP).GT.0) WRITE(*,*) 
     *    'Will use PB''s tables for Hydrogen line profiles'
C
      END IF  ! ENDS PART ONLY DONE FIRST TIME
C
C  Interpolation section
C
      LINE=LINDEX(NLOW,NUP)
C
C  Check if line is in the table
C
      IF(LINE.EQ.0) THEN
        HLIN=-1.
        RETURN
      END IF
C
C  Temperature and perturber density interpolation (no extrapolation)
C
      IF(TEMP.LT.T(1).OR.TEMP.GT.T(NTEMP)) THEN
        HLIN=-1.
        RETURN
      END IF
      ITEMP=2
   3  IF(T(ITEMP).LT.TEMP.AND.ITEMP.LT.NTEMP) THEN
        ITEMP=ITEMP+1
        GO TO 3
      END IF
      DT=DLOG10(T(ITEMP))-DLOG10(T(ITEMP-1))
      DT1=(DLOG10(TEMP)-DLOG10(T(ITEMP-1)))/DT
      DT2=(DLOG10(T(ITEMP))-DLOG10(TEMP))/DT
C
      IF(NHYD.LT.NH(1).OR.NHYD.GT.NH(NNH)) THEN
        HLIN=-1.
        RETURN
      END IF
      INH=2
   4  IF(NH(INH).LT.NHYD.AND.INH.LT.NNH) THEN
        INH=INH+1
        GO TO 4
      END IF
      DNH=DLOG10(NH(INH))-DLOG10(NH(INH-1))
      DNH1=(DLOG10(NHYD)-DLOG10(NH(INH-1)))/DNH
      DNH2=(DLOG10(NH(INH))-DLOG10(NHYD))/DNH
C
C  WE MUST SET UP TO ALLOW FOR THE FACT THAT THE PROFILE IS NOT ASSUMED 
C  TO BE SYMMETRIC
C
C  Interpolate over alpha. For each IT and IN we consider 3 cases:
C  1. alpha is in the first interval - simple linear interpolation
C  2. alpha is outside the table     - extrapolate using wing parameter
C  3. alpha is in between            - log interpolation
C
c  no wing extrapolation (impact approx invalid), we just take it as zero.
c
      IF ((WAVE-WAVE0).GE.0.d0) THEN
C
C  CASE OF POSITIVE DETUNING
C
      IT=ITEMP-1
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S11=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S11=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   5    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 5
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S11=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S11=DEXP10(S11)
      END IF
      if (S11.ne.0.d0) S11=DLOG10(S11/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S21=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S21=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   6    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 6
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S21=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S21=DEXP10(S21)
      END IF
      if (S21.ne.0.d0) S21=DLOG10(S21/F0(IN,LINE))
C
      IT=ITEMP-1
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S12=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S12=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   7    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 7
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S12=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S12=DEXP10(S12)
      END IF
      if (S12.ne.0.d0) S12=DLOG10(S12/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.PALPHA(2,IT,IN,LINE)) THEN
        DAL=PALPHA(2,IT,IN,LINE)-PALPHA(1,IT,IN,LINE)
        DAL1=(ALP-PALPHA(1,IT,IN,LINE))/DAL
        DAL2=(PALPHA(2,IT,IN,LINE)-ALP)/DAL
        S22=DAL2*PSPROF(1,IT,IN,LINE)+DAL1*PSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.PALPHA(PNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S22=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
   8    IF(DLOG10(PALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.PNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 8
        END IF
        DAL=DLOG10(PALPHA(K,IT,IN,LINE))-DLOG10(PALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(PALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(PALPHA(K,IT,IN,LINE))-ALP)/DAL
        S22=DAL2*DLOG10(PSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(PSPROF(K,  IT,IN,LINE))
        S22=DEXP10(S22)
      END IF
      if (S22.ne.0.d0) S22=DLOG10(S22/F0(IN,LINE))
C
      ELSE
C
C  NOW THE CASE OF NEGATIVE ALPHA
C
      IT=ITEMP-1
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S11=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S11=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  15    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 15
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S11=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S11=DEXP10(S11)
      END IF
      if (S11.ne.0.d0) S11=DLOG10(S11/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH-1
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S21=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S21=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  16    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 16
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S21=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S21=DEXP10(S21)
      END IF
      if (S21.ne.0.d0) S21=DLOG10(S21/F0(IN,LINE))
C
      IT=ITEMP-1
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S12=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S12=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  17    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 17
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S12=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S12=DEXP10(S12)
      END IF
      if (S12.ne.0.d0) S12=DLOG10(S12/F0(IN,LINE))
C
      IT=ITEMP
      IN=INH
      ALP=DABS(WAVE-WAVE0)/F0(IN,LINE)
      IF(ALP.LT.MALPHA(2,IT,IN,LINE)) THEN
        DAL=MALPHA(2,IT,IN,LINE)-MALPHA(1,IT,IN,LINE)
        DAL1=(ALP-MALPHA(1,IT,IN,LINE))/DAL
        DAL2=(MALPHA(2,IT,IN,LINE)-ALP)/DAL
        S22=DAL2*MSPROF(1,IT,IN,LINE)+DAL1*MSPROF(2,IT,IN,LINE)
      ELSE IF(ALP.GT.MALPHA(MNALPHA(IT,IN,LINE),IT,IN,LINE)) THEN
         S22=0.d0
      ELSE
        ALP=DLOG10(ALP)
        K=2
  18    IF(DLOG10(MALPHA(K,IT,IN,LINE)).LT.ALP.AND.
     *     K.LT.MNALPHA(IT,IN,LINE)) THEN
          K=K+1
          GO TO 18
        END IF
        DAL=DLOG10(MALPHA(K,IT,IN,LINE))-DLOG10(MALPHA(K-1,IT,IN,LINE))
        DAL1=(ALP-DLOG10(MALPHA(K-1,IT,IN,LINE)))/DAL
        DAL2=(DLOG10(MALPHA(K,IT,IN,LINE))-ALP)/DAL
        S22=DAL2*DLOG10(MSPROF(K-1,IT,IN,LINE))+
     +      DAL1*DLOG10(MSPROF(K,  IT,IN,LINE))
        S22=DEXP10(S22)
      END IF
      if (S22.ne.0.d0) S22=DLOG10(S22/F0(IN,LINE))
C
      END IF
C
C  Interpolate over T and NH
C
      if (S11.eq.0.d0.and.S21.eq.0.d0.and.S12.eq.0.d0.and.S22.eq.0.d0)
     * then
         HLIN=0.d0
      else
      HLIN=DEXP10(DT2*DNH2*S11+DT1*DNH2*S21+DT2*DNH1*S12+DT1*DNH1*S22)
      end if
      RETURN
C
C  I/O error processing
C
   9  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'ERROR reading binary file '//HFILE(1:I)//' in HGRID'
      STOP
  10  I=INDEX(HFILE,' ')-1
      IF(I.LE.0) I=LEN(HFILE)
      WRITE(*,*) 'EOF found while reading binary file '//HFILE(1:I)//
     *           ' in HGRID'
      STOP
      END


      SUBROUTINE HSELF_PDWIDTH(LINE,TEMP,H1FRC,WIDTH)
C
C   Computes the Lorentzian WIDTH self broadening 
C   of the Hydrogen Balmer Lines in the p-d approximation
C   Output in usual rad /s cm^3
C   
C   Here LINE = 1,2,3... for Halpha, Hbeta, Hgamma ...
C
      IMPLICIT NONE
      REAL*4 H1FRC
      REAL*4 C,PI,KB,M0,A0,SIGMAH(4),ALPHAH(4),
     *       X,GAMMAF,GX,VBAR,SIGMA,TEMP
      REAL*8 GVW, WIDTH
      INTEGER LINE
C
      PARAMETER (C=2.997925E+18,PI=3.14159265)
      PARAMETER (KB=1.380658E-23,M0=1.660540E-27)
      PARAMETER (A0=5.29177249E-11)  
      DATA SIGMAH/ 1180., 2320., 4208., 0.0/
      DATA ALPHAH/ 0.677, 0.455, 0.380, 0.0/
C  
      X = 2. - ALPHAH(LINE) * .5
      GX = X - 1.0
      GAMMAF = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX
     ;           )*GX)*GX)*GX)*GX
      SIGMA = SIGMAH(LINE) * A0 * A0
      GVW = (4./PI)**(ALPHAH(LINE)*0.5)*GAMMAF*1.E4*SIGMA
      VBAR = SQRT(8.*KB*TEMP/PI/M0 * 2./1.008)
      GVW = GVW * ((VBAR/1.E4)**(1.-ALPHAH(LINE)))
      GVW = GVW * H1FRC * 1.E6
C 
      WIDTH=GVW
      RETURN
C
      END

      SUBROUTINE HE_PDWIDTH(LINE,TEMP,HE1FRC,HEWID)
C 
C  Compute approximate broadening by He collisions (p-d approx)
C
C  Output in usual rad /s cm^3
C   
C  27/4/00 PB
C
      IMPLICIT NONE
      REAL*4 PI,C,HE1FRC,GAMMAF
      REAL*4 SIGMAH(4),ALPHAH(4),X,GX,A0,SIGMA,VBAR,KB,M0,TEMP
      REAL*4 HEWID,GVW
      INTEGER LINE
      PARAMETER (C=2.997925E+18,PI=3.14159265,A0=5.29177249E-11)    
      PARAMETER (KB=1.380658E-23,M0=1.660540E-27)
      DATA SIGMAH/ 834., 1998., 3837., 5700. /
      DATA ALPHAH/ .280, .330, .260, .250 /
C
      X = 2. - ALPHAH(LINE) * .5
      GX = X - 1.0
      GAMMAF = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX
     ;           )*GX)*GX)*GX)*GX
      SIGMA = SIGMAH(LINE) * A0 * A0
      GVW = (4./PI)**(ALPHAH(LINE)*0.5)*GAMMAF*1.E4*SIGMA
      VBAR = SQRT(8.*KB*TEMP/PI/M0 * (1./1.008 + 1./4.003))
      HEWID = GVW * ((VBAR/1.E4)**(1.-ALPHAH(LINE))) 
      HEWID = HEWID * 0.625            !polarisability difference term
      HEWID = HEWID * HE1FRC * 1.E6    !in cgs
C
      RETURN
      END

      SUBROUTINE RAD_WIDTH(LINE,WIDTH)
C
C  Returns the radiation damping width (hwhm) in rad/s cm^3
C
C  27/4/00 PB
C
      INTEGER LINE
      REAL*4 WIDTH,GAMRAD(4)
      DATA GAMRAD/ 5.698, 5.001, 4.815, 4.740 / 
C
      WIDTH = GAMRAD(LINE) * 1.e8 / 2.
C
      RETURN
      END


      SUBROUTINE ALI_GRIEM(LINE,WAVE0,H1FRC,WIDTH)
C
C   Computes the Resonance Width (AA) of the Hydrogen Balmer Lines
C   Halfwidth given H1FRC perturbers per cm^3 
C
C   NB: assumes 2p level only! more correct than adding the 3p or 4p etc.
C   since 2p-nd dominates
C   so LINE is not used actually
C
C  GAMMA = lambda^2/(2 pi c) * 1.15 * e^2/(mc) f12 lambda12 N(H)   fuhrmann
C  GAMMA = lambda^2/(2 pi c) * 1.108 * e^2/(mc) f12 lambda12 N(H)   ali and griem
C
      IMPLICIT NONE
      REAL*4 H1FRC
      REAL*8 WAVE0,WIDTH,FW
      INTEGER LINE
C
C     WAVE0 in AA, WIDTH in AA  
C
C   the two different leading numbers makes only a minute difference
C     FW=2.6037D-27*WAVE0*WAVE0   !balmer lower level (fuhrmann 2 pi)
C
      FW=2.516D-27*WAVE0*WAVE0   !balmer lower level (ali and griem 1.92 pi)
      FW=FW*H1FRC
C      
      WIDTH=FW/2.d0
C
      RETURN
C
      END


C
C  THE FOLLOWING IS KURUCZ'S SUBROUTINE
C
C  SLIGHTLY ALTERED BY PB & REGNER TRAMPEDACH
C  UNIT 25->1 NEEDED TO WORK WITH MULTI
C
      SUBROUTINE VCS(PR,XNE,T,DEL,II,N,M,HVCSFILE)     
C     CALCULATES VIDAL COOPER AND SMITH PROFILES FOR FIRST FOUR BALMER 
C     LINES.  THEY ARE RETURNED IN ANGSTROM UNITS.                     
C     ASSUMES RAW PROFILES ARE IN ALPHA UNITS.                      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PR(40),DEL(40)                                       
C     DIMENSION ALPHA(40),PRALPH(40)                                  
      DIMENSION PRALPH(40)                                   
      DIMENSION SVCS(6,17,40,4),ALPHA0(4)  
      CHARACTER*80 HVCSFILE                                
      DATA SVCS(1,1,1,1)/0./,ALPHA0/-3.,-3.,-3.,-3./   
      SAVE SVCS
C
      EXP10(X)=EXP(2.30258509299405E0*X)                       
      IF(SVCS(1,1,1,1).NE.0.)GO TO 3
C
C     READ IN VCS ARRAYS                                          
      OPEN(1,FILE=HVCSFILE,FORM='FORMATTED',STATUS='OLD')
      READ(1,10)
      DO 20 LINE=1,4
        READ(1,10)
        READ(1,10)(((SVCS(I,J,K,LINE),J=1,17),I=1,6),K=1,40)
   10   FORMAT (10F8.4)
   20 CONTINUE
      CLOSE(1)
C
    3 LINE=M-N                                  
C     TEMPERATURE AND ELECTRON DENSITY INTERPOLATION        
      AT=LOG10(T)                                                              
      BTEMP=(AT-3.3979411)/.3010300+1.                     
      ITEMP=BTEMP                                       
C     ITEMP=(LOG10(T)-LOG10(2500.))/LOG10(2.)+1.                             
      ITEMP=MAX(MIN(ITEMP,5),1)
      WTTEMP=BTEMP-ITEMP                                   
      ANE=LOG10(XNE)                                                           
C     TO AVOID BAD EXTRAPOLATION AT LOW DENSITY
      ANE=MAX(10.D0,ANE)
      BNE=(ANE-10.)/.5+1.                  
      INE=BNE                      
      INE=MAX(MIN(INE,16),1)
      WTXNE=BNE-INE                             
      DO 21 I=1,40                                                
   21 PRALPH(I)=(1.-WTXNE)*(1.-WTTEMP)*SVCS(ITEMP  ,INE  ,I,LINE)+
     1               (1.-WTXNE)*WTTEMP*SVCS(ITEMP+1,INE  ,I,LINE)+           
     2               WTXNE*(1.-WTTEMP)*SVCS(ITEMP  ,INE+1,I,LINE)+
     3                    WTXNE*WTTEMP*SVCS(ITEMP+1,INE+1,I,LINE)   
C     NOW ALPHA INTERPOLATION                                             
      FO=1.25E-9*XNE**.66666667                                    
      DO 50 I=1,II                                         
      IF (DEL(I).NE.0.) GO TO 25                            
      PR(I)=EXP10(PRALPH(1))/FO                               
      GO TO 50                     
   25 AALP=LOG10(ABS(DEL(I))/FO)                                               
      BALP=(AALP-ALPHA0(LINE))/.2+1.          
      IALP=BALP                                       
      IF (IALP.GT.1) GO TO 30                 
      PR(I)=EXP10(PRALPH(1))/FO           
      GO TO 50                                                  
   30 IF (IALP.GT.39) GO TO 40                
      WT=BALP-IALP                        
      PR(I)=EXP10((1.-WT)*PRALPH(IALP)+WT*PRALPH(IALP+1))/FO    
      GO TO 50                                        
   40 PR(I)=EXP10(PRALPH(40)+2.5*(ALPHA0(LINE)+7.8-AALP))/FO      
   50 CONTINUE
      RETURN                       
      END  


