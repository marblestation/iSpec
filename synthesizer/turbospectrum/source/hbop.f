C***********************************************************************
C  Hydrogen bound opacity.  Designed to give the total neutral H opacity 
C  (bound-bound and bound free) at a given wavelength for a given list of 
C  lines at specified plasma conditions.  Occupation probability formalism
C  is used. Appropriate for applications like model stellar atmosphere 
C  calculations.
C
C  Also contains the code for normalised line profiles HLINOP which can 
C  be used separately.
C
C  Paul Barklem, Nik Piskunov, Uppsala August 2003. 
C
C  A number of parts from or based on code by Deane Peterson and Bob 
C  Kurucz. 
C
C  Some contributions and significant input from Kjell Eriksson.
C  
C  Modified by B. Plez to include detailed treatment of lower Balmer lines
C  April 2, 2019
C
C  We would very much appreciate bug reports to: barklem@astro.uu.se 
C
C  Table of Contents:
C  ------------------
C 
C  REAL HBOP(WAVE,N,NLO,NUP,WAVEH,NH,NHE,NE,T,DOP,NPOP,NL,TOTAL,CONTIN)
C  REAL FUNCTION HLINOP(WAVE,NBLO,NBUP,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPLE)
C  REAL FUNCTION HPROFL(N,M,WAVE,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPH)
C  REAL FUNCTION STARK1(N,M,WAVE,WAVEH,T,XNE)
C  REAL FUNCTION HFNM(N,M)
C  REAL FUNCTION VCSE1F(X)
C  REAL FUNCTION SOFBET(B,P,N,M)
C  REAL*8 FUNCTION AIRVAC(W)
C  REAL*8 FUNCTION VACAIR(W)     
C  REAL FUNCTION WCALC(NH, NE, NHW, NS, TEMP)	
C  REAL FUNCTION HBF(WAVE,NLO)
C  SUBROUTINE LOGINT(X,Y,N,XI,YI)
C  SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C
C  Variables to be set / aware of before use:
C  ------------------------------------------
C
C  HBOP  : NLEVELS, NBF
C  HLINOP: CUTNEXT, REFLECT
C  HPROFL: RAYLCUT
C
C***********************************************************************

      SUBROUTINE HBOP(WAVE, N, NLO, NUP, WAVEH, NH, NHE, NE, T, DOP,
     *                  NPOP, NL, TOTAL, CONTIN, contonly)
C
C  Returns the total absorption coefficient due to H bound levels at WAVE, 
C  for a given line list employing the occupation probability formalism 
C  of Daeppen, Anderson & Milhalas ApJ 1987 (DAM).
C
C  The extension of the occupation probability formalism to non-LTE by 
C  Hubeny, Hummer and Lanz A&A 1994 (HHL) is implemented (though see 
C  comments below).
C
C  To use the code in LTE:
C         set NL = 0 
C         NPOP is then not used, population computed internally
C  To use the code in non-LTE:
C         set NL = number of states to treat in non-LTE
C         set NPOP = array of number densities for each state
C         states above NL will be treated in LTE if needed
C                            
C  Note that Hubeny, Hummer and Lanz 1994 have extended the occupation 
C  probability formalism to non-LTE, though in a somewhat unclear manner 
C  in my opinion.  The main point to be aware of is that w(HHL) does not 
C  equal w(DAM).  The first is a conditional probability and the second 
C  is unconditional.  The problem is that the definition in sect 2.1 of HHL 
C  and derivation for w in HHL Appendix A is for the unconditional 
C  probability w(DAM).  In the rest of the paper it's conditional. The paper 
C  doesn't really make this difference clear IMO. 
C
C  Thus w(j)(HHL) = P(j|i) = P(i & j)/P(i) = w(j)/w(i) for j>i
C                                          = 1         for j<i
C  where w(m) is as computed in App A which is as defined in DAM
C  and sect 2.1 of HHL.  Since HHL assume w(i)=1 this leads to no problems.
C  
C 
C  wave   = interrogation wavelength in A               real*8
C  n      = number of lines                             integer
C  nlo    = array of lower principal quantum numbers    int array
C  nup    = array of upper principal quantum numbers    int array
C  waveh  = array of central wavelengths                real*8 array 
C  t      = kinetic temperature in K
C  ne     = electron number density in cm-3
C  nh     = number density of H I in cm-3
C  nhe    = number density of He I in cm-3
C  dop    = reduced Doppler width delta_lambda / lambda_0
C         = reduced Doppler width delta_nu / nu_0
C  npop   = number density of each level in cm-3
C  nl     = number of levels for which populations are pre-specified 
C           -- higher levels are assumed in LTE
C  total  = returns the total (line + continuous) absorption coefficient
C  contin = returns the continuous absorption coefficient
C  contonly = computes only contin if true (BPz 03/04-2019)    logical
C
C  Important parameters
C  
C  nlevels = number of levels summed in partition functions and computed
C            populations and occupation probabilities. Must be larger or 
C            equal to the highest level in the line list and NL.
C  nbf     = highest level for which the bound-free component is 
C            accounted for (limited to 15 at present). 
C
C  Paul Barklem, Uppsala, August 2003
C
      IMPLICIT NONE
      INTEGER N, NLO(*), NUP(*), NLEVELS, I, J, NBF, NL, FF,FN
      PARAMETER (NLEVELS = 100)  
      PARAMETER (NBF = 6)  
      REAL NH, NHE, NE, T, CONTIN, DOP, TOTAL
      REAL NHL, NHEL, NEL, TL
      REAL*8 W(NLEVELS), G(NLEVELS), WGE(NLEVELS)
      REAL NLTE(NLEVELS), NPOP(*), NP(NLEVELS)
      REAL Z, H, C, HC, K, KT, LINE, SIGMA, CHI, SF, PROF
      REAL IONH, X, HFNM, FNM(NLEVELS,NLEVELS), HLINOP, HBF
      REAL TS, TF
      REAL*8 WAVE, WAVEH(*), REDCUT
      REAL*8 EHYD(NLEVELS), CONTH(NLEVELS), WCALC, D, WSTAR, TAPER
      LOGICAL FIRST,contonly
      SAVE 
      DATA FIRST/.TRUE./
      PARAMETER (H=6.62607E-27, C=2.9979E10, K=1.38065E-16)
      PARAMETER (IONH=2.17991e-11) 
      PARAMETER (HC=H*C)
      PARAMETER (SF=0.0265) ! pi*e*e/m/c [cgs]  sigma = sf * f * norm-prof
C
      IF (FIRST) THEN
C
C  Compute energy levels in the H atom, and store f-values
C
        EHYD(1) =      0.000D0
        EHYD(2) =  82259.105D0
        EHYD(3) =  97492.302D0
        EHYD(4) = 102823.893D0
        EHYD(5) = 105291.651D0
        EHYD(6) = 106632.160D0
        EHYD(7) = 107440.444D0
        EHYD(8) = 107965.051D0
        DO 1 I = 9, NLEVELS
 1      EHYD(I) = 109678.764D0 - 109677.576D0/(I*I)
        DO 2 I = 1, NLEVELS
 2      CONTH(I) = 109678.764D0 - EHYD(I)
        FIRST = .FALSE.
        DO 4 I = 1, NLEVELS-1
        DO 3 J = I+1, NLEVELS
 3      FNM(I,J) = HFNM(I,J)
 4      CONTINUE
      ENDIF 
C
C  Compute partition function and populations in LTE.
C  Need LTE populations even in non-LTE case
C  for calculating stimulated free-bound emission
C
C  Create population array NP for all levels 
C
      IF ((ABS(NH-NHL).LT.1E-3*NH).AND.(ABS(NHE-NHEL).LT.1E-3*NHE).AND.
     +   (ABS(NE-NEL).LT.1E-3*NE).AND.(ABS(T-TL).LT.1E-3*T))
     *   GOTO 16
      Z = 0.
      KT = K*T
      DO 10 I = 1, NLEVELS
      W(I) = WCALC(NH, NE, NHE, FLOAT(I), T)
      G(I) = 2.*I*I
      WGE(I) = W(I)*G(I)*DEXP(-HC*EHYD(I)/KT)
 10   Z = Z + WGE(I)
      DO 15 I = 1, NLEVELS
 15   NLTE(I) = NH*WGE(I)/Z  ! LTE populations
      NHL = NH
      NHEL = NHE
      NEL = NE
      TL = T
 16   CONTINUE
C
      DO 20 I = 1, NLEVELS
      IF (I .LE. NL) THEN 
        NP(I) = NPOP(I)        ! non-LTE populations for states below NL
      ELSE 
        NP(I) = NLTE(I)        ! LTE populations otherwise
      ENDIF
 20   CONTINUE

C
C  Compute line opacity components
C
      LINE = 0.
!
      if (.not. contonly) then
!
       DO 30 I = 1, N

        if (nlo(i).eq.2 .and. nup(i).le.5) then
! special treatment for lower Balmer lines
! BPz 2/04-2019
          call hlinprof(wave,waveh(i),t,ne,nlo(i),nup(i),
     &                   nh, nhe, dop,prof)
          prof=prof*wave*wave/c/1.e8
          sigma = sf * fnm(nlo(i),nup(i)) * prof
          chi=0.
          IF (W(NLO(I)).GT.0.)  ! populations would be zero also
     +       CHI = SIGMA * ( W(NUP(I))/W(NLO(I))*NP(NLO(I)) -
     -                 NP(NUP(I))*G(NLO(I))/G(NUP(I)) )  !stim em term
        else
! all other lines
          SIGMA = SF * FNM(NLO(I),NUP(I)) *
     *     HLINOP(WAVE,NLO(I),NUP(I),WAVEH(I),T,NE,NP(1),NHE,DOP)
          CHI = 0.
          IF (W(NLO(I)).GT.0.)  ! populations would be zero also
     +       CHI = SIGMA * ( W(NUP(I))/W(NLO(I))*NP(NLO(I)) -
     -                 NP(NUP(I))*G(NLO(I))/G(NUP(I)) )  !stim em term
        endif

 30    LINE = LINE + CHI
!
      endif
!
C
C  Compute the continuous components - sum over NBF lowest states
C     
      CONTIN = 0.
      DO 40 I = 1, NBF
C
C  Compute dissolved fraction for the wavelength for this lower level
C
C  D does not decline fast enough to overcome the bound-free extrapolation
C  and this can lead to large dissolved contibutions very far from 
C  the bound-free jump, particularly from the Lyman continuum in cool stars.
C  This is certainly unphysical as the extrapolation will be invalid 
C  for such cases (which corresponds basically to H2+ if the perturber 
C  is a proton).  I apply an arbitrary tapering redward of the alpha line.
C 
CC      CHI = 0.     
CC      IF (WAVE.GT.(1.d8/CONTH(I))) THEN 
      D = 1.
      TAPER = 1.
      X=1./(FLOAT(I)*FLOAT(I)) - HC/(WAVE*1.E-8*IONH)
      IF(X .GT. 0.) THEN
        X = 1./SQRT(X)
        WSTAR = WCALC(NH, NE, NHE, X, T)
        D = 1.0-(WSTAR/W(I))
        REDCUT = 1.D8/(EHYD(I+1)-EHYD(I))
        IF (WAVE .GT. REDCUT) TAPER = DEXP((REDCUT-WAVE)*1.D-2)
      ENDIF
      CHI = (NP(I) - NLTE(I)*DEXP(-HC/(WAVE*1E-8)/KT)) * 
     *       D * HBF(WAVE, I) * TAPER
CC      ENDIF
 40   CONTIN = CONTIN + CHI
C
C  Sum up and we're done
C
      TOTAL = CONTIN + LINE
C
      RETURN
      END

C***********************************************************************
      REAL FUNCTION HLINOP(WAVE,NBLO,NBUP,WAVEH,T,XNE,H1FRC,HE1FRC,
     *                     DOPPLE)
C
C Hydrogen line opacity
C
C hlinop = normalised line profile 
C          normalised to unity with integration over frequency
C wave   = interrogation wavelength in A           real*8
C nblo   = lower level principal quantum number
C nbup   = upper level principal quantum number
C waveh  = hydrogen line center wavelength in A    real*8
C t      = kinetic temperature in K
C xne    = electron number density in cm-3
C h1frc  = number density of H I in cm-3
C he1frc = number density of He I in cm-3
C dopple = reduced Doppler width delta_lambda / lambda_0
C        = reduced Doppler width delta_nu / nu_0
C
C  This wrapper program ensures lines are not used at unreasonably 
C  large detunings. 
C
      REAL*8 WAVEH,WAVE,WCON,WREFLECT
      REAL*8 REDCUT,BLUECUT
      REAL*8 EHYD(100),CONTH(100),WCUT(100)
      REAL*8 VACAIR
      LOGICAL FIRST,REFLECT,CUTNEXT
      SAVE
      DATA FIRST/.TRUE./
C
C  Set CUTNEXT true to cut lines at the centre of the next line. 
C  Otherwise they are cut at the series continuum edge in the blue and
C  at the red wing cut specified for alpha lines below. 
C
      DATA CUTNEXT/.FALSE./
C
C  Set REFLECT true to reflect profiles off the cutoffs (following 
C  Däppen et al 1987, ApJ 319, 195). The shape of the far wings should 
C  *not* be relied on in detail as state mixing in quasimolecular states
C  will determine potentials and therefore line shape, the opacity 
C  should perhaps be there somewhere (in the form of a satellite), and 
C  thus this *might* estimate this effect.
C
      DATA REFLECT/.FALSE./
C
C  Lyman series is treated in vacuum. The rest is converted to air.
C
      IF (NBLO.EQ.1) THEN
        IFVAC = 1
      ELSE
        IFVAC = 0
      END IF
C
      IF (FIRST) THEN
C
C  Compute energy levels in the H atom.
C
        EHYD(1) =      0.000D0
        EHYD(2) =  82259.105D0
        EHYD(3) =  97492.302D0
        EHYD(4) = 102823.893D0
        EHYD(5) = 105291.651D0
        EHYD(6) = 106632.160D0
        EHYD(7) = 107440.444D0
        EHYD(8) = 107965.051D0
        DO 1 I = 9, 100
 1      EHYD(I) = 109678.764D0 - 109677.576D0/I**2
        DO 2 I = 1, 100
 2      CONTH(I) = 109678.764D0 - EHYD(I)
C
C  Red cutoff wavelengths in Angstroms.
C  Arbitrarily chosen to be the same energy below the upper state of the
C  alpha line of the series as the next state is above it.
C
C  For Lyman alpha the H2 and H2+ satellites are well studied (Allard 
C  et al's papers ) and these are included, so just choose a suitable 
C  large number.
C 
C  For Balmer alpha the position of the main H2+ satellite is known 
C  (Kielkopf et al 2002, Eur Phys J D 18, 51) for densities 
C  N_H+ < 10^20 cm^-3 but the detailed line shape is not yet accounted 
C  for.
C
        WCUT(1) = 3647.D0
        WCUT(2) = 8650.D0
        DO 3 I = 3, 98
 3      WCUT(I) = 1.D8/((EHYD(I+1)-EHYD(I))-(EHYD(I+2)-EHYD(I+1)))
C
        FIRST = .FALSE.
      END IF
C
C  Compute series limit
C
        WCON = 1.D8/CONTH(NBLO)
        IF (IFVAC.EQ.0) WCON = VACAIR(WCON)
C
C  Either (depending on CUTNEXT setting) cut off at centre of next line
C  or: 
C  For red wings:
C    Cut at a predefined cutoff, and arbitrarily taper off exponentially 
C    from the cut (to avoid discontinuities in spectra).
C  For blue wings:
C    Cut at appropriate continuum edge.  
C
C  Reflect opacity in cutoff if REFLECT set.
C
        IF (WAVE.GE.WAVEH) THEN
C
C  Red wing
C
           HLINOP = 0.0
           REDCUT = WCUT(NBLO)
           IF (CUTNEXT.AND.(NBUP-NBLO.NE.1))
     ;               REDCUT = 1.D8/(EHYD(NBUP-1)-EHYD(NBLO))
           IF (IFVAC.EQ.0) REDCUT = VACAIR(REDCUT)
           IF (WAVE.LE.REDCUT) THEN
              HLINOP = HPROFL(NBLO,NBUP,WAVE,WAVEH,T,XNE,H1FRC,HE1FRC,
     ;                     DOPPLE)
              IF (REFLECT.AND.(NBUP-NBLO.NE.1)) THEN
                 WREFLECT = REDCUT+(REDCUT-WAVE)
                HLINOP = HLINOP + HPROFL(NBLO,NBUP,WREFLECT,WAVEH,T,XNE,
     ;                    H1FRC,HE1FRC,DOPPLE)
              ENDIF
           ELSE 
              IF (NBUP-NBLO.EQ.1) THEN
                 HLINOP = HPROFL(NBLO,NBUP,WAVE,WAVEH,T,XNE,H1FRC,
     ;                           HE1FRC,DOPPLE)
                 HLINOP = HLINOP * DEXP((REDCUT-WAVE)*1.D-2)
              ELSE
                 HLINOP = 0.
              ENDIF
           ENDIF
        ELSE
C
C  Blue wing
C
           HLINOP = 0.0
           IF (CUTNEXT) THEN 
              BLUECUT = 1.D8/(EHYD(NBUP+1)-EHYD(NBLO))
           ELSE
              BLUECUT = 1.D8/CONTH(NBLO)
           ENDIF
           IF (IFVAC.EQ.0) BLUECUT = VACAIR(BLUECUT)
           IF (WAVE.LT.BLUECUT) RETURN
           HLINOP = HPROFL(NBLO,NBUP,WAVE,WAVEH,T,XNE,H1FRC,HE1FRC,
     ;                    DOPPLE)
           IF (REFLECT) THEN
               WREFLECT = BLUECUT-(WAVE-BLUECUT)  
               IF (WREFLECT.GT.0.d0) HLINOP = HLINOP + 
     ;        HPROFL(NBLO,NBUP,WREFLECT,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPLE)
           ENDIF
        END IF
C
      RETURN
      END

C***********************************************************************
      REAL FUNCTION HPROFL(N,M,WAVE,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPH)
C  
C hprofl = normalised line profile 
C          normalised to unity wrt integration over frequency
C n      = lower level principal quantum number
C m      = upper level principal quantum number
C wave   = interrogation wavelength in A           real*8
C waveh  = hydrogen line center wavelength in A    real*8
C t      = kinetic temperature in K
C xne    = electron number density in cm-3
C h1frc  = number density of H I in cm-3
C he1frc = number density of He I in cm-3
C dopph  = reduced Doppler width delta_lambda / lambda_0
C        = reduced Doppler width delta_nu / nu_0
C
C  Based on code by Deane Peterson and Bob Kurucz
C
      REAL*8 DELW,WAVE,WAVEH,DOP,D,FREQ,FREQNM,RAYLCUT,WAVE4000
      REAL*8 FINEST(14),FINSWT(14)
      REAL*8 XN2,F,FO,HTOTAL,FREQSQ
      DIMENSION STCOMP(5,4),STALPH(34),
     1          ISTAL(4),LNGHAL(4),STWTAL(34),
     2          STCPWT(5,4),LNCOMP(4)
      DIMENSION ASUM(100),ASUMLYMAN(100),XKNMTB(4,3)
      DIMENSION SIGMA(3),ALF(3),BALFAC(3)
      REAL*8 LYMANH2PLUS(111),LYMANH2(91)
      REAL*8 CLIGHT
      LOGICAL LYMANALF
      SAVE
C
C  Einstein A-value sums for H lines
C
      DATA ASUM/
     1 0.000E+00, 4.696E+08, 9.980E+07, 3.017E+07, 1.155E+07, 5.189E+06,
     2 2.616E+06, 1.437E+06, 8.444E+05, 5.234E+05, 3.389E+05, 2.275E+05,
     3 1.575E+05, 1.120E+05, 8.142E+04, 6.040E+04, 4.560E+04, 3.496E+04,
     4 2.719E+04, 2.141E+04, 1.711E+04, 1.377E+04, 1.119E+04, 9.166E+03,
     5 7.572E+03, 6.341E+03, 5.338E+03, 4.523E+03, 3.854E+03, 3.302E+03,
     6 2.844E+03, 2.460E+03, 2.138E+03, 1.866E+03, 1.635E+03, 1.438E+03,
     7 1.269E+03, 1.124E+03, 9.983E+02, 8.894E+02, 7.947E+02, 7.120E+02,
     8 6.396E+02, 5.759E+02, 5.198E+02, 4.703E+02, 4.263E+02, 3.873E+02,
     9 3.526E+02, 3.215E+02, 2.938E+02, 2.689E+02, 2.465E+02, 2.264E+02,
     A 2.082E+02, 1.918E+02, 1.769E+02, 1.634E+02, 1.512E+02, 1.400E+02,
     1 1.298E+02, 1.206E+02, 1.121E+02, 1.043E+02, 9.720E+01, 9.066E+01,
     2 8.465E+01, 7.912E+01, 7.403E+01, 6.933E+01, 6.498E+01, 6.097E+01,
     3 5.725E+01, 5.381E+01, 5.061E+01, 4.765E+01, 4.489E+01, 4.232E+01,
     4 3.994E+01, 3.771E+01, 3.563E+01, 3.369E+01, 3.188E+01, 3.019E+01,
     5 2.860E+01, 2.712E+01, 2.572E+01, 2.442E+01, 2.319E+01, 2.204E+01,
     6 2.096E+01, 1.994E+01, 1.898E+01, 1.808E+01, 1.722E+01, 1.642E+01,
     7 1.566E+01, 1.495E+01, 1.427E+01, 1.363E+01/
C
C  For Lyman lines only the s-p transition is allowed.
C
      DATA ASUMLYMAN/
     1 0.000E+00, 6.265E+08, 1.897E+08, 8.126E+07, 4.203E+07, 2.450E+07,
     2 1.236E+07, 8.249E+06, 5.782E+06, 4.208E+06, 3.158E+06, 2.430E+06,
     3 1.910E+06, 1.567E+06, 1.274E+06, 1.050E+06, 8.752E+05, 7.373E+05,
     4 6.269E+05, 5.375E+05, 4.643E+05, 4.038E+05, 3.534E+05, 3.111E+05,
     5 2.752E+05, 2.447E+05, 2.185E+05, 1.959E+05, 1.763E+05, 1.593E+05,
     6 1.443E+05, 1.312E+05, 1.197E+05, 1.094E+05, 1.003E+05, 9.216E+04,
     7 8.489E+04, 7.836E+04, 7.249E+04, 6.719E+04, 6.239E+04, 5.804E+04,
     8 5.408E+04, 5.048E+04, 4.719E+04, 4.418E+04, 4.142E+04, 3.888E+04,
     9 3.655E+04, 3.440E+04, 3.242E+04, 3.058E+04, 2.888E+04, 2.731E+04,
     A 2.585E+04, 2.449E+04, 2.322E+04, 2.204E+04, 2.094E+04, 1.991E+04,
     1 1.894E+04, 1.804E+04, 1.720E+04, 1.640E+04, 1.566E+04, 1.496E+04,
     2 1.430E+04, 1.368E+04, 1.309E+04, 1.254E+04, 1.201E+04, 1.152E+04,
     3 1.105E+04, 1.061E+04, 1.019E+04, 9.796E+03, 9.419E+03, 9.061E+03,
     4 8.721E+03, 8.398E+03, 8.091E+03, 7.799E+03, 7.520E+03, 7.255E+03,
     5 7.002E+03, 6.760E+03, 6.530E+03, 6.310E+03, 6.100E+03, 5.898E+03,
     6 5.706E+03, 5.522E+03, 5.346E+03, 5.177E+03, 5.015E+03, 4.860E+03,
     7 4.711E+03, 4.569E+03, 4.432E+03, 4.300E+03/
C
      DATA N1/0/, M1/0/
C
C  Fine structure components for alpha lines in FREQ*10**-7
C
      DATA STALPH/ -730.,  370.,  188.,  515.,  327.,  619., -772.,
     1             -473., -369.,  120.,  256.,  162.,  285., -161.,
     2             -38.3,  6.82, -174., -147., -101., -77.5,   55.,
     3              126.,   75.,  139.,  -60.,   3.7,   27.,  -69., 
     4              -42.,  -18.,  -5.5,  -9.1,  -33.,  -24./
C
C  Alpha component weights
C
      DATA STWTAL/1.,2.,1.,2.,1.,2.,1.,2.,3.,1.,2.,1.,2.,1.,4.,6.,1.,
     1            2.,3.,4.,1.,2.,1.,2.,1.,4.,6.,1.,7.,6.,4.,4.,4.,5./
      DATA ISTAL /1, 3, 10, 21/
      DATA LNGHAL/2, 7, 11, 14/
C
C  Fine structure for M.EQ.INFINITY IN FREQ*10**-7
C
      DATA STCOMP/   0.,    0.,    0.,    0.,    0.,
     2             468.,  576., -522.,    0.,    0.,
     3             260.,  290.,  -33., -140.,    0.,
     4             140.,  150.,   18.,  -27.,   -51./
C
C  Weights for fine structure components
C
      DATA STCPWT/1., 0., 0., 0., 0.,
     2            1., 1., 2., 0., 0.,
     3            1., 1., 4., 3., 0.,
     4            1., 1., 4., 6., 4./
      DATA LNCOMP/1, 3, 4, 5/
C
      DATA XKNMTB/0.0001716, 0.0090190, 0.1001000, 0.5820000,
     1            0.0005235, 0.0177200, 0.1710000, 0.8660000,
     2            0.0008912, 0.0250700, 0.2230000, 1.0200000/
C
C  Lyman alpha red wing profiles following Allard et al (1998, A&A 335,
C  1124). Taken from ATLAS12 code (see Castelli & Kurucz 2001, A&A 372,
C  260)
C  
C  Store log10(I(d omega)) for N(H) = 1e14 cm^-3.
C
C     LYMAN ALPHA QUASI H2+ PROFILE
C     DELTA WAVENO =  -15000+100*(N-1) N=1,111   UP TO -4000
C
      DATA LYMANH2PLUS/
     1 -15.14, -15.06, -14.97, -14.88, -14.80, -14.71, -14.62, -14.53,
     2 -14.44, -14.36, -14.27, -14.18, -14.09, -14.01, -13.92, -13.83,
     3 -13.74, -13.65, -13.57, -13.48, -13.39, -13.30, -13.21, -13.13,
     4 -13.04, -12.95, -12.86, -12.77, -12.69, -12.60, -12.51, -12.40,
     5 -12.29, -12.15, -12.02, -11.90, -11.76, -11.63, -11.53, -11.41,
     6 -11.30, -11.22, -11.15, -11.09, -11.07, -11.06, -11.07, -11.09,
     7 -11.12, -11.16, -11.19, -11.21, -11.24, -11.27, -11.30, -11.33,
     8 -11.36, -11.39, -11.42, -11.45, -11.48, -11.48, -11.48, -11.48,
     9 -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48, -11.48,
     A -11.48, -11.48, -11.48, -11.48, -11.41, -11.40, -11.39, -11.38,
     1 -11.37, -11.36, -11.35, -11.34, -11.33, -11.32, -11.30, -11.29,
     2 -11.28, -11.27, -11.27, -11.27, -11.26, -11.25, -11.24, -11.23,
     3 -11.22, -11.21, -11.20, -11.19, -11.18, -11.17, -11.15, -11.14,
     4 -11.13, -11.12, -11.11, -11.10, -11.09, -11.08, -11.07/
C
C     LYMAN ALPHA QUASI H2 PROFILE
C     DELTA WAVENO = -22000+200*(N-1)  N=1,91  -4000
C
      DATA LYMANH2/
     1 -13.43, -13.32, -13.21, -13.10, -12.98, -12.86, -12.79, -12.72,
     2 -12.65, -12.58, -12.51, -12.47, -12.45, -12.45, -12.48, -12.51,
     3 -12.53, -12.56, -12.59, -12.62, -12.65, -12.69, -12.73, -12.77,
     4 -12.81, -12.85, -12.87, -12.89, -12.90, -12.90, -12.90, -12.90,
     5 -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     6 -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     7 -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     8 -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.90,
     9 -12.90, -12.90, -12.90, -12.90, -12.90, -12.90, -12.89, -12.88,
     A -12.87, -12.86, -12.85, -12.84, -12.83, -12.81, -12.80, -12.79,
     1 -12.78, -12.76, -12.74, -12.72, -12.70, -12.68, -12.65, -12.62,
     2 -12.59, -12.56, -12.53/
C
      PARAMETER (PI = 3.14159265359, SQRTPI = 1.77245385)
      PARAMETER (CLIGHT = 2.9979258E18)
      PARAMETER (CLIGHTCM = 2.99792458E10)
C
C  Most model atmosphere codes include Rayleigh scattering by H atoms 
C  elsewhere, eg. quantum mechanical calculations. This parameter cuts
C  the Lyman alpha natural absorption at this chosen point.  
C
      PARAMETER (RAYLCUT = 1240.D0) ! in Angstroms
c      PARAMETER (RAYLCUT = 3000.D0) ! in Angstroms
C
C  Data for self-broadening from calculations of Barklem, Piskunov and 
C  O'Mara (2000, A&A 363, 1091).
C  SIGMA is in m^2.
C  BALFAC = (4/Pi)**(alf/2) Gamma(2-2/alf) the leading factor for
C  computing the width.
C
      DATA SIGMA /3.304E-18, 6.497E-18, 1.178E-17/
      DATA ALF   /    0.677,     0.455,     0.380/
      DATA BALFAC/  0.97875,   0.97659,   0.97794/
C
C  Precompute variables depending on current physical conditions
C
      T43 = (T/10000.)**0.3
      T3NHE = T43*HE1FRC
      FO = 1.25E-9*XNE**0.66667 ! Holtsmark normal field strength
C
C  Convert wavelengths to frequencies and compute detunings
C
      DELW = WAVE-WAVEH
      FREQ = CLIGHT/WAVE
      FREQNM = CLIGHT/WAVEH
      DEL = FREQ-FREQNM
C
C  Variables dependent on line - compute first time only
C
      IF((N.NE.N1).OR.(M.NE.M1)) THEN
         N1 = N
         M1 = M
         MMN = M-N
         XN = N
         XN2 = XN*XN
         XM = M
         XM2 = XM*XM
         XMN2 = XM2*XN2
         XM2MN2 = XM2-XN2
         GNM = XM2MN2/XMN2
         IF ((MMN.LE.3).AND.(N.LE.4)) THEN
            XKNM = XKNMTB(N,MMN)
         ELSE
            XKNM = 5.5E-5/GNM*XMN2/(1.+.13/FLOAT(MMN))
         END IF
C
C  Lyman alpha wings require special treatment
C
         LYMANALF = .FALSE.
         IF ((N.EQ.1).AND.(M.EQ.2)) LYMANALF = .TRUE.
C
C  Natural damping taken from tables 
C  dnu/ nu_0 = 1/2 * 1/nu_0 * 1/2 pi * Gamma
C
         IF (N.EQ.1) THEN
            RADAMP = ASUMLYMAN(M)
         ELSE
            RADAMP = ASUM(N)+ASUM(M)
         ENDIF
         RADAMP = RADAMP/FREQNM/(4*PI)
C
C  Resonance broadening following Ali & Griem (1966, Phys Rev 144, 366).
C  RESONT is dnu/nu per unit H density. Only the lower state is included 
C  (p-d approx for Balmer lines).  For N > 3 p-states play a very small 
C  role, and van der Waals might even dominate, there is however no 
C  available theory for this.  The lower state resonance broadening is
C  used as a guess.
C
         IF (N.NE.1) THEN 
            RESONT = HFNM(1,N)/(1.-1./ XN2)
         ELSE
            RESONT = HFNM(1,M)/(1.-1./ XM2)
         ENDIF
         RESONT = RESONT * 2.07E-24/GNM
         VDW = 4.45E-26/GNM*(XM2*(7.*XM2+5.))**0.4
         STARK = 1.6678E-18*FREQNM*XKNM
C
C  Approximately include fine structure.  Exact pattern for alpha lines, 
C  M infinite pattern used for others.
C
         IF (N.GT.4 .OR. M.GT.10) THEN
            IFINS = 1
            FINEST(1) = 0.
            FINSWT(1) = 1.
         ELSE IF (MMN.GT.1) THEN
            IFINS = LNCOMP(N)
            DO 1 I = 1,IFINS
            FINEST(I) = STCOMP(I, N)*1.D7
            FINSWT(I) = STCPWT(I, N)/XN2
   1        CONTINUE
         ELSE
C
C  eg: Ly alpha IFINS=2, IPOS=1, FINEST=-7.3E9,3.7E9, FINSWT=1/3, 2/3
C
            IFINS = LNGHAL(N)
            IPOS = ISTAL(N)
            DO 2 I = 1,IFINS
            K = IPOS-1+I
            FINEST(I) = STALPH(K)*1.D7
            FINSWT(I) = STWTAL(K)/XN2/3.D0
   2        CONTINUE
         END IF
      END IF
C
C  Now compute the profile for the given physical conditions.
C
C  Firstly half-widths: DOPPH, HWSTK, HWLOR, are dnu/nu. DOP is the 
C  doppler width dnu.
C
      HWSTK = STARK*FO
      DOP = DOPPH*FREQNM
C      
C  For low Balmer lines use Barklem, Piskunov & O'Mara (2000, A&A 363, 
C  1091) instead of Ali & Griem for self-broadening. 1E6 factor is due 
C  to SI -> CGS.
C
      IF ((N.EQ.2) .AND. (M.LE.5)) THEN
         VBARH = SQRT(42008.*T)
         RESONT = BALFAC(M-N) * 1.E4*SIGMA(M-N) 
     *                        * (VBARH/1.E4)**(1.-ALF(M-N))
         RESONT = RESONT/FREQNM*1.E6/2./PI
      ENDIF
C
      HWLOR = RESONT*H1FRC + VDW*T3NHE + RADAMP
C
C  *Approximate* convolution is achieved by two cases:
C   (a) close to the core -> take the dominant mechanism.
C   (b) outside the core -> add all mechanisms.
C
C  Mechanisms are divided into three:
C   (1) Doppler core with fine structure.
C   (2) (Lorentz) damping parts - self, natural and helium.
C   (3) Stark broadening.
C

C
C  Determine the largest half-width 
C    NWID=1, Doppler
C        =2, Lorentz
C        =3, Stark
C
      IF (DOPPH.GE.HWLOR.AND.DOPPH.GE.HWSTK) THEN
        NWID = 1
        HFWID = DOPPH*FREQNM
      ELSE IF(HWLOR.GE.DOPPH.AND.HWLOR.GE.HWSTK) THEN
        NWID = 2
        HFWID = HWLOR*FREQNM
      ELSE
        NWID = 3
        HFWID = HWSTK*FREQNM
      END IF
C
C  HTOTAL is the profile.
C
      HTOTAL=0.D0
C
C  Case 1: Doppler core including fine structure 
C
      IF ((ABS(DEL).GT.HFWID).OR.(NWID.EQ.1)) THEN
          DO 3 I = 1,IFINS
          D = DABS(FREQ-FREQNM-FINEST(I))/DOP
          IF(D.LE.7.) HTOTAL = HTOTAL + EXP(-D*D)*FINSWT(I)/(SQRTPI*DOP)
   3      CONTINUE
      ENDIF
C  
C  Case 2: Self-broadening section - includes natural and vdW due to 
C  helium for convenience.
C
      IF ((ABS(DEL).GT.HFWID).OR.(NWID.EQ.2)) THEN
C         
C  All lines except Lyman alpha:
C
C  The full classical expression for the line shape is used (eqn 4-73 in
C  Aller). This extends the profile to large detuning, but note the 
C  impact approximation is not valid.  
C
C  The classical damping constant is frequency dependent. Assume similar 
C  dependence of quantum mechanical damping (cf. eqn 4-116 Aller).
C
         IF (.NOT.LYMANALF) THEN
            HHW = HWLOR*FREQNM 
            HHW = HHW * FREQ*FREQ/(FREQNM*FREQNM)    
            FREQSQ = FREQ*FREQ
            XSELF = 4.*FREQSQ*HHW/PI/
     ;              ((DEL*(FREQNM+FREQ))**2. + FREQSQ*4.*HHW**2.)
C
C  Lyman alpha:
C
         ELSE
C
C  Natural broadening 
C
C  This is cut at some specified detuning in the red.
C
            IF (WAVE.LT.RAYLCUT) THEN 
               HHW = RADAMP*FREQNM   
               HHW = HHW * FREQ*FREQ/(FREQNM*FREQNM)    
               FREQSQ = FREQ*FREQ
               HTOTAL = HTOTAL + 4.*FREQSQ*HHW/PI/
     ;             ((DEL*(FREQNM+FREQ))**2. + FREQSQ*4.*HHW**2.)
            ENDIF
C
C  Self broadening of Lyman alpha:
C 
C  Red wing of Lyman Alpha before the satellite region and blue wing 
C  (ie. redward of ~1190 A).
C  
C  The impact approximation breaks down quickly here validity is ~ 1 A 
C  at T~5000 K (Lortet & Roueff 1969, A&A 3, 462) but we use Ali & Griem 
C  anyway (~ Barklem et al for 2p).
C
C  A factor 1.13 is needed to match this region to the far red wing 
C  (below) at -4000 cm^-1 (62 A) detuning. A linear transition is used.
C
            IF (FREQ.GT.(82259.105-4000.)*CLIGHTCM) THEN
               XFAC = 1.+0.13*(DABS(DELW)/62.)
               IF (FREQ.GT.FREQNM) XFAC = 1.0
               HHW = (XFAC*RESONT*H1FRC + VDW*T3NHE) * FREQNM
               XSELF = HHW/PI/(DEL*DEL+HHW*HHW)
C
C  Far red wing of Lyman Alpha in satellite region:
C
C  Self-broadening following Kurucz's implementation of Allard et al 
C  (A&A 335, 1124, 1998), see Castelli & Kurucz (A&A 372, 260 2001) for 
C  more details.
C
C  Tables LYMANH2() are on a delta wavenumber grid spacing of 200 
C  waves/cm -22000, -21800, ...., -4000 waves/cm (91 points). The 
C  profiles store log10(I(d omega)) for N(H) = 1e14 cm^-3. Assumed 
C  insensitive to T, and linear scaling with N(H).  Note I(d freq) = 
C  I(d omega)/c 
C
            ELSE IF (FREQ.GT.20000.*CLIGHTCM) THEN 
               SPACING = 200.*CLIGHTCM
               FREQ22000 = (82259.105-22000.)*CLIGHTCM
C
C  If redward of last point (1660 -> 5000 Angstrom) extrapolate,
C  otherwise (1278 -> 1660 Angstrom) linear (in log10 profiles) 
C  interpolation 
C  
               IF (FREQ.LT.FREQ22000) THEN
                  XLYMANH2 = (LYMANH2(2)-LYMANH2(1))/SPACING*
     *                         (FREQ-FREQ22000)+LYMANH2(1)
               ELSE
                  ICUT = (FREQ-FREQ22000)/SPACING
                  GRIDFREQ = ICUT*SPACING+FREQ22000
                  XLYMANH2 = (LYMANH2(ICUT+2)-LYMANH2(ICUT+1))/SPACING*
     *                         (FREQ-GRIDFREQ)+LYMANH2(ICUT+1)
               ENDIF
               XLYMANH2 = (10.**(XLYMANH2-14.))*H1FRC/CLIGHTCM 
               XSELF = XLYMANH2
            ELSE
               XSELF = 0.
            END IF
C
C  End Lyman alpha section
C
         ENDIF
C
         HTOTAL = HTOTAL + XSELF
C
C  End self-broadening section
C
      ENDIF
C
C  Case 3: Stark section
C
      IF ((ABS(DEL).GT.HFWID).OR.(NWID.EQ.3)) THEN
C
C  Stark wings due to ions and electrons.  If Lyman alpha this is for 
C  before (blueward of) the satellite region.
C
          IF ((.NOT.(LYMANALF)).OR.
     ;             (FREQ.GT.(82259.105-4000.)*CLIGHTCM)) THEN        
              XSTARK = STARK1(N,M,WAVE,WAVEH,T,XNE)
C
C  If Lyman alpha we match the static ion part to the Allard et al 
C  value at -4000 cm^-1 detuning on red wing.
C
              IF ((LYMANALF).AND.(FREQ.LT.(82259.105)*CLIGHTCM)) THEN
                WAVE4000 = 1.D8/(82259.105-4000.)
                XHOLT4000 = 0.5 * STARK1(N,M,WAVE4000,WAVEH,T,XNE)
                IF (XHOLT4000.GT.0.D0) THEN
                  SAT4000 = (10.**(-11.07-14.))/CLIGHTCM*XNE
                  XFAC = 1.+(SAT4000/XHOLT4000-1.)*DABS(DELW)/62.
                  XSTARK = XSTARK * 0.5 * (1. + XFAC)
                ENDIF
              ENDIF
C
C  Stark wings due to protons in the satellite region:
C
C  Following Kurucz's implementation of Allard et al (A&A 335, 1124, 
C  1998), see Castelli & Kurucz (A&A 372, 260 2001) for more details.
C
C  Tables LYMANH2PLUS() are on a delta wavenumber grid spacing of 100 
C  waves/cm -15000, -14900, ...., -4000 waves/cm (111 points). 
C
C  The profiles store log10(I(d omega)) for N(H+) = 1e14 cm^-3. Assumed 
C  insensitive to T, and linear scaling with N(H+). Note I(d freq) = 
C  I(d omega)/c and we assume N(H+) = N(e-).
C
         ELSE IF (FREQ.GT.20000.*CLIGHTCM) THEN 
            SPACING=100.*CLIGHTCM
            FREQ15000=(82259.105-15000.)*CLIGHTCM
C
C  If redward of last point (1487 -> 5000 Angstroem) extrapolate,
C  otherwise (1278 -> 1487 Angstroem) interpolation. 
C  
            IF (FREQ.LT.FREQ15000) THEN
               XLYMANH2PLUS = (LYMANH2PLUS(2)-LYMANH2PLUS(1))/SPACING*
     ;                      (FREQ-FREQ15000)+LYMANH2PLUS(1)
            ELSE
               ICUT = (FREQ-FREQ15000)/SPACING
               GRIDFREQ = ICUT*SPACING+FREQ15000
               XLYMANH2PLUS = (LYMANH2PLUS(ICUT+2)-LYMANH2PLUS(ICUT+1))
     ;                      /SPACING*(FREQ-GRIDFREQ)+LYMANH2PLUS(ICUT+1)
            ENDIF
            XLYMANH2PLUS = (10.**(XLYMANH2PLUS-14.))*XNE/CLIGHTCM 

C
C  There may be a contribution due to electrons also at long range. 
C  The static Holtsmark profile (1/2 the total STARK1) is used noting 
C  that the correction for quantum effects suggested by Stehle (1994, 
C  A&AS 104, 509 eqn 7) has been applied already.
C
            XSTARK = XLYMANH2PLUS + 0.5 * STARK1(N,M,WAVE,WAVEH,T,XNE)  
         ELSE
            XSTARK = 0.
         END IF
         HTOTAL = HTOTAL + XSTARK
C
C  End Stark section
C
      END IF
C
      HPROFL = HTOTAL
C
      RETURN
      END

C***********************************************************************
      FUNCTION STARK1(N,M,WAVE,WAVEH,T,XNE)
C
C  Returns the Stark broadened line profile.  The code follows Griem's 
C  theories (mostly 1960 ApJ, 132, 883) with corrections to approximate 
C  the Vidal, Cooper & Smith (1973, ApJS 25, 37) profiles.
C
C  Area normalised to unity with frequency.
C
C  by Deane Peterson & Bob Kurucz.
C  (adapted, corrected and comments added by PB)
C
      REAL*8 WAVE,WAVEH,DELW,DEL,F,FO,CLIGHT,FREQ,FREQNM
      REAL*4 K
      DIMENSION Y1WTM(2,2),XKNMTB(4,3)
      LOGICAL LYMANALF
      SAVE
C
C  Knm constants as defined by Griem (1960, ApJ 132, 883) for the long 
C  range Holtsmark profile (due to ions only). Lyman and Balmer series 
C  are from VCS, higher series from elsewhere.
C
      DATA XKNMTB/0.0001716, 0.0090190, 0.1001000, 0.5820000,
     1            0.0005235, 0.0177200, 0.1710000, 0.8660000,
     2            0.0008912, 0.0250700, 0.2230000, 1.0200000/
C
      DATA Y1WTM/1.E18, 1.E17, 1.E16, 1.E14/
      DATA N1/0/, M1/0/
C
      PARAMETER (CLIGHT = 2.9979258E18)
      PARAMETER (PI = 3.14159265359, SQRTPI = 1.77245385)
      PARAMETER (H = 6.62618E-27)  !Planck in cgs
      PARAMETER (K = 1.38066E-16)  !Boltzmann in cgs
C
C  Variables depending on conditions
C
      T4 = T/10000.
      T43 = T4**0.3
      T3NHE = T43*HE1FRC
      XNE16 = XNE**0.1666667
      PP = XNE16*0.08989/SQRT(T) ! the shielding parameter 
      FO = XNE16**4*1.25E-9      ! Holtsmark normal field strength
      Y1B = 2./(1.+0.012/T*SQRT(XNE/T))
      Y1S = T43/XNE16
      C1D = FO*78940./ T
      C2D = FO**2/5.96E-23/XNE
      GCON1 = 0.2+0.09*SQRT(T4)/(1.+XNE/1.E13)
      GCON2 = 0.2/(1.+XNE/1.E15)
C
      DELW = WAVE-WAVEH
      FREQNM = CLIGHT/WAVEH
      FREQ = CLIGHT/WAVE
      DEL = FREQ-FREQNM
C
C  Variables dependent on line - compute first time only
C
      IF((N.NE.N1).OR.(M.NE.M1)) THEN  
         N1 = N
         M1 = M
         MMN = M-N
         XN = N
         XN2 = XN*XN
         XM = M
         XM2 = XM*XM
         XMN2 = XM2*XN2
         XM2MN2 = XM2-XN2
         GNM = XM2MN2/XMN2
C
C  Knm constants not tabulated from approximate asymptotic expression 
C  (Griem 1960 eqn 33) where 1./(1.+.13/FLOAT(MMN)) appears to be a 
C  correction factor to match to the tables.
C
         IF ((MMN.LE.3).AND.(N.LE.4)) THEN
            XKNM = XKNMTB(N,MMN)
         ELSE
            XKNM = 5.5E-5/GNM*XMN2/(1.+.13/FLOAT(MMN))
         END IF
C
C  Some weighting factors which relate to y1, which is the velocity at 
C  which the minimum impact parameter (where second order perturbation 
C  theory breaks down) and the Lewis cutoff (the limit of validity of 
C  the impact approximation) are the same.
C
         IF(M.EQ.2) THEN
            Y1NUM = 550.
         ELSE IF (M.EQ.3) THEN
            Y1NUM = 380.
         ELSE
            Y1NUM = 320.
         END IF
         IF (MMN.LE.2 .AND. N.LE.2) THEN
            Y1WHT = Y1WTM(N,MMN)
         ELSE IF (MMN.LE.3) THEN
            Y1WHT = 1.E14
         ELSE
            Y1WHT = 1.E13
         END IF
C
         C1CON = XKNM/WAVEH*GNM*XM2MN2
         C2CON = (XKNM/WAVEH)**2
C
      ENDIF
C
C  Compute line profile
C
C  PRQS is the quasistatic ion contribution
C  FNS  is the quasistatic electron contribution rel to PRQS
C  F    is the impact electron contribution
C
C  First compute the width of the impact electron profile roughly Griem
C  (1967, ApJ 147, 1092) eqn for w.
C
      WTY1 = 1./(1.+XNE/Y1WHT)
      Y1SCAL = Y1NUM*Y1S*WTY1+Y1B*(1.-WTY1)
      C1 = C1D*C1CON*Y1SCAL
      C2 = C2D*C2CON
      G1 = 6.77*SQRT(C1)
      BETA = DABS(DELW)/FO/XKNM
      Y1 = C1*BETA
      Y2 = C2*BETA**2
      IF ((Y2.LE.1.E-4).AND.(Y1.LE.1.E-5)) THEN
         GAM = G1*AMAX1(0.,0.2114+LOG(SQRT(C2)/C1))*(1.-GCON1-GCON2)
      ELSE
         GAM = G1*(0.5*EXP(-MIN(80.,Y1))+VCSE1F(Y1)-0.5*VCSE1F(Y2))*
     *            (1.-GCON1/(1.+(90.*Y1)**3)-GCON2/(1.+2000.*Y1))
         IF (GAM.LE.1.E-20) GAM = 0.
      END IF
C
C  Compute individual quasistatic and impact profiles.
C
      PRQS = SOFBET(BETA,PP,N,M)
      IF (GAM.GT.0.) THEN
         F = GAM/PI/(GAM*GAM+BETA*BETA)
      ELSE
         F = 0.D0
      ENDIF
C
C  Fraction of electrons which count as quasistatic. A fit to eqn 8 
C  (2nd term) of Griem (1967, ApJ 147, 1092).
C
      P1 = (0.9*Y1)**2
      FNS = (P1+0.03*SQRT(Y1))/(P1+1.)
C
C  DBETA (=dBeta/dfreq) changes the area normalisation. 
C  DSQRT(WAVE/WAVEH) corrects the long range part to dfreq**-5/2
C  asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).
C
      DBETA = CLIGHT/FREQ/FREQ/XKNM/FO
      STARK1 = (PRQS*(1.+FNS)+F)*DBETA * DSQRT(WAVE/WAVEH)
C
C  The red wing is multiplied by the Boltzmann factor to roughly account
C  for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume 
C  absorption case.  If emission do for DEL.GT.0.
C
      IF (DEL.LT.0.d0) STARK1 = STARK1 * DEXP(-DABS(H*DEL)/K/T)
C
      END


C***********************************************************************
      FUNCTION HFNM(N,M)
C
C  HFNM calculates hydrogen oscillator strengths
C
C  From Kurucz codes.
C
      SAVE
      DATA NSTR/0/,MSTR/0/
C
      HFNM=0.
      IF(M.LE.N) RETURN
      IF(N.NE.NSTR) THEN
        XN=N
        GINF=0.2027/XN**0.71
        GCA=0.124/XN
        FKN=XN*1.9603
        WTC=0.45-2.4/XN**3*(XN-1.)
        NSTR=N
        XM=M
        XMN=M-N
        FK=FKN*(XM/(XMN*(XM+XN)))**3
        XMN12=XMN**1.2
        WT=(XMN12-1.)/(XMN12+WTC)
        FNM=FK*(1.-WT*GINF-(0.222+GCA/XM)*(1.-WT))
        MSTR=M
      ELSE IF(M.NE.MSTR) THEN
        XM=M
        XMN=M-N
        FK=FKN*(XM/(XMN*(XM+XN)))**3
        XMN12=XMN**1.2
        WT=(XMN12-1.)/(XMN12+WTC)
        FNM=FK*(1.-WT*GINF-(0.222+GCA/XM)*(1.-WT))
        MSTR=M
      END IF
      HFNM=FNM
C
      RETURN
      END



C***********************************************************************
      FUNCTION VCSE1F(X)
C
C  E1 function calculator for VCS approximation. It's rough, but 
C  arranged to be fast. X must be >=0.
C
C  From Kurucz codes.
C
      VCSE1F=0.0
      IF(X.LE.0.0) RETURN
      IF(X.LE.0.01) THEN
        VCSE1F=-LOG(X)-0.577215+X
      ELSE IF(X.LE.1.0) THEN
        VCSE1F=-LOG(X)-0.57721566+X*(0.99999193+X*(-0.24991055+
     +                            X*(0.05519968+X*(-0.00976004+
     +                            X*0.00107857))))
      ELSE IF(X.LE.30.) THEN
        VCSE1F=(X*(X+2.334733)+0.25062)/(X*(X+3.330657)+
     +         1.681534)/X*EXP(-X)
      END IF
C
      RETURN
      END


C***********************************************************************
      FUNCTION SOFBET(B,P,N,M)
C
C  Calculates S(BETA,P) for Hydrogen lines ie. the Holtsmark profile for
C  quasistatic charged particles.  The alpha and beta lines of the first
C  three series are explicitly included. All other cases use the H18 
C  profile. Profiles are normalised to full oscillator strength. Method 
C  is based on Griem (1960, ApJ 132, 883).
C
C  By Deane Peterson and Bob Kurucz.
C
C  STORAGE FOR CORRECTIONS (P,BETA,IND),(P,IND),(P,IND)
C
      DIMENSION PROPBM(5,15,7),C(5,7),D(5,7)
      DIMENSION PP(5),BETA(15)
      DIMENSION PROB1(75),PROB2(75),PROB3(75),PROB4(75),PROB5(75)
      DIMENSION PROB6(75),PROB7(75)
      DIMENSION C1(5),C2(5),C3(5),C4(5),C5(5),C6(5),C7(5)
      DIMENSION D1(5),D2(5),D3(5),D4(5),D5(5),D6(5),D7(5)
      EQUIVALENCE (PROPBM(1,1,1),PROB1(1)),(PROPBM(1,1,2),PROB2(1))
      EQUIVALENCE (PROPBM(1,1,3),PROB3(1)),(PROPBM(1,1,4),PROB4(1))
      EQUIVALENCE (PROPBM(1,1,5),PROB5(1)),(PROPBM(1,1,6),PROB6(1))
      EQUIVALENCE (PROPBM(1,1,7),PROB7(1))
      EQUIVALENCE (C(1,1),C1(1)),(C(1,2),C2(1)),(C(1,3),C3(1))
      EQUIVALENCE (C(1,4),C4(1)),(C(1,5),C5(1)),(C(1,6),C6(1))
      EQUIVALENCE (C(1,7),C7(1))
      EQUIVALENCE (D(1,1),D1(1)),(D(1,2),D2(1)),(D(1,3),D3(1))
      EQUIVALENCE (D(1,4),D4(1)),(D(1,5),D5(1)),(D(1,6),D6(1))
      EQUIVALENCE (D(1,7),D7(1))
      SAVE PROPBM,C,D,PP,BETA
C
C  Lyman alpha
C
      DATA PROB1/
     1-.980,-.967,-.948,-.918,-.873,-.968,-.949,-.921,-.879,-.821,
     2-.950,-.922,-.883,-.830,-.764,-.922,-.881,-.830,-.770,-.706,
     3-.877,-.823,-.763,-.706,-.660,-.806,-.741,-.682,-.640,-.625,
     4-.691,-.628,-.588,-.577,-.599,-.511,-.482,-.484,-.514,-.568,
     5-.265,-.318,-.382,-.455,-.531,-.013,-.167,-.292,-.394,-.478,
     6 .166,-.056,-.216,-.332,-.415, .251, .035,-.122,-.237,-.320,
     7 .221, .059,-.068,-.168,-.247, .160, .055,-.037,-.118,-.189,
     8 .110, .043,-.022,-.085,-.147/
      DATA C1 /-18.396, 84.674,-96.273,  3.927, 55.191/
      DATA D1 / 11.801,  9.079, -0.651,-11.071,-26.545/
C
C  Lyman beta
C
      DATA PROB2/
     1-.242, .060, .379, .671, .894, .022, .314, .569, .746, .818,
     2 .273, .473, .605, .651, .607, .432, .484, .489, .442, .343,
     3 .434, .366, .294, .204, .091, .304, .184, .079,-.025,-.135,
     4 .167, .035,-.082,-.189,-.290, .085,-.061,-.183,-.287,-.374,
     5 .032,-.127,-.249,-.344,-.418,-.024,-.167,-.275,-.357,-.420,
     6-.061,-.170,-.257,-.327,-.384,-.047,-.124,-.192,-.252,-.306,
     7-.043,-.092,-.142,-.190,-.238,-.038,-.070,-.107,-.146,-.187,
     8-.030,-.049,-.075,-.106,-.140/
      DATA C2 / 95.740, 18.489, 14.902, 24.466, 42.456/
      DATA D2 / -6.665, -7.136,-10.605,-15.882,-23.632/
C
C  Balmer alpha
C
      DATA PROB3/
     1-.484,-.336,-.206,-.111,-.058,-.364,-.264,-.192,-.154,-.144,
     2-.299,-.268,-.250,-.244,-.246,-.319,-.333,-.337,-.336,-.337,
     3-.397,-.414,-.415,-.413,-.420,-.456,-.455,-.451,-.456,-.478,
     4-.446,-.441,-.446,-.469,-.512,-.358,-.381,-.415,-.463,-.522,
     5-.214,-.288,-.360,-.432,-.503,-.063,-.196,-.304,-.394,-.468,
     6 .063,-.108,-.237,-.334,-.409, .151,-.019,-.148,-.245,-.319,
     7 .149, .016,-.091,-.177,-.246, .115, .023,-.056,-.126,-.189,
     8 .078, .021,-.036,-.091,-.145/
      DATA C3 /-25.088,145.882,-50.165,  7.902, 51.003/
      DATA D3 /  7.872,  5.592, -2.716,-12.180,-25.661/
C
C  Balmer beta
C
      DATA PROB4/
     1-.082, .163, .417, .649, .829, .096, .316, .515, .660, .729,
     2 .242, .393, .505, .556, .534, .320, .373, .394, .369, .290,
     3 .308, .274, .226, .152, .048, .232, .141, .052,-.046,-.154,
     4 .148, .020,-.094,-.200,-.299, .083,-.070,-.195,-.299,-.385,
     5 .031,-.130,-.253,-.348,-.422,-.023,-.167,-.276,-.359,-.423,
     6-.053,-.165,-.254,-.326,-.384,-.038,-.119,-.190,-.251,-.306,
     7-.034,-.088,-.140,-.190,-.239,-.032,-.066,-.103,-.144,-.186,
     8-.027,-.048,-.075,-.106,-.142/
      DATA C4 / 93.783, 10.066,  9.224, 20.685, 40.136/
      DATA D4 / -5.918, -6.501,-10.130,-15.588,-23.570/
C
C  Paschen alpha
C
      DATA PROB5/
     1-.819,-.759,-.689,-.612,-.529,-.770,-.707,-.638,-.567,-.498,
     2-.721,-.659,-.595,-.537,-.488,-.671,-.617,-.566,-.524,-.497,
     3-.622,-.582,-.547,-.523,-.516,-.570,-.545,-.526,-.521,-.537,
     4-.503,-.495,-.496,-.514,-.551,-.397,-.418,-.448,-.492,-.547,
     5-.246,-.315,-.384,-.453,-.522,-.080,-.210,-.316,-.406,-.481,
     6 .068,-.107,-.239,-.340,-.418, .177,-.006,-.143,-.246,-.324,
     7 .184, .035,-.082,-.174,-.249, .146, .042,-.046,-.123,-.190,
     8 .103, .036,-.027,-.088,-.146/
      DATA C5 /-19.819, 94.981,-79.606,  3.159, 52.106/
      DATA D5 / 10.938,  8.028, -1.267,-11.375,-26.047/
C
C  Paschen beta
C
      DATA PROB6/
     1-.073, .169, .415, .636, .809, .102, .311, .499, .639, .710,
     2 .232, .372, .479, .531, .514, .294, .349, .374, .354, .279,
     3 .278, .253, .212, .142, .040, .215, .130, .044,-.051,-.158,
     4 .141, .015,-.097,-.202,-.300, .080,-.072,-.196,-.299,-.385,
     5 .029,-.130,-.252,-.347,-.421,-.022,-.166,-.275,-.359,-.423,
     6-.050,-.164,-.253,-.325,-.384,-.035,-.118,-.189,-.252,-.306,
     7-.032,-.087,-.139,-.190,-.240,-.029,-.064,-.102,-.143,-.185,
     8-.025,-.046,-.074,-.106,-.142/
      DATA C6 /111.107, 11.910,  9.857, 21.371, 41.006/
      DATA D6 / -5.899, -6.381,-10.044,-15.574,-23.644/
C
C  Balmer 18
C
      DATA PROB7/
     1 .005, .128, .260, .389, .504, .004, .109, .220, .318, .389,
     2-.007, .079, .162, .222, .244,-.018, .041, .089, .106, .080,
     3-.026,-.003, .003,-.023,-.086,-.025,-.048,-.087,-.148,-.234,
     4-.008,-.085,-.165,-.251,-.343, .018,-.111,-.223,-.321,-.407,
     5 .032,-.130,-.255,-.354,-.431, .014,-.148,-.269,-.359,-.427,
     6-.005,-.140,-.243,-.323,-.386, .005,-.095,-.178,-.248,-.307,
     7-.002,-.068,-.129,-.187,-.241,-.007,-.049,-.094,-.139,-.186,
     8-.010,-.036,-.067,-.103,-.143/
      DATA C7 /511.318,  1.532,  4.044, 19.266, 41.812/
      DATA D7 / -6.070, -4.528, -8.759,-14.984,-23.956/
      DATA PP/0.,.2,.4,.6,.8/
      DATA BETA/1.,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.310,
     1          7.943,10.,12.59,15.85,19.95,25.12/
C
      IF(B.GT.500.) THEN
C
C  Very large B
C
        B2=B*B
        SOFBET=(1.5/SQRT(B)+27./B2)/B2
        RETURN
      END IF
C
C  Other cases
C
      CORR=1.
      B2=B*B
      SB=SQRT(B)
      INDX=7
      MMN=M-N
      IF(N.LE.3.AND.MMN.LE.2) INDX=2*(N-1)+MMN
C
C  Determine relevant Debye range
C
      IM=MIN(INT(5.*P)+1,4)
      IP=IM+1
      WTPP=5.*(P-PP(IM))
      WTPM=1.-WTPP
      IF(B.LE.25.12) THEN
C
        JP=2
   1    IF(B.GT.BETA(JP).AND.JP.LT.15) THEN
          JP=JP+1
          GO TO 1
        END IF
        JM=JP-1
C
        WTBP=(B-BETA(JM))/(BETA(JP)-BETA(JM))
        WTBM=1.-WTBP
        CBP=PROPBM(IP,JP,INDX)*WTPP+PROPBM(IM,JP,INDX)*WTPM
        CBM=PROPBM(IP,JM,INDX)*WTPP+PROPBM(IM,JM,INDX)*WTPM
        CORR=1.+CBP*WTBP+CBM*WTBM
C
C  Get approximate profile for the inner part
C
        PR1=0.
        PR2=0.
        WT=AMAX1(MIN(0.5*(10.-B),1.),0.)
        IF(B.LE.10.) PR1=8./(83.+(2.+0.95*B2)*B)
        IF(B.GE.8.)  PR2=(1.5/SB+27./B2)/B2
        SOFBET=(PR1*WT+PR2*(1.-WT))*CORR
      ELSE
C
C  Asymptotic part for medium B's
C
        CC=C(IP,INDX)*WTPP+C(IM,INDX)*WTPM
        DD=D(IP,INDX)*WTPP+D(IM,INDX)*WTPM
        CORR=1.+DD/(CC+B*SB)
        SOFBET=(1.5/SB+27./B2)/B2*CORR
      END IF
C
      RETURN
      END


C***********************************************************************
      REAL*8 FUNCTION AIRVAC(W)
      REAL*8 W, WAVEN, WNEW
C
C  W is air wavelength in Angstroms. WAVEN is air wavenumber which is 
C  usually good enough, must iterate for exact solution
C
      WAVEN = 1.D8 / W
      WNEW = W * (1.0000834213D0 + 2406030.D0 / (1.30D10 - WAVEN**2) +
     +       15997.D0 / (3.89D9 - WAVEN**2))
      WAVEN = 1.D8 / WNEW
      WNEW = W * (1.0000834213D0 + 2406030.D0 / (1.30D10 - WAVEN**2) +
     +       15997.D0 / (3.89D9 - WAVEN**2))
      WAVEN = 1.D8 / WNEW
      AIRVAC = W * (1.0000834213D0 + 2406030.D0 / (1.30D10 - WAVEN**2) +
     +         15997.D0 / (3.89D9 - WAVEN**2))
C
      RETURN
      END


C***********************************************************************
      REAL*8 FUNCTION VACAIR(W)
      REAL*8 W, WAVEN
C
C  W is vacuum wavelength in Angstroms
C
      WAVEN = 1.D8 / W
      VACAIR = W / (1.0000834213D0 + 2406030.D0 / (1.30D10 - WAVEN**2) +
     +       15997.D0 / (3.89D9 - WAVEN**2))
C
      RETURN
      END


     

C***********************************************************************
 
      REAL*8 FUNCTION WCALC(NH, NE, NHE, NS, TEMP)
C  	
C  Computes the occupation probability for a H atom in 
C  state with effective principle quantum number NS in a plasma
C  enviroment with NH, NE, NHE number densities of H, ions,
C  and He atoms respectively.  This code assumes the perturbing
C  neutral H and He are in the ground state, (noting the hard
C  sphere approximation used is quite crude anyway) and that ions 
C  are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
C
C  See eqn 4.71 Hummer & Milhalas (1988) 
C  Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
C  which is a fit to detailed calculation including correlations,
C  thus the temperature dependence
C
C  Sizes of neutral atoms adopted are sqrt(<r^2>)
C  
C  Coded by Paul Barklem and Kjell Eriksson, Aug 2003
C 	
      IMPLICIT NONE
      REAL NH, NE, NHE, NS, TEMP
      REAL*8 IONH, A0, E, PI, K, CHI, RIH, NEUTR, ION, NS2, NS4 
      REAL*8 F, A, BETAC, X, WNEUTR, WION, X1, X2      
      PARAMETER (IONH=2.17991E-11) 
      PARAMETER (A0=5.29177E-9)
      PARAMETER (E=4.803207E-10) 
      PARAMETER (PI=3.1415927)
C
      NS2 = NS*NS
      NS4 = NS2*NS2
C  
C  Neutral perturbers
C
      CHI=IONH/(NS2)
      RIH=SQRT(2.5*NS4 + 0.5*NS2)*A0
      X1=RIH + 1.73*A0
      X2=RIH + 1.02*A0
      NEUTR=NH*X1*X1*X1 + NHE*X2*X2*X2
      WNEUTR = DEXP(-4.18879 * NEUTR)  ! 4.18879 is 4*pi/3
C
C  Charged perturbers
C      
      K=1.
      IF(NS .gt. 3.) THEN
C       K=5.33333333*(NS/(NS+1.))**2 * (NS + 7./6.)/(NS2+NS+0.5)
        K=5.33333333*NS/(NS+1.)/(NS+1.)
      ENDIF
c      ION  = NE*16.*(E*E/CHI/DSQRT(K))**3
c      WION = DEXP(-4.*PI/3. * (ION))
      IF ((NE.GT.10.).AND.(TEMP.GT.10.)) THEN  ! just in case!!! 
        A = 0.09 * DEXP(0.16667d0*LOG(NE)) / SQRT(TEMP)
        X = DEXP(3.15*DLOG(1.+A))
        BETAC = 8.3E14 * EXP(-0.66667*LOG(NE)) * K / NS4
        F = 0.1402d0*X*BETAC*BETAC*BETAC /
     /             (1.+0.1285*X*BETAC*DSQRT(BETAC))
        WION = F/(1.d0+F)
      ELSE 
        WION = 1.0d0     
      ENDIF 
C      
      WCALC = WION * WNEUTR 
      RETURN
      END
     

C***********************************************************************

      REAL FUNCTION HBF(WAVE, NLO)
C
C  Returns the bound-free cross section due to H in level NLO
C  at wavelength WAVE (in Angstroms).  
C
C  It gives an extrapolation redward of usual bound free jumps
C  so that we can use the occupation probability formalism
C  for "broadening" the jump due to nearby particles
C
C  Data (Gaunt factors - Karzas & Latter) from MARCS code
C
      INTEGER N(15), NLO
      REAL L1(12), G1(12), L2(12), G2(12), L3(14), G3(14)
      REAL L4(14), G4(14), L5(14), G5(14), L6(14), G6(14)
      REAL L7(15), G7(15), L8(15), G8(15), L9(6), G9(6)
      REAL L10(16), G10(16), L11(6), G11(6), L12(6), G12(6)
      REAL L13(6), G13(6), L14(6), G14(6), L15(6), G15(6)
      REAL G, W
      REAL*8 WAVE
C
      DATA N/ 12, 12, 14, 14, 14, 14, 15, 15, 6, 16, 6, 6, 6, 6, 6/
C
      DATA L1/   182.,   241.,   300.,   356.,   408.,   456.,
     *           538.,   631.,   729.,   821.,   877.,   912./
      DATA G1/  .9939,  .9948,  .9856,  .9721,  .9571,  .9423,
     *          .9157,  .8850,  .8531,  .8246,  .8076,  .7973/
      DATA L2/   215.,   301.,   398.,   503.,   614.,   729.,
     *           965.,  1313.,  1824.,  2525.,  3144.,  3647./
      DATA G2/  1.043,  1.058,  1.063,  1.061,  1.056,  1.049,
     *          1.033,  1.009,  .9755,  .9342,  .9011,  .8761/
      DATA L3/   222.,   316.,   424.,   545.,   677.,   821.,
     *          1132.,  1641.,  2525.,  4103.,  6034.,  6933.,
     *          7529.,  8207./
      DATA G3/  1.053,  1.072,  1.081,  1.083,  1.081,  1.078,
     *          1.068,  1.051,  1.024,  .9854,  .9458,  .9293,
     *          .9189,  .9075/
      DATA L4/   224.,   321.,   434.,   561.,   703.,   859.,
     *          1205.,  1799.,  2918.,  5252.,  8895., 10998.,
     *         12576., 14588./
      DATA G4/  1.056,  1.077,  1.087,  1.091,  1.091,  1.089,
     *          1.082,  1.069,  1.047,  1.014,  .9736,  .9542,
     *          .9408,  .9247/
      DATA L5/   226.,   324.,   438.,   569.,   715.,   877.,
     *          1241.,  1882.,  3144.,  6034., 11397., 15095.,
     *         18235., 22794./
      DATA G5/  1.058,  1.080,  1.090,  1.095,  1.095,  1.094,
     *          1.089,  1.078,  1.060,  1.030,  .9925,  .9719,
     *          .9563,  .9358/
      DATA L6/   226.,   325.,   441.,   573.,   722.,   887.,
     *          1262.,  1931.,  3282.,  6564., 13448., 18916.,
     *         24120., 32915./
      DATA G6/  1.059,  1.081,  1.092,  1.097,  1.098,  1.097,
     *          1.092,  1.083,  1.067,  1.041,  1.006,  .9851,
     *          .9682,  .9436/
      DATA L7/   227.,   326.,   442.,   576.,   726.,   894.,
     *          1275.,  1961.,  3372.,  6933., 15095., 22347.,
     *         29992., 36616., 44694./
      DATA G7/  1.060,  1.082,  1.093,  1.098,  1.099,  1.099,
     *          1.095,  1.086,  1.071,  1.047,  1.015,  .9952,
     *          .9776,  .9641,  .9494/
      DATA L8/   227.,   326.,   443.,   578.,   729.,   897.,
     *          1284.,  1982.,  3433.,  7196., 16398., 25326.,
     *         35615., 45361., 58446./
      DATA G8/  1.060,  1.082,  1.094,  1.099,  1.100,  1.100,
     *          1.096,  1.088,  1.074,  1.052,  1.022,  1.003,
     *          .9853,  .9708,  .9540/
      DATA L9/   227.,   901.,  7383., 40703., 67537., 74126./
      DATA G9/  1.060,  1.101,  1.056,  .9917,  .9632,  .9576/
      DATA L10/  227.,   327.,   445.,   580.,   732.,   903., 
     *          1294.,  2006.,  3507.,  7529., 18235., 29992.,
     *         45588., 63316., 72940., 91175./
      DATA G10/ 1.060,  1.083,  1.095,  1.100,  1.102,  1.101,
     *          1.098,  1.091,  1.078,  1.058,  1.032,  1.014,
     *          .9970,  .9814,  .9737,  .9606/
      DATA L11/   82.,   905.,  7636., 49822., 96995.,111189./
      DATA G11/ .9379,  1.102,  1.060,  1.002,  .9703,  .9632/
      DATA L12/   82.,   905.,  7720., 53950.,112562.,132138./
      DATA G12/ .9379,  1.102,  1.062,  1.005,  .9732,  .9653/ 
      DATA L13/   82.,   906.,  7793., 57343.,130250.,154534./
      DATA G13/ .9380,  1.102,  1.063,  1.008,  .9758,  .9672/
      DATA L14/   82.,   907.,  7846., 60381.,147056.,178775./
      DATA G14/ .9380,  1.103,  1.064,  1.011,  .9780,  .9689/
      DATA L15/   82.,   908.,  7887., 63316.,162813.,207216./
      DATA G15/ .9380,  1.103,  1.064,  1.014,  .9803,  .9703/
C
      W = WAVE
      G = 0.
      IF (NLO.EQ.1) CALL LOGINT(L1,G1,N(1),W,G)
      IF (NLO.EQ.2) CALL LOGINT(L2,G2,N(2),W,G)
      IF (NLO.EQ.3) CALL LOGINT(L3,G3,N(3),W,G)
      IF (NLO.EQ.4) CALL LOGINT(L4,G4,N(4),W,G)
      IF (NLO.EQ.5) CALL LOGINT(L5,G5,N(5),W,G)
      IF (NLO.EQ.6) CALL LOGINT(L6,G6,N(6),W,G)
      IF (NLO.EQ.7) CALL LOGINT(L7,G7,N(7),W,G)
      IF (NLO.EQ.8) CALL LOGINT(L8,G8,N(8),W,G) 
      IF (NLO.EQ.9) CALL LOGINT(L9,G9,N(9),W,G)
      IF (NLO.EQ.10) CALL LOGINT(L10,G10,N(10),W,G)
      IF (NLO.EQ.11) CALL LOGINT(L11,G11,N(11),W,G)
      IF (NLO.EQ.12) CALL LOGINT(L12,G12,N(12),W,G)
      IF (NLO.EQ.13) CALL LOGINT(L13,G13,N(13),W,G)
      IF (NLO.EQ.14) CALL LOGINT(L14,G14,N(14),W,G)
      IF (NLO.EQ.15) CALL LOGINT(L15,G15,N(15),W,G)
C
      HBF = 1.044E-26 * G * WAVE**3.d0 / (NLO**5.d0)
      RETURN
      END


C***********************************************************************
      SUBROUTINE LOGINT(X, Y, N, XI, YI)
C
C  Does log interpolation and extrapolation of arrays X and Y
C  of length N.  Returns value YI at XI.
C
      IMPLICIT NONE
      INTEGER N, I, J
      REAL X(*), Y(*), XI, YI
      REAL XL(N), YL(N), XIL, YIL, XLS(3), YLS(3), DY        
C
      J = 0
      XIL =  LOG10(XI)
      DO 10 I = 1, N
      XL(I) = LOG10(X(I))
      YL(I) = LOG10(Y(I))
      IF ((XL(I).LT.XIL) .AND. J.EQ.0) J = I
 10   CONTINUE
      IF (J.EQ.0) J = 1
      IF (J.GT.N-3) J = N-3
C
C  polynomial 3 point interpolation or extrapolation
C     
      XLS(1) = XL(J)
      XLS(2) = XL(J+1)
      XLS(3) = XL(J+2)
      YLS(1) = YL(J)
      YLS(2) = YL(J+1)
      YLS(3) = YL(J+2)
      CALL POLINT(XLS, YLS, 3, XIL, YIL, DY)
      YI = 10.**YIL
C
      RETURN
      END
   

C***********************************************************************
C  polynomial interpolation routine

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      INTEGER I,N,NMAX,M,NS
      REAL XA(*),YA(*),X,Y,DY,DIF,DIFT,HO,HP,W,DEN
      PARAMETER (NMAX=10) 
      DIMENSION C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
c          IF(DEN.EQ.0.d0)PAUSE
          IF(DEN.EQ.0.d0)STOP
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


C***********************************************************************
