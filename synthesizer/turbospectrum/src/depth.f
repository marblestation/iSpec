      SUBROUTINE DEPTH(IEL)
*
*-----------------------------------------------------------------------
*
* THIS SUBROUTINE COMPUTES THE POPULATION IN A SPECIFIC LEVEL OF THE
* ATOM UNDER STUDY (N, IN PER GRAM STELLAR MATTER). IT ALSO GIVES THE
* DOPPLER BROADENING VELOCITY (XI0), THE DOPPLER WIDTHS (DNUD) AND THE
* DAMPING PARAMETER A OF THE LINE STUDIED. ALL THIS IS DONE FOR ALL
* DEPTHS IN A GIVEN MODEL ATMOSPHERE.
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
* IF CONSTANT DOPPLER BROADENING VELOCTY (DBVCON) IS PUT IN VIA
* COMMON /POP/, DNUD WILL BE INDEPENDENT OF DEPTH.
*
*   XITE IS MICROTURB IN CM/S
*   XL LINE WAVELENGTH IN CM
*   MA ATOMIC MASS IN GRAMS
*   CHI IS FIRST IONIZATION ENERGY IN EV
*   CHI2 IS SECOND IONIZATION ENERGY IN EV
*   CHIE IS THE LINES LOWER EXCITATION ENERGY IN EV
*   G IS ALWAYS=1 WHILE F IS G*F
*   IDAMP MUST NOT BE=1 THEN A BECOMES = 0 --> only Gaussian profile
*   FDAMP IS FUDGE FACTOR TO CHANGE VAN DER WAALS DAMPING
* if FDAMP < 0, it is instead log(gamma(vdW) at 10000K) as read from, e.g.,
* the VALD3 database (BPz, 05/04-2013)
*   GAMRAD IS RADIATION DAMPING in s^-1. gamrad=sum(1/tau_life(levels)),
*            gamrad/(2.pi.c)=deltanu (cm-1); gamrad/(2.pi)=deltanu (Hz)
*   ALOGC6 IS UNSOELDS (DIFFERENCE BETWEEN THE 2 LEVELS)
*   IJON IS THE IONIZATIN STATE OF THE LINE SPECIES (1=NEUTRAL)
*
* ROUTINE MODIFIED FOR CANARY CODE
*
*-----------------------------------------------------------------------
*
      INCLUDE 'spectrum.inc'
*
      character*1 recipe
      REAL N,MU,MA,MH,M
      REAL NTOT,NTT
      REAL MUM,gamvdw10000
      doubleprecision ionpot
      DIMENSION Q(NDP),theta(ndp),xi0(ndp)
*
      COMMON/POP/ N(NDP),A(NDP),DNUD(NDP),STIM(NDP),ANJON(NDP),DBVCON
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),
     &              NDEPTH
      COMMON/ATOM/ XL,MA,CHI,CHI2,chi3,CHIE,G,IDAMP,FDAMP,
     &             GAMRAD,ALOGC6,IJON
      COMMON/CONST/ BOLTZ,MH,H,C,E,M
      COMMON/CQ/ Q1(NDP),Q2(NDP),Q3(NDP),AQ1(3),AQ2(3),TLIM1,
     &           TLIM2,Q1LIM,Q2LIM,AQ3(3),TLIM3,Q3LIM,TQA(3)
      COMMON/CNUMB/ NTOT(NDP),ntt(ndp),fpartition(ndp),
     &        PH(NDP),HEH,phe(ndp),ph2(ndp)
      COMMON/ATMOL/ NAT,NMOL
      COMMON/PLDATA/ DUMA,DUMB,DUMC,IDUM,DUMD,DUME,IPRESS,GAMPRE
      logical oldpart
      common/oldpart/ oldpart
      real sigmacross,velexp
      common/damping/sigmacross,velexp,recipe
*
      save xi0,theta,q,sigma,constb
*
      EX10(X) = EXP(X/0.43429)
      FYRAPI = 16.*ATAN(1.)
      CONSTA = H*C/(BOLTZ*XL)
      if (.not.oldpart) then
        redmas=1./((1./MH)+(1./MA))
        constb=8.*boltz*4.0/fyrapi/redmas
      endif

      if (recipe.eq.'W') then
* fdamp contains gamma van der Waals at 10000K instead of the fudge factor 
* for Unsoeld
        gamvdw10000=fdamp
*        recipe='W'
      endif
*
* to save time: if same species (iel and ijon): 
*  store q(i),ntt(i) and goto after 20 (move up q=q1...)
      DO I=1,NDEPTH
*
        if (oldpart) goto 21

        TEMP = T(I)
        XITE=XI(I)
        THETA(i) = 5040./TEMP
*
* ntt < 0 means element not present in chemical equilibrium. 
* ntt computed here. Otherwise uses ntt from eqmol (as well as partf, chi).
* When ntt computed here, only neutral, first and second ion considered.

        q(i)=fpartition(i)
        if (ntt(i).ge.0.) goto 30

        NTT(i) = NTOT(I)
        CALL INP3(TQA,AQ1,T(I),Q1(I))
        Q1(I)=EXP(Q1(I))
        IF(TEMP.LE.TLIM1) Q1(I)=Q1LIM
*
* Irwin used for atoms and molecules. Some molecules might not
* be included in his data. ZrO and HF are from Sauval.
* See partf.f for details.  BPz 12/10-95
*  beware also that we are not able to search for species
*  identified by Irwin like: 900.00. Our getinfospecies is
*  fooled by 00 in atoms.
*
        call partf(IEL,1,T(I),1,Q1(I),ionpot)
        chi=ionpot
        IF(IEL.GT.NAT) GOTO 20
        CALL INP3(TQA,AQ2,T(I),Q2(I))
        Q2(I)=EXP(Q2(I))
        IF(TEMP.LE.TLIM2) Q2(I)=Q2LIM
        call partf(IEL,2,T(I),1,Q2(I),ionpot)
        chi2=ionpot
        CALL INP3(TQA,AQ3,T(I),Q3(I))
        Q3(I)=EXP(Q3(I))
        IF(TEMP.LE.TLIM3) Q3(I)=Q3LIM
        call partf(IEL,3,T(I),1,Q3(I),ionpot)
        chi3=ionpot
        SAHA = 2.5*ALOG10(TEMP) - ALOG10(PE(I)) - 0.1761
        SAHA1 = SAHA + ALOG10(Q2(I)/Q1(I)) - THETA(i)*CHI
        SAHA2 = SAHA + ALOG10(Q3(I)/Q2(I)) - THETA(i)*CHI2
        SAHA1 = EX10(SAHA1)
        SAHA2 = EX10(SAHA2)
        IF(IJON.NE.1) GOTO 4
        NTT(i) = NTT(i)/(1. + SAHA1*(1. + SAHA2))
        GOTO 20
    4   IF(IJON.NE.2) GOTO 2
        NTT(i) = NTT(i)*SAHA1/(1. + SAHA1*(1. + SAHA2))
        GOTO 20
    2   PRINT 900, IJON
  900   FORMAT(' STOP|    IJON =',I2,'  IN SUBR. DEPTH **********')
        STOP
   20   CONTINUE
        IF(IJON.EQ.1) Q(I) = Q1(I)
        IF(IJON.EQ.2) Q(I) = Q2(I)
   30   continue
        XI0(i) = SQRT(2.*BOLTZ*TEMP/MA + XITE*XITE)
*
   21   continue
*
        DNUD(I) = XI0(i)/XL
        IF(DBVCON.GT.1.E-6) DNUD(I)=DBVCON/(XL*1.E-5)
*
* NTOT(I) = NO OF CURRENT NUCLEI / GRAM *-MATTER
* NTT(i) = NO OF ATOMS IN THIS IONISATION STATE/ GRAM *-MATTER
* N(I) = NO OF ATOMS IN THIS EXCITATION STATE / GRAM *-MATTER
* ANJON(I) =  RELATIVE NO OF ATOMS IN THIS EXCITATION STATE / ALL
*             EXCITATION AND IONISATION STATES
*
        N(I) = G*NTT(i)*EX10(-CHIE*THETA(i))/Q(I)
        ANJON(I) = N(I)/NTOT(I)
        STIM(I) = 1. - EXP(-CONSTA/T(i))

        if (IDAMP.eq.1) then 
          gamma=0.
        else if (recipe.eq.'R') then
* Only radiative damping (if gamrad present in input line list)
* BPz 12/04-2012
          gamma = gamrad
        else if (recipe.eq.'T') then
* A Tsuji like recipe with Gamma6 /(2c) =0.1cm-1 *sqrt(273/T) * Pg/1atm
* BPz 21/07-1998
          gam6 = 0.1 *2. * 3.e10 *sqrt(273./T(i))* Pg(i)/1.01325e6
          gamma = gam6 +gamrad
        else if (recipe.eq.'U') then
* Classical van der Waals broadening 
*  (Theta inst of T explains factor 0.7log5040)
          ALGH = 8.6735 + .7*ALOG10(THETA(i)) + ALOG10(PH(I)) + 
     &          .4*ALOGC6
          GAMH = EX10(ALGH)*1.E08
cccc          GAM6 = GAMH*(1. + HEH*0.41336)*FDAMP
          GAM6 = GAMH*(1. + phe(i)/ph(i)*0.41336 + 
     &                      ph2(i)/ph(i)*0.85      ) * FDAMP
          GAMMA = GAMRAD + GAM6
cc          if (iel.eq.2) then
cc            if (i.eq.1) print*,' depth,n_spec_low,gamrad,gam6,gamma,a'
cc            print*,i,n(i),gamrad,gam6,gamma,a(i)
cc          endif
        else if (recipe.eq.'W') then
* van der Waals broadening with gamma(vdW) at 10000K from line list 
* BPz 05/04-2013
          gam6= 10.**gamvdw10000*(T(i)/10000.)**0.3*
     &          ph(i)/boltz/T(i)
          gam6=gam6*(1. + phe(i)/ph(i)*0.41336 + 
     &                           ph2(i)/ph(i)*0.85           )
          GAMMA = GAMRAD + GAM6
        else if (recipe.eq.'A'.or.recipe.eq.'B'.or.recipe.eq.'C'.or.
     &           recipe.eq.'S'.or.recipe.eq.'D') then
* Anstee&O'Mara 1995, MNRAS 276, 859, Barklem&O'Mara 1997: 
*  (their w=HWHM!), add helium . BPz adds H2, 29/07-1997
* sigmacross (atomic units, multiply by Bohr radius area=a0^2 to get cm2)
          sigma=sigmacross*2.8002e-17
          relvel=sqrt(constb*t(i))
          gam6anstee=(16./fyrapi)**(velexp/2.)*
     &                exp(gammln(2.-velexp/2.))*
     &                relvel*sigma*(relvel/1.e6)**(-velexp)
*         print *,'w/N=',gam6anstee,
*     &        ' HalfWHM/H atom number density in rad/s cm3'
          gam6anstee=gam6anstee*2.*(PH(I)/boltz/t(i))
* polarizability H=6.70e-25; He=2.07e-25 see Allen73 Sect36
* 0.41 is polarizability(He/H)*velocity ratio in the expression for gamma6
          gam6anstee=gam6anstee*(1. + phe(i)/ph(i)*0.41336 + 
     &                           ph2(i)/ph(i)*0.85           )
*         print *,gam6anstee,' FullWHM in rad/s'
          gamma=gamrad+gam6anstee
        else
          print*,'DEPTH. WARNING!!'
          print*,'depth: recipe= ',recipe,' not implemented!'
          stop
        endif
        a(i) = gamma/(fyrapi*dnud(i))
      enddo
*
      oldpart=.true.
      RETURN
      END
* From BE's depth 29/07-1997. Maybe useful sometime?
**********************************************************************************
* Quadratic Stark broadening for all lines! log C4=-13 typical?
*     alogc4=-12.5
*     gam4=ex10(18.92 + 2.*alogc4/3. + alog10(pe(i)) -
*    &          5.*alog10(t(i))/6.)
*     gam6=gam6 + gam4
**********************************************************************************
      FUNCTION gammln(xx) 
* returns ln gamma(xx) The Gamma function, Numerical recipes, 2nd ed, p 206
      REAL gammln,xx 
      INTEGER j 
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6) 
      SAVE cof,stp 
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, 
     *-.5395239384953d-5,2.5066282746310005d0/ 
      x=xx 
      y=x 
      tmp=x+5.5d0 
      tmp=(x+0.5d0)*log(tmp)-tmp 
      ser=1.000000000190015d0 
      do 11 j=1,6 
        y=y+1.d0 
        ser=ser+cof(j)/y 
11    continue 
      gammln=tmp+log(stp*ser/x) 
      return 
      END 
C  (C) Copr. 1986-92 Numerical Recipes Software 1"n26@1kN. 
* see link ~/smallcalc/Numerical.recipes
*
