C
      SUBROUTINE DETABS(J,JP,NTP,IOUTR)
C
* New version including improvements of CIA, now tabulated into jonabs.dat
* as well as new sources: MgIff, SiIff, and MeIff (other metals)
* OH and CH continuous opac from Kurucz are also in jonabs now.
* from Uppsala version of KE+BE
*      BPz 19/09-2007
*     MARCS v07.3
*
* included Borysow H2-H2 and H2-He opacities. logical switch
* borysow allows choice between these and older ones (Linsky + Patch?)
* BPz 6/10/95.
* Updated 90.03.38 together with jonabs.dat /Bengt Gustafsson
* Data taken from Mathiesen, 1984, Publ. Ser. No. 1, Oslo ...
* parametrized in NDP and included to spherical code:
* VIKTR has been changed to VIKTR**(-1) and H2 scattering added.
*      11-may-90 BPz
* New OP and IRON data for FeI and FeII by Dan Kiselman 98-08-26
C        THIS ROUTINE GIVES THE DETAILS OF THE ABSORPTION MECHANISMS.
C        CHANGES IN THE ABSORPTION-COEFFICIENT PROGRAM ARE EXPECTED TO
C        BE CONFINED TO THE TABLES AND TO THIS ROUTINE.
C        DETABS HAS TWO PURPOSES
C        1. JP=0   DETERMINATION OF WAVELENGTH-INDEPENDENT FACTORS (DEP.ON
C                  T, PE AND THE COMPONENT) STORED IN FAKT.
C        2. JP= THE ACTUAL WAVELENGTH NUMBER.
C                  MULTIPLICATION OF AB, COMPUTED IN SUBROUTINE ABSKO,
C                  BY WAVELENGTH-DEPENDENT FACTORS. SUMMATION OF THE TOTAL
C                  ABSORPTION AND SCATTERING COEFFICIENTS ( SUMABS AND
C                  SUMSCA ).
C
C        N O T E .  BEFORE A CALL ON DETABS FOR PURPOSE 1, SUBROUTINE
C        JON MUST HAVE BEEN CALLED.
C
C        IF J IS LESS THAN OR EQUAL TO ZERO, THE WEIGHT FOR A ROSSELAND MEAN
C        WILL BE COMPUTED AND STORED IN VIKTR (THE WEIGHT BEING 1/VIKTR).
C        NTP IS THE ARRAY INDEX OF THE T-PE POINT.
C        IF IOUTR IS GREATER THAN ZERO AT A CALL WITH JP GREATER THAN ZERO
C        (PART TWO OF THE ROUTINE), DETAILS OF THE ABSORPTION COEFFICIENTS
C        ARE PRINTED. IF IOUTR IS GREATER THAN ONE, A TABLE HEADING IS ALSO
C        PRINTED.
C
C
C        CONTENTS OF COMMON/CI5/, COMMUNICATING PHYSICAL INFORMATION FROM
C        SUBROUTINE JON.
C             ABUND  ABUNDANCES
C             ANJON  FRACTIONS OF IONIZATION
C             H      QUANTUM NUMBER OF THE HIGHEST EXISTING HYDROGENIC LEVEL
C             PART   PARTITION FUNCTIONS
C             DXI    DECREASE OF IONIZATION ENERGY OF HYDROGEN IN ELECTRON VOLTS
C             F1     N(HI)/N(H)
C             F2     N(HII)/N(H)
C             F3     N(H-)/N(H)
C             F4     N(H2+)/N(H)
C             F5     N(H2)/N(H)
C             XKHM   'DISSOCIATION CONSTANT' OF H-
C             XMH    MASS OF THE HYDROGEN ATOM IN GRAMS
C             XMY    GRAMS OF STELLAR MATTER/GRAMS OF HYDROGEN
C
C        DIMENSIONS NECESSARY
C        ABUND(NEL),ANJON(NEL,MAX(NJ)),ELS(NT),H(5),HREST(NT),PART(NEL,MAX(NJ)),
C        PROV(NKOMP)
C        THE DIMENSIONS ARE LOWER LIMITS. DIMENSIONS IN COMMON /CA5/ ARE
C        COMMENTED ON IN SUBROUTINE ABSKO.
C        NEL IS THE NUMBER OF CHEMICAL ELEMENTS INITIATED IN SUBROUTINE INJON
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUTRAL
C             STAGE, FOR ELEMENT I
C        NKOMP IS THE NUMBER OF COMPONENTS, NOT INCLUDING THOSE ADDED BY
C             ANALYTICAL EXPRESSIONS AFTER STATEMENT NO. 13.
C        NT   IS THE NUMBER OF TEMPERATURES-ELECTRON PRESSURES GIVEN AT THE
C             CALL OF SUBROUTINE ABSKO.
C
C
*
      implicit none
*
      include 'spectrum.inc'
      include 'tsuji.par'
*
      real    ELS(ndp),HREST(ndp),FAKRAY(ndp),h2ray(ndp),rayh2,rayhe
      real    KT,LAMBDA3,heray(ndp)
      real    hn,addf,hmin,expa,teta,rayh,expj,stim,xray
      real    xniv,xm2,xm3,hnh,xray2,umc,teta31,xfakh
      real*8  rhokt
      real    omega,propac,h2pres,phtva(ndp),phel(ndp),hepres,
     &        ph2h(ndp),phhe(ndp),h2hpres,hhepres,presmetion
      integer ioutr,jp,kp,komp,nniv,j,m,ntp,n1,i
* common block declarations:
      real    ABUND(16),ANJON(16,5),H(5),PART(16,5),DXI,F1,F2,F3,F4     /ci5/
      real    F5,XKHM,XMH,XMY                                           /ci5/
      COMMON /CI5/ ABUND,ANJON,H,PART,DXI,F1,F2,F3,F4,F5,XKHM,XMH,XMY   /ci5/
*
      real    AB(mkomp),FAKT(mkomp),PE(ndp),T(ndp),XLA(20),XLA3(20),RO        /ca5/
      real    SUMABS,SUMSCA,VIKTR                                       /ca5/
      integer ISET,NLB                                                  /ca5/
      COMMON /CA5/AB,FAKT,PE,T,XLA,XLA3,RO,SUMABS,SUMSCA,VIKTR,ISET,NLB /ca5/
*
      integer IREAD,IWRIT                                               /utput/
      COMMON /UTPUT/ IREAD,IWRIT                                        /utput/
*
      real    PROV(mkomp)                                                  /carc4/
      integer nprova,nprovs,nprov                                       /carc4/
      COMMON /CARC4/ PROV,nprova,nprovs,nprov                           /carc4/
*
      real    ALES,BLES                                                 /ldopac/
      COMMON /LDOPAC/ ALES,BLES                                         /ldopac/

      real rho,xmytsuji,ejontsuji
      doubleprecision partryck,xmettryck,xiontryck,presneutral,
     &     presion,presion2,presion3
      character*20 nametryck
      common/rhotsu/rho,xmytsuji,ejontsuji
      common/fullequilibrium/partryck(ndp,maxmol),
     &       xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     &       nametryck(maxmol)
      common/orderedpress/presneutral(ndp,100),presion(ndp,100),
     &                    presion2(ndp,100),presion3(ndp,100)
*
      CHARACTER*8 SOURCE,ABNAME
      COMMON /CHAR/ ABNAME(mkomp),SOURCE(mkomp)
      integer nkomp,kompla,kompr,komps
      real abkof
      COMMON/CA2/ABKOF(nabdim),KOMPLA(mkomp*20),KOMPR,KOMPS,NKOMP
      CHARACTER*8 NHMIN,NH2PR,NHEPR,NELS,NHRAY,NH2RAY,nh,NHERAY
      DATA NH2PR/'H2+H2'/,NHEPR/'He+H2'/,NELS/'e-sc'/,NHRAY/'H-sc'/,
     & NH2RAY/'H2sc'/,NHMIN/'H-'/,nh/' H'/,NHERAY/'He-sc'/
      LOGICAL FIRST,borysow
      DATA FIRST/.TRUE./
      data borysow/.true./
C
C SAVE ABSORPTION COMPONENT NAMES THE FIRST TIME DETABS IS CALLED
C
      if (FIRST) then
**********************************************
***        open(100,file='dumpopac.dat')
**********************************************
        ABNAME(1)=NHMIN
        ABNAME(2)=nh
        ABNAME(3)=ABNAME(17)
        DO 2 KOMP=20,NKOMP
    2   ABNAME(KOMP-16)=ABNAME(KOMP)
* this renaming is just to avoid printing all the separate H levels
        NPROVA=NKOMP-16+4
        ABNAME(NPROVA-3)=NH2PR
        ABNAME(NPROVA-2)=NHEPR
        abname(nprova-1)='H2+H'
        abname(nprova)='He+H'
        NPROVS=4
        NPROV=NPROVA+NPROVS
        ABNAME(NPROV-3)=NELS
        ABNAME(NPROV-2)=NHRAY
        ABNAME(NPROV-1)=NH2RAY
        ABNAME(NPROV)=NHERAY
        DO KOMP=1,NPROV
           PRINT*,ABNAME(KOMP)
        ENDDO
        FIRST=.FALSE.
      endif
*
      TETA=5040./T(NTP)
      HN=1./(XMH*XMY)
      IF(JP.GT.0) goto 7
C
C        1. COMPUTATION OF WAVELENGTH-INDEPENDENT QUANTITIES
C
      HNH=F1*HN
      rhokt=rho*1.38066d-16*t(ntp)
      kt=1.38066e-16*t(ntp)
***************************************************************************
*** alternative, strictement equivalente:
***      hnh=presneutral(ntp,1)/(rho*1.38066e-16*t(ntp))
*** xmytsuji est mieux que xmy
***      hn=1./(xmh*xmytsuji)
* ON POURRAIT CHANGER PART EN G0 G1 etc...
***************************************************************************
C        H-
* data in jonabs.dat are 10 times larger than in the past,
* in case you would wonder about this 1.e-18 factor for H-
      fakt(1)=(1.d-18/rhokt)*partryck(ntp,1)
      if (nametryck(1).ne.'H -') then 
        print*,nametryck(1), 'should be H-'
        stop 'PROBLEM in detabs!'
      endif
      fakt(19)=pe(ntp)*presneutral(ntp,1)*(1.d-26/rhokt)
C        HI
      TETA31=31.30364*TETA
      xfakh=(2.0898d-26/rhokt)*presneutral(ntp,1)/part(1,1)
      NNIV=15
      XNIV=15.
      IF(H(1).LT.XNIV)NNIV=INT(H(1))
      DO3 M=1,NNIV
        XM2=M*M
        XM3=XM2*FLOAT(M)
    3 FAKT(M+1)=XFAKH*EXP(-TETA31*(1.-1./XM2))/XM3
      FAKT(NNIV+1)=FAKT(NNIV+1)*AMIN1(H(1)-NNIV,1.)
      IF(NNIV.GE.15) goto 6
    4 N1=NNIV+1
      DO5 M=N1,15
    5 FAKT(M+1)=0.
C
C        FREE-FREE HI ABSORPTION
    6 UMC=2.302585*DXI*TETA
      EXPJ=XFAKH*EXP(-TETA31+UMC)/(2.*TETA31)
      ADDF=EXP(TETA31/((FLOAT(NNIV)+0.5)**2)-UMC)-1.
      IF(H(1).LT.XNIV+0.5)ADDF=0.
      FAKT(18)=EXPJ
      HREST(NTP)=EXPJ*ADDF
C
C        H+H
* note ro == rho (set in jon)
      fakt(20)=(presneutral(ntp,1)*(1.d-25/rhokt))*
     &         (presneutral(ntp,1)*(1.e-25/kt))
C        H2+
      fakt(21)=(presneutral(ntp,1)*(1.d-20/rhokt))*
     &         (presion(ntp,1)*(1.e-20/kt))
C        H2- ff
      fakt(22)=pe(ntp)*partryck(ntp,2)/rhokt
      if (nametryck(2).ne.'H H') then 
        print*,nametryck(2), 'should be H2'
        stop 'PROBLEM in detabs!'
      endif
C        He I
      fakt(23)=presneutral(ntp,2)/(rhokt*part(2,1))
* presneutral/(rho*k*T) = no of He I per gram stellar matter
C        He I ff
      fakt(24)=pe(ntp)*(1.e-20/kt)*presion(ntp,2)*(1.d-20/rhokt)
C        He- (stimulated emission included in table)
      fakt(25)=pe(ntp)*presneutral(ntp,2)*(1.d-26/rhokt)
C        C I
      fakt(26)=presneutral(ntp,6)/(rhokt*part(3,1))
C        C II
      fakt(27)=presion(ntp,6)/(rhokt*part(3,2))
C        C I ff
      fakt(28)=(pe(ntp)/kt)*presion(ntp,6)*(1.d-40/rhokt)
C        C II ff
      fakt(29)=(pe(ntp)/kt)*presion2(ntp,6)*(1.d-40/rhokt)
C        C- (stimulated emission included in table)
      fakt(30)=pe(ntp)*presneutral(ntp,6)*(1.d-27/rhokt)
C        N I
      fakt(31)=presneutral(ntp,7)/(rhokt*part(4,1))
C        N II
      fakt(32)=presion(ntp,7)/(rhokt*part(4,2))
C        N- (stimulated emission included in table)
      fakt(33)=pe(ntp)*presneutral(ntp,7)*(1.d-27/rhokt)
C        O I
      fakt(34)=presneutral(ntp,8)/(rhokt*part(5,1))
C        O II
      fakt(35)=presion(ntp,8)/(rhokt*part(5,2))
C        O- (stimulated emission included in table)
      fakt(36)=pe(ntp)*presneutral(ntp,8)*(1.d-26/rhokt)
C        CO-ff (stimulated emission included in table) BPz 1/10-1996
      fakt(37)=pe(ntp)*partryck(ntp,7)*(1.d-26/rhokt)
      if (nametryck(7).ne.'C O ') then 
        print*,nametryck(7), 'should be CO'
        stop 'PROBLEM in detabs!'
      endif
C        H2O-ff (stimulated emission included in table) BPz 30/09-1996
      fakt(38)=pe(ntp)*partryck(ntp,4)*(1.d-26/rhokt)
      if (nametryck(4).ne.'H O H ') then 
        print*,nametryck(4), 'should be H2O'
        stop 'PROBLEM in detabs!'
      endif
C        Mg I
      fakt(39)=presneutral(ntp,12)/(rhokt*part(8,1))
C        Mg II
      fakt(40)=presion(ntp,12)/(rhokt*part(8,2))
C        Al I
      fakt(41)=presneutral(ntp,13)/(rhokt*part(9,1))
C        Al II
      fakt(42)=presion(ntp,13)/(rhokt*part(9,2))
C        Si I
      fakt(43)=presneutral(ntp,14)/(rhokt*part(10,1))
C        Si II
      fakt(44)=presion(ntp,14)/(rhokt*part(10,2))
C        Ca I
      fakt(45)=presneutral(ntp,20)/(rhokt*part(13,1))
C        Ca II
      fakt(46)=presion(ntp,20)/(rhokt*part(13,2))
C        FE I
      fakt(17)=presneutral(ntp,26)/(rhokt*part(15,1))
C        FE II
      fakt(47)=presion(ntp,26)/(rhokt*part(15,2))
* presion/(rho*k*T) = no of Fe II per gram stellar matter
c        CH
      fakt(48)=partryck(ntp,6)*13.02/(6.023d23*rhokt)
      if (nametryck(6).ne.'C H') then 
        print*,nametryck(6), 'should be CH'
        stop 'PROBLEM in detabs!'
      endif
c        OH
      fakt(49)=partryck(ntp,5)*17.01/(6.023d23*rhokt)
      if (nametryck(5).ne.'O H') then 
        print*,nametryck(5), 'should be OH'
        stop 'PROBLEM in detabs!'
      endif
C        Mg I ff
      fakt(50)=(pe(ntp)/kt)*(presion(ntp,12)*(1.d-40/rhokt))
C        Si I ff
      fakt(51)=(pe(ntp)/kt)*(presion(ntp,14)*(1.d-40/rhokt))
C        Me I ff for the sum of all other metals; Al, Li-B, N-Na, P-U
C                using Peach G. 1970, Mem RAS 73,1 hydrogenic approx.
      presmetion=0.0
      do i=3,5
        if(presion(ntp,i).gt.0.0) presmetion=presmetion+presion(ntp,i)
      enddo
      do i=7,11
        if(presion(ntp,i).gt.0.0) presmetion=presmetion+presion(ntp,i)
      enddo
      if(presion(ntp,13).gt.0.0) presmetion=presmetion+presion(ntp,13)
      do i=15,92
        if(presion(ntp,i).gt.0.0) presmetion=presmetion+presion(ntp,i)
      enddo
      fakt(52)=(pe(ntp)/kt)*(presmetion*(1.d-40/rhokt))
* presmetion is similar to the MgII pressure almost throughout the Sun
**
C        ELECTRON SCATTERING
      ELS(NTP)=4.8206E-9*PE(NTP)/(T(NTP)*RO)
**
* added for cool stars -> collision induced absorption (CIA)
cc      ph2=f5*hn*ro
cc      ph2=ph2*1.38e-16*0.987e-6*273.
cc      phtva(ntp)=ph2*ph2/ro
      phtva(ntp)=(partryck(ntp,2)*273./t(ntp)*0.987e-6)*
     &           (partryck(ntp,2)*273./t(ntp)*0.987e-6)/rho
      if (nametryck(2).ne.'H H') then 
        print*,nametryck(2), 'should be H2'
        stop 'PROBLEM in detabs!'
      endif
cc      PHEL(NTP)=(ABUND(2)/ABUND(1))*RO*HN/1.008
cc      PHEL(NTP)=PHEL(NTP)*1.38E-16*0.987E-6*273.
cc      PHEL(NTP)=PHEL(NTP)*PH2/RO
      phel(ntp)=(presneutral(ntp,2)*273./t(ntp)*0.987e-6)*
     &           (partryck(ntp,2)*273./t(ntp)*0.987e-6)/rho
      if (nametryck(2).ne.'H H') then 
        print*,nametryck(2), 'should be H2'
        stop 'PROBLEM in detabs!'
      endif
      ph2h(ntp)=(presneutral(ntp,1)/t(ntp)*2.6945e-04)*
     &             (partryck(ntp,2)/t(ntp)*2.6945e-04)/rho
      if (nametryck(2).ne.'H H') then 
        print*,nametryck(2), 'should be H2'
        stop 'PROBLEM in detabs!'
      endif
      phhe(ntp)=(presneutral(ntp,1)/t(ntp)*2.6945e-04)*
     &          (presneutral(ntp,2)/t(ntp)*2.6945e-04)/rho
**
C        RAYLEIGH SCATTERING
cc      FAKRAY(NTP)=HNH*2./PART(1,1)
      fakray(ntp)=2.*presneutral(ntp,1)/(rhokt*part(1,1))
cc      H2RAY(NTP)=F5*HN
      h2ray(ntp)=partryck(ntp,2)/rhokt
      if (nametryck(2).ne.'H H') then 
        print*,nametryck(2), 'should be H2'
        stop 'PROBLEM in detabs!'
      endif
      heray(ntp)=presneutral(ntp,2)/rhokt

cc      do komp=1,nkomp
cc        print*,'detabs fakt(',komp,')',fakt(komp)
cc      enddo
*
      RETURN
*
C        N O T E . APART FROM VECTORS HREST AND ELS, NONE OF THE
C        TEMPERATURE- OR PRESSURE-DEPENDENT VARIABLES DEFINED ABOVE CAN
C        GENERALLY BE USED AT THE NEXT VISIT BELOW.
C        ANY SET OF FACTORS WHICH IS WANTED SHOULD BE STORED IN AN ARRAY WITH
C        DIMENSION = NT, LIKE HREST AND ELS, OR IN FAKT, WHERE THE DATA FOR
C        FURTHER USE ARE STORED IN SUBR. ABSKO.
C
C        2. WAVELENGTH-DEPENDENT FACTORS. SUMMATION.
C        CORRECTION FOR STIMULATED EMISSION
    7 continue
      EXPA=EXP(-28556.*TETA/XLA(JP))
   11 STIM=1.-EXPA
      LAMBDA3=(XLA(JP)/911.27)**3
C
C        ABSORPTION
      SUMABS=0.
C        H I
      DO 12 KOMP=2,16
        SUMABS=SUMABS+AB(KOMP)
   12 CONTINUE
      sumabs=sumabs+ab(18)
      SUMABS=(SUMABS+HREST(NTP))*XLA3(JP)
      PROV(2)=SUMABS
C        H-
      HMIN=AB(1)+AB(19)/STIM
      SUMABS=SUMABS+HMIN
      PROV(1)=HMIN
C Fe I bf     NEW: OP data from DK 98-08-23 revised later. Still needs revision 
C   (smoothing resonance and not remove them) BPz 19/09-2007
      SUMABS=SUMABS+AB(17)
      PROV(3)=AB(17)
c stimulated emission already included in table in jonabs.dat for negative ions
c He-
      AB(25)=AB(25)/STIM
c C-
      AB(30)=AB(30)/STIM
c N-
      AB(33)=AB(33)/STIM
c O- 
      AB(36)=AB(36)/STIM
c CO-
      AB(37)=AB(37)/stim
c H2O-
      AB(38)=AB(38)/stim
c similar also for CH, OH:
      AB(48)=AB(48)/stim
      AB(49)=AB(49)/stim
c He I ff
      AB(24)=AB(24)*LAMBDA3
c C I ff
      AB(28)=AB(28)*LAMBDA3
c C II ff
      AB(29)=AB(29)*LAMBDA3
c Mg I ff
      AB(50)=AB(50)*LAMBDA3
c Si I ff
      AB(51)=AB(51)*LAMBDA3
c Me I ff (the remaining metals)
      AB(52)=AB(52)*LAMBDA3
C 
      DO 13 KOMP=20,NKOMP
        SUMABS=SUMABS+AB(KOMP)
        PROV(KOMP-16)=AB(KOMP)
   13 CONTINUE

***********************************
***      if (xla(j).gt.1100.and.xla(j).lt.2000..and.ab(20).gt.0.e0) then
***        write(100,*) xla(jp),T(ntp),ab(20),'Doyle H+H'
***      endif
***********************************
C
C        HERE FURTHER ABSORPTION MECHANISMS, GIVEN IN TABLES OR BY
C        ANALYTICAL EXPRESSIONS, MAY BE ADDED.
*
* from old detabs (C-stars) H2-H2 and H2-He pressure induced opacity
* 16/01-1996: Aleksandra Borysow me signale que stim deja inclu dans CIA.
      OMEGA=1./XLA(JP)*1.E+8
      if (borysow) then
cc        call h2borysowopac(omega,t(ntp),propac)
        call CIAh2h2(omega,t(ntp),propac)
      else
        CALL H2OPAC(OMEGA,T(NTP),PROPAC)
      endif
      H2PRES=PROPAC*PHTVA(NTP)/stim
      PROV(NPROVA-3)=H2PRES
      if (borysow) then
cc        call heborysowopac(omega,t(ntp),propac)
        call CIAh2he(omega,t(ntp),propac)
      else
        CALL HEOPAC(OMEGA,T(NTP),PROPAC)
      endif
      HEPRES=PROPAC*PHEL(NTP)/stim
      PROV(NPROVA-2)=HEPRES
      if(borysow) then
        call CIAh2h(omega,t(ntp),propac)
      else
        stop '1 in detabs'
      endif
      h2hpres=propac*ph2h(ntp)/stim
      prov(nprova-1)=h2hpres
      if(borysow) then
        call CIAhhe(omega,t(ntp),propac)
      else
        stop '2 in detabs'
      endif
      hhepres=propac*phhe(ntp)/stim
      prov(nprova)=hhepres
*
      SUMABS=SUMABS+H2PRES+HEPRES+h2hpres+hhepres
      SUMABS=SUMABS*STIM
C
C        SCATTERING
      XRAY=AMAX1(XLA(JP),1026.)
      XRAY2=1./(XRAY*XRAY)
c Changed by BPz 08/10-2007 to stop where Paul BArklem's Hlinop.f starts.
c He accounts for Rayleigh scattering in the far red wings out to that cut 
c (but it is added into the pure absorption coefficient, as all his H I opacity). 
      if (xla(jp).ge.1400.) then
        RAYH=XRAY2*XRAY2*(5.799E-13+XRAY2*(1.422E-6+XRAY2*2.784))*
     &     FAKRAY(NTP)
      else
        RAYH=0.0
      endif
*****************
* Scattering by HeI in cm**2/atom, from Bues & Wehrse, 1976, A&A 51, 461
*  introduced 5/02-2008 by BPz
*
      RAYHe=0.66520e-24*4.*(500./xla(jp))**4*heray(ntp)
*****************
* add h2 from old detabs
      RAYH2=XRAY2*XRAY2*(8.14E-13+XRAY2*(1.28E-6+XRAY2*1.61))*H2RAY(NTP)
      SUMSCA=ELS(NTP)+RAYH+RAYH2+RAYHe
      PROV(NPROV-3)=ELS(NTP)
      PROV(NPROV-2)=RAYH
      PROV(NPROV-1)=RAYH2
      PROV(NPROV)=RAYHe
***********************************
***      if (xla(j).gt.1100.and.xla(j).lt.2000.) then
***        write(100,*) xla(jp),T(ntp),sumabs,sumsca,'tot abs & scatt'
***        if (RAYH.gt.0.e0) then
***          write(100,*) xla(jp),T(ntp),RAYH,'rayleigh H '
***        endif
***        if (RAYH2.gt.0.e0) then
***          write(100,*) xla(jp),T(ntp),RAYH2,'rayleigh H2 '
***        endif
***      endif
***********************************
C
      IF(J.GT.0) goto 15
C
C        WEIGHT FOR A ROSSELAND MEAN
   14 VIKTR=EXPA/(STIM*STIM*(XLA3(JP)*1E-3)**2)
   15 CONTINUE
C
      IF(IOUTR-1)23,21,20
C
C        **** PRINT-OUT ****
   20 WRITE(IWRIT,200)XLA(JP),(ABNAME(KP),KP=1,NPROV)
  200 FORMAT('0WAVEL.=',F7.0,'    ABS       SCAT  ',36A6)
   21 DO22 KP=1,nprova
   22 PROV(KP)=PROV(KP)/SUMABS*STIM
      DO24 KP=1,NPROVS
   24 PROV(NPROVA+KP)=PROV(NPROVA+KP)/SUMSCA
      WRITE(IWRIT,201) T(NTP),sumabs,sumsca,(PROV(KP),KP=1,nprov)
  201 FORMAT(3H T=,F7.1,5X,2E10.4,36F6.4)
   23 CONTINUE


      RETURN
      END
