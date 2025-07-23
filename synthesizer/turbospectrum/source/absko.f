C
      SUBROUTINE ABSKO(NEWT,NT,TSKAL,PESKAL,ISETA,J,ABSK,SPRID)
C
C        THE ROUTINE ADMINISTERS THE COMPUTATION OF ABSORPTION
C        COEFFICIENTS. IT CALLS THE ROUTINES, GIVING THE PROPER THERMO-
C        DYNAMIC INFORMATION ( J O N ) , THE DETAILS OF THE ABSORPTION
C        MECHANISMS ( D E T A B S ) AND THE FACTORS FOR THE INTERPOLATION
C        IN T  ( T A B S ) . IT CHOOSES (IF NECESSARY READS) THE RIGHT SET
C        OF ABSORPTION-COEFFICIENT DATA (ISETA), STATEMENT NO. 5 AND MAKES
C        THE INTERPOLATION IN T, STATEMENTS NOS. 10-18, AND THE SUMMATION
C        OF A ROSSELAND MEAN, IF INDICATED BY J = 0, STATEMENTS NOS. 25-28
C
C        NEWT SHOULD BE GT 1 THE FIRST TIME THIS ROUTINE IS USED,
C                       EQ 1 WHEN A NEW SET OF T-PE IS USED,
C                       EQ 0 OTHERWISE.
C
C        NT IS THE NUMBER OF T-PE POINTS. THE TEMPERATURES T SHOULD BE EX-
C        PRESSED IN KELVIN, THE ELECTRON PRESSURES PE IN DYNES PER CM2
C        ISETA IS THE WAVELENGTH-SET NUMBER, J THE WAVELENGTH NUMBER IN THAT
C        SET. J EQUAL TO ZERO INDICATES THAT A ROSSELAND MEAN IS WANTED.
C        THIS MEAN IS COMPUTED USING THE WAVELENGTH POINTS OF THE ACTUAL
C        SET (ISETA) AND THE QUADRATURE WEIGHTS GIVEN IN ROSW.
C        IN ABSK AND SPRID THE ABSORPTION AND SCATTERING COEFFICIENTS PER GRAM
C        OF STELLAR MATTER ARE STORED.
C
C        DIMENSIONS NECESSARY
C        AB(NKOMP),ABSK(1),FAKT(NKOMP),FAKTP(IFADIM),NTPO(NTO),PE(NT),PESKAL(1),
C        ROSW(MAX(NL)),SPRID(1),SUMW(NT),T(NT),TSKAL(1),XLA(MAX(NL)),
C        XLA3(MAX(NL))
C        THE DIMENSIONS ARE LOWER LIMITS.
C        DIMENSIONS OF ARRAYS IN COMMONS /CA2/,/CA3/ AND /CFIL/ ARE COMMENTED
C        ON IN SUBROUTINE INABS, THOSE OF ARRAYS IN COMMON /CA4/ IN SUBROUTINE
C        TABS.
C        NKOMP IS THE NUMBER OF 'COMPONENTS'
C        NL(I) IS THE NUMBER OF WAVELENGTHS IN WAVELENGTH SET I
C        NT IS THE NUMBER OF T-PE POINTS SUPPLIED IN TSKAL AND PESKAL
C        NTO IS THE NUMBER OF POINTS IN THOSE SCALES FOR WHICH A DETAILED
C              PRINT-OUT IS WANTED.
      implicit none
      include 'spectrum.inc'
      include 'tsuji.par'
C
      integer newt,iseta,isetp,nt,ifak,kfak,jp,kp,ntp,ioutr,j,komp
      integer ireadp,nabkof,nkompl,j1,j2,iu,position,nop,np,lunit
      real pg,dum,delsum

      logical first
      real TSKAL(NT),PESKAL(NT),ABSK(NT),SPRID(NT)
      real FAKTP(ifadim)
      real SUMW(NT)
      real tioabs(ndp),h2oabs(ndp)
      integer iread, iwrit
      COMMON/UTPUT/IREAD,IWRIT
      real abkof
      integer kompla,kompr,komps,nkomp
      COMMON/CA2/ABKOF(nabdim),KOMPLA(mkomp*20),KOMPR,KOMPS,NKOMP
      integer ilogta,null
      COMMON/CA3/ILOGTA(mkomp),NULL
      real afak
      integer nofak,nplats
      COMMON/CA4/AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM)
      real ab,fakt,pe,t,xla,xla3,ro,sumabs,sumsca,viktr
      integer iset,nlb
      COMMON/CA5/AB(mkomp),FAKT(mkomp),PE(NDP),T(NDP),XLA(20),XLA3(20),
     &           RO,
     &           SUMABS,SUMSCA,VIKTR,ISET,NLB
      integer ireset,islask,ireat
      COMMON/CFIL/IRESET(numbset),ISLASK,IREAT
      integer nto,ntpo
      COMMON/COUTR/NTO,NTPO(10)
      real rosw
      COMMON/CROS/ROSW(20)
      real f1p,f3p,f4p,f5p,hnic,presmo
      COMMON /CARC3/ F1P,F3P,F4P,F5P,HNIC,PRESMO(30)
      real prov
      integer nprova,nprovs,nprov
      COMMON /CARC4/ PROV(mkomp),NPROVA,NPROVS,NPROV
      real ptio,rosav,poxg1
      COMMON /TIO/PTIO(NDP),ROsav(NDP),POXG1(NDP)
      integer ielem,ion
      real tmolim,molh
      COMMON/CI4/ IELEM(16),ION(16,5),TMOLIM,MOLH
      real eh,fe,fh,fhe,fc,fce,fn,fne,fo,foe,fk,fke,fs,fse
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
      real rotest,prh2o
      COMMON /DENSTY/ ROTEST(NDP),PRH2O(NDP)

      doubleprecision partryck,xmettryck,xiontryck
      character*20 nametryck
      common/fullequilibrium/partryck(ndp,maxmol),
     &       xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),
     &       nametryck(maxmol)

      doubleprecision presneutral,presion,presion2,presion3
      common/orderedpress/presneutral(ndp,100),presion(ndp,100),
     &                    presion2(ndp,100),presion3(ndp,100)
      character species*20, comment*100,path*128
      integer i,nline,ioniz,pathlen,NHtbl_unit
      doubleprecision xlambda
      real nh1,ne,nhe1,nhop,opnh
! 150 is the max allowed number of H lines
      doubleprecision hlambda(150)
      real xlo(150),xup(150),gf(150),npop(150),b_departure(150)
      real cont,total,dopple
      integer nlo(150),nup(150)
      character*9 lname(150)
      logical contonly,firsth,lineonly
      data contonly /.true./, firsth /.true./, lineonly /.false./
!

      data first/.true./

      save first
C
C
      if (first) then
        newt=2
        first=.false.

! LTE for hydrogen. Can be changed later.
        b_departure=1.0

      endif
      ISET=ISETA
      IF(NEWT.GT.1)ISETP=-1
      IF(NEWT.EQ.0)GO TO 5
C
C        FACTORS ONLY DEPENDENT ON T-PE
C
      CALL TABS(NT,TSKAL)
      IFAK=1
      KFAK=1
      JP=0
      KP=1
C
C        LOOP OVER THE T-PE POINTS ('THE FIRST NTP-LOOP')
      DO4 NTP=1,NT
      T(NTP)=TSKAL(NTP)
      PE(NTP)=PESKAL(NTP)
C        IS PRINT-OUT WANTED FOR T-PE POINT NO. NTP
      IOUTR=0
      IF(KP.GT.NTO)GO TO 3
      IF(NTP.EQ.NTPO(KP))GO TO 1
      GO TO 3
    1 IOUTR=1
      KP=KP+1
    3 CONTINUE
C
      CALL JON(T(NTP),PE(NTP),1,PG,RO,DUM,IOUTR,ntp)
! save ro for HI absorption calculation
      rosav(ntp)=ro
cc      print*,'absko back from jon, calling detabs'
      CALL DETABS(J,0,NTP,IOUTR)
cc      print*,'absko back from detabs, nkomp, kompr,komps',
cc     &      nkomp,kompr,komps
C
C        WE STORE THE FAKT ARRAY, MADE IN JON-DETABS IN LONGER ARRAYS NAMELY
C                  IN AFAK FOR TEMPERATURE-INDEPENDENT COMPONENTS
C                  IN FAKTP FOR TEMPERATURE-DEPENDENT ONES
      DO2 KOMP=1,KOMPR
      AFAK(KFAK)=FAKT(KOMP)
cc      print*,'absko fakt(',komp,')',fakt(komp)
    2 KFAK=KFAK+1
      DO4 KOMP=KOMPS,NKOMP
      FAKTP(IFAK)=FAKT(KOMP)
cc      print*,'absko fakt(',komp,')',fakt(komp)
      KFAK=KFAK+NOFAK(IFAK)
    4 IFAK=IFAK+1
C        END OF 'THE FIRST NTP-LOOP'
C
C        READING  OF A NEW WAVELENGTH SET IF INDICATED BY ISET
    5 IF(ISET.EQ.ISETP)GO TO 6
      IREADP=IRESET(ISET)
   51 READ(IREADP,END=52)ISETP,NLB,XLA,XLA3,NABKOF,ABKOF,NKOMPL,KOMPLA
      GO TO 5
   52 REWIND IREADP
      GO TO 51
C        ROSSELAND MEAN OR NOT
    6 IF(J.GT.0)GO TO 9
    7 J1=1
      J2=NLB
      DO8 NTP=1,NT
      SUMW(NTP)=0.
    8 ABSK(NTP)=0.
      GO TO 10
    9 J1=J
      J2=J
C
C        INTERPOLATION IN T
C        LOOP OVER ALL THE WAVELENGTHS IN CASE OF ROSSELAND MEAN. THIS
C        LOOP ENDS IN STATEMENT NO. 26
   10 CONTINUE
      DO26 JP=J1,J2
      KFAK=1
      IFAK=1
      KP=1
C
C        LOOP OVER THE T-PE POINTS ('THE SECOND NTP-LOOP')
      DO26 NTP=1,NT
C
C        IS PRINT-OUT WANTED FOR T-PE POINT NO. NTP
      IOUTR=0
      IF(KP.GT.NTO)GO TO 93
      IF(NTP.EQ.NTPO(KP))GO TO 92
      GO TO 93
   92 IOUTR=1
      KP=KP+1
      IF(KP.EQ.2)IOUTR=2
   93 CONTINUE
      IU=JP
C
C        COMPONENTS WITH ABSORPTION COEFFICIENTS INDEPENDENT OF THE
C        TEMPERATURE
C
      DO14 KOMP=1,KOMPR
      IF(KOMPLA(IU).LE.0)GO TO 12
C        THE VECTOR KOMPLA IS DETERMINED IN SUBROUTINE INABS.
C        KOMPLA GREATER THAN ZERO GIVES THE INDEX IN ABKOF, WHERE THE TABLE FOR
C        THIS COMPONENT AND WAVELENGTH BEGINS.
C        KOMPLA LESS THAN OR EQUAL TO ZERO INDICATES THAT THE ACTUAL ABSORPTION
C        COEFFICIENT FOR THIS COMPONENT AND WAVELENGTH IS ZERO, AS FOUND IN SUB-
C        ROUTINE INABS.
C
   11 position=KOMPLA(IU)
      AB(KOMP)=AFAK(KFAK)*ABKOF(position)
      GO TO 13
   12 AB(KOMP)=0.
   13 KFAK=KFAK+1
   14 IU=IU+NLB
C
C        COMPONENTS WITH T-DEPENDENT ABSORPTION COEFFICIENTS
      DO19 KOMP=KOMPS,NKOMP
      NOP=NOFAK(IFAK)
      IF(NOP.EQ.0)GO TO 17
      IF(KOMPLA(IU).LE.0)GO TO 17
   15 position=NPLATS(IFAK)-1+KOMPLA(IU)
C        THE VECTOR NPLATS IS DETERMINED BY SUBROUTINE TABS. IT GIVES THE ARRAY
C        INDEX OF THE TEMPERATURE AT WHICH THE INTERPOLATION IN ABKOF
C        BEGINS. NOFAK, GIVING INFORMATION ON THE T-INTERPOLATION AND
C        POSSIBLY INDICATING THAT AB=0 (NOFAK=0) IS ALSO DETERMINED BY TABS
C
C        INTERPOLATION
      DELSUM=0.
      DO16 NP=1,NOP
      DELSUM=DELSUM+AFAK(KFAK)*ABKOF(position)
      KFAK=KFAK+1
   16 position=position+1
C
C        HAS THE INTERPOLATION BEEN MADE ON THE LOGARITHM
cc      print*,'absko ifak delsum exp?',ifak,delsum,ilogta(komp)
      IF(ILOGTA(KOMP).GT.0)DELSUM=EXP(DELSUM)
cc      print*,'absko ifak delsum ',ifak,komp,delsum
C        MULTIPLICATION BY FACTOR FROM JON-DETABS
cc      print*,'absko fakt        ',ifak,faktp(ifak)
      DELSUM=DELSUM*FAKTP(IFAK)
cc      print*,'absko fakt*delsum ',ifak,faktp(ifak),delsum
      IF(DELSUM.GE.0.)GO TO 162
C      print*,'osabsko.f faktp(',ifak,')= ',faktp(ifak)
C
C        A NEGATIVE INTERPOLATION RESULT
  161 IF(NULL.GT.0)WRITE(IWRIT,200)KOMP,DELSUM,JP,ISET,T(NTP)
  200 FORMAT(4H AB(,I4,11H) NEGATIVE=,E12.4,5X,17HFOR WAVELENGTH NO,I5,
     *5X,6HSET NO,I5,5X,2HT=,F10.4,'  AND THEREFORE PUT =0 ***ABSKO***')
      AB(KOMP)=0.
      GO TO 18
  162 AB(KOMP)=DELSUM
      GO TO 18
   17 AB(KOMP)=0.
      KFAK=KFAK+NOP
   18 IU=IU+NLB
   19 IFAK=IFAK+1
C
C        WE MULTIPLY BY WAVELENGTH-DEPENDENT  FACTORS AND ADD UP. THIS IS
C        DONE IN DETABS.
cc      print*,'absko, calling detabs 2nd time'
      CALL DETABS(J,JP,NTP,IOUTR)
cc      print*,'absko, back from detabs'
!
! Must add HI bf absorption with improved treatment from Barklem&Piskunov
! BPz 03/04-2019
!
      if (firsth) then
!       read file
        lunit=77
        open(lunit,file='DATA/Hlinedata', status='old')
        read(lunit,*) species,ioniz,nline
        read(lunit,*) comment
*        print*,species,ioniz,nline,comment
        if (species(1:6) /= '01.000' .or. ioniz /= 1) then
          print*, 'wrong H line data file!'
          print*,species,ioniz
          stop 'ERROR!'
        endif
c babsma uses air wavelengths, and so does hbop.f, except for Lyman series
c
        if (nline > 150 ) stop 'increase nline dimension in babsma!'
        do i=1,nline
!         oscillator strengths are computed in hbop.f
          read(lunit,*) hlambda(i),nlo(i),nup(i),xlo(i),xup(i),gf(i),
     &                lname(i)
        enddo
        close(lunit)
        firsth=.false.
      endif
*
      ne=pe(ntp)/(t(ntp)*1.38066e-16)
      nh1=sngl(presneutral(ntp,1))/(t(ntp)*1.38066e-16)
      nhe1=sngl(presneutral(ntp,2))/(t(ntp)*1.38066e-16)
! dopple not used for continuum
      dopple=0.0
      xlambda=dble(xla(jp))

! COMPUTE in LTE. CAN BE CHANGED LATER.

      call hbop(xlambda,nline,nlo,nup,hlambda,
     &         nh1,nhe1,ne,t(ntp),dopple,npop,
     &         b_departure,0,1.,1.,1.,1.,
     &         total,cont,contonly,lineonly,.false.,.false.)

      sumabs=sumabs+cont/rosav(ntp)
!
! now HI bf absorption is included !

! add NH continuous absorption from Stancil. BPz 24/07-2019
      path='DATA/'
      pathlen=index(path,' ')-1
      NHtbl_unit=75
      opnh=nhop(xlambda,t(ntp),path,pathlen,NHtbl_unit)
! opnh is absorption in A^2/molecule
      if (nametryck(13).ne.'N H') then
        print*,nametryck(13), 'should be NH'
        stop 'PROBLEM in absko!'
      endif
! convert to cm^2/g of star
      sumabs = sumabs + opnh * 1.d-16 *
     &    partryck(ntp,13) / (t(ntp)*1.38066d-16) / rosav(ntp)
! NH cont included !
C
      if (J.gt.0) then
        ABSK(NTP)=SUMABS
        SPRID(NTP)=SUMSCA
      else
C
C        SUMMATION TO GET A ROSSELAND MEAN
ccc   25 CONTINUE
C        if(j.le.0) sumabs=sumabs + tioabs(ntp) + h2oabs(ntp)
        IF(J.EQ.0) ABSK(NTP)=ABSK(NTP)+ROSW(JP)*VIKTR/(SUMABS+SUMSCA)
        IF(J.LT.0) ABSK(NTP)=ABSK(NTP)+ROSW(JP)*VIKTR/SUMABS
        SUMW(NTP)=SUMW(NTP)+ROSW(JP)*VIKTR
      endif
! end of NTP loop
   26 CONTINUE
C
C        END OF 'THE SECOND NTP-LOOP'
C
      if (J.le.0) then
        DO NTP=1,NT
          SPRID(NTP)=0.
          ABSK(NTP)=SUMW(NTP)/ABSK(NTP)
        enddo
      endif
C
      RETURN
      END
