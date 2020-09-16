      PROGRAM eqwidt
*
*-----------------------------------------------------------------------
*
* Main program for lte line calculations. 
* Calculates equivalent widths for a line list and 
* iterate on abundance or computes curve fo growth for each line
* in the list. Line list format identical with bsyn.
* Based on Bsyn (turbospectrum) by BPz 06/07-2000
*
*-----------------------------------------------------------------------
*
      INCLUDE 'spectrum.inc'
      include 'tsuji.par'
*
      parameter (maxim=1000)
      INTEGER    TLUNIT,   LUNIT,alunit
      PARAMETER (TLUNIT=13,alunit=16)
      INTEGER    SPUNIT
      PARAMETER (SPUNIT=15)
      INTEGER    MAXNL,maxfil
      PARAMETER (MAXNL=200,maxfil=100)
      parameter (nmemol=16)
*
      CHARACTER*256 DETOUT,INATOM,INMOD,INLINE,INSPEC,OUTFIL,inabun,
     &              filterfil
      character*80 filttitle,comment
      character*50 MCODE

      real scattfrac,absfrac,mindnud
      real obs(3),abufact(3),eqwidth(3)
      REAL N,MA,MH,M,L,MUM,NTOT,MAM,MABUND(16),ntt
      REAL XIH,XIHM,XKHM,HJONH,HJONC,HJONN,HJONO,XNECNO
      REAL ABUVC(NDP),ABUVN(NDP),ABUVO(NDP),ABUVS(NDP),ABUVK(NDP)
      doubleprecision ionpot
      DIMENSION
     &          JLEV(NDP/5+1),EMOL(NMEMOL,NDP),BPLAN(NDP),
     &          X(NDP),S(NDP),plez(ndp),contop(ndp)
      doubleprecision XL1,XL2,DEL,XLMARG,XL1L,XL2R,XLBOFF,XLB,step
      doubleprecision DLAMB0,DOPPLC,DXLAMB,xlb_vshifted(ndp),lshift(ndp)
      CHARACTER*20 LELE
      real newvoigt
      logical formatlong

      COMMON/POP/ N(NDP),A(NDP),DNUD(NDP),STIM(NDP),QUO(NDP),DBVCON
      COMMON/ATOM/ XL,MA,CHI,CHI2,chi3,CHIE,G,IDAMP,FDAMP,
     &             GAMRAD,ALOGC6,ION
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),NTAU
      logical hydrovelo,firstiter,strongflag,outofrange,notconvflag(3),
     &        infoonly

      real velocity
      common/velo/velocity(ndp),hydrovelo
      COMMON/CONST/ BOLTZ,MH,H,C,E,M
      COMMON/CQ/ Q1(NDP),Q2(NDP),Q3(NDP),AQ1(3),AQ2(3),TLIM1,TLIM2,
     &           Q1LIM,Q2LIM,AQ3(3),TLIM3,Q3LIM,TQA(3)
      COMMON/TAUC/ TAU(NDP),DTAULN(NDP),JTAU
      COMMON/ROSSC/ ROSS(NDP),cross(ndp)
      COMMON/CWAVES/ XLS
      COMMON/CNUMB/ NTOT(NDP),ntt(ndp),fpartition(ndp),
     &              PH(NDP),HEH,phe(ndp),ph2(ndp)
      COMMON/ATMOL/ NAT,NMOL
      COMMON/MODID/ MCODE
      COMMON/PIECES/ XL1,XL2,DEL,EPS,NMY,NLBLDU,IINT,XMYC,IWEAK
      COMMON/UTPUT/ IREAD,IWRIT
      COMMON/CI5/ MABUND,ANJON(16,5),DUMT(94)
      COMMON/CHECK/ ABURC,ABURN,ABURO,ABURS,ABURK
* COMMONs for the carbon MOL()
      dimension presmo(30)
      COMMON/CMOL1/ EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
      COMMON/CMOL2/ NNMOL,PK(30)
*
* Special for spherical 
*
      COMMON /RHOC/RHO(NDP)
      COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP)
      COMMON /CSTYR/MIHAL  /CTAUM/TAUM
      logical spherical,limbdark,multidump,skiprelim
*
* extension for large number of wavelengths and lines (monster II)
      character*256 linefil(maxfil),mongofil,filprint
      integer     isotope(5),atom(5)
      character*26   species,blabla

      doubleprecision xlambda
      common/large/ xlambda(lpoint),maxlam,abso(ndp,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)
      logical oldpart,Ames,scan2001,oldscan,Barber
      common/oldpart/oldpart
      common/h2ochoice/Ames,scan2001,oldscan,Barber
*
      common/count/icount1,icount2,icount3,icount4
*
* common for continuum lambdas. Nlcont lambdas at
* xlambda(jlcont(1:nlcont)). 
      common/continuum/nlcont,jlcont(lpoint)
      real xx(lpoint),yyy(lpoint),zzz(lpoint)

      doubleprecision xlp
      common/babcont/xlp(20*numbset),nlq

* common for damping recipe
      real sigmacross,velexp,xlbr
      character*1 levlo,levup,recipe
      common/damping/sigmacross,velexp,recipe

*
* common for partial pressures
      logical tsuswitch,tsuji,chck
      doubleprecision parptsuji,presneutral,presion,presion2,presion3,
     &                partryck,xmettryck,xiontryck
      common /tsuji/ tsuji,tsuswitch,nattsuji,nmotsuji,
     &               parptsuji(maxim+400)
      character*128 filmet,filmol
      character*20 nametryck
      common /filetsuji/ filmet,filmol
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),nametryck(maxmol)
      common /orderedpress/ presneutral(ndp,100),presion(ndp,100),
     &                      presion2(ndp,100),presion3(ndp,100)
      real rhotsuji,xmytsuji,ejontsuji
      common/rhotsu/ rhotsuji,xmytsuji,ejontsuji

      logical dattsuji,datspherical,datlimbdark,databfind,abfind,
     &        datmultidump,datxifix,datmrxf,dathydrovelo,datpureLTE
      integer datnoffil,datncore,datmaxfil,datmihal,datiint
      real    isoch(1000),isochfact(1000),datisoch(1000),
     &        datisochfact(1000)
      real    datxmyc,datscattfrac
      character*128 datfilmet,datfilmol,datfilwavel
      character*256 datlinefil(maxfil),datdetout,
     &          datinatom,datinmod,datinabun,datcontinopac,datinpmod,
     &          datinspec,datoutfil,datmongofil,datfilterfil
      doubleprecision   datxl1,datxl2,datdel,datxlmarg,datxlboff
      common/inputdata/datmaxfil,dattsuji,datfilmet,datfilmol,
     &                 datnoffil,datlinefil,
     &                 datspherical,datmihal,dattaum,datncore,
     &                 datdiflog,datdetout,datinatom,
     &                 datinmod,datinabun,datinspec,
     &                 datoutfil,datmongofil,databch(100),
     &                 datlimbdark,datfilterfil,
     &                 datoverall,databfind,
     &                 datmultidump,datisoch,datisochfact,
     &                 dathelium,datalpha,datrabund,datsabund,
     &                 datxifix,datxic,datmrxf,datinpmod,datcontinopac,
     &                 datfilwavel,dathydrovelo,
     &                 datxl1,datxl2,datdel,datxlmarg,datxlboff,
     &                 datiint,datxmyc,datscattfrac,datpureLTE

      real amass(92,0:250),abund(92),fixabund(92),
     &         isotopfrac(92,0:250)
      real overall,alpha,helium,rabund,sabund
      character*2 aname(92)

      common/refabundances/ abund,amass,aname,isotopfrac


      common/filter/limbdark,ifilt,filtlam(1000),filttrans(1000),
     &              filterfil,filttitle
      common/abundch/abch(100)
      real absave(100),symmfactor

******************************************
      integer version
      data version /191/
******************************************
      data nat/92/
      logical newformat
      character oneline*256

ccc      external commn_handler

      tsuswitch =.false.
      Ames=.false.
      scan2001=.false.
      oldscan=.false.
      Barber=.false.
      lunit =23
      chck=.false.
      do k=1,100
        abch(k)=-99.9
      enddo
*
*
      datmaxfil=maxfil
*
      print*
      print*,'***********************'
      print10,version*0.1
10    format(' * EQWIDT version ',f4.1,' *')
      print*,'***********************'
      print*
      call input

*  fraction of the line opacity to count as scattering.
*  Remaining is counted in absorption    BPz 27/09-2002
      scattfrac=datscattfrac
      absfrac=1.0-scattfrac
      if (scattfrac.gt.0.) then
        print*
        print*,' WARNING!!!!! ', scattfrac,' of the line opacity',
     &     ' counted as scattering!!!!!'
        print*
      endif
      tsuji=dattsuji
      filmet=datfilmet
      filmol=datfilmol
      noffil=datnoffil
      do i=1,noffil
        linefil(i)=datlinefil(i)
      enddo
      hydrovelo=dathydrovelo
      spherical=datspherical
      mihal=datmihal
      taum=dattaum
      ncore=datncore
      diflog=datdiflog
      detout=datdetout
      inatom=datinatom
      inmod=datinmod
      inabun=datinabun
      inspec=datinspec
      outfil=datoutfil
      mongofil=datmongofil
      limbdark=datlimbdark
      multidump=datmultidump
      filterfil=datfilterfil
      overall=datoverall
      alpha=datalpha
      rabund=datrabund
      sabund=datsabund
      helium=dathelium
      abfind=databfind
      xl1=datxl1
      xl2=datxl2
      del=datdel
      xlmarg=datxlmarg
      xlboff=datxlboff
      iint=datiint
      xmyc=datxmyc

      print*,' EQWIDT CHECK: xl1= ',xl1
      print*,' EQWIDT CHECK: xl2= ',xl2
      print*,' EQWIDT CHECK: del= ',del
      print*,' EQWIDT CHECK: iint=', iint
      print*,' EQWIDT CHECK: xmyc=',xmyc

      do i=1,92
        fixabund(i)=databch(i)
      enddo
      do i=1,1000
        isoch(i)=datisoch(i)
        isochfact(i)=datisochfact(i)
      enddo

      print*,tsuji,filmet(1:index(filmet,' ')),
     &       filmol(1:index(filmol,' ')),
     &       noffil,(linefil(i),i=1,noffil),spherical,
     & mihal,taum,ncore,diflog,detout(1:index(detout,' ')),
     & inatom(1:index(inatom,' ')),inmod(1:index(inmod,' ')),
     & outfil(1:index(outfil,' '))
      do i=1,92
        if (fixabund(i).gt.-90.0) print*,'ab. changed: elt ',
     &      i,fixabund(i)
      enddo
*
* Tsuji molecular equilibrium?
*
      if (tsuji) then
       print*,'   Full molecular equilibrium, using:'
       print*,'     ',filmet(1:index(filmet,' '))
       print*,'     ',filmol(1:index(filmol,' '))
       tsuswitch=.true.
      else
       print*,'   Carbon Marcs molecular equilibrium'
      endif
*
* which kind of file for lines
*
      if (noffil.gt.maxfil) stop 'number of files too large. See maxfil'
      lunit=tlunit
*
* Read info for spherical
*
      print*,' '
      if (spherical) then
        print*,' Transfer treated in the spherical scheme'
      else
        print*,' Transfer treated in the plane parallel approximation'
      endif
      print*,' '
*
      OPEN(UNIT=7,FILE=DETOUT,STATUS='UNKNOWN')
      call clock
cc      OPEN(UNIT=12,FILE=INMOD,STATUS='OLD',recl=4*2*200)
      OPEN(UNIT=12,FILE=INMOD,STATUS='OLD')
cc      open(unit=16,file=inabun,status='old')
cc      OPEN(UNIT=27,STATUS='SCRATCH',FORM='FORMATTED')
ccc      OPEN(UNIT=20,FILE=OUTFIL,STATUS='UNKNOWN',FORM='UNFORMATTED')
cc      OPEN(UNIT=20,FILE=OUTFIL,STATUS='UNKNOWN',FORM='UNFORMATTED',
cc     &     RECL=412)
cc      OPEN(UNIT=66,FILE='sphlimb',STATUS='UNKNOWN',FORM='UNFORMATTED',
cc     &     RECL=412)
      OPEN(UNIT=23,STATUS='SCRATCH',FORM='FORMATTED')
*
* Initiate
*
      open(20,file=outfil,status='unknown')
      if (abfind) then
        print*,'WARNING! You may consider to rerun with the abundances'
        print*,'WARNING! found by eqwidt, in order to have the correct'
        print*,'WARNING! molecular equilibrium.'
        write(20,460)
        write(20,461)
      endif

      IREAD=5
      IWRIT=6
      IP=1
*
      BOLTZ=1.38066E-16
      MH=1.6735E-24
      M=9.1091E-28
      C=2.9979E10
      E=4.80298E-10
      H=6.6256E-27
      constant=SQRT(4.*ATAN(1.))*E*E/M/C
      IDAMP=2
      IELP=0
      IONP=0
*
* NALLIN is total no of lines used; NREJCT total no rejected
*
      NALLIN=0
      NREJCT=0
*
* DBVCON depthconstant doppler broadening velocity, if zero then micro-
* turbulence is used.
*
      DBVCON=0.0
*
* get abundances and scale them by overall. Then 
* if appropriate replace by fixabund.
* for molecular equilibrium calculation and damping
*
      call makeabund(overall,alpha,helium,rabund,sabund,fixabund,
     &                  abund,amass,aname,isotopfrac)

ccc      OPEN(UNIT=11,FILE=INATOM,STATUS='OLD')
      print*,'metallicity changed by ',overall,' dex'
      do i=2,92
        if (abund(i).lt.-28.) then
          abund(i)=0.0
        else
          abund(i)=10**(abund(i)-abund(1))
        endif
      enddo
      abund(1)=1.00
* change isotopic mixture if appropriate
      do i=1,1000
        if (isoch(i).gt.0.) then
          write(blabla,'(f6.3,20x)') isoch(i)
          print*,blabla
          call getinfospecies(blabla,iel,natom,atom,isotope)
          if (natom.gt.1)  then
            print*,'eqwidt: error in isotope specification',isoch(i)
            stop
          endif
          print*,'ISOTOPIC fraction changed from ',
     &            isotopfrac(atom(1),isotope(1))
          isotopfrac(atom(1),isotope(1))=isochfact(i)
          print*,'to ',isotopfrac(atom(1),isotope(1)),' for element ',
     &            atom(1)
        endif
      enddo
*
      mabund(1)=abund(1)
      mabund(2)=abund(2)
      mabund(3)=abund(6)
      mabund(4)=abund(7)
      mabund(5)=abund(8)
      mabund(6)=abund(10)
      mabund(7)=abund(11)
      mabund(8)=abund(12)
      mabund(9)=abund(13)
      mabund(10)=abund(14)
      mabund(11)=abund(16)
      mabund(12)=abund(19)
      mabund(13)=abund(20)
      mabund(14)=abund(24)
      mabund(15)=abund(26)
      mabund(16)=abund(28)
      HEH=MABUND(2)
      ABURC=MABUND(3)
      ABURN=MABUND(4)
      ABURO=MABUND(5)
      ABURS=MABUND(11)
      ABURK=MABUND(10)
*
* Calculate grams of hydrogen/grams of stellar matter (ABUNDH)
*
      WGT=1.008+4.003*HEH+12.01*ABURC+14.01*ABURN+16.00*ABURO+
     &    20.183*MABUND(6)+22.9898*MABUND(7)+24.312*MABUND(8)+
     &    26.9815*MABUND(9)+28.086*MABUND(10)+32.064*MABUND(11)+
     &    39.102*MABUND(12)+40.08*MABUND(13)+51.996*MABUND(14)+
     &    55.847*MABUND(15)+58.71*MABUND(16)
      ABUNDH=1.008/WGT
      print*,'old abundh:',abundh
*
* Start reading data on lines to be calculated
*
* IP>=1 gives lots of printout, IP=0 less, IP=2 gives some on term.
* The spacing of points is proportional to DOPPLC
* NL is the number of points
*
cc      call cstrip(alunit,27)
cc      READ(27,101) IP,EPS
      ip=0
      eps=0.0001
      iweak=0
      nmy=6
*
      if (iint.gt.0.and.spherical) then
        print*,'Spherical Eqwidt cannot be run in intensity mode.'
        Stop
cc        print*,'WARNING! Running in flux mode this time.'
cc        iint=0
      endif

      IF(ABS(XLMARG).LT.1.E-6) XLMARG=5.0000
      XLM=(XL1+XL2)/2.
      XL1L=XL1-XLMARG
      XL2R=XL2+XLMARG
      do 400 j=1,lpoint
       do 401 k=1,ndp
        abso(k,j)=0.0
        absos(k,j)=0.0
        absocont(k,j)=0.0
        absoscont(k,j)=0.0
401    continue
400   continue
       
*
* EPS tells which small l/kappa is neglected.
* XLM is a characteristic wavelength
* XL1L is the left consideration limit for lines
* XL2R is the right consideration limit for lines
*
cc      IF(EPS.LE.0.0) EPS=0.0001
cc      READ(27,102) XITE5
      xite5=0.0
      xlboff=0.0
cc      READ(27,109) XLBOFF
      if (xlboff.lt.1.e-6) xlboff=0.0
*
* XITE5  is depth independent microturbulence parameter in
* km/s. If this is zero or unspecified, XI is taken from
* model file, read by subr. READMO.
* XLBOFF is the wavelength offset for each line (Angstroms)
*
* Next read model atm. and continuous abs. coeffs. at XLM
* readmo also interpolates the continuous opacities to all the
* xlambda(1:maxlam) . They are put into absocont and absoscont. 
* the wavelengths that
* were used in babsma for the calculation of the cont. opacities
* are put into xlp.
*
* We provide a single fake xlm lambda for continuum.
      maxlam=-200
      xlm=(xl1+xl2)*0.5e0
      CALL READMO(XLM,X,S)
*
* X  is kappa/stndop,
* S  is sigma/stndop
* at all depths of the model, at wavelength xlm.
* the complete description of the continuum absorption is in 
* absocont (pure absorption) and absoscont (scattering) at 
* wavelengths xlambda(1:maxlam). They are also 
* divided by stndop. They are computed in readmo.
* stndop is called ross in this routine and is not ross.
*
      WRITE(6,237) MCODE(1:lenstr(mcode))
!      IF(IP.EQ.0) WRITE(7,303)
      NLEV=0
      DO 25 K=1,NTAU,5
        NLEV=NLEV+1
        JLEV(NLEV)=K
   25 CONTINUE
      JLEV0=NTAU
      IF(XITE5.GT.0) THEN
        XITE=XITE5*1.E5
        DO 318 K=1,NTAU
          XI(K)=XITE
318     CONTINUE
      END IF
* Calculate number densities of molecules, put in EMOL()
*
      tp=1.e6
      pep=-1.
      print*,'eqwidt; k, pgmod, pg_calc, romod, ro_calc'
      do 43 k=1,ntau
        if ((abs((t(k)-tp)/t(k)).lt.3.e-2).and.
     &      (abs((pe(k)-pep)/pe(k)).lt.0.6)) then
          skiprelim=.true.
        else
          skiprelim=.false.
        endif
        tp=t(k)
        pep=pe(k)
        call eqmol_pe(t(k),pg(k),pgpg,pe(k),1.,1.,k,niter,skiprelim)
c        print*,'eqmol_pe calculated ',niter,' iterations'
c        print*,k,pg(k),pgpg,ro(k),rhotsuji

ccc      write(*,'(i3,15e10.3,/,3x,15e10.3)') k,presmo
*
        PH(K)=presneutral(k,1)
        phe(k)=presneutral(k,2)
        ph2(k)=partryck(k,2)
43    continue
      abundh=1./xmytsuji
      print*,'new abundh:',abundh
*
*
* molecules in presmo/partryck:
*        I=1 H-, 2 H2, 3 H2+, 4 H2O, 5 OH, 6 CH, 7 CO, 8 CN, 9 C2, 10 N2,
*         11 O2, 12 NO, 13 NH, 14 C2H2, 15 HCN, 16 C2H, 17 -, 18 HS
*         19 SIH, 20 C3H, 21 C3, 22 CS, 23 SIC, 24 SIC2, 25 NS
*         26 SIN, 27 SIO, 28 SO, 29 S2, 30 SIS etc etc (see takemolec.f)
*
* molecules present in atomda are: 
*  CN, CH, C2, H2O, OH, NH, CO, N2, O2, H2, TiO, MgH,ZrO, HF
*                   et VO (12/2/96)
*                   et CaH (1/4/96) etc etc
*
* Here great loop over all elements and lines starts
*
* first, loop over line files
*

      do 98 ifil=1,noffil
*
        inline=linefil(ifil)
        OPEN(UNIT=13,FILE=INLINE,STATUS='old')
        call clock
        print*,' starting scan of linelist'
    1   CONTINUE
* security:
        ibadc6=0
        iel=-1
        ion=-1
        nline=0
*
        ielp=0
        read(lunit,*,end=9874) species,ion,nline
        read(lunit,*) comment
        call getinfospecies(species,iel,natom,atom,isotope)
* find out if Ames H2O or Joergensen's
        if (iel.eq.10108) then
          do ic=1,77
            if (comment(ic:ic+3).eq.'Ames') then
              Ames=.true.
              print*,'H2O line list from Ames. Using their partf.'
            endif
          enddo
          do ic=1,75
            if (comment(ic:ic+5).eq.'Barber') then
              Barber=.true.
              print*,'H2O line list from Barber.  Using their partf.'
            endif
          enddo
          do ic=1,73
            if (comment(ic:ic+7).eq.'scan2001') then
              scan2001=.true.
              print*,'H2O line list from SCAN2001. Using their partf.'
            endif
          enddo
          do ic=1,74
            if (comment(ic:ic+6).eq.'oldscan') then
              oldscan=.true.
              print*,'H2O line list from old SCAN. Using their partf.'
            endif
          enddo
        endif

* H I lines with Stark broadening. Special treatment.
        if (iel.eq.1) then
          print*,'eqwidt cannot treat H I lines '
          print*,' in particular the normalizing with the continuum '
          print*,' opacity is not made at the right wavelength'
          print*,' Also, iteration on abundance is not done'
          stop
ccc          lele='H '
ccc          print 1234,species,lele,iel,ion,(isotope(nn),nn=1,natom)
ccc          print*, 'nlines ', nline
ccccc          call Hlineadd(lunit,nline,xlboff)
ccc          call hydropac(lunit,xlboff)
ccc          goto 9874
        endif

*default fdamp (may be changed for each element further down)
****        if (iel.gt.92) then
****          fdamp=0.
****        else
****          fdamp=2.0
****        endif
        g=1.
ccc        raddmp=0.0
* new element. we must compute more in depth
        oldpart=.false.
        symmfactor=1.e0
        IF(IEL.LE.0) GOTO 9874
*
* IEL  is identification of the species (e.g. 3 = Li, 822 = TiO, 10108 = H2O)
*      defined in getinfospecies. 
*
* ION  is the stage of ionization (=1 for neutral atoms or
*      molecules, 2 for singly ionized species etc)
*
*The abundance is by default the
* Anders and Grevesse abundance stored in makeabund.f)
* scaled by overall, helium, etc (see options in input.f). 
* This value may be overwritten by using ABCHANGE in input file.
*
* We don't read the atomic/molec data table anymore (atomda): abunp not used.
* Partition functions
* come now from partf and molecpartf. Chi, chi2, chi3 come from partf.f for 
* atoms (Irwin data tables).
* CHI is not needed for molecules anymore! as we don't want to use the 
* Unsoeld recipe for vanderWaals broadening. 
* Finally the mass is computed here.
*
ccc        CALL ATOMDA(IEL,LELE,CHI,CHI2,MAM,ABUNP)

        if (iel.le.92) then
          lele=aname(iel)
        else
          call getlele(iel,ion,lele)
        endif
        mam=0.
        do i=1,natom
* compute mass for isotopomer. If isotope is not specified (isotope=0), 
* then standard mix is assumed, with mass=amass(atom(i),0)
          mam=mam+amass(atom(i),isotope(i))
        enddo
        MA=MAM*1.6603E-24
        print*,'Mass used for lines of ',iel,lele,' is ', mam

        print 1234,species,lele,iel,ion,(isotope(nn),nn=1,natom)
        print*, 'nlines ', nline
1234    format('species: ',a17,1x,a20,' iel: ',i8,' ion: ',i2,
     &         ' isotopes: ',5i3)
        print*,comment
        if (iel.le.92) then
          ABUL=abund(iel)
          abunp=abund(iel)
        else
* abunp is used when searching for abundance from equivalent width.
* we arbitrarily set it to zero for molecules that don't have "abundances".
* This avoids "inf" in output
          abunp=1.e-12
        endif
* Test for molecular list format
* partly allows backward compatibility for pre-v14.1 format molecular line lists
* can compute eqw for molecules, but not iterate on "abundance" as there is no eqw in the input
* So, abfind should be .false.
        if (iel.gt.92) then
          read(lunit,'(a)') oneline
          backspace(lunit)
          read(oneline,*,err=11,end=11) xlb,chie,gfelog,fdamp,gu,
     &                raddmp,levlo,levup,obseqw,eqwerror
          newformat=.true.
          goto 12
11        newformat=.false.
12        continue 
        else
          newformat=.true.
        endif
*
* Start wavelength loop
*
      ILINE=0
      NALLIN=NALLIN+NLINE
*
* NLINE is the number of lines of the element IEL.
*
***************************************************************
*
*
***************************************************************
*
* Big jump to 64 from far below, line loop.
*
   64 continue
      dgfe=0.
*
* Loop for reading the lines
*
   50 CONTINUE
ccc        IF(IP.LE.0) WRITE(7,303)
ccc        IF(IP.LE.0) WRITE(7,302) LELE,ION,XLB,CHIE,XITE5,FDAMP,GFELOG
      IF(ILINE.EQ.NLINE) then
        print*,iline,' considered for element ',iel,ion
        print*,ibadc6,' lines rejected because of negative c6'
        GOTO 1
      ENDIF
*
* Now read the line data
*
* XLB= wavelength
* CHIE= excitation pot (in eV) of lower level
* GFELOG= log(gf)
* GU the upper statistical weight for the line, is only of
*    importance for damping and if raddmp is not 0.
* FDAMP  is a fudge factor to increase damping constant
* RADDMP externaly calculated radiation damping (if needed)
* f = f-value (as g=1., f is in fact gf-value).
*
* warning! xlb is real*8
*
* new format for molecules, identical to that for atoms, starting with v14.1
      if (newformat) then
        read(lunit,*) xlb,chie,gfelog,fdamp,gu,raddmp,levlo,levup,
     &                obseqw,eqwerror
      else
* allows backward compatibility for older format molecular line lists
        read(lunit,*) xlb,chie,gfelog,fdamp,gu,raddmp
      endif
      outofrange=.false.
*
* first, reject line if outside the wavelength range of babsma.
      if (xlb.le.xl1.or.xlb.ge.xl2) then
cc        print*,' line out of range for continuous opacity. Dropped.'
        outofrange=.true.
        abufact(2)=1.0
        NREJCT=NREJCT+1
        goto 2345
      endif
*
cc      if (iel.le.nat) then
      velexp=0.0
      sigmacross=0.0
      if (abfind) then
* we determine abundance for obseqwidt and obs+-error
        ieqmin=1
        ieqmax=3
      else
* otherwise we only do one computation
        abufact(1)=1.0
        abufact(3)=1.0
        ieqmin=2
        ieqmax=2
      endif
      if (obseqw.gt.0.) then
        obs(1)=max(obseqw-eqwerror,0.1)
        obs(2)=obseqw
        obs(3)=obseqw+eqwerror
        if (eqwerror.le.0.0) then
* do only one computation for that one, at nominal abundance.
          ieqmin=2
          ieqmax=2
        endif
      else 
        obs(1)=0.0
        obs(2)=0.0
        obs(3)=0.0
      endif
*
      f=10**(gfelog)
      xlb=xlb+xlboff

      if (.not.abfind.or.obseqw.gt.0.0) then
*
* if abfind and obseqw=0. we skip the line and do not output anyhting
* see tag "123" below
* Otherwise, 
*   start line calculations
*
      if (IEL.ne.IELP) then 
*
* Calculate abundance of molecule/atom per gram stellar matter
*
        symmfactor=1.e0
        infoonly=.false.
        if (iel.gt.nat) then
          molindex=0
          infoonly=.true.
          call takemolec(1,infoonly,lele,molindex)
          if (lele.eq.'C C                 ') then
            if ((isotope(1).eq.12.and.isotope(2).eq.13).or.
     &          (isotope(2).eq.12.and.isotope(1).eq.13)) then
              symmfactor=2.
              print*,' SYMMFACTOR = 2 for 12C13C'
            else if ((isotope(1).eq.12.and.isotope(2).eq.12).or.
     &               (isotope(1).eq.13.and.isotope(2).eq.13)) then
              symmfactor=1.
            else if (isotope(1).eq.0.and.isotope(2).eq.0) then
              symmfactor=1.
            else
              stop 'Eqwidt: Problem with C2 isotopic mix!!'
            endif
          endif
          if (molindex.eq.0) then
           print*,'eqwidt: molecular species not implemented in atomda',
     &           lele
            stop
          endif
          do k=1,ntau
            ntot(k)=partryck(k,molindex)/boltz/t(k)/ro(k)*symmfactor
            ntt(k)=partryck(k,molindex)/boltz/t(k)/ro(k)*symmfactor
          enddo
          call partffordepth(ntau,t,lele,fpartition)
        else 
          do k=1,ntau
            if (presneutral(k,iel).ge.0.) then
              ntot(k)=(presneutral(k,iel)+presion(k,iel)+
     &                 presion2(k,iel)+presion3(k,iel))/
     &               boltz/t(k)/ro(k)
              if (ion.eq.1) then
                ntt(k)=presneutral(k,iel)/boltz/t(k)/ro(k)
              else if (ion.eq.2) then
                ntt(k)=presion(k,iel)/boltz/t(k)/ro(k)
              else if (ion.eq.3) then
                ntt(k)=presion2(k,iel)/boltz/t(k)/ro(k)
              else if (ion.eq.4) then
                ntt(k)=presion3(k,iel)/boltz/t(k)/ro(k)
              endif
              call partf(iel,1,t(k),1,fpartition(k),ionpot)
              chi=ionpot
              call partf(iel,2,t(k),1,fpartition(k),ionpot)
              chi2=ionpot
              call partf(iel,3,t(k),1,fpartition(k),ionpot)
              chi3=ionpot
              call partf(iel,ion,t(k),1,fpartition(k),ionpot)
            else
              if (k.eq.1) then
                print*,'element not present in chemical equilibrium',
     &          ' adopted abundance: ',log10(abunp)+12.
              endif
              ntot(k)=abunp*abundh/mh
              ntt(k)=-1.0
            endif
          enddo
        endif

        do k=1,ntau

* prepare line shift vs. depth, from velocity in model (in cm/s)
* velocity in model should be positive outwards.
* If xlb_vshifted is the wavelength in the observer's frame at which
* the line position is shifted [(lambda_0-lambda)/lambda_0=v/c], 
* xlb_vshifted=xlb*lshift, with lshift calculated here:

          lshift(k)=1.d0-velocity(k)/2.99792458d10

          do i=1,natom
            ntot(k)=ntot(k)*isotopfrac(atom(i),isotope(i))
            ntt(k)=ntt(k)*isotopfrac(atom(i),isotope(i))
            if (ntot(k).eq.0.) then
              print*,'eqwidt. WARNING!, ntot=0 for species: ',lele
              print*,'atom=',atom(i),' isotope=',isotope(i)
              print*,'isotopfrac =',isotopfrac(atom(i),isotope(i))
            endif
          enddo
***          ntot(k)=max(ntot(k),1.e-30)
        enddo
 
      endif
*
* Print line information
*
ccc      IF(IP.GE.1) WRITE(7,990)IEL,ION,NMY,FDAMP,IINT,IMY,MAM,ABUL,
ccc     &                       CHI,XITE
      CHIU=CHIE+3.40*3647./XLB
      if (iel.le.nat) then
        XIONP=CHI
        IF(ION.EQ.2) XIONP=CHI2
        IF(ION.EQ.3) XIONP=CHI3
        IF(ION.GT.3) THEN
          PRINT *,
     &       '*****************************************************'
          PRINT *,
     &       'Error in Eqwidt, ION.GT.3 which has not been foreseen'
          PRINT *,
     &       '*****************************************************'
          STOP '***** Eqwidt *****'
        END IF
      endif
*
* Calculate damping parameters
*
      if (fdamp.ge.20.) then
*
* BPz 02/06-2014
* 1) use ABO theory (Anstee, Barklem, O'Mara) for collisional damping with H,
* with data taken from line list: fdamp contains sigma.alpha.
* This number is available starting with VALD3 version of the VALD database.
* See : http://www.astro.uu.se/~barklem/howto.html
* 2) if (1) not available check if something can be computed in the anstee.f
* routine
* 3) if (2) not available, check in linelist for a gamma6 at 10000K
* 4) if nothing else worked, comput Unsoeld approximation.
*
        sigmacross=int(fdamp)
        velexp=fdamp-int(fdamp)
        recipe='S'
*
      else IF (FDAMP.GT.0..and.fdamp.lt.20.) THEN
* We may use Unsoeld theory with fudge factor fdamp, and prepare for it.
        XXXXX=ION**2*(1./(XIONP-CHIU)**2-1./(XIONP-CHIE)**2)
      END IF

      if (RADDMP.NE.0.) then
*       Use radiative damping data in line list if available
        GAMRAD=RADDMP
      else
*       default recipe for radiative damping
        GAMRAD=6.669E15*G*F/(GU*XLB**2)
      endif

      XL=XLB*1.D-8
*
* check whether there are quantum mechanical damping data for this line
*  (atoms only)
* Unsoeld recipe for atomic lines is default, but Barklem et al.'s treatment 
* may be set in anstee.f
* we should not use the Unsoeld recipe for collisional broadening of molecular lines!
*
      if (iel.le.nat) then
        xlbr=xlb
        idamp=2
        if (fdamp.lt.20.) then
*
* We call anstee only if we don't have quantum mechanical collisional data
* in the line list
*
          call anstee(iel,ion,xlbr,chie,xionp,sigmacross,velexp,levlo,
     &              levup,recipe)
        endif
        if (recipe.eq.'U') then
* if the transition collisional broadening is not handled by anstee.f (e.g.
* for x ->x transitions), we have recipe='U' (Unsoeld approximation)
* However we may have gamma van der Waals data in the line list, which we then use

          if (fdamp.lt.0.) then
* fdamp contains log(gammavdW at 10000K) instead of the fudge factor for Unsoeld 
* recipe. BPz 08/04-2013

            recipe='W'
          else if (xxxxx.le.0.) then
* skip the line
            ibadc6=ibadc6+1
            NALLIN=NALLIN-1
            NLINE=NLINE-1
            NREJCT=NREJCT+1
            goto 50
          else
* Unsoeld approximation
            ALOGC6=ALOG10(XXXXX)-29.7278
          endif
          ALOGC6=ALOG10(XXXXX)-29.7278
        endif
      else
        idamp=2
ccc        recipe='U'
ccc        recipe='T'
        recipe='R'
* we cannot use Unsoeld recipe for molecules !
* 'R' is for pure radiative damping
      endif
*
* Calculate occupation numbers
*
   54 CONTINUE
      CALL DEPTH(IEL)
      IELP=IEL
      IONP=ION
cc      IF(IP.GE.1) WRITE(7,230)
cc      DBVK=XL*1.E-05
cc      DO 6 JJ=1,NTAU,10
cc        DBV=DBVK*DNUD(JJ)
cc        IF(IEL.LE.NAT)  ANTEL=N(JJ)/(ABUND/MH/MUM(JJ))
cc        IF(IEL.GT.NAT) ANTEL=N(JJ)/NTOT(JJ)
cc        PARTK=Q1(JJ)
cc        IF(ION.EQ.2) PARTK=Q2(JJ)
cc        IF(IP.GE.1) WRITE(7,231)JJ,TAU(JJ),QUO(JJ),ANTEL,A(JJ),
cc     &                         PARTK,DBV,STIM(JJ),XC(JJ),S(JJ)
cc  6   CONTINUE
      DLAMB0=DOPPLC*DNUD(JLEV0)*XL**2/C
*
* Constants
*
      CALF=constant*F
*
cc      IF(FDAMP.GT.0..AND.IP.GE.1) WRITE(7,267) GAMRAD,ALOGC6
cc      IF(IP.GE.1) WRITE(7,265) LELE,ABUL,ABUND,NTAU
      mindnud=1.e30
      do 111 j=1,ntau
       plez(j)=n(j)*stim(j)/dnud(j)/ross(j)
       xlb_vshifted(j)=xlb*lshift(j)
       mindnud=min(mindnud,dnud(j))
111   continue
*
* Start wavelength loop for this line
*
* Reference opacity for line contribution
      xlmb=xlb
      call readmo(xlmb,x,s)
      do j=1,ntau
cc        contop(j)=x(j)+s(j)
        contop(j)=x(j)
      enddo

      strongflag=.false.

! adaptive step for profile calculation, given in fraction of Doppler width.
! minimum value is del in [AA] given in input.

      step=min(del,mindnud*xl*xlb/c*0.5)
      print*,'lambda',xlb,' step for calculation ',step,'A'

      do ieq=ieqmin,ieqmax
! ieq=2 is for nominal abundance, ieq=1 for lower error bar, ieq=3 for upper error bar

        notconvflag(ieq)=.false.
        firstiter=.true.
        iterabu=1
        if ((ieq.eq.1).or.(ieqmin.eq.ieqmax)) then
          abufact(ieq)=1.0
        else
          abufact(ieq)=obs(ieq)/eqwidth(ieq-1)/1000.*abufact(ieq-1)
        endif
        abufactold=2.0

9876    continue
        lmin=1
        lmax=1
        do iloop=int(lpoint/2)+1,lpoint,1
          xkmax=0.
          i=iloop-int(lpoint/2)-1
          xlambda(iloop)=xlb+float(i)*step
          if (xlambda(iloop).gt.xl2) then 
            lmax=iloop-1
            goto 155
          endif
          vt=(float(i)*step)*1.d-8
          vt=c*vt/xl**2
          do j=1,ntau
            v=vt/dnud(j)
c           CALL VOIGT(A(j),V,HVOIGT)
            hvoigt=newvoigt(a(j),v)
            l=calf*hvoigt*plez(j)
            l=l*abufact(ieq)
            abso(j,iloop)=l*absfrac
            absos(j,iloop)=l*scattfrac
            xkmax=max(xkmax,l/contop(j))
          enddo
          lmax=iloop
          if (xkmax.lt.eps.and.iloop.gt.int(lpoint/2)+1) goto 15
        enddo
155     continue
        print*,' element: ', lele, ion
        print*,' the line may still contribute at ',xl2
        print*,' central lambda: ',xlmb
        print*,' abufact: ',abufact(ieq)
        print*,' obseqw: ',obs(ieq),' eqwidth: ',eqwidth(ieq)*1000.
        strongflag=.true.

15      CONTINUE
* and now the other side of the profile
        do iloop=int(lpoint/2),1,-1
          xkmax=0.0
          i=int(lpoint/2)+1-iloop
          xlambda(iloop)=xlb-float(i)*step
          if (xlambda(iloop).lt.xl1) then 
            lmin=iloop+1
            goto 165
          endif
          vt=(float(i)*step)*1.d-8
          vt=c*vt/xl**2
          do j=1,ntau
            v=vt/dnud(j)
c           CALL VOIGT(A(j),V,HVOIGT)
            hvoigt=newvoigt(a(j),v)
            l=calf*hvoigt*plez(j)
            l=l*abufact(ieq)
            abso(j,iloop)=l*absfrac
            absos(j,iloop)=l*scattfrac
            xkmax=max(xkmax,l/contop(j))
            lmin=iloop
          enddo
          if (xkmax.lt.eps) goto 16
        enddo
165     continue
        print*,' element: ', lele, ion
        print*,' the line may still contribute at ',xl1
        print*,' central lambda: ',xlmb
        print*,' abufact: ',abufact(ieq)
        print*,' obseqw: ',obs(ieq),' eqwidth: ',eqwidth(ieq)*1000.
        strongflag=.true.
  
16      CONTINUE
*
cc        print*, 'EQWIDT: Line at ', sngl(xlb)
cc        print*, 'EQWIDT: lambda min and max ', 
cc       &       lmin,lmax,sngl(xlambda(lmin)),sngl(xlambda(lmax))
*
* Continuum points for continuum flux calculations
* call readmo again to get continuum absorption and scattering coefficients
* in abso and absos (common "large") at the wavelengths of the line profile.
* The line opacity is in abso(lmin:lmax),and absos(lmin:lmax)
* The continuum opacity, after interpolation in readmo, is in absocont(1:maxlam)
* and absoscont(1:maxlam) with maxlam=lmax-lmin+1
* 
        maxlam=lmax-lmin+1
        do i=1,maxlam
          xlambda(i)=xlambda(lmin+i-1)
        enddo
        call readmo(xlm,x,s)
*
* X  is kappa/stndop,
* S  is sigma/stndop
* stndop is called ross in the routine readmo but is not ross.
*
*
* We shift everything between (1:maxlam), and set lmin=1, lmax=maxlam
* we also add the continuum extinction to the line extinction
*
        do i=1,maxlam
          do k=1,ntau
            abso(k,i)=abso(k,lmin+i-1)+absocont(k,i)
            absos(k,i)=absos(k,lmin+i-1)+absoscont(k,i)
          enddo
        enddo
        lmin=1
        lmax=maxlam

        eqwidth(ieq)=0.0
        if (lmax.gt.lmin) then
          if (spherical) then
            call eqwidtb(lmin,lmax,step,eqwidth(ieq))
          else
            call eqwidtbplatt(lmin,lmax,step,eqwidth(ieq))
          endif
        endif
*
* iterate on abundance; epsabu = precision of determination in dex
*
        epsabu=0.001
        epsabu=log(10**epsabu)

        if (abfind) then
          if (obs(ieq).gt.0.e0) then
ccc          if (abs(log(obs(ieq)/eqwidth(ieq)/1000.)).gt.epsabu) then
            if (abs(log(abufact(ieq)/abufactold)).gt.epsabu.and.
     &        iterabu.le.14) then
              if (firstiter) then
                abufactold=abufact(ieq)
                abufact(ieq)=obs(ieq)/eqwidth(ieq)/1000.*abufact(ieq)
                eqwidthold=eqwidth(ieq)
                firstiter=.false.
              else 
                dlnabufact=log(abufact(ieq)/abufactold)/
     &                log(eqwidth(ieq)/eqwidthold)*
     &                  log(obs(ieq)/eqwidth(ieq)/1000.)

                if (iterabu.eq.8) then
c we are likely stuck in an infinite loop
c we try to get out of it with a small kick
                  dlnabufact=dlnabufact+1.1
                endif
  
                abufactold=abufact(ieq)
                abufact(ieq)=abufactold*exp(dlnabufact)
                eqwidthold=eqwidth(ieq)
                iterabu=iterabu+1
              endif
cc            print*, iterabu,abufactold,
cc     &              obs(ieq),eqwidth(ieq)*1000.,exp(dlnabufact)
cc            print*
              goto 9876
            else if (iterabu.gt.14) then
c obviously the kick at iteration 8 was not enough...
              print*,' WARNING!! abundance not converged for line:'
              print*,LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                  obs(2),
     &                  log10(abunp)+12.+log10(abufact(2))
* fake abufact
              abufact(ieq)=1.0
              notconvflag(ieq)=.true.
cc            print*, iterabu,abufactold,
cc     &          obs(ieq),eqwidth(ieq)*1000.,exp(dlnabufact)
cc            print*
            endif
          else
* observed equivalent width is zero mAA
            abufact(ieq)=1.e-30
          endif
        endif
* end of the ieq loop
      enddo

* Write line data
*
      gfelog=log10(f)
      if (gfelog.lt.-9.99) then
        formatlong=.true.
      else 
        formatlong=.false.
      endif
      if (strongflag) then
        if (formatlong) then
          WRITE(20,455) LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                  obs(2),eqwerror,log10(abufact(2)),
     &                (log10(abunp)+12.+log10(abufact(ieq)),ieq=1,3)
        else
          WRITE(20,457) LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                  obs(2),eqwerror,log10(abufact(2)),
     &                (log10(abunp)+12.+log10(abufact(ieq)),ieq=1,3)
        endif
        strongflag=.false.
      else if (notconvflag(1).or.notconvflag(2).or.notconvflag(3)) then
        if (formatlong) then
          WRITE(20,4572) LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                  obs(2),eqwerror,log10(abufact(2)),
     &                (log10(abunp)+12.+log10(abufact(ieq)),ieq=1,3)
        else
          WRITE(20,4571) LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                  obs(2),eqwerror,log10(abufact(2)),
     &                (log10(abunp)+12.+log10(abufact(ieq)),ieq=1,3)
        endif
C Note that we could try to recover 1 or 2 abundances if 
C they were converged. Here we throw away everything, even if only
C one of the cases did not work.
        do ieq=1,3
          notconvflag(ieq)=.false.
        enddo
      else
        if (eqwerror.le.0.0 .or. .not.abfind) then
* Case : find abundance, only for nominal measured eqw. No error
* on measured eqw.
* Or do not find abundance. Print same format, with all abundances the same.
          abufact(1)=abufact(2)
          abufact(3)=abufact(2)
        endif
        if (formatlong) then
          WRITE(20,4561) LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                obs(2),eqwerror,log10(abufact(2)),
     &                (log10(abunp)+12.+log10(abufact(ieq)),ieq=1,3)
        else
          WRITE(20,456) LELE,ION,XLB,CHIE,GFELOG,eqwidth(2)*1000.,
     &                obs(2),eqwerror,log10(abufact(2)),
     &                (log10(abunp)+12.+log10(abufact(ieq)),ieq=1,3)
        endif
      endif
* We arrive here, if obseqw=0. and abfind. We have skipped the calculations
* and we don't print anything
123   endif

*
* End of line calculation
*
2345  continue
      iline=iline+1
*
      goto 64
*
* End model loop
*
* end of line lists loop
9874  close(lunit)
98    continue
* 
      WRITE(6,214) NREJCT,XL1L,XL2R,XLM
 214  FORMAT(1X,I8,' LINES WERE REJECTED, ONLY LINES BETWEEN',F10.3,
     &       ' AND',F10.3,' A CONSIDERED',/,'  CENTRAL WAVELENGTH=',
     &       F10.3,' A')
*
  100 FORMAT(4X,I1,6X,I3)
  101 FORMAT(3X,I1,5X,F6.4)
  102 FORMAT(6X,F4.0)
  103 FORMAT(4X,I3,5X,I1,7X,F6.2)
 1030 format(a3,17x,1x,e15.8,f9.3,1x,f6.3)
 1031 format(a3,7x,i3,7x,1x,e15.8,f9.3,1x,f6.3)
 1032 format(a3,i3,i3,1x,i3,7x,1x,e15.8,f9.3,1x,f6.3)
  107 FORMAT(7X,I3)
  109 FORMAT(7X,F6.3)
  137 FORMAT(7X,F10.0)
  212 FORMAT(4X,I1,7X,I1)
  213 FORMAT(' HE=',F6.2,'  C=',F6.2,'  N=',F6.2,'  O=',F6.2)
  230 FORMAT('0JLEV',3X,'  TAU ',8X,'ANJON',6X,'ANTEL',6X,'DAMP',7X,
     &       'PART',7X,'DBV',8X,'STIM',7X,'X',11X,'S')
  231 FORMAT(I5,1P10E11.3)
  235 FORMAT(2X,'CORRECTIONS FOR STIMULATED EMISSION'/20X,6F15.4)
  236 FORMAT(' ***STOP IN eqwidt***.NTOT(',I2,').LT.0.0, IEL=',
     & I3,/,' ION=',I3,' ABUND=',E10.3,' ILINE=',I4,' XLB=',F9.2)
  237 FORMAT(' MODEL IDENTIFICATION=',A,'; MAIN PROGRAMME eqwidt')
  240 FORMAT(' **DATA FOR LAMBDA',F10.3,'  FROM ',A2,I1,
     &       ' WRITTEN ON UNIT 14')
  265 FORMAT('0THE ABUNDANCE OF ',A3,' IS',F6.2,
     &       ' (NO OF FREE NUCLEI PER HYDROGEN: ',1PE9.2,
     &       ' AT DEPTHPOINT',I3,' )')
  266 FORMAT('0'/,' CHEMICAL COMPOSITION'//16(2X,A4)//16F6.2)
  267 FORMAT('0THE DAMPING WAS COMPUTED USING GAMRAD =',1PE9.2,
     &       ' AND LOG C6 =',0PF7.2)
  301 FORMAT()
 3301 FORMAT(D10.3,4F10.3,2E10.2)
  300 FORMAT(A2)
  302 FORMAT(1X,A2,I1,F6.0,F5.2,2F4.1,F6.3)
  303 FORMAT(' EL   XLB  CHIE   X  FD   GFE')
  460 format('       lambda    Exc  log(gf)           obs.eqw+- error',
     &'   delta    lower    abund    upper')
  461 format('                                                       ',
     &'   abund.   abund             abund')
 4572 format(a3,1x,i1,1x,f13.3,1x,f5.2,1x,f6.2,2(1x,f8.3),' +- ',f6.3,
     &       1x,f6.2,3(1x,f8.3),' abu not converged ')
 4571 format(a3,1x,i1,1x,f13.3,1x,f5.2,1x,f6.3,2(1x,f8.3),' +- ',f6.3,
     &       1x,f6.2,3(1x,f8.3),' abu not converged ')
 455  format(a3,1x,i1,1x,f13.3,1x,f5.2,1x,f6.2,2(1x,f8.3),' +- ',f6.3,
     &       1x,f6.2,3(1x,f8.3),' strong? wings cut?')
 457  format(a3,1x,i1,1x,f13.3,1x,f5.2,1x,f6.3,2(1x,f8.3),' +- ',f6.3,
     &       1x,f6.2,3(1x,f8.3),' strong? wings cut?')
 4561 format(a3,1x,i1,1x,f13.3,1x,f5.2,1x,f6.2,2(1x,f8.3),' +- ',f6.3,
     &       1x,f6.2,3(1x,f8.3))
  456 format(a3,1x,i1,1x,f13.3,1x,f5.2,1x,f6.3,2(1x,f8.3),' +- ',f6.3,
     &       1x,f6.2,3(1x,f8.3))
  990 FORMAT(' INPUT PARAMETERS:'/' IEL=',I3,
     & ' ION=',I1,' NMY=',I1,' FDAMP=',F3.1,' IINT=',I1,
     & ' IMY=',I1,' MA=',F6.2,' ABUND=',F4.2,' CHI=',F5.2,' XITE=',
     & 1PE8.2)
*
      END
