      PROGRAM BSYN
*
*-----------------------------------------------------------------------
*
* Main program for lte line calculations. The program calculates
* several lines from one element.
* Author: Kjell Eriksson, Astron. Obs., S-75120 Uppsala, Sweden.
* Modified by: Ake Nordlund, Nordita, Blegdamsvej 17, DK-2100 Copenhagen
*              Denmark.
* 76.03.30
* The Canarian version by B. Gustafsson, Uppsala 82.01.22
*
* Modified by                  K.Eriksson and O. Morell  Uppsala 1984.09
* ECNO --> MOL, F77 + kosmetics               O. Morell  Uppsala 1986.09
* New input file format +
* Double precision & more F77                 O. Morell  Uppsala 1986.10
* Export version  1988-03-24           ********* Olof Morell *** Uppsala
* Small improvements on NL & dimensions       O. Morell  Uppsala  890530
*
* Major upgrade: new jon and detabs blocks from MARCS. major revisions
*    in main bsyn below. BPz 29/10-1996
*
* New jon/eqmol_pe_lu/die_pe_lu . Solves chemical equilibrium with 
* matrix inversion.  BPz 09/07-1999
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

      real scattfrac,absfrac
      REAL N,MA,MH,M,L,MUM,NTOT,MAM,MABUND(16),ntt
      REAL XIH,XIHM,XKHM,HJONH,HJONC,HJONN,HJONO,XNECNO
      REAL ABUVC(NDP),ABUVN(NDP),ABUVO(NDP),ABUVS(NDP),ABUVK(NDP)
      doubleprecision ionpot
      DIMENSION
     &          JLEV(NDP/5+1),EMOL(NMEMOL,NDP),BPLAN(NDP),
     &          X(NDP),S(NDP),plez(ndp),contop(ndp)
      doubleprecision XL1,XL2,DEL,XLMARG,XL1L,XL2R,XLBOFF,XLB
      doubleprecision xl1l_tmp, xl2r_tmp !SBC
      integer not_discard !SBC
      doubleprecision xlb_vshifted(ndp),lshift(ndp)
      CHARACTER*20 LELE
      real newvoigt

      COMMON/POP/ N(NDP),A(NDP),DNUD(NDP),STIM(NDP),QUO(NDP),DBVCON
      COMMON/ATOM/ XL,MA,CHI,CHI2,chi3,CHIE,G,IDAMP,FDAMP,
     &             GAMRAD,ALOGC6,ION
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),NTAU
      logical hydrovelo,debug,infoonly,computeIplus
      real velocity
      common/velo/velocity(ndp),hydrovelo,computeIplus
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

      doubleprecision xlambda,dist
      common/large/ xlambda(lpoint),maxlam,ABSO(NDP,lpoint),
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

      logical dattsuji,datspherical,datlimbdark,databfind,
     &        datmultidump,datxifix,datmrxf,dathydrovelo,datpureLTE,
     &        pureLTE
      integer datnoffil,datncore,datmaxfil,datmihal,datiint
      real    isoch(1000),isochfact(1000),datisoch(1000),
     &        datisochfact(1000)
      real    datxmyc,datscattfrac
      character*128 datfilmet,datfilmol,datfilwavel
      character*256 datlinefil(maxfil),datdetout,
     &          datinatom,datinmod,datinabun,datcontinopac,datinpmod,
     &          datinspec,datoutfil,datmongofil,datfilterfil
      doubleprecision  datxl1,datxl2,datdel,datxlmarg,datxlboff
      integer nsegments_lambda_min, nsegments_lambda_max !SBC
      doubleprecision segments_lambda_min(400), segments_lambda_max(400)
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
     &                 datiint,datxmyc,datscattfrac,datpureLTE,
     &                 nsegments_lambda_min,nsegments_lambda_max,
     &                 segments_lambda_min, segments_lambda_max

      real amass(92,0:250),abund(92),fixabund(92),
     &         isotopfrac(92,0:250)
      real overall,alpha,helium,rabund,sabund
      character*2 aname(92)

      common/refabundances/ abund,amass,aname,isotopfrac


      common/filter/limbdark,ifilt,filtlam(1000),filttrans(1000),
     &              filterfil,filttitle
      common/abundch/abch(100)
      real absave(100),symmfactor
***********************************************
      integer version
      data version /191/
***********************************************
      data debug/.false./
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
      datmaxfil=maxfil
*
      print*
      print*,'**************************************************'
      print10,version*0.1
10    format(' * BSYN version ',f4.1,'                              *')
      print*,'* Warning, now Flux is so that F=sigma.Teff^4 !  *'
      print*,'* consistent with MARCS.                         *'
      print*,'* Versions 10.1 and lower give f=sigma.Teff^4/pi *'
      print*,'**************************************************'
      print*
 
      call input

* fraction of the line opacity to count as scattering.
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
      computeIplus=.false.
      if (hydrovelo) then
* We make sure that every layer will have an intercepting ray.
* The forward computation of flux in extended atmospheres is not
* so clever and produces wrong results if the source function
* is not well sampled.
        diflog=1.0001
        computeIplus=.true.
      endif
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
      xl1=datxl1
      xl2=datxl2
      del=datdel
      xlmarg=datxlmarg
      xlboff=datxlboff
      iint=datiint
      xmyc=datxmyc
      do i=1,92
        fixabund(i)=databch(i)
      enddo
      do i=1,1000
        isoch(i)=datisoch(i)
        isochfact(i)=datisochfact(i)
      enddo

      print*,tsuji,filmet(1:index(filmet,' ')),
     &       filmol(1:index(filmol,' '))
      print*,noffil,' line lists:'
      do i=1,noffil
          print*,linefil(i)(1:index(linefil(i),' '))
      enddo
      print*,spherical,mihal,taum,ncore,diflog
      print*,detout(1:index(detout,' '))
      print*,inatom(1:index(inatom,' '))
      print*,inmod(1:index(inmod,' '))
      print*,outfil(1:index(outfil,' '))
      do i=1,92
        if (fixabund(i).gt.-90.0) print*,'ab. changed: elt ',
     &      i,fixabund(i)
      enddo
      if (limbdark) then
cc* reads filter transmission, opens output files, 
cc* 46 contains spectrum in interval and 47 contains center-to-limb profile
cc* filter in common filter. Is interpolated weightlimb
cc        ifiltunit=46
cc        open(ifiltunit,file=filterfil,status='old')
cc        call readfilt(ifiltunit)
cc        close(ifiltunit)
cc        open(47,file=outfil,status='unknown')
        open(46,file=outfil(1:lenstr(outfil))//'.spec',
     &       status='unknown')
      else if (multidump) then
* opacities (kappa cont, sigma cont and kappa line) for MULTI input.
        open(46,file=outfil(1:lenstr(outfil))//'.multi',
     &        status='unknown',form='unformatted')
CCCconvert='big_endian')
      else
* spectrum file in ascii format (2 columns: lambda,flux)
        open(46,file=outfil,status='unknown')
      endif
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
cc      OPEN(UNIT=15,FILE=INSPEC,STATUS='OLD')
cc      open(unit=16,file=inabun,status='old')
cc      OPEN(UNIT=27,STATUS='SCRATCH',FORM='FORMATTED')
ccc      OPEN(UNIT=20,FILE=OUTFIL,STATUS='UNKNOWN',FORM='UNFORMATTED')
cc      OPEN(UNIT=20,FILE=OUTFIL,STATUS='UNKNOWN',FORM='UNFORMATTED',
cc     &     RECL=412)
      OPEN(UNIT=23,STATUS='SCRATCH',FORM='FORMATTED')
*
* Initiate
*
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

cc      OPEN(UNIT=11,FILE=INATOM,STATUS='OLD')
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
            print*,'bsyn: error in isotope specification',isoch(i)
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
*
cccc      call cstrip(alunit,27)
cccc      READ(27,101) IP,EPS
      IP=0
      EPS=0.0001
      NMY=6
      iweak=0
*
* limbdarkening?
      if (limbdark) then
cc        xl1=filtlam(1)
cc        xl2=filtlam(ifilt)
        iint=1
      endif
!      if (iint.gt.0) then
!      OPEN(UNIT=66,FILE='sphlimb',STATUS='UNKNOWN',FORM='UNFORMATTED')
!      endif
      if (hydrovelo) then
* this allows line shifts from the first lambda to the last of the list 
* (at least for v/c*lambda <5A.
        xlmarg=5.0
        if (iint.gt.0.or.limbdark) then
          print*
          print*,' TEST VERSION WITH VELOCITY FIELD'
          print*,' CENTRAL INTENSITY MAYBE OK, BUT NOT FLUX !!'
ccc          print*,' ERROR! cannot compute intensities and limbdarkening'
ccc          print*,' ERROR! with velocity shifts. Not yet...'
          print*
ccc          stop
        endif
        if (.not.spherical) then
          print*
          print*,' ERROR! cannot compute PP flux with velocities'
          print*,' ERROR! Use spherical'
          print*
        endif
      endif
      IF(ABS(XLMARG).LT.1.E-6) XLMARG=5.0000
      XLM=(XL1+XL2)/2.
      XL1L=XL1-XLMARG
      XL2R=XL2+XLMARG
      maxlam=int((xl2-xl1)/del)+1
      if (maxlam.gt.lpoint) stop 'bsyn: too many wavelengths'
      do 400 j=1,min(maxlam+100,lpoint)
* concerning this min() see readmo, set up of the continuum opacities
       xlambda(J)=XL1+FLOAT(J-1)*DEL
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
      IF(EPS.LE.0.0) EPS=0.0001
      XITE5=0.0
      XLBOFF=0.0
cc      READ(27,102) XITE5
cc      READ(27,109) XLBOFF
cc      if (xlboff.lt.1.e-6) xlboff=0.0
*
* XITE5  is depth independent microturbulence parameter in
* km/s. If this is zero or unspecified, XI is taken from
* model file, read by subr. READMO.
* XLBOFF is the wavelength offset for each line (Angstroms)
*
* Next read model atm. and continuous abs. coeffs. at XLM
* readmo also interpolates the continuous opacities to all the
* xlambda(1:maxlam) . 
* They are put into absocont and absoscont. the wavelengths that
* were used in babsma for the calculation of the cont. opacities
* are put into xlp.
*
      CALL READMO(XLM,X,S)
      if ( (rr(1).eq.0. .or. abs((rr(1)-rr(2))/rr(1)).lt.1.e-8) 
     &                 .and.spherical) then
        print*,
     &   'WARNING! transfer should be treated PP, as the model is PP!'
        stop ' COMPUTATION HALTED!'
**        print*,'WARNING! Swaping to PP transfer!'
**        spherical=.false.
      endif
*
* initiate total extinction
      do j=1,maxlam
        do K=1,ntau
          absos(k,j)=absoscont(k,j)
          abso(k,j)=absocont(k,j)
        enddo
      enddo
*
* X  is kappa/stndop,
* S  is sigma/stndop
*    at each depth for a series of nlcont lambdas (xlmb).
* the complete description of the continuum absorption is in 
* absocont (pure absorption) and absoscont (scattering). They are also 
* divided by stndop. They are computed in readmo.
* stndop is called ross in this routine and is not ross.
*
      WRITE(6,237) MCODE(1:lenstr(mcode))
      IF(IP.EQ.0) WRITE(7,303)
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
      print*,'Bsyn; k, pgmod, pg_calc, romod, ro_calc'

c      skiprelim=.false.
c      print*,'calling eqmol for intializing'
c      print*,t(20),pg(20),pe(20)
c      call eqmol_pe(t(20),pg(20),pgpg,pe(20),
c     &      1.,1.,k,niter,skiprelim)

* Test Plez 11-May-2018
      do 43 k=ntau,1,-1
*      do 43 k=1,ntau
* end of test

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
c        print*,'soon opening file :' , inline(1:index(inline,' '))
        OPEN(UNIT=13,FILE=INLINE,STATUS='old')
cc        print*,'opened file '
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
        print*,species
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
              print*,'H2O list from Barber et al. Using their partf.'
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
          lele='H '
          print 1234,species,lele,iel,ion,(isotope(nn),nn=1,natom)
          print*, 'nlines ', nline
cc          call Hlineadd(lunit,nline,xlboff)
          call hydropac(lunit,xlboff)
          goto 9874
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
        print*,'after getlele: ',iel,lele
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
        endif
* With Tsuji equilibrium, abund is not needed.
*
* Test for molecular list format
* allows backward compatibility for pre-v14.1 format molecular line lists
        if (iel.gt.92) then
          read(lunit,'(a)') oneline
          backspace(lunit)
          read(oneline,*,err=11,end=11) xlb,chie,gfelog,fdamp,gu,raddmp,
     &                levlo,levup
          newformat=.true.
          goto 12
11        newformat=.false.
12        continue 
        else
          newformat=.true.
        endif
          
* Start wavelength loop
*
      ILINE=0
      NALLIN=NALLIN+NLINE
*
* NLINE is the number of lines of the element IEL.
*
***************************************************************
*
      iatouca=0
*
***************************************************************
*
* Big jump to 64 from far below, line loop.
*
   64 CONTINUE
      DGFE=0.
*
* Loop for reading the lines
*
   50 CONTINUE
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
          read(lunit,*) xlb,chie,gfelog,fdamp,gu,raddmp,levlo,levup
          !print*,'SBC*********** ',xlb,chie,gfelog,fdamp,gu,raddmp
          !print*,'SBC----------- ',levlo,levup
        else
* allows backward compatibility for older format molecular line lists
          read(lunit,*) xlb,chie,gfelog,fdamp,gu,raddmp
        endif
*
      if (xlb.lt.xl1l.or.xlb.gt.xl2r) then
        NALLIN=NALLIN-1
        NLINE=NLINE-1
        NREJCT=NREJCT+1
        GOTO 50
      endif
        ! SBC:
        if (nsegments_lambda_min.gt.0) then
            not_discard = 0
            do i = 1,nsegments_lambda_min 
              xl1l_tmp = segments_lambda_min(i)-XLMARG
              xl2r_tmp = segments_lambda_max(i)+XLMARG
              if (xlb.gt.xl1l_tmp.and.xlb.lt.xl2r_tmp) then
                   not_discard = 1
                   exit
              endif
            enddo
            if (not_discard.eq.0) then
                NALLIN=NALLIN-1
                NLINE=NLINE-1
                NREJCT=NREJCT+1
                GOTO 50
            endif
        endif
      f=10**(gfelog)
      xlb=xlb+xlboff
*
* Start line calculations
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
              stop 'Bsyn: Problem with C2 isotopic mix!!'
            endif
* For TiO line scattering **********************************
cc          else if (lele.eq.'TiO                 ') then
cc            scattfrac=1.0
cc            absfrac=1.0-scattfrac
cc            if (scattfrac.gt.0.) then
cc              print*
cc              print*,' WARNING!!!!! ', scattfrac,' of the line opacity',
cc     &           ' counted as scattering!!!!!'
cc              print*
cc            endif
* For TiO line scattering **********************************
          endif
          if (molindex.eq.0) then 
            print*,'bsyn: molecular species not implemented in atomda',
     &           lele
            stop
          endif
          do k=1,ntau
            ntot(k)=partryck(k,molindex)/boltz/t(k)/ro(k)*symmfactor
            ntt(k)=partryck(k,molindex)/boltz/t(k)/ro(k)*symmfactor
          enddo
          call partffordepth(ntau,t,lele,fpartition)
        else 
* For K I line scattering **********************************
          if (iel.eq.19.and.ion.eq.1) then
            scattfrac=0.0
            absfrac=1.0-scattfrac
          else
            scattfrac=0.0
            absfrac=1.0-scattfrac
          endif
          if (scattfrac.gt.0.) then
            print*, 'Element ',iel,' Ion',ion
            print*,' WARNING!!!!! ', scattfrac,' of the line opacity',
     &         ' counted as scattering!!!!!'
            print*
          endif
* For K I line scattering **********************************
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

* we do not use this anymore. It could be reactivated easily,
* to compute tau-scales with velocity.
* We now compute line profiles with velocity-fields in bsynb.f
* BPz 08/08-2001
*
* prepare line shift vs. depth, from velocity in model (in cm/s)
* velocity in model should be positive outwards.
* If xlb_vshifted is the wavelength in the observer's frame at which
* the line position is shifted [(lambda_0-lambda)/lambda_0=v/c], 
* xlb_vshifted=xlb*lshift, with lshift calculated here:

ccc          lshift(k)=1.d0-velocity(k)/2.99792458d10
          lshift(k)=1.d0

          do i=1,natom
            ntot(k)=ntot(k)*isotopfrac(atom(i),isotope(i))
            ntt(k)=ntt(k)*isotopfrac(atom(i),isotope(i))
            if (ntot(k).eq.0.) then
              print*,'Bsyn. WARNING!, ntot=0 for species: ',lele
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
         PRINT *,'***************************************************'
         PRINT *,'Error in BSYN, ION.GT.3 which has not been foreseen'
         PRINT *,'***************************************************'
         STOP '***** BSYN *****'
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
cc      DLAMB0=DOPPLC*DNUD(JLEV0)*XL**2/C
*
* Constants
*
      CALF=constant*F
*
cc      IF(FDAMP.GT.0..AND.IP.GE.1) WRITE(7,267) GAMRAD,ALOGC6
cc      IF(IP.GE.1) WRITE(7,265) LELE,ABUL,ABUND,NTAU
      do 111 j=1,ntau
       plez(j)=n(j)*stim(j)/dnud(j)/ross(j)
       xlb_vshifted(j)=xlb*lshift(j)
111   continue
*
* Start wavelength loop for this line
*
ccc      print*,xlb
      zap=(xlb-xl1)/del
      IF (zap.le.0.) THEN
* treatment of the lines lying between xl1l and xl1
       do k=1,ntau
          contop(k)=absocont(k,1)
       enddo 
       do 222 i=1,maxlam
         xkmax=0.
         do 333 j=1,ntau
          vt=(xlambda(i)-xlb_vshifted(j))*1.d-8
          vt=c*vt/xl**2
           v=vt/dnud(j)
cc           CALL VOIGT(A(j),V,HVOIGT)
           hvoigt=newvoigt(a(j),v)
           l=calf*hvoigt*plez(j)
           ABSO(j,i)=ABSO(j,i)+l*absfrac
           ABSOS(j,i)=ABSOS(j,i)+l*scattfrac
* we compare the line absorption to the continuum x at the rightmost
* lambda of the interval to set the limit of inclusion of this line 
* kappa. We take this continuum x just for convenience.
           xkmax=max(xkmax,l/contop(j))
333      continue
         if (xkmax.lt.eps) goto 255
222    continue
*
      ELSE
*
* now lines lying in the [xl1,xl2] interval
* lindex in xlambda of the closest wavelength > to the wavelength of 
* the line
      lindex=int(zap)+2
      lindex=min(maxlam+1,lindex)
      do k=1,ntau
         contop(k)=absocont(k,min(maxlam,lindex))
      enddo 
      iii=0
      do 2 i=lindex,maxlam
       xkmax=0.
       iii=iii+1
       do 3 j=1,ntau
        vt=(xlambda(i)-xlb_vshifted(j))*1.d-8
        vt=c*vt/xl**2
        v=vt/dnud(j)
cc        CALL VOIGT(A(j),V,HVOIGT)
        hvoigt=newvoigt(a(j),v)
        l=calf*hvoigt*plez(j)
        ABSO(j,i)=ABSO(j,i)+l*absfrac
        ABSOS(j,i)=ABSOS(j,i)+l*scattfrac
        xkmax=max(xkmax,l/contop(j))
*
ccc        if(iprint.eq.1.and.j.eq.1) then
ccc         print*,xlambda(i),hvoigt,xkmax,a(1),v
ccc        endif
*
3      continue
       if (xkmax.lt.eps) goto 15
2     continue
* we get out without being under eps. 2possibilities:
* 1) the line lies towards the end of the wavelength array
* 2) the line encompasses the whole array
      do 422 i=lindex-1,1,-1
       xkmax=0.
       do 433 j=1,ntau
        vt=(xlb_vshifted(j)-xlambda(i))*1.d-8
        vt=c*vt/xl**2
        v=vt/dnud(j)
cc        CALL VOIGT(A(j),V,HVOIGT)
        hvoigt=newvoigt(a(j),v)
        l=calf*hvoigt*plez(j)
        ABSO(j,i)=ABSO(j,i)+l*absfrac
        ABSOS(j,i)=ABSOS(j,i)+l*scattfrac
        xkmax=max(xkmax,l/contop(j))
*
ccc        if(iprint.eq.1.and.j.eq.1) then
ccc         print*,xlambda(i),hvoigt,xkmax,a(1),v
ccc        endif
*
433    continue
       if (xkmax.lt.eps) goto 255
422   continue
      goto 255
* here is the normal continuation
15    CONTINUE
* and now the other side of the profile
*
********************** check n(contributing lines) ******************
*
      if (iii.gt.1) iatouca=iatouca+1
*
*********************************************************************
      istart=max(1,lindex-iii-2)
ccc      print*,istart,lindex-1
      do 22 i=istart,lindex-1
       do 33 j=1,ntau
        vt=(xlb_vshifted(j)-xlambda(i))*1.d-8
        vt=c*vt/xl**2
        v=vt/dnud(j)
cc        CALL VOIGT(A(j),V,HVOIGT)
        hvoigt=newvoigt(a(j),v)
        l=calf*hvoigt*plez(j)
        ABSO(j,i)=ABSO(j,i)+l*absfrac
        ABSOS(j,i)=ABSOS(j,i)+l*scattfrac
*
ccc        if(iprint.eq.1.and.j.eq.1) then
ccc         print*,xlambda(i),hvoigt,xkmax,a(1),v
ccc        endif
*
33     continue
22    continue
*
      ENDIF
*
255   continue
*
* End of line calculation
*
      ILINE=ILINE+1
*
      iannonce=mod(iline,30000)
      if (iannonce.eq.0) then
        call clock
        print*,iline,' lines done'
      endif
      GOTO 64
*
* End model loop
*
* end of line lists loop
******************************************************
*
      print*,'nombre de raies contribuant (environ) : ',iatouca,
     &       ' n. total : ', iline
*
******************************************************
9874  close(lunit)
98    continue
* 
      WRITE(6,214) NREJCT,XL1L,XL2R,XLM
 214  FORMAT(1X,I8,' LINES WERE REJECTED, ONLY LINES BETWEEN',F10.3,
     &       ' AND',F10.3,' A CONSIDERED',/,'  CENTRAL WAVELENGTH=',
     &       F10.3,' A')
*
      call clock
      print*, 'now solve'
      print*,nallin,' lines in the interval'
ccc      print*,icount1,' slow voigt',icount2,' nordlund',
ccc     & icount3,' quick',
ccc     &   icount4,' look-up table'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c        print*,'CHECK kappa !!'
c        do j=1,maxlam
c          if (abs(xlambda(j)-5240.41d0).lt.1.d-4.or.
c     &        abs(xlambda(j)-5242.49d0).lt.1.d-4.or.
c     &        abs(xlambda(j)-5241.62d0).lt.1.d-4)  then
c            print*,'lambda=',xlambda(j)
c            print*,'  K  log(chi) log(kappa) log(sigma)'
c            do k=1,ntau
c              print*,k,T(k),log10(ro(k)),
c     &          log10((abso(k,j)+absos(k,j))*ross(k)),
c     &          log10(abso(k,j)*ross(k)),
c     &          log10(absos(k,j)*ross(k))
c            enddo
c          endif
c        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if (.not.multidump) then
        if (spherical) then
          CALL BSYNB(NALLIN)
        else
          CALL BSYNBplatt(NALLIN)
        endif
      else
* dump opacities for MULTI input. 20 juin 1994.

        do ieliel=1,83
          if (ieliel.ne.61) then
            absave(ieliel) = abund(ieliel)
          endif
        enddo

        write(46) '''* Output from TurboCanary'''
        write(46) '''*'''
        write(46) '''* Model atmosphere'''
        write(46) '''* Linelists'''
        write(46) '''* Atomic abundances '''
        write(46) 
     &  '''* Ntau, Nlam, lambda start, lambda end, delta lambda'''
        write(46) '''* Lambda for tau-scale'''
        write(46) '''* Tau-scale and microturbulence'''
        write(46) '''*'''
        write(46) '''* Continuous absorption coefficient'''
        write(46) '''* ( (Kcont(k,l), k=1,ntau), l=1,nlam)'''
        write(46) '''* Line absorption coefficient same format'''
        write(46) 
     &  '''* Continuous scattering coefficient same format'''
        write(46) '''*'''

        write(46) '''* ',inmod(1:index(inmod,' ')),''''
        do i=1,noffil
          filprint=linefil(i)
          write(46) '''* ',filprint(1:index(filprint,' ')),''''
        enddo
        write(46) (ieliel,log10(absave(ieliel))+12.,ieliel=1,60),
     &              (ieliel,log10(absave(ieliel))+12.,ieliel=62,83)
        write(46) ntau,maxlam,xl1,xl2,del
        write(46) xls
        write(46) (tau(k),k=1,ntau),(xi(k),k=1,ntau)
        write(46) ((absocont(k,j)*ross(k),k=1,ntau),j=1,maxlam)
        write(46) (((abso(k,j)-absocont(k,j))*ross(k),k=1,ntau),
     &                                                  j=1,maxlam)
        write(46) ((absos(k,j)*ross(k),k=1,ntau),j=1,maxlam)
      endif
*
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
  236 FORMAT(' ***STOP IN BSYN***.NTOT(',I2,').LT.0.0, IEL=',
     & I3,/,' ION=',I3,' ABUND=',E10.3,' ILINE=',I4,' XLB=',F9.2)
  237 FORMAT(' MODEL IDENTIFICATION=',A,'; MAIN PROGRAMME BSYN')
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
  990 FORMAT(' INPUT PARAMETERS:'/' IEL=',I3,
     & ' ION=',I1,' NMY=',I1,' FDAMP=',F3.1,' IINT=',I1,
     & ' IMY=',I1,' MA=',F6.2,' ABUND=',F4.2,' CHI=',F5.2,' XITE=',
     & 1PE8.2)
*
      END
