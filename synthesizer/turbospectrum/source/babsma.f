      PROGRAM BABSMA
*
*-----------------------------------------------------------------------
*
* THIS IS A MAIN PROGRAM FOR MAKING CONTINUOUS ABSORPTION COEFF.
* FOR THE CANARY STAFF. MADE BY BG IN JANUARY 1982.
*
* Export version  1988-03-24  ********* Olof Morell *** Uppsala
*
*-----------------------------------------------------------------------
*
      INCLUDE 'spectrum.inc'
      INCLUDE 'tsuji.par'
*
* the dimension 200 corresponds to 20*10 (cf. xl)
*
      parameter (maxim=1000)
      
      doubleprecision rrr(ndp)

      DIMENSION TAU(NDP),PGL(NDP),XI(NDP),X(20*numbset),
     &          S(20*numbset),TAUS(NDP),taur(ndp)
      DIMENSION XLP(20*numbset),RHO(NDP),DUMDUM(NDP),RR(NDP),drr(ndp)
      real dumdumm(ndp)
      real Teff,metallicity,kboltz,mh
      dimension rhobow(ndp),coldens(ndp),xe(ndp),ghoefner(ndp)
      dimension comparison(30),iidum(16),xlr(20)
      real abskk(ndp),spridd(ndp),pgk(ndp)
      character*8 abname,source
      character*9 key
      character*1 bla
      character*80 firstline
      dimension abname(mkomp),source(mkomp)
      LOGICAL MRXF,XIFIX
      character*50 blabla
      CHARACTER*1024 DETOUT,OUTMOD,outmod2
      CHARACTER*50 MOCODE,marcsformat
*
* from bsyn; this ensures that continuum opacities are computed for all
* wavelengths of the bsyn calculation. Babsma must be called before each 
* bsyn call.
      integer    spunit
      parameter (spunit=44)
      doubleprecision xl1,xl2,del,xlmarg,xl1l,xl2r
      common/pieces/ xl1,xl2,del,eeps,nmy,nlbldu,iint,xmyc,iweak
      doubleprecision xlambda(lpoint)
cc      common/large/ xlambda,maxlam,ABSO(NDP,lpoint),
cc     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)

********
      COMMON/CFIL/   IRESET(numbset),ISLASK,IREAT
      COMMON/UTPUT/  IREAD,IWRIT
      COMMON/CXLSET/ NSET,NL(numbset),XL(20,numbset)
      COMMON/COUTR/  NTO,NTPO(10)
      real T(ndp),pe(ndp),ro
      common/CA5/ ab(mkomp), fakt(mkomp),ppe(ndp),tt(ndp),xla(20),
     &            xla3(20),ro,sumabs,sumsca,viktr,iset,nlb
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(30)
* common molecules contains the molecular pressures from a marcs model.
* Do not use! May be used for check if you wish.
      common/molecules/nmol,pressure(ndp,50)
* common for partial pressures
      logical tsuswitch,tsuji,changeab,exist
      doubleprecision parptsuji
      common /tsuji/ tsuji,tsuswitch,nattsuji,nmotsuji,
     &               parptsuji(maxim+400)
      character*128 filmet,filmol,fileab
      common /filetsuji/ filmet,filmol
      character*20 nametryck(maxmol)
      doubleprecision partryck,xmettryck,xiontryck
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),nametryck
      common/rhotsu/rhorho,xmytsuji,ejontsuji
      common/chabu/ changeab,fileab
*
      common/abundch/abch(100)
      real xmass(maxim+400), atmass(100)
      doubleprecision molweight,ndensity
      common/density/ndensity,molweight,xmass,atmass
*
      real kaprefmass(ndp),rhox(ndp),first,tottau
      real datscattfrac
      character*256 inspec

* input from input.f

      parameter (maxfil=100)
      integer datnoffil,datncore,datmmaxfil,datmihal
      character*128 datfilmet,datfilmol,datfilwavel
      character*256 datlinefil,datdetout,datinatom,datinmod,
     &              datinabun,datinspec,datoutfil,datmongofil,
     &              datfilterfil,datcontinopac,datinpmod
      logical dattsuji,datspherical,datlimbdark,databfind,
     &        datmultidump,datxifix,datmrxf,dathydrovelo,
     &        datpureLTE,pureLTE
      integer datiint
      real    datisoch(1000),datisochfact(1000),datxmyc
      doubleprecision  datxl1,datxl2,datdel,datxlmarg,datxlboff
      common/inputdata/datmmaxfil,dattsuji,datfilmet,datfilmol,
     &                 datnoffil,datlinefil(maxfil),datspherical,
     &                 datmihal,dattaum,datncore,
     &                 datdiflog,datdetout,datinatom,
     &                 datinmod,datinabun,datinspec,datoutfil,
     &                 datmongofil,databch(100),datlimbdark,
     &                 datfilterfil,datoverall,databfind,datmultidump,
     &                 datisoch,datisochfact,dathelium,datalpha,
     &                 datrabund,datsabund,datxifix,datxic,datmrxf,
     &                 datinpmod,datcontinopac,datfilwavel,dathydrovelo,
     &                 datxl1,datxl2,datdel,datxlmarg,datxlboff,
     &                 datiint,datxmyc,datscattfrac,datpureLTE
*
      real amass(92,0:250),abund(92),fixabund(92),
     &         isotopfrac(92,0:250)
      real overall,alpha,helium,rabund,sabund
      character*2 aname(92)

      common/refabundances/ abund,amass,aname,isotopfrac

      logical hydrovelo
      real velocity
      common/velo/velocity(ndp),hydrovelo

      real*8 vacair
      real*8 edge(6)
      data edge/ 912.00, 3648.04, 8209.26, 14591.99, 22800.22, 32923.98/ ! vacuum edges from hydropac.f
cc      data edge/ 912.00, 3647.98, 8209.26, 14591.99, 22800.22, 32923.98/ ! vacuum edges from marcs.f
cc      data edge /912., 3647., 8207., 14588., 22794., 32915./
cc* air edges from continuous opacity file.
      integer version
      data version /191/
      data kboltz / 1.3806e-16 /
      data mh / 1.660e-24 /

      print*
      print*,'***********************'
      print10,version*0.1
10    format(' * BABSMA version ',f4.1,' *')
      print*,'***********************'
      print*

      gravl=4.44
      intryc=0
      scale=0.
* 
c bsyn uses air wavelengths:
      do k=2,6
        edge(k)=vacair(edge(k))
        print*,'Babsma.f HI k,edge(k)',k,edge(k)
      enddo

      datmmaxfil=maxfil
      tsuji=.false.
      tsuswitch=.false.
      do k=1,100
        abch(k)=-99.9
      enddo
       
      do k=1,ndp
        velocity(k)=0.
      enddo

      IREAD=11
      IWRIT=7
      IREAT=11
      ISLASK=15
*
      DO 1 K=1,numbset
        IRESET(K)=16
    1 CONTINUE  
      IMOD=12
      IMODUT=13
      IMODUT2=42
*
      call input
*
* XL content:
* 20*(numbset-1) wavelengths allowed.
* Starting at set #2. Set #1 contains std lambda for tau-scale (read in model)
* BPz 30/08-1999 new version using same wavelength set as bsyn.
*   This creates too large a contopac file. I take every AA instead. BPz 1 hour later.
*
      xl1=datxl1
      xl2=datxl2
      del=datdel
      xmyc=datxmyc
      iweak=0
      xlmarg=datxlmarg
      iint=datiint
      pureLTE=datpureLTE

      IF(ABS(XLMARG).LT.1.E-6) XLMARG=5.0000
      XLM=(XL1+XL2)/2.
      XL1L=XL1-XLMARG
      XL2R=XL2+XLMARG
************ !!! ******************
      DEL=max(del,1.0d0)
************ !!! ******************
      maxlam=int((xl2-xl1)/del)+1
      if (maxlam.gt.lpoint) stop 'babsma: too many wavelengths'
      jmax=min(maxlam+100,lpoint)
* concerning this min() see readmo, set up of the continuum opacities
      do j=1,jmax
        xlambda(J)=XL1+dble(J-1)*DEL
cc        do k=1,ndp
cc          abso(k,j)=0.0
cc          absos(k,j)=0.0
cc          absocont(k,j)=0.0
cc        enddo
      enddo

* Try to avoid erroneous interpolation of continuum opacity near H I edges
* producing spurious spikes in spectra. BPz 16/05-2007 

      do k=1,6
        if (xlambda(1).lt.edge(k)) then
          if (xlambda(jmax).gt.edge(k)) then
* there is an HI edge inside the interval
            do j=2,jmax-1
c              if (xlambda(j).gt.edge(k)) then
              if (xlambda(j).ge.edge(k)) then
                xlambda(j)=edge(k)
                xlambda(j-1)=edge(k)-min(datdel,0.9d0)
                xlambda(j+1)=edge(k)+min(datdel,0.9d0)
                print*,'Found edge: ',edge(k),xlambda(j-1),xlambda(j),
     &                  xlambda(j+1)
                if (xlambda(j-1).le.xlambda(j-2)) then
                  print*,'edge, lambda cont: ',edge(k),xlambda(j-2),
     &                 xlambda(j-1),xlambda(j)
                  stop 'problem with continuum wavelength'
                endif
                if (j.lt.jmax.and.xlambda(j+1).le.xlambda(j)) then
                  print*,'edge, lambda cont: ',edge(k),xlambda(j-1),
     &                 xlambda(j),xlambda(j+1)
                  stop 'problem with continuum wavelength'
                endif
                goto 111
              endif
            enddo
111         continue
          endif
        endif
      enddo
********************************
      np=0
      nnl=0
      j=0
      do k=2,numbset
* read x sets 'k' of 20 wavelength points 'kk' (AA)
        do kk=1,20
          j=j+1
          xl(kk,k)=xlambda(j)
          nnl=nnl+1
          if (nnl.gt.(numbset-1)*20) then
            stop 'babsma, too many wavelengths for continuum. ERROR!'
          endif
          if (j.eq.jmax) then
             goto 248
          endif
        enddo
      enddo
  248 continue
      if (j.lt.jmax) stop 
     &    'babsma: too few wavelengths for continuum. Increase numbset!'
      np=int(nnl/20)
      do k=2,np+1
        nl(k)=20
      enddo
      nlast=nnl-20*np
      if (nlast.gt.0) then
        np=np+1
        nl(np+1)=nlast
      endif
      print *,nnl,' wavelengths for continuum between: ',xl(1,2),
     &        ' and ',xl(nl(np+1),np+1)
*
      hydrovelo=dathydrovelo
      detout=datdetout
      outmod=datinmod
      outmod2=outmod(1:lenstr(outmod))//'.mod'
      tsuji=dattsuji
      mrxf=datmrxf
      xifix=datxifix
      xic=datxic
      filmol=datfilmol
      overall=datoverall
      alpha=datalpha
      rabund=datrabund
      sabund=datsabund
      helium=dathelium
      do i=1,92
        fixabund(i)=databch(i)
      enddo

      do i=1,92
        if (fixabund(i).gt.-90.0) print*,'ab. changed: elt ',
     &      i,fixabund(i)
      enddo
*
* get abundances and scale them by overall. Then 
* if appropriate replace by fixabund.
* for molecular equilibrium calculation
*
      call makeabund(overall,alpha,helium,rabund,sabund,fixabund,
     &                  abund,amass,aname,isotopfrac)

      print*,'metallicity changed by ',overall,' dex'
*
*
* JON DATA AND ABSKOEFF DATA TO BE READ FROM UNIT 11
* PRINTOUT ON UNIT 7
* PRELIMINARY FILES ON UNIT 15 AND 16
* MODEL READ FROM UNIT 12, NEW MODEL DATA WRITTEN ON UNIT 13
*
      OPEN(UNIT=IWRIT,FILE=DETOUT,STATUS='UNKNOWN')
c,recl=4*2*200)
      OPEN(UNIT=IREAT,FILE=datcontinopac,STATUS='OLD')
      OPEN(UNIT=IMODUT,FILE=OUTMOD,STATUS='UNKNOWN')
      OPEN(UNIT=IMODUT2,FILE=outmod2,STATUS='UNKNOWN')
c,recl=4*2*200)
      OPEN(UNIT=ISLASK,STATUS='SCRATCH',FORM='UNFORMATTED')
c,RECL=8192)
      OPEN(UNIT=16,STATUS='SCRATCH',FORM='UNFORMATTED')
c,RECL=16384)
ccc      OPEN(UNIT=34,FILE='pression.dat',STATUS='UNKNOWN')

**** input parameter read. Now read model.

      if (MRXF) then

* check if MARCS model is a binary or ascii file.
* If it is an ascii file, we assume it is one of the 
* models of the new UPPSALA grid (2004-2005)
*   BPz 12/01-2005
*
        inquire(file=datinpmod,exist=exist)
        if (.not.exist) then
          print*,'ERROR! Model file does not exist'
          stop
        endif
        marcsformat='unknown'
* try if model file is binary
        open(unit=imod,file=datinpmod,form='unformatted')
        print*,'opened unformatted'
        read(imod,err=3003,iostat=ios) dum
        print*,'read record 1'
        print*,'iostat=',ios
        if (ios.ne.0) then
           goto 3003
        endif
        read(imod,err=3003,iostat=ios) idum,dum,dum
        print*,'read record 2'
        print*,'iostat=',ios
        if (ios.ne.0) then
           goto 3003
        endif
        marcsformat='binary'
        goto 3002
* try if model file is ascii
3003    close(imod)
        open(unit=imod,file=datinpmod,form='formatted')
        read(imod,'(a)',err=3001) mocode
        print*,'mocode',mocode
        marcsformat='ascii'
        goto 3002
3001    print*,'Error! Unable to determine MARCS model format!'
        stop
3002    print*,'MARCS model format is ',marcsformat(1:10)
        close(imod)
      
        if (marcsformat(1:6).eq.'binary') then
*
* Marcs model reading. if any question see oslistmo.f
          open(unit=imod,file=datinpmod,status='old',form='unformatted')
CCCC     &       convert='big_endian')
C     &     ,RECL=152600)
*
          read(imod) dum
          read(imod) idum,dum,dum
          read(imod) dum,dum,dum,dum,dum,dum,dum,jdum,jdum,jdum,jdum,
     &           jdum,jdum,jdum,jdum,jdum,jdum,nel,
     &           (dumdum(i),i=1,nel),dum,dum,dum,dum,dum
*
* Some MARCS binary models (e.g. RCrB from KE) may be with RR in 
* double precision. We try it. If it does not work we try again in
* single precision
*      BPz 17/09-2007
*
****** THIS DOES NOT WORK ANYMORE WITH THE NEW INTEL COMPILER????  BPz 26/03-2018
* Error condition does not happen with single precision RR() models
******** REMOVED !!
*          read(imod,err=999) ntau,jdum,dum,ddum,dddum,(rrr(i),i=1,ntau)
*          do i=1,ntau
*            rr(i)=sngl(rrr(i))
*          enddo
*          goto 998
*999       backspace(imod)
*********** END OF REMOVED

          read(imod) ntau,jdum,dum,ddum,dddum,(rr(i),i=1,ntau)
998       continue
          if (ntau.gt.ndp) then
            print*,' ndp = ',ndp,'       ntau = ', ntau
            stop 'ndp too small!'
          endif
cccc          read(imod) ntau,(dumdumm(i),i=1,ntau),(dumdum(i),i=1,ntau)
cccc          print*,'reading ntau again ',ntau
          read(imod) 
          pi=3.14159
          do k=1,ntau
* tau is tau (xls). The 2nd dum is tauross.
* nmol+1 is TiO
            read(imod) dum,taur(k),tau(k),dum,t(k),pe(k),pgl(k),dum,dum,
     &       dum,rho(k),dum,dum,dum,dum,dum,dum,dum,dum,dum,nmol,
     &       (pressure(k,i),i=1,nmol+1)
          enddo
          read(imod)(iidum(i),i=1,nel),nnlp,(xlr(i),i=1,nnlp),
     &     jdum,idum,idum,(abname(kp),source(kp),kp=1,jdum)
*
          xls=xlr(nnlp)
          Print*, ' Lambda standard = ',xls
*
          if(.not.xifix) then
***            read 105,(xi(k),k=1,ntau)
            print*,' xifix=.false., but where is microturbulent',
     &              'velocity given??'
            stop 'error'
          else
            do k=1,ntau
              xi(k)=xic
            enddo
          endif
        else if (marcsformat(1:5).eq.'ascii') then
* MARCS model in ascii format
          open(unit=imod,file=datinpmod,status='old')
          read(imod,'(a)') mocode
          print*,mocode
          if (mocode(1:1).eq.'s'.or.mocode(1:1).eq.'p') then
            print*,'This model seems to be aan ascii MARCS model'
          else
            print*,' This model may not be a MARCS model!'
          endif

          read(imod,'(a)') blabla
          do while (blabla(14:19).ne.'Radius'.and.
     &              blabla(19:24).ne.'radius')
            read(imod,'(a)') blabla
          enddo
          read(blabla,*) radius
          if (radius.ge.2.) then
            print*,' this model is SPHERICALLY SYMMETRIC'
          else if (radius.lt.2.) then
            print*,' this model is PLANE PARALLEL'
          endif

          do while (blabla.ne.'Model structure')
            read(imod,'(a)') blabla
          enddo
          backspace(imod)
          backspace(imod)
          read(imod,*) ntau
          if (ntau.gt.ndp) then
            print*,'ntau = ',ntau,' greater than ndp = ',ndp
            print*,'Increase NDP'
            stop
          endif
          read(imod,*)
          read(imod,*)
          do k=1,ntau
* the first dum is tauross, the second is Prad, the third is Pturb
            read(imod,*) idum,dum,tau(k),rr(k),T(k),
     &                   pe(k),pgl(k),dum,dum
            tau(k)=10.**tau(k)
            rr(k)=radius-rr(k)
          enddo
          xls=5000.
          print*,'Beware !! Lambda standard assumed to be 5000 A'
          print*,' Make sure it is consistent with the model'
          if(.not.xifix) then
            print*,' xifix=.false., but where is microturbulent',
     &              'velocity given??'
            stop 'error'
          else
            do k=1,ntau
              xi(k)=xic
            enddo
          endif
        else
          print*,'ERROR in MARCS format model'
          stop
        endif
      else
        print*,' This is an ascii model '
        open(unit=imod,file=datinpmod,status='old')
        read(imod,*,err=765) mocode,ntau,xls,gravl,intryc,scale
        if (ntau.gt.ndp) then
            print*,' ndp = ',ndp,'       ntau = ', ntau
            stop 'ndp too small!'
        endif
        print*, mocode,ntau,xls,gravl,intryc,scale
*
* MOCODE IS A 50 LETTERs IDENTIFICATOR FOR THE MODEL
* NTAU IS THE NUMBER OF DEPTH POINTS IN THE MODEL
* DEPTH SCALE AT WAVELENGTH XLS (IN AANGSTROEM)
* GRAVL IS LOGARITHMIC GRAVITY
* INTRYC GT O IF PRESSURE INTEGRATION IS WANTED
*
        goto 764
765     continue
* We try if this might be an ATLAS model
        backspace(imod)
        read(imod,'(a4)',err=762) mocode
        if (mocode.eq.'TEFF') then
          do while (mocode.ne.'READ')
            read(imod,'(a4)') mocode
          enddo
          backspace(imod)
          read(imod,763) ntau
763       format(10x,i3)
          if (ntau.gt.ndp) then
            print*,' ntau = ',ntau,' > ndp = ',ndp
            stop 'increase ndp!'
          endif
        endif
        mocode='KURUCZ'
        xls=5000.
        gravl=0.0
        intryc=0
        scale=0.
        goto 764
762     stop 'COULD NOT READ MODEL ATMOSPHERE FILE !!!'
764     continue

***********************************************
        if (mocode(1:3).eq.'sph') then
          do k=1,ntau
            read(imod,*) tau(k),t(k),pe(k),pgl(k),xi(k),rr(k)
            tau(k)=10.**tau(k)
            t(k)=(1.000+scale)*t(k)
            pe(k)=10.**pe(k)
            pgl(k)=10.**pgl(k)
          enddo

***********************************************
        else if (mocode(1:4).eq.'bowe'.or.mocode(1:4).eq.'BOWE') then
          mocode(1:4)='bowe'
          do k=1,12
            read(imod,*)
          enddo
          do k=ntau,1,-1
* BOwen's models are for increasing radius.
cccLuttermoser            read(imod,1963) iii,rr(k),velocity(k),t(k),rhobow(k),
cccLuttermoser     &                 coldens(k),xe(k)
cccLuttermoser1963        format(i3,6e12.0)
            read(imod,*) iii,rr(k),drr(k),t(k),rhobow(k),
     &                     velocity(k),xe(k)
cc1963        format(i3,2x,6(e12.0))
* xe is ne/ntot * 0.908. Pending confirmation for Bowen's new models. (1996)
* Bowen's says this xe should not be used, esp. not outside hot regions.
* guess pg and pe
            pgl(k)=rhobow(k)*1.38e-16*T(k)/1.26/1.67e-24
            pe(k)=xe(k)*0.908*pgl(k)
          enddo
***********************************************
        else if (mocode(1:4).eq.'alva') then
* modeles de Rodrigo sans echelle de profondeur optique.
* Initially from R. Alvarez. Models without tau scale.
* These models have T, Pgas, R, and either a depth dependent microturbulence, 
* or a radial velocity field. The velocity field is used for models with a
* wind extension. Calculations of flux spectra do not work in that case in 
* the v12.1 and previous versions of the code. 
          read (imod,'(a)') firstline
          read(firstline,*,err=99,end=99) t(1),pgl(1),rr(1),xi(1)
          print*,'checking firstline: ' ,t(1),pgl(1),rr(1),xi(1)
          backspace(imod)
          if (hydrovelo) then
            print*, 'Reading a model with radial velocity field'
          else
            print*,
     &       ' Reading a model with depth dependent microturbulence'
          endif
          do k=1,ntau
            if (hydrovelo) then
              read(imod,*) t(k),pgl(k),rr(k),velocity(k)
              xi(k)=xic
              print*, k,t(k),pgl(k),rr(k),xi(k),velocity(k)
            else
              read(imod,*) t(k),pgl(k),rr(k),xi(k)
              print*, k,t(k),pgl(k),rr(k),xi(k)
            endif
* check definition of rr. between tau,t,pe points? 10/12-1996
            pe(k)=1.e-4
            if (t(k).lt.2000.) then
* This estimate of Pe is crucial for the convergence of the molecular equilibrium
* in the low temperature regime. In case of non convergence of the equilibrium, 
* one may try to lower Pe through an increase of the exponent. I change this
* exponent from 20 to 30 today. BPz 15/06-2012. 
              pe(k)=1.e-4*pgl(k)*(t(k)/3000.)**30
            endif
          enddo
          goto 98
 99       backspace (imod)
          print*,' Reading a model without microturbulence or ',
     &       'radial velocity field'
          do k=1,ntau
            xifix=.true.
            read(imod,*) t(k),pgl(k),rr(k)
            print*, k,t(k),pgl(k),rr(k)
            pe(k)=1.e-4
            if (t(k).lt.2000.) then
              pe(k)=1.e-4*pgl(k)*(t(k)/3000.)**30
            endif
          enddo
 98       continue
          do k=1,ntau-1
            drr(k)=rr(k)-rr(k+1)
          enddo
          drr(ntau)=drr(ntau-1)
***********************************************
        else if (mocode(1:7).eq.'Stagger') then
* Stagger average <3D> model containing: depth, T, rho, averaged on constant tau surfaces
* added by BPz 16/04-2018
          read(imod,*) Teff,gravl,metallicity,ntau
          do k=1,ntau
            read(imod,*) rr(k),T(k),rhobow(k)
            pe(k)=1.e-10                                                  ! guess value only
            pgl(k)=rhobow(k)*kboltz*t(k)/1.3/mh                           ! guess value only
          enddo
          intryc=0
          scale=0.

***********************************************
        else if (mocode(1:7).eq.'Hoefner') then
* S. hoefner models, with velocity. BPz 04/03-2002
* read header!
          bla='#'
          do while (bla.eq.'#')
            read(imod,'(a)') bla
          enddo
          backspace(imod)
          do k=1,ntau
            read(imod,*) rr(k),rhobow(k),T(k),
     &                   pgl(k),tauhoefner,kappahoefner,tdust,k3hoefner,
     &                   velocity(k),
     &                   tradhoefner
            pe(k)=1.e-10
* test !!!!!!!
            print*
            print*,' WARNING !!!!!!! TEST !!!!!! velocity=5km/s!!!!!!'
            print*
            velocity(k)=5.e5
          enddo
* extrapolate inwards the models that have a tau_max too small.
          kk=1
          k=ntau+kk
          do while ((T(k-1).le.6000.).and.(ntau+kk.le.ndp))
            Tstep=min(400.,(t(ntau)-t(ntau-1))*2.**(float(kk)/2.))
            t(k)=Tstep+t(k-1)
            scalefactor=Tstep/(t(ntau)-t(ntau-1))
            rhobow(k)=(rhobow(ntau)-rhobow(ntau-1))*scalefactor+
     &                  rhobow(k-1)
            pgl(k)=(log(pgl(ntau))-log(pgl(ntau-1)))*scalefactor+
     &                  log(pgl(k-1))
            pgl(k)=exp(pgl(k))
            pe(k)=1.e-10
            velocity(k)=velocity(ntau)
            print*,' k',k,' T ',t(k)
            kk=kk+1
            k=ntau+kk
          enddo
          ntauinput=ntau
          ntau=ntau+kk-1
* the geometrical depths are computed using rho and the hydrostatic equation
*
          print*,'WARNING!!!! model extrapolated inwards',
     &              ' using the hydrostatic approximation!!!'
          do k=2,ntauinput-1
* this computed "hydrostatic" gravity is defined at integer tau-points
* (i.e. where T, P, rho etc are defined)
            ghoefner(k)=(pgl(k+1)-pgl(k-1))/2./(rr(k+1)-rr(k))/rhobow(k)
          enddo
          gr2=0.
* We estimate the gravity from an average of the 4 inner points.
          do k=ntauinput-4,ntauinput-1
            gr2=gr2+ghoefner(k)*((rr(k)+rr(k+1))*0.5)**2
          enddo
          gr2=gr2/4.
          do k=ntauinput+1,ntau
* we solve for the radius, using gr2:
            aequa=(pgl(k)-pgl(k-2))/(gr2*rhobow(k-1))/8.
            bequa=rr(k-1)
            gammaequa=bequa**2 + bequa/aequa
            betaequa=2.*bequa - 1./aequa
            rr(k)=(-betaequa + sqrt(betaequa**2-4.*gammaequa))/2.
          enddo
          do k=ntauinput,ntau-1
            ghoefner(k)=(pgl(k+1)-pgl(k-1))/2./(rr(k+1)-rr(k))/rhobow(k)
          enddo
* check model:
          print*,'Hoefner model. CHECK extrapolation'
          print*,'   R       T       ro        Pg      gstatic*r^2'
          do k=1,ntau-1
            print*,rr(k),t(k),rhobow(k),pgl(k),
     &             ghoefner(k)*((rr(k)+rr(k+1))*0.5)**2
          enddo
          k=ntau
          print*,rr(k),t(k),rhobow(k),pgl(k),
     &             ghoefner(k)*rr(k)**2
*
* SH models have velocity at r points, and rho,T, P in between.
* T, P, rho from one line of input describe conditions between
* the r of that line the r at next line of input. Models
* starting from the outer layers and going inwards.
          do k=1,ntau-1
            drr(k)=rr(k)-rr(k+1)
          enddo
          drr(ntau)=drr(ntau-1)

***********************************************
        else if (mocode(1:6).eq.'KURUCZ') then
*
* Kurucz models. Reading + rinteg taken from moog. 06/04-2001 BPz+ST (Sivarani)
*
* Modified 16/07-2015 by BPz. rinteg moved down, to compute the tau-scale 
* directly at the reference wavelength without using the Rosseland scale.
* 
          do k=1,ntau
            read (imod,*) rhox(k),t(k),pgl(k),pe(k),kaprefmass(k)
            pe(k)=pe(k)*T(k)*1.38054e-16
            print*,'reading: ', k, rhox(k),t(k),pgl(k),pe(k),
     &              kaprefmass(k)
          enddo
        else
          DO 11 K=1,NTAU
            READ(IMOD,*) TAU(K),T(K),PE(K),PGL(K),XI(K)
            TAU(K)=10.**TAU(K)
            T(K)=(1.000+SCALE)*T(K)
            PE(K)=10.**PE(K)
            PGL(K)=10.**PGL(K)
   11     CONTINUE
        endif
        if (xifix) then
          do K=1,NTAU
            XI(K)=XIC
          enddo
        endif
*
      endif
*
* MODELS CHEMICAL COMPOSITION AND VARIOUS DATA FOR IONIZATION
* EQUILIBRIUM ARE READ BY INJON
* A MARCS MODEL IS REWINDED IN INJON  !!!!!!!!
*
      IO=0
      CALL INJON(IO,MRXF)
**********
c      print*,'injon done'
**********
      XL(1,1)=XLS
      NL(1)=1
*
* THE STANDARD WAVELENTH XLS IS IN THE FIRST WAVELENGTH SET.
*
* NP IS THE NUMBER OF WAVELENGTH SETS TO BE READ SUBSEQUENTLY,
* NL THE NUMBER OF WAVLENGTHS IN EACH SUCH SET.
*
      NSET=NP+1
      NLQ=0
      IPP=1
      DO 14 K=1,NP
        NLP=NL(K+1)
        NLQ=NLP+NLQ
        DO 14 I=1,NLP
          XLP(IPP)=XL(I,K+1)
          IPP=IPP+1
   14 CONTINUE
*
* NOW, INITIATE ABSORPTION COEFFICIENT TABLES BY CALLING INABS
*
      CALL INABS(IO)
**********
      print*,'inabs done'
**********
*
* IF YOU WANT LESS PRINTOUT FROM INABS, CHANGE IWRIT TO A DUMMY
* UNIT AT THIS CALL (AND BACK AGAIN AFTER) OR SPECIFY IO TO 0.
*
      NTO=0
      IREADP=IRESET(1)
      REWIND IREADP
*
      GRAV=10.**(GRAVL)
*
* TIME TO START MODEL LOOP WITH CALCULATIONS FOR EACH DEPTH
* OF IONIZATION EQ. AND OF ABS. COEFF.
*
      IF(MRXF) THEN
        WRITE(IMODUT,200) 'MRXF',NTAU,XLS
        WRITE(IWRIT,300) 'MRXF',NTAU,XLS
      ELSE
        WRITE(IMODUT,200) MOCODE(1:lenstr(mocode)),NTAU,XLS
        WRITE(IMODUT2,2000) MOCODE(1:lenstr(mocode)),NTAU,XLS
        WRITE(IWRIT,300) MOCODE(1:lenstr(mocode)),NTAU,XLS
      ENDIF
      WRITE(IMODUT,203) NLQ
      WRITE(IMODUT,204) (XLP(IPP),IPP=1,NLQ)
ccc      IF(NLQ.GT.NDP) PRINT 303,NLQ
ccc      IF(NLQ.GT.NDP) STOP
      IF(INTRYC.GT.0) CALL TRYCK(GRAV,NTAU,TAU,T,PE,PGL)
      NEWT=2
      if (tsuji) tsuswitch=.true.
******
* BPz 16/07-2015
* modified for Kurucz models. In order to compute directly the 
* tau-scale for the reference wavelength xls, we store the reference
* opacity in the kaprefmass array, that originally contained the 
* Kurucz model rosseland opacity.
      if (mocode(1:6).eq.'KURUCZ') then
        CALL ABSKO(NEWT,ntau,T,PE,1,1,ABSKK,SPRIDD)
        do k=1,ntau
          kaprefmass(k)=ABSKk(k)+SPRIDd(k)
        enddo
        newt=1
        first = rhox(1)*kaprefmass(1)
        tottau = rinteg(rhox,kaprefmass,tau,ntau,first)
        tau(1) = first
        print*,'Kurucz model. Computed tau-scale at lambda= ',xls,'A'
        print*,'tau(1)=',tau(1)
        do k=2,ntau
          tau(k) = tau(k-1) + tau(k)
          print*,'tau(',k,')=',tau(k)
        enddo
        do k=1,ntau
cc          kapref(k) = kaprefmass(k)*rho(k)
          print 222, k,log10(tau(k)),tau(k),T(k), log10(pgl(k)),
     &                pgl(k),log10(pe(k))
222       format(i2,x,f5.2,x,1pe11.4,2x,0pf7.1,2x,f6.3,x,1pe11.4,
     &             2x,0pf6.3)
        enddo
      endif

******   Model read, start calculations

* attempt to ease convergence at low T / BPz 15/05-2018
      if (mocode(1:4).eq.'alva') then
        if (t(1).lt.1000.) then
          do k=ntau,1,-1
            pg=pgl(k)
            if (k.ne.ntau) then
* better guess of pe from previous depth point
              pe(k)=pe(k+1)*pgl(k)/pgl(k+1)
            endif
            call jon(t(k),pe(k),1,pg,ro,dum,io,k)
            print888,t(k),pe(k),pg,pgl(k)
888         format('first jon call T, Pe, Pg Pgin ',f6.0,3(1x,1pe9.3))
* must iterate on pe to get right pg in model.
* try better guess input pe:
            pein=pe(k)*pgl(k)/pg
            call pemake(t(k),pein,pgl(k),pe(k))
            call jon(t(k),pe(k),1,pg,ro,dum,io,k)
* pg should be within eps of pgl (cf. pemake).
            print889,t(k),pe(k),pg,pgl(k)
889         format('after jon iter T, Pe, Pg Pgin ',f6.0,3(1x,1pe9.3))
          enddo
        endif
      endif

* end of attempt

      DO 25 K=1,NTAU
        if (mrxf) then
* guess input pg to help eqmol_pe.
          pg=pgl(k)
          CALL JON(T(K),PE(K),1,PG,RO,DUM,IO,k)
**********
c      print*,'jon done'
**********
          print*,'MARCS model. romod= ',rho(k),'rocalc= ',ro
          RO=RHO(K)
        else
          if (mocode(1:4).ne.'BOWE'.and.mocode(1:4).ne.'bowe'.and.
     &        mocode(1:4).ne.'alva'.and.mocode(1:7).ne.'Hoefner'.and.
     &        mocode(1:7).ne.'Stagger') then
* guess input pg to help eqmol_pe.
            pg=pgl(k)
            call jon(T(K),PE(K),1,PG,RO,DUM,IO,k)
            if (abs(pg/pgl(k)-1.).gt.0.03) then
              print*,'WARNING ! calculated gas pressure difference > 3%'
              print*,'model rocalc=',ro,' pgcalc=',pg,
     &               ' pgmod=',pgl(k)
            endif
          else
* for funny models, like miras
            if (k.ne.1) then
* better guess of pe from previous depth point
ccccc	      if (mocode.ne.'bowe'.or.mocode.ne.'BOWE') then
              if (mocode(1:4).eq.'alva') then
                pe(k)=pe(k-1)*pgl(k)/pgl(k-1)
              else
                pe(k)=pe(k-1)*rhobow(k)/rhobow(k-1)
              endif
            endif
            iteration=0
1964        continue
C
* guess input pg to help eqmol_pe.
            pg=pgl(k)
            print*,'babsma, alva, call jon'
            call jon(T(K),pe(k),-1,PG,RO,DUM,IO,k)
            print*,'return from jon, pg=',pg
ccc unecessary ?
C try to improve convergence
            pefirst=pe(k)
            do while (pg.lt.0..and.pe(k).gt.1.e-30)
c not converged in jon
              pe(k)=pe(k)/10.
              print*,'babsma, trying with lower Pe:',pe
              call jon(T(K),pe(k),-1,PG,RO,DUM,IO,k)
            enddo
            if (pg.lt.0.) pe(k)=pefirst
            do while (pg.lt.0..and.pe(k).lt.1000.)
c not converged in jon
              pe(k)=pe(k)*10.
              print*,'babsma, trying with lower Pe:',pe
              call jon(T(K),pe(k),-1,PG,RO,DUM,IO,k)
            enddo
            if (pg.lt.0.) then
              print*,'giving up in babsma'
              stop 'Too many unsuccessful tries.'
            endif
ccc
* must iterate on pe here to get right ro in Bowen, Hoefner or Stagger models.
* pgl/ro==kT/muamu 
            if (mocode(1:7).eq.'Hoefner'.or.
     &          mocode(1:7).eq.'Stagger') then
              ratiobow=rhobow(k)/ro
              print*,'rho_input/rho ',k,ratiobow
              pein=pe(k)
              call pemakero(t(k),pein,rhobow(k),pe(k))
              call jon(T(K),pe(k),1,PG,RO,DUM,IO,k)
              ratiobow=rhobow(k)/ro
              print*,'rho_input/rho ',k,ratiobow
              
              print*,'k  rho_input    ro_calc     pe     ',
     &              'pg_input   pg'
              print*,k,rhobow(k),ro,pe(k),pgl(k),pg

            else if (mocode(1:4).eq.'bowe'.or.
     &               mocode(1:4).eq.'BOWE') then
c              ratiobow=rhobow(k)/ro*molweight/1.26
              ratiobow=rhobow(k)/ro
              print*,'rhobow/ro',k,ratiobow
              pgl(k)=pgl(k)*ratiobow*1.26/molweight
              pe(k)=xe(k)*0.908*pgl(k)
* for Bowen's new models (10/12-1996)
              pein=pe(k)
              call pemakero(t(k),pein,rhobow(k),pe(k))
              CALL JON(T(K),pe(k),1,PG,RO,DUM,IO,k)
              ratiobow=rhobow(k)/ro
              print*,'rhobow/ro',k,ratiobow
              print*,'k     rhobow      ro_calc     pe     ',
     &              'xebow     xe_calc   pg'
1965          print*,k,rhobow(k),ro,pe(k),xe(k),pe(k)/0.908/pgl(k),
     &                 pgl(k)
            else if (mocode(1:4).eq.'alva') then
* must iterate on pe to get right pg in Rodrigo's models.
* try better guess input pe:
              pein=pe(k)*pgl(k)/pg
              print*,'babsma, improving on guess Pe = ',pein
              call pemake(t(k),pein,pgl(k),pe(k))
              print*,'Alva format. Converged'
              print*,'Pe=',pe(k)
              CALL JON(T(K),pe(k),1,PG,RO,DUM,IO,k)
* pg should be within eps of pgl (cf. pemake).
            endif

          ENDIF
        endif
*
* COMPARE PRESSURES, TABULATED/CALCULATED
*
        PCOMP=PGL(K)/PG
        WRITE(IWRIT,301) TAU(K),T(K),PE(K),PG,PCOMP,RO,XI(K)
*
ccc        print*,'here calculated/marcs calculated'
ccc        do 1111 i=1,30
ccc         if (i.eq.17) goto 1111
ccc         comparison(i)=presmo(i)/pressure(k,i)
ccc1111    continue
ccc        print 1112,(comparison(i),i=1,10)
ccc        print 1112,(comparison(i),i=11,20)
ccc        print 1112,(comparison(i),i=21,30)
ccc        print*,' '
ccc        print 1113,presmo(24),pressure(k,24)
ccc        print*,' '
ccc1112    format(10(1x,f6.4))
ccc1113    format(2(2x,e10.4))
*
ccc      write(34,'(i3,15e10.3,/,3x,15e10.3)') k,presmo
*
* CALCULATED STANDARD OPaCITY
*

        CALL ABSKO(NEWT,1,T(K),PE(K),1,1,ABSK,SPRID)
        STNDOP=ABSK+SPRID
        if (mocode(1:7).eq.'Hoefner') then
*
* we place tau points at T, P, rho points.
* Then R should be at the same level, but in Hoefner's model,
* it isn't. So we recompute the R scale at integer tau points.
*
          if(k.eq.1) then
            tau(k)=drr(k)*stndop*ro*0.5
            rr(k)=rr(k)-drr(k)*0.5
            stndopprev=stndop
            roprev=ro
          else
            tau(k)=tau(k-1)+
     &             0.5*(drr(k)*ro*stndop+drr(k-1)*roprev*stndopprev)
            rr(k)=rr(k-1)-0.5*(drr(k)+drr(k-1))
            stndopprev=stndop
            roprev=ro
          endif
        else if (mocode(1:4).eq.'bowe'.or.mocode(1:4).eq.'BOWE'.or.
     &           mocode(1:4).eq.'alva') then
          if(k.eq.1) then
            tau(1)=drr(1)*stndop*ro*0.5
            stndopprev=stndop
            roprev=ro
          else
ccc            tau(k)=tau(k-1)+(rr(k-1)-rr(k))*stndop*ro
* corrected 10/12-1996 BPz
* at T(k), pe(k) pg(k) etc are defined at tau(k). r(k) is halfway between
* tau(k) and tau(k-1). This is the way MARCS works (maybe not... It seems 
* Z is midway between tau points but R is at tau points.). YES. That 
* last statement is correct R is at integer tau-points. BPz 15/11-2002
* Maybe Alva models are
* not constructed this way. Bowen's models are.
**            tau(k)=tau(k-1) + ( (drr(k-1)*roprev*stndopprev) +
**     &                (drr(k)*stndop*ro) )*0.5
** corrected 14/11-2007 BPz
            tau(k)=tau(k-1) + ( roprev*stndopprev +
     &                stndop*ro )*0.5*drr(k-1)
            stndopprev=stndop
            roprev=ro
          endif
************   Stagger models
        else if (mocode(1:7).eq.'Stagger') then
* compute tau-scale
          if (k.eq.1) then
            print*,'Stagger model'
            print*,' r           tau'
            tau(1)=stndop*(rr(2)-rr(1))*rhobow(1)/2.
            stndopprev=stndop
            print*,rr(1),tau(1)
          else
            tau(k)=tau(k-1)+(rr(k)-rr(k-1))*
     &         (stndop*rhobow(k)+stndopprev*rhobow(k-1))*0.5
            stndopprev=stndop
            print*,rr(k),tau(k)
          endif
        endif
************
   22   NEWT=0
c        print*,'check ro'
        write(6,66) k,tau(k),t(k),pe(k),pg,ro
66      format('tau,T,Pe,Pg,ro ',i3,1x,1pe10.3,1x,0pf8.1,3(1x,1pe10.3))
*******************************************************************
        if (mocode(1:7).eq.'Hoefner'.or.mocode(1:4).eq.'alva'.or.
     &      mocode(1:4).eq.'bowe'.or.mocode(1:7).eq.'Stagger') then
*
* write out model at babsma input format.
* It can be injected in MARCS35 for a global OS spectrum run (newmod=8)
*  BPz 14/11-2002
*
          write(imodut2,123) log10(tau(k)),T(k),log10(pe(k)),
     &                        log10(pg),xi(k),rr(k)
123       format(f8.4,1x,f8.1,2(1x,f8.4),1x,f5.2,1x,1pe15.8)
        endif
*******************************************************************
        if (hydrovelo) then
          WRITE(IMODUT,211) RR(K),TAU(K),T(K),PE(K),PG,RO,
     &                      XI(K),STNDOP,velocity(k)
          if (mocode(1:6).eq.'KURUCZ') then
            stop 'incompatible options: Kurucz and hydrovelo'
          endif
        else
          WRITE(IMODUT,201) RR(K),TAU(K),T(K),PE(K),PG,RO,
     &                      XI(K),STNDOP
        endif
*
* FINALLY, THE ABS.COEFF. RATIOS
*
        if (k.eq.1) then
          if (nametryck(6).ne.'C H  ') then
            print*,'Babsma, error in molecule. Not CH',nametryck(6)
            stop
          endif
          if (nametryck(5).ne.'O H  ') then
            print*,'Babsma, error in molecule. Not OH',nametryck(5)
            stop
          endif
        endif
        J0=1
        DO 24 I=2,NSET
          NLP=NL(I)
          DO 23 J=1,NLP
            CALL ABSKO(NEWT,1,T(K),PE(K),I,J,ABSK,SPRID)
            X(J0)=ABSK/STNDOP
            S(J0)=SPRID/STNDOP
            if (pureLTE) then
              x(j0)=x(j0)+s(j0)
              s(j0)=0.
            endif
            J0=J0+1
   23     CONTINUE
*
   24   CONTINUE
        WRITE(IMODUT,202)(X(J),S(J),J=1,NLQ)
        WRITE(IWRIT,302)(X(J),S(J),J=1,NLQ)
        NEWT=1

* save gas pressure
        pgk(k)=pg

   25 CONTINUE
      tsuswitch=.false.
      if (mocode(1:6).eq.'KURUCZ') then  
        WRITE(*,*) 'KURUCZ model converted to:' 
        print*,'log(tau',int(xls),'), T,  log(Pe),  log(Pg), Xi'
        do k=1,ntau
          WRITE(*,*) log10(TAU(K)),T(K),log10(PE(K)),log10(pgk(k)),
     &                      XI(K)
        enddo
      else      
        print*,' CONTROLE: k, R, tau, log10tau, T'
        do k=1,ntau
          print*,k,rr(k),tau(k),log10(tau(k)),t(k)
        enddo
      endif
cc      print*,'test G-band: k, logtau, logPH  logPH+  logPCH,  logPCN'
cc      do k=1,ntau
cc        print*,k,log10(tau(k)),log10(sngl(xmettryck(k,1))),
cc     &      log10(sngl(xiontryck(k,1))),
cc     &      log10(sngl(partryck(k,6))),log10(sngl(partryck(k,8)))
cc      enddo
*
  100 FORMAT(A4,I3,F6.0,F5.2,I2,F6.3)
  101 FORMAT(5F10.0)
  102 FORMAT(I5,I10,10I5)
  103 FORMAT(10F8.0)
  104 FORMAT(11X,L1)
  105 FORMAT(10F6.2)
  109 FORMAT(5F10.0)
  200 FORMAT(1X,'''',A,'''',I5,F10.2)
 2000 FORMAT(1X,'''',A,'''',1x,I5,1x,F10.2,' 0. 0 0.')
  201 FORMAT(1X,1P8E14.7)
  211 FORMAT(1X,1P8E14.7,1PE15.7)
  202 FORMAT(1X,1P6E11.4)
  203 FORMAT(1X,I5)
  204 FORMAT(1X,10F11.3)
  205 FORMAT(15X,F6.0,9X,F5.2,3X,A4,6X,I2)
  206 FORMAT(I4,3X,2E10.3,F6.0,4E10.3)
  300 FORMAT(' MODEL CODE=',A,5X,'  NTAU=',I3,5X,'DEPTH SCALE AT',
     &       F10.2,' AANGSTROEMS')
  301 FORMAT(' TAU=',F8.5,' T=',F6.0,' PE=',1PE10.2,' PG=',E10.2,
     &       ' PRATIO=',0PF6.3,' RO=',1PE10.2,' XI=',0PF6.3)
  302 FORMAT(' ',10F8.4)
  303 FORMAT('*** ERROR IN BABSMA ***, NLQ=',I5,' GREATER THAN
     &       DIM. FOR XS,SS,XLP')
*
      END
C***********************************************************************

