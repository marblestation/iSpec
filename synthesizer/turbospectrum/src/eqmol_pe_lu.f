*
      subroutine eqmol_pe(tt,pgin,pgas,pe,xih,xihm,kk,niter,skiprelim)
*     
* this version for computing pgas from T,Pe with LU-decompsition of
* the matrix in die_pe_lu. BPz 09/07-1999
* double precision
* It uses partition functions for atomic and molecular species as input
* from which it computes a consistent set of equilibrium constants.
* 4 ionisation stages (neutral, 1st, 2nd, 3rd) are included for atoms and
* more can easily be added.
* BPz 28/08-1995

C        GG2=N(HII)/N(HI), GC=N(CII)/N(CI) ETC.
C        XIH = THE IONIZATION ENERGY OF HYDROGEN
C        XIHM= THE 'DISSOCIATION ENERGY' OF H-
C        XKHM= THE 'DISSOCIATION CONSTANT' OF H-

C        THE ROUTINE GIVES FH, FC, FO, FN, FE. FH=P(HI)/PH, FC=P(CI)/PH ETC.,
C        WHERE PH=NH*KT (NH IS THE NUMBER OF HYDROGEN NUCLEI PER CM3).
C        THE SUBSCRIPT IN PK(I)  HAS THE FOLLOWING MEANING
C        I=1 H-, 2 H2, 3 H2+, 4 H2O, 5 OH, 6 CH, 7 CO, 8 CN, 9 C2, 10 N2,
C         11 O2, 12 NO, 13 NH, 14 C2H2, 15 HCN, 16 C2H, 17 -, 18 HS
C         19 SIH, 20 C3H, 21 C3, 22 CS, 23 SIC, 24 SIC2, 25 NS
C         26 SIN, 27 SIO, 28 SO, 29 S2, 30 SIS

c							b. plez   7/7/87
C							modif Pe- 6/11/87
C                                                       modif    22/11/88
* for a change of delta(D0) of the dissociation energy (in eV) of a molecule,
* change the c(2) coefficient in log(Kp) expansion by -delta(D0)
*
* VERSION FOR INCLUSION IN THE MARCS CODE. Uses only one depth point at a time
*                    BPz 010890
* 
* BPz 290595
* common block plez contains the names of all the molecules. variable names
* is just a copy of mol. This is done to avoid carrying around the common
* block COMFH1. This is used in takemolec.

        implicit none
        include 'spectrum.inc'
        integer maxim
        parameter (maxim=1000)
	character*128 MET,MOLEC
        common/files/molec,met
	REAL KPLOG,IPI,fictpres_c,fictpres_h,xlog,pmoll,pionl,
     &       plog,fplog,pglog,theta,xmytsuji,ipii,ipiii,exponent(maxim)
        real fictpres_c_noco,fract_c_in_atom_noco,fract_c_in_ion,
     &  fract_c_in_atom,fract_c_in_ion_noco,fract_c_in_mol,
     &  fract_c_in_mol_noco,seuil,cclogi,pdfpl,pelog,avo
        character*20   names(maxim),molinquire
        integer tselem(100),kk,iii,niter,indexanswer
        logical switer,found(maxim),molkeep(5),molekeep,converge,
     &          infoonly,skiprelim
        common /plez/ names,tselem,fictpres_h
        integer i,j,m,mmax(maxim),mmaxj,natom(5,maxim),nelem(5,maxim),
     &          nelemx(100),nmetal,nmol,nimax,jjj,
     &          k,nelemi
        character*20   MOL(maxim)
* really real, not supposed to become doubleprecision
        real tem,pg,pe,pgin,tt,pgas,g0(100),g1(100),g2(100),g3(100)
* may become dbleprec
        doubleprecision IP(100),KP(100),uiidui(100),eps,fp(100),
     &       ppmol(maxim),apm(maxim),c(maxim,5),ccomp(100),p(100),
     &       econst,exp10,x,ppk(30),ipp(100),ippp(100),d00(maxim),
     &       qmol(maxim),reducedmass15(maxim)
        character*20 molcode(maxim)
        COMMON/COMFH1/C,NELEM,NATOM,MMAX,PPMOL,d00,qmol,
     &                APM,MOL,IP,ipp,ippp,g0,g1,g2,g3,
     &                CCOMP,exponent,reducedmass15,
     &                UIIDUI,P,FP,KP,eps,NELEMX,
     &                NIMAX,NMETAL,NMOL,switer,molcode,elem
        integer natomm(5),nelemm(5),ig1,ig0,nel,iel(16),nj(16)
        real fl2(5),parco(45),parq(180),shxij(5),tparf(4),xiong(16,5),
     &       eev,enamn,sumh,xkbol,summ
        common/ci1/ fl2,parco,parq,shxij,tparf,xiong,eev,enamn,sumh,
     &              xkbol,nj,iel,nel,summ
        doubleprecision presneutral,presion,presion2,presion3
        common/orderedpress/presneutral(ndp,100),presion(ndp,100),
     &                      presion2(ndp,100),presion3(ndp,100)
        character*2 atominclude(100)
        common/species/atominclude

        character*2 elem(100), elemnt(100), elemxi
        real CCLOG(100),molenergy
        real xmass(maxim+400),atmass(100)
        doubleprecision ndensity,molweight
        common /density/ndensity,molweight,xmass,atmass
        real xkh2,xkh2p,deh2,deh2p,xih,xihm,gg2,deh2pnodis,deh2nodis

* commons from injon
* declarations should not become dbleprec!
      real abund,anjon,h,part,dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy
      COMMON/CI5/ABUND(16),ANJON(16,5),H(5),PART(16,5),DXI,
     & F1,F2,F3,F4,F5,XKHM,XMH,XMY
      real eh,fe,fh,fhe,fc,fce,fn,fne,fo,foe,fk,fke,fs,fse
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
      real pk
      integer jonnmol
      COMMON /CMOL2/jonnmol,PK(30)
      real eh2,eh2p,ehm,ehj,eh2o,eoh,ech,eco,ecn,ec2,en2,eo2,eno,enh
      real DIS(10)
      real absc,abti,abv,abmn,abco
      common /auxabund/ absc,abti,abv,abmn,abco
* common for partial pressures
      logical tsuswitch,tsuji
      integer nattsuji,nmotsuji
      doubleprecision parptsuji
      common /tsuji/ tsuji,tsuswitch,nattsuji,nmotsuji,
     &               parptsuji(maxim+400)
      character*128 filmet,filmol
      real factmet,rho,ejontsuji
      common /filetsuji/ filmet,filmol
      common/rhotsu/rho,xmytsuji,ejontsuji
      common/zzzz/factmet

      real inputabund(92),inputamass(92,0:250),isotopfrac(92,0:250)
      character*2 inputaname(92)
      common/refabundances/ inputabund,inputamass,inputaname,isotopfrac

 
      logical first,test,testh
      DATA ELEMNT(99),ELEMNT(100)/'E-','H*'/
      data DIS/9.50,4.38,3.47,11.11,7.90,6.12,9.76,5.12,6.51,3.21/
      data first/.true./, infoonly/.false./
      save first,dis,elemnt
*
        exp10(x)=exp(2.302585*x)
	ECONST=4.342945E-1
	AVO=0.602217E+24
*
** do NOT use xmass(99) for electrons, it is reserved for a molecule! **
*
        atmass(99)=5.4858e-4
C
* number of atoms, max iterations, convergence criterium, full or half
* pressure correction (see die_pe.f)
        nimax=2000
        eps=1.e-4
* should not become false. Convergence is then too slow.
        switer=.true.
C
* read atomic data file

      if (first) then
	WRITE(6,6102)
*
cc        open(UNIT=26,FILE=filmet,STATUS='OLD')  ! suppressed by BPz 29/02-2012
cc* file 26 opened in input.f, as scratch, contains a list of atomic
cc* species to include in chemical equilibrium. One column, character*2.
cc        rewind(26)
        nmetal=0
	do i=1,100
cccc* name, number of e-, ionisation pot, partition fct for neutral and
cccc* first ion, log of abundance.
cccc          READ (26,*,end=99) ELEMXI,NELEMI,IPI,IG0,IG1,CCLOGI,xmass(i),
cccc     &                ipii,ipiii
cc          read(26,*,end=99) elemxi  ! suppressed by BPz 29/02-2012
         elemxi=atominclude(i)
         if (elemxi.ne.'  ') then
          nelemi=0
          do iii=1,92
            if(elemxi.eq.inputaname(iii)) then
              if (inputabund(iii).gt.0.) then
                nelemi=iii
                goto 111
              else
* skip atomic species with zero abundance (flagged with -29. in
* makeabund.f, changed to 0 . in bsyn or in injon).
                goto 112
              endif
            endif
          enddo
111       continue
          if (nelemi.eq.0) then
            print*,'Eqmol_pe; element not found: ',elemxi
            stop
          endif
            
          nmetal=nmetal+1
          NELEMX(nmetal)=NELEMI
          tselem(nmetal)=nelemi
          ELEM(nmetal)=ELEMXI
* here we assume mass of standard mix
          atmass(nelemi)=inputamass(nelemi,0)
          xmass(nmetal)=atmass(nelemi)
          ELEMNT(NELEMI)=ELEMXI
* ionization potentials are set in die_pe.
cc          IP(NELEMI)=IPI
cc          ipp(nelemi)=ipii
cc          ippp(nelemi)=ipiii
          cclogi=log10(inputabund(nelemi))
          CCLOG(NELEMI)=CCLOGI
ckeep for debug
	  WRITE(6,6103) ELEMXI,NELEMI,CCLOGI+12.
c 
112       continue
         endif
        enddo
99      continue
        do i=1,nmetal
          xmass(i+nmetal)=xmass(i)
          xmass(i+2*nmetal)=xmass(i)
          xmass(i+3*nmetal)=xmass(i)
        enddo
        nattsuji=nmetal
cc        close(26)
*
* normalization to H abundance added 11/03/93, and computation of 
* gram of */g of H
        xmytsuji=0.0
        if (nelemx(1).ne.1) then
          print*,'Eqmol_pe; WARNING: first element in list is not H.'
          print*,'    It is: ',nelemx(1)
        endif
        ccomp(1)=exp10(dble(cclog(1)))
        do i=2,nmetal
          nelemi=nelemx(i)
          ccomp(nelemi)=exp10(dble(cclog(nelemi)))
          ccomp(nelemi)=ccomp(nelemi)/ccomp(1)
        enddo
        ccomp(1)=1.0
cc        print*,'Chemical composition from atomic list'
        do i=1,nmetal
          nelemi=nelemx(i)
cc          print*,nelemi,elemnt(nelemi),log10(sngl(ccomp(nelemi)))+12.
          xmytsuji=xmytsuji+atmass(nelemi)*ccomp(nelemi)
        enddo
        xmytsuji=xmytsuji/(ccomp(1)*atmass(1))
* not initialized if not called injon before. Skippping.
ccc        enamn=eev/(xmytsuji*xmh)
c
* reads molecular data file
	J=0
        
	open(UNIT=26,FILE=filmol,STATUS='OLD')
1010   	J=J+1
1011    READ (26,*,end=1014) MOL(J),(C(J,K),K=1,5),MMAX(J),
     &    (NELEMM(M),NATOMM(M),M=1,mmax(j))
        molekeep=.true.
        do m=1,mmax(j)
          molkeep(m)=.false.
          if (nelemm(m).eq.99) then
            molkeep(m)=.true.
          endif
          do iii=1,nmetal
            if (nelemm(m).eq.nelemx(iii)) then
              molkeep(m)=.true.
            endif
          enddo
          molekeep=molekeep.and.molkeep(m)
        enddo
        if (.not.molekeep) then
          print*,' Eqmol_pe; rejecting : ',mol(j)
          goto 1011
        endif
****************************** debug**************************
ccc        if (j.eq.1) open(33,file='outputtest',status='unknown')
ccc	write (33,5012) MOL(J),(C(J,K),K=1,5),MMAX(J),
ccc     &    (NELEMM(M),NATOMM(M),M=1,4)
****************************** debug**************************

cccc      print*,'eqmol reading molec',j,mol(j)
        
	MMAXJ=MMAX(J)
	IF(MMAXJ.EQ.0) GOTO 1014
	DO M=1,MMAXJ
          NELEM(M,J)=NELEMM(M)
          NATOM(M,J)=NATOMM(M)
        enddo
1110    GOTO 1010
C
1014    NMOL=J-1
        close(26)
        if (4*nmetal+nmol.gt.maxim+400) stop 
     &           'eqmo: dim=maxim+400 too small!'
	DO I=1,NMETAL
          NELEMI=NELEMX(I)
          P(NELEMI)=1.0E-20
        enddo
        first=.false.
      endif

c
* tt(ndp) model temperatures to which equilibrium has to be computed
* ppe(ndp) same for e- pressure
* pgg(ndp) same for gas pressure
* nh(ndp)  no need
        THETA=5040./tt
	TEM=tt
	PG=Pgin
c
cc        print*,'eqmol_pe: t, pe ',tem,pe
cc        print*,'eqmol before molecpartf, nmol',nmol
        call molecpartf(tem,found)
cc        print*,'eqmol after molecpartf, nmol',nmol
        nmotsuji=nmol
	call die_pe(tem,pe,pg,found,converge,niter,skiprelim)
c
        if (.not.converge) then
          pgas=-1.0
          return
        endif

        pgas=pg
	PgLOG=LOG10(Pg)
        pelog=log10(pe)
C
C ndensity is number density and molweight is the mean molecular weight
c  in amu.
        molweight=0.0
        do i=1,4*nmetal+nmol
          parptsuji(i)=max(parptsuji(i),1.d-99)
          molweight=molweight+parptsuji(i)*xmass(i)
        enddo
        ndensity=pg/1.38066e-16/tem
        molweight=(molweight+pe*atmass(99))/1.38066e-16/tem
        molweight=molweight/ndensity

C compute fictitious pressure of H and C.

	DO J=1,NMOL
        test=.false.
        if (test) then

        do jjj=1,mmax(j)
          if (nelem(jjj,j).eq.0) then
          print 1111,nelem(jjj,j),
     &         mol(j),parptsuji(j+4*nmetal)*natom(jjj,j)/
     &        (pgas*exp10(dble(cclog(nelem(jjj,j))))),
     &           log10(parptsuji(j+4*nmetal))
1111      format(i3,2x,a8,e10.2,x,f6.2)
         endif
        enddo

        endif
**************************************
        enddo
        do i=1,maxim
          names(i)=mol(i)
        enddo
**************************************
* print out fraction of C (or H, if testh=.t.) locked in each species.
* only species with containing a fraction larger than seuil are printed.
* BPz 1/6/95
*
* must have test=.true. for computation of fictpres_h, always needed!!
      test=.true.
      testh=.false.
      if (test) then

      seuil=0.

ccc      open(30,file='fractions_c_in_species.dat',status='unknown')
ccc      fictpres_c=0.0
      fictpres_h=0.0
      do j=1,nmetal
        if (nelemx(j).eq.1) then
          fictpres_h=parptsuji(j)+parptsuji(j+nmetal)+
     &          parptsuji(j+2*nmetal)+parptsuji(j+3*nmetal)
ccc        else if (nelemx(j).eq.22) then
ccc          fictpres_c=parptsuji(j)+parptsuji(j+nmetal)+
ccc     &          parptsuji(j+2*nmetal)+parptsuji(j+3*nmetal)
        endif
      enddo
      do j=1,nmol
        do jjj=1,mmax(j)
ccc          if (nelem(jjj,j).eq.22) then
ccc            fictpres_c = fictpres_c + 
ccc     &                   parptsuji(j+4*nmetal)*natom(jjj,j)
ccc          endif
          if (nelem(jjj,j).eq.1) then
            fictpres_h = fictpres_h + 
     &                   parptsuji(j+4*nmetal)*natom(jjj,j)
          endif
        enddo
      enddo

ccc      do j=1,nmetal
ccc        if (nelemx(j).eq.22) then
ccc          i=nelemx(j)
ccc          fract_c_in_atom=parptsuji(j)/fictpres_c
ccc          fract_c_in_ion=parptsuji(j+nmetal)/fictpres_c
ccc          if (fract_c_in_atom.ge.seuil) then
ccc          write(30,1112) char(39),elemnt(i),char(39),
ccc     &               log10(fract_c_in_atom),
ccc     &           log10(parptsuji(j)), tt,log10(pgas)
ccc          endif
ccc          if (fract_c_in_ion.ge.seuil) then
ccc          write(30, 1113) char(39),elemnt(i),char(39),
ccc     &               log10(fract_c_in_ion),
ccc     &           log10(parptsuji(j+nmetal)), tt,log10(pgas)
ccc          endif
ccc        endif
ccc      enddo
ccc      do j=1,nmol
ccc        do jjj=1,mmax(j)
ccc          if (nelem(jjj,j).eq.22) then
ccc            fract_c_in_mol = 
ccc     &          parptsuji(j+4*nmetal)*natom(jjj,j)/fictpres_c
ccc            if(fract_c_in_mol.ge.seuil) then
ccc            write(30, 1112) char(39),
ccc     &         mol(j),char(39),log10(fract_c_in_mol),
ccc     &           log10(parptsuji(j+4*nmetal)), tt,log10(pgas)
1112        format(a1,a8,a1,x,f6.2,x,f6.2,3x,f6.0,x,f7.3)
1113        format(a1,a6,'+ ',a1,x,f6.2,x,f6.2,3x,f6.0,x,f7.3)
ccc            endif
ccc          endif
ccc        enddo
ccc      enddo
      if (testh) then
      do j=1,nmetal
        if (nelemx(j).eq.1) then
          i=nelemx(j)
          write(20, 1112) char(39),elemnt(i),char(39),
     &               log10(parptsuji(j)/fictpres_h),0.0,
     &           log10(parptsuji(j)), tt,log10(pgas)
          write(20, 1113) char(39),elemnt(i),char(39),
     &               log10(parptsuji(j+nmetal)/fictpres_h),0.0,
     &           log10(parptsuji(j+nmetal)), tt,log10(pgas)
        endif
      enddo
      do j=1,nmol
        do jjj=1,mmax(j)
          if (nelem(jjj,j).eq.1) then
            write(20, 1112) char(39),
     &         mol(j),char(39),log10(parptsuji(j+4*nmetal)*natom(jjj,j)/
     &                fictpres_h),0.0,
     &           log10(parptsuji(j+4*nmetal)), tt,log10(pgas)
          endif
        enddo
      enddo
      endif

      endif
*
* calculates quantities needed by jon (were computed in moleq before)
*
      call takemolec(kk,infoonly,molinquire,indexanswer)
      rho=0.
      do 9874 i=1,4*nmetal+nmol
        rho=rho+parptsuji(i)*xmass(i)
9874  continue
      rho=rho+pe*atmass(99)
      rho=rho*1.2123e-8/tt
* 1.2123e-8 == mH/k = amu*1.00797/k
cc        print*,'t,pe,ndensity,molweight ro',
cc     &           tem,pe,ndensity,molweight,rho
* pk

      call molfys(tt,xkh2,xkh2p,deh2,deh2p,deh2nodis,deh2pnodis)
*
      jonnmol=30
      do j=1,nmol
        if (mol(j).eq.'H -                 ') then
          ppk(1)=apm(j)
        else if (mol(j).eq.'H H                ') then
          ppk(2)=apm(j)
        else if (mol(j).eq.'H H +              ') then
cc corrected 4/12-2007: 
cc apm(H2+) is not Kp(H2+). pk should be Pe/Kp.   BPz
cc          ppk(3)=apm(j)
          do jjj=1,nmetal
            if (nelemx(jjj).eq.1) then
              ppk(3)=apm(j)*pe*parptsuji(jjj+nmetal)/parptsuji(jjj)
            endif
          enddo                                         
        else if (mol(j).eq.'H O H              ') then
          ppk(4)=apm(j)
        else if (mol(j).eq.'O H                ') then
          ppk(5)=apm(j)
        else if (mol(j).eq.'C H                ') then
          ppk(6)=apm(j)
        else if (mol(j).eq.'C O                ') then
          ppk(7)=apm(j)
        else if (mol(j).eq.'C N                ') then
          ppk(8)=apm(j)
        else if (mol(j).eq.'C C                ') then
          ppk(9)=apm(j)
        else if (mol(j).eq.'N N                ') then
          ppk(10)=apm(j)
        else if (mol(j).eq.'O O                ') then
          ppk(11)=apm(j)
        else if (mol(j).eq.'N O                ') then
          ppk(12)=apm(j)
        else if (mol(j).eq.'N H                ') then
          ppk(13)=apm(j)
        else if (mol(j).eq.'C C H H            ') then
          ppk(14)=apm(j)
        else if (mol(j).eq.'H C N              ') then
          ppk(15)=apm(j)
        else if (mol(j).eq.'C C H              ') then
          ppk(16)=apm(j)
        else if (mol(j).eq.'H S                ') then
          ppk(18)=apm(j)
        else if (mol(j).eq.'SiH                ') then
          ppk(19)=apm(j)
        else if (mol(j).eq.'C C C H            ') then
          ppk(20)=apm(j)
        else if (mol(j).eq.'C C C              ') then
          ppk(21)=apm(j)
        else if (mol(j).eq.'C S                ') then
          ppk(22)=apm(j)
        else if (mol(j).eq.'SiC                ') then
          ppk(23)=apm(j)
        else if (mol(j).eq.'SiC C              ') then
          ppk(24)=apm(j)
        else if (mol(j).eq.'N S                ') then
          ppk(25)=apm(j)
        else if (mol(j).eq.'SiN                ') then
          ppk(26)=apm(j)
        else if (mol(j).eq.'SiO                ') then
          ppk(27)=apm(j)
        else if (mol(j).eq.'S O                ') then
          ppk(28)=apm(j)
        else if (mol(j).eq.'S S                ') then
          ppk(29)=apm(j)
        else if (mol(j).eq.'SiS                ') then
          ppk(30)=apm(j)
        endif
      enddo
      ppk(17)=1.
      do j=1,jonnmol
        ppk(j)=pe/ppk(j)
      enddo
      PPK(4)=PE*PPK(4)
      PPK(14)=PE*PE*PPK(14)
      PPK(15)=PE*PPK(15)
      PPK(16)=PE*PPK(16)
      PPK(20)=PE*PE*PPK(20)
      PPK(21)=PE*PPK(21)
      PPK(24)=PE*PPK(24)
      do j=1,jonnmol
* pick up single precision values
        pk(j)=ppk(j)
      enddo
* fe, fh, etc
      do i=1,nmetal
        j=nelemx(i)
        if (j.eq.1) then
          fh=parptsuji(i)/fictpres_h
          gg2=parptsuji(i+nmetal)/parptsuji(i)
        else if (j.eq.6) then
          fc=parptsuji(i)/fictpres_h
        else if (j.eq.7) then
          fn=parptsuji(i)/fictpres_h
        else if (j.eq.8) then
          fo=parptsuji(i)/fictpres_h
        else if (j.eq.14) then
          fk=parptsuji(i)/fictpres_h
        else if (j.eq.16) then
          fs=parptsuji(i)/fictpres_h
        endif
      enddo
      fe=pe/fictpres_h
 
      fhe=fh/fe
      fce=fc/fe
      fne=fn/fe
      foe=fo/fe
      fke=fk/fe
      fse=fs/fe
      f1=fh
      f2=gg2*fh
      f3=fh*ppk(1)
      f4=fh*fhe*gg2*ppk(3)
      f5=fh*fhe*ppk(2)
C
C COMPUTATION OF THE INNER ENERGY. DEH2 AND DEH2P ARE THE SUM OF
C ROTATION AND VIBRATION ENERGIES (IN EV PER MOLECULE) FOR
C H2 AND H2+. DIS(I) IS THE DISSOCIATION ENERGY FOR THE MOLECULE (I+3)
C IN THE LIST OF MOLECULES (VALUES ARE FROM TSUJI). FOR THESE MOLECULES
C THE ROTATION AND VIBRATION ENERGIES ARE NEGLECTED.
C
c      EH2=(-2.*XIH+DEH2)*FHE*FH*PPK(2)
c      EH2P=(DEH2P-XIH)*FHE*FH*GG2*PPK(3)
c      EHM=-(XIHM+XIH)*FH*PPK(1)
c      EHJ=-XIH*FH
CCC  old and wrong (dicovered 19/09-96)    EH2O=-(2.*XIH+DIS(1))*FHE*FH*FO*PPK(4)
c      EH2O=-(2.*XIH+DIS(1))*FHE*FHE*FO*PPK(4)
c      EOH=-(XIH+DIS(2))*FOE*FH*PPK(5)
c      ECH=-(XIH+DIS(3))*FHE*FC*PPK(6)
c      ECO=-DIS(4)*FCE*FO*PPK(7)
c      ECN=-DIS(5)*FCE*FN*PPK(8)
c      EC2=-DIS(6)*FCE*FC*PPK(9)
c      EN2=-DIS(7)*FNE*FN*PPK(10)
c      EO2=-DIS(8)*FOE*FO*PPK(11)
c      ENO=-DIS(9)*FNE*FO*PPK(12)
c      ENH=-(DIS(10)+XIH)*FNE*FH*PPK(13)
c      EH=EH2+EH2P+EHM+EHJ+EH2O+EOH+ECH+ECO+ECN+EC2+EN2+EO2+ENO+ENH
c      print*,'old EH = ',eh
C                            NOTE THAT ENERGIES INCLUDE ONLY MOLECULES 1-13
**********************************************************
      if (abs(xkbol/1.38e-16-1.).gt.0.1) goto 1963
* In that case we are calling eqmol without prior call to injon/jon.
* The common block CI1 is not initialized. For bsyn 31/10-96 BPz

      eh=-xih*presneutral(kk,1)
      do j=1,nmol
        if (mol(j).eq.'H H                ') then
          molenergy=(deh2nodis-(2.*xih+d00(j)))*
     &               parptsuji(j+4*nmetal)
        else if (mol(j).eq.'H H +              ') then
          molenergy=(deh2pnodis-(2.*xih+d00(j)))*
     &               parptsuji(j+4*nmetal)
        else
          molenergy=-d00(j)
          mmaxj=mmax(j)
          do m=1,mmaxj
            if (nelem(m,j).eq.1) then
              molenergy=molenergy-xih*natom(m,j)
            endif
          enddo
          molenergy=molenergy*parptsuji(j+4*nmetal)
        endif
        eh = eh + molenergy
      enddo
cc      print*,'enamn rho k t',enamn, rho,xkbol,tt,eh,enamn*xkbol*rho*tt
      eh=eh*eev/(enamn*rho*xkbol*tt)
cc      print*,'new EH = ',eh
* we try to compute a better ejon (ionisation energy of all elements except
* H. We account for molecule formation (?). BPz 20/09-1996
      ejontsuji=0.
      do i=1,nmetal
        nelemi=nelemx(i)
        if (nelemi.ne.1) then
          ejontsuji=(ip(nelemi)-dxi)*presion(kk,nelemi)+
     &       (ipp(nelemi)-2.*dxi)*presion2(kk,nelemi)+
     &       (ippp(nelemi)-3.*dxi)*presion3(kk,nelemi)+ejontsuji
        endif
      enddo
      ejontsuji=ejontsuji*eev/(enamn*rho*xkbol*tt)
1963  continue
**************************************
C
C------------formats--------------------------------------------------------
c
  51	FORMAT(7E11.4)
 501    FORMAT (I4)
 500	FORMAT (3(1X,I5))
 600	FORMAT (7(1X,I4,1X,A4,1X))
 620	FORMAT (7(1X,I4,1X,A4,'+',1X))
 650	FORMAT (7(1X,I4,1X,A8,1X))
5000    FORMAT (2I5,2F10.5,I10)
5001    FORMAT (A4,I6,F10.3,2I5,F10.3,F10.0)

CC format to write tsuji's data in file readable in free format

5012    FORMAT ('''',A8,'''',x,1pE11.5,x,4(1pE12.5,x),I1,4(I3,I4))

5021    FORMAT (F10.3,E12.5,E12.6)
5030    FORMAT (A)
5031    FORMAT (1X,A)
6031    FORMAT(1H1,20A4/)
6091    FORMAT (/,10X,'LOG PG=',F8.4,20X,'LOG PE=',F8.4,20X,'PE=',E13.6
     2  /,10X,'THETA =',F8.4,20X,'TEMP. =',F8.0,20X,'PROF. =',E14.6/)
cc6102    FORMAT (1H0, ' ELEMENT  ATOMIC NUMBER       I.P.        G(0)   G
cc     1(1)    LOGN/NH')
cc6103    FORMAT(1H,5X,A4,8X,I5,3X,F10.3,5X,2I5,3X,F10.5)
6102    FORMAT (' ELEMENT  ATOMIC NUMBER   LOG(ABUNDANCE)')
6103    FORMAT(5X,A2,8X,I5,3X,F10.3)
6300    FORMAT(1H0,' EQUILIBRIUM PARTIAL PRESSURES OF THE GASEOUS',
     &    ' ATOMS'  ///)
6301    FORMAT(1H0,'ELEMENT',3X,'LOG N(E)/N(H)',4X,'P(E)',6X,'LOG P(E)'
     1,2X,'LOG P(A)',2X,'LOG P(A)/P(E)'/)
6302    FORMAT(1H,1X,A4,1X,I2,4X,F10.3,1X,E11.4,2F10.3,4X,F10.3)
6307    FORMAT(1H0,/      ' P(E) **** FICTITIOUS PRESSURE OF THE NUCLE
     1US OF THE ELEMENT',/    ' P(A) ****PARTIAL PRESSURE OF THE MONA
     2TOMIC GAS OF THE ELEMENT')
6092    FORMAT (1H1,3(5X,'MOLECULE   LOG P   LOG P/PG  LOG KP')//)
6495    FORMAT ( 3(8X,A4,'+',3F9.3))
6496    FORMAT ( 3(9X,A4,3F9.3))
6992    FORMAT (///,3(9X,'ION     LOG P   LOGP/PG  LOG KP ')//)
7000    FORMAT (21(A9,I3))
C
1100    return
        END
