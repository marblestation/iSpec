      subroutine takemolec(kk,infoonly,molinquire,indexanswer)
*
* to add new molecules, just fill in with name as in Tsuji molecular file
* starting with the first '        ' at th ened of data block molinpresmo.
*    BPz 10/10-95
*
* this routine is to be used after a call to jon,
*  -if eqmol has been called also-, in order to get the pressures
* one needs placed in the common 'fullequilibrium'.
* This is the routine to change if one wants more or/and other
* molecular pressure to be kept.
* The values in the common fullequilibrium can then be used for
* computation of opacities, printing, etc.
* kk is the depth point to be adressed.
* 020890 BPlez
*
      include 'spectrum.inc'
      include 'tsuji.par'
      integer maxim
      parameter (maxim=1000)
      logical tsuswitch,tsuji,first
      character*20 molinpresmo(maxmol),nametryck(maxmol)
      doubleprecision parptsuji,xmettryck,xiontryck,partryck
      doubleprecision presneutral,presion,presion2,presion3
      common /tsuji/ tsuji,tsuswitch,nattsuji,nmotsuji,
     &             parptsuji(maxim+400)
      common /fullequilibrium/ partryck(ndp,maxmol),
     &  xmettryck(ndp,maxmet),xiontryck(ndp,maxmet),nametryck
      common /orderedpress/ presneutral(ndp,100),presion(ndp,100),
     &       presion2(ndp,100),presion3(ndp,100)
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO(30)
      common/ci1/fl2(5),parco(45),parq(180),shxij(5),tparf(4),
     &           xiong(16,5),eev,enamn,sumh,xkbol,nj(16),iel(16),
     &           nel,summ
* common from eqmol
      logical switer,infoonly
      character*20 mol(maxim),molinquire
      integer nelem(5,maxim),natom(5,maxim),mmax(maxim),nelemx(100),
     &          nmol,
     &        nmetal,nimax,marcsnelemx(20),marcsnelemj,j,atindex(20),
     &        molindex(maxmol),index,indexanswer
      real exponent(maxim),g0(100),g1(100),g2(100),g3(100),sumpress
      doubleprecision IP(100),KP(100),uiidui(100),eps,fp(100),
     &       ppmol(maxim),apm(maxim),c(maxim,5),p(100),ccomp(100),
     &       ipp(100),ippp(100),d00(maxim),qmol(maxim),
     &       reducedmass15(maxim)
      character*20 molcode(maxim)
      COMMON/COMFH1/C,NELEM,NATOM,MMAX,PPMOL,d00,qmol,APM,MOL,IP,
     &              ipp,ippp,g0,g1,g2,g3,CCOMP,exponent,reducedmass15,
     &              UIIDUI,P,FP,KP,eps,NELEMX,NIMAX,NMETAL,NMOL,switer,
     &              molcode
      real abund(16),anjon(16,5),h(5),part(16,5),dxi,f1,f2,f3,
     &           f4,f5,xkhm,xmh,xmy
      common/ci5/abund,anjon,h,part,dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy

      save atindex,molindex,first

      data molinpresmo/'H -                 ','H H                 ',
     &                 'H H +               ','H O H               ',
     &                 'O H                 ','C H                 ',
     &                 'C O                 ','C N                 ',
     &                 'C C                 ','N N                 ',
     &                 'O O                 ','N O                 ',
     &                 'N H                 ','C C H H             ',
     &                 'H C N               ','C C H               ',
     &                 '                    ','H S                 ',
     &                 'SiH                 ','C C C H             ',
     &                 'C C C               ','C S                 ',
     &                 'SiC                 ','SiC C               ',
     &                 'N S                 ','SiN                 ',
     &                 'SiO                 ','S O                 ',
     &                 'S S                 ','SiS                 ',
     &                 'TiO                 ','V O                 ',
     &                 'ZrO                 ','MgH                 ',
     &                 'CaH                 ','H F                 ',
     &                 'SiO                 ','H Cl                ',
     &                 'FeH                 ','SiH                 ',
     &                 'N O                 ','C H H H H           ',
     &                 'AlH                 ','CrH                 ',
     &                 'LaO                 ','TiH                 ',
     &                 'Y O                 ','NaH                 ',
     &                 '                    ','                    ',
     &                 '                    ','                    ',
     &                 '                    ','                    ',
     &                 '                    ','                    ',
     &                 '                    ','                    ',
     &                 '                    ','                    '/
      data first/.true./
      data marcsnelemx/ 1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 19,
     &                 20, 24, 26, 28, 21, 22, 23, 25/
      data molindex/maxmol*0/


* store the metals and ions
* The 16 first are indexed like in jon:
* H, He, C, N, O, Ne, Na, Mg, Al, Si, S, K, Ca, Cr, Fe, Ni
* 1   2  3  4  5   6   7   8   9  10 11 12  13  14  15  16
* and then:
* Sc, Ti, V, Mn
* 17  18 19  20
*
cc      if(nattsuji.gt.maxmet) stop 'takemolec: maxmet too small'

      if (infoonly) then
* get index for molecular pressure in partryck for molecule molinquire
*
        indexanswer=0
        do j=1,maxmol
          if (molinquire.eq.molinpresmo(j)) then
            indexanswer=j
          endif
        enddo
        if (indexanswer.eq.0) then
          print*,'takemolec infoonly: molec. not implemented',molinquire
          print*,' See and complete data block molinpresmo!'
cccc          print*,' complete also atomdata table with ID and IP'
          stop
        endif
        return
      endif

      if (first) then
        do i=1,100
          do j=1,ndp
            presneutral(j,i)=-1.0
            presion(j,i)=-1.0
            presion2(j,i)=-1.0
            presion3(j,i)=-1.0
          enddo
        enddo
* find index for atomic and molecular pressures
        do i=1,nmetal
          nelemi=nelemx(i)
          do j=1,20
            marcsnelemj=marcsnelemx(j)
            if (nelemi.eq.marcsnelemj) then
              atindex(j)=i
            endif
          enddo
        enddo
* the indexes from 1 to 30 for partryck correspond to
* the ones in presmo (common carc3). See data block molinpresmo.
* 31 is TiO, 32 is VO, 33 is ZrO, 34 is MgH., 35 is CaH , 36 is HF
*  ( compatible with new bsyn !!!)
*********************************************************************
        do i=1,maxmol
          if (molinpresmo(i).ne.'                    ') then
            j=1
            do while (molinpresmo(i).ne.mol(j).and.j.lt.nmol)
              j=j+1
            enddo
            if (molinpresmo(i).eq.mol(j)) then
              molindex(i)=4*nmetal+j
            else 
* this is when we reach the last molecule in Tsuji's list
* without finding the right molecule. we flag with negative pressure.
              molindex(i)=4*nmetal+nmol+2
            endif
          else
* partial pressure=0.0
            molindex(i)=4*nmetal+nmol+1
          endif
        enddo
        if (4*nmetal+nmol+2.gt.maxim+400) stop 
     &          'takemolec; maxim dim. too small! ERROR!'
        if (nmol+2.gt.maxim) stop 
     &          'takemolec; maxim dim. too small! ERROR!'
        parptsuji(4*nmetal+nmol+1) = 1.e-31
        parptsuji(4*nmetal+nmol+2) = -1.0
        mol(nmol+1)='zero pressure'
        mol(nmol+2)='not existing'
      endif
*
      do j=1,20
        xmettryck(kk,j)=parptsuji(atindex(j))
        xiontryck(kk,j)=parptsuji(atindex(j)+nmetal)
      enddo
* Set partition functions and ionisation fractions to the values used 
* in eqmol_pe/die_pe. Part(16,5) and anjon(16,5) are used in detabs for 
* the computation of continuous opacities.

cc      print*,'takemolec old part, new part for atoms for depth',kk

      do j=1,16
        if (first) then
          if ((nj(j).eq.4.and.g3(marcsnelemx(j)).eq.0.).or.
     &        (nj(j).gt.4) ) print*,
     &    'WARNING, takemolec: number of ionisation stages considered',
     &     ' for',marcsnelemx(j),
     &     ' lower than requested in input file jonabs.dat!'
        endif
        index=atindex(j)
        nelemi=nelemx(index)
cc        print*,nelemi,marcsnelemx(j),' =?'
        marcsnelemj=marcsnelemx(j)

cc        print*,part(j,1),part(j,2),part(j,3),part(j,4)
cc        print*,g0(marcsnelemj),g1(marcsnelemj),g2(marcsnelemj),
cc     &         g3(marcsnelemj)

        part(j,1)=g0(marcsnelemj)
        part(j,2)=g1(marcsnelemj)
        part(j,3)=g2(marcsnelemj)
        part(j,4)=g3(marcsnelemj)
        sumpress=parptsuji(index)+
     &           parptsuji(index+nmetal)+
     &           parptsuji(index+2*nmetal)+
     &           parptsuji(index+3*nmetal)
        anjon(j,1)=parptsuji(index)/sumpress
        anjon(j,2)=parptsuji(index+nmetal)/sumpress
        anjon(j,3)=parptsuji(index+2*nmetal)/sumpress
        anjon(j,4)=parptsuji(index+3*nmetal)/sumpress
      enddo
*
* now, the pressures indexed with the atomic number
* The arrays are initially set to -1.0 to allow selection later in bsyn of the
* elements not treated in eqmol. NOTE: does not work for
* presion3 not treated in eqmol, and which is 0.0 after next loop.
      do i=1,nmetal
        nelemi=nelemx(i)
        presneutral(kk,nelemi)=parptsuji(i)
        presion(kk,nelemi)=parptsuji(i+nmetal)
        presion2(kk,nelemi)=parptsuji(i+2*nmetal)
        presion3(kk,nelemi)=parptsuji(i+3*nmetal)
      enddo
* and the molecules
      do i=1,maxmol
        partryck(kk,i)=parptsuji(molindex(i))
        if (first) then
          nametryck(i)=mol(molindex(i)-4*nmetal)
        endif
      enddo
      first=.false.
***************************** debug*******************
ccc      if (kk.eq.1) print*,'same?',(partryck(kk,i),i=1,35)
***************************** debug*******************
      return
      end

