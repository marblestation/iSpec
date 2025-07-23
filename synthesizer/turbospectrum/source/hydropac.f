      subroutine hydropac(lunit,xlboff,nlte,
     &                    nlte_species,maxlevel,modnlevel,b_departure,
     &                    modenergy,modg,modion,modid)
!
! include NLTE capability with departure coefficients as input. BPz 17/11-2020
!
! new version for use with hbop.f  ! H line + continuum bf opacity
! Bpz 02/02-2019
! BPz: inclusion of nlte case. Necessary regardless of case for HI lines, as 
! we need emissivity (stored in source_function), and not only extinction coefficient.

      implicit none

      include 'spectrum.inc'

      character comment*100,species*20
      integer ion,nline,ntau,l,k,i,
     &        nmy,nlbldu,iint,iweak,lunit,maxlam,
     ,        jlcont, nlcont
      doubleprecision xl1,xl2,del,xlboff,xlambda
      doubleprecision source_function,source_f

! 150 is the max allowed number of H lines
      doubleprecision hlambda(150)
      real xlo(150),xup(150),gf(150),npop(150)
      real cont,total,contrib,planckfct
      integer nlo(150),nup(150)
      character*9 lname(150)
!
      real pe,t,pg,xi,mum,ro,eps,xmyc,
     &     ne(ndp),fact,
     &     nh1(ndp),nhe1(ndp),
     &     hckt(ndp),abso,absos,absocont,
     &     dopple(ndp),xkapr,cross,absoscont
      logical  contonly
      logical  lineonly
! NLTE
      logical nlte,nlte_species
      character*40 modid(*)
      real xlsingle,bpl,expcorr,ediff
      integer modnlevel,modion(*),modnlev,nlevlist,maxlevel
      real b_departure(ndp,0:maxlevel),modenergy(*),modg(*),bd(1000)
      real*8 planck,clight,electron,boltzmann
      real kayser2eV,ev2kayser,hck,kbcgs
*      
      common /atmos/ t(ndp),pe(ndp),pg(ndp),xi(ndp),mum(ndp),ro(ndp),
     &               ntau
      doubleprecision presneutral,presion,presion2,presion3
      common /orderedpress/ presneutral(ndp,100),presion(ndp,100),
     &                      presion2(ndp,100),presion3(ndp,100)
      common /pieces/ xl1,xl2,del,eps,nmy,nlbldu,iint,xmyc,iweak
      common/large/ xlambda(lpoint),source_function(ndp,lpoint),
     & maxlam,abso(ndp,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)
      common/rossc/ xkapr(ndp),cross(ndp)
      common/continuum/nlcont,jlcont(lpoint)
      
      data clight /2.99792458d8/, planck /6.62607015d-34/
      data boltzmann /1.380649d-23/, electron /1.602176634d-19/

!      kayser2ev = 0.0001239854563934909
!      ev2kayser = 8065.4621040899665
      kbcgs=sngl(boltzmann*1.d7)
      kayser2eV = sngl(planck * clight / electron * 1.d2)
      ev2kayser=1./kayser2ev
      hck=sngl(planck * clight / boltzmann * 1.d10)     ! ready for lambda in AA !
! 
      if (nlte_species) then
        do i=1,modnlevel-1
          if (modion(i).ne.1) then
            print*,'hydropac: problem with HI levels'
            stop 'stop in hydropac'
          endif
        enddo
        if (modion(modnlevel).ne.2) then
          print*,'HII level missing in model atom !'
          stop 'stop in hydropac'
        endif
      endif
      modnlev=max(modnlevel-1,0)   ! modnlev=0 if LTE, or =number of HI levels if NLTE

! set populations to 0. (computed in hbop)
      npop=0.0
! compute lines only 
      contonly=.false.
      lineonly=.true.
! read line list
      rewind (lunit)
      read(lunit,*) species,ion,nline
      read(lunit,*) comment
*      print*,species,ion,nline,comment
      if (species(1:6).ne.'01.000'.or.ion.ne.1) then
        print*, 'wrong H line data file!'
        print*,species,ion
        stop 'ERROR!'
      endif
c bsyn uses air wavelengths, and so does hbop.f, except for Lyman series
c      
      nlevlist=0
      do i=1,nline
!       oscillator strengths are computed in hbop.f
        read(lunit,*) hlambda(i),nlo(i),nup(i),xlo(i),xup(i),gf(i),
     &                lname(i)
        nlevlist=max(nlevlist,nup(i))               ! number of levels in line list
! check levels in line list vs levels in model atom !
        if (modnlevel.gt.0) then
          if (nlo(i).le.modnlev) then
            ediff = abs(xlo(i)*ev2kayser - modenergy(nlo(i)))/
     &               modenergy(nlo(i))
            if (ediff.gt.1.e-3) then
              print*,'H level identification problem'
              print*,'model atom ', modenergy(nlo(i)), 
     &               'line list',xlo(i),xlo(i)*ev2kayser
              stop 'hydropac, stopping'
            endif
          endif
          if (nup(i).le.modnlev) then
            ediff = abs(xup(i)*ev2kayser - modenergy(nup(i)))/
     &               modenergy(nup(i))
            if (ediff.gt.1.e-3) then
              print*,'H level identification problem'
              print*,'model atom ', modenergy(nup(i)), 
     &               'line list',xup(i),xup(i)*ev2kayser
              stop 'hydropac, stopping'
            endif
          endif
        endif
!
      enddo
*
! add opacity
      do k=1,ntau
        ne(k)=pe(k)/(t(k)*kbcgs)
        nh1(k)=sngl(presneutral(k,1))/(t(k)*kbcgs)
        nhe1(k)=sngl(presneutral(k,2))/(t(k)*kbcgs)
* hckt=hc/kT ready for lambda in AA.
        hckt(k)=hck/T(k)
        dopple(k)=sqrt( xi(k)**2 * 1. + 
     &                 2.*kbcgs*t(k)/1.6738e-24) /
     &                 sngl(clight*1.d2)
        fact = 1.0 / ro(k) / xkapr(k)
! departure coefficients for that depth
        if (nlte_species) then
          do i=1,modnlev
            bd(i)=b_departure(k,i)
          enddo
          do i=max(modnlev+1,1),nlevlist
            bd(i)=1.0
          enddo
        else
          bd=1.0
        endif
! loop over wavelengths
        do l=1,maxlam
          xlsingle=sngl(xlambda(l))
          expcorr = exp(hckt(k)/xlsingle)
          planckfct=bpl(T(k),xlsingle)
          if (nlte) source_f = source_function(k,l)
          call hbop(xlambda(l),nline,nlo,nup,hlambda,
     &     nh1(k),nhe1(k),ne(k),t(k),dopple(k),npop,
     &     bd,modnlev,expcorr,fact,source_f,
     &     planckfct,
     &     total,cont,contonly,lineonly,nlte,nlte_species)
! HI bf is already included in babsma.f
          contrib = (total - cont) * fact
          abso(k,l) = abso(k,l) + contrib
! save modified source function
          if (nlte) source_function(k,l) = source_f
        enddo
      enddo

      return
      end
