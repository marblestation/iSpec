      subroutine hydropac(lunit,xlboff)
!
! new version for use with hbop.f  ! H line + continuum bf opacity
! Bpz 02/02-2019

      implicit none

      include 'spectrum.inc'

      character comment*100,species*20
      integer ion,nline,ntau,iter,l,k,i,lstart,lstartp,
     &        nmy,nlbldu,iint,iweak,lunit,maxlam,ll,
     ,        jlcont, nlcont
      doubleprecision ionpot,wave,xl1,xl2,del,xlboff,xlambda

! 150 is the max allowed number of H lines
      doubleprecision hlambda(150)
      real xlo(150),xup(150),gf(150),npop(150)
      real cont,total
      integer nlo(150),nup(150)
      character*9 lname(150)
!
      real fpart(ndp),pe,t,pg,xi,mum,ro,eps,xmyc,
     &     ne(ndp),
     &     nh1(ndp),nhe1(ndp),hnorm,hnormnostim(ndp),stim,
     &     diff,diffp,theta(ndp),hckt(ndp),abso,absos,absocont,
     &     dopple(ndp),contrib,xkapr,cross,hlinop,h1bfg,alpha,
     ,     ee(6),d0(ndp),d(lpoint,ndp),absoscont
      logical  lymanalpha, usedam,notfound,contonly
*      
      common /atmos/ t(ndp),pe(ndp),pg(ndp),xi(ndp),mum(ndp),ro(ndp),
     &               ntau
      doubleprecision presneutral,presion,presion2,presion3, h1bfgc(30)    
      common /orderedpress/ presneutral(ndp,100),presion(ndp,100),
     &                      presion2(ndp,100),presion3(ndp,100)
      common /pieces/ xl1,xl2,del,eps,nmy,nlbldu,iint,xmyc,iweak
      common /large/ xlambda(lpoint),maxlam,abso(ndp,lpoint),
     &               absos(ndp,lpoint),absocont(ndp,lpoint),
     &               absoscont(ndp,lpoint)
      common/rossc/ xkapr(ndp),cross(ndp)
      common/continuum/nlcont,jlcont(lpoint)
      
      character*50 zMC
      integer znlc ,  ii
      real    zxlm,zBP(ndp),zXC(ndp),zS(ndp),zXI(ndp)

! set populations to 0. (computed as LTE pop in hbop)
      npop=0.0
! compute lines and continuum
      contonly=.false.
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
      do i=1,nline
!       oscillator strengths are computed in hbop.f
        read(lunit,*) hlambda(i),nlo(i),nup(i),xlo(i),xup(i),gf(i),
     &                lname(i)
      enddo
*
      do k=1,ntau
        ne(k)=pe(k)/(t(k)*1.38066e-16)
        nh1(k)=sngl(presneutral(k,1))/(t(k)*1.38066e-16)
        nhe1(k)=sngl(presneutral(k,2))/(t(k)*1.38066e-16)
* hckt=hc/kT ready for lambda in AA.
        hckt(k)=2.99792458e10*6.626075e-27/1.380658e-16*1.e8/t(k)
        dopple(k)=sqrt( xi(k)**2 * 1. + 
     &                 2.*1.380658e-16*t(k)/1.6738e-24) /
     &                 2.99792458e10
      enddo
! add opacity
      do k=1,ntau
        do l=1,maxlam
*	  print*,' l,xlam ',l,xlambda(l)
          call hbop(xlambda(l),nline,nlo,nup,hlambda,
     &         nh1(k),nhe1(k),ne(k),t(k),dopple(k),npop,0,total,cont,
     &         contonly)
! HI bf is already included in babsma.f
          abso(k,l) = abso(k,l) + (total - cont)/xkapr(k)/ro(k)
!!!!          absocont(k,l) = absocont(k,l) + cont/xkapr(k)/ro(k)

!!!!!!!!!          if(contrib/absocont(k,l) .le. eps) goto 5    !  in order to implement a condition for not calling hbop, we need to start from line centers 
        enddo      
   5    continue
      enddo

      return
      end
