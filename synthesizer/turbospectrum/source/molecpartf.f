      subroutine molecpartf(temp,found)
*
* reads in Irwin and Sauval molecular data, select the relevant
* molecules and computes partition functions. 
*  BPz 21/08-1996
*
      implicit none
      integer maxim
      parameter (maxim=1000)
      integer i,j,n,m,mmaxj,nelemj,natomj,jhm,jj,k
      real temp,exponent(maxim),g0(100),g1(100),g2(100),g3(100)
      doubleprecision qcoeff(10,maxim),qc(10),qq,logt
      character mol(maxim)*20,molname*20,ligne*10
      logical first,switer,found(maxim),ilenreste,sfound(maxim)
      integer nelem(5,maxim),natom(5,maxim),mmax(maxim),nimax,
     &        nelemx(100),
     &        nmetal,nmol,idumb,nqcoeff(maxim),nq
* may become dbleprec
      doubleprecision ip(100),kp(100),uiidui(100),eps,fp(100),
     &     ppmol(maxim),apm(maxim),c(maxim,5),ccomp(100),p(100),
     &     ipp(100),ippp(100),d00(maxim),qmol(maxim),
     &     dissoc,reducedmass15(maxim),ndensity,molweight,d00hm
      character*20 molcod,molcode(maxim)
       common/comfh1/c,nelem,natom,mmax,ppmol,d00,qmol,apm,mol,ip,
     &            ipp,ippp,g0,g1,g2,g3,ccomp,exponent,reducedmass15,
     &            uiidui,p,fp,kp,eps,nelemx,nimax,nmetal,nmol,switer,
     &            molcode,elem
      character*2 elem(100)
      integer natomm(5),nelemm(5)
      real xmass(maxim+400),atmass(100)
      common /density/ndensity,molweight,xmass,atmass
      real abund,anjon,h,part,dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy
      COMMON/CI5/ABUND(16),ANJON(16,5),H(5),PART(16,5),DXI,
     & F1,F2,F3,F4,F5,XKHM,XMH,XMY

      data first /.true./

      save qcoeff,nqcoeff,first
 
      if (first) then
        do j=1,nmol
          found(j)=.false.
        enddo
* read molecular data and sort according to list mol
        open(55,file=
     &    'DATA/IRWIN_molecules_v15.1.dat',
     &       status='old')
123     read(55,'(a)') ligne
        if (ligne(1:1).eq.'*') then
          goto 123
        else
          backspace(55)
        endif
        do i=1,100000
          read(55,*,end=99) idumb, molcod,molname,dissoc,nq,
     &               (qc(n),n=1,nq),
     &               mmaxj,(nelemm(m),natomm(m),m=1,mmaxj)
          if (nq.gt.10) then
            print*, 'molecpartf: nq=',nq,' too large!'
            stop
          endif
          j=1
          ilenreste=.true.
          do while (ilenreste)
            if (molname.eq.mol(j)) then
              if (dissoc.ne.0..or.
     &             c(j,1).ne.0..or.c(j,2).ne.0..or.c(j,3).ne.0..or.
     &              c(j,4).ne.0..or.c(j,5).ne.0.) then
* either Irwin's data is complete (Q + D00) or if D00=0. (unknown) there
* is available Tsuji data (Kp is defined, i.e. the coefficients are not 0.)
                found(j)=.true.
                molcode(j)=molcod
                d00(j)=-dissoc
                mmax(j)=mmaxj
                do m=1,mmaxj
                  nelem(m,j)=nelemm(m)
                  natom(m,j)=natomm(m)
                enddo
                do n=1,nq
                  qcoeff(n,j)=qc(n)
                enddo
                nqcoeff(j)=nq
              else
* there is no Tsuji data and D00 is unknown in Irwin's file.
* Suppress this molecule and rearrange list.
              print*,'MOLECPARTF. suppressing:',j,' ',mol(j),
     &                '(unknown D00)'
                do jj=j,nmol-1
                  mol(jj)=mol(jj+1)
                  molcode(jj)=molcode(jj+1)
                  d00(jj)=d00(jj+1)
                  found(jj)=found(jj+1)
                  nqcoeff(jj)=nqcoeff(jj+1)
                  do n=1,nqcoeff(jj)
                    qcoeff(n,jj)=qcoeff(n,jj+1)
                  enddo
                  do k=1,5
                    c(jj,k)=c(jj+1,k)
                  enddo
                  mmax(jj)=mmax(jj+1)
                  do m=1,mmax(jj)
                    nelem(m,jj)=nelem(m,jj+1)
                    natom(m,jj)=natom(m,jj+1)
                  enddo
                enddo
                nmol=nmol-1
              endif
            endif
            j=j+1
            if (j.gt.nmol) ilenreste=.false.
          enddo
        enddo
99      continue
        do j=1,nmol
          sfound(j)=found(j)
        enddo
        do j=1,nmol
          if (mol(j).eq.'H -                 ') then
            print*,'H- electron affinity will be changed by -2dxi'
            d00hm=d00(j)
            jhm=j
          endif
          xmass(j+4*nmetal)=0.
          exponent(j)=-1
          reducedmass15(j)=1.
          mmaxj=mmax(j)
          do m=1,mmaxj
            nelemj=nelem(m,j)
            natomj=natom(m,j)
            xmass(j+4*nmetal)=xmass(j+4*nmetal)+
     &                        natomj*atmass(nelemj)
            exponent(j)=exponent(j)+natomj
            reducedmass15(j)=reducedmass15(j)*atmass(nelemj)**natomj
          enddo
          reducedmass15(j)=(reducedmass15(j)/xmass(j+4*nmetal))**1.5
          if (.not.found(j)) then
cc            print*,'Molecule ',mol(j),' not found in data file!!'
            if (first) then
              print*,' data may be added in IRWIN_molecule file'
              print*,' to update older data from Tsuji'
              first=.false.
            endif
* note that molecules not appearing in IRWIN_molecules file do not
* have a molcode. This makes them unsuitable for use in getlele.f
* and hence in line calculations.
            print*,'using Tsuji''s old data for ',mol(j)
            d00(j)=0.00
            nqcoeff(j)=1
            qcoeff(1,j)= -1.d300
          endif
        enddo
        close(55)
        first=.false.
      endif

* first correct e-affinity of H- ; note that unlike other molecules, d00
* does not contain the original energy.

      d00(jhm)=d00hm-2.*dxi

      logt=log(temp)
      do j=1,nmol
        qq=qcoeff(nqcoeff(j),j)
        do n=nqcoeff(j)-1,1,-1
          qq=qq*logt + qcoeff(n,j)
        enddo
        qmol(j)=exp(qq)
        found(j)=sfound(j)
      enddo

      return
      end
