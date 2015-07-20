      subroutine getlele(iel,ion,lele)
*
      implicit none
      integer j,iel,maxim,ion
      integer*8 ielprop,ionprop
      parameter (maxim=1000)
      character*20 lele,molchar,ielchar,ionchar
      real exponent(maxim),g0(100),g1(100),g2(100),g3(100)
      character mol(maxim)*20,elem(100)*2
      logical switer,found(maxim),first
      integer nelem(5,maxim),natom(5,maxim),mmax(maxim),nimax,
     &        nelemx(100),nmetal,nmol
* may become dbleprec
      doubleprecision ip(100),kp(100),uiidui(100),eps,fp(100),
     &     ppmol(maxim),apm(maxim),c(maxim,5),ccomp(100),p(100),
     &     ipp(100),ippp(100),d00(maxim),qmol(maxim),
     &     reducedmass15(maxim)
      character*20     molcode(maxim)
       common/comfh1/c,nelem,natom,mmax,ppmol,d00,qmol,apm,mol,ip,
     &            ipp,ippp,g0,g1,g2,g3,ccomp,exponent,reducedmass15,
     &            uiidui,p,fp,kp,eps,nelemx,nimax,nmetal,nmol,switer,
     &            molcode,elem
       data first/.true./

* One call to molecpartf to initialize arrays (some molecules sorted out).
      if (first) then
        call molecpartf(2000.,found)
        first=.false.
      endif

      lele='                    '
      if (iel.gt.92) then
        do j=1,nmol
          molchar=molcode(j)
          if (molchar.ne.'                    ') then
            ielchar=molchar(1:index(molchar,'.')-1)
            ionchar=molchar(index(molchar,'.')+1:index(molchar,'.')+2)
            read(ielchar,*) ielprop
            read(ionchar,*) ionprop
            ionprop=ionprop+1
            if (ielprop.eq.iel.and.ionprop.eq.ion) then
              lele=mol(j)
              print*,'getlele ',lele
              return
            endif
          endif
        enddo
      else
        print*,
     &   'getlele does not provide identification for atomic species'
        print*,'ERROR in getlele'
        stop
      endif
      if (lele.eq.'                    ') then
        print*,'ERROR ! getlele: molecule not found:',iel
        stop 'check Irwin molecular datafile and line list!'
      endif

      return
      end
