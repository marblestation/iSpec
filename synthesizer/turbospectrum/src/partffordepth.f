      subroutine partffordepth(ntau,temp,molname,fpartition)
*
      implicit none
      integer j,ntau,k,molindex,maxim
      parameter (maxim=1000)
      real temp(ntau),exponent(maxim),g0(100),g1(100),g2(100),g3(100),
     &          fpartition(ntau)
      doubleprecision qcoeff(10,maxim),qc(10),qq,logt
      doubleprecision dpartf,dtemp
      character mol(maxim)*20,molname*20
      logical switer,found(maxim),Ames,first,scan2001,oldscan,Barber
      integer nelem(5,maxim),natom(5,maxim),mmax(maxim),nimax,
     &        nelemx(100),
     &        nmetal,nmol,idumb,nqcoeff(maxim),nq
* may become dbleprec
      doubleprecision ip(100),kp(100),uiidui(100),eps,fp(100),
     &     ppmol(maxim),apm(maxim),c(maxim,5),ccomp(100),p(100),
     &     ipp(100),ippp(100),d00(maxim),qmol(maxim),
     &     dissoc,reducedmass15(maxim),ndensity,molweight,d00hm
      character*20     molcode(maxim)
       character*2 elem(100)
       common/comfh1/c,nelem,natom,mmax,ppmol,d00,qmol,apm,mol,ip,
     &            ipp,ippp,g0,g1,g2,g3,ccomp,exponent,reducedmass15,
     &            uiidui,p,fp,kp,eps,nelemx,nimax,nmetal,nmol,switer,
     &            molcode,elem
       common/h2ochoice/Ames,scan2001,oldscan,Barber
       data first/.true./


      if (molname.eq.'H O H               ') then 
        if (Ames) then
          print*,'partffordepth: using Ames partf.'
          do k=1,ntau
            dtemp=temp(k)
            call partfAmesH2O(dtemp,dpartf)
            fpartition(k)=dpartf
          enddo
        else if (scan2001) then
          print*,'Using SCAN 2001 partition function for H2O'
          call h2opartf2001(ntau,temp,fpartition)
        else if (oldscan) then
          print*,'partffordepth: using uffe partf.'
          call h2opartf(ntau,temp,fpartition)
        else if (Barber) then
          do k=1,ntau
            dtemp=temp(k)
            call partfBarberH2O(dtemp,dpartf)
            fpartition(k)=dpartf
          enddo
        else
          print*,'No partition function for the H2O line list!'
          stop
        endif
*
        return
*
      endif

* One call to molecpartf to initiate arrays (some molecules sorted out).
      if (first) then
        call molecpartf(temp(1),found)
        first=.false.
      endif

      do j=1,nmol
        if (molname.eq.mol(j)) then
          molindex=j
          goto 1
        endif
      enddo
1     continue

      do k=1,ntau
        call molecpartf(temp(k),found)
        if (.not.found(molindex)) stop 
     &            'partffordepth: molecule not found'
        fpartition(k)=qmol(molindex)
      enddo

      return
      end
