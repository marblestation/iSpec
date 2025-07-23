	SUBROUTINE funcv(fvec,fmin,pgin,pg)

C SUBROUTINe of die_pe.
*   computes fvec, for Newton iteration in die_pe.
*            fmin is 0.5*fvec*fvec
C 
C 03/06-1999 BPz
C
        implicit none

        integer maxim
        parameter (maxim=1000)
        integer i,j,m,mmax(maxim),mmaxj,natom(5,maxim),nelem(5,maxim),
     &          natomj,nelemj,nelemx(100),nmetal,nmol,nimax,ntry,niter,
     &          k,km5,nelemi,n,nelmolec(maxim),natomi,nhemolec(maxim)
        integer indx(100),ii,jj
        character*20   MOL(maxim)
        doubleprecision fvec(100),fmin
        character*20 molcode(maxim)
* really real, not supposed to become doubleprecision
        real spe,tem,spg,g0(100),g1(100),g2(100),g3(100)
        real exponent(maxim)
* may become dbleprec
	doubleprecision IP(100),KP(100),uiidui(100),eps,tem25,
     &       delta,deltae,pg,pgold,fpold(100),fifi,epsp,
     &       ppmol(maxim),apm(maxim),p(100),fp(100),c(maxim,5),
     &       ccomp(100),maxerror,pstart(100),fp99back,
     &       smallest,epsdie,econst,pgin,dhh,t,pmoljl,dfp99dp1,
     &       heh,aplogj,pph,phh,pgin5,reducedmass15(maxim),
     &       factor,pmolj,atomj,prev(100),pold(100),dfdp(100),
     &       kp1(100),kp2(100),kp3(100),kpe(100),ipp(100),ippp(100),
     &       d00(maxim),qmol(maxim),ipppp(100)
        logical switer,readjust,first,elmolec(maxim),converge,
     &          hemolec(maxim),molions
	COMMON/COMFH1/C,NELEM,NATOM,MMAX,PPMOL,d00,qmol,APM,MOL,IP,
     &              ipp,ippp,g0,g1,g2,g3,CCOMP,exponent,reducedmass15,
     &              UIIDUI,P,FP,KP,EPS,nelemx,nimax,NMETAL,NMOL,switer,
     &              molcode,elem
        character*2 elem(100)
        double precision  pe
        common/funco/kpe,dfdp,pe,elmolec,molions
*
      pgin5=pgin*1.d5
      molions=.true.
      readjust=.false.
1     continue

      pg=0.0
      fp(99)=0.

      do i=1,nmetal
        nelemi=nelemx(i)

* compute part of fictitious pressure fp(nelemi) due to p(nelemi), dfdp=dfp/dp
* compute also contribution to Pgas and to electron pressure

        dfdp(nelemi)=1.+kp(nelemi)
        fp(nelemi)=p(nelemi)*dfdp(nelemi)
        pg=pg+fp(nelemi)
        fp(99)=fp(99)-p(nelemi)*kpe(nelemi)
      enddo

      do J=1,nmol
* log(Pmol) = -logKp + sum(ni * log Pi)
        mmaxj=mmax(j)
        pmolj=1./apm(j)
        do m=1,mmaxj
          nelemj=nelem(m,j)
          natomj=natom(m,j)
          pmolj=pmolj*p(nelemj)**natomj

* add contribution of molecules to the fictitious pressure.
* fp(99) is -pe_calc until the end of the loop, where it changes sign and
* becomes pe_calc. 

        enddo
        do m=1,mmaxj
          nelemj=nelem(m,j)
          natomj=natom(m,j)
          if (molions.or.nelemj.ne.99) then
cc          if (molions) then
            fp(nelemj)=fp(nelemj)+pmolj*natomj
          else
********************
* artificially suppress contribution of molecular ions to 
* electron and atomic partial pressures.
********************
cc            print*,'funcv: electron contrib ',mol(j),pmolj
******
* try rescaling p(nelemi) of the negative contribution to fp99 that
* is larger than fp99 itself, if there is one. It would be better
* to search the largest negative contribution, rescale it, recompute fp99,
* and if negative again, search the next largest contribution, etc until
* fp99>0.
*****
            if (pmolj*natomj.gt.-(fp99back/2.)) then
              readjust=.true.
            endif
          endif
        enddo
        if (readjust) then
          do m=1,mmaxj
            nelemj=nelem(m,j)
            if (nelemj.ne.99) then
              natomj=natom(m,j)
*              print*,'funcv readjust 99',mol(j),(pmolj),
*     &         (fp99back),nelemj,(p(nelemj))
              p(nelemj)=p(nelemj)*1.d-5
            endif
          enddo
          readjust=.false.
          goto 1
        endif
        pg=pg+pmolj
        ppmol(j)=pmolj

        if (pmolj.gt.pgin5) then
*          print*,'funcv: P(',mol(j),')/Pgas5 =',sngl(pmolj/pgin5),
*     &      'readjusting'
          do m=1,mmaxj
              nelemj=nelem(m,j)
              natomj=natom(m,j)
              if (nelemj.ne.99) p(nelemj)=p(nelemj)*1.d-5
          enddo
          goto 1
        endif
*
      enddo
      fp(99)=-fp(99)
      fp99back=fp(99)
      pg=pg+fp(99)

******************************************************************************
      if (fp(99).le.0.d0) then
* try without molecular ions first
        if (.not.molions) then
          print*,'funcv: fp99=',fp(99)
          print*,'funcv: negative fp99 without molecular ions!'
          print*,'funcv: individual contributions to electron pressure:'
          do i=1,nmetal
             print*,nelemx(i),p(nelemx(i))*kpe(nelemx(i))
          enddo
          STOP 'funcv: negative fp99 without molecular ions!'
        else
          molions=.false.
*          print*,'funcv: fp99=',fp(99),
*     &         'trying without molecular ions'
        endif
        goto 1
      endif
******************************************************************************
*
* now, we have the fictitious pressures, and we need to compute
* fvec(i)=ln(fictp(i))-ln(fictp(He))-ln(abi), with abi=ccomp(i)/ccomp(He),
* the abundance of i relative to He. We take He, because it is
* involved in few molecules and hopefully the convergence should be 
* more stable than with H as a reference. 
* Fvec(i) should be 0, and the Newton-Raphson iteration is trying
* to achieve that, using the fjac(i,j) to correct the error.
* We need the Jacobian fjac(i,j)=dfvec(i)/dlnpj
*
* Equation (2) is for electron pressure, as there is no equation
* for He.
* fvec(2)= ln(pe_calc/pe)
* and fjac(2,j)= dlnpe_calc/dlnpj
*
* the correction is :  P_corrected(i)=p_prev(i)*exp(pp(i))
*

      do i=1,nmetal
        nelemi=nelemx(i)
cc            print*,' atom ',nelemi,' fictP ',fp(nelemi),' abund/He',
cc     &        ccomp(nelemi)/ccomp(2),' fp/fpHe ',fp(nelemi)/fp(2)

          fvec(i)=log(fp(nelemi))-log(fp(2))+
     &            log(ccomp(2)/ccomp(nelemi))
      enddo
* instead of Helium we set the equation for electronic pressure.
* max (fp99, 1.d.-100) in order to avoid log(<0)

      fvec(2)=log(max(fp(99),1.d-100)/pe)

CC the above solution crashes when fp(99) becomes negative...
CC because the perturbation in fvec to compute fjac has all (2,j)
CC components null.
CC
CC But that one does not work either (roundoff problem in die....)
CC            fvec(2)=fp(99)/pe-1.d0

      fmin=0.
      do i=1,nmetal
        fmin=fmin+fvec(i)**2
      enddo
      fmin=fmin*0.5

************
**        print*,'********************** FVEC ***************************'
**        print*,fvec
**        print*,'********************** FVEC ***************************'
************

      return
      end
