      SUBROUTINE DIE_pe(TEM,SPE,spg,found,converge,niter,skiprelim)

C DIE_pe SUBROUTINE DE EQMOL_pe
C 
C based on cleaned-up die.f
C uses T,Pe as input instead of T,Pg.
C
C 28/08-1996 BPz
C
        implicit none

        integer maxim
        parameter (maxim=1000)
        integer i,j,m,mmax(maxim),mmaxj,natom(5,maxim),nelem(5,maxim),
     &          natomj,nelemj,nelemx(100),nmetal,nmol,nimax,ntry,niter,
     &          k,km5,nelemi,n,nelmolec(maxim),natomi,nhemolec(maxim)
        integer indx(100),ii,jj
        character*20   MOL(maxim)
        real xmass(maxim+400), atmass(100),exponent(maxim)
        doubleprecision ndensity, molweight,tsuapm,fmin,fminold,
     &                  fmin2,fminold2,alam,alam2,tmplam,disc,slope,
     &                  rhs1,rhs2,a,b,g(100),sum,stpmax,alamin,
     &                  maxcorr,mincorr,perturb
        doubleprecision dd,pp(100),fvec(100),fjac(100,100),errx,errf,
     &                  fvecold(100),plncorr(100),pback,pgdumm,fmindumm
* really real, not supposed to become doubleprecision
        real spe,tem,spg,g0(100),g1(100),g2(100),g3(100)
* may become dbleprec
        doubleprecision IP(100),KP(100),uiidui(100),eps,tem25,
     &       delta,deltae,pg,pgold,fpold(100),fifi,epsp,pgindumm,
     &       ppmol(maxim),apm(maxim),p(100),fp(100),c(maxim,5),
     &       ccomp(100),maxerror,pstart(100),epsf,fp99back,
     &       smallest,epsdie,econst,pgin,dhh,t,pmoljl,dfp99dp1,
     &       heh,aplogj,pph,phh,pgin5,reducedmass15(maxim),
     &       factor,pmolj,atomj,prev(100),pold(100),dfdp(100),
     &       kp1(100),kp2(100),kp3(100),kpe(100),ipp(100),ippp(100),
     &       d00(maxim),qmol(maxim),ipppp(100)
        character*20 molcode(maxim)
        logical switer,readjust,first,elmolec(maxim),converge,
     &          hemolec(maxim),molions,singular,skiprelim
        common /density/ndensity,molweight,xmass,atmass
        COMMON/COMFH1/C,NELEM,NATOM,MMAX,PPMOL,d00,qmol,APM,MOL,IP,
     &              ipp,ippp,g0,g1,g2,g3,CCOMP,exponent,reducedmass15,
     &              UIIDUI,P,FP,KP,EPS,nelemx,nimax,NMETAL,NMOL,switer,
     &              molcode,elem
        character*2 elem(100)
        double precision  error(100),pe,xx,bpzlog10,pep
        double precision konst,qprod,konsttem25
        real abund,anjon,h,part,dxi,f1,f2,f3,f4,f5,xkhm,xmh,xmy
        common /ci5/abund(16),anjon(16,5),h(5),part(16,5),dxi,
     &              f1,f2,f3,f4,f5,xkhm,xmh,xmy
        logical tsuswitch,tsuji,found(maxim)
        integer nattsuji,nmotsuji
        doubleprecision parptsuji
        common /tsuji/ tsuji,tsuswitch,nattsuji,nmotsuji,
     &                 parptsuji(maxim+400)
        common/funco/ kpe,dfdp,pe,elmolec,molions
        data epsp /0.05/
* epsp is the convergence criterium for a rough equilibrium, without Pe
* being checked. It only insures that the individual atomic pressures 
* are not changing by a larger amount anymore.
        data first/.true./
        save first
* this save statement may not be complete. I have not checked if anything else
* needs saving.
*
      konst=25947.172
*konst=(2 pi amu / h^2)^1.5 * k^2.5
*
      if (first) then
        smallest = tiny(parptsuji)
      endif

      ECONST=4.342945E-1
***** TEST : skip preliminary equilibrium iteration.
cc      if (tem.lt.4000.) then
cc        epsdie=1.e30
cc      else
        epsdie=0.1
cc      endif

      converge=.true.
      pe=spe
      p(99)=pe
      T=5040.0/TEM
      pgin=spg
      if (pgin.le.0.) pgin=1.
C
C heh=helium/hydrogene ratio by number
      HEH=CCOMP(2)/CCOMP(1)
C
* backup kp(H2) 
      DHH=(((0.1196952E-02*T-0.2125713E-01)*T+0.1545253E+00)*T
     &    -0.5161452E+01)*T+0.1277356E+02
      DHH=EXP(DHH/ECONST)
C
C  EVALUATION OF THE IONIZATION CONSTANTS
      TEM25=TEM**2*SQRT(TEM)
      konsttem25=konst*tem25
      g0(99)=2.

* security. This must be fulfilled in order for the Newton-Raphson
* iteration to work .

      if (nelemx(2).ne.2) stop 
     & 'STOP in die_pe !! Helium MUST be second in the list of atoms!!'

****

      DO I=1,NMETAL
        NELEMI = NELEMX(I)
*
* calculation of the partition functions following Irwin (1981)
        call partf(nelemi,1,tem,1,g0(nelemi),ip(nelemi))
        call partf(nelemi,2,tem,1,g1(nelemi),ipp(nelemi))
        call partf(nelemi,3,tem,1,g2(nelemi),ippp(nelemi))
        call partf(nelemi,4,tem,1,g3(nelemi),ipppp(nelemi))
        uiidui(nelemi)=g1(nelemi)/g0(nelemi)*0.6667654
*
*****
* We include neutral, 1st, 2nd and 3rd ion.
* kp1 = pe * p+ / p
* kp2 = pe * p++ / p+
* kp3 = pe * p+++ / p++
* kp is  (p+ +  p++ +  p+++) / p
* kpe is (p+ + 2p++ + 3p+++) / p
*        (ie the electron contribution of the species / p)
*****
        if (nelemi.eq.1) then
          KP1(NELEMI) =UIIDUI(NELEMI)*TEM25*
     &                 EXP((dxi-IP(NELEMI))*T/ECONST)
          kp2(NELEMI) =g2(nelemi)/g1(nelemi)*0.6667654*TEM25*
     &                 EXP((2.*dxi-ipp(NELEMI))*T/ECONST)
          kp3(NELEMI) =0.0
        else
          KP1(NELEMI) =UIIDUI(NELEMI)*TEM25*
     &                 EXP((dxi-IP(NELEMI))*T/ECONST)
          kp2(NELEMI) =g2(nelemi)/g1(nelemi)*0.6667654*TEM25*
     &                 EXP((2.*dxi-ipp(NELEMI))*T/ECONST)
          kp3(NELEMI) =g3(nelemi)/g2(nelemi)*0.6667654*TEM25*
     &                 EXP((3.*dxi-ippp(NELEMI))*T/ECONST)
        endif
        kp(nelemi)=kp1(nelemi)/pe *
     &             (1. + kp2(nelemi)/pe *
     &                   (1. + kp3(nelemi)/pe))
        kpe(nelemi)=kp1(nelemi)/pe *
     &             (1. + 2.* kp2(nelemi)/pe *
     &                   (1. + 1.5 * kp3(nelemi)/pe))
      enddo
****************************************************
C evaluation of kp(mol) ==apm(mol)
      do J=1,NMOL
* if we have partition function data, use it.
        if (found(j)) then
          if (d00(j).eq.0.0) then
* but not if D00 is missing , instead we use kp from old Tsuji's data
            tsuapm=c(j,5)
            do k=4,1,-1
              tsuapm=tsuapm*t+c(j,k)
            enddo
            apm(j)=exp(tsuapm/econst)
            if (first) then
              print*,'WARNING ! D00(',mol(j),')=',d00(j)
              print*,'using old Tsuji''s data!'
            endif
          else
            qprod=1./qmol(j)
            mmaxj=mmax(j)
            do m=1,mmaxj
              nelemj=nelem(m,j)
              natomj=natom(m,j)
              qprod=qprod*g0(nelemj)**natomj
            enddo
            APM(J)=(konsttem25)**exponent(j)*qprod*
     &              reducedmass15(j)*exp(-d00(j)*t/econst)
          endif
        else
* We use kp from old Tsuji's data
          tsuapm=c(j,5)
          do k=4,1,-1
            tsuapm=tsuapm*t+c(j,k)
          enddo
          apm(j)=exp(tsuapm/econst)
        endif
      enddo

      if (.not.skiprelim) then
****************************************************
*
*   INCLUDE CO and C2H2 and ????
*
*********************************************************
C
C   PRELIMINARY VALUE OF Pg and fph assuming Pg=PH+PHH+PPH+Pe+PHe 
C use first pg=pgin (guess from outside, or pg=1. then iterate on
C atoms and H2 only, with low accuracy (10%) to get a reasonable pg.
C as fictP(H) = PH+2*PHH+PPH, we get:
C PG = PH+PHH+PPH + pe + HEH*(PH+2.0*PHH+PPH)
* approximate fph for starting:
      ntry=0
1     if (tem.lt.2000.) then
* at low T, pg approx = ph2 + phe = fph/2. + fph*heh
        fp(1)=pgin/(0.5+heh)
      else
        fp(1)=pgin/(1.+heh)
      endif
      p(1)=dhh/4.* (sqrt((1.+kp(1))**2+8.*fp(1)/dhh) - 
     &  (1.+kp(1)))
      do niter=1,100
        pph=p(1)*kp(1)
        phh=p(1)**2/dhh
        pg=p(1)+phh+pph+pe
        fp(99)=pph
        fp(1)=p(1)+pph+2.*phh
        dfp99dp1=kp(1)*ccomp(1)
        dfdp(1)=(1.+kp(1))+4.*p(1)/dhh
ccc        print*,p(1),phh,pph,dfdp(1)
C
C evaluation of the fictitious pressure of each element, fp, of its
C pressure (neutral) of the error on fp(nelemi)/fph (should be equal
C to ccomp), and of the derivative of fp/fph relative to p(nelemi).
C these are only first estimates. In the iteration loop below, values
C from the previous iteration are used.
C p(nelemi) estimated assuming only neutral and first ion present
C (no molecules nor higher ionisation stages).
C higher ionisation stages easy to implement.
        do I=2,NMETAL
          NELEMI=NELEMX(I)
          FP(NELEMI)=CCOMP(NELEMI)*fp(1)
          p(nelemi)=fp(nelemi)/(1.+kp(nelemi))
          fp(99)=fp(99)+fp(nelemi)*kpe(nelemi)/(1.+kp(nelemi))
          dfdp(nelemi)=1.+kp(nelemi)
          dfp99dp1=dfp99dp1+kpe(nelemi)*ccomp(nelemi)*
     &             dfdp(1)/dfdp(nelemi)
          pg=pg+fp(nelemi)
        enddo
        error(99)=pe-fp(99)
ccc        print*,'error',error(99), dfp99dp1
        if (abs(error(99)/pe).lt.epsdie) goto 100
        p(1)=p(1)+error(99)/dfp99dp1
ccc	WRITE(6,6109) TEM,PG,fp(1),P(1)
6109    FORMAT('TEM=',F10.2,1X,'PG=',E13.5,1X,
     &  'FPH=',E12.3,1X,'PH=',E12.3)
      enddo
100   continue
      if (niter.ge.100) then
**** TEST ***********************************************
cc        print*,'WARNING PRELIMINARY ATOMIC ',
cc     &  'EQUILIBRIUM NOT CONVERGED AFTER 100 iterations'
*********************************************************
        if (tem.lt.2000.) then
          fp(1)=pgin/(0.5+heh)
        else
          fp(1)=pgin/(1.+heh)
        endif
        do i=2,nmetal
          nelemi=nelemx(i)
          fp(nelemi)=ccomp(nelemi)*fp(1)
        enddo
        p(1)=dhh/4.* (sqrt((1.+kp(1))**2+8.*fp(1)/dhh) -
     &    (1.+kp(1)))
        pgin=1.d10
      else
        pgin=pg
        pold(1)=p(1)
*        print*,'PRELIMINARY EQUILIBRIUM converged after',niter,
*     &         'iterations'
*        print*,'atoms and H2 only: '
*        print*,'t pe pg ph ph2 rel. error on pe ',
*     &    tem, sngl(pe), sngl(pg), sngl(p(1)), sngl(phh), 
*     &    sngl(error(99)/pe)
      endif
 
*
* try to avoid negative fp(99). Scale down inital pressures.
********************
      if (tem.le.2.) then
        do i=1,nmetal
          nelemi=nelemx(i)
          if (nelemi.ne.1.and.nelemi.ne.2) then
c          if (nelemi.ne.2) then
            p(nelemi)=fp(nelemi)*exp(-2.*t/econst)
            if (tem.le.1500.) p(nelemi)=fp(nelemi)*exp(-5.*t/econst)
            if (tem.le.900.) p(nelemi)=fp(nelemi)*exp(-10.*t/econst)
          endif
          p(nelemi)=max(p(nelemi),smallest)
          pstart(nelemi)=p(nelemi)
        enddo
      endif
      P(100)=PPH
      pgold=pgin

      else
*******************
* we end up here if we skipped the preliminary equilibrium
        pgin=pg
        do i=1,nmetal
          p(nelemx(i))=p(nelemx(i))*(pe/pep)
        enddo
      endif
C
C    RUSSELL EQUATIONS

      if (first) then
* find out molecules with electrons or with He
        do j=1,nmol
          elmolec(j)=.false.
          hemolec(j)=.false.
          mmaxj=mmax(j)
          do m=1,mmaxj
            nelemj=nelem(m,j)
            natomj=natom(m,j)
            if (nelemj.eq.99) then
              elmolec(j)=.true.
              nelmolec(j)=natomj
*              print*,'DIE: molecule: ',mol(j),' will be used for Pe'
            else if (nelemj.eq.2) then
              hemolec(j)=.true.
              nhemolec(j)=natomj
*              print*,'DIE: molecule: ',mol(j),' will be used for He'
            endif
          enddo
        enddo
        first=.false.
      endif
*
* Start of the Newton-Raphson iteration on the neutral atomic pressures

      epsf=1.e-10
      pgindumm=1.d200

      do niter=1,nimax
        delta=0.0
        pg=0.0
        fp(99)=0.

        do i=1,100
          do j=1,100
            fjac(j,i)=0.
          enddo
        enddo

* compute fvec, fmin and pg from the P().
* pgin is used as a pressure limitator (cf. funcv.f)

        call funcv(fvec,fmin,pgin,pg)

cc        print*,'**********************************'
cc        print*,'fvec ',fvec
cc        print*,'**********************************'

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
* compute dfvec/dlnPi numerically.
*
        do i=1,nmetal
          fvecold(i)=fvec(i)
        enddo
        do j=1,nmetal
          pback=p(nelemx(j))
          perturb=abs(log(pback))*1.d-4
          if (perturb.eq.0.) perturb=1.d-4
          p(nelemx(j))=exp(log(pback)+perturb)
          perturb=log(p(nelemx(j)))-log(pback)
          call funcv(fvec,fmindumm,pgindumm,pgdumm)
          if (fp(99).lt.0.) then
*            print*,'die_pe is going to die... fp99=',fp(99)
          endif
cc             print*,'delta fvec',sngl(fvec-fvecold)
          do i=1,nmetal
            fjac(i,j)=(fvec(i)-fvecold(i))/perturb
          enddo
* if fp(99)<0, the security in funcv works and avoids crash there,
* but the result is that fjac computed in the loop is identically
* 0 for (2,j).  That's why fvec(2) is recomputed without the contribution
* from the molecular ions in funcv.
          p(nelemx(j))=pback
        enddo
        do i=1,nmetal
          fvec(i)=fvecold(i)
        enddo
*        print*,'die fvec2',tem,sngl(pe),sngl(fp(99)),sngl(fvec(2))
*********************************************************************
        goto 1234
        if (.not.molions) then
* compute fjac(2,i) analytically.
*
* fvec(2)= ln(pe_calc/pe)
* and fjac(2,j)= dlnpe_calc/dlnpj
*
        do i=1,nmetal
          nelemi=nelemx(i)
          fjac(2,i)=kpe(nelemi)
cc        print*,'fjac ',2,i,fjac(2,i)
        enddo

        do j=1,nmol
          if (elmolec(j)) then
            mmaxj=mmax(j)
            do m=1,mmaxj
              nelemj=nelem(m,j)
              natomj=natom(m,j)
              if (nelemj.ne.99) then
                do i=1,nmetal
                  if (nelemx(i).eq.nelemj) jj=i
                enddo
                fjac(2,jj)=fjac(2,jj)-
     &                nelmolec(j)*natomj*ppmol(j)/p(nelemj)
cc        print*,'fjac ',mol(j),2,jj,fjac(2,jj)
              endif
            enddo
          endif
        enddo

        do i=1,nmetal
          fjac(2,i)=fjac(2,i)*p(nelemx(i))/fp(99)
cc        print*,'fjac ',2,i,sngl(fjac(2,i)),' final'
        enddo

        endif
1234    continue
*********************************************************************

cc        print*,'Pe, pg ' ,fp(99),pg

        deltae=fp(99)/pe -1.

        errf=0.
        do i=1,nmetal
          errf=errf+abs(fvec(i))
        enddo
        if (errf.le.epsf) goto 6054
*
* epsf is much smaller than eps (convergence criterium for individual pressure),
*  as the fvec may be converged to eps, but the 
* p(i) may not be as well converged.... ?? hmmm?

        do i=1,nmetal
          pp(i)=-fvec(i)
        enddo

        do i=1,nmetal
          sum=0.
          do j=1,nmetal
            sum=sum+fjac(j,i)*fvec(j)
          enddo
          g(i)=sum
          pold(nelemx(i))=p(nelemx(i))
        enddo
        fminold=fmin

***************************************************

        call ludcmp(fjac,nmetal,100,indx,dd,singular)
        if (singular) then
          converge=.false.
          return
        endif
        call lubksb(fjac,nmetal,100,indx,pp)
***************************************************

        slope =0.
        do i=1,nmetal
          slope=slope+g(i)*pp(i)
        enddo
        if (slope.gt.0.) then
           print*,'roundoff problem in die_pe!! '
           print*,'die_pe WARNING, unconverged T=',
     &             tem,' Pe=',sngl(pe)
           converge=.false.
           return
        endif
        alam=1.0
cc        print*,'slope ',slope

***************************************************
        errx=0
        maxcorr=-1.d30
        mincorr=1.d30
        do i=1,nmetal
          nelemi=nelemx(i)
          if (abs(pp(i)).lt.5) then
          errx=errx+abs(exp(pp(i))-1.)
cc          print*,'iter ',niter,'  atom ',i, nelemi,' P ', 
cc     &       sngl(p(nelemi)),
cc     &          'corr factor ',sngl(exp(pp(i)))
          else
            errx=1.e10
          endif
          maxcorr=max(maxcorr,pp(i))
          mincorr=min(mincorr,pp(i))
        enddo
*
* limit correction to avoid overflows
* stpmax=11. represents approximately a change by a factor 1.e5 in P 
*
        stpmax=30.
        do i=1,nmetal
          pp(i)=min(stpmax,max(pp(i),-stpmax))
        enddo

cc        print*,' AFTER SCALING OF CORRECTIONS:'
        do i=1,nmetal
          nelemi=nelemx(i)
cc          print*,'iter ',niter,'  atom ',i, nelemi,' P ', 
cc     &       sngl(p(nelemi)),
cc     &          'corr factor ',sngl(exp(pp(i)))
        enddo
***************************************************
*        print*,'initial errx ',errx
*        print*,'initial f:', fmin

* BPz: just a guess. In Num. Rec. a computation of alamin is done...
        alamin=1.e-3
* ... like this:
        maxcorr=0.
        do i=1,nmetal
          maxcorr=max(abs(exp(pp(i)))/max(abs(p(nelemx(i))),1.d0),
     &                maxcorr)
        enddo
        alamin=eps/maxcorr
*        print*,'alamin= ',alamin

124     continue
        do i=1,nmetal
          nelemi=nelemx(i)
          p(nelemi)=pold(nelemi)*exp(alam*pp(i))
          plncorr(i)=alam*pp(i)
        enddo
        call funcv(fvec,fmin,pgindumm,pg)
*        print*,'alam ',alam,'fmin ',fmin
        if (alam.lt.alamin) then
*
* may not be converged, we try one more iteration.
* It may end up here if in asecondary minimum....
* we just try the full (damped) correction ...
*
*          print*,' SECONDARY MIN? '
          alam=1.
          do i=1,nmetal
            nelemi=nelemx(i)
*first: damp correction
            pp(i)=min(20.d0,max(pp(i),-20.d0))
cc            if (molions) then
              p(nelemi)=pold(nelemi)*exp(alam*pp(i))
              plncorr(i)=alam*pp(i)
cc            else 
cc* desesperate move... BPz 28/06-1999
cc              p(nelemi)=pold(nelemi)/1000.
cc            endif
          enddo
          goto 123
        else if (fmin.lt.fminold+1.d-4*alam*slope) then
          goto 123
        else
          if (alam.eq.1.) then
            tmplam=-slope/(2.*(fmin-fminold-slope)) 
          else 
            rhs1=fmin-fminold-alam*slope 
            rhs2=fmin2-fminold2-alam2*slope 
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2) 
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2) 
            if(a.eq.0.)then 
              tmplam=-slope/(2.*b) 
            else 
              disc=b*b-3.*a*slope 
              if(disc.lt.0.) stop 'roundoff problem in lnsrch' 
              tmplam=(-b+sqrt(disc))/(3.*a) 
            endif 
            if(tmplam.gt..5*alam)tmplam=.5*alam 
          endif 
        endif 
        alam2=alam 
        fmin2=fmin
        fminold2=fminold 
        alam=max(tmplam,.1*alam) 
        goto 124

123     continue
cc        pgin=pg
***************************************************

        deltae=fp(99)/pe -1.

*        print*,' errf : ' , errf
*        print*,' errx : ' , errx
        
*        print*,'iter ',niter,'      Pe ',(fp(99)),
*     &         '      Pg  ',(pg)
*        print*,'iter ',niter,' deltaPe ',(deltae),
*     &   ' deltaP ',errx,'errf',errf
        if (errx.le.eps) then
          if (molions) then
            goto 6054
          else
            print*,'die_pe: converged without mol. ions.'
            print*,'iterating further with mol ions',tem,sngl(pe)
          endif
        endif
      enddo
6054  continue

ccc        do i=1,nmetal
ccc          nelemi = nelemx(i)
ccc          p(nelemi)=max(p(nelemi),smallest)
ccc        enddo

      if (NITER.ge.NIMAX.and.(errx.gt.eps.or.errf.gt.epsf))  then
* not converged 
        WRITE(6,6055) NIMAX
6055    FORMAT('*P(nelemi) DO NOT CONVERGE AFTER ',
     &         I4,' ITERATIONS')
        print*,'T, Pe ',Tem,sPe
        print*,'error on Press., on Pe:',errx,deltae,
     &             ' convergence criterium ',eps
        converge=.false.
        return
cc        stop
      endif
* everything converged
1051  continue
      if (.not.molions) then
* Spurious convergence without molecular ions.
        print*,'die_pe WARNING. Converged without molecular ions'
        print*,'continuing with unconverged equilibrium'
* Try to improve on this, BPz, and find a better solution than:
        converge=.false.
        return
      endif

* recompute partial pressures to be on the safe side.
      call funcv(fvec,fmin,pgindumm,pg)
      spg=pg
* store partial pressures in parptsuji array
      do i=1,nmetal
        nelemi=nelemx(i)
        parptsuji(i)=p(nelemi)
        parptsuji(i+nmetal)=p(nelemi)*kp1(nelemi)/pe
        parptsuji(i+2*nmetal)=parptsuji(i+nmetal)*kp2(nelemi)/pe
        parptsuji(i+3*nmetal)=parptsuji(i+2*nmetal)*kp3(nelemi)/pe
ccc        print*,'check die',nelemi,parptsuji(i),parptsuji(i+nmetal),
ccc     &   parptsuji(i+2*nmetal),parptsuji(i+3*nmetal)
      enddo
      do i=1,nmol
        parptsuji(i+4*nmetal)=ppmol(i)
ccc        print*,'check die',mol(i),ppmol(i)
      enddo
c      if (skiprelim) then
c      print*,'Die_pe_lu skipped prelim. Conv after ',niter,'iterations'
c      else
c      print*,'Die_pe_lu with prelim.    Conv after ',niter,'iterations'
c      endif
**      print*,niter,'iterations '
**,' deltaPe ',sngl(deltae),
cc     &                           ' deltaPat ',sngl(delta)

      errf=0.
      do i=1,nmetal
        errf=errf+abs(fvec(i))
      enddo
*      print*,'T Pe, pg ' ,tem,sngl(fp(99)),spg,' wanted Pe:',sngl(pe)
*      print*,'ERRF',sngl(errf),' ERRX',sngl(errx),' DPE',sngl(deltae),
*     &      sngl(tem),sngl(pe)
**      print*,'FVEC ',sngl(fvec)
*      do i=1,nmetal
*        if (kpe(nelemx(i))*p(nelemx(i)).ge.pe*0.01) print*,
*     &      'e-CONTRIB ',nelemx(i),kpe(nelemx(i))*p(nelemx(i))/pe
*      enddo
*      do j=1,nmol
*        if (elmolec(j)) then
*          if(abs(nelmolec(j)*ppmol(j)/pe).ge.0.01) then
*          print*,'e-CONTRIB ',mol(j),
*     &                -nelmolec(j)*ppmol(j)/pe
*          endif
*        endif
*      enddo

* make sure we save the electron pressure for next call.
       pep=spe

      RETURN
      END
