      subroutine hydropac(lunit,xlboff)
*
* Adapted from MARCS35 routine. BPz 30/06-1998. Redone 2003-11-14/KE
*  This version is from Paul Barkleem.

* further modified by BPz 10/09-2004 for compatibility with Montpellier
* version of bsyn. Only slight remaining problem: the continuous
* contribution of H I dissolved states is properly included in the
* continuum opacity for computation of spectrum, but is not included
* in absocont, which is used for dump of opacity for multi for example.
* maybe it is now. 071204
*
* BPz 27/02-2015: 

* compute H I line opacity     BPz 3/06-1998
*
      implicit none

      include 'spectrum.inc'

      character comment*100,species*20,lname*9
      integer ion,nline,ntau,iter,nlo,nup,l,k,i,lstart,lstartp,
     &        nmy,nlbldu,iint,iweak,lunit,maxlam,ll,
     ,        jlcont, nlcont
      doubleprecision ionpot,hlambda,wave,xl1,xl2,del,xlboff,xlambda
      real*8  vacair,edge(6)
      real fpart(ndp),pe,t,xlo,pg,xi,mum,ro,eps,xmyc,
     &     xup,gf,ne(ndp),
     &     nh1(ndp),nhe1(ndp),hnorm,hnormnostim(ndp),stim,
     &     diff,diffp,theta(ndp),hckt(ndp),abso,absos,absocont,
     &     dopple(ndp),contrib,xkapr,cross,hlinop,h1bfg,alpha,
     ,     ee(6),d0(ndp),d(lpoint,ndp)
      logical  lymanalpha, usedam,notfound
*      
      common /atmos/ t(ndp),pe(ndp),pg(ndp),xi(ndp),mum(ndp),ro(ndp),
     &               ntau
      doubleprecision presneutral,presion,presion2,presion3, h1bfgc(30)    
      common /orderedpress/ presneutral(ndp,100),presion(ndp,100),
     &                      presion2(ndp,100),presion3(ndp,100)
      common /pieces/ xl1,xl2,del,eps,nmy,nlbldu,iint,xmyc,iweak
      common /large/ xlambda(lpoint),maxlam,abso(ndp,lpoint),
     &               absos(ndp,lpoint),absocont(ndp,lpoint)
      common/rossc/ xkapr(ndp),cross(ndp)
      common/continuum/nlcont,jlcont(lpoint)
      
      character*50 zMC
      integer znlc ,  ii
      real    zxlm,zBP(ndp),zXC(ndp),zS(ndp),zXI(ndp)
      
           
      data  ee/ 0.,  10.199, 12.088, 12.749, 13.055, 13.221 /
c In babsma I use slighty different edges from continuous opacity file (air wavelength)
c     data  edge/ 912.00, 3647.00, 8207.00, 14588.00, 22794.00, 32915.00/
c This is not anymore the case, since we use the same jonabs_vac.dat file
c in version 7.3, and transform them to vacuum, just as in marcs. 
c   BPz 21/01-2008
c vacuum edge used in MARCS :     
cc      data edge/ 912.00, 3647.98, 8209.26, 14591.99, 22800.22, 32923.98/
      data edge/ 912.00, 3648.04, 8209.26, 14591.99, 22800.22, 32923.98/  
      data  h1bfgc/1.0711223,  -3.0033216e-04,  
     ;             1.0648009,  -5.1846584e-05,
     ;             1.0478152,  -1.7104666e-05,
     ;             1.0443061,  -8.2075951e-06,
     ;             1.0421281,  -4.6750806e-06,
     ;             1.0402398,  -2.9468910e-06,
     ;             1.0347312,  -1.9145532e-06,
     ;             1.0335221,  -1.3661146e-06,
     ;             1.0335879,  -1.0323393e-06,
     ;             1.0282216,  -7.4292798e-07,
     ;             1.0337342,  -6.4224310e-07,
     ;             1.0325113,  -5.1576699e-07,
     ;             1.0321166,  -4.2484186e-07,
     ;             1.0324288,  -3.6090637e-07,
     ;             1.0329178,  -3.0957602e-07/  
*
c bsyn uses air wavelengths:
      do k=2,6
        edge(k)=vacair(edge(k))
*	print*,'k,edge(k)',k,edge(k)
      enddo

      do 2 k=2,6

* BPz 061204: potential problem: catches only k=2, if in 
* lambda interval. Other edges ignored, and no
* provision if more than one edge.

ccccc        if(xl1.gt.edge(k) .or. xl2.lt.edge(k)) goto 3
        if(xl1.lt.edge(k) .and. xl2.gt.edge(k)) then
          ll=0
          diffp=1.e30
          do l=1,maxlam
            diff=abs(xlambda(l)-edge(k))
            if(diff.ge.diffp) then
              ll=l-1
              goto 1
            endif
            diffp=diff
          enddo
    1     continue
c
c          print*,'Hlinop, changed edge ',k,edge(k),xlambda(ll)
c          print*,'nearby lambdas',xlambda(ll-1),xlambda(ll+1)
          edge(k)=xlambda(ll)		
        endif
    2 continue
    3 continue
c      
      rewind (lunit)
      read(lunit,*) species,ion,nline
      read(lunit,*) comment
*      print*,species,ion,nline,comment
      if (species(1:6).ne.'01.000'.or.ion.ne.1) then
        print*, 'wrong H line data file!'
        print*,species,ion
        stop 'ERROR!'
      endif
*
***change nline   temporarily!
*      nline=37
*      
***bort:!
* efter 1 continue:
*        do l=2,nlcont-1
*	  if(jlcont(l).eq.ll) goto 111
*	enddo
*        print*,"Hydropac: check:ll,l,jlcont(l,l-1),l",ll,l,jlcont(l),
*     ,         jlcont(l-1),xlambda(jlcont(l-1))
*	stop 'Hydropac: ll-value not found?!'
*  111   continue
*        print*,"Hydropac: check:ll,l,jlcont(l,l-1),l",ll,l,jlcont(l),
*     ,         jlcont(l-1),xlambda(jlcont(l-1))	
*******************
*
      call partf(1,1,t,ntau,fpart,ionpot)

      do k=1,ntau
*        print*,"k,t(k)",k,t(k)
        theta(k)=5040./t(k)
        ne(k)=pe(k)/(t(k)*1.38066e-16)
        nh1(k)=sngl(presneutral(k,1))/(t(k)*1.38066e-16)
        nhe1(k)=sngl(presneutral(k,2))/(t(k)*1.38066e-16)
* hckt=hc/kT ready for lambda in AA.
c        hckt(k)=1.43877e8/t(k)
        hckt(k)=2.99792458e10*6.626075e-27/1.380658e-16*1.e8/t(k)
        dopple(k)=sqrt( xi(k)**2 * 1. + 
     &                 2.*1.380658e-16*t(k)/1.6738e-24) /
     &                 2.99792458e10
      enddo

c
cccc
c  First add the continuous part of the H opacity from n=2,3,.,6 to dissolved
c  states following Daeppen et al. 1987. See subr. dcalc
c
cccc
c Note added 20030902: Mainly by luck the loop below starts at i=2. In PB's
c hbop package the corresponding part goes from the ground state (& up to i=15);
c this caused problems since the 'dissolved n=1 state' contributes 'too much'
c opacity at 'long' wavelengths (> 2000? AA). Whether this is due to a too 
c large extrapolated sigma-bf for n=1 or to a 'd' that does not decrease fast
c enough (with wavelength) is presently not known. /KE
c
      do i=2,6
c        print*,' i edge(i)',i,edge(i)
        do k=1,ntau
*	  print*,'k',k
*          if(k.eq.30)
*     ,	  print*,k,xlambda(1),maxlam,nh1(k),ne(k),nhe1(k),d(1,k),t(k)
          call dcalc(i,xlambda,maxlam,nh1(k),ne(k),nhe1(k),d(1,k),t(k))     
          do l=1,maxlam
*	    print*,' l,xlam ',l,xlambda(l)
	    if(xlambda(l) .ge. edge(i)) then
	      stim = 1. - exp(-hckt(k)/xlambda(l))
              h1bfg = h1bfgc(2*i-1) + h1bfgc(2*i)*xlambda(l)
 	      alpha = 1.044e-26 * h1bfg * xlambda(l)**3 / float(i**5)
CC added float() above and below on 09/05-2007.    BPz.
              contrib = nh1(k)* 2.*float(i)**2/fpart(k)*
     &               10.**(-ee(i)*theta(k))
     ,               / ro(k) * d(l,k) * alpha * stim  / xkapr(k) 
c above division by xkapr added 030829...... 
c    
              abso(k,l)=abso(k,l) + contrib
              absocont(k,l)=absocont(k,l) + contrib
*	      
cc              if (contrib.le.0.) print*,'H',xlambda(l),contrib
*      if(k.eq.30) print*," DAMc",i,k,l,xlambda(l),d(l,k),contrib,
*     ,   contrib/abso(k,l)
***
* now check if this lambda-depth point is in the 'continuum'array, if so
* then add to the continuous x also this contribution
*
              notfound=.true.
              ii=1
              do while (notfound.and.ii.le.nlcont)
                if (jlcont(ii).eq.l) then
                  read(14,rec=ii) zMC,znlc,zxlm,zBP,zXC,zS,zXI
                  zXC(k)=zXC(k)+contrib
                  write(14,rec=ii) zMC,znlc,zxlm,zBP,zXC,zS,zXI
cc	          print*,"updated XC,depth",k,zxlm,xlambda(l),contrib
                  notfound=.false.
                endif
                ii=ii+1
              enddo
* BPz 27/02-2015: abso changed to absocont in following statement
              if(contrib/absocont(k,l) .le. eps) goto 5
	    endif
          enddo      
   5      continue
	enddo      
      enddo
*
cc      do l=1,maxlam
cc        write(88,*) xlambda(l),abso(30,l),absos(30,l),absocont(30,l)
cc      enddo      
cc      write(88,*) "*****"
*
      lstartp=1
      do 20 i=1,nline
        read(lunit,*) hlambda,nlo,nup,xlo,xup,gf,lname
cc	if(nlo.gt.1)  hlambda=vacair(hlambda)
*        print*,'Hydropac; line ',i,' of ',nline,' ',nlo,nup,hlambda
**********************************
* removed by BPz 24/01-2012
* this skipped Lyman series
*	if(nlo.eq.1) goto 20
**********************************
ccc
c Determine the fraction (d) of H b-b opacity that is lost by absorption to
c dissolved states. Left of the line is (1-d). 
c See subr. dcalc  and Daeppen, Anderson & Mihalas 1987
c
     	usedam =.false.
     	if(nlo .ge. 2 .and. nlo .le. 6) usedam=.true.
     	if(usedam) then
c    	    calculate d at line-centre. (1-d)*bb is added further down. 
     	  do k=1,ntau
     	    call dcalc(nlo,hlambda,1,nh1(k),ne(k),nhe1(k),d0(k),
     ,                 t(k))	
     	  enddo   
     	endif
        lymanalpha=.false.
	if(nlo.eq.1 .and. nup.eq.2) lymanalpha=.true.

        gf=10.**gf
*
* in bsyn we normally use airwavlengths
*	
*        if (hlambda.gt.2000.) hlambda=airvac(hlambda)
* we search for the OS point nearer the line center.
        if (hlambda.le.xl1) then
          lstart=1
        else if (hlambda.gt.xl2) then
          lstart=maxlam+1
        else
          diffp=1.e30
          do l=lstartp,maxlam
            diff=abs(hlambda-xlambda(l))
            if (diff.gt.diffp) then
              lstart=l-1
* make sure everything is fine (this test has a large margin)
              if (diff.gt.(xlambda(min(maxlam,lstart+1))-
     &                     xlambda(max(1,lstart-1)))) then
                print*,'Problem in hydropac!'
                print*,' hydropac: hlambda',sngl(hlambda),
     &            'not nearest OS:',xlambda(lstart),lstart
                print*,'xlb-1,xlb,xlb+1:',xlambda(max(1,lstart-1)),
     &                  xlambda(lstart),xlambda(min(maxlam,lstart+1))
                stop 'ERROR!'
              endif
              goto 10
            endif
            diffp=diff
          enddo
        endif
10      continue
        lstartp=lstart 
	
*        print*,'H opac: H lambda, xlb start:',sngl(hlambda),lstart,
*     .    lstartp,xlambda(lstart)
*	if(usedam) print*," DAM ",i,sngl(hlambda),d0(30)  
        do k=1,ntau
* numerical factor 0.014... is sqrt(pi)e2/mec
c          hnormnostim(k)=0.014974*gf*nh1(k)/fpart(k)*
c     &                   10.**(-xlo*theta(k))/rho(k)

          hnormnostim(k)=0.026540045*gf*nh1(k)/fpart(k)*
     .       10.**(-xlo*theta(k))/ro(k)
     
*          if(nlo.eq.2     .and. k.eq.5) then
*            print*,"Hydropac:t,ne,nh1,nhe,dopp",t(k),ne(k),nh1(k),
*     .        nhe1(k),dopple(k),hnormnostim(k)
*          endif
          do l=lstart,maxlam
	    ll=l
cctest: 	do l=nlb-100,nlb
            wave=xlambda(l)
            stim=1.-exp(-hckt(k)/xlambda(l))
            hnorm=hnormnostim(k)*stim
cc	    
            contrib=HLINOP(wave,nlo,nup,hlambda,t(k),ne(k),
     *                nh1(k),nhe1(k),dopple(k))  /xkapr(k)
*          If( k.eq.31 .and. abs(wave-hlambda).lt.10. 
*     *        .and. nlo.eq.2 .and. nup.eq.3 )
*     *       print*,wave,contrib,hnorm
     
c     *     print 2224,contrib,hnorm,stim

            if(usedam) then
	      contrib=contrib*(1. - d0(k))
*	      if(k.eq.30 .and. l.eq.lstart)
*     , print*," DAMl",contrib,contrib*hnorm/abso(k,l),(1. - d0(k))
	    endif
c    
            contrib=contrib*hnorm
	    
cc              if (contrib.le.0.) print*,'H',xlambda(l),contrib

            abso(k,l)=abso(k,l) + contrib
c	    if()then
c	      write(59,   )
c	    endif

*       If( k.eq.30 )
*     . Print2222,l,wave,contrib,abso(k,l), contrib/abso(k,l)
 2222 format("l,wave,contrb,osopx,ratio",i6,f9.1,1p,3e12.3)
 2223 format("l,wave,c*ntrb,osopx,ratio",i6,f9.1,1p,3e12.3)
 2224 format("  hlinop,hnorm,stim      ",15x,1p,3e12.3)
 
* BPz 27/02-2015: abso changed to absocont in following statement
            if (contrib/absocont(k,l).le.eps) goto 11
*
* we jump out if contribution less than 10-4 of the continuum at that depth.
* Before calling this routine we have only added continuum opacities in osop
*
          enddo
11        continue

cc          print*,k,sngl(hlambda),'computed out to: ',xlb(ll),ll
	  ll=0
cc	    Print*,"In next loop we start at l= ",lstart-1

          do l=lstart-1,1,-1
            if(l.le.0) Print*," In the loop...;l= ",l
	    ll=l
            wave=xlambda(l)
            stim=1.-exp(-hckt(k)/xlambda(l))
            hnorm=hnormnostim(k)*stim
            contrib=HLINOP(wave,nlo,nup,hlambda,t(k),ne(k),
     *                nh1(k),nhe1(k),dopple(k))  /xkapr(k)
            if(usedam) contrib=contrib*(1. - d0(k))
 	    contrib=contrib*hnorm     

cc              if (contrib.le.0.) print*,'H',xlambda(l),contrib

            abso(k,l)=abso(k,l) + contrib
	    
*      If( k.eq.30)
*     . Print2223,l,wave,contrib,abso(k,l), contrib/abso(k,l)	
    
* BPz 27/02-2015: abso changed to absocont in following statement
            if (contrib/absocont(k,l).le.eps) goto 12
          enddo
12        continue

c          If(ll.gt.0) then
c	    print*,k,sngl(hlambda),'computed*out to: ',xlb(ll),ll
c	  else
c	    print*,k,sngl(hlambda),': no lambda-points considered'
c	  endif

        enddo
c	print*," hydropac; after big depthloop.. "
*
cc      do l=1,maxlam
cc        write(88,*) xlambda(l),abso(30,l),absos(30,l),absocont(30,l)
cc      enddo      
cc      write(88,*) "*****"
*
  20    continue

      close(lunit)

      return
      end

      subroutine dcalc(ni, wave, nwav, nH, ne, nHe, d, temp)
c      
c  dcalc calculates the D-function of Daeppen, Anderson & Mihalas 1987
c  (ApJ 319, 195). Accounts for opacity from bound to dissolved states in H.
c  Coded by Paul Barklem. Translated to Fortran 20030204/KE
c       
      implicit none
      real     h, c, ionH, nH, ne, nHe, d(*),
     ,         x, xx, xi, temp
cc this was wrong!  BPz 16/05-2007     ,         x, xx, xi, wi, wstar, wcalc, temp
      real*8   wave(*),wi,wstar,wcalc
      integer  ni, ij, nwav
      data h/6.626e-27/, c/2.997925e10/, ionH/2.17991e-11/ 
c      
      do ij=1,nwav
	d(ij)=1.
        x=1./(float(ni)**2) - h*c/sngl((wave(ij))*1.e-8*ionH)
	if(x .gt. 0.) then
	  xx=1./sqrt(x)
	  xi=float(ni)
	  wi=wcalc(nH, ne, nHe, xi, temp)
	  wstar=wcalc(nH, ne, nHe, xx, temp)
	  d(ij)=1.-(wstar/wi)
	endif	
      enddo
      return	
      end
      
      REAL*8 FUNCTION WCALC(NH, NE, NHE, NS, TEMP)	
C  	
C  Computes the occupation probability for a H atom in 
C  state with effective principle quantum number NS in a plasma
C  enviroment with NH, NE, NHE number densities of H, ions,
C  and He atoms respectively.  This code assumes the perturbing
C  neutral H and He are in the ground state, (noting the hard
C  sphere approximation used is quite crude anyway) and that ions 
C  are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
C
C  See eqn 4.71 Hummer & Milhalas (1988) 
C  Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
C  which is a fit to detailed calculation including correlations,
C  thus the temperature dependence
C
C  Sizes of neutral atoms adopted are sqrt(<r^2>)
C  
C  Coded by Paul Barklem and Kjell Eriksson, Aug 2003
C 	
      IMPLICIT NONE
      REAL NH, NE, NHE, NS, TEMP
      REAL*8 IONH, A0, E, PI, K, CHI, RIH, NEUTR, ION 
      REAL*8 F, A, BETAC, X, WNEUTR, WION      
      DATA IONH/2.17991E-11/ 
      DATA A0  /5.29177E-9/
      DATA E   /4.803207E-10/ 
      DATA PI  /3.1415927/
C  
C  Neutral perturbers
C
      CHI=IONH/(NS*NS)
      RIH=SQRT(2.5*NS**4 + 0.5*NS*NS)*A0
      NEUTR=NH*(RIH + 1.73*A0)**3 + NHE*(RIH + 1.02*A0)**3
      WNEUTR = EXP(-4.*PI/3. * (NEUTR))
C
C  Charged perturbers
C      
      K=1.
      IF(NS .gt. 3.) THEN
        K=16./3.*(NS/(NS+1.))**2 * (NS + 7./6.)/(NS*NS+NS+0.5)
      ENDIF
c      ION  = NE*16.*(E**2/CHI/DSQRT(K))**3
c      WION = DEXP(-4.*PI/3. * (ION))
      IF ((NE.GT.10.).AND.(TEMP.GT.10.)) THEN  ! just in case!!! 
        A = 0.09 * NE**0.16667 / SQRT(TEMP)
        X = (1. + A)**3.15
        BETAC = 8.3E14 * NE**(-0.667) * K / (NS**4.d0)
       F = 0.1402d0*X*BETAC*BETAC*BETAC/(1.+0.1285*X*BETAC*DSQRT(BETAC))
        WION = F/(1.d0+F)
      ELSE 
        WION = 1.0d0     
      ENDIF 
C      
      WCALC = WION * WNEUTR 
      RETURN
      END
