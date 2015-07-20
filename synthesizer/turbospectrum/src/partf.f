      subroutine partf(jatom,ion,temp,ntemp,u,ionpot)
*
* Correction for ThII, ThIII, UII and UIII introduced by BPz 20/09-2000
* yields partition functions with polynomial data from
* ref. Irwin, A.W., 1981, ApJ Suppl. 45, 621. Updated to read new
* input file on 21/08-1996 BPz
* ln u(temp)=sum(a(i)*(ln(temp))**i) 0<=a<=5
* Input:
* jatom = element number in periodic table
*       NOT implemented for molecules; see ref. in ref.
* ion   = 1 for neutral, 2 for once ionized and 3 for twice ionized
* temp  = array with ntemp values of the temperature
* ntemp = number of temperatures for which partf. is calculated
* Output:
* u     = partf. (linear scale) for iat,ion at the ntemp temperatures
*
      implicit none
*
      logical uexist(4,92),dummy
      integer ntemp,i,j,k,jatom,ion,iread,iaa,idumb
      real temp(ntemp),u(ntemp),spec,ipot,ip(0:4,92),tt,tt1,tt2
      real partitu2,partitu3,partitth2,partitth3
      doubleprecision a(0:5,1:4,1:92),ulog,t,aa(0:5),ionpot,ulog1,ulog2
      character*20 specname
      character*1 dum
*
      data iread /0/
      data ip/460*0./
      data uexist/368*.false./
*
      save iread,a,ip,uexist
*
      if(iread.ne.1) then
* read data if first call:
        
        open(67,file= 
     &    'DATA/IRWIN_atoms_v07.3.dat',
     &    status='old')
        dummy=.true.
        do while (dummy)
          read(67,'(a)') dum
          if (dum.ne.'#') then
            dummy=.false.
            backspace(67)
          endif
        enddo
        do while (.true.)
          read(67,*,end=99) idumb,spec,specname,ipot,iaa,aa
          j=int(spec)
          i=int((spec-j)*101.)+1
          if (i.lt.5)  then
* ip(i,j) contains ionisation pot. for ion. stage i+1. ip(0,j) contains 0,
* ip(1,j) contains IP for I to II etc.
            ip(i-1,j)=ipot
            uexist(i,j)=.true.
            do k=0,5
              a(k,i,j)=aa(k)
            enddo
          endif
        enddo
99      close(67)
        iread=1
        print *,'Atomic partition function data read'
      endif

* some ions not defined
      if (.not.uexist(ion,jatom)) then
        do i=1,ntemp
          u(i)=0.0
        enddo
      else
        do 30 i=1,ntemp
          tt=temp(i)
          if (temp(i).le.16000.) then
            if(temp(i).lt.1000.) then
              tt=1000.
              print*,'WARNING: atomic partf; temp<1000 K, using Q(1000)'
            endif
            t=dlog(dble(tt))
            ulog=   a(0,ion,jatom)+
     &           t*(a(1,ion,jatom)+
     &           t*(a(2,ion,jatom)+
     &           t*(a(3,ion,jatom)+
     &           t*(a(4,ion,jatom)+
     &           t*(a(5,ion,jatom))))))
            u(i)=real(dexp(ulog))
          else if(temp(i).gt.16000.) then
            print*,'WARNING: atomic partf; temp=',temp(i), 
     &             ' extrapolating ',
     &               'from Q(10000) and Q(16000)'
            tt1=10000.
            t=dlog(dble(tt1))
            ulog1=   a(0,ion,jatom)+
     &           t*(a(1,ion,jatom)+
     &           t*(a(2,ion,jatom)+
     &           t*(a(3,ion,jatom)+
     &           t*(a(4,ion,jatom)+
     &           t*(a(5,ion,jatom))))))
            tt2=16000.
            t=dlog(dble(tt2))
            ulog2=   a(0,ion,jatom)+
     &           t*(a(1,ion,jatom)+
     &           t*(a(2,ion,jatom)+
     &           t*(a(3,ion,jatom)+
     &           t*(a(4,ion,jatom)+
     &           t*(a(5,ion,jatom))))))
            ulog=(ulog2-ulog1)*(tt-tt2)/(tt2-tt1) + ulog2
            u(i)=real(dexp(ulog))
          endif
   30   continue
        if (jatom.eq.92.and.ion.eq.2) then
C correction for UII    BPz 20/09-2000
          do i=1,ntemp
            u(i)=partitu2(temp(i))
          enddo
        else if (jatom.eq.92.and.ion.eq.3) then
C correction for UIII    BPz 20/09-2000
          do i=1,ntemp
            u(i)=partitu3(temp(i))
          enddo
        else if (jatom.eq.90.and.ion.eq.2) then
C correction for ThII    BPz 20/09-2000 from Holweger
          do i=1,ntemp
            u(i)=partitth2(temp(i))
          enddo
        else if (jatom.eq.90.and.ion.eq.3) then
C correction for ThIII    BPz 20/09-2000 from Holweger
          do i=1,ntemp
            u(i)=partitth3(temp(i))
          enddo
        endif
      endif
      ionpot=ip(ion,jatom)

      return
*
      end
**************************************************************************8

	FUNCTION PARTITU2(T)

C 	calcule la fonction de partition de UII
C	use 'tableu2.dat'	
	IMPLICIT NONE
	INTEGER*4 k, FLAG,nlevel
	REAL*4 C2, T, U, E(100), G(100),partitu2
	DATA  FLAG /0/ 
	SAVE
	IF(FLAG.EQ.0)THEN
           OPEN(70, file=
CC'/b1/plez/BIGGRID/DATA/tableU2.dat',
     &    'DATA/tableU2.dat',
     &            status='old')
           READ(70, '(A)')
           nlevel=0
           DO k=1,100 	
             READ(70,*,end=99) E(k), g(k)
             nlevel=k
           ENDDO
99         continue
	   CLOSE (70)
           FLAG=1
	ENDIF
	C2=1.43877
	U=0.
	DO K=1,nlevel
		U=U+G(K)*EXP(-C2*E(K)/T)
	ENDDO
	partitu2=U
C
	END
**************************************************************************8

	FUNCTION PARTITU3(T)

C 	calcule la fonction de partition de UIII
C	use 'tableU3.dat'	
C from Cayrel 20/09-2000
	IMPLICIT NONE
	INTEGER*4 k, FLAG, nlevel
	REAL*4 C2, T, U, E(100), G(100),partitu3
	DATA  FLAG /0/ 
	SAVE
	IF(FLAG.EQ.0)THEN
           OPEN(70, file=
CC'/b1/plez/BIGGRID/DATA/tableU3.dat',
     &    'DATA/tableU3.dat',
     &            status='old')
           READ(70, '(A)')
           nlevel=0
           Do k=1,100
             READ(70,*,end=99) E(k), g(k)
             nlevel=k
           ENDDO
99         continue
           CLOSE (70)
           FLAG=1
	ENDIF
	C2=1.43877
	U=0.
	DO K=1,nlevel
		U=U+G(K)*EXP(-C2*E(K)/T)
	ENDDO
	partitu3=U

C
	END
**************************************************************************8
	
	function partitth2(temp)
*
* Inspired from : SUBROUTINE IONDIS(TEMP,PELE,MODE,IPRINT,PGAS,DENSNC)
* Which was received From: supas059@astrophysik.uni-kiel.de (Hartmut Holweger)
* Date: Fri, 15 Sep 2000 11:28:13 +0200
* To: plez@graal.univ-montp2.fr
* Subject: iondis.f

      tn=-5040./temp*log(10.)

C     THORIUM II
C     ZALUBAS(1968),ZALUBAS&CORLISS(1974),KLINKENBERG(1950)
      partitth2=4.+10.0*EXP(0.21*TN)+20.0*EXP(0.53*TN)+24.0*EXP(0.78*TN)
     1 +22.0*EXP(0.89*TN)+22.0*EXP(1.04*TN)+52.0*EXP(1.17*TN)
     2 +60.0*EXP(1.33*TN)+102.0*EXP(1.64*TN)+104.0*EXP(1.90*TN)
     3 +88.0*EXP(2.19*TN)+58.0*EXP(2.40*TN)+128.0*EXP(2.62*TN)
     4 +138.0*EXP(2.94*TN)+140.0*EXP(3.25*TN)+128.0*EXP(3.56*TN)
     5 +344.0*EXP(4.03*TN)+596.0*EXP(4.65*TN)+46.0*EXP(5.27*TN)
     6 +878.0*EXP(6.4*TN)

      end
**************************************************************************8

	function partitth3(temp)
*
* Inspired from : SUBROUTINE IONDIS(TEMP,PELE,MODE,IPRINT,PGAS,DENSNC)
* Which was received From: supas059@astrophysik.uni-kiel.de (Hartmut Holweger)
* Date: Fri, 15 Sep 2000 11:28:13 +0200
* To: plez@graal.univ-montp2.fr
* Subject: iondis.f

      tn=-5040./temp*log(10.)

C     THORIUM III
C     ZALUBAS(1968),ZALUBAS&CORLISS(1974),KLINKENBERG(1950)
      partitth3=14.+5.0*EXP(0.06*TN)+21.0*EXP(0.37*TN)+41.0*EXP(0.58*TN)
     1 +35.0*EXP(0.84*TN)+44.0*EXP(1.06*TN)+42.0*EXP(1.35*TN)
     2 +26.0*EXP(1.85*TN)+60.0*EXP(2.47*TN)+31.0*EXP(3.57*TN)
     3 +41.0*EXP(4.51*TN)+84.0*EXP(5.66*TN)+37.0*EXP(6.79*TN)

      end
