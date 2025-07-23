C
      SUBROUTINE JON(T,PE,IEPRO,PG,RO,E,IOUTR,kk)
C
C
C        THIS ROUTINE COMPUTES IONIZATION EQUILIBRIA FOR A GIVEN TEMPERATURE 
C        (T, EXPRESSED IN KELVIN) AND A GIVEN ELECTRON PRESSURE (PE, IN
C        DYNES PER CM2). THE FRACTIONS OF IONIZATION ARE PUT IN THE ANJON VECTOR
C        AND THE PARTITION FUNCTIONS ARE PUT IN PART. 
C*************************************
C        BPz (20/02-2012) IEPRO is not used for this anymore. The role of IEPRO
C        to tell whether one stops on the absence of convergence in eqmol_pe_lu (iepro>0)
C        or if one goes back with another try after modifying the input guess. (iepro<0)
C*************************************
C C      IF IEPRO IS GREATER THAN
C C       ZERO, THE GAS PRESSURE (PG,IN DYNES PER CM2), DENSITY (RO, IN GRAMS
C C       PER CM3) AND INNER ENERGY (E, IN ERGS PER GRAM) ARE ALSO EVALUATED.
C        N O T E . RADIATION PRESSURE IS NOT INCLUDED IN E.
C
C        THE ENERGIES OF IONIZATION ARE REDUCED BY DXI, FOLLOWING BASCHEK ET 
C        AL., ABH. HAMB. VIII, 26 EQ. (10). THESE REDUCTIONS ARE ALSO MADE IN
C        THE COMPUTATION OF E.
C        THE ENERGY OF DISSOCIATION FOR H- HAS BEEN REDUCED BY 2*DXI, FOLLOWING
C        TARAFDAR AND VARDYA, THIRD HARV. SMITHS. CONF., PAGE 143. THE FORMATION
C        OF MOLECULES IS CONSIDERED FOR T LESS THAN TMOLIM.
C
C        IF IOUTR IS GREATER THAN ZERO, A DETAILED PRINT-OUT WILL BE GIVEN.
C
C
C        THE FUNCTION  QTRAV AND SUBROUTINE MOLEQ ARE CALLED.
C        THEY CALL QAS AND MOLFYS RESPECTIVELY.
C
C        DIMENSIONS NECESSARY
C        A(5),DQ(4),F(MAX(NJ)),PFAK(MAX(NJ)),RFAK(JMAX)
C        DIMENSIONS OF ARRAYS IN COMMONS /CI1/,/CI4/,/CI5/ AND /CI6/ ARE
C        COMMENTED ON IN SUBROUTINE INJON.
C        JMAX IS THE TOTAL NUMBER OF STAGES OF IONIZATION, INCLUDING NEUTRAL
C             ATOMS.
C        NJ(I) IS THE NUMBER OF STAGES OF IONIZATION, INCLUDING THE NEUTRAL
C             STAGE, FOR ELEMENT I.
C
C
      include 'spectrum.inc'
      integer maxim
      parameter (maxim=1000)
C
      DIMENSION DQ(4),F(5),PFAK(5),RFAK(45)
      DIMENSION PRESMO(30)
      COMMON/CI1/FL2(5),PARCO(45),PARQ(180),SHXIJ(5),TPARF(4),
     *XIONG(16,5),EEV,ENAMN,SUMH,XKBOL,NJ(16),IEL(16),NEL,SUMM
      COMMON/CI4/ IELEM(16),ION(16,5),TMOLIM,MOLH
      COMMON/CI5/ABUND(16),ANJON(16,5),H(5),PART(16,5),DXI,F1,F2,F3,F4,
     *F5,XKHM,XMH,XMY
      COMMON/CI6/IQFIX(16,5),NQTEMP,TP
      COMMON/CI7/A(5),PFISH,ITP
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/RABELL/XXRHO(NDP),XYRHO
      COMMON/CI8/YYPG,YYRHO,YYE
      COMMON/CMOL1/EH,FE,FH,FHE,FC,FCE,FN,FNE,FO,FOE,FK,FKE,FS,FSE
      COMMON/CMOL2/NMOL,PK(30)
      COMMON/CARC3/F1P,F3P,F4P,F5P,HNIC,PRESMO
      character*20 names(maxim)
      integer tselem,niter
      common/plez/names,tselem(100),fictpres_h
      real    phydro
      common /cphydro/ phydro
* phydro is the total hydrogen pressure
* common for partial pressures
      logical tsuswitch,tsuji,skiprelim
      doubleprecision parptsuji
      common /tsuji/ tsuji,tsuswitch,nattsuji,nmotsuji,
     &               parptsuji(maxim+400)
      real rho,ejontsuji
      common/rhotsu/ rho,xmytsuji,ejontsuji
      real xmass( maxim+400),atmass(100)
      doubleprecision ndensity,molweight
      common /density/ndensity,molweight,xmass,atmass
      real pgpgpg
      logical chck
      data chck/.false./ ,pgpgpg/-1.0/
      data kkp/-1/
      data pep/-1.e10/
C
C
C STATEMENT FUNCTION FOR 10.**
      EXP10(X)=EXP(2.302585*X)
C
*****************************************************************************************
* Babsma, jon. PG is guess PG in input, replaced by computed Pg at return
*
cc      print*,'entering jon with iepro = ',iepro
cc      print*,'jon, T, Pe, pgin = ',t,pe,pg
      pgpgpg=pg
*****************************************************************************************
      ITP=1
C Is it the same T,pe point?
      if (kk.eq.kkp) then
       if ((abs((t-tp)/t).lt.1.e-8).and.(abs((pe-pep)/pe).lt.1.e-8)) 
     &   then
cc        print*,'CALLING jon with same T, P. ',
cc     &      'Exiting jon without further calculations'
cc        print*,' JON T=',T,' Pe=',pe,' kk=',kk,' Previous:',tp,pep,kkp
* but first make sure we return correct values. The compiler, even
* with -K, does not save variables in call statements.
        ro=rosave
        pg=pgsave
        e=esave
        return
       endif
      endif
C
C IS T=THE TEMPERATURE OF THE PRECEDING CALL ?
      IF(ABS((T-TP)/T).LT.1.E-8)GO TO 53
   51 ITP=0
C
C SOME QUANTITIES, ONLY DEPENDENT ON T
      TETA=5040./T
      TETA25=1.202E9/(TETA*TETA*SQRT(TETA))
      DO52 J=1,5
   52 A(J)=FL2(J)*TETA
C A=ALFA(BASCHEK ET AL., CITED ABOVE)
C
      IF(NQTEMP.EQ.0)GO TO 53
C
C PREPARATION FOR INTERPOLATION OF PARTITION FUNCTIONS IN T
      DQ(1)=1.
      DQ(2)=T-TPARF(1)
      DQ(3)=DQ(2)*(T-TPARF(2))
      DQ(4)=      DQ(3)*(T-TPARF(3))
C
C        SOME QUANTITIES ALSO DEPENDENT ON PE
C        THE PFAK FACTORS ARE USED IN THE SAHA EQUATION. H(J) IS THE
C        QUANTUM NUMBER OF THE CUT OF THE PARTITION FUNCTIONS (ACCORDING
C        TO BASCHEK ET AL., CITED ABOVE) FOR J-1 TIMES IONIZED ATOMS. H IS
C        USED IN QAS.
C
C        XNEL= THE ELECTRON (NUMBER) DENSITY (PER CM3)
C        PFISH= P(FISCHEL AND SPARKS, ASTROPHYS. J. 164, 359 (1971)) IS USED IN
C        FUNCTION QAS.
C
   53 DXI=4.98E-4*TETA*SQRT(PE)
      DUM=TETA25/PE
      DIM=EXP10(DXI*TETA)
      PFAK(1)=DIM*DUM
      SQDXI=1./SQRT(DXI)
      H(1)=SHXIJ(1)*SQDXI
      DO54 J=2,5
      PFAK(J)=PFAK(J-1)*DIM
   54 H(J)=SHXIJ(J)*SQDXI
      XNEL=PE/(XKBOL*T)
      PFISH=4.2E3/XNEL**0.166666667
C
C        PARTITION FUNCTIONS AND IONIZATION EQUILIBRIA
C
      XNECNO=0.
      XNENH=0.
      EJON=0.
      JA=1
C
C        BEGINNING OF LOOP OVER ELEMENTS ('THE I-LOOP').
      DO24 I=1,NEL
      NJP=NJ(I)
C
C        SHOULD ELEMENT NO. I BE CONSIDERED
      IF(IELEM(I).GT.0)GO TO 9
      DO 55 J=1,NJP
      ANJON(I,J)=0.
      PART(I,J)=0.
   55 CONTINUE
      GO TO 23
C
C        BEGINNING OF LOOP OVER STAGES OF IONIZATION ('THE J-LOOP')
    9 DO19 J=1,NJP
      JM1=J-1
C
C        SHOULD STAGE OF IONIZATION NO. J BE CONSIDERED
      IF(ION(I,J).GT.0)GO TO 10
      ANJON(I,J)=0.
      PART(I,J)=0.
      GO TO 18
C
C        WHICH KIND OF PARTITION FUNCTION SHOULD BE COMPUTED
C
   10 IF(IQFIX(I,J)-1)14,11,13
   11 IF(T.LT.TPARF(1).OR.T.GT.TPARF(4))GO TO 13
      PARTP=PART(I,J)
      IF(ITP.GT.0)GO TO 15
C
C        PARTITION FUNCTIONS TO BE INTERPOLATED IN T
      JPARF=(JA-1)*4+1
      PARTP=0.
      DO12 IP=1,4
      PARTP=PARTP+PARQ(JPARF)*DQ(IP)
   12 JPARF=JPARF+1
      GO TO 15
C
C        PARTITION FUNCTIONS FOLLOWING TRAVING ET AL., ABH. HAMB. VIII,1 (1966)
   13 PARTP=QTRAV(TETA,H(J),J,JA)
      GO TO 15
C
C        THE PARTITION FUNCTION IS CONSTANT
   14 PARTP=PARCO(JA)
   15 PART(I,J)=PARTP
C
C        IONIZATION EQUILIBRIA AND TOTAL NUMBER OF ELECTRONS
C
      IF(J.LE.1)GO TO 19
      IF(ITP.GT.0)GO TO 17
      RFAK(JA)=EXP10(-XIONG(I,JM1)*TETA)
   17 F(JM1)=PFAK(JM1)*RFAK(JA)*PARTP/PART(I,J-1)
      GO TO 19
   18 IF(J.GT.1)F(JM1)=0.
   19 JA=JA+1
C        END OF 'THE J-LOOP'
C
      FIL=1.
      DO20 J=2,NJP
      LL=NJP-J+1
   20 FIL=1.+F(LL)*FIL
      ANJON(I,1)=1./FIL
      XNEN=0.
      DO21 J=2,NJP
      JM1=J-1
      ANJON(I,J)=ANJON(I,JM1)*F(JM1)
      IF(I.LE.1)GO TO 24
      FLJM1=JM1
   21 XNEN=ANJON(I,J)*FLJM1+XNEN
      IF(I.GT.2.AND.I.LT.6) XNECNO=XNECNO+XNEN*ABUND(I)
      XNENH=XNEN*ABUND(I)+XNENH
C        XNENH=NUMBER OF ELECTRONS FROM ELEMENTS OTHER THAN HYDROGEN (Q IN
C        MIHALAS, METH. COMP. PHYS. 7, 1 (1967), EQ. (35))
C        XNECNO=NUMBER OF ELECTRONS FROM ELEMENTS OTHER THAN H, C, N, O
C
C
C        COMPUTATION OF THE ENERGY OF IONIZATION (EJON). HYDROGEN IS NOT
C        INCLUDED.
C
      XERG=0.
C        XERG= THE ENERGY OF IONIZATION PER ATOM (IN ELECTRON VOLTS)
C
      DO22 J=2,NJP
      JM1=J-1
      FLJM1=JM1
   22 XERG=ANJON(I,J)*(XIONG(I,JM1)-DXI*FLJM1)+XERG
      EJON=XERG*ABUND(I)+EJON
      GO TO 24
   23 JA=JA+NJP
   24 CONTINUE
C        END OF 'THE I-LOOP'
C
C
      XNECNO=XNENH-XNECNO
ccc      TP=T
ccccccc      IF(IEPRO.LE.0)GO TO 71  ! removed by BPz 20/06-2012
C
C        COMP. OF PRESSURE, DENSITY AND INNER ENERGY
C
      XIH=XIONG(1,1)-DXI
      XIHM=0.747-2.*DXI
C        XIH AND XIHM ARE THE ENERGIES OF IONIZATION FOR H AND H- RESPECTIVELY
C        (IN ELECTRON VOLTS).
C
      XKHM=TETA25*2.*EXP10(-TETA*XIHM)

C        XKHM = THE 'DISSOCIATION CONSTANT' FOR H-.
C
      HJONH=ANJON(1,2)/ANJON(1,1)
      IF(T.GT.TMOLIM)GO TO 42
      IF(MOLH.LE.0) GOTO 45
C
C        FORMATION OF MOLECULES. ONLY H2 AND H2+
   41 CALL MOLEQ(T,PE,HJONH,XIH,XKHM,XIHM,XNENH,F1,F2,F3,F4,F5,FE,FSUM,
     *   EH)
      FEPE=PE/FE
      F1P=F1*FEPE
      F3P=F3*FEPE
      F4P=F4*FEPE
      F5P=F5*FEPE
      phydro=fsum*pe/fe
      GO TO 43
C        FORMATION OF MOLECULES COMPOSED OF H,C,N,O
* to be removed or changed???
   45 IF(ANJON(3,1).LE.0..OR.ANJON(4,1).LE.0..OR.ANJON(5,1).LE.0.)
     * GOTO 41
      HJONC=ANJON(3,2)/ANJON(3,1)
      HJONN=ANJON(4,2)/ANJON(4,1)
      HJONO=ANJON(5,2)/ANJON(5,1)
      ABUC=ABUND(3)/ABUND(1)
      ABUN=ABUND(4)/ABUND(1)
      ABUO=ABUND(5)/ABUND(1)
** might try a better guess, if Pg already known.
      if (pgpgpg.le.0.) then
        pgin= 1.
      else
        pgin=pgpgpg
      endif
*
* is the T, P point close to previous one?

*****      if ((abs((t-tp)/t).lt.1.e-1).and.(abs((pe-pep)/pe).lt.0.7)) 
      if ((abs((t-tp)/t).lt.1.e-1).and.(abs((pe-pep)/pe).lt.1.0)) 
     &    then
* try saving time in die_pe  BPz 08/07-1999
        skiprelim=.true.
      else
        skiprelim=.false.
      endif
      call eqmol_pe(t,pgin,pg,pe,xih,xihm,kk,niter,skiprelim)

      if (pg.le.0.) then
        print*,' Jon: eqmol not converged after ',niter,' iterations!'
        print*,'  T=',t,' Pe=',pe
        if (iepro.gt.0) then
* we stop on non-convergence condition
          stop 'ERROR !!'
        else if (iepro.lt.0) then
* we go back 
          print*,'jon: We try again with other initial guess values'
          return
        else
          print*,'Jon, bad condition on iepro'
          stop 'ERROR !!'
        endif
      endif
CCCC      CALL MOL(T,PE,HJONH,HJONC,HJONN,HJONO,ABUC,ABUO,ABUN,XIH,XKHM,XIHM
CCCC     *,XNECNO,F1,F2,F3,F4,F5)
      xmy=xmytsuji
      enamn=eev/(xmytsuji*xmh)
      ejon=ejontsuji
      SUMPMO=0.
      PRESMO(1)=FHE*PK(1)
      PRESMO(2)=FHE*FHE*PK(2)
      PRESMO(3)=FHE*FHE*HJONH*PK(3)
      PRESMO(4)=FHE*FHE*FOE*PK(4)
      PRESMO(5)=FHE*FOE*PK(5)
      PRESMO(6)=FHE*FCE*PK(6)
      PRESMO(7)=FCE*FOE*PK(7)
      PRESMO(8)=FCE*FNE*PK(8)
      PRESMO(9)=FCE*FCE*PK(9)
      PRESMO(10)=FNE*FNE*PK(10)
      PRESMO(11)=FOE*FOE*PK(11)
      PRESMO(12)=FNE*FOE*PK(12)
      PRESMO(13)=FNE*FHE*PK(13)
      PRESMO(14)=FCE*FCE*FHE*FHE*PK(14)
      PRESMO(15)=FHE*FCE*FNE*PK(15)
      PRESMO(16)=FCE*FCE*FHE*PK(16)
      PRESMO(17)=0.0
      PRESMO(18)=FHE*FSE*PK(18)
      PRESMO(19)=FKE*FHE*PK(19)
      PRESMO(20)=FCE*FCE*FCE*FHE*PK(20)
      PRESMO(21)=FCE*FCE*FCE*PK(21)
      PRESMO(22)=FCE*FSE*PK(22)
      PRESMO(23)=FKE*FCE*PK(23)
      PRESMO(24)=FKE*FCE*FCE*PK(24)
      PRESMO(25)=FNE*FSE*PK(25)
      PRESMO(26)=FKE*FNE*PK(26)
      PRESMO(27)=FKE*FOE*PK(27)
      PRESMO(28)=FSE*FOE*PK(28)
      PRESMO(29)=FSE*FSE*PK(29)
      PRESMO(30)=FKE*FSE*PK(30)
      DO 30 I=1,NMOL
      PRESMO(I)=PRESMO(I)*PE
   30 SUMPMO=SUMPMO+PRESMO(I)
      SUMPA=PE*(FHE+FCE+FNE+FOE)
      SUMPI=PE*(FHE*HJONH+FCE*HJONC+FNE*HJONN+FOE*HJONO)
      HNIC=PE*FHE
      HPNIC=HNIC*HJONH
CCCC      PG=PE+SUMPMO+SUMPA+SUMPI+PE*SUMM/FE
      pgpgpg=pg
      phydro=hnic/f1
      ro=rho
*
ccc      print*,'jon, T,Pe: ',t,pe
ccc      print*,'jon: old ejon:',ejon
ccc      print*,'jon: new ejon:',ejontsuji
*
***************************
      GOTO 46
C
C        NO MOLECULES
   42 F2=ANJON(1,2)
      FE=XNENH+F2
      F1=ANJON(1,1)
      F3=0.
      F4=0.
      F5=0.
      FSUM=1.
      EH=-XIH*F1
      phydro=fsum*pe/fe
   43 PG=PE*(1.+(FSUM+SUMH)/FE)
      RO=PE*XMY*(XMH/XKBOL)/(FE*T)
***** Partial pressures for T>TMOLIM
      HNIC=PE*F1/fe
      DO 31 I=1,30
        PRESMO(I)=0.0
   31 continue
***** END
   46 continue
      XYRHO=RO
      E=1.5*PG/RO+(EH+EJON)*ENAMN
ccc      print*,'jon pg/ro eh ejon ejontsuji e ',
ccc     &       1.5*pg/ro,eh*enamn,ejon*enamn,ejontsuji*enamn,e
      YYPG=PG
      YYRHO=RO
      YYE=E
C
      IF(IOUTR.LE.0)GO TO 71
C
C        **** PRINT-OUT ****
C
      WRITE(IWRIT,204)T,PE,PG,RO,E
      WRITE(IWRIT,201)
      WRITE(IWRIT,202)
      DO93 I=1,NEL
      NJP=NJ(I)
      WRITE(IWRIT,203)IEL(I),ABUND(I),(ANJON(I,J),J=1,NJP)
      WRITE(IWRIT,207)(PART(I,J),J=1,NJP)
   93 CONTINUE
      IF(T.GT.TMOLIM)GO TO 44
      IF(MOLH.LE.0)GOTO 47
      WRITE(IWRIT,205)F1P,F3P,F5P,F4P
      GO TO 71
   47 CONTINUE
      WRITE(IWRIT,208)HNIC,(PRESMO(I),I=1,13)
      GOTO 71
   44 WRITE(IWRIT,206)
C
   71 CONTINUE
c
      pep=pe
      tp=t
      pgsave=pg
      rosave=ro
      esave=e
      kkp=kk

      RETURN
  201 FORMAT(1H0,'ELEMENT  ABUNDANCE  IONIZATION FRACTIONS',17X,
     *'PARTITION FUNCTIONS')
  202 FORMAT(1H ,23X,1HI,7X,2HII,6X,3HIII,5X,2HIV,12X,1HI,9X,2HII,8X,
     *3HIII,7X,2HIV)
  203 FORMAT(6H      ,A2,E12.4,4F8.4)
  204 FORMAT(3H0T=,F7.1,5X,3HPE=,E12.4,5X,3HPG=,E12.4,5X,3HRO=,E12.4,
     *5X,2HE=,E12.4)
  205 FORMAT(1H0,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',7X,'H2+'/1X,
     *4(1PE9.2))
  206 FORMAT(1H0,'NO MOLECULES CONSIDERED')
  207 FORMAT(1H+,56X,4E10.3)
  208 FORMAT(1H0,'PARTIAL PRESSURES'/4X,'H',8X,'H-',7X,'H2',7X,'H2+',6X,
     *'H2O',6X,'OH',7X,'CH',7X,'CO',7X,'CN',7X,'C2',7X,'N2',7X,'O2',7X,
     *'NO',7X,'NH'/1X,14(1PE9.2))
      END
