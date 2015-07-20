C
      SUBROUTINE INABS(IOUTS)
C
C        THIS ROUTINE  READS ABSORPTION COEFFICIENT TABLES AND INTER/EXTRA-
C        POLATES THEM TO OUR WAVELENGTHS GIVEN IN XL. THE INTERPOLATION IS
C        PERFORMED SEPARATELY FOR EACH WAVELENGTH SET.
C
C        NKOMP IS THE NUMBER OF COMPONENTS IN THE FULL TABLE.
C        NEXTL SHOULD BE GREATER THAN ZERO IF A PRINT-OUT IS WANTED ON EXTRA-
C              POLATION IN WAVELENGTH,
C        NUTZL IF PRINT-OUT IS WANTED WHEN WE PUT THE COEFFICIENT =0 OUTSIDE THE
C              WAVELENGTH REGION OF THE TABLES.
C        NEXTT AND NUTZL ARE THE CORRESPONDING QUANTITIES ON INTERPOLATION IN
C              T, MADE IN SUBROUTINE TABS.
C        NULL  SHOULD BE GREATER THAN ZERO IF A PRINT-OUT IS WANTED (FROM SUB-
C              ROUTINE ABSKO) WHEN A COEFFICIENT IS FOUND TO BE LESS THAN ZERO
C              ON INTERPOLATION IN T AND THEREFORE PUT EQUAL TO ZERO.
C
C        FOR EACH COMPONENT THE FOLLOWING PARAMETERS MUST BE SPECIFIED
C        ABNAME IS THE NAME OF, OR A SYMBOL FOR, THE ABSORPTION MECHANISM.
C        SOURCE INDICATES THE SOURCE OR REFERENCE OF THE DATA
C
C        1. PARAMETERS FOR THE WAVELENGTH INTERPOLATION.
C          ILOGL SHOULD BE GREATER THAN ZERO IF INTERPOLATION IN WAVELENGTH IS
C                TO BE PERFORMED ON THE LOGARITHMIC ABSORPTION COEFFICIENTS
C                (WITH SUBSEQUENT EXPONENTIATION OF THE RESULTS - HERE IF ILOGT
C                IS EQUAL TO ZERO OR IN SUBROUTINE ABSKO IF ILOGT IS GREATER
C                THAN ZERO). OTHERWISE INTERPOLATION IN WAVELENGTH IS MADE
C                DIRECTLY ON THE ABSORPTION COEFFICIENTS THEMSELVES.
C          KVADL SHOULD BE GREATER THAN ZERO IF QUADRATIC INTERPOLATION IN
C                WAVELENGTH IS WANTED. OTHERWISE INTERPOLATION WILL BE LINEAR
C          MINEX SHOULD BE GT 0 IF LINEAR EXTRAPOLATION (INSTEAD OF PUTTING THE
C                COEFFICIENT = 0) IS WANTED TOWARDS SHORTER WAVELENGTHS.
C          MAXEX, CORRESPONDING TOWARDS LONGER WAVELENGTHS.
C          NLATB IS THE NUMBER OF WAVELENGTH POINTS OF THE ABSORPTION COEFFI-
C                CIENT TABLE TO BE READ.
C          XLATB ARE THOSE WAVELENGTHS. THEY SHOULD BE GIVEN IN INCREASING ORDER
C
C        2. PARAMETERS FOR THE TEMPERATURE INTERPOLATION.
C          ILOGT, KVADT, MINET, MAXET AND NTETB ARE THE T-INTERPOLATION
C                ANALOGUES TO ILOGL-NLATB.
C          ITETA IS PUT GREATER THAN ZERO WHEN TETA VALUES (TETA=5040./T) ARE
C                GIVEN IN XTET INSTEAD OF TEMPERATURES.
C          XTET ARE THE TEMPERATURE (TETA) VALUES OF THE ABSORPTION
C                COEFFICIENT TABLE TO BE READ. THE XTET VALUES SHOULD BE GIVEN
C                IN INCREASING ORDER AND EQUIDISTANTLY, HOWEVER (IELMAX-1)
C                CHANGES OF THE INTERVAL ARE ALLOWED. THE PROGRAM CHECKS  THAT
C                THIS NUMBER IS NOT EXCEEDED.
C        XKAP IS THE ABSORPTION COEFFICIENT TABLE FOR THE ACTUAL COMPONENT. THE
C                WAVELENGTHS INCREASES MORE RAPIDLY THAN T (TETA).
C
C        THE TABLES FOR  T E M P E R A T U R E - I N D E P E N D E N T
C        C O M P O N E N T S  S H O U L D  B E  P U T  F I R S T .
C           THE RESULTING TABLE IS PUT IN ABKOF. HERE T (TETA) INCREASES MORE
C        RAPIDLY THAN XLA, WHICH INCREASES MORE RAPIDLY THAN KOMP. IF THE RESULT
C        OF THE INTERPOLATION IS ZERO FOR A CERTAIN XLA(J) AND KOMP, THIS IS NOT
C        PUT IN ABKOF. INSTEAD A NOTE IS MADE IN KOMPLA (KOMPLA(NLB*(KOMP-
C        1)+J) IS PUT EQUAL TO ZERO). OTHERWISE THE KOMPLA VALUE TELLS WHERE IN
C        ABKOF THE TABLE FOR THE COMPONENT KOMP AND THE WAVELENGTH J BEGINS.
C
C        A DETAILED PRINT-OUT IS GIVEN IF IOUTS IS GREATER THAN ZERO.
C
C
C        DIMENSIONS NECESSARY
C        ABKOF(NABDIM),ABNAME(NKOMP),DELT(NKOMP,IELMAX),IDEL(NKOMP),
C        IDISKV(MAX(NLATB)),ILOGTA(NKOMP),IRESET(NSET),ISVIT(NKOMP),ITETA(NKOMP)
C        KOMPLA(MAX(NL)*NKOMP),KVADT(NKOMP),MAXET(NKOMP),MINET(NKOMP),
C        NL(NSET),NTAET(NKOMP),NTM(NKOMP,IELMAX),SOURCE(NKOMP),
C        TBOLT(NKOMP,IELMAX),XKAP(MAX(NLATB),MAX(NTETB)),XL(MAX(NL),NSET)
C        XLA(MAX(NL)),XLA3(MAX(NL)),XLATB(MAX(NLATB)),XTET(MAX(NTETB)),
C        XTETP(MAX(NTETB))
C
C        THE DIMENSIONS ARE LOWER LIMITS
C        IELMAX IS THE MAXIMUM NUMBER OF DIFFERENT T INTERVALS (GIVEN BELOW) IN
C               ANY ABSORPTION COEFFICIENT TABLE.
C        NABDIM IS THE DIMENSION OF THE ABKOF ARRAY (GIVEN BELOW).
C        NKOMP IS THE NUMBER OF 'COMPONENTS', I.E. EQUAL TO THE NUMBER OF
C               DIFFERENT ABSORPTION COEFFICIENT TABLES TO BE READ.
C        NL(I)  IS THE NUMBER OF WAVELENGTHS IN THE WAVELENGTH SET I.
C        NLATB(KOMP) IS THE NUMBER OF WAVELENGTH POINTS IN THE TABLE TO BE READ
C               FOR THE COMPONENT KOMP.
C        NSET   IS THE NUMBER OF WAVELENGTH SETS.
C        NTETB  IS THE NUMBER OF TEMPERATURE POINTS IN THE TABLE FOR THE COM-
C               PONENT BEING CONSIDERED.
C
C
      include 'spectrum.inc'

      CHARACTER*8 ABNAME,SOURCE
      DIMENSION IDISKV(40),XLATB(40),XTET(30),NTAET(mkomp),XKAP(40,30),
     *XLA3(20),XLA(20),XTETP(30)
      COMMON/CARC4/PROV(mkomp),NDUM(3)
      COMMON /CHAR/ ABNAME(mkomp),SOURCE(mkomp)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA1/DELT(mkomp,2),TBOT(mkomp,2),IDEL(mkomp),ISVIT(mkomp),
     &           ITETA(mkomp),KVADT(mkomp),MAXET(mkomp),MINET(mkomp),
     &           NTM(mkomp,2),NEXTT,NUTZT
      COMMON/CA2/ABKOF(nabdim),KOMPLA(mkomp*20),KOMPR,KOMPS,NKOMP
      COMMON/CA3/ILOGTA(mkomp),NULL
      COMMON/CFIL/IRESET(numbset),ISLASK,IREAT
      COMMON/CXLSET/NSET,NL(numbset),XL(20,numbset)
      real*8 vacair,lambda

      logical quadratic, linear
C
C        IELMAX IS THE MAXIMUM NUMBER OF DIFFERENT T INTERVALS IN THE XKAP-
C        TABLE. THE DIMENSIONS OF TBOT, DELT AND NTM ARE AFFECTED BY THIS NUMBER
      IELMAX=2
C        THE DIMENSION OF THE ABKOF ARRAY

      do L=1,30
        XTETP(L)=0.
      enddo
      if (IOUTS.GT.0) WRITE(IWRIT,229)
C
      READ(IREAT,101) NKOMP,NEXTL,NUTZL,NEXTT,NUTZT,NULL
C
cc           print *,'inabs,    NKOMP=',nkomp
C
      KOMPR=0
      REWIND ISLASK
C
C        LOOP OVER COMPONENTS STARTS (THE 'FIRST KOMP-LOOP')
      do KOMP=1,NKOMP
        READ(IREAT,105) ABNAME(KOMP),SOURCE(KOMP)
        READ(IREAT,102) ILOGL,KVADL,MINEX,MAXEX,NLATB
        READ(IREAT,103) (XLATB(J),J=1,NLATB)
C
C        WE FIND THE DISCONTINUITIES IN WAVELENGTH
C        A DISCONTINUITY IN A TABLE IS DEFINED BY TWO WAVELENGTH POINTS
C        WITHIN LESS THAN TWO ANGSTROEMS.
        IDISK=0
        IDISKV(1)=0
        DO J=2,NLATB
          IDISKV(J)=0
          IF((XLATB(J)-XLATB(J-1)).lt.2.) then
            IDISKV(J-1)=1
            IDISKV(J)=1
            IDISK=1
          endif
        enddo
*
* shift wavelengths to air wavelengths, from the vacuum values of the MARCS jonabs file
*
        do j=1,nlatb
          if (xlatb(j).gt.2000.) then
            lambda=dble(xlatb(j))
            lambda=vacair(lambda)
            xlatb(j)=sngl(lambda)
          endif
        enddo
C
C        CONTINUE READING
        READ(IREAT,102) ILOGT,KVADT(KOMP),MINET(KOMP),MAXET(KOMP),NTETB,
     &     ITETA(KOMP)
C
cc          print*, 'inabs    NTETB= ',ntetb
C
        ILOGTA(KOMP)=ILOGT
        if (NTETB.GT.1) then
          READ(IREAT,103)(XTET(L),L=1,NTETB)
        else
          KOMPR=KOMPR+1
        endif
C
C        FINALLY THE ABSORPTION COEFFICIENT TABLE IS READ
        do K=1,NTETB
          READ(IREAT,104)(XKAP(JJ,K),JJ=1,NLATB)
        enddo
C
C        WE TAKE THE LOGARITHMS BEFORE THE WAVELENGTH INTERPOLATION
C        IF ILOGL IS GREATER THAN ZERO.
        if (ILOGL.ge.1) then
          do K=1,NTETB
            do JJ=1,NLATB
              if (XKAP(JJ,K).le.0.) then
C
C A COEFFICIENT FOR WHICH THE LOGARITHM SHOULD BE TAKEN IS ZERO
                WRITE(IWRIT,207)JJ,K,XKAP(JJ,K),KOMP
                XKAP(JJ,K)=1.E-37
              endif
              XKAP(JJ,K)=ALOG(XKAP(JJ,K))
            enddo
          enddo
        endif
C
C        PREPARATION OF THE T-INTERPOLATION IN SUBROUTINE TABS
C
C        WE FIND OUT WHETHER ISVIT(KOMP) CAN BE CHOSEN GREATER THAN ZERO. THIS
C        IS THE CASE IF THE T SCALE AND MAXET, MINET AND KVADT ARE IDENTICAL
C        WITH THOSE OF THE PREVIOUS COMPONENT. IF ISVIT IS GREATER THAN ZERO
C        THE TIME SPENT IN SUBR. TABS WILL BE DECREASED.
        ISVIT(KOMP)=0
        if (NTETB.gt.1) then
          if ((NTETB.eq.NTETBP).and.(MAXET(KOMP).eq.MAXETP).and.
     &     (MINET(KOMP).eq.MINETP).and.(KVADT(KOMP).eq.KVADTP)) then
            ISVIT(KOMP)=1
          endif
          do L=1,NTETB
            IF(XTET(L).NE.XTETP(L)) isvit(komp)=0
          enddo
C
C          WE REMEMBER TEMPERATURES ETC. FOR NEXT COMPONENT
          do L=1,NTETB
            XTETP(L)=XTET(L)
          enddo
          NTETBP=NTETB
          MAXETP=MAXET(KOMP)
          MINETP=MINET(KOMP)
          KVADTP=KVADT(KOMP)
C
C          WE FIND THE INTERVALS IN THE T (TETA) SCALE
          TBOT(KOMP,1)=XTET(1)
          DELT(KOMP,1)=XTET(2)-XTET(1)
          NTM(KOMP,1)=1
          IDEL(KOMP)=1
          if (NTETB.gt.2) then
C
            J=1
            LF=1
            do L=3,NTETB
*
*  PRobleme ici avec le decompte des LF ?????
*
              DIFF=XTET(L)-XTET(L-1)
              if (ABS(1.-DIFF/DELT(KOMP,J)).ge.1.E-4) then
                J=J+1
                if (J.GT.IELMAX) then
C TOO MANY DIFFERENT INTERVALS IN THE T-TABLE FOR THIS COMPONENT
                  WRITE(IWRIT,203)KOMP,IELMAX
                  WRITE(IWRIT,206)(XTET(LL),LL=1,NTETB)
                  STOP 'ERROR INABS 1: too many components in T-table!'
                endif
                TBOT(KOMP,J)=XTET(L-1)
                DELT(KOMP,J)=DIFF
                NTM(KOMP,J-1)=LF
                LF=0
              endif
              LF=LF+1
            enddo
            NTM(KOMP,J)=LF
            IDEL(KOMP)=J
          endif
        endif
* added 20/05-2008 by BPz
        ntetbp=ntetb
C
        NTAET(KOMP)=NTETB
C        ALL DATA NECESSARY BELOW FOR THIS COMPONENT ARE STORED ON UNIT
C        ISLASK
        WRITE(ISLASK)KVADL,MINEX,MAXEX,NLATB,ILOGL,IDISK,(IDISKV(J),J=1,
     &    NLATB),(XLATB(J),J=1,NLATB),NTETB,ILOGT,(XTET(L),L=1,NTETB),
     &    ((XKAP(JJ,K),JJ=1,NLATB),K=1,NTETB)
C
C        **** PRINT-OUT ****
        if (IOUTS.gt.0) then
          WRITE(IWRIT,211)KOMP,ABNAME(KOMP),SOURCE(KOMP),ABNAME(KOMP)
          WRITE(IWRIT,212)
          WRITE(IWRIT,213)XLATB(1),XLATB(NLATB)
          IF(MINEX.EQ.0)WRITE(IWRIT,214)
          IF(MINEX.GT.0)WRITE(IWRIT,215)
          IF(KVADL.EQ.0)WRITE(IWRIT,216)
          IF(KVADL.GT.0)WRITE(IWRIT,217)
          IF(ILOGL.EQ.0)WRITE(IWRIT,218)
          IF(ILOGL.GT.0)WRITE(IWRIT,219)
          IF(MAXEX.EQ.0)WRITE(IWRIT,220)
          IF(MAXEX.GT.0)WRITE(IWRIT,221)
          IF(IDISK.GT.0)WRITE(IWRIT,222)
          if (NTETB.eq.0.or.ntetb.eq.1) then
            WRITE(IWRIT,230)
          else 
            WRITE(IWRIT,223)
            WRITE(IWRIT,213)XTET(1),XTET(NTETB)
            IF(MINET(KOMP).EQ.0)WRITE(IWRIT,214)
            IF(MINET(KOMP).GT.0)WRITE(IWRIT,215)
            IF(KVADT(KOMP).EQ.0)WRITE(IWRIT,216)
            IF(KVADT(KOMP).GT.0)WRITE(IWRIT,217)
            IF(ILOGTA(KOMP).EQ.0)WRITE(IWRIT,218)
            IF(ILOGTA(KOMP).GT.0)WRITE(IWRIT,219)
            IF(MAXET(KOMP).EQ.0)WRITE(IWRIT,220)
            IF(MAXET(KOMP).GT.0)WRITE(IWRIT,221)
            IF(ISVIT(KOMP).GT.0)WRITE(IWRIT,224)
            WRITE(IWRIT,*)
          endif
        endif
      enddo
C        END OF 'THE FIRST KOMP-LOOP'
C
      KOMPS=KOMPR+1
C
C
C        WE BUILD THE ABKOF ARRAY. INTERPOLATION IN WAVELENGTH.
C
C        LOOP OVER WAVELENGTH SETS ('THE ISET-LOOP')
      do ISET=1,NSET
        REWIND ISLASK
        NLB=NL(ISET)
        do J=1,NLB
          XLA(J)=XL(J,ISET)
          XLA3(J)=XLA(J)**3
        enddo
        INDEX=1
C
C        LOOP OVER COMPONENTS STARTS ('THE SECOND KOMP-LOOP')
        do KOMP=1,NKOMP
          READ(ISLASK) KVADL,MINEX,MAXEX,NLATB,ILOGL,IDISK,(IDISKV(J),
     &                 J=1,NLATB),(XLATB(J),J=1,NLATB),NTETB,ILOGT,
     &                 (XTET(L),L=1,NTETB),
     &                 ((XKAP(JJ,K),JJ=1,NLATB),K=1,NTETB)
          JI=1
C
C LOOP OVER WAVELENGTHS ('THE J-LOOP') STARTS
          do J=1,NLB
C SEARCHING IN WAVELENGTH
            IU=NLB*(KOMP-1)+J
            KOMPLA(IU)=INDEX
            if (XLA(J).lt.XLATB(1)) then
              ihelp=1
              lambi=1
            else if (XLA(J).ge.XLATB(nlatb)) then
              ihelp=nlatb
              lambi=nlatb
            else
              do JJ=1,NLATB-1
                if (XLA(J).ge.XLATB(jj).and.xla(j).lt.xlatb(jj+1)) then
                  lambi=jj
                  ihelp=jj+1
                endif
              enddo
            endif
            linear=.false.
            quadratic=.false.
            if (IHELP.le.1) then
              if (MINEX.gt.0) then
                IF(NEXTL.GT.0)WRITE(IWRIT,201)KOMP,XLA(J)
                linear=.true.
              else
C too small a wavelength - outside the table
C abs. coeff. is put = zero
                KOMPLA(IU)=0
                if (NUTZL.GT.0)  WRITE(IWRIT,202)KOMP,XLA(J)
              endif
            else if (lambi.eq.nlatb) then
              if (maxex.le.0) then
C too great a wavelength - outside the table
C abs. coeff. is put = zero
                KOMPLA(IU)=0
                if (NUTZL.GT.0)  WRITE(IWRIT,202)KOMP,XLA(J)
              else
                lambi=lambi-1
                linear=.true.
              endif
            else if (lambi.lt.nlatb) then
              if (kvadl.lt.1) then
                linear=.true.
              else if (kvadl.ge.1) then
                lambi=min(nlatb-2,lambi)
                if ((idisk.le.0).or.(idiskv(lambi+1).le.0)) then
                  quadratic=.true.
                else if (xla(j).gt.xlatb(lambi+1)) then
                  lambi=lambi+1
                  if ((IDISKV(LAMBI+1).GT.0).or.(LAMBI+1.EQ.NLATB)) then
                    linear=.true.
                  else
                    quadratic=.true.
                  endif
                else if ((idiskv(lambi).gt.0).or.(lambi.eq.1)) then
                  linear=.true.
                else 
                  lambi=lambi-1
                  quadratic=.true.
                endif
              endif
            endif
          
            if (linear) then
C LINEAR INTER- AND EXTRAPOLATION
              A2=(XLA(J)-XLATB(LAMBI))/(XLATB(LAMBI+1)-XLATB(LAMBI))
              A1=1.-A2
              do K=1,NTETB
                ABKOF(INDEX)=A1*XKAP(LAMBI,K)+A2*XKAP(LAMBI+1,K)
                INDEX=INDEX+1
              enddo
            else if (quadratic) then
C QUADRATIC INTERPOLATION
              DXX1=XLA(J)-XLATB(LAMBI)
              DXX2=XLA(J)-XLATB(LAMBI+1)
              DXX3=XLA(J)-XLATB(LAMBI+2)
              DX21=XLATB(LAMBI+1)-XLATB(LAMBI)
              DX32=XLATB(LAMBI+2)-XLATB(LAMBI+1)
              DX31=XLATB(LAMBI+2)-XLATB(LAMBI)
              A1=DXX2*DXX3/(DX21*DX31)
              A2=DXX1*DXX3/(DX21*DX32)
              A3=DXX1*DXX2/(DX31*DX32)
              do K=1,NTETB
                ABKOF(INDEX)=A1*XKAP(LAMBI,K)-A2*XKAP(LAMBI+1,K)+A3*
     &                        XKAP(LAMBI+2,K)
                INDEX=INDEX+1
              enddo
            endif
      
            if (linear.or.quadratic) then
              if (ilogl.ge.1) then
                if (ilogt.le.0) then
C LOGARITHMIC INTERPOLATION ONLY IN WAVELENGTH
                  LIP=INDEX-NTETB
                  LAP=INDEX-1
                  do LL=LIP,LAP
                    ABKOF(LL)=EXP(ABKOF(LL))
                  enddo
                endif
              else
                if (ILOGT.gt.0) then
C WE TAKE THE LOGARITHM BEFORE THE T INTERPOLATION IF ILOGT GT 0
                  LIP=INDEX-NTETB
                  LAP=INDEX-1
                  do LL=LIP,LAP
                    if (ABKOF(LL).le.0.) then
C IMPOSSIBLE TO TAKE THE LOGARITHM OF A NEGATIVE COEFFICIENT
                      LUS=LL-LIP+1
                      WRITE(IWRIT,208) LL,ABKOF(LL),KOMP,J,ISET,LUS
                      ABKOF(LL)=1.E-37
                    endif
                    ABKOF(LL)=ALOG(ABKOF(LL))
                  enddo
                endif
              endif
            endif
          enddo
C END OF 'THE J-LOOP'
        enddo
C END OF 'THE SECOND KOMP-LOOP'
C
C WRITE THE DATA OF THE SET ISET ON UNIT IRESET(ISET)
        NABKOF=INDEX-1
        NKOMPL=IU
        IREADP=IRESET(ISET)
        WRITE(IREADP)ISET,NLB,XLA,XLA3,NABKOF,ABKOF,NKOMPL,KOMPLA
C
cc ?? BPz 13/07-2001. next 2 statements were active in babsma, not in marcs.
        END FILE IREADP
        BACKSPACE IREADP
C
C        CHECK DIMENSION OF ABKOF
        if (IOUTS.GT.0) WRITE(IWRIT,204)NABKOF,ISET
        if (NABKOF.gt.NABDIM) then
C TOO SMALL DIMENSION FOR ABKOF
          WRITE(IWRIT,205)NABDIM
          STOP ' ERROR! INABS 2: abkof dimension too small!!'
        endif
      enddo
C
C        END OF 'THE ISET-LOOP'
C
      do ISET=1,NSET
        IREADP=IRESET(ISET)
        REWIND IREADP
      enddo
      if  (IOUTS.gt.0) then
C **** PRINT-OUT ****
C ON WAVELENGTH SETS AND FILES
        WRITE(IWRIT,225)IREAT,ISLASK
        WRITE(IWRIT,226)
        do M=1,NSET
          NP=NL(M)
          WRITE(IWRIT,227)M,IRESET(M)
          WRITE(IWRIT,228)(XL(J,M),J=1,NP)
        enddo
        WRITE(IWRIT,*)
      endif

      RETURN

  101 FORMAT(8X,I2,5(9X,I1))
  102 FORMAT(4(9X,I1),8X,I2,9X,I1)
  103 FORMAT(6F10.0)
  104 FORMAT(6E10.3)
  105 FORMAT(2A8)
  201 FORMAT(' EXTRAPOLATION FOR COMPONENT',I5,5X,'AND WAVELENGTH=',
     & F10.3,5X,'***INABS***')
  202 FORMAT(' ABS.COEFF. PUT =0 AT WAVELENGTH-INTER/EXTRAPOLATION FOR
     *COMPONENT',I5,5X,'AND WAVELENGTH=',F10.3,5X,'***INABS***')
  203 FORMAT(' TOO MANY DIFFERENT INTERVALS IN THE T-(TETA-)TABLE FOR CO
     *MPONENT ',I5,5X,'MAX IS',I5,5X,'***INABS***')
  204 FORMAT(' NECESSARY DIMENSION FOR ABKOF=',I5,5X,'IN SET',I5,5X,'***
     *INABS***')
  205 FORMAT(' DIMENSION ALLOWED =',I5,5X,'TOO SMALL     ***INABS***')
  206 FORMAT(' XTET=',10E12.4)
  207 FORMAT(' XKAP(',I2,',',I2,')=',E12.5,' PUT = 1.E-37 BEFORE LOG HAS
     *BEEN TAKEN.   ***INABS***')
  208 FORMAT(' ABKOF(',I4,')=',E12.5,' FOR COMPONENT ',I2,' WAVELENGTH',
     *I3,' SET ',I2,' XTET NR ',I2,' ABKOF PUT=1.E-37   ***INABS***')
  211 FORMAT('0************  COMPONENT NO',I3,', ',A8,' SOURCE  ',A8,'
     *****************',5X,A8)
  212 FORMAT('0     I N T E R P O L A T I O N  I N  W A V E L E N G T H'
     *)
  213 FORMAT(15X,F10.2,25X,F10.2)
  214 FORMAT('   KAPPA=0   - ')
  215 FORMAT(' LIN. EXTRAP.- ')
  216 FORMAT(25X,' -LIN. INT. ')
  217 FORMAT(25X,' -QUAD. INT.')
  218 FORMAT(37X,'  KAPPA    - ')
  219 FORMAT(37X,' LOG(KAPPA)- ')
  220 FORMAT(60X,' -   KAPPA=0')
  221 FORMAT(60X,' -LIN. EXTRAP.')
  222 FORMAT(' DISCONTINUITIES PRESENT.')
  223 FORMAT('0     I N T E R P O L A T I O N  I N  T  ( T E T A )')
  224 FORMAT('0 T SCALE ETC. IDENTICAL WITH PRECEEDING COMPONENT')
  225 FORMAT('F I L E S  U S E D  B Y  T H E  A B S - B L O C K'//'
     * INITIAL FILE ',I3,', PRELIMINARY FILE',I3)
  226 FORMAT('SET           WAVELENGTHS',81X,'FILE')
  227 FORMAT (I2,105X,I2)
  228 FORMAT(5X,10F10.2)
  229 FORMAT('D A T A  F R O M  S U B R O U T I N E  I N A B S')
  230 FORMAT('     N O  T-  ( T E T A - ) D E P E N D E N C E')
      END
