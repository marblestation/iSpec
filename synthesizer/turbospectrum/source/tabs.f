C
      SUBROUTINE TABS(NT,T)
C
C        THIS ROUTINE COMPUTES  FACTORS FOR INTERPOLATION IN T (TETA IF
C        ITETA(KOMP) IS GREATER THAN ZERO) IN THE ABKOF TABLE, INITIATED BY
C        SUBROUTINE INABS. CONCERNING THE OTHER CONTROL INTEGERS, SEE INABS.
C        THE RESULTING FACTORS ARE PUT IN AFAK. THE NUMBER OF FACTORS FOR
C        THE COMPONENT KOMP AT TEMPERATURE T(NTP) IS GIVEN IN
C        NOFAK((NKOMP-KOMPR)*(NTP-1)+KOMP-KOMPR). HERE KOMPR IS THE NUMBER
C        OF COMPONENTS WITH T-INDEP. COEFFICIENTS. NOFAK=0 MEANS THAT THE
C        ABSORPTION COEFFICIENT SHOULD BE =0. NPLATS (INDEX AS FOR NOFAK)
C        GIVES THE ARRAY INDEX OF THE TEMPERATURE POINT AT WHICH THE
C        INTERPOLATION IN ABKOF SHOULD START.
C
C        NT=NUMBER OF TEMPERATURES
C        T= ARRAY OF TEMPERATURES
C
C        DIMENSIONS NECESSARY
C        AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM),T(1)
C        THE DIMENSIONS ARE LOWER LIMITS. DIMENSIONS OF ARRAYS IN COMMON /CA1/
C        AND /CA2/ ARE COMMENTED ON IN SUBROUTINE INABS.
C        IFADIM SHOULD BE AT LEAST =(NKOMP-KOMPR)*NT, WHERE NKOMP IS THE NUMBER
C               OF COMPONENTS, KOMPR THE NUMBER OF TEMPERATURE-INDEPENDENT
C               COMPONENTS AND NT THE NUMBER OF TEMPERATURE POINTS (IN THE PARA-
C               METER LIST).
C        KFADIM SHOULD BE AT LEAST =KOMPR*NT+(NKOMP-KOMPR)*NT*NUM, WHERE NUM IS
C               BETWEEN 2 AND 3 AND DEPENDENT ON THE TYPE OF TEMPERATURE
C               INTERPOLATION USED.
C
C
      include 'spectrum.inc'
C
C      PARAMETER (IFADIM=1000,KFADIM=4000)
      DIMENSION T(NT)
      COMMON/UTPUT/IREAD,IWRIT
      COMMON/CA1/DELT(mkomp,2),TBOT(mkomp,2),IDEL(mkomp),ISVIT(mkomp),
     &           ITETA(mkomp),KVADT(mkomp),MAXET(mkomp),MINET(mkomp),
     &           NTM(mkomp,2),NEXTT,NUTZT
      COMMON/CA2/ABKOF(nabdim),KOMPLA(mkomp*20),KOMPR,KOMPS,NKOMP
      COMMON/CA4/AFAK(KFADIM),NOFAK(IFADIM),NPLATS(IFADIM)
      logical goon
C
      IFAK=1
      KFAK=1
      NSVIT=1
C        THIS IS JUST A DUMMY STATEMENT TO GIVE NSVIT A FORMAL VALUE
C
      do NTP=1,NT
        TP=T(NTP)
        KFAK=KFAK+KOMPR
        do KOMP=KOMPS,NKOMP
          if (ISVIT(KOMP).le.0) then
            if (ITETA(KOMP).gt.0) then
              TS=5040./T(NTP)
            else
              TS=T(NTP)
            endif
C SEARCHING
            if ((TS-TBOT(KOMP,1)).lt.0.) then
              if (MINET(KOMP).LE.0) then
C No extrapolation downwards
                nsvit=3
              else
C EXTRAPOLATION DOWNWARDS
                if (NEXTT.GT.0) WRITE(IWRIT,200)TS,KOMP
                INTA=1
                AP=(TS-TBOT(KOMP,1))/DELT(KOMP,1)
                IP=0
                A2=AP-FLOAT(IP)
                A1=1.-A2
                nsvit=2
              endif
C SEARCHING CONTINUES
            else
              INTAP=1
              IDP=IDEL(KOMP)
              goon=.true.
              I=1
              do while (goon.and.I.le.idp)
                AP=(TS-TBOT(KOMP,I))/DELT(KOMP,I)
                IP=INT(AP)
                INTA=IP+INTAP
                INAP=NTM(KOMP,I)-1+INTAP
C inta is the index just before the T we want.
C inap is the index of the beginning of the last possible T-interval with this spacing.
                if (INTA.LE.INAP) then
                  goon=.false.
                  if (KVADT(KOMP).LE.0.or.NTM(KOMP,I).eq.1) then
                    A2=AP-FLOAT(IP)
                    A1=1.-A2
                    nsvit=2
                  else
C QUADRATIC INTERPOLATION
                    if (INTA.ge.INAP) then
C This happens if our T is in the last T-interval for this spacing.
C We take the T-point one step below to allow quadratic interpolation.
C The temperature we interpolate to is then in the second T-interval.
C We cannot allow this if there is only one interval with that spacing.
C This used to cause an erroneous interpolation in the past (inta=0, 
C and IP=-1, or a quadratic interpolation using points with different
C spacings). This is why I have added ".or.NTM(KOMP,I).eq.1" above.
C This forces linear interpolation.
C       BPz 12/07-2001
                      INTA=INTA-1
                      IP=IP-1
                    endif
C
                    DXX1=AP-FLOAT(IP)
                    DXX2=DXX1-1.
                    DXX3=DXX1-2.
                    A1=DXX2*DXX3*0.5
                    A2=-DXX1*DXX3
                    A3=DXX1*DXX2*0.5
                    nsvit=1
                  endif
                else
                  INTAP=INAP+1
                endif
                I=I+1
              enddo
              if (goon) then
                if (MAXET(KOMP).LE.0) then
                  nsvit=3
                else 
                  if (NEXTT.GT.0) WRITE(IWRIT,200)TS,KOMP
C
C        EXTRAPOLATION UPWARDS
                  INTA=INAP
                  IP=NTM(KOMP,IDP)-1
                  A2=AP-FLOAT(IP)
                  A1=1.-A2
                  nsvit=2
                endif
              endif
            endif
          endif
C
          if (nsvit.eq.1) then

            AFAK(KFAK)=A1
            AFAK(KFAK+1)=A2
            AFAK(KFAK+2)=A3
            NPLATS(IFAK)=INTA
            NOFAK(IFAK)=3
            IFAK=IFAK+1
            KFAK=KFAK+3
            NSVIT=1

          else if (nsvit.eq.2) then

            AFAK(KFAK)=A1
            AFAK(KFAK+1)=A2
            NPLATS(IFAK)=INTA
            NOFAK(IFAK)=2
            IFAK=IFAK+1
            KFAK=KFAK+2
            NSVIT=2
C
          else if (nsvit.eq.3) then
C OUTSIDE TABLE. ABS.COEFF. SHOULD BE = 0
            IF(NUTZT.GT.0)WRITE(IWRIT,201)TS,KOMP
            NOFAK(IFAK)=0
            IFAK=IFAK+1
            NSVIT=3
          endif
C
          if (KFAK.GT.KFADIM) then
             WRITE(IWRIT,202)KFAK,KFADIM,NT
             STOP 'TABS 1'
          endif
          if (IFAK.GT.IFADIM+1) then
             WRITE(IWRIT,203)IFAK,IFADIM,NT
             STOP 'TABS 2'
          endif
        enddo
      enddo
C
      RETURN
C
  200 FORMAT('EXTRAPOLATION IN TABS, T (TETA)=',E12.5,5X,
     &       ' COMPONENT NO ',I5)
  201 FORMAT('ZERO IN TABS, T (TETA)=',E12.5,5X,' COMPONENT NO',I5)
  202 FORMAT(' KFAK=',I5,5X,' GT KFADIM=',I5,5X,'IN TABS, NT=',I5)
  203 FORMAT(' IFAK=',I5,5X,' GT IFADIM=',I5,5X,'IN TABS, NT=',I5)
      END
