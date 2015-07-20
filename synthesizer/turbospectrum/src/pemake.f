C
      SUBROUTINE PEMAKE(T,PE,PG,PEX)
C
C 'PEMAKE-R' IS COMPATIBLE WITH 'PEMAKE' AND SOMEWHAT FASTER. A MODIFIED
C REGULA-FALSI PROCEDURE IS USED ON THE LOG-LOG PG-PE RELATION.
C 76.03.08  *NORD*
C YOU FEED IT WITH T AND PG AND A GUESS PE. YOU GET PEX
C
      DATA IT,N,EPS/0,20,1.E-3/
C
C START
c guess input Pg to help jon
      FA=pg
      CALL JON(T,PE,-1,FA,RO,E,0,1)

C try to improve convergence
      pefirst=pe
      do while (fa.lt.0..and.pe.gt.1.e-30)
c not converged in jon
        pe=pe/1.e1
        print*,'pemake, trying with lower Pe:',pe
        CALL JON(T,PE,-1,FA,RO,E,0,1)
      enddo
      if (fa.lt.0.) pe=pefirst
      do while (fa.lt.0..and.pe.lt.1000.)
c not converged in jon
        pe=pe*1.e1
        print*,'pemake, trying with higher Pe:',pe
        CALL JON(T,PE,-1,FA,RO,E,0,1)
      enddo
      if (fa.lt.0.) then
        print*,'giving up in pemake' 
        stop 'Too many unsuccessful tries.'
      endif

      A=ALOG(PE)
      PEX=PE
c
C******WRITE(7,40) T,PG,PE,FA
40    FORMAT(' T,PG,PE,PGP=',4E10.4)
      IT=IT+1
      FA=ALOG(FA/PG)
      IF(ABS(FA).LT.EPS) GOTO 101
      B=A-0.69*FA
      PEX=EXP(B)
C ONE PEMAKE ITERATION, CF. PEMAKE
C input pg to help jon
      FB=FA*PEX/PE
      CALL JON(T,PEX,1,FB,RO,E,0,1)
C******WRITE(7,40) T,PG,PE,FB
      IT=IT+1
      FB=ALOG(FB/PG)
      IF(ABS(FB).LT.EPS) GOTO 101
      X=B
C
C LOOP
      fx=fb
      DO 100 I=1,N
      XOLD=X
C
C INTERPOLATE TO FIND NEW X
      X=A-(B-A)/(FB-FA)*FA
      PEX=EXP(X)
      IF(ABS(X-XOLD).LT.EPS) GOTO 101
      FX=pex/exp(xold)*fx
      CALL JON(T,PEX,1,FX,RO,E,0,1)
C******WRITE(7,40) T,PG,PEX,FX
      IT=IT+1
      FX=ALOG(FX/PG)
C
C CHECK IF A OR B CLOSEST TO X
      IF(ABS(A-X).LT.ABS(B-X)) GOTO 102
      A=X
      FA=FX
      GOTO 100
102   B=X
      FB=FX
C
C END OF LOOP
100   CONTINUE
      WRITE(7,51) N,T,PE,PG,A,B,FA,FB,EPS
51    FORMAT('0***PEMAKE, MAX ITER.: N,T,PE,PG,A,B,FA,FB,EPS=',
     * /,1X,I2,8E11.4)
      RETURN
C
C NORMAL END
101   CONTINUE
      RETURN
C
C COUNT ENTRY
      ENTRY PECNT
      WRITE(7,52) IT
52    FORMAT('0TOTAL NUMBER OF CALLS TO JON FROM PEMAKE-R =',I5)
      RETURN
      END
