      SUBROUTINE eqwidtbplatt(lmin,lmax,eqwidth)
*
*-----------------------------------------------------------------------
*
* eqwidtbplatt computes equivalent width for eqwidt, in PP geometry.
*  06/07-2000  BPz
*
*-----------------------------------------------------------------------
*
      include 'spectrum.inc'
*
      INTEGER    MAXNL
      PARAMETER (MAXNL=200)
      INTEGER    NLBL
      PARAMETER (NLBL=100)
*
      CHARACTER*20  LELE(100)
      character*50 mcode
      DIMENSION ION(100),SRXLAL(2*MAXNL),
     &          GFELOG(100),ABUND(100),ETAK(2*MAXNL),ETAR(ndp*2*MAXNL),
     &          ETAD(ndp),XC(ndp),CHIE(100),UTREC(103)
      dimension fluxme(nlbl)
      real mum,eqwidth
      doubleprecision XL1,XL2
      doubleprecision XLBEG,XLEND,DEL,XLB(100),XL(100),XLAL(2*MAXNL)
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),
     &  nnNTAU
*
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON/CSURF/ HSURF,Y1(NRAYS)
      COMMON/CANGLE/ NMY,XMY(nrays),XMY2(nrays),WMY(nrays)
      COMMON/TAUC/ TAU(ndp),DTAULN(ndp),NTAU
      COMMON/PIECES/ XL1,XL2,DEL,EPS,NMX,NLBLDU,IINT,XMYC,IWEAK
*
* extension for large number of wavelengths and lines (monster II)
      doubleprecision xlambda
      common/large/ xlambda(lpoint),maxlam,ABSO(NDP,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint)
*
      dimension fcfc(lpoint),y1cy1c(lpoint),
     & xlm(lpoint),jlcont(lpoint)

      logical findtau1,hydrovelo,first
      real velocity
      common/velo/velocity(ndp),hydrovelo
*
* NLBLDU is a dummy, real val. of NLBL is set as parameter
*
      DATA first/.true./
*
      if (first) then
* Initiate angle quadrature points and number of wavelengths
* per block (NLBL)
*
        NMY=NMX
        CALL GAUSI(NMY,0.,1.,WMY,XMY)
        DO 1 I=1,NMY
          XMY2(I)=XMY(I)*XMY(I)
    1   CONTINUE
*
* Initiate mode of calculation
* IINT =1  : intensity at MY=XMYC
* IINT =0  : flux
* XL1      : wavelength where syntetic spectrum starts
* XL2      : wavelength where syntetic spectrum stops
* DEL      : wavelength step
* IWEAK =1 : weak line approximation for L/KAPPA le. EPS
* note  XL2.gt.XL1
*
        iweak=0
        IF(IINT.GT.0) THEN
          NMY=NMY+1
          XMY(NMY)=XMYC
          WRITE(7,200) XMYC,XL1,XL2,DEL
        END IF
        IF(IINT.EQ.0) WRITE(7,201) XL1,XL2,DEL
        IF(IWEAK.GT.0) WRITE(7,202) EPS
        first=.false.
      endif
      profold=0.e0
      eqwidth=0.e0
*
* Read model atmosphere
*
      REWIND 14
      READ(14) MCODE,nlcont,xlm(1),BPLAN,XC,S,XI
      REWIND 14
      WRITE(7,203) MCODE(1:lenstr(mcode))
*
* Continuum calculations:
*
      do 1963 jc=1,nlcont
        READ(14) MCODE,idum,xlm(jc),BPLAN,XC,S,XI
        DO 9 K=1,NTAU
          X(K)=XC(K)
    9   CONTINUE
        CALL traneqplatt(0)
        Y1CY1C(jc)=Y1(NMY)
        FCFC(jc)=4.*HSURF
        IF(IINT.LE.0) WRITE(7,204) fcFC(jc),xlm(jc)
        IF(IINT.GT.0) WRITE(7,205) Y1Cy1c(jc),xlm(jc)
1963  continue
*
*
        DO 39 j=lmin,lmax
          xlsingle=xlambda(j)
          IF(IWEAK.LE.0.OR.IINT.LE.0) THEN
            DO 30 K=1,NTAU
* the continuum opacity is already included in abso
ccc              X(K)=XC(K)+ABSO(K,J)
              X(K)=ABSO(K,J)
              s(k)=absos(k,j)
              BPLAN(k)=BPL(T(k),xlsingle)
   30       CONTINUE
            CALL traneqplatt(0)
* interpolate continuum flux
            do 3691 jc=2,nlcont
              if(xlm(jc)-xlsingle.gt.0.) then
                jjc=jc-1
                goto 3692
              endif
3691        continue
3692        continue
	    if (nlcont.gt.1) then
              jjc=min(jjc,nlcont-1)
              fc=(fcfc(jjc+1)-fcfc(jjc))/
     &           (xlm(jjc+1)-xlm(jjc))*(xlsingle-xlm(jjc)) + fcfc(jjc)
              y1c=(y1cy1c(jjc+1)-y1cy1c(jjc))/
     &           (xlm(jjc+1)-xlm(jjc))*(xlsingle-xlm(jjc)) + y1cy1c(jjc)
	    else
              fc=fcfc(1)
              y1c=y1cy1c(1)
            endif
            PRF=4.*HSURF/FC
            IF(IINT.GT.0) then
              PRF=Y1(NMY)/Y1C
            endif
            PROF=1.-PRF
          ELSE
            DO 31 K=1,NTAU
              ETAD(K)=ABSO(K,J)
              IF(ETAD(K).GT.EPS) GOTO 32
   31       CONTINUE
            CALL TRANW(NTAU,TAU,XMYC,BPLAN,XC,ETAD,DELI)
            PROF=DELI/Y1C
            GOTO 39
   32       DO 33 K=1,NTAU
              X(K)=ABSO(K,J)
              s(k)=absos(k,j)
              BPLAN(k)=BPL(T(k),xlsingle)
   33       CONTINUE
            CALL traneqplatt(0)
* interpolate continuum flux
            do 3693 jc=2,nlcont
              if(xlm(jc)-xlsingle.gt.0.) then
                jjc=jc-1
                goto 3694
              endif
3693        continue
3694        continue
	    if (nlcont.gt.1) then
              jjc=min(jjc,nlcont-1)
              fc=(fcfc(jjc+1)-fcfc(jjc))/
     &           (xlm(jjc+1)-xlm(jjc))*(xlsingle-xlm(jjc)) + fcfc(jjc)
              y1c=(y1cy1c(jjc+1)-y1cy1c(jjc))/
     &           (xlm(jjc+1)-xlm(jjc))*(xlsingle-xlm(jjc)) + y1cy1c(jjc)
            else
              fc=fcfc(1)
              y1c=y1cy1c(1)
	    endif
            PRF=4.*HSURF/FC
cc            fluxme(j)=hsurf*4.
            IF(IINT.GT.0) then
              PRF=Y1(NMY)/Y1C
              fluxme(j)=y1(nmy)
            endif
            PROF=1.-PRF
          END IF
          eqwidth=eqwidth+(prof+profold)*del/2.
          profold=prof
*
39      CONTINUE
ccc      print 1111,eqwidth*1000.
1111    format(' eqwidth: ',f9.1,' mA')
        call clock
      RETURN
*
  100 FORMAT(4X,I1,6X,I3)
  207 FORMAT(' SPECTRUM CALCULATED FOR THE FOLLOWING ',I3,' LINES'
     &      /' ELEMENT LAMBDA XI(EV) LOG(GF) LOG(ABUND) LINE NO')
  208 FORMAT('   ',A2,I2,F9.2,F5.2,1X,F6.2,3X,F7.2,4X,I5)
  200 FORMAT(' INTENSITY SPECTRUM AT MU=',F6.2,' BETWEEN',F10.3,
     &      ' AND',F10.3,' A WITH STEP',F6.3,' A')
  201 FORMAT(' FLUX SPECTRUM BETWEEN',F10.3,' A AND ',F10.3,' A WITH',
     &       ' STEP',F6.3,' A')
  202 FORMAT(' ***WEAK LINE APPROX FOR L/KAPPA L.T.',F7.4)
  203 FORMAT(' MODEL ATMOSPHERE:',A)
  204 FORMAT(' CONTINUUM FLUX=',E12.4, ' AT ',f10.2,'AA')
  205 FORMAT(' CONTINUUM INTENSITY=',E12.4,' at ',f10.2,'AA')
  206 FORMAT(1X,10F8.4)
 2341 FORMAT(1X,'Block number: ',I8)
*
      END
