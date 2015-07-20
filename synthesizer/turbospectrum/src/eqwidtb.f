      SUBROUTINE eqwidtb(lmin,lmax,eqwidth)
*
*-----------------------------------------------------------------------
*
* eqwidtb computes equivalent width for eqwidt (in spherical geometry)
*      06/07-2000 BPz
*
*-----------------------------------------------------------------------
*
      INCLUDE 'spectrum.inc'
*
      INTEGER    MAXNL
      PARAMETER (MAXNL=200)
      INTEGER    NLBL
      PARAMETER (NLBL=100)
      parameter (imax=1000)
*
      real pos(imax),intens(imax),mum 
      character*50 mcode
      CHARACTER*20  LELE(100)
      DIMENSION ION(100),SRXLAL(2*MAXNL),
     &          GFELOG(100),ABUND(100),ETAK(2*MAXNL),ETAR(NDP*2*MAXNL),
     &          ETAD(NDP),XC(NDP),CHIE(100),UTREC(103)
      real eqwidth
      doubleprecision XL1,XL2
      doubleprecision XLBEG,XLEND,DEL,XLB(100),XL(100),XLAL(2*MAXNL)
      COMMON/ATMOS/ T(NDP),PE(NDP),PG(NDP),XI(NDP),MUM(NDP),RO(NDP),
     & nnNTAU
*
* Special for spherical
*
      COMMON /CTRAN/X(NDP),S(NDP),BPLAN(NDP),XJ(NDP),HFLUX(NDP),XK(NDP)
     & ,FJ(NDP),SOURCE(NDP),TAUS(NDP),DTAUS(NDP),JTAU0,JTAU1,ISCAT
      COMMON /TAUC/TAU(NDP),DTAULN(NDP),NTAU   
      COMMON /ROSSC/ROSS(NDP),CDUMM(NDP) 
      COMMON /RHOC/RHO(NDP)
      COMMON /CSPHER/NCORE,DIFLOG,RADIUS,RR(NDP)
      COMMON /CSTYR/MIHAL  /CTAUM/TAUM
      COMMON /SPACE2/ SPACEDUM1(NDP*7+NDP*NRAYS*5+NRAYS*3+1),
     &                MMU(NDP),SPACEDUM2(NDP*2),PFEAU(NRAYS,NDP),
     &                XMU(NRAYS,NDP)
      COMMON /TRDBUG/IDEBUG
      common /limbdk/ pos,intens,totintens,tottrans
      LOGICAL DEBUG
      DIMENSION FLUXME(nlbl)
      DIMENSION Y1C(NRAYS),Y1L(NRAYS,100),XMUC(NRAYS)
      dimension fcfc(lpoint),y1cy1c(nrays,lpoint),
     & xmucxmuc(nrays,lpoint),mmucmmuc(lpoint),xlm(lpoint),
     & jlcont(lpoint)
*
      COMMON/CSURF/ HSURF,Y1(NRAYS)
CCC      COMMON/CANGLE/ NMY,XMY(6),XMY2(6),WMY(6)
      COMMON/PIECES/ XL1,XL2,DEL,EPS,NMX,NLBLDU,IINT,XMYC,IWEAK
*
* extension for large number of wavelengths and lines (monster II)
      doubleprecision xlambda
      common/large/ xlambda(lpoint),maxlam,ABSO(NDP,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint)
*
      character*256 filterfil
      character*80 filttitle
      logical limbdark,first,findtau1
      common/filter/limbdark,ifilt,filtlam(1000),filttrans(1000),
     &              filterfil,filttitle
      logical hydrovelo
      real velocity
      common/velo/velocity(ndp),hydrovelo
*
*
* NLBLDU is a dummy, real val. of NLBL is set as parameter
*
      DATA first/.true./
      data debug/.false./
      PI=3.141593
*
* Initiate mode of calculation
* IINT =1  : intensity at all mu-points
* IINT =0  : flux
* XL1      : wavelength where synthetic spectrum starts
* XL2      : wavelength where synthetic spectrum stops
* DEL      : wavelength step
* IWEAK =1 : weak line approximation for L/KAPPA le. EPS
* note  XL2.gt.XL1
*
      if (first) then
        iweak=0
        IF(IINT.EQ.0) WRITE(7,201) XL1,XL2,DEL
        IF(IWEAK.GT.0) WRITE(7,202) EPS
        first=.false.
      endif
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
        if (debug) then
          print*,'eqwidtb: calling traneq for continuum',jc,xlm(jc)
          do k=1,ntau,5
            print*,(x(kk),s(kk),kk=k,k+4)
          enddo
        endif
        CALL TRANEQ
*
* Spherical fluxes
*
C
C FLUX TO PRINT
        HFLUX1C=4.*PI*HSURF*(RR(1)/RADIUS)**2
ccc        print*,hflux1c
        HFLUX2C=4.*PI*HFLUX(NTAU)*(RR(NTAU)/RADIUS)**2
        FCFC(JC)=HFLUX1C/PI
C        GFLUX1C=HFLUX1/PI*XL(J)**2/CLIGHT
C        GFLUX2C=HFLUX2/PI*XL(J)**2/CLIGHT
C        FFLUX1C=-2.5*ALOG10(AMAX1(GFLUX1,1E-20))
C        FFLUX2C=-2.5*ALOG10(AMAX1(GFLUX2,1E-20))
        DO 900 MM=1,MMU(1)
          Y1CY1C(MM,jc)=Y1(MM)
          XMUCXMUC(MM,jc)=XMU(MM,1)
 900    CONTINUE
        MMUCMMUC(jc)=MMU(1)
C
*
* End of spherical fluxes
*
        if (debug) print*,jc,xlm(jc),fcfc(jc)
        IF(IINT.LE.0) WRITE(7,204) fcFC(jc),xlm(jc)
CC      IF(IINT.GT.0) WRITE(7,205) Y1C
1963  continue
*
* initialize
      eqwidth=0.e0
      profold=0.e0
*
* Now solve transfer equation in the wavelength range around the line
*
      DO 39 j=lmin,lmax
        xlsingle=xlambda(j)
        DO 30 K=1,NTAU
* the continuum opacity is already included in abso
ccc          X(K)=XC(K)+ABSO(K,J)
          X(K)=ABSO(K,J)
          s(k)=absos(k,j)
          BPLAN(k)=BPL(T(k),xlsingle)
   30   CONTINUE
        CALL TRANEQ
*
* Spherical fluxes
*
C
C FLUX TO PRINT
        HFLUX1=4.*PI*HSURF*(RR(1)/RADIUS)**2
        HFLUX2=4.*PI*HFLUX(NTAU)*(RR(NTAU)/RADIUS)**2
ccc        FLUXME(J)=HFLUX1/PI
C        GFLUX1=HFLUX1/PI*XL(J)**2/CLIGHT
C        GFLUX2=HFLUX2/PI*XL(J)**2/CLIGHT
C        FFLUX1=-2.5*ALOG10(AMAX1(GFLUX1,1E-20))
C        FFLUX2=-2.5*ALOG10(AMAX1(GFLUX2,1E-20))
C
* interpolate continuum flux
        do 3691 jc=2,nlcont
          if(xlm(jc)-xlsingle.gt.0.) then
            jjc=jc-1
            goto 3692
          endif
3691    continue
3692    continue
        if (nlcont.gt.1) then
          jjc=min(jjc,nlcont-1)
          fc=(fcfc(jjc+1)-fcfc(jjc))/
     &       (xlm(jjc+1)-xlm(jjc))*(xlsingle-xlm(jjc)) + fcfc(jjc)
        else
          fc=fcfc(1)
        endif
        PRF=hflux1/pi/FC
        if (debug) print*,j,xlsingle,hflux1/pi,prf
        IF(IINT.GT.0) THEN
          DO 901 MM=1,MMU(1)
            Y1L(MM,J)=Y1(MM)
 901      CONTINUE
        ENDIF
*
* End of spherical fluxes
*
        PROF=1.-PRF
        eqwidth=eqwidth+(prof+profold)*del/2.
        profold=prof
   39 CONTINUE
*
ccc      print 1111,eqwidth*1000.
1111  format(' eqwidth: ',f9.1,' mA'/'warning! flux calculation only!')
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
  205 FORMAT(' CONTINUUM INTENSITY=',E12.4)
  206 FORMAT(1X,10F8.4)
 2341 FORMAT(1X,'Block number: ',I8)
*
      END
