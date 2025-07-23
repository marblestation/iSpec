      SUBROUTINE eqwidtb(lmin,lmax,step,eqwidth)
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
      parameter (imax=1000)
*
      real pos(imax),intens(imax),mum 
      character*50 mcode
      real eqwidth
      doubleprecision XL1,XL2
      doubleprecision DEL,step
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
      DIMENSION Y1C(NRAYS),Y1L(NRAYS,lpoint),XMUC(NRAYS)
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
     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)
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
* NLBLDU is a dummy
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
        IF(IINT.EQ.0) WRITE(6,201) XL1,XL2,DEL
        IF(IWEAK.GT.0) WRITE(6,202) EPS
        first=.false.
      endif
*
* initialize
      eqwidth=0.e0
      profold=0.e0
*
* Now solve transfer equation in the wavelength range around the line
*
      do j=lmin,lmax
        xlsingle=xlambda(j)
*
* Continuum calculations:
*
        do k=1,ntau
          x(k)=absocont(k,j)
          s(k)=absoscont(k,j)
          bplan(k)=bpl(t(k),xlsingle)
        enddo
        call traneq
*
* Spherical continuum fluxes
*
        hflux1c=4.*pi*hsurf*(rr(1)/radius)**2
        hflux2c=4.*pi*hflux(ntau)*(rr(ntau)/radius)**2
        fcfc(j)=hflux1c/pi
        do mm=1,mmu(1)
          y1cy1c(mm,J)=y1(mm)
          xmucxmuc(mm,J)=xmu(mm,1)
        enddo
        mmucmmuc(J)=mmu(1)
* Line flux
        do k=1,ntau
          x(k)=abso(k,j)
          s(k)=absos(k,j)
        enddo
        call traneq
*
* Spherical line fluxes
*
        hflux1=4.*pi*hsurf*(rr(1)/radius)**2
        hflux2=4.*pi*hflux(ntau)*(rr(ntau)/radius)**2
* normalised flux
        prf=hflux1/pi/fcfc(j)
        if (debug) print*,j,xlsingle,hflux1/pi,prf
        if(iint.gt.0) then
          do mm=1,mmu(1)
            y1l(mm,j)=y1(mm)
          enddo
        endif
*
* End of spherical fluxes
*
        prof=1.-prf
        eqwidth=eqwidth+(prof+profold)*step/2.
        profold=prof
      enddo
* add small contribution past last point of the profile.
      eqwidth=eqwidth+profold*step/2.
      if (profold.gt.0.0001) then
        print*,' WARNING! last point of calculated profile has depth ',
     &         profold
      endif
*
ccc      print 1111,eqwidth*1000.
1111  format(' eqwidth: ',f9.1,' mA'/'warning! flux calculation only!')
      call clock
      return
*
  100 FORMAT(4X,I1,6X,I3)
  207 FORMAT(' SPECTRUM CALCULATED FOR THE FOLLOWING ',I3,' LINES'
     &      /' ELEMENT LAMBDA XI(EV) LOG(GF) LOG(ABUND) LINE NO')
  208 FORMAT('   ',A2,I2,F9.2,F5.2,1X,F6.2,3X,F7.2,4X,I5)
  200 FORMAT(' INTENSITY SPECTRUM AT MU=',F6.2,' BETWEEN',F10.3,
     &      ' AND',F10.3,' A WITH STEP',F6.3,' A')
  201 FORMAT(' FLUX SPECTRUM BETWEEN',F10.3,' A AND ',F10.3,' A WITH',
     &       'max STEP',F6.3,' A (or 0.3 x Doppler width)')
  202 FORMAT(' ***WEAK LINE APPROX FOR L/KAPPA L.T.',F7.4)
  203 FORMAT(' MODEL ATMOSPHERE:',A)
  204 FORMAT(' CONTINUUM FLUX=',E12.4, ' AT ',f10.2,'AA')
  205 FORMAT(' CONTINUUM INTENSITY=',E12.4)
  206 FORMAT(1X,10F8.4)
 2341 FORMAT(1X,'Block number: ',I8)
*
      END
