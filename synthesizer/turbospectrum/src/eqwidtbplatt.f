      SUBROUTINE eqwidtbplatt(lmin,lmax,step,eqwidth)
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
      character*50 mcode
      real  maxetad,etad(ndp),xc(ndp)
      dimension fluxme(lpoint)
      real mum,eqwidth
      doubleprecision XL1,XL2
      doubleprecision DEL,step
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
     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)
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
* Initiate angle quadrature points
*
        NMY=NMX
        CALL GAUSI(NMY,0.,1.,WMY,XMY)
        DO I=1,NMY
          XMY2(I)=XMY(I)*XMY(I)
        enddo
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
      do j=lmin,lmax
        xlsingle=xlambda(j)
*
* continuum calculations:
*
        do k=1,ntau
          x(k)=absocont(k,j)
          s(k)=absoscont(k,j)
          bplan(k)=bpl(t(k),xlsingle)
        enddo
        call traneqplatt(0)
        y1cy1c(j)=y1(nmy)
        fcfc(j)=4.*hsurf
*            
* line calculations
*
        if (iweak.le.0.or.iint.le.0) then
          do k=1,ntau
            x(k)=abso(k,j)
            s(k)=absos(k,j)
          enddo
          call traneqplatt(0)
          prf=4.*hsurf/fcfc(j)
          if(iint.gt.0) then
            prf=y1(nmy)/y1cy1c(j)
          endif
          prof=1.-prf
        else
          maxetad=0.
          do k=1,ntau
            etad(k)=abso(k,j)
            maxetad=max(maxetad,etad(k))
          enddo
          if (maxetad.le.eps) then
            call tranw(ntau,tau,xmyc,bplan,xc,etad,deli)
            prof=deli/y1cy1c(j)
          else
            do k=1,ntau
              x(k)=abso(k,j)
              s(k)=absos(k,j)
            enddo
            call traneqplatt(0)
            prf=4.*hsurf/fcfc(j)
            if(iint.gt.0) then
              prf=y1(nmy)/y1cy1c(j)
              fluxme(j)=y1(nmy)
            endif
            prof=1.-prf
          endif
        endif
        eqwidth=eqwidth+(prof+profold)*step/2.
        profold=prof
*
      enddo
* add small contribution past last point of the profile.
      eqwidth=eqwidth+profold*step/2.
      if (profold.gt.0.0001) then
        print*,' WARNING! last point of calculated profile has depth ',
     &         profold
      endif

ccc      print 1111,eqwidth*1000.
1111    format(' eqwidth: ',f9.1,' mA')
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
     &       ' STEP',F6.3,' A')
  202 FORMAT(' ***WEAK LINE APPROX FOR L/KAPPA L.T.',F7.4)
  203 FORMAT(' MODEL ATMOSPHERE:',A)
  204 FORMAT(' CONTINUUM FLUX=',E12.4, ' AT ',f10.2,'AA')
  205 FORMAT(' CONTINUUM INTENSITY=',E12.4,' at ',f10.2,'AA')
  206 FORMAT(1X,10F8.4)
 2341 FORMAT(1X,'Block number: ',I8)
*
      END
