C >>> YCK
C Fall 2000 Written for the standard set of the YY
C Jan 2002  Modified for the alpha enhanced set of the YY
C May 2002  Modified for the interpolation in [a/Fe]
C           modified again for multiple age input
C           modified once again for multiple file outputs
C Dec 2002  modified for the introduction of pseudo-eep
C               numage has increased from 35 to 37
C May 2003      numage has increased from 37 to 41
C <<< YCK
C WARNING     THIS ROUTINE ASSUMES numage=41 AT EACH ISOCHRONE
C             and the set of the 41 ages should be the same.
C numfile = the number of the isochrone files + 1
C
C The names of the two input files are fixed to be 'YY.nml' and 'YY.age'
C
C This program reads data from YY.nml file.
C when the read-in Age=0., the age is read-in from YY.age 
C The name of the source isochrone tables (yy00l.* or yy00g.*)
C should be assigned at YYiso(*).
CCCCCCCCCC   YY.nml   CCCCCCCCCCCCC
C $INPUT
C AFe=0.45       ! [a/Fe] This should be between 0 and 0.6
C targetZ=-0.03  ! if negative, read in the target [Fe/H]
C FeH=-1.5       ! [Fe/H] Be careful not to make Z > 0.08
C Age= 9.        ! age in Gyr: if negative, interactive input
C                !             if 0., read 'YY.age'
C YYout=' '      ! output file name: if ' ', interactive input
C                ! if 'each', the output file name will be asked 
C                !   for each age, and stored in seperate files.
C $END
C $ISET1         ! solar mixture set
C NYYiso=11
C YYiso(1)='a0o2/yy00l.x76997z00001'
C YYiso(2)='a0o2/yy00l.x7697z0001'
C YYiso(3)='a0o2/yy00l.x7688z0004'
C YYiso(4)='a0o2/yy00l.x767z001'
C YYiso(5)='a0o2/yy00l.x758z004'
C YYiso(6)='a0o2/yy00l.x749z007'
C YYiso(7)='a0o2/yy00l.x74z01'
C YYiso(8)='a0o2/yy00l.x71z02'
C YYiso(9)='a0o2/yy00l.x65z04'
C YYiso(10)='a0o2/yy00l.x59z06'
C YYiso(11)='a0o2/yy00l.x53z08'
C $END
C $ISET2         ! twice alpha set
C NYYiso=11
C YYiso(1)='a2o2/yy00l.x76997z00001a2o2'
C YYiso(2)='a2o2/yy00l.x7697z0001a2o2'
C YYiso(3)='a2o2/yy00l.x7688z0004a2o2'
C YYiso(4)='a2o2/yy00l.x767z001a2o2'
C YYiso(5)='a2o2/yy00l.x758z004a2o2'
C YYiso(6)='a2o2/yy00l.x749z007a2o2'
C YYiso(7)='a2o2/yy00l.x74z01a2o2'
C YYiso(8)='a2o2/yy00l.x71z02a2o2'
C YYiso(9)='a2o2/yy00l.x65z04a2o2'
C YYiso(10)='a2o2/yy00l.x59z06a2o2'
C YYiso(11)='a2o2/yy00l.x53z08a2o2'
C $END
C $ISET4         ! four times alpha set
C NYYiso=11
C YYiso(1)='a4o2/yy00l.x76997z00001a4o2'
C YYiso(2)='a4o2/yy00l.x7697z0001a4o2'
C YYiso(3)='a4o2/yy00l.x7688z0004a4o2'
C YYiso(4)='a4o2/yy00l.x767z001a4o2'
C YYiso(5)='a4o2/yy00l.x758z004a4o2'
C YYiso(6)='a4o2/yy00l.x749z007a4o2'
C YYiso(7)='a4o2/yy00l.x74z01a4o2'
C YYiso(8)='a4o2/yy00l.x71z02a4o2'
C YYiso(9)='a4o2/yy00l.x65z04a4o2'
C YYiso(10)='a4o2/yy00l.x59z06a4o2'
C YYiso(11)='a4o2/yy00l.x53z08a4o2'
C $END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program YYinbetween
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: IYY=4
      integer, parameter :: UNIN=7
      integer, parameter :: UNOUT=9
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
C
      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
C
      character*30 YYiso(numfile),YYout,YYout2
      data YYiso/numfile*' '/
      data YYout/' '/
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      data NYYiso,targetZ,FeH,Age/1,-1.0,-1.0,12.0/
      NAMELIST /INPUT/AFe,targetZ,FeH,Age,YYout
      NAMELIST /ISET1/NYYiso,YYiso
      NAMELIST /ISET2/NYYiso,YYiso
      NAMELIST /ISET4/NYYiso,YYiso
C
      open(UNIT=IYY,file='YY.nml',status='old')
      READ(UNIT=IYY, NML=INPUT)
      close(IYY)
         trgtZ0=targetZ
         FeH0=FeH

      nforyi=0
      if(YYout.ne."each")then
      nforyi=1
      YYout2=YYout
      if(YYout2.eq." ")then
      write(6,'(1x,A,$)')' output file name ? '
      read(5,'(a)')YYout2
      endif
      open(UNIT=UNOUT,file=YYout2,status='new')
      endif

      if(Age.gt.0.0d0)then
        agenew=Age
      elseif(Age.eq.0.0d0)then
        open(UNIT=UNIN,file='YY.age',status='old')
        read(UNIN,*,end=900,err=900)agenew
      else
        write(6,'(A,$)')' New Age (input a negative value to stop) '
        read(5,*,end=900,err=900)agenew
        if(agenew.le.0.0e0)go to 900
      endif
      goto 401
  600 continue
      if(Age.lt.0.0d0)then
        write(6,'(A,$)')' New Age (input a negative value to stop) '
        read(5,*,end=900,err=900)agenew
        if(agenew.le.0.0e0)go to 900
      else
        read(UNIN,*,end=900,err=900)agenew
      endif
  401 continue
C
      if(YYout.eq."each")then
      write(6,'(1x,A,$)')' output file name ? '
      read(5,'(a)')YYout2
      open(UNIT=UNOUT,file=YYout2,status='new')
      endif
C
      open(UNIT=IYY,file='YY.nml',status='old')

      if(AFe.eq.0.0)then
         READ(UNIT=IYY, NML=ISET1)
         Nalp=1
         call a_each(targetZ,FeH,YYiso,agenew,NYYiso,Nalp)
         call fin_write(Nalp,nforyi)
       elseif(AFe.eq.0.3)then
         READ(UNIT=IYY, NML=ISET2)
         Nalp=1
         call a_each(targetZ,FeH,YYiso,agenew,NYYiso,Nalp)
         call fin_write(Nalp,nforyi)
       elseif(AFe.eq.0.6)then
         READ(UNIT=IYY, NML=ISET4)
         Nalp=1
         call a_each(targetZ,FeH,YYiso,agenew,NYYiso,Nalp)
         call fin_write(Nalp,nforyi)
       else
         READ(UNIT=IYY, NML=ISET1)
         Nalp=1
         call a_each(targetZ,FeH,YYiso,agenew,NYYiso,Nalp)
CD         call fin_write(Nalp,nforyi)
c
         READ(UNIT=IYY, NML=ISET2)
         Nalp=2
         targetZ=trgtZ0
         FeH=FeH0
         call a_each(targetZ,FeH,YYiso,agenew,NYYiso,Nalp)
CD         call fin_write(Nalp,nforyi)
c
         READ(UNIT=IYY, NML=ISET4)
         Nalp=3
         targetZ=trgtZ0
         FeH=FeH0
         call a_each(targetZ,FeH,YYiso,agenew,NYYiso,Nalp)
CD        call fin_write(Nalp,nforyi)
c
         targetZ=trgtZ0
         FeH=FeH0
         call intpol_alp(FeH,targetZ,AFe)
         Nalp=1
         call fin_write(Nalp,nforyi)
      endif
C
      close(IYY)

      if(YYout.eq."each")then
      close(UNOUT)
      write(6,*)' File created ==> ',YYout2
      endif

      if(Age.gt.0.0d0)goto 900
      go to 600
  900 continue
      if(Age.eq.0.0d0)close(UNIN)
      if(YYout.ne."each")then
      close(UNOUT)
      write(6,*)' File created ==> ',YYout2
      endif
      stop 'All done...'
      end

      subroutine a_each(targetZ,FeH,YYiso,Age,NYYiso,Nalp)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
C
      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
C
      character*30 YYiso(numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass

C Read tables
      do 10 ifile=1,NYYiso
      call yyread(YYiso(ifile), ifile)
   10 continue
      if(NYYISO.le.1)go to 100
C Calculate the [Fe/H] and Z
      call feh2z(FeH,targetZ)
C Z interpolation
      Znew=targetZ
      ifile=NYYiso*0.5
      call findz(YYzarr,NYYiso,Znew,ifile)
C in case
      if( DABS((YYzarr(ifile)-targetZ)/YYzarr(ifile)).lt.0.5d-7)then
       goto 100
      endif
C
      ifile=ifile-1
      if(ifile.ge.(NYYiso-3))ifile=NYYiso-3
      if(ifile.le.1)ifile=1
      write(6,'(A)')'  Z interpolation         V'
      write(6,'(5f10.6)')YYzarr(ifile),YYzarr(ifile+1),Znew,
     +          YYzarr(ifile+2),YYzarr(ifile+3)
      ofile=NYYiso+1
      Ynew=dydz*(Znew-zp)+yp
      YYzarr(ofile)=Znew
      YYyarr(ofile)=Ynew
      call intpol_z(ifile,ofile)
C the output stored in (NYYiso+1)
C YCK      call writeyy(ofile)
      ifile=ofile
C
  100 continue
C Age interpolation
      if(NYYiso.le.1)then
         ifile=1
      endif

        agenew=Age
      call findz(YYage(1,ifile),numage,agenew,iage)
C for linear
      if(iage.ge.(numage-1))iage=numage-1
      if(iage.le.1)iage=1
      write(6,'(A)')'  Age interpolation       V'
      write(6,'(12x,5f10.5)')YYage(iage,ifile),agenew,
     +                       YYage(iage+1,ifile)
C
      call intpol_age(ifile,iage,agenew,Nalp)
  900 continue
      return
      end
C
      subroutine writeyy(ifile,YYisoz)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: IYY=14
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
C
      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
      parameter (FeHa2=-0.217D0)
      parameter (FeHa4=-0.470D0)
C
      character*30 YYisoz
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      real*8 YYage(numage,numfile),YYz(numfile),YYy(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYz,YYy,iYYmass
C
      open(IYY,file=YYisoz,status='new')
      write(IYY,830)YYz(ifile),YYy(ifile),YYos,YYalp,YYfeh,YYalpfe
  830 format(2HZ=,f8.6,3H Y=,f8.6,4H OS=,f4.2,6H l/Hp=,f8.6,8H [Fe/H]=,
     +f9.3,12H [Alpha/Fe]=,f5.2)
      write(IYY,'(A)')header
      do 890 iage=1,numage
C      write(6,*)YYage(iage,ifile),iYYmass(iage,ifile)
      write(IYY,840)YYage(iage,ifile),iYYmass(iage,ifile)
  840 format(9Hage(Gyr)=,f6.3,i4,7H points)
      do 100 imass=1,iYYmass(iage,ifile)
      write(IYY,chfmt)(YYarray(i,imass,iage,ifile),i=1,nump)
  100 continue
      write(IYY,*)' '
  890 continue
      close(IYY)
      write(6,*)' File created -->  ',YYisoz
      return
      end
      subroutine yyread(YYisofl,ifile)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: IYY=24
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
      character*30 YYisofl
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYz(numfile),YYy(numfile)
      integer iYYmass(numage,numfile),ofile
      common /workspace/YYarray,YYage,YYz,YYy,iYYmass
c for nomalization
      integer meep,kipnm
      real*8 eep(numpoint),xeep(numpoint),dteep,tlkip,blkip,slope
      real*8 temar(numparam,numpoint),xxeep
      pollin(x1,y1,x2,y2,x)=
     +  (x-x2)*y1/(x1-x2) +(x-x1)*y2/(x2-x1)
      open(IYY,file=YYisofl,status='old')
      read(IYY,810)YYz(ifile),YYy(ifile),YYos,YYalp,YYfeh,YYalpfe
  810 format(2x,f8.6,3x,f8.6,4x,f4.2,6x,f8.6,8x,f9.6,12x,f5.2)
      read(IYY,'(A)')header
      if ( header(73:75) .eq. 'V-J' )then
        nump=numparam
        chfmt(1)='(f10.7,3f8.4,10f7.3,1p3e12.4)'
      else
        nump=numparam-5
        chfmt(1)='(f10.7,3f8.4,5f7.3,1p3e12.4)'
      endif
      iage=0
 901  continue
      read(IYY,'(9x,f6.3,i4)',end=910)age,nmass
      iage=iage+1
      YYage(iage,ifile)=age
      iYYmass(iage,ifile)=nmass
      do 100 imass=1,nmass
C#       m,t,l,g,mv,ub,bv,vr,vi,vj,vh,vk,vl,vm,n1,n135,n3
      read(IYY,chfmt,err=901)(YYarray(i,imass,iage,ifile),i=1,nump)
C >>> YCK modified for the eep
c--------------yi----------------
       if(imass.le.1)then
        tlkip=YYarray(2,imass,iage,ifile)
        blkip=YYarray(3,imass,iage,ifile)
        eep(imass)=1.0d0
       else
        eep(imass)=eep(imass-1)
     +           + dsqrt(1.0d2*(YYarray(2,imass,iage,ifile)-tlkip)
     +                  *(YYarray(2,imass,iage,ifile)-tlkip)
     +                  +(YYarray(3,imass,iage,ifile)-blkip)
     +                  *(YYarray(3,imass,iage,ifile)-blkip))
        tlkip=YYarray(2,imass,iage,ifile)
        blkip=YYarray(3,imass,iage,ifile)
       endif
c--------------yi----------------
C        eep(imass)=YYarray(1,imass,iage,ifile)
        xeep(imass)=YYarray(2,imass,iage,ifile)
        if(nmass.lt.numpoint)kipnm=nmass
C <<< YCK  Now we have eep(140)
  100 continue
C >>> YCK modified for the eep
         slope=3.d-2
c to find the turnoff mass to utilize it as an anchor point
        if(nmass.ne.kipnm)then
c         call eepset(eep,xeep,nmass,100,slope) 
         call eepset(eep,xeep,nmass,kipnm,slope) 
c         call eepset(eep,xeep,nmass,40,slope) 
        endif

         do im=1,nmass
         xim=dfloat(im-1)*slope
         xeep(im)= datan(xim) *
     +      (eep(nmass)-eep(1))/(datan(slope*dfloat(nmass-1)))
     +      +eep(1)
         enddo
c
       meep=1
       do 2003 im=1,nmass
       xxeep=xeep(im)
       call findz(eep,nmass,xxeep,meep)
       if(meep.ge.(nmass-1))meep=nmass-1
       if(meep.le.1)meep=1
       do  20031 ii=1,nump
       temar(ii,im)= pollin(
     +  eep(meep),YYarray(ii,meep,iage,ifile),
     +  eep(meep+1),YYarray(ii,meep+1,iage,ifile),
     +  xxeep )
20031  continue
 2003  continue
c normalized
       do 20032 im=1,nmass
       do 20033 ii=1,nump
        YYarray(ii,im,iage,ifile)=temar(ii,im)
20033  continue
20032  continue
cc <<< YCK  Done!
      read(IYY,'()')
      goto 901
 910  continue
      close(IYY)
      nage=iage
      return
      end
      subroutine findz(AX,NX,X,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AX(NX)
C FIND THE 'M'
      IF(M.LT.1.OR.M.GT.NX)M=1
      KIP=M
CCC      DO 21 K=1,NX
        IF(X.LT.AX(KIP))THEN
           DO 211 IYC=KIP-1,1,-1
             IF(AX(IYC).LE.X)THEN
               KIP=IYC
               GOTO 213
             ENDIF
  211      CONTINUE
c           WRITE(6,*)' SMALLER THAN THE RANGE : extrapolation '
C           WRITE(ISHORT,*)' ERROR findex: SMALLER THAN THE RANGE '
C           WRITE(6,*) 1,AX(1),X
           KIP=1
        ELSE
           DO 212 IYC=KIP,NX-1
             IF(AX(IYC+1).GT.X)THEN
               KIP=IYC
               GOTO 213
             ENDIF
  212      CONTINUE
               KIP = NX 
c         WRITE(6,*)' LARGER THAN THE RANGE : extrapolation '
C         WRITE(ISHORT,*)' ERROR findex: LARGER THAN THE RANGE '
C         WRITE(6,*)NX,AX(NX),X
        ENDIF
  213   M=KIP
      RETURN
      END
      subroutine intpol_z(ifile,ofile)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
      character*30 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numfile)
      integer iYYmass(numage,numfile),ofile
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      cubeint(x1,y1,x2,y2,x3,y3,x4,y4,x)=
     +  (x-x2)*(x-x3)*(x-x4)*y1/((x1-x2)*(x1-x3)*(x1-x4))
     + +(x-x1)*(x-x3)*(x-x4)*y2/((x2-x1)*(x2-x3)*(x2-x4))
     + +(x-x1)*(x-x2)*(x-x4)*y3/((x3-x1)*(x3-x2)*(x3-x4))
     + +(x-x1)*(x-x2)*(x-x3)*y4/((x4-x1)*(x4-x2)*(x4-x3))

      do 300 iage=1,numage
        iYYm=iYYmass(iage,ifile+1)
       if((YYage(iage,ifile).ne.YYage(iage,ifile+1)) .or.
     +    (YYage(iage,ifile).ne.YYage(iage,ifile+2)) .or.
     +    (YYage(iage,ifile).ne.YYage(iage,ifile+3)))then
        write(6,*)'Warning. Need age interpolation'
       write(6,*)YYage(iage,ifile),YYage(iage,ifile+1),
     +          YYage(iage,ifile+2),YYage(iage,ifile+3)
       stop 'error'
       endif
      YYage(iage,ofile)=YYage(iage,ifile)
      iYYmass(iage,ofile)=iYYm
      do 200 ipoint=1,iYYm
      do 100 iparam=1,nump
      YYarray(iparam,ipoint,iage,ofile)=
     +cubeint(dlog(YYzarr(ifile)),YYarray(iparam,ipoint,iage,ifile),
     +        dlog(YYzarr(ifile+1)),YYarray(iparam,ipoint,iage,ifile+1),
     +        dlog(YYzarr(ifile+2)),YYarray(iparam,ipoint,iage,ifile+2),
     +        dlog(YYzarr(ifile+3)),YYarray(iparam,ipoint,iage,ifile+3),
     +        dlog(YYzarr(ofile)))
  100 continue
  200 continue
  300 continue

      return
      end
      subroutine intpol_age(ifile,iage,agenew,Nalp)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
      real*8 YYnew(numparam)
      character*30 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numfile)
      integer iYYmass(numage,numfile),ofile
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      common /final/AYYfin(8,3),YYfin(numparam,numpoint,3)

      pollin(x1,y1,x2,y2,x)=
     +  (x-x2)*y1/(x1-x2) +(x-x1)*y2/(x2-x1)

      iYYm=iYYmass(iage+1,ifile)
      if(iYYmass(iage,ifile).lt.iYYm)iage=iage+1
c
       AYYfin(1,Nalp)=YYzarr(ifile)
       AYYfin(2,Nalp)=YYyarr(ifile)
       AYYfin(3,Nalp)=YYos
       AYYfin(4,Nalp)=YYalp
       AYYfin(5,Nalp)=YYfeh
       AYYfin(6,Nalp)=YYalpfe
       AYYfin(7,Nalp)=agenew
       AYYfin(8,Nalp)=iYYm
c
      do 200 ipoint=1,iYYm
      do 100 iparam=1,nump
C output
      YYfin(iparam,ipoint,Nalp)=
     +pollin(YYage(iage,ifile),YYarray(iparam,ipoint,iage,ifile),
     +        YYage(iage+1,ifile),YYarray(iparam,ipoint,iage+1,ifile),
     +        agenew)
  100 continue
  200 continue
      return
      end

      subroutine feh2z(FeH,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
      parameter (FeHa2=-0.217D0)
      parameter (FeHa4=-0.470D0)
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt

      quad(x1,y1,x2,y2,x3,y3,x)=y1*(x2-x)*(x3-x)/((x2-x1)*(x3-x1))
     +                         +y2*(x1-x)*(x3-x)/((x1-x2)*(x3-x2))
     +                         +y3*(x1-x)*(x2-x)/((x1-x3)*(x2-x3))

C If the input is in [Fe/H]
      if(Z.lt.0.0d0)then 
       FeH0=FeH
     +      -quad(0.d0,0.d0,0.3d0,FeHa2,0.6d0,FeHa4,YYalpfe)
C       if(abs(YYalpfe-0.3d0).le.0.1d0)FeH0=FeH-FeHa2
C       if(abs(YYalpfe-0.6d0).le.0.1d0)FeH0=FeH-FeHa4
      ZovX=(10.0d0**FeH0)*Zsun/Xsun
      Z=ZovX*(1.0d0+dydz*zp-yp)/(1.0d0+ZovX*(1.0d0+dydz))
C      write(6,*)Z
C
      else
      Znew=Z
C If the input is in Z
      Ynew=dydz*(Znew-zp)+yp
      Xnew=1.0d0 - Ynew - Znew
      FeH=dlog10(Znew*Xsun/(Zsun*Xnew))
     +      +quad(0.d0,0.d0,0.3d0,FeHa2,0.6d0,FeHa4,YYalpfe)
C       if(abs(YYalpfe-0.3d0).le.0.1d0)FeH=FeH+FeHa2
C       if(abs(YYalpfe-0.6d0).le.0.1d0)FeH=FeH+FeHa4
C      write(6,*)FeH
      endif
      write(6,9)YYalpfe,FeH,Z
    9 format(1x,'New mixture ',
     +       ' [a/Fe] = ',f7.3,' [Fe/H] = ',f7.3,' Z =',f9.5)
      YYfeh=FeH
C      write(6,*) YYfeh, ' in feh'
      return
      end
      subroutine fin_write(Nalp,nforyi)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, parameter :: UNOUT=9
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)
      real*8 YYnew(numparam)
      character*30 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numfile)
      integer iYYmass(numage,numfile),ofile
      real*8 mass(numpoint), logt(numpoint),logl(numpoint)
      real*8 massup, massdo,num(3,numpoint)
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      common /final/AYYfin(8,3),YYfin(numparam,numpoint,3)
      data icounter/0/
      save icounter
      icounter=icounter+1

CCCCCCCCCCCCCCCC
      if(nforyi.le.0)then
      write(UNOUT,830)(AYYfin(i,Nalp),i=1,6)
  830 format(3H#Z=,f8.6,3H Y=,f8.6,4H OS=,f4.2,6H l/Hp=,f8.6,8H [Fe/H]=,
     +f9.3,12H [Alpha/Fe]=,f5.2)
      write(UNOUT,'(1h#,A)')header
      iYYm=AYYfin(8,Nalp)
      write(UNOUT,840)AYYfin(7,Nalp),iYYm
  840 format(10H#age(Gyr)=,f6.3,i4,7H points)
      else
C-------
      if(icounter.le.1)then
      write(UNOUT,831)(AYYfin(i,Nalp),i=1,6)
  831 format(2HZ=,f8.6,3H Y=,f8.6,4H OS=,f4.2,6H l/Hp=,f8.6,8H [Fe/H]=,
     +f9.3,12H [Alpha/Fe]=,f5.2)
      write(UNOUT,'(A)')header
      else
      write(UNOUT,'(1h )')
      endif
      iYYm=AYYfin(8,Nalp)
      write(UNOUT,841)AYYfin(7,Nalp),iYYm
  841 format(9Hage(Gyr)=,f6.3,i4,7H points)
C------
      endif
CCCCCCCCCCCCCCCC
c Apr 2003 >>>
c input : 1 mass(iYYm) 2 logt(iYYm) 3 logl(iYYm)
c output : num(3,iYYm)
      do 3000 ipoint=1,iYYm
       mass(ipoint)=YYfin(1,ipoint,Nalp)
       logt(ipoint)=YYfin(2,ipoint,Nalp)
       logl(ipoint)=YYfin(3,ipoint,Nalp)
 3000 continue
      call imf(mass,logt,logl,iYYm,num)
c The output
      do 200 ipoint=1,iYYm
      write(UNOUT,chfmt)(YYfin(ii,ipoint,Nalp),ii=1,nump-3),
     +   (num(iii,ipoint),iii=1,3)
  200 continue
c <<< Apr 2003
      return
      end
      
      subroutine intpol_alp(FeH,targetZ,AFe)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (numfile=21)
      parameter (numage=41)
      parameter (numpoint=140)
      parameter (numparam=17)

      parameter (dydz=2.0E0)
      parameter (yp=0.23E0)
      parameter (zp=0.0E0)
      parameter (Zsun=0.0181D0)
      parameter (Xsun=0.7148705D0)
      real*8 YYnew(numparam)
      character*30 YYiso(numfile)
      real*8 YYarray(numparam,numpoint,numage,numfile)
      real*8 YYage(numage,numfile),YYzarr(numfile),YYyarr(numfile)
      integer iYYmass(numage,numfile),ofile
      character*138 header,chfmt(1)
      common /head/YYos,YYalp,YYfeh,YYalpfe,nump,header,chfmt
      common /workspace/YYarray,YYage,YYzarr,YYyarr,iYYmass
      common /final/AYYfin(8,3),YYfin(numparam,numpoint,3)

      quad(x1,y1,x2,y2,x3,y3,x)=y1*(x2-x)*(x3-x)/((x2-x1)*(x3-x1))
     +                         +y2*(x1-x)*(x3-x)/((x1-x2)*(x3-x2))
     +                         +y3*(x1-x)*(x2-x)/((x1-x3)*(x2-x3))

      write(6,'(A)')'  [a/Fe] interpolation - quadratic '
      write(6,'(A,f10.6)')'   0.00000   0.30000   0.60000 ',AFe

      YYalpfe=AFe
      call feh2z(FeH,targetZ)
       AYYfin(1,1)=targetZ
       AYYfin(2,1)=dydz*(targetZ-zp)+yp
       AYYfin(3,1)=YYos
       AYYfin(4,1)=YYalp
       AYYfin(5,1)=FeH
       AYYfin(6,1)=AFe
      iYYm=AYYfin(8,1)
c
      do 200 ipoint=1,iYYm
      do 100 iparam=1,nump
      YYfin(iparam,ipoint,1)=
     +       quad(0.0d0,YYfin(iparam,ipoint,1),
     +            0.3d0,YYfin(iparam,ipoint,2),
     +            0.6d0,YYfin(iparam,ipoint,3),AFe)
  100 continue
  200 continue
      return
      end
      subroutine eepset(x,y,nstep,igrd,slope)
c      parameter (nstep=140)
      real*8 x(nstep),y(nstep),xp,xpmean,ybar,ymax
      real*8 tomass, totemp,anchor,xhi,xlo,slope,eeptag
c To find the turnoff mass
      ikip=1
      ymax=0.0
      do i=1,nstep
c      write(1,*)x(i),y(i)
      if(y(i).ge.ymax)then
       ikip=i
       ymax=y(i)
      endif
      enddo
c      write(6,*)ikip,x(ikip),ymax
c
      call ToffM(x(ikip-1),y(ikip-1),xp,ybar,1,*99)
      if(xp.lt.x(ikip))call ToffM(x(ikip-2),y(ikip-2),xp,ybar,2,*99)
      if(xp.lt.x(ikip).and.xp.gt.x(ikip+1))then
        write(6,*)' possibly incorrect '
      endif
      go to 88
  99  continue
      xp=x(ikip)
  88  tomass=xp
      totemp=ybar
c      write(6,*)ikip, xp, ybar, (xp-x(1))/(x(nstep)-x(1))
c Got the turnoff mass
c To find the slope
      anchor=(tomass-x(1))/(x(nstep)-x(1))
      xhi=5.0
      xlo=0.005
      slope=0.5*(xhi+xlo)
 2003 continue
      eeptag=atan(slope*float(igrd-1))/atan(slope*float(nstep-1))
      if(abs(anchor-eeptag).le.0.5e-7)goto 2004
      if(eeptag.ge.anchor)then
       xhi=slope
      else
       xlo=slope
      endif
      slope=0.5*(xhi+xlo)
      if(abs((xhi-xlo)/slope).le.0.5d-12)goto 2004
      goto 2003
 2004 continue
      eeptag=atan(slope*float(igrd-1))/atan(slope*float(nstep-1))
c      write(6,*)' slope ',slope,eeptag,anchor
c Got the slope for the eep setup
      return
      end
      subroutine ToffM(x,yy,xp,ybar,iorder,*)
      parameter (n=4)
      real*8 x(n), y(n),yy(n),a,b,c,s,xp,xm,xbar,ybar
      integer k,l,np,m
      do k=1,n
       y(k)=yy(k)
c      write(6,*)x(k),y(k)
      enddo
      do 10 k=1,n-1
      do 20 l=1,(n-k)
        y(l)=(y(l+1)-y(l))/(x(l+k)-x(l))
   20 continue
   10 continue
      a=3.0d0*y(1)
      b=-2.0d0*y(1)*(x(4)+x(3)+x(2))+2.0d0*y(2)
      c=y(1)*(x(3)*x(4)+x(2)*x(4)+x(2)*x(3))-y(2)*(x(4)+x(3))+y(3)
      s=sqrt(b*b-4.0d0*a*c)
c      write(6,*)a,b,c,s
      xp=(-b+s)/(2.d0*a)
      xm=(-b-s)/(2.d0*a)
c      write(6,*)xp,xm
      if(xp.ge.x(iorder).and.xp.le.x(3))then
        xbar=xp
      else if(xm.ge.x(iorder).and.xm.le.x(3))then
       xbar=xm
      else
       if (iorder.ge.2)then
         return 1
       else
         write(6,*)' failed to find the turnoff mass'
         write(6,*)xp,xm
         stop
       endif
      endif
c      xp=(-2.0d0*c) /(b+s)
c      xm=(-2.0d0*c) /(b-s)
      ybar=y(1)
      do 100 m=2,n
       ybar=ybar*(xbar-x(m))+y(m)
  100 continue  
c      write(6,*)ybar
      xp=xbar
      return
      end
      subroutine imf(mass,logt,logl,nmass,num)
C variable definitions.....
      integer nmass
      real*8 alpha(3),ximf(3)
      real*8 mass(nmass), logt(nmass),logl(nmass)
      real*8 massup, massdo,num(3,nmass)

C     alpha: normalization constant for x(IMF)=-1, 1.35(Salpeter), & 3.
C     1000 stars within 0.5-1 Msun range.

      data ximf/1.0d0,-1.35d0,-3.0d0/

C     alpha/(-ximf)

      data alpha/2000.0000d0,-645.5273d0,-142.8571d0/

C compute number density
C     massdo: M_down_cut   massup: M_up_cut

         do 150 j=1,nmass

            if (j.eq.nmass) then 
               massup=mass(j)
            else
               massup=mass(j)+0.5d0*(mass(j+1)-mass(j))
            endif

            if (j.eq.1) then 
               massdo=mass(j)
            else
               massdo=mass(j)-0.5d0*(mass(j)-mass(j-1))
            endif

            do 120 i=1,3
               num(i,j)=alpha(i)*(massup**ximf(i)-massdo**ximf(i))
 120        continue
 150        continue

      return
      end
