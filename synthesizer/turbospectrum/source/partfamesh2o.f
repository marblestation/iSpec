****************************************************************************

      subroutine partfAmesH2O(temp,partf0)
      
      implicit real*8 (a-h,o-z)
c
c     read files containing ro-vibrational levels and computes 
c     partition function
c     unit 1 contains a list of files containing the ro-vibrational
c     levels
c
c     idrv is number of ro-vibrational levels
c     idj is number of j blocks
c     p,q, or r branch.
c
      character*80 fname                                                8d22s96
      character*1 a(80)                                                 8d22s96
      parameter (idrv=180000,idj=60,igrd=8000)
      dimension erov(idrv),iqnrv(5,idrv),iadd(idj,2,2),
     $          pop(idrv),numbb(idj,2,2),
     $          idig(4),ipack(2),idig1(8,10),                           8d22s96
     $          idig2(80)                                               8d22s96
      common/levels/erov,iadd
      equivalence (pack,ipack)                                          8d22s96
      equivalence (idig1,idig2)                                         8d23s96
      dimension spec(igrd,2),popj(idj)
       tocm=1d0/4.556335d-6                                              5d23s91
       xk=3.16681d-6
       conv=tocm*16.1941d5/6.0221367d23
       todebye=2.5417662d0                                               9d16s94
       todebye2=todebye**2
c
c     temp is temperature in kelvin to evaluate partion function
c
      tokm=6.0221367d18
      print*,'temperature = ',temp
      beta=1d0/(xk*temp)
c
c     read in ro-vibrational states and compute partition function
c
      do 1 i=1,idj
       iadd(i,1,1)=0
       iadd(i,2,1)=0
       iadd(i,1,2)=0
       iadd(i,2,2)=0
       numbb(i,1,1)=0                                                   8d22s96
       numbb(i,2,1)=0                                                   8d22s96
       numbb(i,1,2)=0                                                   8d22s96
       numbb(i,2,2)=0                                                   8d22s96
    1 continue
      do 14 i=1,idj
       popj(i)=0d0
   14 continue
      nrv=1
      partf=0d0
      jmax=-1
      jmin=idj                                                          8d22s96
      ehigh=0d0
cc  200 continue                                                          8d22s96
cc      read(1,*,end=201)fname                                            8d22s96
cc      fname='/b1/plez/H2O/OSMAKE/E161'
      fname='DATA/E161'
* this one for 1H1H16O. Change if needed!
*
      print*,' Partition funstion for 1H1H16O from Ames'
      print*,'opening file ',fname
      open(unit=3,file=fname,form='formatted',status='old')             8d22s96
    2 continue
      read(3,3,end=4)jtot,ipar,isym,numb                                8d22s96
    3 format(4i5)
      jmax=max(jmax,jtot)
      jmin=min(jmin,jtot)                                               8d22s96
      if(jtot+1.gt.idj)then
       print*,'too high a jtot. increase idj to ',jtot+1
       stop
      end if
      ipar=mod(ipar,2)+1
      if(isym.eq.0)isym=2
      print7,jtot,ipar,isym,numb
    7 format('jtot = ',i5,' ipar = ',i1,' isym= ',i1,' numb = ',i4)
      iadd(jtot+1,ipar,isym)=nrv
      numbb(jtot+1,ipar,isym)=numb
      if(nrv+numb-1.gt.idrv)then
       print*,'too many ro-vib levels: increase idrv to ',
     $      nrv+numb-1
       stop
      end if
      xj=dfloat(2*jtot+1)
      if(isym.eq.1)xj=xj*3d0
      sumj=0d0
      do 5 i=1,numb
       read(3,6)erov(nrv),(iqnrv(j,nrv),j=1,5)                          8d22s96
       ehigh=max(ehigh,erov(nrv))
       if(partf.eq.0d0)then
        emin=erov(nrv)
       else
        emin=min(emin,erov(nrv))
       end if
       pop(nrv)=xj*exp(-beta*erov(nrv))
       popj(jtot+1)=popj(jtot+1)+pop(nrv)
       sumj=sumj+pop(nrv)
       partf=partf+pop(nrv)
    6  format(1pe21.14,5i4)
       nrv=nrv+1
    5 continue
      sumj=1d0/sumj
      tsumj=0d0
      plast=pop(nrv-1)*sumj
      do 100 i=1,numb
       tsumj=tsumj+pop(nrv-i)
       if(tsumj*sumj.gt.0.01d0)then
        print 103,numb+1-i,numb,plast
  103   format(1x,'99% of pop at root ',i4,' out of ',i4,' roots.',     8d22s96
     $       ' last level has fraction ',1pe10.2)
        go to 2
       end if
  100 continue
      go to 2
    4 continue
      close(unit=3)                                                     8d22s96
cc      go to 200                                                         8d22s96
  201 continue                                                          8d22s96
      nrv=nrv-1
      print*,'total no. ro-vibrational levels = ',nrv,
     $     ' below ',ehigh
      print*,'ground state energy = ',emin
      print*,'partition function = ',partf
      partf0=partf*exp(beta*emin)
      print*,'partition function with energy zero at ground',
     $      ' state energy = ',partf0
      xnorm=1d0/partf
      do 13 i=1,nrv
       pop(i)=pop(i)*xnorm
       erov(i)=erov(i)-emin
   13 continue
      do 15 j=jmin+1,jmax+1                                             8d22s96
       popj(j)=popj(j)*xnorm
       if(mod(j-1,2).eq.0)then                                          8d22s96
        i1a=1                                                           8d22s96
        i1b=2                                                           8d22s96
       else                                                             8d22s96
        i1a=2                                                           8d22s96
        i1b=1                                                           8d22s96
       end if                                                           8d22s96
       print 16,j-1,popj(j),(numbb(j,i1a,i2),numbb(j,i1b,i2),i2=1,2)
   16  format(1x,'fraction in j = ',i3,' is ',f8.6,4i5)
   15 continue

      return
      end
