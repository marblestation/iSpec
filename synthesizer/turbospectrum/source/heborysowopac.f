*
      subroutine heborysowopac(omega,t,propac)
*
* BPz 6/10/95. interpolates H2-He CIAs from A. Borysow.
* iw is the index of the wavenumber just larger than omega
* it is the index of the temperature just smaller than T
*
      implicit none

      integer j,iw,it,maxwave,maxtemp,i,nwave,ntemp
      parameter (maxwave=4500,maxtemp=7)
      real omega,t,propac,omegalast,qt,qw,h2op(maxtemp),x,exp10
      real wave(maxwave),temp(maxtemp),h2he(maxtemp, maxwave)
      character ident*10
      logical first

      save omegalast,iw,h2op,first,wave,temp,h2he

      data omegalast /0.00/
      data iw        /1/
      data first     /.true./


* ATTENTION EXTRAPOLATION AUX OMEGAS EN DEHORS DE LA TABLE
* DE BORYSOW!!  ET extrapolation a T < 1000K et >7000K!!
*
      if (first) then
        open(33,file=
cc     &  '/b1/plez/BIGGRID/OSFILES/H2Heborysow.opac',
     &    'DATA/H2Heborysow.opac',
     &         status='old')
        read(33,*) ident
        read(33,*) nwave
        read(33,*) ntemp
        if (ntemp.gt.maxtemp) stop 'heborysowopac: maxtemp too small'
        if (nwave.gt.maxwave) stop 'heborysowopac: maxwave too small'
        read(33,*) (temp(i),i=1,ntemp)
        do i=1,nwave
          read(33,*) wave(i),(h2he(j,i),j=1,ntemp)
        enddo
        first=.false.
        close(33)
      endif
 
* no extrapolation to higher wavenumbers!!
      if (omega.gt.wave(1)) then
        propac=0.
        return
      endif

      if (omega.gt.omegalast) then
        iw=1
      endif
      if (omega.ne.omegalast) then
        do while ((omega.le.wave(iw+1)).and.(iw.lt.nwave))
          iw=iw+1
        enddo
        iw=min(iw,nwave-1)
        qw=(omega-wave(iw))/(wave(iw+1)-wave(iw))
        do j=1,ntemp
          h2op(j)=qw*h2he(j,iw+1)+(1.-qw)*h2he(j,iw)
        enddo
      endif
      it=1
      do while (it.lt.ntemp.and.t.ge.temp(it+1))
        it=it+1
      enddo
      it=min(it,ntemp-1)
      qt=(t-temp(it))/(temp(it+1)-temp(it))
      propac=qt*log(h2op(it+1))+(1.-qt)*log(h2op(it))
      propac=exp(propac)
      omegalast=omega

      return
      end
