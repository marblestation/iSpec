      subroutine weightlimb(xl,wl)
*
* must interpolate at fixed xmu before summing contributions!
*  note xmu<0 by definition in routine tranfr
*
        include 'spectrum.inc'
*
        parameter (imax=1000)
        character filterfil*256,filttitle*80
        logical first,limbdark
        real intens(imax),x(nrays),y(nrays),yinterp(imax),pos(imax),
     &            postheta(imax),xtheta(nrays)
        common /limbdk/ pos,intens,totintens,tottrans
        common /csurf/ hsurf,y1(nrays)
        common /space2/ spacedum1(ndp*7+ndp*nrays*5+nrays*3+1),
     &                  mmu(ndp),spacedum2(ndp*2),pfeau(nrays,ndp),
     &                  xmu(nrays,ndp)
        common /filter/ limbdark,ifilt,filtlam(1000),filttrans(1000),
     &                  filterfil,filttitle
        save first,postheta,lamcount
        data first/.true./
*
        if (first) then
          lamcount=0
          totintens=0.
          tottrans=0.
          do i=1,imax
            intens(i)=0.
            pos(i)=float(imax+1-i)/float(imax)
* rem: pos equidistant, for easy use in FFT or HT. outermost point
* (i=1) is outermost in x also. (not r=1)
* postheta(i) est le cos theta correspondant a pos(i)
            postheta(i)=sqrt(1.-(float(imax+1-i)/float(imax))**2)
          enddo
          first=.false.
        endif
*
        n=mmu(1)
        do i=1,n
          x(i)=-xmu(n+1-i,1)
          y(i)=y1(n+1-i)
        enddo
* reordered increasing. ycenter = Imu at cos(theta)=1
        ycenter=y(n)
        do i=1,n
          y(i)=y(i)/ycenter
        enddo
* x est maintenant en cos(theta).
* postheta aussi. (tous les deux ordonnes croissant)
* now, interpolate on pos(i)
	do i=1,imax
          call tint (n,x,y,postheta(i),yinterp(i))
        enddo
* interpolate filter transmission
        call tint(ifilt,filtlam,filttrans,xl,trans)
        totintens=totintens+wl*trans*ycenter
        tottrans=tottrans+wl*trans
* then add up contributions 
        do i=1,imax
          intens(i)=intens(i)+wl*ycenter*yinterp(i)*trans
        enddo
        lamcount=lamcount+1

        return
        end
