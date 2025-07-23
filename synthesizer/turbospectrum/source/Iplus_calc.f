      subroutine Iplus_calc(ntau,tau,source,Iplus,optthin)
*
* Calculates Iplus (outwards directed intensities) along
* a ray at all depths in a model atmosphere, by integration of 
* the source function. Input: tau is the optical depth along
* the ray, source is the source function ,
* both provided at ntau depth points.  BPz 19/07-2001
* Revised by BPz 27/06-2012:
* Revised by BPz on 09/04-2020 following test comparison with analytical solutions
* (S=polynomial(tau))
*
* Numbering of layers goes from 1 at the observer inwards. 
* The tau-scale increases away from the observer.
*
* Opthin means that the ray goes throughthe atmosphere.
* THIS VERSION ASSUMES SPHERICAL SYMMETRY
*
      implicit none
      logical optthin
      integer i,ntau
      real Iplus(*),source(*),tau(*),ex(ntau+1),mex,dtau

      do i=2,ntau
!          ex(i)=exp(tau(i-1)-tau(i))    ! dtau is >0, we need exp(-dtau)
        dtau=tau(i)-tau(i-1)
        if (dtau.lt.0.01) then
          mex=dtau*(1.-dtau/2.*(1.-dtau/3.*(1.-dtau/4.)))
          ex(i)=1.-mex
        else
          ex(i)=exp(-dtau)
        endif
      enddo
      if (optthin) then
* optically thin: assumes no incoming radiation above taumax
* BPlez 26/03-2020 replaces this by calcultaing intensity generated between tau=0 and first depth point.
!        Iplus(ntau)=0.
!        Iplus(ntau)=source(ntau)*(1.-exp(-tau(ntau)))
        mex=tau(1)*(1.-tau(1)/2.*(1.-tau(1)/3.*(1.-tau(1)/4.)))
        if (tau(1).gt.0.01) mex=1.-exp(-tau(1))
        Iplus(ntau)=source(ntau)*mex             ! CHECK THIS ! should compute a tau above the last point in the back
      else
* optically thick: diffusion approximation
        Iplus(ntau)=source(ntau)+
     &     (source(ntau)-source(ntau-1))/(tau(ntau)-tau(ntau-1))
      endif
      do i=ntau-1,1,-1
!        Iplus(i)=Iplus(i+1)*ex(i+1)+
!     &           (source(i+1)+source(i))*0.5*
!     &                  (1.-ex(i+1))
! adopt Trapezoidal rule within each dtau interval (i.e. assumes Sexp(-tau) is affine function
         Iplus(i)=Iplus(i+1)*ex(i+1)+
     &            (source(i+1)*ex(i+1)+source(i))*0.5*(tau(i+1)-tau(i))
      enddo

cc! Iplus(1) corrected for emission/absorption from layer between 
cc! tau=tau(1) and tau=0.

cc      Iplus(1)=Iplus(1)*exp(-tau(1))+
cc     &         source(1)*(1.-exp(-tau(1)))

! BPlez: 25/03-2020. Iplus(1) is intensity at tau(1). 
! It must be corrected to tau=0 to get intensity at observer.

      return
      end

