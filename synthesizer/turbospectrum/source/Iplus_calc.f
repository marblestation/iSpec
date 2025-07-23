      subroutine Iplus_calc(ntau,tau,source,Iplus,optthin)
*
* Calculates Iplus (outwards directed intensities) along
* a ray at all depths in a model atmosphere, by integration of 
* the source function. Input: tau is the optical depth along
* the ray, source is the source function ,
* both provided at ntau depth points.  BPz 19/07-2001
* Revised by BPz 27/06-2012:
* Numbering of layers goes from 1 at the observer inwards. 
* The tau-scale increases away from the observer.
*
* Opthin means that the ray goes throughthe atmosphere.
* THIS VERSION ASSUMES SPHERICAL SYMMETRY
*
      implicit none
      logical optthin
      integer i,ntau
      real Iplus(*),source(*),tau(*),ex(5000)

      if (ntau.gt.5000) stop 'Iplus_calc.f : increase ex() dimension!'
      do i=2,ntau
        ex(i)=exp(tau(i-1)-tau(i))
      enddo
      if (optthin) then
* optically thin: assumes no incoming radiation above taumax

        Iplus(ntau)=0.
      else
* optically thick: diffusion approximation
        Iplus(ntau)=source(ntau)+
     &     (source(ntau)-source(ntau-1))/(tau(ntau)-tau(ntau-1))
      endif
      do i=ntau-1,1,-1
        Iplus(i)=Iplus(i+1)*ex(i+1)+
     &           (source(i+1)+source(i))*0.5*
     &                  (1.-ex(i+1))
      enddo

! Iplus(1) corrected for emission/absorption from layer between 
! tau=tau(1) and tau=0.

      Iplus(1)=Iplus(1)*exp(-tau(1))+
     &         source(1)*(1.-exp(-tau(1)))

      return
      end

