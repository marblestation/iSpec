
c******************************************************************************
c  The subroutines needed to calculate the He- b-f and f-f opacities are in 
c  this file.  These are from ATLAS9.
c******************************************************************************





      subroutine opacHeminus
c******************************************************************************
c  This routine computes the He- bound-free and free-free opacities.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      save
      data freq1 /0./
         
      if (freq .ne. freq1) then
         freq1 = freq
         a1 =  3.397d-46 + (-5.216d-31+7.039d-15/freq)/freq
         b1 = -4.116d-42 + ( 1.067d-26+8.135d-11/freq)/freq
         cc =  5.081d-37 + (-8.724d-23-5.659d-08/freq)/freq
      endif

      do i=1,ntau
         aHeminus(i) = (a1*t(i) + b1 +c1/t(i))*ne(i)*
     .                 numdens(2,1,i)/u(2,1,i)
      enddo

      return
      end










