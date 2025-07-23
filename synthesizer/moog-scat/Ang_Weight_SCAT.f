      subroutine ang_weight_scat
c****************************************************************************************************************
c    Performs a Gaussian Quadrature Summation via the calculation of the Gaussian integration points and weights.
c    Currently operates ONLY with odd numbers of rays (mmu = 1, 3, 5,...).  Possible to correct this.
c    Note that with increasing ray number (>= 5), the mu value more closely approaches 1.  The number of rays (mmu)
c    is hard-coded (but this routine allows for mmu to be changed by the user).
c****************************************************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Scat.com'

      mmu = 5.
      aaa = 0.
      bbb = 1.
      call gauss_quad_sum (mmu,aaa,bbb,wtmu,mu)
      xxx = 0
      do i=1,mmu
         xxx = xxx + wtmu(i)
      enddo     
      write(*,*) 'ang_weight done! wtmu ', wtmu

      return
      end
