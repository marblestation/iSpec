
      subroutine setmols
c******************************************************************************
c     This subroutine is invoked at the end of each call to *eqlib*:
c     it transfers some of the output number densities to arrays needed
c     for opacities, etc.
c****************************************************************************** 

      implicit real*8 (a-h,o-z) 
      include 'Atmos.com'
      include 'Mol.com'
      include 'Dummy.com'
      integer nel(7)
      integer append
      data (nel(i),i=1,7)/  1,  2,  6, 12, 13, 14, 26/


c     In the number density array "numdens", the elements denoted by
c     the first subscripts are named in at the ends of the assignment
c     lines; at present these are the only ones needed for continuous 
c     opacities
c     the desired elements are H, He, C, Mg, Al, Si, Fe, along
c     with the H_2 molecule
c     first do the neutrals
      do j=1,7
         do k=1,neq
            if (nel(j) .eq. iorder(k)) then
               do i=1,ntau
                  numdens(j,1,i) = xamol(k,i)
               enddo
               exit
            endif
         enddo
      enddo


c*****then do the ions
      do j=1,7
         ispec10 = nint(dble(10*nel(j)+1))
         do k=1,nmol
            kmol10 = nint(10.*amol(k))
            if (ispec10 .eq. kmol10) then
               do i=1,ntau
                  numdens(j,2,i) = xmol(k,i)
               enddo
               exit
            endif
         enddo
      enddo


c*****finally add in H_2
      do k=1,nmol
         ispec = 101
         if (ispec .eq. nint(amol(k))) then
            do i=1,ntau
               numdens(8,1,i) = xmol(k,i)
            enddo
            exit
         endif
      enddo


c*****compute partitiion functions for H_2O and CO_2; 
      do i=1,ntau
         h2olog = 0.
         co2log = 0.
         if (t(i) .gt. 5000.) then
            uh2o(i) = 1.d8
            uco2(i) = 1.d8
         else
            do j=1,5
               h2olog = h2olog + h2ocoeff(j)*t(i)**(j-1)
               co2log = co2log + co2coeff(j)*t(i)**(j-1)
            enddo
            uh2o(i) = 10**h2olog
            uco2(i) = 10**co2log
         endif
      enddo


c*****transfer H_2O and CO_2 number densities from the molecular 
c     equilibrium output
c     note:  HITRAN partition functions are given only for T < 5000K
      do j=1,nmol
         if (nint(amol(j)) .eq. 10108) ih2o = j
         if (nint(amol(j)) .eq. 60808) ico2 = j
      enddo
      do i=1,ntau
         xnh2o(i) = xmol(ih2o,i)
         xnco2(i) = xmol(ico2,i)
      enddo
         
      
c*****yes, perhaps not the most efficient coding, but it works - I hope
      return

      end

