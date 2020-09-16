
      subroutine partfn (atom,jmark)                                  
c******************************************************************************
c     This routine computes arrays or single partition functions;
c     It uses the data and formulation of Kurucz's ATLAS9 program.
c     However, the partition functions have been updated for several species,
c     and those are calculated in subroutine *newpart*.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Quants.com'


      iat = 10*nint(atom) - 1
      at = dfloat(iat)/10.
      iatom = nint(atom)


c*****compute partition functions for 4 ionization states of an element.
      do k=1,4
         iat = iat + 1
         at = dfloat(iat)/10.
         if (partflag(iatom,k) .gt. 0) then
            do i=1,ntau   
               u(jmark,k,i) = partnew(at,k,i)
            enddo
         else
            do i=1,ntau   
               u(jmark,k,i) = ucalc(at,i)   
            enddo
         endif
      enddo


      return
      end 




