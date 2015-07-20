      subroutine partfBarberH2O(temp,Q)
*
* calculate partition function for H2O Barber.
*  BPz 10/04-2007
*
*

      implicit none
      integer i,j,ns,nr,n,nlevel
      doubleprecision  kthc,q,temp,energy(222000)
      doubleprecision g(222000),ge,go,khc
      logical first
      data first/.true./


      if (first) then
* ge and go are nuclear spin factors (even and odd respectively)
        ge=1.0d0
        go=3.0d0

        khc=1.380658e-16/6.6260755e-27/2.997925e10

        open(3,file='DATA/energy-levels-BarberH2O.txt')
        i=1
        do while (.true.) 
          read(3,*,end=99) n,j,ns,nr,energy(i)
          i=i+1
          if (ns.le.2) then
            g(i)=ge*dfloat(2*j+1)
          else
            g(i)=go*dfloat(2*j+1)
          endif
        enddo
99      continue
        nlevel=i-1
        close(3)
        first=.false.
      endif
*  kT/hc
      kthc=temp*khc
      Q=0.d0

      do i=1,nlevel
        Q=Q+g(i)*exp(-energy(i)/kthc)
      enddo
      return

      end
