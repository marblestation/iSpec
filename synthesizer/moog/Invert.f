
      subroutine invert(order,a,i,nsiz)                                 
c******************************************************************************
c     This routine inverts an *order,order* matrix *a*. answer is *i*       
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      integer order, orde
      real*8 a(nsiz,nsiz), i(nsiz,nsiz)


      do n=1,order                                                   
         do m=1,order                                                   
            i(m,n) = 0.0
         enddo
      enddo
      do j=1,order                                                   
         i(j,j) = 1.0
      enddo
      orde=order-1                                                      
      do k=1,orde
         kk = k+1
         do j=kk,order
            i(j,k) = 0.0
            i(k,j) = 0.0
         enddo
      enddo

c  here begins the inversion                                             
      do n=1,order                                                   
         if (n .eq. order) go to 244
         g = a(n,n)
         n1 = n+1
         ng = n
         do m=n1,order                                                  
            f = a(m,n)
            if (f .gt. g) then
               g = f
               ng = m
            endif
         enddo
         if (ng .eq. n) go to 244
         do k=1,order
            d = i(n,k)
            e = i(ng,k)
            f = a(n,k)
            g = a(ng,k)
            i(ng,k) = d
            i(n,k) = e 
            a(ng,k) = f
            a(n,k) = g
         enddo
244      g = a(n,n)
         if (g .eq. 0.0) then
            write (nf1out,1001)     
            return
         endif
         do k=1,order
            a(n,k) = a(n,k)/g
            i(n,k) = i(n,k)/g
         enddo
         if (n .eq. order) go to 27
         do k=n1,order
            f = -a(k,n)
            do j=1,order
               a(k,j) = a(k,j) + f*a(n,j)
               i(k,j) = i(k,j)+f*i(n,j)
            enddo
         enddo
      enddo

27    do jk=1,orde                                                   
         n = order - jk + 1
         do j=jk,orde
            m = order - j
            f = -a(m,n)
            do k=1,order
            a(m,k) = a(m,k) + f*a(n,k)
            i(m,k) = i(m,k) + f*i(n,k)
            enddo
         enddo
      enddo
      return


c*****format statements
1001  format('WARNING: AN UN-INVERTABLE ARRAY HAS BEEN ENCOUNTERED!')

      end                                                               


