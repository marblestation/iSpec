
      subroutine taukap
c******************************************************************************
c     This routine calculates the line absorption coefficient and the line  
c     opacity at wavelength *wave* for all lines in the spectrum            
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dummy.com'


c*****compute the total line opacity at each depth                          
      do i=1,ntau     
         kapnu(i) = 0.0
         do j=lim1,lim2
            v = 2.997929d10*dabs(wave-wave1(j))/
     .             (wave1(j)*dopp(j,i))            
            kapnu(i) = kapnu(i) + kapnu0(j,i)*voigt(a(j,i),v)
         enddo                                     

         dummy1(i) = tauref(i)*kapnu(i)/(0.4343*kapref(i))
                                                       
c*****do the same for the strong lines
         if (dostrong .gt. 0) then
            do j=nlines+1,nlines+nstrong
               v = 2.997929d10*dabs(wave-wave1(j))/
     .             (wave1(j)*dopp(j,i)) 
               kapnu(i) = kapnu(i) + kapnu0(j,i)*voigt(a(j,i),v)    
            enddo
         endif

         dummy1(i) = tauref(i)*kapnu(i)/(0.4343*kapref(i))
      enddo      

c*****compute the optical depths                                            
      first = tauref(1)*kapnu(1)/kapref(1)
      dummy1(1) = rinteg(xref,dummy1,taunu,ntau,0.)      
      taunu(1) = first
      do i=2,ntau                                                     
         taunu(i) = taunu(i-1) + taunu(i)                  
      enddo


      return                                              
      end                                                

