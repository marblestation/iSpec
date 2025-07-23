
      subroutine discov(amol,i1,i2)
c******************************************************************************
c     This routine returns i2 = the number of times atom i1 appears in      
c     molecule amol                                                         
c******************************************************************************

      implicit real*8 (a-h,o-z)
      integer itest(5)                                                
      data itest /100000000, 1000000, 10000, 100, 1/
      i2 = 0                                                            
      im = amol                                                         
      do i=1,5                                                        
         i3 = im/itest(i)                                              
         if(i3 .eq. i1) i2 = i2 + 1                                    
         im = im - i3*itest(i)                                         
      enddo
      return                                                            


      end                                                               


