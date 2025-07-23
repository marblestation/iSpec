
      subroutine total
c******************************************************************************
c     This routine integrates to get an equivalent width.                   
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      real*8 waveint(5000), depthint(5000)
      equivalence (waveint,dummy3), (depthint,dummy4)


c*****compute the wavelength array 
      ntotal = (sstop - start)/step + 1.3
      if (ntotal .gt. 5000) then
         write (nf1out,1002) ntotal
         write (nf2out,1002) ntotal
         return
      endif
      do i=1,ntotal                                                    
         waveint(i) = start + step*(i-1)    
      enddo


c*****use the RINTEG routine to do an integration
      answer = 1000.*rinteg(waveint,d,depthint,ntotal,0.)         
      write (nf1out,1001) answer, ntotal                               
      write (nf2out,1001) answer, ntotal                               
      w(ncurve) = answer/1000.


c*****Then recheck using Simpson's rule
      ntot = ntotal
      if(ntotal/2*2 - ntotal .eq. 0) ntot = ntotal - 1
      answer = d(1) + 4.*d(2) + d(ntot)
      ntot = ntot - 2
      do i=3,ntot,2                                                   
         answer = answer + 2.*d(i) + 4.*d(i+1)
      enddo
      answer = 1000.*step/3.*answer
      write (nf1out,1006) answer
      write (nf2out,1006) answer
      return


c*****format statements
1001  format (' equivalent width: ',f8.2,' mA.',5x,
     .        5x,'number of points: ',i5)
1002  format (i7,' POINTS ARE TOO MANY FOR THE INTEGRATION')
1006  format(' Simpson rule check on equivalent width =',f8.2)          


      end







