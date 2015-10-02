
c******************************************************************************
c  The subroutines needed to calculate the opacities from scattering by
c  H I, H2, He I, and e- are in this file.  These are from ATLAS9.
c******************************************************************************






      subroutine opacescat
c******************************************************************************
c  This routine computes electron scattering opacities.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com' 
      include 'Kappa.com'
      save
      data modcount/0/

c  compute scattering, but only if there is a new model atmosphere.
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            sigel(i) = 0.6653d-24*ne(i)
         enddo
      endif

      return
      end


      subroutine opacHscat
c******************************************************************************
c  This routine computes H I Rayleigh scattering opacities.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      include 'Linex.com'

      wavetemp = 2.997925d18/dmin1(freq,2.463d15)
      ww = wavetemp**2
      sig = (5.799d-13+1.422d-6/ww+2.784/(ww*ww))/(ww*ww)
      do i=1,ntau
         sigH(i) = sig*2.*numdens(1,1,i)/u(1,1,i)
      enddo

      return
      end


      subroutine opacH2scat
c******************************************************************************
c  This routine computes H2 I Rayleigh scattering opacities.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      include 'Linex.com'

      wavetemp = 2.997925d18/dmin1(freq,2.463d15)
      ww = wavetemp**2
      sig = (8.14d-13+1.28d-6/ww+1.61/(ww*ww))/(ww*ww)
      do i=1,ntau
       sigH2(i) = sig*numdens(8,1,i)
      enddo

      return
      end


      subroutine opacHescat
c******************************************************************************
c  This routine computes He I Rayleigh scattering opacities.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      include 'Linex.com'

      wavetemp = 2.997925d18/dmin1(freq,5.15d15)
      ww = wavetemp**2
      sig = 5.484E-14/ww/ww*(1.+(2.44d5+5.94d10/(ww-2.90d5))/ww)**2
      do i=1,ntau
         sigHe(i) = sig*numdens(2,1,i)/u(2,1,i)
      enddo

      return
      end









