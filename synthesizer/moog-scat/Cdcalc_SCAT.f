      subroutine cdcalc_SCAT (number)
c*****************************************************************************************
c     Calculates the line depth and the corresponding quantity, adepth
c     (which is *similar* to the quantity, cd(i), in MOOG-LTE).
c     The quantity 'adepth' is the depth of the spectral line at a specified 
c     wavelength. Essentially, the integral of this adepth quantity over the wavelength of 
c     the spectral feature is the equivalent width. 
c******************************************************************************************
      implicit real*8 (a-h,o-z)
      real*8  Residual

      include 'Atmos.com'
      include 'Linex.com'
      include 'Scat.com'

      if (number .eq. 1) then
         call sourcefunc_cont_scat
         return
      else  
         call sourcefunc_line_scat
      endif
      Residual = Flux_line/Flux_cont
      adepth   = 1.- Residual
      return
      end

 
