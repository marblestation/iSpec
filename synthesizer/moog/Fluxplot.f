 
      subroutine fluxplot
c******************************************************************************
c     This subroutine creates plots of flux curves
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Dummy.com'
      include 'Pstuff.com'
      real*8 wavep(2000), fluxp(2000)
      real*4 waveplot(2000), flxplt(2000)
      real*4 style(1)


c*****dump the data into working arrays
      rewind nf2out
      i = 1
1     read (nf2out,*,end=10) wavep(i), fluxp(i), waveplot(i), flxplt(i)
      i = i + 1
      go to 1


c*****define the plot boundaries
10    ntot = i - 1
      xlo = int(waveplot(1))
      xhi = int(waveplot(ntot))
      ylo =  1.e+30
      yhi = -1.e+30
      do i=1,ntot
         ylo = amin1(ylo,flxplt(i))
         yhi = amax1(yhi,flxplt(i))
      enddo
      ylo = real(int(ylo+0.5001)) - 0.5
      yhi = real(int(yhi+0.5001)) + 0.5


c*****start the plot via some setup calls
      call sm_location (3500,31000,5000,31000)
      call sm_limits (xlo,xhi,ylo,yhi)
      call findtic (xlo,xhi,bigxtic,smlxtic)
      call findtic (ylo,yhi,bigytic,smlytic)
      call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)
 
 
c*****draw and label the box for the curve-of-growth
      call sm_expand (0.6)
      call sm_lweight (1.4)
      call defcolor (1)
      call sm_box (0,0,0,0)
      call sm_expand (1.0)
      call sm_box (1,2,4,4)
      array = '1/lambda'
      call sm_relocate (0.5*(xlo+xhi),ylo-0.20*(yhi-ylo))
      call sm_putlabel (5,array)
      array = 'log (flux)'
      call sm_relocate (xlo-0.10*(xhi-xlo),0.5*(yhi+ylo))
      call sm_angle (90.)
      call sm_putlabel (5,array)
      call sm_angle (0.)
      call sm_ltype (1)
      call sm_lweight (0.8)
      call sm_grid (0,0)
      call sm_ltype (0)


c*****plot the computed flux curve points
      call sm_expand (1.0)
      call defcolor (2)
      style(1) = 240.7
      call sm_ptype (style,1)
      call sm_points (waveplot,flxplt,ntot)
      call defcolor (1)
      call sm_relocate ((xhi+xlo)/2.0,ylo+0.07*(yhi-ylo))
      call sm_putlabel (5,moditle)


c*****compute total flux and effective temperature; exit normally
      do i=1,ntot
         wavep(i) = 1.0d-8*wavep(i)
      enddo
      first =  fluxp(1)
      fluxtot = rinteg(wavep,fluxp,dummy3,ntot,first)
      teff = (3.14159*fluxtot/5.67d-5)**0.25
      write (smitle,1001) fluxtot, teff
      call defcolor (5)
      call sm_relocate ((xhi+xlo)/2.0,ylo+0.14*(yhi-ylo))
      call sm_putlabel (5,smitle)
      return


c*****format statements
1001  format ('total flux = ' ,1pd12.4, 5x, 'Teff = ', 0pf6.0)
      end




