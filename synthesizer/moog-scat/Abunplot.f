
      subroutine abunplot
c******************************************************************************
c     This routine produces MONGO plots of line abundances versus
c     excitation potentials and equivalent widths
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      real*4 ep(2500),abb(2500),logrw(2500),wavepl(2500)
      equivalence (ep,dummy3(1)),(abb,dummy3(1251)),(logrw,dummy3(2501))
      real*4 style(1),ymed
      character ion*4


c*****dump the data into working arrays
      j = 0
      do l=lim1obs,lim2obs
         if (abundout(l) .ne. 999.99) then
            j = j + 1
            ep(j) = e(l,1)
            abb(j) = abundout(l)
c           logrw(j) = dlog10(wid1comp(l)/wave1(l))
            logrw(j) = dlog10(width(l)/wave1(l))
            wavepl(j) = wave1(l)
         endif
      enddo


c*****find the plot boundaries for the excitation potential plot
      xlo = 20.
      do j=1,kount
         xlo = amin1(xlo,ep(j))
      enddo
      xhi =  0.
      do j=1,kount
         xhi = amax1(xhi,ep(j))
      enddo
      if (xhi-xlo .lt. 5.) then
         xlo = amax1((xlo+xhi)/2.-2.5,-0.2)
         xhi = xlo + 5.0
      else
         xlo = amax1((xlo+xhi)/2.-5.0,-0.2)
         xhi = xlo + 10.
      endif
      ylo = 20.
      do j=1,kount
         ylo = amin1(ylo,abb(j))
      enddo
      yhi = -20.
      do j=1,kount
         yhi = amax1(yhi,abb(j))
      enddo
      if (yhi-ylo .lt. 0.5) then
         ylo = (ylo+yhi)/2. - 0.30
         yhi = ylo + 0.60
      elseif (yhi-ylo .lt. 1.0) then
         ylo = (ylo+yhi)/2. - 0.55
         yhi = ylo + 1.10
      else
         ylo = (ylo+yhi)/2. - 1.10
         yhi = ylo + 2.20
      endif


c*****start the excitation potential plot via plot setup calls
      call sm_defvar ('y_gutter','0.9')
      call sm_window (1,3,1,3,1,3)
      call sm_limits (xlo,xhi,ylo,yhi)
      call findtic (xlo,xhi,bigxtic,smlxtic)
      call findtic (ylo,yhi,bigytic,smlytic)
      call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)


c*****draw and label the box for the excitation potential plot
      call defcolor (1)
      call sm_lweight (4.0)
      call sm_expand (1.2)
      call sm_box (0,0,0,0)
      call sm_lweight (2.0)
      call sm_expand (0.8)
      call sm_box (1,2,4,4)
      array = 'E. P. (eV)'
      call sm_relocate (0.5*(xlo+xhi),ylo-0.20*(yhi-ylo))
      call sm_putlabel (5,array)
      array = 'log eps'
      call sm_relocate (xlo-0.10*(xhi-xlo),0.5*(yhi+ylo))
      call sm_angle (90.)
      call sm_putlabel (5,array)
      call sm_angle (0.)
      

c*****make the excitation potential plot
      call defcolor (2)
      call sm_expand (2.2)
      style(1) = 43.5
      call sm_ptype (style,1) 
      call sm_points (ep,abb,kount)
      call defcolor (4)
      ymed = average
      call sm_lweight (4.0)
      call sm_ltype (2)
      call sm_relocate (xlo,ymed)
      call sm_draw (xhi,ymed)
      if (kount .gt. 2 .and. deltaep .gt. 1.5) then
         call sm_ltype (3)
         call defcolor (3)
         call sm_relocate (xlo,real(xxm1*xlo+xxb1)) 
         call sm_draw (xhi,real(xxm1*xhi+xxb1)) 
      endif
      call sm_ltype (0)
      call sm_lweight (2.0)
      call defcolor (1)
      call sm_relocate (xlo+0.05*(xhi-xlo),ylo+0.15*(yhi-ylo))
      call sm_expand (0.8)
      ich = idint(charge(lim1obs) + 0.1)
      if (ich .eq. 1) then 
         ion = ' I  '
      elseif (ich .eq. 2) then
         ion = ' II '
      elseif (ich .eq. 3) then
         ion = ' III'
      endif
      iatom = idint(atom1(lim1obs))
      write (array,1002) names(iatom),ion
      call sm_relocate ((xhi+xlo)/2.,yhi-0.12*(yhi-ylo))
      call sm_expand (0.8)
      call sm_putlabel (5,array)
   

c*****define the plot limits for the equivalent width plot
      xlo = 5000.
      do j=1,kount
         xlo = amin1(xlo,logrw(j))
      enddo
      xhi =  -5000.
      do j=1,kount
         xhi = amax1(xhi,logrw(j))
      enddo
      xlo = int((xlo-0.01)*2.-1.)/2.
      xhi = int((xhi+0.01)*2.)/2.
      call sm_defvar ('y_gutter','0.9')
      call sm_window (1,3,1,2,1,2)
      call sm_limits (xlo,xhi,ylo,yhi)
      call findtic (xlo,xhi,bigxtic,smlxtic)
      call findtic (ylo,yhi,bigytic,smlytic)
      call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)


c*****draw and label the box for the equivalent width plot
      call defcolor (1)
      call sm_lweight (4.0)
      call sm_expand (1.2)
      call sm_box (0,0,0,0)
      call sm_lweight (2.0)
      call sm_expand (0.8)
      call sm_box (1,2,4,4)
      array = 'log (EW/lambda)'
      call sm_relocate (0.5*(xlo+xhi),ylo-0.20*(yhi-ylo))
      call sm_putlabel (5,array)
      array = 'log eps'
      call sm_relocate (xlo-0.10*(xhi-xlo),0.5*(yhi+ylo))
      call sm_angle (90.)
      call sm_putlabel (5,array)
      call sm_angle (0.)


c*****make the equivalent width plot
      call defcolor (2)
      call sm_expand (2.2)
      style(1) = 43.5
      call sm_ptype (style,1)
      call sm_points (logrw,abb,kount)
      call defcolor (4)
      call sm_lweight (4.0)
      call sm_ltype (2)
      call sm_relocate (xlo,ymed)
      call sm_draw (xhi,ymed)
      if (kount .gt. 2 .and. deltarw .gt. 0.5) then
         call sm_ltype (3)
         call defcolor (3)
         call sm_relocate (xlo,real(xxm2*xlo+xxb2)) 
         call sm_draw (xhi,real(xxm2*xhi+xxb2)) 
      endif
      call sm_ltype (0)
      call sm_lweight (2.0)
      call defcolor (1)
      call sm_expand (0.8)
      call sm_relocate ((xhi+xlo)/2.,yhi-0.12*(yhi-ylo))
      call sm_putlabel (5,moditle(1:56))
      call sm_relocate ((xhi+xlo)/2.,ylo+0.12*(yhi-ylo))
      call sm_putlabel (5,moditle(57:80))


c*****define the plot limits for the wavelength plot
      xlo = 5000000.
      do j=1,kount
         xlo = amin1(xlo,wavepl(j))
      enddo
      xhi = 0.
      do j=1,kount
         xhi = amax1(xhi,wavepl(j))
      enddo
      xlo = int(xlo-50.)
      xhi = int(xhi+50.)
      call sm_defvar ('y_gutter','0.9')
      call sm_window (1,3,1,1,1,1)
      call sm_limits (xlo,xhi,ylo,yhi)
      call findtic (xlo,xhi,bigxtic,smlxtic)
      call findtic (ylo,yhi,bigytic,smlytic)
      call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)


c*****draw and label the box for the wavelength plot
      call defcolor (1)
      call sm_lweight (4.0)
      call sm_expand (1.2)
      call sm_box (0,0,0,0)
      call sm_lweight (2.0)
      call sm_expand (0.8)
      call sm_box (1,2,4,4)
      array = 'lambda (A)'
      call sm_relocate (0.5*(xlo+xhi),ylo-0.20*(yhi-ylo))
      call sm_putlabel (5,array)
      array = 'log eps' 
      call sm_relocate (xlo-0.10*(xhi-xlo),0.5*(yhi+ylo)) 
      call sm_angle (90.) 
      call sm_putlabel (5,array) 
      call sm_angle (0.) 


c*****make the wavelength plot, and exit normally
      call defcolor (2)
      call sm_expand (2.2)
      style(1) = 43.5 
      call sm_ptype (style,1)
      call sm_points (wavepl,abb,kount)
      call defcolor (4)
      call sm_relocate (xlo,ymed)
      call sm_lweight (4.0)
      call sm_ltype (2)
      call sm_draw (xhi,ymed)
      if (kount .gt. 2 .and. deltawv .gt. 500.) then
         call sm_ltype (3)
         call defcolor (3)
         call sm_relocate (xlo,real(xxm3*xlo+xxb3))
         call sm_draw (xhi,real(xxm3*xhi+xxb3))
      endif
      call sm_ltype (0)  
      call sm_lweight (2.0)
      call defcolor (1)
      call sm_relocate ((xhi+xlo)/2.,yhi-0.12*(yhi-ylo))
      call sm_expand (0.8)
      call sm_putlabel (5,linitle)
      return


c*****format statements
1002  format (a2,a4)


      end








