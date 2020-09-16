
      subroutine specplot
c******************************************************************************
c     This routine produces MONGO plots of the syntheses.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Quants.com'
      include 'Factor.com'
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      include 'Multistar.com'
      real*4 style(1)
      real*4 yup,ydown
      integer iflip


c*****for grid syntheses, dump out relevant information to a file
      if (choice .eq. 'g') then
         write (nf6out,3001) syncount
         write (nf6out,3002) obsitle, moditle, linitle, smitle
      endif
 

c*****begin with a default window
      call sm_location (3500,31000,4000,30000)
      call sm_window (1,1,1,1,1,1)
      call sm_limits (0.0,1.0,0.0,1.0)
      call defcolor (1)


c*****write smoothing information at the top of the plot
      call sm_lweight (2.2)
      call sm_expand (0.7)
      call sm_relocate (-0.120,1.015)
      call sm_label (smitle)


      if (isoitle(1:10) .eq. '          ') then
         isoitle(1:16) = 'no isotopic data'
      endif
      call sm_relocate (-0.120,1.075)
      call sm_label (isoitle(1:120))
      if (numiso .gt. 3) then
         call sm_relocate (-0.120,1.045)
         call sm_label (isoitle(121:240))
      endif
      if (control .eq. 'gridplo' .or.
     .    control .eq. 'gridsyn' .or.
     .    control .eq. 'gridend') then
         write (nf6out,3002) isoitle(1:120)
         write (nf6out,3002) isoitle(121:240)
      endif
 

c*****define the real plot limits
      if (xlo .lt. xhi) then
         call sm_limits (xlo,xhi,ylo,yhi)
         iflip = 0
      else
         call sm_limits (xhi,xlo,ylo,yhi)
         iflip = 1
      endif
      call findtic (xlo,xhi,bigxtic,smlxtic)
      call findtic (ylo,yhi,bigytic,smlytic)
      call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)


c*****draw and label the box for the spectra
      call defcolor (1)
      if (whichwin .eq. '1of1') then
         idev = 1
         call sm_window (1,1,1,1,1,1)
      else
         idev = 2
         call sm_defvar ('y_gutter','0.0')
         call sm_window (1,2,1,1,1,1)
      endif
      call sm_lweight (4.0)
      call sm_expand (1.2)
      call sm_box (0,0,0,0)
      call sm_lweight (2.0)
      call sm_expand (0.8)
      call sm_box (1,2,4,4)
      if (iflip .eq. 1) then
         array = 'Wavenumber'
      else
         array = 'Wavelength'
      endif
      call sm_xlabel (array)
      array = 'Rel  Flux'
      call sm_ylabel (array)


c*****plot the synthetic spectra
      call sm_lweight (2.2)
      call sm_expand (0.7)
      do i=1,100
         if (pec(i) .ne. 0) go to 111
      enddo            
111   do j=1,nsyn
         if (choice.eq.'h' .or. choice.eq.'f' .or.
     .       choice.eq.'g') then
            call defcolor (8)
            call sm_ltype (j-1)
         else
            if (smterm(1:3) .eq. 'x11') then
               call defcolor (j+1)
               call sm_ltype (0)
            else
               call defcolor (1)
               call sm_ltype (j-1)
            endif
         endif
         call sm_connect (xsyn,chunk(1,j),kount) 
         if (iflip .eq. 1) then
            call sm_relocate (xhi+0.045*(xlo-xhi),
     .                     ylo+(0.12+0.06*j)*(yhi-ylo))
            call sm_draw (xhi+0.005*(xlo-xhi),
     .                 ylo+(0.12+0.06*j)*(yhi-ylo))
            call sm_relocate (xhi+0.05*(xlo-xhi),
     .                     ylo+(0.12+0.06*j)*(yhi-ylo))
         else
            call sm_relocate (xlo+0.045*(xhi-xlo),
     .                     ylo+(0.12+0.06*j)*(yhi-ylo))
            call sm_draw (xlo+0.005*(xlo-xhi),
     .                 ylo+(0.12+0.06*j)*(yhi-ylo))
            call sm_relocate (xlo+0.05*(xhi-xlo),
     .                     ylo+(0.12+0.06*j)*(yhi-ylo))
         endif
         noff = 80*(j-1)
         call sm_lweight (2.2)
         call sm_label (abitle(noff+1:noff+80))
         if ((control .eq. 'gridplo' .or.
     .        control .eq. 'gridsyn' .or.
     .        control .eq. 'gridend') .and. 
     .        whichwin.eq.'1of1') then
            write (nf6out,3002) abitle(noff+1:noff+80)
         endif
      enddo
      call defcolor (1)
      call sm_ltype (0)

   
c*****plot the observed spectrum
      if (plotopt .eq. 2) then
         call defcolor (1)
         if (choice.eq.'h' .or. choice.eq.'f' .or.
     .       choice.eq.'g') then
            call sm_lweight (4.0)
         else 
            call sm_lweight (2.2)
         endif
         call sm_ltype (0)
         call sm_expand (3.0)
         style(1) = 43.5
         call sm_ptype (style,1)
         mount = lim2obs - lim1obs + 1
         if (mount .lt. 500) then
            call sm_points (xobs(lim1obs),yobs(lim1obs),mount)
         else
            if (histoyes .eq. 1) then
               call sm_histogram (xobs(lim1obs),yobs(lim1obs),mount)
            else
               call sm_connect (xobs(lim1obs),yobs(lim1obs),mount)
            endif
         endif
         call sm_lweight (2.2)
         call sm_expand (0.7)
         if (iflip .eq. 1) then
            call sm_relocate (xhi+0.05*(xlo-xhi),ylo+0.12*(yhi-ylo))
         else
            call sm_relocate (xlo+0.05*(xhi-xlo),ylo+0.12*(yhi-ylo))
         endif
         call sm_label (obsitle)
      endif
      if (iflip .eq. 1) then
         call sm_relocate (xhi+0.05*(xlo-xhi),ylo+0.06*(yhi-ylo))
      else
         call sm_relocate (xlo+0.05*(xhi-xlo),ylo+0.06*(yhi-ylo))
      endif
      call sm_label (moditle)
      if (whichwin.eq.'1of1' .or. plotopt.ne.2) then
         return
      endif


c*****this section of code is executed only if a deviations plot is desired;
c     find the starting and stopping points in the arrays for the deviations
      if (xsyn(kount) .le. xobs(lim1obs)) go to 1000
      if (xsyn(1) .gt. xobs(lim2obs)) go to 1000
      if (xsyn(1) .gt. xobs(lim1obs)) go to 150
      lim3obs = lim1obs
      do k=2,kount
         if (xsyn(k) .gt. xobs(lim3obs)) then
            lim1syn = k - 1
            go to 155
         endif
      enddo
150   lim1syn = 1
      do l=lim1obs,lim2obs
         if (xsyn(lim1syn) .le. xobs(l)) then
            lim3obs = l 
            go to 155
         endif
      enddo
155   if (xsyn(kount) .lt. xobs(lim2obs)) go to 160
      lim4obs = lim2obs
      do k=lim1syn,kount
         if (xsyn(k) .gt. xobs(lim4obs)) then
            lim2syn = k
            go to 165
         endif
      enddo
160   lim2syn = kount
      do l=lim3obs,lim2obs
         if (xsyn(lim2syn) .lt. xobs(l)) then
            lim4obs = l - 1
            go to 165
         endif
      enddo


c  compute the deviations; linear interpolation in the wavelength array
c  of the synthetic spectra is considered sufficient
165   open (222,file='devfile')
      do j=1,nsyn
         lpoint = lim1syn
         devsigma = 0.
         do i=lim3obs,lim4obs
170         if (xsyn(lpoint+1) .lt. xobs(i)) then
               lpoint = lpoint + 1
               go to 170
            endif
            syninterp = (chunk(lpoint+1,j)-chunk(lpoint,j))*
     .         (xobs(i)-xsyn(lpoint))/(xsyn(lpoint+1)-xsyn(lpoint)) +
     .         chunk(lpoint,j)
            dev(i) = yobs(i) - syninterp
            devsigma = devsigma + dev(i)**2
         enddo
         devsigma = dsqrt(devsigma/(lim4obs-lim3obs-1))
         write (222,1222) j
1222     format ('deviation array ', i1)
         write (222,1223) (xobs(i),dev(i),i=lim3obs,lim4obs)
1223     format (f10.3, f8.3)
         
      
c  from first set of deviations, define the plot limits, draw and label box
         if (j .eq. 1) then
            yup = -1000.
            ydown = +1000.
            do i=lim3obs,lim4obs
               yup = amax1(yup,dev(i))
               ydown = amin1(ydown,dev(i))
            enddo
            ydelta = amin1(0.3,amax1(0.05,1.5*(yup-ydown)/2.))
            ydown = -ydelta
            yup = ydelta
            call sm_defvar ('y_gutter','0.0')
            call sm_window (1,2,1,2,1,2)
            call sm_limits (xlo,xhi,ydown,yup)
            call findtic (ydown,yup,bigytic,smlytic)
            call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)
            call sm_lweight (4.0)
            call sm_expand (1.2)
            call defcolor (1)
            call sm_box (0,0,0,0)
            call sm_lweight (2.0)
            call sm_expand (0.8)
            call sm_box (4,2,4,4)
            array = 'Obs - Comp'
            call sm_ylabel (array)
            call sm_relocate (xlo,0.0)
            call sm_draw (xhi,0.0)
         endif


c  plot the array of deviations
         if (choice.eq.'h' .or. choice.eq.'f' .or.
     .       choice.eq.'g') then
            call defcolor (8)
            call sm_ltype (j-1)
         else
            if (smterm(1:3) .eq. 'x11') then
               call defcolor (j+1)
               call sm_ltype (0)
            else
               call defcolor (1)
               call sm_ltype (j-1)
            endif
         endif
         call defcolor (j+1)
         call sm_lweight (2.2)
         call sm_connect (xobs(lim1obs),dev(lim1obs),mount)
         write (array,1021) devsigma
         call sm_relocate(xhi-0.2*(xhi-xlo),
     .                 ydown+(0.10+0.06*j)*(yup-ydown))
         call sm_relocate(xhi-0.24*(xhi-xlo),
     .                 ydown+(0.10+0.06*j)*(yup-ydown))
         call sm_draw(xhi-0.215*(xhi-xlo),
     .                 ydown+(0.10+0.06*j)*(yup-ydown))
         call sm_label (array)
         if (choice .eq. 'g') then
            noff = 80*(j-1)
            write (nf6out,3002) abitle(noff+1:noff+80)
            write (nf6out,3003) devsigma, velsh
         endif
      enddo
      close (unit=222)


c  reset the spectrum plot boundaries before exiting
      if(xlo .lt. xhi) then
         call sm_limits (xlo,xhi,ylo,yhi)
         iflip = 0
      else
         call sm_limits (xhi,xlo,ylo,yhi)
         iflip = 1
      endif
1000  return


c*****format statements
1021  format ('sigma =',f7.4)
3001  format (/'RUN NUMBER:', i4, 65('-'))
3002  format (a80)
3003  format ('sigma =',f8.5, '   lambda shift =',f7.3,' km/s')


      end



