
      subroutine cogplot
c******************************************************************************
c     This subroutine creates plots of curves-of-growth
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'
      real*4 gfplot(3000), rwplot(3000)
      real*4 xfour, style(1)
      character*4 ion


c*****dump the data into working arrays
      do i=1,ncurve
         gfplot(i) = gf1(i)
         rwplot(i) = w(i)
      enddo


c*****define the plot boundaries
      ylo = real(nint(rwlow*10))/10. - 0.1
      yhi = real(nint(rwhigh*10))/10. + 0.1
      do i=1,ncurve
         if (rwplot(i) .gt. rwlow) then
            if (gfplot(i) .gt. 0) then
               xlo = real(int(gfplot(i)*10))/10
            else
               xlo = real(int(gfplot(i)*10))/10 - 0.1
            endif
            go to 10
         endif
      enddo
10    do i=1,ncurve
         if (rwplot(i) .gt. rwhigh) then
            if (gfplot(i) .gt. 0) then
               xhi = real(int(gfplot(i)*10))/10 + 0.1
            else
               xhi = real(int(gfplot(i)*10))/10
            endif
         endif
      enddo


c*****start the plot via some setup calls
      call sm_limits (xlo,xhi,ylo,yhi)
      smlxtic = 0.20
      smlytic = 0.20
      bigxtic = 1.0
      bigytic = 1.0
      call sm_ticksize (smlxtic,bigxtic,smlytic,bigytic)
 
 
c*****draw and label the box for the curve-of-growth
      call sm_expand (0.6)
      call sm_lweight (1.4)
      call defcolor (1)
      call sm_box (0,0,0,0)
      call sm_expand (0.85)
      call sm_box (1,2,4,4)
      array = 'log gf'
      call sm_relocate (0.5*(xlo+xhi),ylo-0.10*(yhi-ylo))
      call sm_putlabel (5,array)
      array = 'log (EW/lambda)'
      call sm_relocate (xlo-0.06*(xhi-xlo),0.5*(yhi+ylo))
      call sm_angle (90.)
      call sm_putlabel (5,array)
      call sm_angle (0.)
      call sm_ltype (1)
      call sm_lweight (0.8)
      call sm_grid (0,0)
      call sm_ltype (0)


c*****plot the computed curve-of-growth points; exit normally
      call sm_expand (1.05)
      call defcolor (2)
      style(1) = 240.7
      call sm_ptype (style,1)
      call sm_points (gfplot,rwplot,ncurve)
      call defcolor (1)
      ich = idint(charge(lim1) + 0.1)
      if (ich .eq. 1) then
         ion = ' I  '
      elseif (ich .eq. 2) then
         ion = ' II '
      elseif (ich .eq. 3) then
         ion = ' III'
      endif
      iatom = idint(atom1(lim1))
      write (array,1002) names(iatom),ion
      call sm_expand (0.75)
      call sm_relocate (xhi-0.20*(xhi-xlo),ylo+0.35*(yhi-ylo))
      call sm_putlabel (5,array)
      call sm_relocate ((xhi+xlo)/2.0,yhi-0.03*(yhi-ylo))
      call sm_putlabel (5,moditle)
      xfour = dlog10(xabund(iatom)) + 12.0
      write (array,1003) xfour
      call sm_relocate (xhi-0.20*(xhi-xlo),ylo+0.27*(yhi-ylo))
      call sm_putlabel (5,array)
      write (array,1005) wave1(lim1)
      call sm_relocate (xhi-0.20*(xhi-xlo),ylo+0.19*(yhi-ylo))
      call sm_putlabel (5,array)
      write (array,1004) e(lim1,1) 
      call sm_relocate (xhi-0.20*(xhi-xlo),ylo+0.11*(yhi-ylo)) 
      call sm_putlabel (5,array)
      return


c*****format statements
1002  format (a2,a4)
1003  format ('log eps =',f7.3)
1004  format ('E.P. (eV) =',f7.3) 
1005  format ('lambda (A) =',f11.3)

      end




