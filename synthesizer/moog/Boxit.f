 
      subroutine boxit
c******************************************************************************
c     This subroutine figures out the point numbers of the wavelength 
c     boundaries for synthetic spectrum plots
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'
      include 'Linex.com'
      include 'Equivs.com'


c*****if only a synthetic spectrum is being plotted, use beginning and
c     end points of the synthetic spectrum instead of those of the observed
c     spectrum
      if (plotopt .lt. 2) then
         do i=1,kount
            if (xlo .le. xsyn(i)) then
               lim1obs = i
               go to 10
            endif
         enddo
         lim1obs = kount
10       do i=lim1obs,kount
            if (xhi .lt. xsyn(i)) then
               lim2obs = i -1
               return
            endif
         enddo
         lim2obs = kount
         return


c*****and here is the same logic when an observed spectrum exists
      else
         do i=1,lount
            if (xlo .le. xobs(i)) then
               lim1obs = i
               go to 20
            endif
         enddo
         lim1obs = lount
20       do i=lim1obs,lount
            if (xhi .lt. xobs(i)) then
               lim2obs = i -1
               return
            endif
         enddo
         lim2obs = lount
         return
      endif

      end




