
      subroutine pltabun
c******************************************************************************
c     This subroutine controls the decisions that are made around the
c     plotting of individual line abundances
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'
      real*8 xnum
      real*4 shortnum


      if (kount .ge. plotopt) then
         if (plotopt .eq. 0) return
10       choice = 'y'
         plotroutine = 'term_port_abun'
         lscreen = maxline -2
         call makeplot (lscreen)
         array = 'WHAT TO DO NEXT ([n]/h/f/r/m/v)? '
         lscreen = lscreen + 2
         nchars = 35
         call getasci (nchars,maxline)
         choice = chinfo(1:1)
         if (choice.eq.'n' .or. nchars.le.0) then
            return
         elseif (choice .eq. 'm') then
            return   
         elseif (choice .eq. 'v') then
            write (array,*) 'What is the new microturbulence (km/s)? '
            nchars = 41
            call getnum (nchars,lscreen,xnum,shortnum)
            do i=1,ntau
               vturb(i) = xnum*1.0e5
            enddo
            write (moditle(57:64),1010) xnum
            lim1line = 0
            lim2line = 0
            lim1obs = 0
            lim2obs = 0
            return   
         elseif (choice .eq. 'h') then
            plotroutine = 'hard_port_abun'
            call makeplot (lscreen)
         elseif (choice .eq. 'r') then
            go to 10
         elseif (choice .eq. 'f') then
            plotroutine = 'file_port_abun'
            call makeplot (lscreen)
         endif
      endif


      return


c*****format statements
1010  format ('vt=',f5.2)
      end


