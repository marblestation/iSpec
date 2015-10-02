 
      subroutine pltcog
c******************************************************************************
c     This subroutine controls the decisions that are made around the
c     plots of curves-of-growth
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'
      real*8 xnum
      real*4 shortnum


c  call up the curve-of-growth plot
      if (plotopt .eq. 0) return
10    choice = 'y'
      plotroutine = 'term_land_cog '
      lscreen = 12
      whichwin = '1of1'
      call makeplot (lscreen)


c  make a hardcopy, write to a postscript file, point at a place on the
c  screen, try a new model atmosphere, or replot?
1     array = 'WHAT TO DO NEXT ([n]/h/f/p/v/m/r)? '
      lscreen = 12
      nchars = 37
      call getasci (nchars,lscreen)
      choice = chinfo(1:1)
      if (choice.eq.'n' .or. nchars.le.0) then
         return
      elseif (choice .eq. 'h') then
         plotroutine = 'hard_land_cog '
         call makeplot (lscreen)
         go to 10
      elseif (choice .eq. 'v') then
         write (array,*) 'What is the new microturbulence (km/s)? '
         nchars = 41
         lscreen = lscreen + 2
         call getnum (nchars,lscreen,xnum,shortnum)
         do i=1,ntau
            vturb(i) = xnum*1.0e5
         enddo
         write (moditle(57:64),1010) xnum
         rewind nf1out
         rewind nf2out
         rewind nfmodel
         rewind nflines
         return
      elseif (choice .eq. 'm') then
         return
      elseif (choice .eq. 'r') then
         go to 10
      elseif (choice .eq. 'p') then
         array = 'MARK THE POSITION WITH THE CURSOR'
         istat=ivcleof(21,1)
         istat=ivwrite(13,3,array,34)
         call drawcurs
         go to 1
      elseif (choice .eq. 'f') then
         plotroutine = 'file_land_cog '
         call makeplot (lscreen)
         go to 10
      endif
      return


c*****format statements
1010  format ('vt=',f5.2)


      end




