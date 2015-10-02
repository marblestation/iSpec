
      subroutine drawcurs
c******************************************************************************
c     This subroutine draws an arrow and writes the user x- and y-positions
c     of the cursor upon a click
c******************************************************************************

      include 'Pstuff.com'
      integer ichr


      call sm_graphics
      if     (whichwin .eq. '1of1') then
         call sm_window (1,1,1,1,1,1)
      elseif (whichwin .eq. '2of2') then
         call sm_defvar ('y_gutter','0.0')
         call sm_window (1,2,1,1,1,1)
      endif
      call sm_curs (xplotpos,yplotpos,ichr)


      call sm_relocate (xplotpos,yplotpos-0.10*(yhi-ylo))
      call sm_draw (xplotpos,yplotpos)
      call sm_draw (xplotpos-0.01*(xhi-xlo),yplotpos-0.03*(yhi-ylo))
      call sm_relocate (xplotpos,yplotpos)
      call sm_draw (xplotpos+0.01*(xhi-xlo),yplotpos-0.03*(yhi-ylo))
      call writenumber (xplotpos)
      call sm_expand (0.6)
      call sm_relocate (xplotpos,yplotpos-0.11*(yhi-ylo))
      call sm_putlabel (5,array)
      call writenumber (yplotpos)
      if     (whichwin(4:4) .eq. '1') then
         call sm_relocate (xplotpos,yplotpos-0.15*(yhi-ylo))
      elseif (whichwin(4:4) .eq. '2') then
         call sm_relocate (xplotpos,yplotpos-0.18*(yhi-ylo))
      elseif (whichwin(4:4) .eq. '3') then
         call sm_relocate (xplotpos,yplotpos-0.21*(yhi-ylo))
      endif
      call sm_putlabel (5,array)


      call sm_gflush
      call sm_alpha


      
      return
      end



