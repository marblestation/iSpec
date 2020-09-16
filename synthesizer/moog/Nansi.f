
       integer function ivcleof(y,x)
c******************************************************************************
c    This routine clears to the end of the screen, beginning with the
c    position row=y, column=x. The cursor is then placed at (y,x)
c******************************************************************************

       include 'Pstuff.com'
       character blank*79
       integer y,x,ypos,xpos
 
       write(blank,1001)
1001   format(79(' '))
       istat = ivwrite(y,x,blank,79)
       xpos = 1
       ypos = y


       do 11 ipos=ypos,maxline
          istat = ivmove(ipos,xpos)
11        istat = ivcleol()
       istat = ivmove(y,x)
c
       ivcleof = 0
       return
       end





       integer function ivwrite(y,x,arr,ccount)
c
c    This routine writes out a string of characters from 'arr',
c    with the first character beginning at row=y, column=x. The length
c    of the string may be at most 79 characters, and this routine
c    will not write on the 80th column or the 25th row of the screen.
c
       include 'Pstuff.com'
       integer y,x,count,ccount
       character arr*(*),dummy*80,string(80)*1,esc*1
       equivalence (dummy,string(1))
c
       count = ccount
       if (y .gt. maxline .or. x .gt. 79) then
          ivwrite = -1
          return
       endif
c
       esc = char(27)
       dummy = arr
       count = min0(80-x,count)
c       
       if (x .lt. 10) then
          if (y .lt. 10) then
             write (*,1007) esc,y,x,(string(i),i=1,count)
1007         format(1x,a1,'[',i1,';',i1,'H',80a1)
          else
             write (*,1006) esc,y,x,(string(i),i=1,count)
1006         format(1x,a1,'[',i2,';',i1,'H',80a1)
          endif
       else 
          if (y .lt. 10) then
             write (*,1005) esc,y,x,(string(i),i=1,count)
1005         format(1x,a1,'[',i1,';',i2,'H',80a1)
          else
             write (*,1004) esc,y,x,(string(i),i=1,count)
1004         format(1x,a1,'[',i2,';',i2,'H',80a1)
          endif
       endif
c
       ivwrite = 0
       return
       end






       integer function ivmove(y,x)
c
c    This routine moves the cursor to position row=y, column=x.
c    It checks to be sure that y does not exceed 'maxline', and x does not
c    exceed 79 (the screen limits).
c
       include 'Pstuff.com'
       integer y,x
       character esc*1
c
       if (y .gt. maxline .or. x .gt. 79) then
          ivmove = -1
          return
       endif
c
       esc = char(27)
c       
       if (x .lt. 10) then
          if (y .lt. 10) then
             write (*,1007) esc,y,x
1007         format(1x,a1,'[',i1,';',i1,'H')
          else
             write (*,1006) esc,y,x
1006         format(1x,a1,'[',i2,';',i1,'H')
          endif
       else 
          if (y .lt. 10) then
             write (*,1005) esc,y,x
1005         format(1x,a1,'[',i1,';',i2,'H')
          else
             write (*,1004) esc,y,x
1004         format(1x,a1,'[',i2,';',i2,'H')
          endif
       endif
c
       ivmove = 0
       return
       end



    


       integer function ivcleol()
c    This routine clears to the end of a line, beginning with the current
c    cursor position.
c
       character*1 esc
c
       esc = char(27)
       write(*,1001) esc
1001   format(1x,a1,'[K')
c
       ivcleol = 0
       return
       end




       integer function ivbold(on)
c    
c    This routine either turns on the bold-face lettering, if on=1,
c    or turns it off, if on=0
c
       integer on
       character*1 esc
c
       esc = char(27)
       if(on.eq.1) then
             write(*,'(1x,a1,a)') esc,'[1m'
       else
             write(*,'(1x,a1,a)') esc,'[0m'
       endif 
c
       ivbold = 0
       return
       end





       integer function ivundl(on)
c
c    This routine either turns on the underlining of text, if on=1,
c    or turns it off.
c
       integer on
       character*1 esc
c
       esc = char(27)
       if(on.eq.1) then
             write(*,'(1x,a1,a)') esc,'[4m'
       else
             write(*,'(1x,a1,a)') esc,'[0m'
       endif 
c
       ivundl = 0
       return
       end

 
