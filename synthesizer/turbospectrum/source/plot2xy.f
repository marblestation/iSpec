      program plot2xy
*
* writes spectrum plot format output in 2 colums to file
*
*               Olle Morell 911203
*                         
      real*8 wstart,step,x
      dimension irec(103), rec(103)
      integer pixel
      equivalence (irec(1), rec(1))
*
      open(10,file='syntspec/synoutf',
     &        status='old', form='unformatted')
      open(11,file='syntspec/synoutf.xy',
     &        status='unknown', form='formatted')
*
      read(10) rec
      wstart = rec(1)
      step = rec(2)
      pixel=0
*
      do while(.true.)
        read(10, end=999) rec
        do j=2, irec(1)+1
          x = wstart+step*pixel
          y = rec(j)
          write(11,1010) x, y
 1010     format(2f10.4)
          pixel = pixel+1
        end do
      end do
*
 999  continue       
      end
