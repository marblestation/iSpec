      subroutine minimax(x,xmin,xmax,numpts)
c*****this routine simply finds the minimum and maximum values of array x

      real*8 x(*)
 
      xmax=x(1)
      xmin=xmax
      do 10 i=1,numpts
         if (x(i).lt.xmin) xmin=x(i)
         if (x(i).gt.xmax) xmax=x(i)
   10 continue

      return
      end


