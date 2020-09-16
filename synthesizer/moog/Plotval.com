      common /plotval/ xleft,right,up,down,oldlft,oldrgt,oldup,olddown,
     .                 xmax,xmin,ixlcnt,iylcnt,realcolors(50)
      
      real*4 xleft,right,up,down,oldlft,oldrgt,oldup,olddown,
     .                 xmax,xmin,realcolors
  
