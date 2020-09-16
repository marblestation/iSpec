
      real*8 function expint(x,n)
c******************************************************************************
c     This routine computes the n'th exponential integral of x
c     ******* notice limits on x for this program!!! ******
c     input - x,  independent variable (-100. .le. x .le. +100.)
c            n,  order of desired exponential integral (1 .le. n .le. 8)
c     output - expint,  the desired result
c     jd does not output ex to save moog call re-write
c              rex,  expf(-x)
c     jd
c     note   returns with e1(0)=0, (not infinity).
c     based on the share routine nuexpi, written by j. w. cooley,
c     courant institute of mathematical sciences, new york university
c     obtained from rudolf loeser
c     general compilation of 1 august 1967.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      real*8 tab(20),xint(7)
      data xint/1.,2.,3.,4.,5.,6.,7./
      data tab /  .2707662555,.2131473101,.1746297218,.1477309984,
     1.1280843565,.1131470205,.1014028126,.0919145454,.0840790292,
     1.0774922515,.0718735405,.0670215610,.0627878642,.0590604044,
     1.0557529077,.0527977953,.0501413386,.0477402600,.0455592945,
     1.0435694088/
      data xsave /0./
c
      if (x.ge.100.) goto 800 
      if (x.le.-100.) goto 800 
      xu=x
      if(xu)603,602,603
  602 rex=1.d+00
      if(n-1)800,800,801
  800 expint=0.
      goto 777
  801 expint=1./xint(n-1)
      goto 777
  603 if(xu-xsave)604,503,604
  604 xsave=xu
      xm=-xu
      emx = dexp(xm)
c
c  select method for computing ei(xm)
c
      if(xm-24.5)501,400,400
  501 if(xm-5.)502,300,300
  502 if(xm+1.)100,200,200
  503 eisave=-arg
      exsave=emx
c
c  now recurse to higher orders
c
      if(n-1)507,507,505
 505  do 506 i=2,n
        eisave=(xu*eisave-exsave)/(-xint(i-1))
  506 continue
  507 expint=eisave
      rex=exsave
  777 return
c
c  ei(xm) for xm .lt. -1.0
c  hastings polynomial approximation
c
  100 arg=((((((xu+8.573328740 )*xu+18.05901697  )*xu+8.634760893 )*xu
     *+.2677737343)/xm)*emx)/((((xu+9.573322345 )*xu+25.63295615  )*xu
     *+21.09965308  )*xu+3.958496923 )
      goto 503
c     ei(xm) for -1. .le. xm .lt. 5.0
c     power series expansion about zero
  200 arg=dlog(abs(xm))
      arg=((((((((((((((((.41159050d-14*xm+.71745406d-13)*xm+.76404637d-
     *12)*xm+.11395905d-10)*xm+.17540077d-9)*xm+.23002666d-8)*xm+.275360
     *18d-7)*xm+.30588626d-6)*xm+.31003842d-5)*xm+.28346991d-4)*xm+.2314
     *8057d-3)*xm+.0016666574)*xm+.010416668)*xm+.055555572)*xm+.25)*xm+
     *.99999999)*xm+.57721566)+arg
      goto 503
c
c  ei(xm) for 5.0 .le. xm .lt. 24.5
c  table look-up and interpolation
c
  300 i=xm+.5
      xzero=i
      xdelta=xzero-xm
      arg=tab(i-4)
      if(xdelta)303,305,303
  303 y=arg
      deltax=xdelta/xzero
      power=1./deltax
      do 304 i=1,7
        power=power*deltax
        y=((y-power/xzero)*xdelta)/xint(i)
        arg=arg+y
        if(abs(y/arg)-1.d-8)305,304,304
  304 continue
  305 arg=emx*arg
      goto 503
c     ei(xm) for 24.5 .le. xm
c     truncated continued fraction
  400 arg=((((xm-15.)*xm+58.)*xm-50.)*emx)/((((xm-16.)*xm+72.)*xm-96.)
     **xm+24.)
      goto 503
      end
