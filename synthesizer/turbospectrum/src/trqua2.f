************************************************************************
      function trqua2 (n,x,f,dx,d)
cc      dimension x(*),f(*),dx(*),d(*)
      dimension x(n),f(n),dx(n),d(n)
*
*  Calculate the integral over f(x), using dx and d as intermediaries
*  holding the x-differences and the derivatives of f.
*
*  The quadrature is performed with a formula for cubic splines, which
*  may be decomposed into a trapeziodal formula and a correction, where
*  the correction term in each interval is 1/12 times dx**2 times the
*  derivative change across the interval.
*
*  We assume f(x) is an even function of x, and that x(n)=0.0
*
************************************************************************
*
      res=0.
      do 191 i=2,n-1
        dx(i)=x(i)-x(i-1)
        d(i)=(f(i+1)-f(i-1))/(x(i+1)-x(i-1))
191     res=res+dx(i)*(f(i)+f(i-1))
      dx(n)=x(n)-x(n-1)
      res=res+dx(n)*(f(n)+f(n-1))
*
*  This is because f is an even function of x at x(n)=0
*
      d(n)=0.
      d(n-1)=2.*(f(n)-f(n-1))/dx(n)
*
*  Extrapolate derivative to x=1
*
      d1=(f(2)-f(1))/dx(2)
      d2=(f(3)-f(2))/dx(3)
      p=-dx(2)/(dx(2)+dx(3))
      d(1)=p*d2+(1.-p)*d1

      res=res*6.
      do i=2,n
        res=res+(d(i-1)-d(i))*dx(i)**2
      enddo
      res=res*0.083333333
      trqua2=res

      end
