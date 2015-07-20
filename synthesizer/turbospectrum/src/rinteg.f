cc      real*8 function rinteg(x,f,fint,n,start)
      function rinteg(x,f,fint,n,start)
c******************************************************************************
c     This routine is from ATLAS6
* Copied from MOOG 06/04-2001
c******************************************************************************

cc      implicit real*8 (a-h,o-z)
      dimension x(5000), f(5000), fint(5000)
      dimension a(5000), b(5000), c(5000)

      call parcoe (f,x,a,b,c,n)
      fint(1) = start 
      rinteg = start
      n1 = n - 1
      do 10 i=1,n1
         fint(i+1)= (a(i)+b(i)/2.*(x(i+1)+x(i))+ 
     .     c(i)/3.*((x(i+1)+x(i))*x(i+1)+x(i)*x(i)))*(x(i+1)-x(i))
10    rinteg = rinteg + fint(i+1)

      return
      end 




      subroutine parcoe(f,x,a,b,c,n)

cc      implicit real*8 (a-h,o-z)
      dimension f(5000), x(5000), a(5000), b(5000), c(5000)

      c(1)=0.
      b(1)=(f(2)-f(1))/(x(2)-x(1))
      a(1)=f(1)-x(1)*b(1)
      n1=n-1
      c(n)=0.
      b(n)=(f(n)-f(n1))/(x(n)-x(n1))
      a(n)=f(n)-x(n)*b(n) 
      if(n.eq.2)return
      do 1 j=2,n1
      j1=j-1
      d=(f(j)-f(j1))/(x(j)-x(j1)) 
      c(j)=f(j+1)/((x(j+1)-x(j))*(x(j+1)-x(j1)))-f(j)/((x(j)-x(j1))*
     1(x(j+1)-x(j)))+f(j1)/((x(j)-x(j1))*(x(j+1)-x(j1)))
      b(j)=d-(x(j)+x(j1))*c(j)
    1 a(j)=f(j1)-x(j1)*d+x(j)*x(j1)*c(j)
      c(2)=0.
      b(2)=(f(3)-f(2))/(x(3)-x(2))
      a(2)=f(2)-x(2)*b(2)
      c(3)=0.
      b(3)=(f(4)-f(3))/(x(4)-x(3))
      a(3)=f(3)-x(3)*b(3)
      if(c(j).eq.0.)go to 2
      j1=j+1
      wt=abs(c(j1))/(abs(c(j1))+abs(c(j)))
      a(j)=a(j1)+wt*(a(j)-a(j1))
      b(j)=b(j1)+wt*(b(j)-b(j1))
      c(j)=c(j1)+wt*(c(j)-c(j1))
    2 continue
      a(n1)=a(n)
      b(n1)=b(n)
      c(n1)=c(n)
      return
      end 



