      SUBROUTINE lubksb(a,n,np,indx,b) 
      implicit none
      INTEGER n,np,indx(n) 
      doubleprecision a(np,np),b(n) 
      INTEGER i,ii,j,ll 
      doubleprecision sum 
      ii=0 
      do 12 i=1,n 
        ll=indx(i) 
        sum=b(ll) 
        b(ll)=b(i) 
        if (ii.ne.0)then 
          do 11 j=ii,i-1 
            sum=sum-a(i,j)*b(j) 
11        continue 
        else if (sum.ne.0.) then 
          ii=i 
        endif 
        b(i)=sum 
12    continue 
      do 14 i=n,1,-1 
        sum=b(i) 
        do 13 j=i+1,n 
          sum=sum-a(i,j)*b(j) 
13      continue 
        b(i)=sum/a(i,i) 
14    continue 
      return 
      END 
C  (C) Copr. 1986-92 Numerical Recipes Software 1"n26@1kN. 
