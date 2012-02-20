#include <stdio.h>

double trapez(y,x,N)
double y[],x[];
int N;
{
 int i;
 double integ;

 integ = 0.0;
 if(N == 0) return(0.0);
 for(i=0;i<N;i++) integ += 0.5*(y[i+1]+y[i])*(x[i+1]-x[i]);
 return(integ);
}
