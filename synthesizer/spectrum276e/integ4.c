#include <math.h>
#include <errno.h>
#include "spectrum.h"
int imin(),imax();

double integ(y,x,N)
double y[],x[];
int N;
{
 int i,j,k;
 double integral = 0.0;
 double a[2],b[2],c[2],am,bm,cm,w;

 if(N == 0) return(0.0);

 if(N == 1) {
   integral = (y[0] + y[1])*(x[1] - x[0])/2.0;
   return(integral);
 }

 if(N == 2) {
   integral = (y[0]+y[1])*(x[1]-x[0])/2.0 + (y[1]+y[2])*(x[2]-x[1])/2.0;
   return(integral);
 }


 for(i=0;i<N;i++) {
   j = imax(i,1);
   k = imin(i+1,N-1);

   if(i == 0) {
     c[0] = y[j+1]/((x[j+1]-x[j])*(x[j+1]-x[j-1])) -
	    y[j]/((x[j] - x[j-1])*(x[j+1]-x[j])) +
	    y[j-1]/((x[j]-x[j-1])*(x[j+1]-x[j-1]));
     b[0] = (y[j]-y[j-1])/(x[j]-x[j-1]) - (x[j]+x[j-1])*c[0];
     a[0] = y[j-1] - x[j-1]*(y[j]-y[j-1])/(x[j]-x[j-1]) + x[j]*x[j-1]*c[0];
   } else {
     c[0] = c[1];
     b[0] = b[1];
     a[0] = a[1];
   }

   c[1] = y[k+1]/((x[k+1]-x[k])*(x[k+1]-x[k-1])) -
	  y[k]/((x[k] - x[k-1])*(x[k+1]-x[k])) +
	  y[k-1]/((x[k]-x[k-1])*(x[k+1]-x[k-1]));
   b[1] = (y[k]-y[k-1])/(x[k]-x[k-1]) - (x[k]+x[k-1])*c[1];
   a[1] = y[k-1] - x[k-1]*(y[k]-y[k-1])/(x[k]-x[k-1]) + x[k]*x[k-1]*c[1];

   if(c[1] == 0.0) w = 0.0;
   else w = fabs(c[1])/(fabs(c[1]) + fabs(c[0]));
   am = w*a[0] + (1.0-w)*a[1];
   bm = w*b[0] + (1.0-w)*b[1];
   cm = w*c[0] + (1.0-w)*c[1];
   if(c[1] != 0.0) {
     if(fabs(cm/c[1]) < 1.0e-14) cm = 0.0;
   }
   integral += (am + bm*(x[i+1] + x[i])/2.0 +
	       cm*(x[i+1]*x[i+1] + x[i+1]*x[i] + x[i]*x[i])/3.0)*(x[i+1] - x[i]);
 }

 return(integral);
}


