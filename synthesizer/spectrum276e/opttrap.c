#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "spectrum.h"

void opttrap(y,x,N,z)
double y[],x[],z[];
int N;
{
 int i;
 double integral;

 integral = 0.0;
 z[0] = 0.0;

 z[1] = (y[0] + y[1])*(x[1] - x[0])/2.0;

 integral = z[1];

 for(i=2;i<=N;i++) {
   integral += (y[i]+y[i-1])*(x[i]-x[i-1])/2.0;
   z[i] = integral;
 }

}


