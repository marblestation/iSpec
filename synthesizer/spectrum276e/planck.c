#include <math.h>
#include <errno.h>

double planck(wave,temp)
double wave,temp;
{
  static double p = 1.19106e+27;
  double p1;
  p1 = pow(wave,5.0)*(exp(1.43879e+08/(wave*temp)) -1.0);
  return(p/p1);
}


