#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spectrum.h"

double autoion(double wave, double center, double gammar, double ashore, double bshore)
{
  double nu,nu0,epsilon;

  nu = 2.99792458e+18/wave;
  nu0 = 2.99792458e+18/center;

  epsilon = 2.0*(nu - nu0)/gammar;

  return((gammar/2.0)*(ashore*epsilon+bshore)/(epsilon*epsilon+1.0));
}
