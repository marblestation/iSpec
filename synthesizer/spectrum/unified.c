#include <math.h>
#include <errno.h>
#include "spectrum.h"
double partfn(double code,double T,double Ne);
double Hprofl();

double unified(n,m,ntau,El,Eh,model,wave)
int n,m,ntau;
double El,Eh,wave;
atmosphere *model;
{
  double kappa = 0;
  double T,kT,U;

  kT = model->kT[ntau];
  U = model->U[ntau];
  kappa = model->NHI[ntau]*2*n*n*exp(-El/kT)*(1.0 - exp((El-Eh)/kT))/U;
  kappa *= Hprofl(wave,n,m,ntau,model);
  return(kappa);
}

