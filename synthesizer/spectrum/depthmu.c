#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double integ();
double planck();
double depthmu();
double idepthmu();
double trapez();

double depthmu(model,wave,Intensity,mu)
atmosphere *model;
double wave,Intensity,mu;
{
  int i;
  double cd[NTAU],Depth,e1,e2,e3,e4;
  extern int Ntau;
  extern int flagp;

  for(i=1;i<Ntau;i++) {
      cd[i] = (2.30256/(mu*Intensity))*(model->tauref[i]*
      model->kappawave[i]/model->kapparef[i])*planck(wave,model->T[i])*
      (exp(-model->tauwave[i]/mu) - (1.0 + model->kapnu[i]/
      model->kappawave[i])*exp(-(model->tauwave[i] + model->taunu[i])/mu));
  }
  cd[0] = 0.0;
  Depth = integ(cd,model->x,Ntau-1);
  if(Depth < 0.0 && flagp == 0) Depth = 0.0;
  if(Depth > 1.0) Depth = 1.0;
  return(Depth);
}


double idepthmu(model,wave,wstart,wend,mu)
atmosphere *model;
double wave,mu,wstart,wend;
{
  int i;
  double cd[NTAU],Depth,e1,e2,e3,e4;
  extern int Ntau;
  extern int flagp;

  for(i=1;i<Ntau;i++) {
    model->kappawave[i] = model->kappawstart[i] +
      (model->kappawend[i] - model->kappawstart[i])*
      (wave - wstart)/(wend - wstart);
    model->tauwave[i] = model->tauwstart[i] +
      (model->tauwend[i] - model->tauwstart[i])*
      (wave - wstart)/(wend - wstart);
  }

  for(i=1;i<Ntau;i++) {
      cd[i] = (2.30256/mu)*(model->tauref[i]*
      model->kappawave[i]/model->kapparef[i])*planck(wave,model->T[i])*
      (exp(-model->tauwave[i]/mu) - (1.0 + model->kapnu[i]/
      model->kappawave[i])*exp(-(model->tauwave[i] + model->taunu[i])/mu));
  }
  cd[0] = 0.0;
  if(wave > 3645.0 && wave < 3750.0) Depth = trapez(cd,model->x,Ntau-1);
  else Depth = integ(cd,model->x,Ntau-1);
  if(Depth < 0.0 && flagp == 0) Depth = 0.0;
  return(Depth);
}


