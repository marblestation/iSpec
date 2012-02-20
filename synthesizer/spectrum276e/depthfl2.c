#include <math.h>
#include <errno.h>
#include "spectrum.h"
double integ();
double trapez();
double planck();
double expint();

double depthflx(model,wave,wstart,wend)
atmosphere *model;
double wave,wstart,wend;
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
      cd[i] = 4.6052*(model->tauref[i]*
      model->kappawave[i]/model->kapparef[i])*planck(wave,model->T[i])*
      (expint(model->tauwave[i],2) - (1.0 + model->kapnu[i]/
      model->kappawave[i])*expint(model->tauwave[i] + model->taunu[i],2));
  }
  cd[0] = 0.0;
  if(wave > 3645.0 && wave < 3750.0) Depth = trapez(cd,model->x,Ntau-1);
  else Depth = integ(cd,model->x,Ntau-1);
  if(Depth < 0.0 && flagp == 0) Depth = 0.0;
  return(Depth);
}


