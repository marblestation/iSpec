#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double integ();
double planck();
double expint();

double depth(model,wave,Flux)
atmosphere *model;
double wave,Flux;
{
  int i;
  double cd[NTAU],Depth,e1,e2,e3,e4;
  extern int Ntau;
  extern int flagp;
  extern int flagC;
  FILE *out;

  for(i=1;i<Ntau;i++) {
      cd[i] = (4.6052/Flux)*(model->tauref[i]*
      model->kappawave[i]/model->kapparef[i])*planck(wave,model->T[i])*
      (expint(model->tauwave[i],2) - (1.0 + model->kapnu[i]/
      model->kappawave[i])*expint(model->tauwave[i] + model->taunu[i],2));
  }
  cd[0] = 0.0;
  if(flagC == 1) {
    out = fopen("cd.out","a");
    for(i=0;i<Ntau;i++) fprintf(out,"%d %e %e\n",i+1,model->tauref[i],cd[i]);
    fclose(out);
  }
  Depth = integ(cd,model->x,Ntau-1);
  if(Depth < 0.0 && flagp == 0) Depth = 0.0;
  if(Depth > 1.0) Depth = 1.0;
  return(Depth);
}


