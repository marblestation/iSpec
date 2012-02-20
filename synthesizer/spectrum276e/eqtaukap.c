#include <math.h>
#include <stdio.h>
#include <errno.h>
#include "spectrum.h"
double integ();
double voigt();
void optdepth();
double autoion(double wave, double center, double gammar, double ashore, double bshore);

void eqtaukap(wave,model,line)
atmosphere *model;
linedata *line;
double wave;
{
  int i,j;
  double v,y[NTAU];
  extern int Ntau;

  for(i=1;i<Ntau;i++) {
    model->kapnu[i] = 0.0;
    if(line[0].ai == 1) { 
      model->kapnu[i] += line[0].capnu[i]*
              autoion(wave,line[0].wave,line[0].gammar,line[0].gammas,
                      line[0].gammaw);
    } else {
      v = 2.997929e+10*(line[0].wave - wave)/(line[0].wave*line[0].dopp[i]);
      model->kapnu[i] += line[0].capnu[i]*voigt(line[0].a[i],v);
    }
  }
  model->kapnu[0] = model->kapnu[1];

  for(i=1;i<Ntau;i++) y[i] = model->tauref[i]*model->kapnu[i]/
			     (0.4343*model->kapparef[i]);
  y[0] = 0;
  optdepth(y,model->x,Ntau-1,model->taunu);
}
