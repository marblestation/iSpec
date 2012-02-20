#include "spectrum.h"
#include <math.h>
#include <errno.h>
double planck();
double expint();
double integ();

double flux(model,wave)
atmosphere *model;
double wave;
{
 int ntau;
 double y[NTAU],Flux;
 extern int Ntau;

 y[0] = 0.0;
 for(ntau=1;ntau<Ntau;ntau++) {
    y[ntau] = 2.0*model->kappawave[ntau]*model->tauref[ntau]*
    planck(wave,model->T[ntau])*expint(model->tauwave[ntau],2)/
    (0.4343*model->kapparef[ntau]);
 }
 Flux = integ(y,model->x,Ntau-1);
 return(Flux);
}


