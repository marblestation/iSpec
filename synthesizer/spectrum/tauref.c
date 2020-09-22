#include "spectrum.h"
#include <math.h>
#include <errno.h>
double opacity();
double integ();
void optdepth();

void tauref(model,wavref)
double wavref;
atmosphere *model;
{
 int ntau;
 double y[NTAU];
 extern int Ntau;

 for(ntau=0;ntau<Ntau;ntau++) {
   model->kapparef[ntau] = opacity(model,wavref,ntau);
   y[ntau] = model->kapparef[ntau]/model->rho[ntau];
 }
 y[0] = y[1];
 model->kapparef[0] = model->kapparef[1];
 model->mass[0] = 0.0;
 optdepth(y,model->mass,Ntau-1,model->tauref);

 for(ntau=1;ntau<Ntau;ntau++) model->x[ntau] = log10(model->tauref[ntau]);
 model->x[0] = model->x[1] - model->mass[1]*
	      (model->x[2] - model->x[1])/(model->mass[2] - model->mass[1]);
}
