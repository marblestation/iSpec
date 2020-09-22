#include "spectrum.h"
#include <math.h>
#include <errno.h>
double opacity();
double integ();
void optdepth();

void tauwave(model,wave)
double wave;
atmosphere *model;
{
 int ntau;
 extern int Ntau;
 double y[NTAU];

 for(ntau=1;ntau<Ntau;ntau++) {
   model->kappawave[ntau] = opacity(model,wave,ntau);
   y[ntau] = model->kappawave[ntau]/model->kapparef[ntau];
 }
 model->kappawave[0] = model->kappawave[1];
 y[0] = y[1];
 optdepth(y,model->tauref,Ntau-1,model->tauwave);
}


