#include "spectrum.h"
#include <math.h>
#include <errno.h>
double opacity();
double integ();
void optdepth();
void opttrap();

void tauflx(model,wstart,wend)
double wstart,wend;
atmosphere *model;
{
 int ntau;
 extern int Ntau;
 double y[NTAU];
 static int flag = 1;

 for(ntau=1;ntau<Ntau;ntau++) {
   model->kappawstart[ntau] = opacity(model,wstart,ntau);
   y[ntau] = model->kappawstart[ntau]/model->kapparef[ntau];
 }
 model->kappawstart[0] = model->kappawstart[1];
 y[0] = y[1];
 if(wstart > 3645.0 && wstart < 3750.0)
   opttrap(y,model->tauref,Ntau-1,model->tauwstart);
 else optdepth(y,model->tauref,Ntau-1,model->tauwstart);

 for(ntau=1;ntau<Ntau;ntau++) {
   model->kappawend[ntau] = opacity(model,wend,ntau);
   y[ntau] = model->kappawend[ntau]/model->kapparef[ntau];
 }
 model->kappawend[0] = model->kappawend[1];
 y[0] = y[1];
 if(wend > 3645.0 && wend < 3750.0)
   opttrap(y,model->tauref,Ntau-1,model->tauwend);
 else optdepth(y,model->tauref,Ntau-1,model->tauwend);

}
