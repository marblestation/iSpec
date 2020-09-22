#include "spectrum.h"
#include <math.h>
#include <errno.h>
double planck();
double expint();
double integ();
double trapez();

void fluxflx(model,wstart,wend,flxstart,flxend)
atmosphere *model;
double wstart,wend,*flxstart,*flxend;
{
 int ntau;
 double y[NTAU];
 extern int Ntau;

 y[0] = 0.0;
 for(ntau=1;ntau<Ntau;ntau++) {
    y[ntau] = 2.0*model->kappawstart[ntau]*model->tauref[ntau]*
    planck(wstart,model->T[ntau])*expint(model->tauwstart[ntau],2)/
    (0.4343*model->kapparef[ntau]);
 }
 if(wstart > 3645.0 && wstart < 3750.0)
   *flxstart = trapez(y,model->x,Ntau-1);
 else *flxstart = integ(y,model->x,Ntau-1);

 for(ntau=1;ntau<Ntau;ntau++) {
    y[ntau] = 2.0*model->kappawend[ntau]*model->tauref[ntau]*
    planck(wend,model->T[ntau])*expint(model->tauwend[ntau],2)/
    (0.4343*model->kapparef[ntau]);
 }
 if(wend > 3645.0 && wend < 3750.0)  *flxend = trapez(y,model->x,Ntau-1);
 else *flxend = integ(y,model->x,Ntau-1);


}


