#include "spectrum.h"
#include <math.h>
#include <errno.h>
double planck();
double integ();
double intensity();
void intenint();
double trapez();

double intensity(model,wave,mu)
atmosphere *model;
double wave,mu;
{
 int ntau;
 double y[NTAU],Intensity;
 extern int Ntau;

 y[0] = 0.0;
 for(ntau=1;ntau<Ntau;ntau++) {
    y[ntau] = model->kappawave[ntau]*model->tauref[ntau]*
    planck(wave,model->T[ntau])*exp(-model->tauwave[ntau]/mu)/
    (0.4343*mu*model->kapparef[ntau]);
 }
 Intensity = integ(y,model->x,Ntau-1);
 return(Intensity);
}

void intenint(model,wstart,wend,intenstart,intenend,mu)
atmosphere *model;
double wstart,wend,mu;
double *intenstart,*intenend;
{
 int ntau;
 double y[NTAU];
 extern int Ntau;

 y[0] = 0.0;
 for(ntau=1;ntau<Ntau;ntau++) {
    y[ntau] = model->kappawstart[ntau]*model->tauref[ntau]*
    planck(wstart,model->T[ntau])*exp(-model->tauwstart[ntau]/mu)/
    (0.4343*mu*model->kapparef[ntau]);
 }
 if(wstart > 3645.0 && wstart < 3750.0)
   *intenstart = trapez(y,model->x,Ntau-1);
 else *intenstart = integ(y,model->x,Ntau-1);

 y[0] = 0.0;
 for(ntau=1;ntau<Ntau;ntau++) {
    y[ntau] = model->kappawend[ntau]*model->tauref[ntau]*
    planck(wend,model->T[ntau])*exp(-model->tauwend[ntau]/mu)/
    (0.4343*mu*model->kapparef[ntau]);
 }
 if(wend > 3645.0 && wend < 3750.0) *intenend = trapez(y,model->x,Ntau-1);
 else  *intenend = integ(y,model->x,Ntau-1);
}
