#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double integ();
double voigt();
void balmer();
void lyman();
void paschen();
void brackett();
void pfund();
void humphreys();
double stronglines();
void optdepth();
void opttrap();
double helium();
double autoion(double wave, double center, double gammar, double ashore, double bshore);

void taukap(wave,model,atom,line,N,strgln,V,He,POP)
atmosphere *model;
atominfo *atom;
linedata *line,*strgln;
pfunc *V;
Helium *He;
population *POP;
double wave;
int N;
{
  int i,j;
  double v,y[NTAU];
  double hkap[NTAU];
  double ckm = 2.997929e+05;
  double delwave;
  extern int Ntau;
  extern int flagc;
  extern float *velgrad;
  extern double mu;

  for(i=0;i<Ntau;i++) {
    hkap[i] = 0.0;
    model->kapnu[i] = 0.0;
  }
  if(flagc == 0) {
   if(model->teff > 3750.0) {
    lyman(hkap,model,wave);
    balmer(hkap,model,wave);
    paschen(hkap,model,wave);
    for(i=1;i<Ntau;i++) model->kapnu[i] = hkap[i];
    if(wave > 14779.0) {
      for(i=0;i<Ntau;i++) hkap[i] = 0.0;
      brackett(hkap,model,wave);
      for(i=1;i<Ntau;i++) model->kapnu[i] += hkap[i];
      if(wave > 23291.0) {
        for(i=0;i<Ntau;i++) hkap[i] = 0.0;
        pfund(hkap,model,wave);
        for(i=1;i<Ntau;i++) model->kapnu[i] += hkap[i];
	if(wave > 33869.0) {
          for(i=0;i<Ntau;i++) hkap[i] = 0.0;
          humphreys(hkap,model,wave);
          for(i=1;i<Ntau;i++) model->kapnu[i] += hkap[i];
	}
      }
    }
   }
    for(i=1;i<Ntau;i++) {
      model->kapnu[i] += helium(model,i,wave,V,He);
      model->kapnu[i] += stronglines(model,atom,i,strgln,wave,V,POP);
      delwave = wave*mu*velgrad[i]/ckm;
      for(j=0;j<N;j++) {
	if(line[j].ai == 1) 
          model->kapnu[i] += line[j].capnu[i]*
              autoion(wave,line[j].wave,line[j].gammar,line[j].gammas,
                      line[j].gammaw);
	else{
        v = 2.997929e+10*(line[j].wave + delwave - wave)/
	    (line[j].wave*line[j].dopp[i]);
        model->kapnu[i] += line[j].capnu[i]*voigt((double)line[j].a[i],v);
	}
      }
    }
  }
  model->kapnu[0] = model->kapnu[1];

  for(i=1;i<Ntau;i++) y[i] = model->tauref[i]*model->kapnu[i]/
			     (0.4343*model->kapparef[i]);
  y[0] = 0;
  if(wave > 3645.0 && wave < 3750.0) opttrap(y,model->x,Ntau-1,model->taunu);
  else optdepth(y,model->x,Ntau-1,model->taunu);
}


