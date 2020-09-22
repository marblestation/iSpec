#include <math.h>
#include "spectrum.h"
#include <errno.h>
double partfn(double code,double T,double Ne);

double al1op(nu,model,j)
double nu;
atmosphere *model;
int j;
{

  double hckt,A,H,waveno,x,u;

  u = partfn(13.0,model->T[j],model->Ne[j]);
  waveno = 3.335641e-11*nu;
  hckt = 1.438832/model->T[j];

  H = 1.0e-30;

  /* 14369.90 A */
  if(waveno >= 6958.993) {
    x = 0.0;
    A = x*14.0*exp(-41319.377*hckt)*model->stim[j];
    H += A;
  }

  /* 12496.15 A */
  if(waveno >= 8002.467) {
    x = 50.0e-18*pow(8002.467/waveno,3.0);
    A = x*6.0*exp(-40275.903*hckt)*model->stim[j];
    H += A;
  }

  /* 10699.50 A */
  if(waveno >= 9346.231) {
    x = 50.0e-18*pow(9346.231/waveno,3.0);
    A = x*10.0*exp(-38932.139*hckt)*model->stim[j];
    H += A;
  }

  /* 9443.801 A */
  if(waveno >= 10588.957) {
    x = 56.7e-18*pow(10588.957/waveno,1.9);
    A = x*2.0*exp(-37689.413*hckt)*model->stim[j];
    H += A;
  }

  /* 6528.264 A */
  if(waveno >= 15318.007) {
    x = 14.5e-18*15318.007/waveno;
    A = x*6.0*exp(-32960.363*hckt)*model->stim[j];
    H += A;
  }

  /* 6312.283 A */
  if(waveno >= 15842.129) {
    x = 47.0e-18*pow(15842.129/waveno,1.83);
    A = x*10.0*exp(-32436.241*hckt)*model->stim[j];
    H += A;
  }

  /* 4360.982 */
  if(waveno >= 22930.614) {
    x = 10.0e-18*pow(22930.614/waveno,2.0);
    A = x*2.0*exp(-25347.756*hckt)*model->stim[j];
    H += A;
  }

  /* 2076.140 A */
  if(waveno >= 48166.309) {
    x = 65.0e-18*pow(48166.309/waveno,5.0);
    A = x*4.0*exp(-112.061*hckt)*model->stim[j];
    H += A;
  }

  /* 2071.321 A */
  if(waveno >= 48278.370) {
    A = x*2.0*1.0*model->stim[j];
    H += A;
  }

  /* 1788.804 A */
    if(waveno >= 55903.260) {
    x = 10.0e-18*pow(55903.260/waveno,2.0);
    A = x*12.0*exp(-29097.11*hckt)*model->stim[j];
    H += A;
  }

  return(H/u);
}



