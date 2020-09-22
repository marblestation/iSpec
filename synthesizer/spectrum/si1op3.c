#include <math.h>
#include <errno.h>
#include "spectrum.h"
double partfn(double code,double T,double Ne);

double si1op(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  double hckt,A,H,waveno,x,u;

  u = partfn(14.0,model->T[j],model->Ne[j]);
  waveno = 3.335641e-11*nu;
  hckt = 1.438832/model->T[j];

  H = 1.0e-30;

  /* 5948.497 A */
  if(waveno >= 16810.969) {
    x = 0.0;
    A = x*9.0*exp(-49128.131*hckt)*model->stim[j];
    H += A;
  }

  /* 5625.043 A */
  if(waveno >= 17777.641) {
    x = 18.0e-18*pow(17777.641/waveno,3.0);
    A = x*15.0*exp(-48161.459*hckt)*model->stim[j];
    H += A;
  }

  /* 5378.946 A */
  if(waveno >= 18587.546) {
    x = 0.0;
    A = x*5.0*exp(-47351.554*hckt)*model->stim[j];
    H += A;
  }

  /* 5360.482 A */
  if(waveno >= 18655.039) {
    x = 0.0;
    A = x*3.0*exp(-47284.061*hckt)*model->stim[j];
    H += A;
  }

  /* 4008.463 A */
  if(waveno >= 24947.216) {
    x = 4.09e-18*pow(24947.216/waveno,2.0);
    A = x*3.0*exp(-40991.884*hckt)*model->stim[j];
    H += A;
  }

  /* 3834.476 A */
  if(waveno >= 26079.180) {
    x = 1.25e-18*pow(26079.180/waveno,2.0);
    A = x*9.0*exp(-39859.920*hckt)*model->stim[j];
    H += A;
  }

  /* 1985.972 A */
  if(waveno >= 50353.180) {
    x = 46.0e-18*pow(50353.180/waveno,5.0);
    x /= 3.0;
    A = x*1.0*exp(-15394.370*hckt)*model->stim[j];
    H += A;
  }

  /* 1974.699 */
  if(waveno >= 50640.630) {
    A *= 2.0 ;
    H += A;
  }

  /* 1682.123 A */
  if(waveno >= 59448.700) {
    x = 35.0e-18*pow(59488.700/waveno,3.0);
    x /= 3.0;
    A = x*5.0*exp(- 6298.850*hckt)*model->stim[j];
    H += A;
  }

  /* 1674.028 A */
  if(waveno >= 59736.150) {
    A *= 2.0 ;
    H += A;
  }

  /* 1576.131 A */
  if(waveno >= 63446.510) {
    x = 18.0e-18*pow(63446.510/waveno,3.0);
    A = x*15.0*exp(-45303.310*hckt)*model->stim[j];
    H += A;
  }

  /* 1526.149 A */
  if(waveno >= 65524.393) {
    x = 37.0e-18;
    if(waveno >= 74000.0) x = x*pow(74000.0/waveno,5.0);
    x /= 3.0;
    A = x*5.0*exp(-223.157*hckt)*model->stim[j];
    H += A;
  }

  /* 1522.755 A */
  if(waveno >= 65670.435) {
    A = x*3.0*exp(-77.115*hckt)*model->stim[j];
    H += A;
  }

  /* 1520.969 A */
  if(waveno >= 65747.550) {
    A = x*1.0*1.0*model->stim[j];
    H += A;
  }

  /* 1519.483 A */
  if(waveno >= 65811.843) {
    x *= 2.0;
    A = x*5.0*exp(-223.157*hckt)*model->stim[j];
    H += A;
  }

  /* 1516.119 A */
  if(waveno >= 65957.885) {
    A = x*3.0*exp(-77.115*hckt)*model->stim[j];
    H += A;
  }

  /* 1514.348 A */
  if(waveno >= 66035.000) {
    A = x*1.0*1.0*model->stim[j];
    H += A;
  }

  /* 1325.842 A */
  if(waveno >= 75423.767) {
    x = 15.0e-18*pow(75423.767/waveno,3.0);
    A = x*5.0*exp(-33326.053*hckt)*model->stim[j];
    H += A;
  }

  return(H/u);
}


