#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "spectrum.h"
double seaton();
double n1op();
double o1op();
double mg2op();
double ca2op();
double si2op();
double partfn(double code,double T,double Ne);
double lukeop();
int imax();
int imin();

double n1op(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  static float c1130[NTAU];
  static float c1020[NTAU];
  static int flag = 0;
  int i;
  float tkev,x1130,x1020,x853;
  double u;
  extern int Ntau;

  u = partfn(7.0,model->T[j],model->Ne[j]);

  if(flag == 0) {
    flag = 1;
    for(i=0;i<Ntau;i++) {
      tkev = 8.617e-05*model->T[i];
      c1130[i] = 6.0*exp(-3.575/tkev);
      c1020[i] = 10.0*exp(-2.384/tkev);
    }
  }
  x1130 = x1020 = x853 = 0.0;
  if(nu >= 3.517915e+15) x853 = seaton(nu,3.517915e+15,1.142e-17,2.0,4.29);
  if(nu >= 2.941534e+15) x1020 = seaton(nu,2.941534e+15,4.41e-18,1.5,3.85);
  if(nu >= 2.653317e+15) x1130 = seaton(nu,2.653317e+15,4.2e-18,1.5,4.34);

  return((4.0*x853 + c1020[j]*x1020 + c1130[j]*x1130)/u);
}

double o1op(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  float x911 = 0.0;
  double u;

  u = partfn(8.0,model->T[j],model->Ne[j]);

  if(nu >= 3.28805e+15) x911 = seaton(nu,3.28805e+15,2.94e-18,1.0,2.66);
  return(9.0*x911/u);
}

double mg2op(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  static float c1169[NTAU];
  static int flag = 0;
  double tkev,x1169,x824,u;
  int i;
  extern int Ntau;

  u = partfn(12.1,model->T[j],model->Ne[j]);

  if(flag == 0) {
    flag = 1;
    for(i=0;i<Ntau;i++) {
      tkev = 8.617e-05*model->T[i];
      c1169[i] = 6.0*exp(-4.43/tkev);
    }
  }
  x1169 = x824 = 0.0;
  if(nu >= 3.635492e+15) x824 = seaton(nu,3.635492e+15,1.40e-19,4.0,6.7);
  if(nu >= 2.564306e+15) x1169 = 5.11e-19*pow(2.564306e+15/nu,3.0);
  return((2.0*x824 + c1169[j]*x1169)/u);
}

double ca2op(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  static float c1218[NTAU],c1420[NTAU];
  double x1420,x1218,x1044,tkev,u;
  int i;
  static int flag = 0;
  extern int Ntau;

  u = partfn(20.1,model->T[j],model->Ne[j]);

  if(flag == 0) {
    flag = 1;
    for(i=0;i<Ntau;i++) {
      tkev = 8.617e-05*model->T[i];
      c1218[i] = 10.0*exp(-1.697/tkev);
      c1420[i] = 6.0*exp(-3.142/tkev);
    }
  }
  x1420 = x1218 = x1044 = 0.0;
  if(nu >= 2.870454e+15) x1044 = 5.4e-20*pow(2.870454e+15/nu,3.0);
  if(nu >= 2.460127e+15) x1218 = 1.64e-17*sqrt(2.460127e+15/nu);
  if(nu >= 2.110779e+15) x1420 = seaton(nu,2.110779e+15,4.13e-18,3.0,0.69);
  return((2.0*x1044 + c1218[j]*x1218 + c1420[j]*x1420)/u);
}

double si2op(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  static float peach[14][6] =
  {{-43.8941,-43.8941,-43.8941,-43.8941,-43.8941,-43.8941},
   {-42.2444,-42.2444,-42.2444,-42.2444,-42.2444,-42.2444},
   {-40.6054,-40.6054,-40.6054,-40.6054,-40.6054,-40.6054},
   {-54.2389,-52.2906,-50.8799,-49.8033,-48.9485,-48.2490},
   {-50.4108,-48.4892,-47.1090,-46.0672,-45.2510,-44.5933},
   {-52.0936,-50.0741,-48.5999,-47.4676,-46.5649,-45.8246},
   {-51.9548,-49.9371,-48.4647,-47.3340,-46.4333,-45.6947},
   {-54.2407,-51.7319,-49.9178,-48.5395,-47.4529,-46.5709},
   {-52.7355,-50.2218,-48.4059,-47.0267,-45.9402,-45.0592},
   {-53.5387,-50.9189,-49.0200,-47.5750,-46.4341,-45.5082},
   {-53.2417,-50.6234,-48.7252,-47.2810,-46.1410,-45.2153},
   {-53.5097,-50.8535,-48.9263,-47.4586,-46.2994,-45.3581},
   {-54.0561,-51.2365,-49.1980,-47.6497,-46.4302,-45.4414},
   {-53.8469,-51.0256,-48.9860,-47.4368,-46.2162,-45.2266}};
  static double freqsi[7] = {4.996541e+15,3.9466738e+15,1.5736321e+15,
   1.5171539e+15,9.237894e+14,8.3825004e+14,7.6869872e+14};
  static double flog[9] = {36.32984,36.14752,35.91165,34.99216,34.95561,
   34.45951,34.36234,34.27572,34.20161};
  static double tlg[6] = {9.21034,9.39266,9.54681,9.68034,9.79813,9.90349};
  int k,n,it;
  int nt;
  double dt,u;
  double d,d1,si2op,x[6];

  u = partfn(14.1,model->T[j],model->Ne[j]);

  n = (int)imax((int)imin(5,(int)(model->T[j]/2000.0)-4),1);
  nt = n;
  dt = (log(model->T[j])-tlg[n-1])/(tlg[n]-tlg[n-1]);
  for(n=1;n<8;n++) {
    if(nu > freqsi[n-1]) break;
  }
  d = (log(nu) - flog[n-1])/(flog[n]-flog[n-1]);
  if(n > 2) n = 2*n-2;
  if(n == 14) n = 13;
  d1 = 1.0-d;
  for(it=0;it<6;it++) x[it] = peach[n][it]*d + peach[n-1][it]*d1;
  n = nt;
  si2op = 6.0*exp(x[n-1]*(1.0-dt) + x[n]*dt);
  return(si2op/u);

}

double lukeop(nu,model,j)
double nu;
atmosphere *model;
int j;
{
  return((n1op(nu,model,j)*model->NNI[j] + o1op(nu,model,j)*model->NOI[j] +
   mg2op(nu,model,j)*model->NMgII[j] + si2op(nu,model,j)*model->NSiII[j] +
   ca2op(nu,model,j)*model->NCaII[j])*model->stim[j]);
}
