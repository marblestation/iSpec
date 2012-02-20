#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "spectrum.h"
double AN();
double AHeIIN();
double AHminus2();
double Hmff();
double he1op_av(double nu,double T,double Ne);
double *dvector();
double dmin();
double poly();
double cutoff();
double pcutoff();
double lcutoff();
float *vector();
float **matrix();
void free_vector();
void free_matrix();
void free_dvector();
double partfn(double code,double T,double Ne);
double Fff();
double FH2();
void lin2();
double coolop();
double lukeop();
double c = 2.99792458e+18;
double h = 6.62618e-27;
double k = 8.617084e-05;
double kh = 1.38054e-16;


double opacity(model,lambda,ntau)
atmosphere *model;
double lambda;
int ntau;
{
  double nu,khbf,khff,khmbf,khmff,kh2,khray,sig,*nh,*nheii,T,NHeI1,kheray,A,B,C;
  double khmbf2,khmff2,Ne;
  double khemff,lam,khebf,kff,kheIff,kheIIff,kh2ray,kheIIbf,lnai,HeIn,UHeI;
  double stim,klotemp,kluke,kall,lam2,lam4;
  static double Cutoff[NTAU];
  static double pCutoff[NTAU];
  static double lCutoff[NTAU];
  static int flag = 0;
  static double cut;
  double cutmax,cutmin;
  extern int Ntau;
  extern int flagO;
  extern FILE *opout;
  int NM;
  int i,n;

  nh = dvector(1,8);
  nheii = dvector(1,9);
  nu = c/lambda;
  T = model->T[ntau];
  Ne = model->Ne[ntau];
  model->stim[ntau] = stim = 1.0 - exp(-h*nu/(kh*T));

  if(flag == 0) {
    for(i=0;i<Ntau;i++) {
      Cutoff[i] = cutoff(model->T[i],model->Ne[i],&NM);
      pCutoff[i] = pcutoff(model->T[i],model->Ne[i],&NM);
      lCutoff[i] = lcutoff(model->T[i],model->Ne[i],&NM);
    }
/*  cutmin = cutmax = Cutoff[0];
    for(i=1;i<Ntau;i++) {
      if(cutmin > Cutoff[i]) cutmin = Cutoff[i];
      if(cutmax < Cutoff[i]) cutmax = Cutoff[i];
    }
    cut = (cutmax+cutmin)/2.0;  */
    flag = 1;
  }

/* Calculate Hydrogen Bound-Free Opacity */

  for(n=1;n<9;n++) nh[n] = model->NHI[ntau]*(2.0*n*n)*
			   exp(-13.595*(1.0 - 1.0/(n*n))/(k*T))/partfn(1.0,T,Ne);

  khbf = 0.0;
  for(n=1;n<9;n++) khbf += nh[n]*AN(n,nu,Cutoff[ntau],pCutoff[ntau],
                                                      lCutoff[ntau]);
  khbf *= stim;
  khbf += (2.815e+29*model->NHI[ntau]/(pow(nu,3.0)*partfn(1.0,T,Ne)))*
	  stim*(k*T/13.595)*
	  (exp(-13.595*(1.0-1.0/81.0)/(k*T)) - exp(-13.595/(k*T)));



/* Calculate Hydrogen Free-Free Opacity */

  kff = model->Ne[ntau]*Fff(nu,T,1.0)*stim;
  khff = model->Np[ntau]*kff;

/* Calculate Helium I Free-Free Opacity */
  kheIff = model->NHeII[ntau]*kff;

/* Calculate Helium II Free-Free Opacity */
  kff = model->Ne[ntau]*Fff(nu,T,2.0)*stim;
  kheIIff = model->NHeIII[ntau]*kff;

/* Calculate H minus Bound-Free Opacity */

/*  khmbf = model->NHminus[ntau]*AHminus(nu)*stim;  */
  khmbf2 = model->NHminus[ntau]*AHminus2(lambda)*stim;

/* Calculate H minus Free-Free Opacity */

/*  khmff = model->N1[ntau]*model->Ne[ntau]*(1.3727e-25 + (4.3748e-10 -
	  2.5993e-7/T)/nu)/nu;  */
  khmff2 = model->NHI[ntau]*model->Ne[ntau]*Hmff(T,lambda);

/* Calculate H2 plus Opacity */

  kh2 =  model->N1[ntau]*model->Np[ntau]*FH2(nu,T)*stim;

/* Calculate Hydrogen Rayleigh scattering */
  lam = c/(double)dmin(nu,2.463e+15);
  khray = model->N1[ntau]*(5.799e-13/pow(lam,4.0) +
	  1.422e-6/pow(lam,6.0) + 2.784/pow(lam,8.0));

/* Calculate H2 Rayleigh Scattering */

  lam = c/(double)dmin(nu,2.922e+15);
  kh2ray = model->NH2[ntau]*(8.14e-13/pow(lam,4.0) + 1.28e-06/pow(lam,6.0) +
	   1.61/pow(lam,8.0));

/* Calculate Helium Bound-Free Opacity */


  khebf = 0.0;
  UHeI = partfn(2.0,T,Ne);

  khebf = model->NHeI[ntau]*he1op_av(nu,T,Ne)/UHeI*stim;

/* Calculate Helium II bound-free opacity */

  for(n=1;n<10;n++) nheii[n] = model->NHeII[ntau]*(2.0*(double)(n*n))*
			   exp(-54.403*(1.0 - 1.0/(n*n))/(k*T))/partfn(2.1,T,Ne);

  kheIIbf = 0.0;
  for(n=1;n<10;n++) kheIIbf += nheii[n]*AHeIIN(n,nu);
  kheIIbf *= stim;
  kheIIbf += (4.504e+30*model->NHeII[ntau]/(pow(nu,3.0)*partfn(2.1,T,Ne)))*
	  stim*(k*T/54.403)*
	  (exp(-54.403*(1.0-1.0/100.0)/(k*T)) - exp(-54.403/(k*T)));

/* Calculate Helium Rayleigh scattering */

  NHeI1 = model->NHeI[ntau]/partfn(2.0,T,Ne);
  lam = c/(double)dmin(nu,5.15e+15);
  lam2 = lam*lam;
  kheray = 5.484e-14/lam2/lam2*
    pow(1.0 + (2.44e+5 + 5.94e+10/(lam2-2.90e+5))/lam2,2.0);
  kheray *= NHeI1;
  /*  kheray = NHeI1*(5.484e-14/pow(lam,4.0))*pow(1.0 + 2.44e+5/(lam*lam) +
      5.94e-10/((lam*lam)*((lam*lam) - 2.90e+5)),2.0); */

/* Calculate HeI minus free-free opacity */

  A = 3.397e-46 + (-5.216e-31 + 7.039e-15/nu)/nu;
  B = -4.116e-42 + (1.067e-26 + 8.135e-11/nu)/nu;
  C = 5.081e-37 + (-8.724e-23 - 5.659e-8/nu)/nu;
  khemff = NHeI1*model->Ne[ntau]*(A*T + B + C/T);

/* Calculate Electron Scattering Opacity */

  sig = 6.653e-25*model->Ne[ntau];

/* Low temperature Opacity due to CI, SiI, MgI and AlI  */

  klotemp = coolop(nu,model,ntau);

/* Intermediate temperature Opacity due to NI, OI, MgII, CaII and SiII */

  kluke = lukeop(nu,model,ntau);

  free_dvector(nh,1,8);
  free_dvector(nheii,1,9);

  if(flagO == 1) 
    fprintf(opout,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            khbf,khff,khmbf2,khmff2,kh2,khray,khebf,kheIff,kheIIff,kheray,
            khemff,kh2ray,kheIIbf,klotemp,kluke,sig);

  return(khbf+khff+khmbf2+khmff2+kh2+khray+khebf+kheIff+kheIIff+kheray+khemff+
	 kh2ray+kheIIbf+klotemp+kluke+sig);
}

double AN(n,nu,cutoff,pcutoff,lcutoff)
int n;
double nu,cutoff,pcutoff,lcutoff;
{
  double An;
  static double A[] = {0,0.9916,1.105,1.101,1.101,1.102,1.0986,1,1};
  static double B[] = {0,2.719e+13,-2.375e+14,-9.863e+13,-5.765e+13,
		       -3.909e+13,-2.704e+13,0,0};
  static double C[] = {0,-2.268e+30,4.077e+28,1.035e+28,4.593e+27,2.371e+27,
		       1.229e+27,0,0};
  static double Qn = 0.0;
   double *wave0,wav;

/*   wave0 = w0-3;  */
   wav = c/nu;

  if(n == 1 && wav > lcutoff) return(0.0);
  if(n == 2 && wav > cutoff) return(0.0);
/*if(n == 2 && wav > 3670.0) return(0.0);  */
/*if(n == 3 && wav > 8203.601) return(0.0); */
  if(n == 3 && wav > pcutoff) return(0.0);
  if(n == 4 && wav > 14810.0) return(0.0);
  if(n == 5 && wav > 23295.0) return(0.0);
  if(n == 6 && wav > 33874.0) return(0.0);
  if(n >= 7 && nu < 3.28805e+15/((double)(n*n))) return(0.0);


  An = (2.815e+29/(pow(nu,3.0)*pow((double)n,5.0)))*(A[n] + (B[n] + C[n]/nu)/nu);
  return(An);
}

double AHeIIN(n,nu)
int n;
double nu;
{
  double An;
  static double A[] = {0.0,0.9916,1.105,1.101,1.101,1.102,1.0986,1.0,1.0};
  static double B[] = {0.0,2.719e+13,-2.375e+14,-9.863e+13,-5.765e+13,
		       -3.909e+13,-2.704e+13,0,0.0,0.0};
  static double C[] = {0.0,-2.268e+30,4.077e+28,1.035e+28,4.593e+27,2.371e+27,
		       1.229e+27,0,0.0,0.0};

  if(nu < 1.31522e+16/(double)(n*n)) return(0.0);

  An = (1.126e+30/(pow(nu,3.0)*pow((double)n,5.0)))*(A[n] + (4.0*B[n] + 16.0*C[n]/nu)/nu);
  return(An);
}

double Fff(nu,T,Z)
double nu,T,Z;
{
  static float x1n[] = {-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0};
  static float x2n[] = {-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5};
  static float yb[11][12] =
		       {{5.53,4.91,4.29,3.64,3.00,2.41,1.87,1.33,0.90,0.55,0.33,0.19},
			{5.49,4.87,4.25,3.61,2.98,2.41,1.89,1.39,0.95,0.58,0.36,0.21},
			{5.46,4.84,4.22,3.59,2.97,2.41,1.91,1.44,1.00,0.62,0.39,0.24},
			{5.43,4.80,4.18,3.56,2.95,2.41,1.93,1.49,1.08,0.70,0.46,0.28},
			{5.40,4.77,4.15,3.54,2.94,2.41,1.95,1.55,1.17,0.85,0.59,0.38},
			{5.25,4.63,4.02,3.41,2.81,2.32,1.90,1.56,1.30,1.01,0.76,0.53},
			{5.00,4.40,3.80,3.22,2.65,2.19,1.80,1.51,1.32,1.15,0.97,0.76},
			{4.69,4.13,3.57,2.97,2.44,2.02,1.68,1.42,1.30,1.18,1.09,0.96},
			{4.48,3.87,3.27,2.70,2.21,1.84,1.52,1.33,1.20,1.15,1.13,1.08},
			{4.16,3.52,2.98,2.45,2.01,1.67,1.41,1.25,1.15,1.11,1.10,1.09},
			{3.85,3.27,2.70,2.20,1.81,1.50,1.30,1.17,1.11,1.08,1.08,1.09}};
/* Errors corrected in column 9 Jan 25, 1995.  See ATLAS7v code */

  float x1,x2,y,dy;
  double F;
  int i,j;
  float **ya;
  void polynn2();

  ya = matrix(0,10,0,11);

  for(i=0;i<11;i++) {
     for(j=0;j<12;j++) ya[i][j] = yb[i][j];
  }


  x2 = log10(h*nu/(kh*T));
  x1 = log10(3.28805e+15*Z*Z*h/(kh*T));

  if(x1 < -3.0) x1 = -3.0;
  if(x1 > 2.0) x1 = 2.0;
  if(x2 < -4.0) x2 = -4.0;
  if(x2 > 1.5) x2 = 1.5;

  polynn2(x1n,x2n,ya,11,12,x1,x2,&y);
  F = (3.6919e+08*Z*Z*y)/(pow(nu,3.0)*sqrt(T));
  free_matrix(ya,0,10,0,11);
  return(F);
}

/* double AHminus(nu)
double nu;
{
   double anu;

   anu = 0.0;
   if(nu >= 2.111e+14) anu = 6.801e-20 + (5.358e-3 + (1.481e+13 +
			     (-5.519e+27 + 4.808e+41/nu)/nu)/nu)/nu;
   if(nu < 2.111e+14 && nu >= 1.8259e+14) anu = 3.695e-16 + (-1.251e-1 +
					       1.052e+13/nu)/nu;
   return(anu);
}  */

double FH2(nu,T)
double nu,T;
{
   double Es,F;

   Es = -7.342e-3 + (-2.409e-15 + (1.028e-30 + (-4.230e-46 + (1.244e-61 -
	1.351e-77*nu)*nu)*nu)*nu)*nu;
   F = exp(-Es/(k*T) - 3.0233e+03 + (3.7797e+02 + (-1.82496e+01 +
       (3.9207e-01 - 3.1672e-03*log(nu))*log(nu))*log(nu))*log(nu));
   return(F);
}

double AHminus2(lambda)
double lambda;
{
/* The following polynomial fits Wishart's (1979) data to better than
0.2 percent, see Gray, D.F. Obs and Analy of Stellar Phot  */
/* Threshold added June 28, 2006 */

  double coef[] = {1.99654,-1.18267e-5,2.64243e-6,-4.40524e-10,
  3.23992e-14,-1.39568e-18,2.78701e-23};
  double hmff;

  if(lambda > 16300.0) return(0.0);

  hmff = 1.0e-18*poly(lambda,6,coef);

  if(hmff < 0.0) return(0.0);
  else return(hmff);
}

double Hmff(T,lambda)
double T,lambda;
{
/* The following expression is a fit to the data of Bell and
   Berrington 1987 */
  double loglam,f0,f1,f2,logth;
  double k = 1.380622e-16;
  double F0[] = {-2.2763,-1.6850,0.76661,-0.0533464};
  double F1[] = {15.2827,-9.2846,1.99381,-0.142631};
  double F2[] = {-197.789,190.266,-67.9775,10.6913,-0.625151};

  loglam = log10(lambda);
  logth = log10(5040.0/T);
  f0 = poly(loglam,3,F0);
  f1 = poly(loglam,3,F1);
  f2 = poly(loglam,4,F2);

  return(1.0e-26*pow(10.0,f0+f1*logth+f2*logth*logth)*k*T);
}















