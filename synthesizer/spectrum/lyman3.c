#include <math.h>
#include <stdio.h>
#include "spectrum.h"
#include <errno.h>
double unified();
int Nmax();
double lcutoff();
double dmax();
double dmin();
int imin();
double lgffactor();

void lyman(hkap,model,wave)
atmosphere *model;
double hkap[];
double wave;
{
   static double w0[130] =
	  {1215.684,1025.734, 972.548, 949.753, 937.814, 930.758,
            926.236, 923.160, 920.973, 919.361, 918.139, 917.191,
            916.439, 915.834, 915.339, 914.929, 914.586, 914.296,
            914.048, 913.836, 913.651, 913.490, 913.349, 913.225,
            913.114, 913.016, 912.928, 912.849, 912.778, 912.713,
            912.655, 912.601, 912.553, 912.508, 912.467, 912.430,
            912.395, 912.363, 912.334, 912.306, 912.281, 912.257,
            912.235, 912.214, 912.194, 912.176, 912.159, 912.143,
            912.128, 912.114, 912.101, 912.088, 912.076, 912.065,
            912.054, 912.044, 912.034, 912.025, 912.017, 912.008,
            912.001, 911.993, 911.986, 911.979, 911.973, 911.966,
            911.961, 911.955, 911.949, 911.944, 911.939, 911.934,
            911.930, 911.925, 911.921, 911.917, 911.913, 911.909,
            911.906, 911.902, 911.899, 911.896, 911.893, 911.890,
            911.887, 911.884, 911.881, 911.878, 911.876, 911.873,
            911.871, 911.869, 911.867, 911.864, 911.862, 911.860,
            911.858, 911.856, 911.855, 911.853, 911.851, 911.849,
            911.848, 911.846, 911.844, 911.843, 911.842, 911.840,
            911.839, 911.837, 911.836, 911.835, 911.834, 911.832,
            911.831, 911.830, 911.829, 911.828, 911.827, 911.826,
            911.825, 911.824, 911.823, 911.822, 911.821, 911.820,
            911.819, 911.818, 911.817, 911.816};
   static double Ehigh[130] =
	  {10.198715,12.087366,12.748393,13.054355,13.220556,13.320770,
           13.385813,13.430406,13.462303,13.485904,13.503854,13.517823,
           13.528907,13.537850,13.545168,13.551233,13.556316,13.560618,
           13.564291,13.567451,13.570191,13.572581,13.574678,13.576529,
           13.578171,13.579633,13.580942,13.582117,13.583177,13.584136,
           13.585007,13.585799,13.586523,13.587186,13.587794,13.588353,
           13.588869,13.589346,13.589787,13.590197,13.590578,13.590932,
           13.591262,13.591571,13.591860,13.592131,13.592384,13.592623,
           13.592847,13.593058,13.593257,13.593445,13.593623,13.593791,
           13.593950,13.594101,13.594244,13.594380,13.594509,13.594632,
           13.594749,13.594860,13.594966,13.595068,13.595165,13.595257,
           13.595346,13.595430,13.595511,13.595589,13.595663,13.595735,
           13.595803,13.595869,13.595932,13.595993,13.596051,13.596107,
           13.596162,13.596214,13.596264,13.596312,13.596359,13.596404,
           13.596448,13.596490,13.596530,13.596570,13.596608,13.596644,
           13.596680,13.596714,13.596747,13.596780,13.596811,13.596841,
           13.596870,13.596899,13.596927,13.596953,13.596979,13.597005,
           13.597029,13.597053,13.597076,13.597099,13.597121,13.597142,
           13.597163,13.597183,13.597202,13.597221,13.597240,13.597258,
           13.597276,13.597293,13.597310,13.597326,13.597342,13.597358,
           13.597373,13.597388,13.597402,13.597416,13.597430,13.597443,
           13.597456,13.597469,13.597482,13.597494};


   extern float **bkap;
   int i,j,m,ifcore;
   double *wave0,*Eh,s,s0s;
   float rad = 500.0;
   float lya_rad = 2430.5;
   double minrad = 0.1;
   double fac = 4.0;
   static double lambda[4] = {1400.0,1400.0,1400.0,1400.0};
   static int nmax[NTAU];
   static double Cutoff[NTAU];
   static double gffac2[150];
   double *gffac;
   static int flag = 0;
   static int flag1 = 0;
   static double cut;
   double cutmax,cutmin;
   extern int Ntau;
   int NM = 129;
   extern memo reset;
   char tmp[10];
   FILE *tst;

   if(reset.lyman == 1) {
     flag = 0;
     reset.lyman = 0;
   }
   gffac = gffac2-2;
   ifcore = 0;
   wave0 = w0-2;
   Eh = Ehigh-2;

   if(wave < 911.6 || wave > 1215.7+lya_rad) return;
   if(flag1 == 0) {
     for(i=0;i<Ntau;i++) {
       Cutoff[i] = lcutoff(model->T[i],model->Ne[i],&NM);
       nmax[i] = NM;
     }
     for(m=2;m<130;m++) gffac[m] = lgffactor(m);
     flag1 = 1;
   }

   for(m=2;m<130;m++) {
     if(fabs(wave - wave0[m]) < 5.0) ifcore = 1;
   }
   if(ifcore == 0 && flag == 0) {
	lambda[1] = floor(wave);
	lambda[0] = lambda[1] - 1.0;
	lambda[2] = lambda[1] + 1.0;
	lambda[3] = lambda[1] + 2.0;
	for(i=0;i<Ntau;i++) {
	  bkap[0][i] = bkap[1][i] = bkap[2][i] = bkap[3][i] = 0.0;
	}
	for(m=2;m<130;m++) {
          rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(1000.0,rad);
	  if(m == 2) rad = lya_rad;
	  if(rad < minrad) rad = minrad; 
	  if(fabs(wave - wave0[m]) < rad) {
	    for(i=0;i<Ntau;i++) {
	      if(m > nmax[i] || wave < Cutoff[i]) continue;
	      for(j=0;j<4;j++) bkap[j][i] += gffac[m]*
			unified(1,m,i,0.000,Eh[m],model,lambda[j]);
	    }
	  }
	}
	flag = 1;
   }


   if(ifcore == 0 && wave >= lambda[2]) {
     for(j=0;j<4;j++) lambda[j] += 1.0;
       for(i=0;i<Ntau;i++) {
	 for(j=0;j<3;j++) bkap[j][i] = bkap[j+1][i];
	 bkap[3][i] = 0.0;
       }
       for(m=2;m<130;m++) {
	 rad = fac*(wave0[m] - wave0[m+1]);
	 rad = dmin(500.0,rad);
	 if(m == 2) rad = lya_rad; 
	 if(rad < minrad) rad = minrad; 
	 if(fabs(wave-wave0[m]) < rad) {
	   for(i=0;i<Ntau;i++) {
	     if(m > nmax[i] || wave < Cutoff[i]) continue;
	     bkap[3][i] += gffac[m]*
			unified(1,m,i,0.000,Eh[m],model,lambda[3]);
	   }
	 }
       }
   }

   if(ifcore == 0) {
     for(i=0;i<Ntau;i++) {
       hkap[i] = bkap[3][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[2]))/6.0 -
	     bkap[0][i]*((wave-lambda[1])*
	     (wave-lambda[2])*(wave-lambda[3]))/6.0 +
	     bkap[1][i]*((wave-lambda[0])*
	     (wave-lambda[2])*(wave-lambda[3]))/2.0 -
	     bkap[2][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[3]))/2.0;
       }
   }
   if(ifcore == 1) {
     for(i=0;i<Ntau;i++) hkap[i] = 0.0;
     for(m=2;m<130;m++) {
       rad = fac*(wave0[m] - wave0[m+1]);
       rad = dmin(500.0,rad);
       if(m == 2) rad = lya_rad;
       if(rad < minrad) rad = minrad; 
       if(fabs(wave - wave0[m]) < rad) {
	     for(i=0;i<Ntau;i++) {
	       if(m > nmax[i] || wave < Cutoff[i]) continue;
	       hkap[i] += gffac[m]*
			     unified(1,m,i,0.000,Eh[m],model,wave);
	     }
       }
     }
     flag = 0;
   }
}


double lgffactor(m)
int m;
{
  return(1.0);
}

double lcutoff(T,Ne,nm)
double Ne,T;
int *nm;
{
  double Eion = 109678.758;
  double E130 = 109672.167;
  double dE;
  double Cutoff;
  double fac = 1.0;

  /* Calculate lowering of ionization potential using Unsold Approx */
/* fac = sqrt(10000.0/T); */

  dE = 5.61e-03*fac*pow(Ne,1.0/3.0);
  *nm = Nmax(Eion - dE);
  dE = dmax(dE,Eion - E130);

  Cutoff = 1.0e+08/(109678.758 - dE);
  return(Cutoff);
}
