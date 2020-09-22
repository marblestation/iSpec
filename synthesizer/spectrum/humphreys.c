#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double unified();
double dmin();

void humphreys(hkap,model,wave)
atmosphere *model;
double hkap[];
double wave;
{
   double w0[] = {123685.04,75004.34,59065.91,51272.49,46712.26,43752.53,
                   41696.49,40197.63,39064.75,38184.02,37483.64,36916.19,
                   36449.22,36059.77,35731.23,35451.26,35210.58,35002.02,
                   34820.03,34660.20,34519.04,34393.70,34281.87,34181.66,
                   34091.49,34010.06,33936.25,33869.14};
   double Ehigh[] = {13.31834,13.38338,13.42796,13.45985,
                     13.48345,13.50139,13.51536,13.52644,13.53538,13.54270,
                     13.54877,13.55385,13.55815,13.56182,13.56498,13.56772,
                     13.57011,13.57221,13.57406,13.57570,13.57716,13.57847,
                     13.57964,13.58070,13.58166,13.58252,13.58333,13.58405};
   extern float **bkap4;
   int i,j,m,ifcore;
   double *wave0,*Eh;
   float rad = 400.0;
   double fac = 4.0;
   static double lambda[4];
   static int flag = 0;
   extern int Ntau;
   extern memo reset;

   if(reset.humphreys == 1) {
     flag = 0;
     reset.humphreys = 0;
   }

   if(wave < 33869 || wave > 125000.0) return;
   ifcore = 0;
   wave0 = w0-7;
   Eh = Ehigh-7;
   for(m=7;m<35;m++) {
     if(fabs(wave - wave0[m]) < 5.0) ifcore = 1;
   }
   if(ifcore == 0 && flag == 0) {
	lambda[1] = floor(wave);
	lambda[0] = lambda[1] - 1.0;
	lambda[2] = lambda[1] + 1.0;
	lambda[3] = lambda[1] + 2.0;
	for(i=0;i<Ntau;i++) {
	  bkap4[0][i] = bkap4[1][i] = bkap4[2][i] = bkap4[3][i] = 0.0;
	}
	for(m=7;m<35;m++) {
          rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(400.0,rad);
	  if(fabs(wave - wave0[m]) < rad) {
	    for(i=0;i<Ntau;i++) {
	      for(j=0;j<4;j++) bkap4[j][i] +=
			unified(6,m,i,13.21815,Eh[m],model,lambda[j]);
	    }
	  }
	}
	flag = 1;
   }


   if(ifcore == 0 && wave >= lambda[2]) {
     for(j=0;j<4;j++) lambda[j] += 1.0;
       for(i=0;i<Ntau;i++) {
	 for(j=0;j<3;j++) bkap4[j][i] = bkap4[j+1][i];
	 bkap4[3][i] = 0.0;
       }
       for(m=7;m<35;m++) {
         rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(400.0,rad);
	 if(fabs(wave-wave0[m]) < rad) {
	   for(i=0;i<Ntau;i++) bkap4[3][i] +=
			unified(6,m,i,13.21815,Eh[m],model,lambda[3]);
	 }
       }
   }

   if(ifcore == 0) {
     for(i=0;i<Ntau;i++) {
       hkap[i] = bkap4[3][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[2]))/6.0 -
	     bkap4[0][i]*((wave-lambda[1])*
	     (wave-lambda[2])*(wave-lambda[3]))/6.0 +
	     bkap4[1][i]*((wave-lambda[0])*
	     (wave-lambda[2])*(wave-lambda[3]))/2.0 -
	     bkap4[2][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[3]))/2.0;
       }
   }
   if(ifcore == 1) {
     for(i=0;i<Ntau;i++) hkap[i] = 0.0;
     for(m=7;m<35;m++) {
       rad = fac*(wave0[m] - wave0[m+1]);
       rad = dmin(400.0,rad);
       if(fabs(wave - wave0[m]) < rad) {
	     for(i=0;i<Ntau;i++) hkap[i] +=
			     unified(6,m,i,13.21815,Eh[m],model,wave);
       }
     }
     flag = 0;
   }
}


