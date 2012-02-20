#include <math.h>
#include <stdio.h>
#include "spectrum.h"
#include <errno.h>
double unified();
double dmin();

void brackett(hkap,model,wave)
atmosphere *model;
double hkap[];
double wave;
{
   double w0[] = {40511.6,26251.5,21655.3,19445.6,18174.1,17362.1,
                  16806.5,16407.2,16109.3,15880.5,15700.7,15556.4,
                  15438.9,15341.8,15260.5,15191.8,15133.2,15082.8,
                  15039.0,15000.9,14967.3,14937.7,14911.4,14888.0,
                  14867.0,14848.1,14831.1,14815.7,14801.6,14788.9};
   double Ehigh[] = {13.05198,13.21815,13.31834,13.38338,13.42796,13.45985,
                     13.48345,13.50139,13.51536,13.52644,13.53538,13.54270,
                     13.54877,13.55385,13.55815,13.56182,13.56498,13.56772,
                     13.57011,13.57221,13.57406,13.57570,13.57716,13.57847,
                     13.57964,13.58070,13.58166,13.58252,13.58333,13.58405};
   extern float **bkap2;
   int i,j,m,ifcore;
   double *wave0,*Eh;
   float rad = 400.0;
   double fac = 4.0;
   static double lambda[4];
   static int flag = 0;
   extern int Ntau;
   extern memo reset;

   if(reset.brackett == 1) {
     flag = 0;
     reset.brackett = 0;
   }

   if(wave < 14779.0 || wave > 42000.0) return;
   ifcore = 0;
   wave0 = w0-5;
   Eh = Ehigh-5;
   for(m=5;m<35;m++) {
     if(fabs(wave - wave0[m]) < 5.0) ifcore = 1;
   }
   if(ifcore == 0 && flag == 0) {
	lambda[1] = floor(wave);
	lambda[0] = lambda[1] - 1.0;
	lambda[2] = lambda[1] + 1.0;
	lambda[3] = lambda[1] + 2.0;
	for(i=0;i<Ntau;i++) {
	  bkap2[0][i] = bkap2[1][i] = bkap2[2][i] = bkap2[3][i] = 0.0;
	}
	for(m=5;m<35;m++) {
          rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(400.0,rad);
	  if(fabs(wave - wave0[m]) < rad) {
	    for(i=0;i<Ntau;i++) {
	      for(j=0;j<4;j++) bkap2[j][i] +=
			unified(4,m,i,12.74607,Eh[m],model,lambda[j]);
	    }
	  }
	}
	flag = 1;
   }


   if(ifcore == 0 && wave >= lambda[2]) {
     for(j=0;j<4;j++) lambda[j] += 1.0;
       for(i=0;i<Ntau;i++) {
	 for(j=0;j<3;j++) bkap2[j][i] = bkap2[j+1][i];
	 bkap2[3][i] = 0.0;
       }
       for(m=5;m<35;m++) {
         rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(400.0,rad);
	 if(fabs(wave-wave0[m]) < rad) {
	   for(i=0;i<Ntau;i++) bkap2[3][i] +=
			unified(4,m,i,12.74607,Eh[m],model,lambda[3]);
	 }
       }
   }

   if(ifcore == 0) {
     for(i=0;i<Ntau;i++) {
       hkap[i] = bkap2[3][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[2]))/6.0 -
	     bkap2[0][i]*((wave-lambda[1])*
	     (wave-lambda[2])*(wave-lambda[3]))/6.0 +
	     bkap2[1][i]*((wave-lambda[0])*
	     (wave-lambda[2])*(wave-lambda[3]))/2.0 -
	     bkap2[2][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[3]))/2.0;
       }
   }
   if(ifcore == 1) {
     for(i=0;i<Ntau;i++) hkap[i] = 0.0;
     for(m=5;m<35;m++) {
       rad = fac*(wave0[m] - wave0[m+1]);
       rad = dmin(400.0,rad);
       if(fabs(wave - wave0[m]) < rad) {
	     for(i=0;i<Ntau;i++) hkap[i] +=
			     unified(4,m,i,12.74607,Eh[m],model,wave);
       }
     }
     flag = 0;
   }
}


