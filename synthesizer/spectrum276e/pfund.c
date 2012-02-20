#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double unified();
double dmin();

void pfund(hkap,model,wave)
atmosphere *model;
double hkap[];
double wave;
{
   double w0[] = {74578.17,46525.04,37395.32,32960.89,30383.70,28722.09,
                  27575.12,26743.98,26119.32,25636.24,25253.98,24945.70,
                  24693.10,24483.29,24306.95,24157.23,24028.93,23918.12,
                  23821.70,23737.26,23662.88,23596.99,23538.36,23485.92,
                  23438.85,23396.42,23358.03,23323.20,23291.48};
   double Ehigh[] = {13.21815,13.31834,13.38338,13.42796,13.45985,
                     13.48345,13.50139,13.51536,13.52644,13.53538,13.54270,
                     13.54877,13.55385,13.55815,13.56182,13.56498,13.56772,
                     13.57011,13.57221,13.57406,13.57570,13.57716,13.57847,
                     13.57964,13.58070,13.58166,13.58252,13.58333,13.58405};
   extern float **bkap3;
   int i,j,m,ifcore;
   double *wave0,*Eh;
   float rad = 250.0;
   double fac = 4.0;
   static double lambda[4];
   static int flag = 0;
   extern int Ntau;
   extern memo reset;

   if(reset.pfund == 1) {
     flag = 0;
     reset.pfund = 0;
   }

   if(wave < 23291.0 || wave > 75500.0) return;
   ifcore = 0;
   wave0 = w0-6;
   Eh = Ehigh-6;
   for(m=6;m<35;m++) {
     if(fabs(wave - wave0[m]) < 5.0) ifcore = 1;
   }
   if(ifcore == 0 && flag == 0) {
	lambda[1] = floor(wave);
	lambda[0] = lambda[1] - 1.0;
	lambda[2] = lambda[1] + 1.0;
	lambda[3] = lambda[1] + 2.0;
	for(i=0;i<Ntau;i++) {
	  bkap3[0][i] = bkap3[1][i] = bkap3[2][i] = bkap3[3][i] = 0.0;
	}
	for(m=6;m<35;m++) {
          rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(250.0,rad);
	  if(fabs(wave - wave0[m]) < rad) {
	    for(i=0;i<Ntau;i++) {
	      for(j=0;j<4;j++) bkap3[j][i] +=
			unified(5,m,i,13.05198,Eh[m],model,lambda[j]);
	    }
	  }
	}
	flag = 1;
   }


   if(ifcore == 0 && wave >= lambda[2]) {
     for(j=0;j<4;j++) lambda[j] += 1.0;
       for(i=0;i<Ntau;i++) {
	 for(j=0;j<3;j++) bkap3[j][i] = bkap3[j+1][i];
	 bkap3[3][i] = 0.0;
       }
       for(m=6;m<35;m++) {
         rad = fac*(wave0[m] - wave0[m+1]);
	  rad = dmin(250.0,rad);
	 if(fabs(wave-wave0[m]) < rad) {
	   for(i=0;i<Ntau;i++) bkap3[3][i] +=
			unified(5,m,i,13.05198,Eh[m],model,lambda[3]);
	 }
       }
   }

   if(ifcore == 0) {
     for(i=0;i<Ntau;i++) {
       hkap[i] = bkap3[3][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[2]))/6.0 -
	     bkap3[0][i]*((wave-lambda[1])*
	     (wave-lambda[2])*(wave-lambda[3]))/6.0 +
	     bkap3[1][i]*((wave-lambda[0])*
	     (wave-lambda[2])*(wave-lambda[3]))/2.0 -
	     bkap3[2][i]*((wave-lambda[0])*
	     (wave-lambda[1])*(wave-lambda[3]))/2.0;
       }
   }
   if(ifcore == 1) {
     for(i=0;i<Ntau;i++) hkap[i] = 0.0;
     for(m=6;m<35;m++) {
       rad = fac*(wave0[m] - wave0[m+1]);
       rad = dmin(250.0,rad);
       if(fabs(wave - wave0[m]) < rad) {
	     for(i=0;i<Ntau;i++) hkap[i] +=
			     unified(5,m,i,13.05198,Eh[m],model,wave);
       }
     }
     flag = 0;
   }
}


