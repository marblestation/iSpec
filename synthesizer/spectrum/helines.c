/* The neutral helium line profiles contained in modules he14a.c
to he617a.c are from calculations by Beauchamp, Wesemael and Bergeron
(see ApJS 108, 559). */ 

#include <stdio.h>
float he3705_14(),he3705_314(),he3705_15(),he3705_315(),he3705_16();
float he3705_316(),he3705_17(),he3705_317(),he3705_617();
float he3819_14(),he3819_314(),he3819_15(),he3819_315(),he3819_16();
float he3819_316(),he3819_17(),he3819_317(),he3819_617();
float he3867_14(),he3867_314(),he3867_15(),he3867_315(),he3867_16();
float he3867_316(),he3867_17(),he3867_317(),he3867_617();
float he3888_14(),he3888_314(),he3888_15(),he3888_315(),he3888_16();
float he3888_316(),he3888_17(),he3888_317(),he3888_617();
float he3964_14(),he3964_314(),he3964_15(),he3964_315(),he3964_16();
float he3964_316(),he3964_17(),he3964_317(),he3964_617();
float he4009_14(),he4009_314(),he4009_15(),he4009_315(),he4009_16();
float he4009_316(),he4009_17(),he4009_317(),he4009_617();
float he4023_14(),he4023_314(),he4023_15(),he4023_315(),he4023_16();
float he4023_316(),he4023_17(),he4023_317(),he4023_617();
float he4026_14(),he4026_314(),he4026_15(),he4026_315(),he4026_16();
float he4026_316(),he4026_17(),he4026_317(),he4026_617();
float he4120_14(),he4120_314(),he4120_15(),he4120_315(),he4120_16();
float he4120_316(),he4120_17(),he4120_317(),he4120_617();
float he4143_14(),he4143_314(),he4143_15(),he4143_315(),he4143_16();
float he4143_316(),he4143_17(),he4143_317(),he4143_617();
float he4168_14(),he4168_314(),he4168_15(),he4168_315(),he4168_16();
float he4168_316(),he4168_17(),he4168_317(),he4168_617();
float he4387_14(),he4387_314(),he4387_15(),he4387_315(),he4387_16();
float he4387_316(),he4387_17(),he4387_317(),he4387_617();
float he4437_14(),he4437_314(),he4437_15(),he4437_315(),he4437_16();
float he4437_316(),he4437_17(),he4437_317(),he4437_617();
float he4471_12(),he4471_13(),he4471_313(),he4471_14(),he4471_314();
float he4471_15(),he4471_315(),he4471_16(),he4471_17();
float he4713_14(),he4713_314(),he4713_15(),he4713_315(),he4713_16();
float he4713_316(),he4713_17(),he4713_317(),he4713_617();
float he4922_13(),he4922_313(),he4922_14(),he4922_314();
float he4922_15(),he4922_315(),he4922_16(),he4922_17();
float he5016_316(),he5016_17(),he5016_317(),he5016_617();
float he5047_14(),he5047_314(),he5047_15(),he5047_315(),he5047_16();
float he5047_316(),he5047_17(),he5047_317(),he5047_617();
float he5875_14(),he5875_314(),he5875_15(),he5875_315(),he5875_16();
float he5875_316(),he5875_17(),he5875_317(),he5875_617();
float he6678_14(),he6678_314(),he6678_15(),he6678_315(),he6678_16();
float he6678_316(),he6678_17(),he6678_317(),he6678_617();

void broadparam();
double jx();

double he3705(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he3705_14,he3705_314,he3705_15,he3705_315,he3705_16,
                      he3705_316,he3705_17,he3705_317,he3705_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;

     
 if(Ne < 1.0e+14) {
   broadparam(3705.00,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 3705.00)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he3705_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;
 
 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }
                 
 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he3819(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he3819_14,he3819_314,he3819_15,he3819_315,he3819_16,
                      he3819_316,he3819_17,he3819_317,he3819_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;
 
     
 if(Ne < 1.0e+14) {
   broadparam(3819.60,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 3819.60)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he3819_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;
 
 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }
                 
 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he3867(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he3867_14,he3867_314,he3867_15,he3867_315,he3867_16,
                      he3867_316,he3867_17,he3867_317,he3867_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;
 
     
 if(Ne < 1.0e+14) {
   broadparam(3867.50,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 3867.50)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he3867_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;
 
 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }
                 
 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he3888(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he3888_14,he3888_314,he3888_15,he3888_315,he3888_16,
                      he3888_316,he3888_17,he3888_317,he3888_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;
 
     
 if(Ne < 1.0e+14) {
   broadparam(3888.65,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 3888.65)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he3888_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;
 
 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }
                 
 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he3964(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he3964_14,he3964_314,he3964_15,he3964_315,he3964_16,
                      he3964_316,he3964_17,he3964_317,he3964_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(3964.73,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 3964.73)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he3964_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4009(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4009_14,he4009_314,he4009_15,he4009_315,he4009_16,
                      he4009_316,he4009_17,he4009_317,he4009_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


if(Ne < 1.0e+14) {
   broadparam(4009.27,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4009.27)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 } 
 if(Ne >= 6.0e+17) return(he4009_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);
 /* printf("4009  k = %d, p1 = %e\n",k,p1);
 printf("4009  k+1 = %d, p2 = %e\n",k+1,p2); */

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4023(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4023_14,he4023_314,he4023_15,he4023_315,he4023_16,
                      he4023_316,he4023_17,he4023_317,he4023_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4023.95,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4023.95)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4023_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4026(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4026_14,he4026_314,he4026_15,he4026_315,he4026_16,
                      he4026_316,he4026_17,he4026_317,he4026_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4026.20,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4026.20)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4026_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4120(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4120_14,he4120_314,he4120_15,he4120_315,he4120_16,
                      he4120_316,he4120_17,he4120_317,he4120_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4120.80,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4120.80)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4120_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4143(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4143_14,he4143_314,he4143_15,he4143_315,he4143_16,
                      he4143_316,he4143_17,he4143_317,he4143_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4143.76,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4143.76)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4143_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4168(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4168_14,he4168_314,he4168_15,he4168_315,he4168_16,
                      he4168_316,he4168_17,he4168_317,he4168_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4168.97,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4168.97)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4168_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4387(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4387_14,he4387_314,he4387_15,he4387_315,he4387_16,
                      he4387_316,he4387_17,he4387_317,he4387_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4387.93,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4387.93)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4387_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}


double he4437(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4437_14,he4437_314,he4437_15,he4437_315,he4437_16,
                      he4437_316,he4437_17,he4437_317,he4437_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4437.55,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4437.55)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4437_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}


double he4471(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4471_12,he4471_13,he4471_313,he4471_14,he4471_314,
		      he4471_15,he4471_315,he4471_16,he4471_17};
 double NE[9] = {1.0e+12,1.0e+13,3.0e+13,1.0e+14,3.0e+14,1.0e+15,3.0e+15,
		 1.0e+16,1.0e+17};
 double p,p1,p2;
 int i,k;

 if(Ne <= 1.0e+12) return(he4471_12(T,wave));
 if(Ne > 1.0e+17) return(he4471_17(T,wave));
 if(T < 5000.0 || T > 40000.0) return(0.0);

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}


double he4713(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he4713_14,he4713_314,he4713_15,he4713_315,he4713_16,
                      he4713_316,he4713_17,he4713_317,he4713_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(4713.20,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4713.20)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he4713_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he4922(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[8]) () = {he4922_13,he4922_313,he4922_14,he4922_314,
		      he4922_15,he4922_315,he4922_16,he4922_17};
 double NE[8] = {1.0e+13,3.0e+13,1.0e+14,3.0e+14,1.0e+15,3.0e+15,
		 1.0e+16,1.0e+17};
 double p,p1,p2,omega,w,dw,a,sig;
 int i,k;

 if(Ne <= 1.0e+13) {
   broadparam(4921.93,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 4921.93)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }

 if(Ne > 1.0e+17) return(he4922_17(T,wave));
 if(T < 5000.0 || T > 40000.0) return(0.0);

 for(i=0;i<7;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}

double he5016(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[4]) () = {he5016_316,he5016_17,he5016_317,he5016_617};
 double NE[8] = {3.0e+16,1.0e+17,3.0e+17,6.0e+17};
 double p,p1,p2,omega,w,dw,a,sig;
 int i,k;

 if(Ne <= 3.0e+16) {
   broadparam(5015.68,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 5015.68)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }

 if(Ne > 6.0e+17) return(he5016_617(T,wave));
 if(T < 5000.0 || T > 40000.0) return(0.0);

 for(i=0;i<3;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}


double he5047(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he5047_14,he5047_314,he5047_15,he5047_315,he5047_16,
                      he5047_316,he5047_17,he5047_317,he5047_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(5047.74,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 5047.74)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he5047_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}


double he5875(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he5875_14,he5875_314,he5875_15,he5875_315,he5875_16,
                      he5875_316,he5875_17,he5875_317,he5875_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(5875.70,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 5875.70)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he5875_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}


double he6678(wave,T,Ne)
double wave,T,Ne;
{
 float (*pf[9]) () = {he6678_14,he6678_314,he6678_15,he6678_315,he6678_16,
                      he6678_316,he6678_17,he6678_317,he6678_617};
 double NE[9] = {1.0e+14,3.0e+14,1.0e+15,3.0e+15,1.0e+16,3.0e+16,1.0e+17,
                 3.0e+17,6.0e+17};
 int i,k;
 double p,p1,p2,w,dw,a,sig,omega;


 if(Ne < 1.0e+14) {
   broadparam(6678.15,Ne,T,&w,&dw,&a,&sig);
   omega = (wave - 6678.15)/w + dw;
   p = jx(omega,a,sig)/w;
   return(p);
 }
 if(Ne >= 6.0e+17) return(he6678_617(T,wave));
 if(T < 10000.0) T = 10000.0;
 if(T > 40000.0) T = 39999.9;

 for(i=0;i<8;i++) {
   if(Ne <= NE[i+1]) {
     k = i;
     break;
   }
 }

 p1 = pf[k] (T,wave);
 p2 = pf[k+1] (T,wave);

 p = p1 + (p2 - p1)*(Ne - NE[k])/(NE[k+1] - NE[k]);
 return(p);
}
