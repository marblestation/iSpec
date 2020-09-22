#include "spectrum.h"
#include <math.h>
#include <errno.h>
#include <string.h>
struct strong {
  double wave;
  double code;
  int iso;
  double El;
  double Eu;
  double loggf;
  double radius;
  double dampfac;
  char T[5];
  double SA;
  int flag;
};

void pop();
void broad(atmosphere *model,linedata *line,int N,double sig,
           double alp,double fac);
void capnu();
double voigt();
double abund();
int approx();

double stronglines(model,atom,n,strgln,wave,V,POP)
atmosphere *model;
atominfo *atom;
pfunc *V;
population *POP;
int n;
linedata *strgln;
double wave;
{
 static struct strong slines[] = {
   { 2585.876,26.1,  0,       0, 38660,-0.190,150.0,1.00,"AO",167.225,0},
   { 2598.370,26.1,  0,     385, 38859,-0.100,150.0,1.00,"AO",167.222,0},
   { 2599.396,26.1,  0,       0, 38459, 0.350,250.0,1.00,"AO",166.228,0},
   { 2607.088,26.1,  0,     668, 39013,-0.170,150.0,1.00,"AO",167.224,0},
   { 2611.874,26.1,  0,     385, 38660,-0.050,150.0,1.00,"AO",167.226,0},
   { 2613.825,26.1,  0,     863, 39109,-0.390,150.0,1.00,"AO",167.223,0},
   { 2617.618,26.1,  0,     668, 38859,-0.570,150.0,1.00,"AO",167.222,0},
   { 2625.668,26.1,  0,     385, 38459,-0.460,150.0,1.00,"AO",166.230,0},
   { 2628.294,26.1,  0,     977, 39013,-0.450,150.0,1.00,"AO",167.223,0},
   { 2631.048,26.1,  0,     863, 38859,-0.320,150.0,1.00,"AO",167.222,0},
   { 2631.324,26.1,  0,     668, 38660,-0.300,150.0,1.00,"AO",167.227,0},
   { 2795.528,12.1,  0,     0.0, 35761, 0.085,500.0,0.75,"AO",246.357,0},
   { 2802.705,12.1,  0,     0.0, 35669,-0.218,500.0,0.75,"AO",246.357,0},
   { 2852.126,12.0,  0,     0.0, 35051, 0.255,3000.0,0.75,"01",000.000,0},
   { 3933.66, 20.1,  0,     0.0, 25414, 0.140,400.0,0.75,"AO",234.223,0},
   { 3968.47, 20.1,  0,     0.0, 25192,-0.162,400.0,0.75,"AO",234.223,0},
   { 4077.715,38.1,  0,     0.0, 24517, 0.158, 15.0,0.80,"AO",262.208,0},
   { 4226.730,20.0,  0,     0.0, 23652, 0.243,3000.0,1.00,"01",000.000,0},
   { 4246.830,21.1,  0,  2541.0, 26081, 0.322, 15.0,1.00,"99",000.000,0},
   { 5889.951,11.0,  0,     0.0, 16973, 0.108,1000.0,0.75,"01",000.000,0},
   { 5895.924,11.0,  0,     0.0, 16956,-0.194,1000.0,0.75,"01",000.000,0},
   { 7664.899,19.0,  0,     0.0, 13043, 0.127,1000.0,1.00,"01",000.000,0},
   { 7698.965,19.0,  0,     0.0, 12985,-0.176,1000.0,1.00,"01",000.000,0}};
 double v,kappa,dl,code;
 int i,j,l;
 int icode = 1;
 extern int Ntau;
 extern memo reset;

 kappa = 0.0;

 if(reset.strong == 1) {
   for(i=0;i<NSTRG;i++) slines[i].flag = 0;
   reset.strong = 0;
 }

 for(i=0;i<NSTRG;i++) {
   if(fabs(wave - slines[i].wave) <= slines[i].radius) {
     if(slines[i].flag == 0) {
       slines[i].flag = 1;
       strgln[i].wave = slines[i].wave;
       code = strgln[i].code = slines[i].code;
       if(approx(code,floor(code),0.001) == 1) icode = 0;
       else if(approx(code,floor(code)+0.1,0.001) == 1) icode = 1;
       else if(approx(code,floor(code)+0.2,0.001) == 1) icode = 2;
       else if(approx(code,floor(code)+0.3,0.001) == 1) icode = 3;
       else icode = 0;
       for(j=0;j<NATOM;j++) {
	 if((int)floor(code) == atom[j].code) {
	   l = j;
	   break;
	 }
       }
       strgln[i].atomass = atom[l].amass;
       strgln[i].abund = abund(atom,(int)floor(code));
       strgln[i].chi1 = atom[l].I1;
       strgln[i].chi2 = atom[l].I2;
       strgln[i].chi3 = atom[l].I3;
       strgln[i].chi4 = atom[l].I4;
       if(icode == 0) strgln[i].chi = atom[l].I1;
       else if(icode == 1) strgln[i].chi = atom[l].I2;
       else if(icode == 2) strgln[i].chi = atom[l].I3;
       else if(icode == 3) strgln[i].chi = atom[l].I4;
       else strgln[i].chi = atom[l].I2;
       strgln[i].El = slines[i].El*1.23981e-04;
       strgln[i].Eu = slines[i].Eu*1.23981e-04;
       strgln[i].gf = pow(10.0,slines[i].loggf);
       strgln[i].sig = floor(slines[i].SA);
       strgln[i].alp = slines[i].SA - floor(slines[i].SA);
       strgln[i].fac = slines[i].dampfac;
       strcpy(strgln[i].T,slines[i].T);
       for(j=0;j<Ntau;j++) strgln[i].xnum[j] = strgln[i].a[j] =
			   strgln[i].dopp[j] = 0.0;
       pop(strgln,i,model,V,POP);
       broad(model,strgln,i,strgln[i].sig,strgln[i].alp,strgln[i].fac);
       capnu(strgln,i,model);
     }
     v = 2.997929e+10*(strgln[i].wave - wave)/
	 (strgln[i].wave*strgln[i].dopp[n]);
     /* The following factor involving dl takes into account the quasi- */
     /* static contribution to van der Waal's broadening.  I have applied */
     /* it only to the red-wing, as this gives the best fit to the solar */
     /* Ca II K and H profiles, for reasons which I do not yet understand */
     
     dl = wave - strgln[i].wave;
     if(dl < 0.0) dl = 0.0;
     /* Note: factor set to zero, March 13, 2001 */
     dl = 0.0; 
     kappa += strgln[i].capnu[n]*voigt((double)strgln[i].a[n],v)*
     (1.0+dl/strgln[i].dlg[n]);  
   }
 }

 return(kappa);
}

