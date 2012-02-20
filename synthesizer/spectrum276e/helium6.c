#include <math.h>
#include "spectrum.h"
double heprof();
double he3705();
double he3819();
double he3867();
double he3888();
double he3964();
double he4009();
double he4023();
double he4026();
double he4120();
double he4143();
double he4168();
double he4387();
double he4437();
double he4471();
double he4713();
double he4922();
double he5016();
double he5047();
double he5875();
double he6678();
double pfunction();
void getparam();
void broadparam();

double helium(model,n,wave,V,He)
atmosphere *model;
int n;
double wave;
pfunc *V;
Helium *He;
{
  double kappa = 0.0;
  double c = 2.997924562e+18;
  double kap,kT,U;
  int i;
  extern memo reset;
  static int flag[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			 -1,-1,-1,-1,-1,-1,-1};

  if(reset.helium == 1) {
    for(i=0;i<27;i++) flag[i] = -1;
    for(i=0;i<NHE;i++) He[i].flag = 0;
    reset.helium = 0;
  }

  if(model->teff < 7501.0) return(0.0);
  kT = model->kT[n];
  U =  pfunction(V,2.0,n);
  
  for(i=0;i<NHE;i++) {
    if(He[i].flag == 1 && He[i].end < wave) He[i].flag = 0;
  }

  if(fabs(wave - 2829.07) <= 20.0) {
     if(flag[0] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[0] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(2829.07,flag[0],20.0,model,He,V,19.8146,0.0182);
     }
     kappa += He[flag[0]].kap[n]*wave*wave*heprof(2829.07,wave,n,He,flag[0]);
  }

  if(fabs(wave - 2945.10) <= 20.0) {
     if(flag[1] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[1] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(2945.10,flag[1],20.0,model,He,V,19.8146,0.0343);
     }
     kappa += He[flag[1]].kap[n]*wave*wave*heprof(2945.10,wave,n,He,flag[1]);
  }


  if(fabs(wave - 3187.74) <= 20.0) {
     if(flag[2] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[2] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3187.74,flag[2],20.0,model,He,V,19.8146,0.0693);
     }
     kappa += He[flag[2]].kap[n]*wave*wave*heprof(3187.74,wave,n,He,flag[2]);
  }

  if(fabs(wave - 3354.55) <= 20.0) {
     if(flag[3] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[3] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3354.55,flag[3],20.0,model,He,V,20.6106,0.0066);
     }
     kappa += He[flag[3]].kap[n]*wave*wave*heprof(3354.55,wave,n,He,flag[3]);
  }

  if(fabs(wave - 3447.59) <= 20.0) {
     if(flag[4] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[4] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3447.59,flag[4],20.0,model,He,V,20.6106,0.0128);
     }
     kappa += He[flag[4]].kap[n]*wave*wave*heprof(3447.59,wave,n,He,flag[4]);
  }

  if(fabs(wave - 3613.64) <= 20.0) {
     if(flag[5] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[5] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3613.64,flag[5],20.0,model,He,V,20.6106,0.0221);
     }
     kappa += He[flag[5]].kap[n]*wave*wave*heprof(3613.64,wave,n,He,flag[5]);
  }

  if(fabs(wave - 3652.00) <= 20.0) {
     if(flag[6] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[6] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3652.00,flag[6],20.0,model,He,V,20.9588,0.0065);
     }
     kappa += He[flag[6]].kap[n]*wave*wave*heprof(3652.00,wave,n,He,flag[6]);
  }

  if(fabs(wave - 3705.00) <= 80.0) {
	kap = model->NHeI[n]*9.0*exp(-20.9558/kT)*model->stim[n]/U;
	kap *= 2.65386e-2*0.0152*wave*wave*he3705(wave,model->T[n],model->Ne[n])/c;
	kappa += kap;
  }

  if(fabs(wave - 3819.60) <= 200.0) {
	kap = model->NHeI[n]*9.0*exp(-20.9558/kT)*model->stim[n]/U;
	kap *= 2.65386e-2*0.0215*wave*wave*he3819(wave,model->T[n],model->Ne[n])/c;
	kappa += kap;
  }

  if(fabs(wave - 3867.50) <= 200.0) {
	kap = model->NHeI[n]*9.0*exp(-20.9558/kT)*model->stim[n]/U;
	kap *= 2.65386e-2*0.00176*wave*wave*he3867(wave,model->T[n],model->Ne[n])/c;
	kappa += kap;
  }

  if(wave >= 3688.65 && wave <= 4388.65) {
        kap = model->NHeI[n]*3.0*exp(-19.8158/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.06446*wave*wave*he3888(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

if(fabs(wave - 3926.53) <= 20.0) {
     if(flag[7] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[7] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3926.53,flag[7],20.0,model,He,V,21.2127,0.0225);
     }
     kappa += He[flag[7]].kap[n]*wave*wave*heprof(3926.53,wave,n,He,flag[7]);
  }  

  if(fabs(wave - 3935.91) <= 20.0) {
     if(flag[8] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[8] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(3935.91,flag[8],20.0,model,He,V,21.2127,0.001668);
     }
     kappa += He[flag[8]].kap[n]*wave*wave*heprof(3935.91,wave,n,He,flag[8]);
  }

  if(fabs(wave - 3964.73) <= 200.0) {
        kap = model->NHeI[n]*1.0*exp(-20.6118/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.0507*wave*wave*he3964(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(wave >= 3929.27 && wave <= 4069.27) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.0112*wave*wave*he4009(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }  

  if(fabs(wave - 4023.95) <= 200.0) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*8.81e-04*wave*wave*he4023(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4026.20) <= 200.0) {
        kap = model->NHeI[n]*9.0*exp(-20.9601/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.0474*wave*wave*he4026(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4120.80) <= 200.0) {
        kap = model->NHeI[n]*9.0*exp(-20.9601/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.00365*wave*wave*he4120(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4143.76) <= 200.0) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.0213*wave*wave*he4143(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4168.97) <= 200.0) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.00153*wave*wave*he4168(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4387.93) <= 200.0) {
	kap = model->NHeI[n]*3.0*exp(-21.2127/kT)*model->stim[n]/U;
	kap *= 2.65386e-2*0.0436*wave*wave*he4387(wave,model->T[n],model->Ne[n])/c;
	kappa += kap;
  }

  if(fabs(wave - 4437.55) <= 200.0) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.00308*wave*wave*he4437(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4471.48) <= 50.0) {
     kap = model->NHeI[n]*9.0*exp(-20.9588/kT)*model->stim[n]/U;
     kap *= 2.65386e-2*0.125*wave*wave*he4471(wave,model->T[n],model->Ne[n])/c;
     kappa += kap;
  }

  if(wave >= 4433.2 && wave <= 4913.2) {
        kap = model->NHeI[n]*9.0*exp(-20.9601/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.0118*wave*wave*he4713(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 4921.93) <= 50.0) {
     kap = model->NHeI[n]*3.0*exp(-21.2127/kT)*model->stim[n]/U;
     kap *= 2.65386e-2*0.122*wave*wave*he4922(wave,model->T[n],model->Ne[n])/c;
     kappa += kap;
  }

  if(fabs(wave - 5015.68) <= 50.0) {
     kap = model->NHeI[n]*1.0*exp(-20.6106/kT)*model->stim[n]/U;
     kap *= 2.65386e-2*0.1514*wave*wave*he5016(wave,model->T[n],model->Ne[n])/c;
     kappa += kap;
  }

  if(fabs(wave - 5047.74) <= 200.0) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.00834*wave*wave*he5047(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }  

  if(fabs(wave - 5875.7) <= 300.0) {
	kap = model->NHeI[n]*9.0*exp(-20.9588/kT)*model->stim[n]/U;
	kap *= 2.65386e-2*0.609*wave*wave*he5875(wave,model->T[n],model->Ne[n])/c;
	kappa += kap;
  }

  if(fabs(wave - 6678.15) <= 200.0) {
        kap = model->NHeI[n]*3.0*exp(-21.2139/kT)*model->stim[n]/U;
        kap *= 2.65386e-2*0.711*wave*wave*he6678(wave,model->T[n],model->Ne[n])/c;
        kappa += kap;
  }

  if(fabs(wave - 7065.19) <= 20.0) {
     if(flag[9] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[9] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(7065.19,flag[9],20.0,model,He,V,20.9588,1.1072);
     }
     kappa += He[flag[9]].kap[n]*wave*wave*heprof(7065.19,wave,n,He,flag[9]);
  }

  if(fabs(wave - 7281.35) <= 20.0) {
     if(flag[10] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[10] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
       getparam(7281.35,flag[10],20.0,model,He,V,21.2127,0.144);
     }
     kappa += He[flag[10]].kap[n]*wave*wave*heprof(7281.35,wave,n,He,flag[10]);
  }

  if(fabs(wave - 8361.77) <= 20.0) {
     if(flag[11] == -1) {
       for(i=0;i<NHE;i++) {
	 if(He[i].flag == 0) {
	   flag[11] = i;
	   He[i].flag = 1;
	   break;
	 }
       }
     getparam(8361.77,flag[11],20.0,model,He,V,22.7128,0.0068);
     }
     kappa += He[flag[11]].kap[n]*wave*wave*heprof(8361.77,wave,n,He,flag[11]);
  }


  return(kappa);
}

void getparam(lambda,flag,radius,model,He,V,El,gf)
double lambda,radius,El,gf;
int flag;
atmosphere *model;
Helium *He;
pfunc *V;
{
  int i;
  double w,dw,a,sig,U,kT;
  double c = 2.997924562e+18;
  extern int Ntau;

  He[flag].linecenter = lambda;
  He[flag].end = lambda+radius;

  for(i=0;i<Ntau;i++) {
    kT = model->kT[i];
    U =  pfunction(V,2.0,i);
    broadparam(lambda,model->Ne[i],model->T[i],&w,&dw,&a,&sig);
    He[flag].w[i] = w;
    He[flag].dw[i] = dw;
    He[flag].a[i] = a;
    He[flag].sig[i] = sig;
    He[flag].kap[i] = 2.65386e-2*model->NHeI[i]*gf*
			exp(-El/kT)*model->stim[i]/(U*c);
  }
}
