#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "spectrum.h"
double partfn(double code,double T,double Ne);
double *dvector();
double **dmatrix();
double abund();
double Xi();
void gaussj();
int imin();
void free_dmatrix();
void free_dvector();
void printhotDensity();

void hotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw)
{
   double T,Ne,NA,U,P,ntot,T32,kT,ne2,ne3,ne4,ne5,rD,DE;
   double phe2,phe3,ph,ph2p,ph2,phm,ps2,ps3,ps4,ps5;
   double pc2,pc3,pc4,pc5,pna2,pna3,pna4,pna5,pmg2,pmg3,pmg4,pmg5;
   double pal2,pal3,pal4,pal5,psi2,psi3,psi4,psi5,pk2,pk3,pk4,pk5;
   double pca2,pca3,pca4,pca5,pfe2,pfe3,pfe4,pfe5,pn2,pn3,pn4,pn5;
   double po2,po3,po4,po5,pti2,pti3;
   double ncI,nnI,noI,nnaI,nmgI,nalI,nsiI,nsI,nkI,ncaI,nfeI,ntiI;
   double ac,an,ao,ana,amg,aal,asi,as,ak,aca,afe,ati;
   double k = 8.617084e-05;
   double kerg = 1.38054e-16;
   double delta = 1.0;
   double **df,**f,*n;
   double neorig;
   double rho;
   int i,j;
   extern int Ntau;
   extern int flagmgh;
   extern int flagP;

/* n[1] = nA, n[2] = nHI, n[3] = nHeI, n[4] = nCI, n[5] = nNI */
/* n[6] = nOI, n[7] = ne */

   n = dvector(1,7);
   f = dmatrix(1,7,1,1);
   df = dmatrix(1,7,1,7);

   ac = abund(atom,6);
   an = abund(atom,7);
   ao = abund(atom,8);
   ana = abund(atom,11);
   amg = abund(atom,12);
   aal = abund(atom,13);
   asi = abund(atom,14);
   as = abund(atom,16);
   ak  = abund(atom,19);
   aca = abund(atom,20);
   afe = abund(atom,26);
   ati = abund(atom,22);

   if(flagw == 1) printf("Completed Level: ");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      ntot = model->P[i]/(kerg*T);
      n[1] = ntot - model->Ne[i];
      n[2] = ah*n[1] - model->Ne[i];
      n[3] = ahe*n[1];
      n[4] = ac*n[1];
      n[5] = an*n[1];
      n[6] = ao*n[1];
      n[7] = model->Ne[i];
      neorig = Ne = model->Ne[i];
      U = partfn(1.0,T,Ne);


      ph = 4.830e+15*T32*exp(-13.595/kT)/U;
      phm = 2.07e-16*exp(0.7552/kT)/(U*T32);
      ph2 = Xi(101.0,T,Ne);
      phe2 = 4.830e+15*T32*(partfn(2.1,T,Ne)/partfn(2.0,T,Ne))*exp(-24.581/kT);
      phe3 = 4.830e+15*T32*exp(-54.403/kT)/partfn(2.1,T,Ne);
      pc2 = 4.830e+15*T32*(partfn(6.1,T,Ne)/partfn(6.0,T,Ne))*exp(-11.260/kT);
      pc3 = 4.830e+15*T32*(partfn(6.2,T,Ne)/partfn(6.1,T,Ne))*exp(-24.383/kT);
      pc4 = 4.830e+15*T32*(partfn(6.3,T,Ne)/partfn(6.2,T,Ne))*exp(-47.8878/kT);
      pc5 = 4.830e+15*T32*(partfn(6.4,T,Ne)/partfn(6.3,T,Ne))*exp(-64.4939/kT);
      pn2 = 4.830e+15*T32*(partfn(7.1,T,Ne)/partfn(7.0,T,Ne))*exp(-14.534/kT);
      pn3 = 4.830e+15*T32*(partfn(7.2,T,Ne)/partfn(7.1,T,Ne))*exp(-29.601/kT);
      pn4 = 4.830e+15*T32*(partfn(7.3,T,Ne)/partfn(7.2,T,Ne))*exp(-47.4492/kT);
      pn5 = 4.830e+15*T32*(partfn(7.4,T,Ne)/partfn(7.3,T,Ne))*exp(-77.4735/kT);
      po2 = 4.830e+15*T32*(partfn(8.1,T,Ne)/partfn(8.0,T,Ne))*exp(-13.618/kT);
      po3 = 4.830e+15*T32*(partfn(8.2,T,Ne)/partfn(8.1,T,Ne))*exp(-35.116/kT);
      po4 = 4.830e+15*T32*(partfn(8.3,T,Ne)/partfn(8.2,T,Ne))*exp(-54.9355/kT);
      po5 = 4.830e+15*T32*(partfn(8.4,T,Ne)/partfn(8.3,T,Ne))*exp(-77.4135/kT);
      pna2 = 4.830e+15*T32*(partfn(11.1,T,Ne)/partfn(11.0,T,Ne))*exp(-5.139/kT);
      pna3 = 4.830e+15*T32*(partfn(11.2,T,Ne)/partfn(11.1,T,Ne))*exp(-47.286/kT);
      pna4 = 4.830e+15*T32*(partfn(11.3,T,Ne)/partfn(11.2,T,Ne))*exp(-71.62/kT);
      pna5 = 4.830e+15*T32*(partfn(11.4,T,Ne)/partfn(11.3,T,Ne))*exp(-98.91/kT);
      pmg2 = 4.830e+15*T32*(partfn(12.1,T,Ne)/partfn(12.0,T,Ne))*exp(-7.646/kT);
      pmg3 = 4.830e+15*T32*(partfn(12.2,T,Ne)/partfn(12.1,T,Ne))*exp(-15.035/kT);
      pmg4 = 4.830e+15*T32*(partfn(12.3,T,Ne)/partfn(12.2,T,Ne))*exp(-80.1437/kT);
      pmg5 = 4.830e+15*T32*(partfn(12.4,T,Ne)/partfn(12.3,T,Ne))*exp(-109.2655/kT);
      pal2 = 4.830e+15*T32*(partfn(13.1,T,Ne)/partfn(13.0,T,Ne))*exp(-5.986/kT);
      pal3 = 4.830e+15*T32*(partfn(13.2,T,Ne)/partfn(13.1,T,Ne))*exp(-18.828/kT);
      pal4 = 4.830e+15*T32*(partfn(13.3,T,Ne)/partfn(13.2,T,Ne))*exp(-28.4477/kT);
      pal5 = 4.830e+15*T32*(partfn(13.4,T,Ne)/partfn(13.3,T,Ne))*exp(-119.992/kT);
      psi2 = 4.830e+15*T32*(partfn(14.1,T,Ne)/partfn(14.0,T,Ne))*exp(-8.151/kT);
      psi3 = 4.830e+15*T32*(partfn(14.2,T,Ne)/partfn(14.1,T,Ne))*exp(-16.345/kT);
      psi4 = 4.830e+15*T32*(partfn(14.3,T,Ne)/partfn(14.2,T,Ne))*exp(-33.4930/kT);
      psi5 = 4.830e+15*T32*(partfn(14.4,T,Ne)/partfn(14.3,T,Ne))*exp(-45.1418/kT);
      ps2 = 4.830e+15*T32*(partfn(16.1,T,Ne)/partfn(16.0,T,Ne))*exp(-10.36/kT);
      ps3 = 4.830e+15*T32*(partfn(16.2,T,Ne)/partfn(16.1,T,Ne))*exp(-23.3379/kT);
      ps4 = 4.830e+15*T32*(partfn(16.3,T,Ne)/partfn(16.2,T,Ne))*exp(-34.79/kT);
      ps5 = 4.830e+15*T32*(partfn(16.4,T,Ne)/partfn(16.3,T,Ne))*exp(-47.2219/kT);
      pk2 = 4.830e+15*T32*(partfn(19.1,T,Ne)/partfn(19.0,T,Ne))*exp(-4.341/kT);
      pk3 = 4.830e+15*T32*(partfn(19.2,T,Ne)/partfn(19.1,T,Ne))*exp(-31.625/kT);
      pk4 = 4.830e+15*T32*(partfn(19.3,T,Ne)/partfn(19.2,T,Ne))*exp(-45.8060/kT);
      pk5 = 4.830e+15*T32*(partfn(19.4,T,Ne)/partfn(19.3,T,Ne))*exp(-60.9134/kT);
      pca2 = 4.830e+15*T32*(partfn(20.1,T,Ne)/partfn(20.0,T,Ne))*exp(-6.113/kT);
      pca3 = 4.830e+15*T32*(partfn(20.2,T,Ne)/partfn(20.1,T,Ne))*exp(-11.871/kT);
      pca4 = 4.830e+15*T32*(partfn(20.3,T,Ne)/partfn(20.2,T,Ne))*exp(-50.9131/kT);
      pca5 = 4.830e+15*T32*(partfn(20.4,T,Ne)/partfn(20.3,T,Ne))*exp(-67.27/kT);
      pfe2 = 4.830e+15*T32*(partfn(26.1,T,Ne)/partfn(26.0,T,Ne))*exp(-7.870/kT);
      pfe3 = 4.830e+15*T32*(partfn(26.2,T,Ne)/partfn(26.1,T,Ne))*exp(-16.18/kT);
      pfe4 = 4.830e+15*T32*(partfn(26.3,T,Ne)/partfn(26.2,T,Ne))*exp(-30.6514/kT);
      pfe5 = 4.830e+15*T32*(partfn(26.4,T,Ne)/partfn(26.3,T,Ne))*exp(-54.8010/kT);
      pti2 = 4.830e+15*T32*(partfn(22.1,T,Ne)/partfn(22.0,T,Ne))*exp(-6.83/kT);
      pti3 = 4.830e+15*T32*(partfn(22.2,T,Ne)/partfn(22.1,T,Ne))*exp(-13.63/kT);

      while(delta > 0.0001) {
	 ne2 = n[7]*n[7];
	 ne3 = ne2*n[7];
	 ne4 = ne2*ne2;
	 ne5 = ne3*ne2;
	 nnaI = ana*n[1]/(1.0 + pna2/n[7] + pna2*pna3/ne2 + pna2*pna3*pna4/ne3 +
                pna2*pna3*pna4*pna5/ne4);
	 nmgI = amg*n[1]/(1.0 + pmg2/n[7] + pmg2*pmg3/ne2 + pmg2*pmg3*pmg4/ne3 +
                pmg2*pmg3*pmg4*pmg5/ne4);
	 nalI = aal*n[1]/(1.0 + pal2/n[7] + pal2*pal3/ne2 + pal2*pal3*pal4/ne3 +
                pal2*pal3*pal4*pal5/ne4);
	 nsiI = asi*n[1]/(1.0 + psi2/n[7] + psi2*psi3/ne2 + psi2*psi3*psi4/ne3 +
                psi2*psi3*psi4*psi5/ne4);
         nsI = as*n[1]/(1.0 + ps2/n[7] + ps2*ps3/ne2 + ps2*ps3*ps4/ne3 +
                ps2*ps3*ps4*ps5/ne4);
	 nkI = ak*n[1]/(1.0 + pk2/n[7] + pk2*pk3/ne2 + pk2*pk3*pk4/ne3 +
                pk2*pk3*pk4*pk5/ne4);
	 ncaI = aca*n[1]/(1.0 + pca2/n[7] + pca2*pca3/ne2 + pca2*pca3*pca4/ne3 +
                pca2*pca3*pca4*pca5/ne4);
	 nfeI = afe*n[1]/(1.0 + pfe2/n[7] + pfe2*pfe3/ne2 + pfe2*pfe3*pfe4/ne3 +
                pfe2*pfe3*pfe4*pfe5/ne4);
	 ntiI = ati*n[1]/(1.0 + pti2/n[7] + pti2*pti3/ne2);
	 f[1][1] =  n[2] + n[2]*ph/n[7] + n[2]*n[7]*phm + 
                    n[3] + n[3]*phe2/n[7] + n[3]*phe2*phe3/ne2 + 
                    n[4] + n[4]*pc2/n[7] + n[4]*pc2*pc3/ne2 + 
                    n[4]*pc2*pc3*pc4/ne3 + n[4]*pc2*pc3*pc4*pc5/ne4 +
                    n[5] + n[5]*pn2/n[7] + n[5]*pn2*pn3/ne2 +
                    n[5]*pn2*pn3*pn4/ne3 + n[5]*pn2*pn3*pn4*pn5/ne4 + 
                    n[6] + n[6]*po2/n[7] + n[6]*po2*po3/ne2 + 
                    n[6]*po2*po3*po4/ne3 + n[6]*po2*po3*po4*po5/ne4 +
                    n[7] - ntot;
	 f[2][1] =  n[2] + n[2]*ph/n[7] + n[2]*n[7]*phm  - ah*n[1];
	 f[3][1] =  n[3] + n[3]*phe2/n[7] + n[3]*phe2*phe3/ne2 - ahe*n[1];
	 f[4][1] =  n[4] + n[4]*pc2/n[7] + n[4]*pc2*pc3/ne2 +
                    n[4]*pc2*pc3*pc4/ne3 + n[4]*pc2*pc3*pc4*pc5/ne4  
                    - ac*n[1];
	 f[5][1] =  n[5] + n[5]*pn2/n[7] + n[5]*pn2*pn3/ne2 + 
                    n[5]*pn2*pn3*pn4/ne3 + n[5]*pn2*pn3*pn4*pn5/ne4
                    - an*n[1];
	 f[6][1] =  n[6] + n[6]*po2/n[7] + n[6]*po2*po3/ne2 +
                    n[6]*po2*po3*po4/ne3 + n[6]*po2*po3*po4*po5/ne4 
                    - ao*n[1];
	 f[7][1] =  n[2]*ph/n[7] + n[3]*phe2/n[7] + 2.0*n[3]*phe2*phe3/ne2 +
		    n[4]*pc2/n[7] + 2.0*n[4]*pc2*pc3/ne2 +
                    3.0*n[4]*pc2*pc3*pc4/ne3 + 4.0*n[4]*pc2*pc3*pc4*pc5/ne4 +
		    n[5]*pn2/n[7] + 2.0*n[5]*pn2*pn3/ne2 +
	            3.0*n[5]*pn2*pn3*pn4/ne3 + 4.0*n[5]*pn2*pn3*pn4*pn5/ne4 +
		    n[6]*po2/n[7] + 2.0*n[6]*po2*po3/ne2 +
                    3.0*n[6]*po2*po3*po4/ne3 + 4.0*n[6]*po2*po3*po4*po5/ne4 +
		    nnaI*pna2/n[7] + 2.0*nnaI*pna2*pna3/ne2 +
                    3.0*nnaI*pna2*pna3*pna4/ne3 + 4.0*nnaI*pna2*pna3*pna4*pna5/ne4 +
		    nmgI*pmg2/n[7] + 2.0*nmgI*pmg2*pmg3/ne2 + 
                    3.0*nmgI*pmg2*pmg3*pmg4/ne3 + 4.0*nmgI*pmg2*pmg3*pmg4*pmg5/ne4 +
		    nalI*pal2/n[7] + 2.0*nalI*pal2*pal3/ne2 +
                    3.0*nalI*pal2*pal3*pal4/ne3 + 4.0*nalI*pal2*pal3*pal4*pal5/ne4 +
		    nsiI*psi2/n[7] + 2.0*nsiI*psi2*psi3/ne2 +
                    3.0*nsiI*psi2*psi3*psi4/ne3 + 4.0*nsiI*psi2*psi3*psi4*psi5/ne4 +
                    nsI*ps2/n[7] + 2.0*nsI*ps2*ps3/ne2 +
                    3.0*nsI*ps2*ps3*ps4/ne3 + 4.0*nsI*ps2*ps3*ps4*ps5/ne4 +
		    nkI*pk2/n[7] + 2.0*nkI*pk2*pk3/ne2 +
                    3.0*nkI*pk2*pk3*pk4/ne3 + 4.0*nkI*pk2*pk3*pk4*pk5/ne4 +
		    ncaI*pca2/n[7] + 2.0*ncaI*pca2*pca3/ne2 +
                    3.0*ncaI*pca2*pca3*pca4/ne3 + 4.0*ncaI*pca2*pca3*pca4*pca5/ne4 +
		    nfeI*pfe2/n[7] + 2.0*nfeI*pfe2*pfe3/ne2 +
                    3.0*nfeI*pfe2*pfe3*pfe4/ne3 + 4.0*nfeI*pfe2*pfe3*pfe4*pfe5/ne4 - 
                    n[2]*n[7]*phm - n[7];

	 /* completed to here */


	 df[1][1] = 0.0;
	 df[1][2] = 1.0 + ph/n[7] + n[7]*phm;
	 df[1][3] = 1.0 + phe2/n[7] + phe2*phe3/ne2;
	 df[1][4] = 1.0 + pc2/n[7] + pc2*pc3/ne2 + pc2*pc3*pc4/ne3 +
                    pc2*pc3*pc4*pc5/ne4;
	 df[1][5] = 1.0 + pn2/n[7] + pn2*pn3/ne2 + pn2*pn3*pn4/ne3 +
                    pn2*pn3*pn4*pn5/ne4;
	 df[1][6] = 1.0 + po2/n[7] + po2*po3/ne2 + po2*po3*po4/ne3 +
                    po2*po3*po4*po5/ne4;
	 df[1][7] = -n[2]*ph/ne2 + n[2]*phm -
		    n[3]*phe2/ne2 - 2.0*n[3]*phe2*phe3/ne3 -
		    n[4]*pc2/ne2 - 2.0*n[4]*pc2*pc3/ne3 -
                    3.0*n[4]*pc2*pc3*pc4/ne4 - 4.0*pc2*pc3*pc4*pc5/ne5 -
		    n[5]*pn2/ne2 - 2.0*n[5]*pn2*pn3/ne3 -
                    3.0*pn2*pn3*pn4/ne4 - 4.0*pn2*pn3*pn4*pn5/ne5 -
		    n[6]*po2/ne2 - 2.0*n[6]*po2*po3/ne3 -
                    n[6]*po2*po3*po4/ne4 - n[6]*po2*po3*po4*po5/ne5  
                    + 1.0;
	 df[2][1] = -ah;
	 df[2][2] = 1.0 + ph/n[7] + n[7]*phm;
	 df[2][3] = 0.0;
	 df[2][4] = 0.0;
	 df[2][5] = 0.0;
	 df[2][6] = 0.0;
	 df[2][7] = -n[2]*ph/ne2 + n[2]*phm;
	 df[3][1] = -ahe;
	 df[3][2] = 0.0;
	 df[3][3] = 1.0 + phe2/n[7] + phe2*phe3/ne2;
	 df[3][4] = df[3][5] = df[3][6] = 0.0;
	 df[3][7] = -n[3]*phe2/ne2 - 2.0*n[3]*phe2*phe3/ne3;
	 df[4][1] = -ac;
	 df[4][2] = 0.0;
	 df[4][3] = 0.0;
	 df[4][4] = 1.0 + pc2/n[7] + pc2*pc3/ne2 + pc2*pc3*pc4/ne3 +
                    pc2*pc3*pc4*pc5/ne4;
	 df[4][5] = 0.0;
	 df[4][6] = 0.0;
	 df[4][7] = -n[4]*pc2/ne2 - 2.0*n[4]*pc2*pc3/ne3 -
                    3.0*n[4]*pc2*pc3*pc4/ne4 - 4.0*n[4]*pc2*pc3*pc4*pc5/ne5;
	 df[5][1] = -an;
	 df[5][2] = 0.0;
	 df[5][3] = 0.0;
	 df[5][4] = 0.0;
	 df[5][5] = 1.0 + pn2/n[7] + pn2*pn3/ne2 + pn2*pn3*pn4/ne3 +
                    pn2*pn3*pn4*pn5/ne4;
	 df[5][6] = 0.0;
	 df[5][7] = -n[5]*pn2/ne2 - 2.0*n[5]*pn2*pn3/ne3 -
	            3.0*n[5]*pn2*pn3*pn4/ne4 - 4.0*n[5]*pn2*pn3*pn4*pn5/ne5;
	 df[6][1] = -ao;
	 df[6][2] = 0.0;
	 df[6][3] = 0.0;
	 df[6][4] = 0.0;
	 df[6][5] = 0.0;
	 df[6][6] = 1.0 + po2/n[7] + po2*po3/ne2 + po2*po3*po4/ne3 +
                    po2*po3*po4*po5/ne4;
	 df[6][7] = -n[6]*po2/ne2 - 2.0*n[6]*po2*po3/ne3 -
	            3.0*n[6]*po2*po3*po4/ne4 - 4.0*n[6]*po2*po3*po4*po5/ne5;
	 df[7][1] = 0.0;
	 df[7][2] = ph/n[7] - n[7]*phm;
	 df[7][3] = phe2/n[7] + 2.0*phe2*phe3/ne2;
	 df[7][4] = pc2/n[7] + 2.0*pc2*pc3/ne2 + 3.0*pc2*pc3*pc4/ne3 +
                    4.0*pc2*pc3*pc4*pc5/ne4;
	 df[7][5] = pn2/n[7] + 2.0*pn2*pn3/ne2 + 3.0*pn2*pn3*pn4/ne3 +
                    4.0*pn2*pn3*pn4*pn5/ne4;
	 df[7][6] = po2/n[7] + 2.0*po2*po3/ne2 + 3.0*po2*po3*po4/ne3 +
                    4.0*po2*po3*po4*po5/ne4;

	 df[7][7] = -n[2]*ph/ne2 - n[3]*phe2/ne2 -
		    4.0*n[3]*phe2*phe3/ne3 -
		    n[4]*pc2/ne2 - 4.0*n[4]*pc2*pc3/ne3 -
                    9.0*n[4]*pc2*pc3*pc4/ne4 - 16.0*pc2*pc3*pc4*pc5/ne5 -
		    n[5]*pn2/ne2 - 4.0*n[5]*pn2*pn3/ne3 -
                    9.0*n[4]*pn2*pn3*pn4/ne4 - 16.0*pn2*pn3*pn4*pn5/ne5 -
		    n[6]*po2/ne2 - 4.0*n[6]*po2*po3/ne3 -
                    9.0*n[4]*po2*po3*po4/ne4 - 16.0*po2*po3*po4*po5/ne5 -
		    nnaI*pna2/ne2 - 4.0*nnaI*pna2*pna3/ne3 -
                    9.0*nnaI*pna2*pna3*pna4/ne4 - 16.0*nnaI*pna2*pna3*pna4*pna5/ne5 -
		    nmgI*pmg2/ne2 - 4.0*nmgI*pmg2*pmg3/ne3 -
                    9.0*nmgI*pmg2*pmg3*pmg4/ne4 - 16.0*nmgI*pmg2*pmg3*pmg4*pmg5/ne5 -
		    nalI*pal2/ne2 - 4.0*nalI*pal2*pal3/ne3 -
                    9.0*nalI*pal2*pal3*pal4/ne4 - 16.0*nalI*pal2*pal3*pal4*pal5/ne5 -
		    nsiI*psi2/ne2 - 4.0*nsiI*psi2*psi3/ne3 -
                    9.0*nsiI*psi2*psi3*psi4/ne4 - 16.0*nsiI*psi2*psi3*psi4*psi5/ne5 -
                    nsI*ps2/ne2 - 4.0*nsI*ps2*ps3/ne3 -
                    9.0*nsI*ps2*ps3*ps4/ne4 - 16.0*nsI*ps2*ps3*ps4*ps5/ne5 -
		    nkI*pk2/ne2 - 4.0*nkI*pk2*pk3/ne3 -
                    9.0*nkI*pk2*pk3*pk4/ne4 - 16.0*nkI*pk2*pk3*pk4*pk5/ne5 -
		    ncaI*pca2/ne2 - 4.0*ncaI*pca2*pca3/ne3 -
                    9.0*ncaI*pca2*pca3*pca4/ne4 - 16.0*ncaI*pca2*pca3*pca4*pca5/ne5 -
		    nfeI*pfe2/ne2 - 4.0*nfeI*pfe2*pfe3/ne3 -
                    9.0*nfeI*pfe2*pfe3*pfe4/ne4 - 16.0*nfeI*pfe2*pfe3*pfe4*pfe5/ne5 -
		    n[2]*phm - 1.0;

	 gaussj(df,7,f,1);
	 delta = 0.0;

	 for(j=1;j<=7;j++) {
	    n[j] -= f[j][1];
	    delta += fabs(f[j][1]/n[j]);
	 }
	 delta /= 7.0;
      }
      delta = 1.0;
      model->Ne[i] = n[7];
      model->NHI[i] = n[2];
      model->N1[i] = 2.0*model->NHI[i]/U;
      model->N2[i] = 8.0*model->NHI[i]*exp(-1.1835e+5/T)/U;
      model->N3[i] = 18.0*model->NHI[i]*exp(-1.4027e+5/T)/U;
      model->NHminus[i] = model->NHI[i]*model->Ne[i]*phm;
      model->NH2[i] = model->NHI[i]*model->NHI[i]*ph2;
      model->NHeI[i] = n[3];
      model->NA[i] = n[1];
      model->Np[i] = n[2]*ph/n[7];
      model->NHeII[i] = n[3]*phe2/n[7];
      model->NHeIII[i] = n[3]*phe2*phe3/(n[7]*n[7]);
      model->NCI[i] = n[4];
      model->NNI[i] = n[5];
      model->NOI[i] = n[6];
      model->NMgI[i] = nmgI;
      model->NMgII[i] = nmgI*pmg2/n[7];
      model->NAlI[i] = nalI;
      model->NSiI[i] = nsiI;
      model->NSiII[i] = nsiI*psi2/n[7];
      model->NCaI[i] = ncaI;
      model->NCaII[i] = ncaI*pca2/n[7];
      model->NFeI[i] = nfeI;
      model->NTiI[i] = ntiI;
      model->NCH[i] = 0.0;
      model->NOH[i] = 0.0;
      model->NMgH[i] = 0.0;
      model->NH2O[i] = 0.0;
      model->rho[i] = 1.66054e-24*(1.008*ah + 4.003*ahe + 12.01*ac + 14.01*an +
	    16.00*ao + 23.00*ana + 24.32*amg + 26.97*aal + 28.06*asi +
	    39.10*ak + 40.08*aca + 55.85*afe)*model->NA[i];
      model->nmax[i] = (int)ceil(2.324*pow(model->T[i]/
			    (kerg*model->Ne[i]),0.25));
      model->nmax[i] = (int)imin(40,model->nmax[i]);
      if(flagw == 1) printf(" %d ",i);
   }
   free_dmatrix(f,1,7,1,1);
   free_dmatrix(df,1,7,1,7);
   free_dvector(n,1,7);
   printf("\n");
   if(flagP == 1) printhotDensity(model);
}

void printhotDensity(model)
     atmosphere *model;
{
  FILE *out;
  double T,Ne,NA,U,P,ntot,T32,kT,ne2,ne3,rD,DE;
  double C,N,O,H;
  double phe2,phe3,ph,phm;
  double pc2,pc3,pc4;
  double pn2,pn3,pn4,po2,po3,po4;
  double pca2,pca3,pfe2,pfe3,pti2,pti3;  
  double k = 8.617084e-05;
  double kerg = 1.38054e-16;
  int i;
  extern int Ntau;

  out = fopen("density.out","w");

  Ne = model->Ne[i];

  fprintf(out,"n  N(e-)     N(H)      N(H+)     N(H-)     N(He)     N(He+)    N(He++)\n");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      U = partfn(1.0,T,Ne);

      ph = 4.830e+15*T32*exp(-13.595/kT)/U;
      phm = 2.07e-16*exp(0.7552/kT)/(U*T32);
      phe2 = 4.830e+15*T32*(partfn(2.1,T,Ne)/partfn(2.0,T,Ne))*exp(-24.581/kT);
      phe3 = 4.830e+15*T32*exp(-54.403/kT)/partfn(2.1,T,Ne);
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              model->Ne[i],model->NHI[i],model->NHI[i]*ph/model->Ne[i],model->NHminus[i],
              model->NHeI[i],model->NHeI[i]*phe2/model->Ne[i],model->NHeI[i]*phe2*phe3/(model->Ne[i]*model->Ne[i]));
   }

   fprintf(out,"\n");
   fprintf(out,"n  N(C)      N(C+)     N(C++)    N(C+3)   N(N)      N(N+)     N(N++)    N(N+3)    N(O)      N(O+)     N(O++)     N(O+3)\n");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      U = partfn(1.0,T,Ne);
      Ne = model->Ne[i];
      ne2 = Ne*Ne;
      ne3 = ne2*Ne;
      pc2 = 4.830e+15*T32*(partfn(6.1,T,Ne)/partfn(6.0,T,Ne))*exp(-11.260/kT);
      pc3 = 4.830e+15*T32*(partfn(6.2,T,Ne)/partfn(6.1,T,Ne))*exp(-24.383/kT);
      pc4 = 4.830e+15*T32*(partfn(6.3,T,Ne)/partfn(6.2,T,Ne))*exp(-47.8878/kT);
      pn2 = 4.830e+15*T32*(partfn(7.1,T,Ne)/partfn(7.0,T,Ne))*exp(-14.534/kT);
      pn3 = 4.830e+15*T32*(partfn(7.2,T,Ne)/partfn(7.1,T,Ne))*exp(-29.601/kT);
      pn4 = 4.830e+15*T32*(partfn(7.3,T,Ne)/partfn(7.2,T,Ne))*exp(-47.4492/kT);
      po2 = 4.830e+15*T32*(partfn(8.1,T,Ne)/partfn(8.0,T,Ne))*exp(-13.618/kT);
      po3 = 4.830e+15*T32*(partfn(8.2,T,Ne)/partfn(8.1,T,Ne))*exp(-35.116/kT);
      po4 = 4.830e+15*T32*(partfn(8.3,T,Ne)/partfn(8.2,T,Ne))*exp(-54.9355/kT);
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              model->NCI[i],model->NCI[i]*pc2/Ne,model->NCI[i]*pc2*pc3/ne2,model->NCI[i]*pc2*pc3*pc4/ne3,
              model->NNI[i],model->NNI[i]*pn2/Ne,model->NNI[i]*pn2*pn3/ne2,model->NNI[i]*pn2*pn3*pn4/ne3,
              model->NOI[i],model->NOI[i]*po2/Ne,model->NOI[i]*po2*po3/ne2,model->NOI[i]*po2*po3*po4/ne3);
   }

   fprintf(out,"\n");

   fprintf(out,"n  N(Fe)     N(Fe+)    N(Fe++)   N(Ti)     N(Ti+)    N(Ti++)   N(H(n=3))\n");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      U = partfn(1.0,T,Ne);
      Ne = model->Ne[i];
      ne2 = Ne*Ne;
      pfe2 = 4.830e+15*T32*(partfn(26.1,T,Ne)/partfn(26.0,T,Ne))*exp(-7.870/kT);
      pfe3 = 4.830e+15*T32*(partfn(26.2,T,Ne)/partfn(26.1,T,Ne))*exp(-16.18/kT);
      pti2 = 4.830e+15*T32*(partfn(22.1,T,Ne)/partfn(22.0,T,Ne))*exp(-6.83/kT);
      pti3 = 4.830e+15*T32*(partfn(22.2,T,Ne)/partfn(22.1,T,Ne))*exp(-13.63/kT);
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              model->NFeI[i],model->NFeI[i]*pfe2/Ne,model->NFeI[i]*pfe2*pfe3/ne2,
              model->NTiI[i],model->NTiI[i]*pti2/Ne,model->NTiI[i]*pti2*pti3/ne2,
              model->N3[i]);
   }

   fprintf(out,"\n");
   fclose(out);
}
