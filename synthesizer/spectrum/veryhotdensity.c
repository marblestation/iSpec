#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "spectrum.h"
double partfn(double code,double T,double Ne);
double abund();
int imin();
void printhotDensity();

/* This routine assumes the electron density given in the model is correct.
   There is no iteration.  Use for temperatures > 25,000K */

void veryhotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw)
{
   double T,Ne,NA,U,P,ntot,T32,kT,ne2,ne3,ne4,ne5,rD,DE;
   double phe2,phe3,ph,ph2p,ph2,phm,ps2,ps3,ps4,ps5;
   double pc2,pc3,pc4,pc5,pna2,pna3,pna4,pna5,pmg2,pmg3,pmg4,pmg5;
   double pal2,pal3,pal4,pal5,psi2,psi3,psi4,psi5,pk2,pk3,pk4,pk5;
   double pca2,pca3,pca4,pca5,pfe2,pfe3,pfe4,pfe5,pn2,pn3,pn4,pn5;
   double po2,po3,po4,po5,pti2,pti3;
   double ncI,nnI,noI,nnaI,nmgI,nalI,nsiI,nsI,nkI,ncaI,nfeI,ntiI;
   double nhI,nheI,nA;
   double ac,an,ao,ana,amg,aal,asi,as,ak,aca,afe,ati;
   double k = 8.617084e-05;
   double kerg = 1.38054e-16;
   double delta = 1.0;
   double **df,**f,*n;
   double neorig;
   double rho;
   int i,j;
   extern int Ntau;
      extern int flagP;

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

   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;     
      ntot = model->P[i]/(kerg*T);
      nA = ntot - model->Ne[i];
      Ne = model->Ne[i];
      ne2 = Ne*Ne;
      ne3 = ne2*Ne;
      ne4 = ne2*ne2;
      ne5 = ne2*ne3;
      U = partfn(1.0,T,Ne);

      ph = 4.830e+15*T32*exp(-13.595/kT)/U;
      phm = 2.07e-16*exp(0.7552/kT)/(U*T32);
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

      nhI = ah*nA/(1.0 + ph/Ne);
      nheI = ahe*nA/(1.0 + phe2/Ne + phe2*phe3/ne2);
      ncI = ac*nA/(1.0 + pc2/Ne + pc2*pc3/ne2 + pc2*pc3*pc4/ne3 +
                pc2*pc3*pc4*pc5/ne4);
      nnI = an*nA/(1.0 + pn2/Ne + pn2*pn3/ne2 + pn2*pn3*pn4/ne3 +
                pn2*pn3*pn4*pn5/ne4);
      noI = ao*nA/(1.0 + po2/Ne + po2*po3/ne2 + po2*po3*po4/ne3 +
                po2*po3*po4*po5/ne4);
      nnaI = ana*nA/(1.0 + pna2/Ne + pna2*pna3/ne2 + pna2*pna3*pna4/ne3 +
                pna2*pna3*pna4*pna5/ne4);
      nmgI = amg*nA/(1.0 + pmg2/Ne + pmg2*pmg3/ne2 + pmg2*pmg3*pmg4/ne3 +
                pmg2*pmg3*pmg4*pmg5/ne4);
      nalI = aal*nA/(1.0 + pal2/Ne + pal2*pal3/ne2 + pal2*pal3*pal4/ne3 +
                pal2*pal3*pal4*pal5/ne4);
      nsiI = asi*nA/(1.0 + psi2/Ne + psi2*psi3/ne2 + psi2*psi3*psi4/ne3 +
                psi2*psi3*psi4*psi5/ne4);
      nsI = as*nA/(1.0 + ps2/Ne + ps2*ps3/ne2 + ps2*ps3*ps4/ne3 +
                ps2*ps3*ps4*ps5/ne4);
      nkI = ak*nA/(1.0 + pk2/Ne + pk2*pk3/ne2 + pk2*pk3*pk4/ne3 +
                pk2*pk3*pk4*pk5/ne4);
      ncaI = aca*nA/(1.0 + pca2/Ne + pca2*pca3/ne2 + pca2*pca3*pca4/ne3 +
                pca2*pca3*pca4*pca5/ne4);
      nfeI = afe*nA/(1.0 + pfe2/Ne + pfe2*pfe3/ne2 + pfe2*pfe3*pfe4/ne3 +
                pfe2*pfe3*pfe4*pfe5/ne4);
      ntiI = ati*nA/(1.0 + pti2/Ne + pti2*pti3/ne2);

      model->NHI[i] = nhI;
      model->N1[i] = 2.0*model->NHI[i]/U;
      model->N2[i] = 8.0*model->NHI[i]*exp(-1.1835e+5/T)/U;
      model->N3[i] = 18.0*model->NHI[i]*exp(-1.4027e+5/T)/U;
      model->NHminus[i] = model->NHI[i]*model->Ne[i]*phm;
      model->NH2[i] = 0.0;
      model->NHeI[i] = nheI;
      model->NA[i] = nA;
      model->Np[i] = nhI*ph/Ne;
      model->NHeII[i] = nheI*phe2/Ne;
      model->NHeIII[i] = nheI*phe2*phe3/ne2;
      model->NCI[i] = ncI;
      model->NNI[i] = nnI;
      model->NOI[i] = noI;
      model->NMgI[i] = nmgI;
      model->NMgII[i] = nmgI*pmg2/Ne;
      model->NAlI[i] = nalI;
      model->NSiI[i] = nsiI;
      model->NSiII[i] = nsiI*psi2/Ne;
      model->NCaI[i] = ncaI;
      model->NCaII[i] = ncaI*pca2/Ne;
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
   printf("\n");
   if(flagP == 1) printhotDensity(model);
}
