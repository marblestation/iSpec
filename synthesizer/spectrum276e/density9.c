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
void printDensity();
/* Jan 19, 1997: Common metal hydrides, such as SiH, MgH, CaH, as well
as TiO have been added to the equilibrium computations for these elements,
but Si, Mg and Ca have not yet been added to the equilibrium computations
for CNO */
/* Aug 18, 2008: Equilibrium calculation including AlH and AlO added */
/* March 17, 2009: Equilibrium calculation including ZrO added */


void Density(model,atom,ah,ahe,flagw)
atmosphere *model;
atominfo *atom;
double ah,ahe;
int flagw;
{
   double T,Ne,NA,U,P,ntot,T32,kT,ne2,ne3,rD,DE;
   double phe2,phe3,ph,ph2p,ph2,phm;
   double pc2,pc3,pna2,pna3,pmg2,pmg3,pal2,pal3,psi2,psi3,pk2,pk3;
   double pca2,pca3,pfe2,pfe3,pn2,pn3,po2,po3,pti2,pti3,pzr2,pzr3;
   double ps2,ps3;
   double ECH,EC2,ECN,ECO,ENO,EN2,EO2,EOH,ENH;
   double EMgH,EMgO,ESiH,ESiC,ESiO,ECaH,ECaO,ETiN,ETiO,EZrN,EZrO,EAlH,EAlO;
   double EH2O;
   double ncI,nnI,noI,nnaI,nmgI,nalI,nsiI,nsI,nkI,ncaI,nfeI,ntiI,nzrI;
   double ac,an,ao,ana,amg,aal,asi,ak,aca,afe,ati,as,azr;
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
   azr = abund(atom,40);

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
      ph2 = Xi(101.0,T,Ne);
      ph2p = 2.525e-05*exp(-10.945/kT)/(U*U);
      phm = 2.07e-16*exp(0.7552/kT)/(U*T32);
      phe2 = 4.830e+15*T32*(partfn(2.1,T,Ne)/partfn(2.0,T,Ne))*exp(-24.581/kT);
      phe3 = 4.830e+15*T32*exp(-54.403/kT)/partfn(2.1,T,Ne);
      pc2 = 4.830e+15*T32*(partfn(6.1,T,Ne)/partfn(6.0,T,Ne))*exp(-11.260/kT);
      pc3 = 4.830e+15*T32*(partfn(6.2,T,Ne)/partfn(6.1,T,Ne))*exp(-24.383/kT);
      pn2 = 4.830e+15*T32*(partfn(7.1,T,Ne)/partfn(7.0,T,Ne))*exp(-14.534/kT);
      pn3 = 4.830e+15*T32*(partfn(7.2,T,Ne)/partfn(7.1,T,Ne))*exp(-29.601/kT);
      po2 = 4.830e+15*T32*(partfn(8.1,T,Ne)/partfn(8.0,T,Ne))*exp(-13.618/kT);
      po3 = 4.830e+15*T32*(partfn(8.2,T,Ne)/partfn(8.1,T,Ne))*exp(-35.116/kT);
      pna2 = 4.830e+15*T32*(partfn(11.1,T,Ne)/partfn(11.0,T,Ne))*exp(-5.139/kT);
      pna3 = 4.830e+15*T32*(partfn(11.2,T,Ne)/partfn(11.1,T,Ne))*exp(-47.286/kT);
      pmg2 = 4.830e+15*T32*(partfn(12.1,T,Ne)/partfn(12.0,T,Ne))*exp(-7.646/kT);
      pmg3 = 4.830e+15*T32*(partfn(12.2,T,Ne)/partfn(12.1,T,Ne))*exp(-15.035/kT);
      pal2 = 4.830e+15*T32*(partfn(13.1,T,Ne)/partfn(13.0,T,Ne))*exp(-5.986/kT);
      pal3 = 4.830e+15*T32*(partfn(13.2,T,Ne)/partfn(13.1,T,Ne))*exp(-18.828/kT);
      psi2 = 4.830e+15*T32*(partfn(14.1,T,Ne)/partfn(14.0,T,Ne))*exp(-8.151/kT);
      psi3 = 4.830e+15*T32*(partfn(14.2,T,Ne)/partfn(14.1,T,Ne))*exp(-16.345/kT);
      ps2 = 4.830e+15*T32*(partfn(16.1,T,Ne)/partfn(16.0,T,Ne))*exp(-10.3600/kT);
      ps3 = 4.830e+15*T32*(partfn(16.2,T,Ne)/partfn(16.1,T,Ne))*exp(-23.3379/kT);
      pk2 = 4.830e+15*T32*(partfn(19.1,T,Ne)/partfn(19.0,T,Ne))*exp(-4.341/kT);
      pk3 = 4.830e+15*T32*(partfn(19.2,T,Ne)/partfn(19.1,T,Ne))*exp(-31.625/kT);
      pca2 = 4.830e+15*T32*(partfn(20.1,T,Ne)/partfn(20.0,T,Ne))*exp(-6.113/kT);
      pca3 = 4.830e+15*T32*(partfn(20.2,T,Ne)/partfn(20.1,T,Ne))*exp(-11.871/kT);
      pfe2 = 4.830e+15*T32*(partfn(26.1,T,Ne)/partfn(26.0,T,Ne))*exp(-7.870/kT);
      pfe3 = 4.830e+15*T32*(partfn(26.2,T,Ne)/partfn(26.1,T,Ne))*exp(-16.18/kT);
      pti2 = 4.830e+15*T32*(partfn(22.1,T,Ne)/partfn(22.0,T,Ne))*exp(-6.83/kT);
      pti3 = 4.830e+15*T32*(partfn(22.2,T,Ne)/partfn(22.1,T,Ne))*exp(-13.63/kT);
      pzr2 = 4.830e+15*T32*(partfn(40.1,T,Ne)/partfn(40.0,T,Ne))*exp(-6.76/kT);
      pzr3 = 4.830e+15*T32*(partfn(40.2,T,Ne)/partfn(40.1,T,Ne))*exp(-14.0/kT);
      ECH = Xi(106.0,T,Ne);
      ENH = Xi(107.0,T,Ne);
      EOH = Xi(108.0,T,Ne);
      EC2 = Xi(606.0,T,Ne);
      ECN = Xi(607.0,T,Ne);
      ECO = Xi(608.0,T,Ne);
      EN2 = Xi(707.0,T,Ne);
      ENO = Xi(708.0,T,Ne);
      EO2 = Xi(808.0,T,Ne);
      if(model->teff < 6500.0) {
	EMgH = Xi(112.0,T,Ne);
	EMgO = Xi(812.0,T,Ne);
	EAlH = Xi(113.0,T,Ne);
	EAlO = Xi(813.0,T,Ne);
	ESiH = Xi(114.0,T,Ne);
	ESiC = Xi(614.0,T,Ne);
	ESiO = Xi(814.0,T,Ne);
	ECaH = Xi(120.0,T,Ne);
	ECaO = Xi(820.0,T,Ne);
	ETiN = Xi(722.0,T,Ne);
	ETiO = Xi(822.0,T,Ne);
	EZrO = Xi(840.0,T,Ne);
	EH2O = Xi(10108.0,T,Ne);
	flagmgh = 1;
      } else EMgH = EMgO = EAlH = EAlO = ESiH = ESiC = ESiO = ECaH = ECaO = 
             ETiN = ETiO = EZrO = EH2O = 0.0;

      while(delta > 0.0001) {
	 ne2 = n[7]*n[7];
	 ne3 = n[7]*n[7]*n[7];
	 nnaI = ana*n[1]/(1.0 + pna2/n[7] + pna2*pna3/ne2);
	 if(model->teff >= 6500.0)
	      nmgI = amg*n[1]/(1.0 + pmg2/n[7] + pmg2*pmg3/ne2);
	 else nmgI = amg*n[1]/(1.0 + pmg2/n[7] + pmg2*pmg3/ne2 + EMgH*n[2] +
		     EMgO*n[6]);
	 if(model->teff >= 6500) 
              nalI = aal*n[1]/(1.0 + pal2/n[7] + pal2*pal3/ne2);
	 else nalI = aal*n[1]/(1.0 + pal2/n[7] + pal2*pal3/ne2 + EAlH*n[2] +
		     EAlO*n[6]);
	 if(model->teff >= 6500.0)
	      nsiI = asi*n[1]/(1.0 + psi2/n[7] + psi2*psi3/ne2);
	 else nsiI = asi*n[1]/(1.0 + psi2/n[7] + psi2*psi3/ne2 + ESiH*n[2] +
		     ESiC*n[4] + ESiO*n[6]);
	 nsI = as*n[1]/(1.0 + ps2/n[7] + ps2*ps3/ne2);
	 nkI = ak*n[1]/(1.0 + pk2/n[7] + pk2*pk3/ne2);
	 if(model->teff >= 6500.0)
	      ncaI = aca*n[1]/(1.0 + pca2/n[7] + pca2*pca3/ne2);
	 else ncaI = aca*n[1]/(1.0 + pca2/n[7] + pca2*pca3/ne2 + ECaH*n[2] +
		     ECaO*n[6]);
	 nfeI = afe*n[1]/(1.0 + pfe2/n[7] + pfe2*pfe3/ne2);
	 if(model->teff >= 6500.0)
	      ntiI = ati*n[1]/(1.0 + pti2/n[7] + pti2*pti3/ne2);
	 else ntiI = ati*n[1]/(1.0 + pti2/n[7] + pti2*pti3/ne2 + ETiN*n[5] +
		     ETiO*n[6]);
	 if(model->teff >= 6500.0)
	      nzrI = azr*n[1]/(1.0 + pzr2/n[7] + pzr2*pzr3/ne2);
	 else nzrI = azr*n[1]/(1.0 + pzr2/n[7] + pzr2*pzr3/ne2 + EZrO*n[6]);
	 f[1][1] =  n[2] + n[2]*ph/n[7] + n[2]*n[2]*ph2 +
		    n[2]*n[2]*ph2p/n[7] + n[2]*n[7]*phm + n[3] +
		    n[3]*phe2/n[7] + n[3]*phe2*phe3/ne2 + n[4] +
		    n[4]*pc2/n[7] + n[4]*pc2*pc3/ne2 + n[5] + n[5]*pn2/n[7] +
		    n[5]*pn2*pn3/ne2 + n[6] + n[6]*po2/n[7] +
		    n[6]*po2*po3/ne2 + n[2]*n[4]*ECH + n[4]*n[4]*EC2 +
		    n[4]*n[5]*ECN + n[4]*n[6]*ECO + n[5]*n[2]*ENH +
		    n[5]*n[6]*ENO + n[5]*n[5]*EN2 + n[6]*n[2]*EOH +
		    n[6]*n[6]*EO2 + n[2]*n[6]*n[2]*EH2O + n[7] - ntot;
	 f[2][1] =  n[2] + n[2]*ph/n[7] + 2.0*n[2]*n[2]*ph2 +
		    2.0*n[2]*n[2]*ph2p/n[7] + n[2]*n[7]*phm + n[2]*n[4]*ECH +
		    n[2]*n[5]*ENH + n[2]*n[6]*EOH + 2.0*n[2]*n[6]*n[2]*EH2O - 
                    ah*n[1];
	 f[3][1] =  n[3] + n[3]*phe2/n[7] + n[3]*phe2*phe3/ne2 - ahe*n[1];
	 f[4][1] =  n[4] + n[4]*pc2/n[7] + n[4]*pc2*pc3/ne2 + n[4]*n[2]*ECH +
		    2.0*n[4]*n[4]*EC2 + n[4]*n[5]*ECN + n[4]*n[6]*ECO -
		    ac*n[1];
	 f[5][1] =  n[5] + n[5]*pn2/n[7] + n[5]*pn2*pn3/ne2 + n[5]*n[2]*ENH +
		    2.0*n[5]*n[5]*EN2 + n[4]*n[5]*ECN + n[5]*n[6]*ENO -
		    an*n[1];
	 f[6][1] =  n[6] + n[6]*po2/n[7] + n[6]*po2*po3/ne2 + n[6]*n[2]*EOH +
		    2.0*n[6]*n[6]*EO2 + n[4]*n[6]*ECO + n[5]*n[6]*ENO + 
                    n[2]*n[6]*n[2]*EH2O - ao*n[1];
	 f[7][1] =  n[2]*ph/n[7] + n[2]*n[2]*ph2p/n[7] +
		    n[3]*phe2/n[7] + 2.0*n[3]*phe2*phe3/ne2 +
		    n[4]*pc2/n[7] + 2.0*n[4]*pc2*pc3/ne2 +
		    n[5]*pn2/n[7] + 2.0*n[5]*pn2*pn3/ne2 +
		    n[6]*po2/n[7] + 2.0*n[6]*po2*po3/ne2 +
		    nnaI*pna2/n[7] + 2.0*nnaI*pna2*pna3/ne2 +
		    nmgI*pmg2/n[7] + 2.0*nmgI*pmg2*pmg3/ne2 +
		    nalI*pal2/n[7] + 2.0*nalI*pal2*pal3/ne2 +
		    nsiI*psi2/n[7] + 2.0*nsiI*psi2*psi3/ne2 +
                    nsI*ps2/n[7] + 2.0*nsI*ps2*ps3/ne2 +
		    nkI*pk2/n[7] + 2.0*nkI*pk2*pk3/ne2 +
		    ncaI*pca2/n[7] + 2.0*ncaI*pca2*pca3/ne2 +
		    nfeI*pfe2/n[7] + 2.0*nfeI*pfe2*pfe3/ne2 -
		    n[2]*n[7]*phm - n[7];

	 df[1][1] = 0.0;
	 df[1][2] = 1.0 + ph/n[7] + 2.0*n[2]*ph2 + 2.0*n[2]*ph2p/n[7] +
		    n[7]*phm + n[4]*ECH + n[5]*ENH + n[6]*EOH + 
                    2.0*n[2]*n[6]*EH2O;
	 df[1][3] = 1.0 + phe2/n[7] + phe2*phe3/ne2;
	 df[1][4] = 1.0 + pc2/n[7] + pc2*pc3/ne2 + n[2]*ECH + 2.0*n[4]*EC2 +
		    n[5]*ECN + n[6]*ECO;
	 df[1][5] = 1.0 + pn2/n[7] + pn2*pn3/ne2 + n[4]*ECN + n[2]*ENH +
		    n[6]*ENO + 2.0*n[5]*EN2;
	 df[1][6] = 1.0 + po2/n[7] + po2*po3/ne2 + n[4]*ECO + n[5]*ENO +
		    n[2]*EOH + 2.0*n[6]*EO2 + n[2]*n[2]*EH2O;
	 df[1][7] = -n[2]*ph/ne2 - n[2]*n[2]*ph2p/ne2 + n[2]*phm -
		    n[3]*phe2/ne2 - 2.0*n[3]*phe2*phe3/ne3 -
		    n[4]*pc2/ne2 - 2.0*n[4]*pc2*pc3/ne3 -
		    n[5]*pn2/ne2 - 2.0*n[5]*pn2*pn3/ne3 -
		    n[6]*po2/ne2 - 2.0*n[6]*po2*po3/ne3 + 1.0;
	 df[2][1] = -ah;
	 df[2][2] = 1.0 + ph/n[7] + 4.0*n[2]*ph2 + 4.0*n[2]*ph2p/n[7] +
		    n[7]*phm + n[4]*ECH + n[5]*ENH + n[6]*EOH + 
                    4.0*n[2]*n[6]*EH2O;
	 df[2][3] = 0.0;
	 df[2][4] = n[2]*ECH;
	 df[2][5] = n[2]*ENH;
	 df[2][6] = n[2]*EOH + 2.0*n[2]*n[2]*EH2O;
	 df[2][7] = -n[2]*ph/ne2 - 2.0*n[2]*n[2]*ph2p/ne2 + n[2]*phm;
	 df[3][1] = -ahe;
	 df[3][2] = 0.0;
	 df[3][3] = 1.0 + phe2/n[7] + phe2*phe3/ne2;
	 df[3][4] = df[3][5] = df[3][6] = 0.0;
	 df[3][7] = -n[3]*phe2/ne2 - 2.0*n[3]*phe2*phe3/ne3;
	 df[4][1] = -ac;
	 df[4][2] = n[4]*ECH;
	 df[4][3] = 0.0;
	 df[4][4] = 1.0 + pc2/n[7] + pc2*pc3/ne2 + n[2]*ECH + 4.0*n[4]*EC2 +
		    n[5]*ECN + n[6]*ECO;
	 df[4][5] = n[4]*ECN;
	 df[4][6] = n[4]*ECO;
	 df[4][7] = -n[4]*pc2/ne2 - 2.0*n[4]*pc2*pc3/ne3;
	 df[5][1] = -an;
	 df[5][2] = n[5]*ENH;
	 df[5][3] = 0.0;
	 df[5][4] = n[5]*ECN;
	 df[5][5] = 1.0 + pn2/n[7] + pn2*pn3/ne2 + n[2]*ENH + 4.0*n[5]*EN2 +
		    n[4]*ECN + n[6]*ENO;
	 df[5][6] = n[5]*ENO;
	 df[5][7] = -n[5]*pn2/ne2 - 2.0*n[5]*pn2*pn3/ne3;
	 df[6][1] = -ao;
	 df[6][2] = n[6]*EOH + 2.0*n[2]*n[6]*EH2O;
	 df[6][3] = 0.0;
	 df[6][4] = n[6]*ECO;
	 df[6][5] = n[6]*ENO;
	 df[6][6] = 1.0 + po2/n[7] + po2*po3/ne2 + n[2]*EOH + 4.0*n[6]*EO2 +
		    n[4]*ECO + n[5]*ENO + n[2]*n[2]*EH2O;
	 df[6][7] = -n[6]*po2/ne2 - 2.0*n[6]*po2*po3/ne3;
	 df[7][1] = 0.0;
	 df[7][2] = ph/n[7] + 2.0*n[2]*ph2p/n[7] - n[7]*phm;
	 df[7][3] = phe2/n[7] + 2.0*phe2*phe3/ne2;
	 df[7][4] = pc2/n[7] + 2.0*pc2*pc3/ne2;
	 df[7][5] = pn2/n[7] + 2.0*pn2*pn3/ne2;
	 df[7][6] = po2/n[7] + 2.0*po2*po3/ne2;
	 df[7][7] = -n[2]*ph/ne2 - n[2]*n[2]*ph2p/ne2 - n[3]*phe2/ne2 -
		    4.0*n[3]*phe2*phe3/ne3 -
		    n[4]*pc2/ne2 - 4.0*n[4]*pc2*pc3/ne3 -
		    n[5]*pn2/ne2 - 4.0*n[5]*pn2*pn3/ne3 -
		    n[6]*po2/ne2 - 4.0*n[6]*po2*po3/ne3 -
		    nnaI*pna2/ne2 - 4.0*nnaI*pna2*pna3/ne3 -
		    nmgI*pmg2/ne2 - 4.0*nmgI*pmg2*pmg3/ne3 -
		    nalI*pal2/ne2 - 4.0*nalI*pal2*pal3/ne3 -
		    nsiI*psi2/ne2 - 4.0*nsiI*psi2*psi3/ne3 -
                    nsI*ps2/ne2 - 4.0*nsI*ps2*ps3/ne3 -
		    nkI*pk2/ne2 - 4.0*nkI*pk2*pk3/ne3 -
		    ncaI*pca2/ne2 - 4.0*ncaI*pca2*pca3/ne3 -
		    nfeI*pfe2/ne2 - 4.0*nfeI*pfe2*pfe3/ne3 -
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
      model->NZrI[i] = nzrI;
      model->NCH[i] = n[2]*n[4]*ECH;
      model->NOH[i] = n[2]*n[6]*EOH;
      model->NMgH[i] = nmgI*n[2]*EMgH;
      model->NH2O[i] = n[2]*n[6]*n[2]*EH2O;
      model->rho[i] = 1.66054e-24*(1.008*ah + 4.003*ahe + 12.01*ac + 14.01*an +
	    16.00*ao + 23.00*ana + 24.32*amg + 26.97*aal + 28.06*asi +
	    32.065*as + 39.10*ak + 40.08*aca + 55.85*afe)*model->NA[i];
      model->nmax[i] = (int)ceil(2.324*pow(model->T[i]/
			    (kerg*model->Ne[i]),0.25));
      model->nmax[i] = (int)imin(40,model->nmax[i]);
      if(flagw == 1) printf(" %d ",i);
   }
   free_dmatrix(f,1,7,1,1);
   free_dmatrix(df,1,7,1,7);
   free_dvector(n,1,7);
   printf("\n");
   if(flagP == 1) printDensity(model);
}



void printDensity(model)
     atmosphere *model;
{
  FILE *out;
  double T,Ne,NA,U,P,ntot,T32,kT,ne2,ne3,rD,DE;
  double C,N,O,H;
  double phe2,phe3,ph,ph2p,ph2,phm;
  double pc2,pc3,pna2,pna3,pmg2,pmg3,pal2,pal3,psi2,psi3,pk2,pk3;
  double pca2,pca3,pfe2,pfe3,pn2,pn3,po2,po3,pti2,pti3;
  double ECH,EC2,ECN,ECO,ENO,EN2,EO2,EOH,ENH,EH2O;
  double EMgH,EMgO,ESiH,ESiC,ESiO,ECaH,ECaO,ETiN,ETiO,EZrN,EZrO;
  double ncI,nnI,noI,nnaI,nmgI,nalI,nsiI,nkI,ncaI,nfeI,ntiI;
  double ac,an,ao,ana,amg,aal,asi,ak,aca,afe,ati;
  double k = 8.617084e-05;
  double kerg = 1.38054e-16;
  int i;
  extern int Ntau;

  out = fopen("density.out","w");

  Ne = model->Ne[i];

  fprintf(out,"n  N(e-)     N(H)      N(H+)     N(H-)     N(H2)     N(He)     N(He+)    N(He++)\n");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      U = partfn(1.0,T,Ne);

      ph = 4.830e+15*T32*exp(-13.595/kT)/U;
      ph2 = Xi(101.0,T,Ne);
      phm = 2.07e-16*exp(0.7552/kT)/(U*T32);
      phe2 = 4.830e+15*T32*(partfn(2.1,T,Ne)/partfn(2.0,T,Ne))*exp(-24.581/kT);
      phe3 = 4.830e+15*T32*exp(-54.403/kT)/partfn(2.1,T,Ne);
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              model->Ne[i],model->NHI[i],model->NHI[i]*ph/model->Ne[i],model->NHminus[i],model->NH2[i],
              model->NHeI[i],model->NHeI[i]*phe2/model->Ne[i],model->NHeI[i]*phe2*phe3/(model->Ne[i]*model->Ne[i]));
   }

   fprintf(out,"\n");
   fprintf(out,"n  N(C)      N(C+)     N(C++)    N(N)      N(N+)     N(N++)    N(O)      N(O+)     N(O++)\n");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      U = partfn(1.0,T,Ne);
      Ne = model->Ne[i];
      ne2 = Ne*Ne;
      pc2 = 4.830e+15*T32*(partfn(6.1,T,Ne)/partfn(6.0,T,Ne))*exp(-11.260/kT);
      pc3 = 4.830e+15*T32*(partfn(6.2,T,Ne)/partfn(6.1,T,Ne))*exp(-24.383/kT);
      pn2 = 4.830e+15*T32*(partfn(7.1,T,Ne)/partfn(7.0,T,Ne))*exp(-14.534/kT);
      pn3 = 4.830e+15*T32*(partfn(7.2,T,Ne)/partfn(7.1,T,Ne))*exp(-29.601/kT);
      po2 = 4.830e+15*T32*(partfn(8.1,T,Ne)/partfn(8.0,T,Ne))*exp(-13.618/kT);
      po3 = 4.830e+15*T32*(partfn(8.2,T,Ne)/partfn(8.1,T,Ne))*exp(-35.116/kT);
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              model->NCI[i],model->NCI[i]*pc2/Ne,model->NCI[i]*pc2*pc3/ne2,
              model->NNI[i],model->NNI[i]*pn2/Ne,model->NNI[i]*pn2*pn3/ne2,
              model->NOI[i],model->NOI[i]*po2/Ne,model->NOI[i]*po2*po3/ne2);
   }

   fprintf(out,"\n");

   fprintf(out,"n  N(Fe)     N(Fe+)    N(Fe++)   N(Ti)     N(Ti+)    N(Ti++)    N(Ca)     N(Ca+)    N(H(n=3))\n");
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
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              model->NFeI[i],model->NFeI[i]*pfe2/Ne,model->NFeI[i]*pfe2*pfe3/ne2,
              model->NTiI[i],model->NTiI[i]*pti2/Ne,model->NTiI[i]*pti2*pti3/ne2,
              model->NCaI[i],model->NCaII[i],model->N3[i]);
   }

   fprintf(out,"\n");

   fprintf(out,"\n  N(CH)     N(NH)     N(OH)     N(C2)     N(CN)     N(CO)     N(N2)     N(NO)     N(O2)    N(H2O)\n");
   for(i=0;i<Ntau;i++) {
      T = model->T[i];
      T32 = pow(T,1.5);
      kT = k*T;
      U = partfn(1.0,T,Ne);
      ECH = Xi(106.0,T,Ne);
      ENH = Xi(107.0,T,Ne);
      EOH = Xi(108.0,T,Ne);
      EC2 = Xi(606.0,T,Ne);
      ECN = Xi(607.0,T,Ne);
      ECO = Xi(608.0,T,Ne);
      EN2 = Xi(707.0,T,Ne);
      ENO = Xi(708.0,T,Ne);
      EO2 = Xi(808.0,T,Ne);
      C = model->NCI[i];
      N = model->NNI[i];
      O = model->NOI[i];
      H = model->NHI[i];
      fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",i+1,
              C*H*ECH,N*H*ENH,O*H*EOH,C*C*EC2,C*N*ECN,C*O*ECO,N*N*EN2,N*O*ENO,O*O*EO2);
   }

   fprintf(out,"\n");
   fprintf(out,"\n  N(Mg I)    N(Mg II)   N(H I)   N(MgH)\n");
   for(i=0;i<Ntau;i++) {
     fprintf(out,"%2d %9.3e %9.3e %9.3e %9.3e\n",i+1,model->NMgI[i],model->NMgII[i],model->NHI[i],model->NMgH[i]);
   }
   fclose(out);
}
