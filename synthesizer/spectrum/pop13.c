#include <math.h>
#include <errno.h>
#include <stdio.h>
#include "spectrum.h"
double pfunction();
int approx();
void popinit();
double molpartfn();

/* Corrections added Nov 12 & 29, 2004 */
/* H2 added Dec 13, 2007 */
/* AlH and AlO added August 2008 */

void pop(line,N,model,V,POP)
linedata *line;
atmosphere *model;
pfunc *V;
population *POP;
int N;
{
  double Q4,ratio,cneutral,c1ion,c2ion,c3ion,mu,D0,Ua,Ub,Na,Nb,theta,psip;
  int i,charge;
  static int l = 2;
  double k = 8.617084e-05;
  double gffac = 1.0;
  extern int Ntau;
  extern double ra1H,ra2H,ra12C,ra13C,ra14N,ra15N,ra16O,ra17O,ra18O;
  extern double ra24Mg,ra25Mg,ra26Mg,ra28Si,ra29Si,ra30Si,ra40Ca,ra42Ca;
  extern double ra43Ca,ra44Ca,ra46Ca,ra48Ca,ra46Ti,ra47Ti,ra48Ti,ra49Ti;
  extern double ra50Ti;
  extern int flagI;

  if(approx(line[N].code,6.0,0.001) == 1) {
    for(i=0;i<Ntau;i++) line[N].xnum[i] = model->NCI[i]*
			exp(-line[N].El/(k*model->T[i]))/pfunction(V,6.0,i);
    return;
  }

  if(approx(line[N].code,7.0,0.001) == 1) {
    for(i=0;i<Ntau;i++) line[N].xnum[i] = model->NNI[i]*
			exp(-line[N].El/(k*model->T[i]))/pfunction(V,7.0,i);
    return;
  }

  if(approx(line[N].code,8.0,0.001) == 1) {
    for(i=0;i<Ntau;i++) line[N].xnum[i] = model->NOI[i]*
			exp(-line[N].El/(k*model->T[i]))/pfunction(V,8.0,i);
    return;
  }


  if(line[N].code < 100.0) {
    cneutral = floor(line[N].code);
    c1ion = cneutral + .1;
    c2ion = cneutral + .2;
    c3ion = cneutral + .3;
    if(approx(line[N].code,cneutral,0.001) == 1) charge = 0;
    else if(approx(line[N].code,c1ion,0.001) == 1) charge = 1;
    else if(approx(line[N].code,c2ion,0.001) == 1) charge = 2;
    else if(approx(line[N].code,c3ion,0.001) == 1) charge = 3;
    else {
	  printf("Could not determine charge in pop\n");
	  printf("Setting charge = 0\n");
	  charge = 0;
    }
    if(approx(POP[l].atom,line[N].code,0.3) != 1) {
      for(i=2;i<NATOM-NMOL;i++) {
	if(approx(POP[i].atom,line[N].code,0.3) == 1) {
	  l = i;
	  break;
	}
      }
    }

    for(i=0;i<Ntau;i++) {
      Q4 = POP[l].Q4[i];
      if(charge == 0) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 exp(-1.1605e+4*line[N].El/model->T[i])/(Q4*pfunction(V,cneutral,i));
      if(charge == 1) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 POP[l].R1[i]*exp(-1.1605e+4*line[N].El/model->T[i])/
	 (Q4*pfunction(V,c1ion,i));
      if(charge == 2) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 POP[l].R1[i]*POP[l].R2[i]*exp(-1.1605e+4*line[N].El/model->T[i])/
	 (Q4*pfunction(V,c2ion,i));
      if(charge == 3) line[N].xnum[i] = line[N].abund*model->NA[i]*
	 POP[l].R1[i]*POP[l].R2[i]*POP[l].R3[i]*
	 exp(-1.1605e+4*line[N].El/model->T[i])/(Q4*pfunction(V,c3ion,i));
    }
    return;
  } else {
    mu = line[N].chi2;
    D0 = line[N].chi1;
    /* Note: In isotope mode we are tacitly assuming that D0 is the
       same for all isotopic molecules.  This isn't correct, but is
       close enough.  As I find the data, I will add it in  
       (Nov 30, 2004)  */
    gffac = line[N].chi3;

    for(i=0;i<Ntau;i++) {
      if(approx(line[N].code,101.0,0.001) == 1) {
        Ua = pfunction(V,1.0,i);
        Ub = pfunction(V,1.0,i);
        Na = model->NHI[i];
        Nb = model->NHI[i];
      }
      if(approx(line[N].code,106.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,6.0,i);
	Na = model->NHI[i];
	Nb = model->NCI[i];
	if(flagI == 1 && line[N].iso != 0) {
	  if(line[N].iso == 2) {
            Na *= ra2H;
            Nb *= ra12C;
	    mu = 1.724561;
	  }
	  if(line[N].iso == 12) {
            Na *= ra1H;
            Nb *= ra12C;
	    mu = 0.929740;
	  }
	  if(line[N].iso == 13) {
	    Na *= ra1H;
            Nb *= ra13C;
	    mu = 0.935332;
	  }
	}
      }
      if(approx(line[N].code,107.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,7.0,i);
	Na = model->NHI[i];
	Nb = model->NNI[i];
	if(flagI == 1 && line[N].iso != 0) {
          if(line[N].iso == 2) {
            Na *= ra2H;
            Nb *= ra14N;
	    mu = 1.760758;
	  }
	  if(line[N].iso == 14) {
	    Na *= ra1H;
            Nb *= ra14N;
	    mu = 0.940160;
	  }
	  if(line[N].iso == 15) {
	    Na *= ra1H;
            Nb *= ra15N;
	    mu = 0.944375;
	  }
	}
      }
      if(approx(line[N].code,108.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,8.0,i);
	Na = model->NHI[i];
	Nb = model->NOI[i];
	if(flagI == 1 && line[N].iso != 0) {
          if(line[N].iso == 2) {
            Na *= ra2H;
            Nb *= ra16O;
	    mu = 1.788767;
	  }
	  if(line[N].iso == 16) {
            Na *= ra1H;
            Nb *= ra16O;
	    mu = 0.948087;
	  }
	  if(line[N].iso == 17) {
	    Na *= ra1H;
            Nb *= ra17O;
	    mu = 0.951418;
	  }
          if(line[N].iso == 18) {
	    Na *= ra1H;
            Nb *= ra18O;
	    mu = 0.954386;
	  }
	}
      }
      if(approx(line[N].code,112.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,12.0,i);
	Na = model->NHI[i];
	Nb = model->NMgI[i];
        if(flagI == 1 && line[N].iso != 0) {
          if(line[N].iso == 2) {
            Na *= ra2H;
	    Nb *= ra24Mg;
	    mu = 1.857987;
	  }
	  if(line[N].iso == 24) {
	    Na *= ra1H;
            Nb *= ra24Mg;
	    mu = 0.967185;
	  }
	  if(line[N].iso == 25) {
	    Na *= ra1H;
            Nb *= ra25Mg;
	    mu = 0.968750;
	  }
          if(line[N].iso == 26) {
	    Na *= ra1H;
            Nb *= ra26Mg;
	    mu = 0.970193;
	  }
	}
      }
      if(approx(line[N].code,113.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,13.0,i);
	Na = model->NHI[i];
	Nb = model->NAlI[i];
      }
      if(approx(line[N].code,114.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,14.0,i);
	Na = model->NHI[i];
	Nb = model->NSiI[i];
        if(flagI == 1 && line[N].iso != 0) {
          if(line[N].iso == 2) {
            Na *= ra2H;
	    Nb *= ra28Si;
	    mu = 1.878753;
	  }
	  if(line[N].iso == 28) {
	    Na *= ra1H;
            Nb *= ra28Si;
	    mu = 0.972782;
	  }
	  if(line[N].iso == 29) {
	    Na *= ra1H;
            Nb *= ra29Si;
	    mu = 0.973950;
	  }
          if(line[N].iso == 30) {
	    Na *= ra1H;
            Nb *= ra30Si;
	    mu = 0.975041;
	  }
	}
      }
      if(approx(line[N].code,120.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = pfunction(V,20.0,i);
	Na = model->NHI[i];
	Nb = model->NCaI[i];
        if(flagI == 1 && line[N].iso != 0) {
          if(line[N].iso == 2) {
            Na *= ra2H;
	    Nb *= ra40Ca;
	    mu = 1.917370;
	  }
	  if(line[N].iso == 40) {
	    Na *= ra1H;
            Nb *= ra40Ca;
	    mu = 0.983034;
	  }
	  if(line[N].iso == 42) {
	    Na *= ra1H;
            Nb *= ra42Ca;
	    mu = 0.984185;
	  }
          if(line[N].iso == 43) {
	    Na *= ra1H;
            Nb *= ra43Ca;
	    mu = 0.984723;
	  }
          if(line[N].iso == 44) {
	    Na *= ra1H;
            Nb *= ra44Ca;
	    mu = 0.985235;
	  }
          if(line[N].iso == 46) {
	    Na *= ra1H;
            Nb *= ra46Ca;
	    mu = 0.986196;
	  }
          if(line[N].iso == 48) {
	    Na *= ra1H;
            Nb *= ra48Ca;
	    mu = 0.987079;
	  }
	}
      }
      if(approx(line[N].code,606.0,0.001) == 1) {
	Ua = pfunction(V,6.0,i);
	Ub = pfunction(V,6.0,i);
	Na = model->NCI[i];
	Nb = model->NCI[i];
	if(flagI == 1 && line[N].iso != 0) {
	  if(line[N].iso == 12) {
	    Na *= ra12C;
	    Nb *= ra12C;
	    mu = 6.000000;
	  }
	  if(line[N].iso == 13) {
	    Na *= 2.0*ra12C;
	    Nb *= ra13C;
	    mu = 6.240773;
	  }
	  if(line[N].iso == 33) {
	    Na *= ra13C;
	    Nb *= ra13C;
	    mu = 6.501678;
	  }
	}
      }
      if(approx(line[N].code,607.0,0.001) == 1) {
	Ua = pfunction(V,6.0,i);
	Ub = pfunction(V,7.0,i);
	Na = model->NCI[i];
	Nb = model->NNI[i];
	if(flagI == 1 && line[N].iso != 0) {
	  if(line[N].iso == 12) {
            Na *= ra12C;
            Nb *= ra14N;
	    mu = 6.462193;
          }
	  if(line[N].iso == 13) {
            Na *= ra13C;
            Nb *= ra14N;
	    mu = 6.742355;
          }
	  if(line[N].iso == 15) {
            Na *= ra12C;
            Nb *= ra15N;
	    mu = 6.666688;
	  }
	}
      }
      if(approx(line[N].code,608.0,0.001) == 1) {
	Ua = pfunction(V,6.0,i);
	Ub = pfunction(V,8.0,i);
	Na = model->NCI[i];
	Nb = model->NOI[i];
        if(flagI == 1 && line[N].iso != 0) {
	  if(line[N].iso == 12) {
            Na *= ra12C;
            Nb *= ra16O;
	    mu = 6.856209;
          }
	  if(line[N].iso == 13) {
            Na *= ra13C;
            Nb *= ra16O;
	    mu = 7.172413;
          }
	  if(line[N].iso == 17) {
            Na *= ra12C;
            Nb *= ra17O;
	    mu = 7.034334;
	  }
          if(line[N].iso == 18) {
            Na *= ra12C;
            Nb *= ra18O;
	    mu = 7.199866;
	  }
	}
      }
      if(approx(line[N].code,813.0,0.001) == 1) {
	Ua = pfunction(V,8.0,i);
	Ub = pfunction(V,13.0,i);
	Na = model->NOI[i];
	Nb = model->NAlI[i];
      }
      if(approx(line[N].code,814.0,0.001) == 1) {
	Ua = pfunction(V,8.0,i);
	Ub = pfunction(V,14.0,i);
	Na = model->NOI[i];
	Nb = model->NSiI[i];
        if(flagI == 1 && line[N].iso != 0) {
          if(line[N].iso == 16) {
            Na *= ra16O;
            Nb *= ra28Si;
	    mu = 10.176707;
          }
	  if(line[N].iso == 17) {
            Na *= ra17O;
            Nb *= ra28Si;
	    mu = 10.574147;
          }
          if(line[N].iso == 18) {
            Na *= ra18O;
            Nb *= ra28Si;
	    mu = 10.952676;
          }
          if(line[N].iso == 29) {
            Na *= ra16O;
            Nb *= ra29Si;
	    mu = 10.306027;
	  }
          if(line[N].iso == 30) {
            Na *= ra16O;
            Nb *= ra30Si;
	    mu = 10.429446;
	  }
	}
      }
      if(approx(line[N].code,822.0,0.001) == 1) {
	Ua = pfunction(V,8.0,i);
	Ub = pfunction(V,22.0,i);
	Na = model->NOI[i];
	Nb = model->NTiI[i];
      }
      if(approx(line[N].code,840.0,0.001) == 1) {
	Ua = pfunction(V,8.0,i);
	Ub = pfunction(V,40.0,i);
	Na = model->NOI[i];
	Nb = model->NZrI[i];
      }
      if(approx(line[N].code,10108.0,0.001) == 1) {
	Ua = pfunction(V,1.0,i);
	Ub = molpartfn(108.0,model->T[i]);
	Na = model->NHI[i];
	Nb = model->NOH[i];
      }

      theta = 5040.0/model->T[i];
      psip = 1.38054e-16*model->T[i]*pow(10.0,D0*theta-13.670)*
	     pow(theta,2.5)/(pow(mu,1.5)*Ua*Ub);

      line[N].xnum[i] = psip*gffac*Na*Nb*exp(-line[N].El/(k*model->T[i]));
    }
    return;
  }

}

void popinit(POP,atom,model,V,flagw)
population *POP;
atominfo *atom;
atmosphere *model;
pfunc *V;
int flagw;
{
  int i,j;
  double Q1,Q2,Q3,R4,T15;
  extern int Ntau;

  if(flagw == 1) {
    printf("\nCalculating Ionization ratios for all atoms at all levels\n");
    printf("Completed atomic number: ");
  }
  for(i=2;i<NATOM-NMOL;i++) {
    POP[i].atom = atom[i].code;
    for(j=0;j<Ntau;j++) {
      Q1 = Q2 = Q3 = 1.0;
      T15 = pow(model->T[j],1.5);
      if(atom[i].maxcharge == 3)  {
	 R4 = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.4,j)/
	 pfunction(V,POP[i].atom+0.3,j))*T15*
	 exp(-1.16e+4*atom[i].I4/model->T[j]);
	 Q1 += R4;
      }
      if(atom[i].maxcharge >= 2)  {
	 POP[i].R3[j] = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.3,j)/
	 pfunction(V,POP[i].atom+0.2,j))*T15*
	 exp(-1.16e+4*atom[i].I3/model->T[j]);
	 Q2 += POP[i].R3[j]*Q1;
      }
      if(atom[i].maxcharge >= 1)  {
	 POP[i].R2[j] = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.2,j)/
	 pfunction(V,POP[i].atom+0.1,j))*T15*
	 exp(-1.16e+4*atom[i].I2/model->T[j]);
	 Q3 += POP[i].R2[j]*Q2;
	 POP[i].R1[j] = (4.830e+15/model->Ne[j])*(pfunction(V,POP[i].atom+0.1,j)/
	 pfunction(V,POP[i].atom,j))*T15*exp(-1.16e+4*atom[i].I1/model->T[j]);
	 POP[i].Q4[j] = 1.0 + POP[i].R1[j]*Q3;
      }
    }
    if(flagw == 1) printf("%2.0f  ",POP[i].atom);
  }
  if(flagw == 1) printf("\n");
}


