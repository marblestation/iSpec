#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "spectrum.h"
int approx();
void getisotope();

void inlin(wave,line,nline,list,nlist,atom,isotope)
double wave;
linedata *line;
linelist *list;
atominfo *atom;
isodata *isotope;
int *nline,nlist;
{
  int i,icode,j,k,l,n;
  double code,code1,El,Eu,loggf,relabund = 1.0;
  static double lastwave = 0.0;
  extern int Ntau;
  extern int flagI;

  for(k=0;k<nlist;k++) {
    if(list[k].flag == 1) continue;
    if(*nline >= NLINE-1) continue;
    if((wave >= list[k].wavel && wave <= list[k].waveh) || k == 0) {
      line[*nline].wave = list[k].wave;
      code = line[*nline].code = list[k].code;
      line[*nline].iso = list[k].iso;
      line[*nline].El = list[k].El*1.23981e-04;
      line[*nline].Eu = list[k].Eu*1.23981e-04;
      line[*nline].gf = pow(10.0,list[k].loggf);
      line[*nline].wavel = list[k].wavel;
      line[*nline].waveh = list[k].waveh;
      line[*nline].fac = list[k].fac;
      line[*nline].alp = list[k].alp;
      line[*nline].sig = list[k].sig;
      line[*nline].gammar = list[k].gammar;
      line[*nline].gammas = list[k].gammas;
      line[*nline].gammaw = list[k].gammaw;
      line[*nline].ai = list[k].ai;
      strcpy(line[*nline].T,list[k].T);
      if(approx(code,floor(code),0.001) == 1) icode = 0;
      else if(approx(code,floor(code)+0.1,0.001) == 1) icode = 1;
      else if(approx(code,floor(code)+0.2,0.001) == 1) icode = 2;
      else if(approx(code,floor(code)+0.3,0.001) == 1) icode = 3;
      else icode = 0;
      i = NATOM+1;
      for(l=0;l<NATOM;l++) {
         if((int)floor(code) == atom[l].code) {
           i = l;
           break;
         }
      } 
      if(i >= NATOM) {
        printf("Warning, SPECTRUM does not support element number %d\n", 
               (int)code);
        continue;
      }
      line[*nline].atomass = atom[i].amass;
      if(flagI == 1 && code < 100.0) {
	relabund = 1.0;
	getisotope(isotope,code,line[*nline].iso,&line[*nline].atomass,
		   &relabund);
	line[*nline].gf *= relabund;
      }
      line[*nline].chi1 = atom[i].I1;
      line[*nline].chi2 = atom[i].I2;
      line[*nline].chi3 = atom[i].I3;
      line[*nline].chi4 = atom[i].I4;
      if(icode == 0) line[*nline].chi = atom[i].I1;
      else if(icode == 1) line[*nline].chi = atom[i].I2;
      else if(icode == 2) line[*nline].chi = atom[i].I3;
      else if(icode == 3) line[*nline].chi = atom[i].I4;
      else line[*nline].chi = atom[i].I2;
      line[*nline].abund = atom[i].abund;
      line[*nline].flag = 0;
      for(j=0;j<Ntau;j++) line[*nline].xnum[j] = line[*nline].a[j] =
			  line[*nline].dopp[j] = 0.0;
      list[k].flag = 1;

      lastwave = wave;
      (*nline)++;
      if(*nline >= NLINE-1) break;
    }
  }

  n = 0;
  while(n < (*nline)-1) {
     if(line[n].waveh < wave) {
       for(i=n;i<(*nline)-1;i++) {
	 line[i].wave = line[i+1].wave;
	 line[i].code = line[i+1].code;
	 line[i].iso = line[i+1].iso;
	 line[i].atomass = line[i+1].atomass;
	 line[i].abund = line[i+1].abund;
	 line[i].chi1 = line[i+1].chi1;
	 line[i].chi2 = line[i+1].chi2;
	 line[i].chi3 = line[i+1].chi3;
	 line[i].chi4 = line[i+1].chi4;
	 line[i].chi = line[i+1].chi;
	 line[i].gf = line[i+1].gf;
	 line[i].El = line[i+1].El;
	 line[i].Eu = line[i+1].Eu;
	 line[i].wavel = line[i+1].wavel;
	 line[i].waveh = line[i+1].waveh;
	 line[i].fac = line[i+1].fac;
	 line[i].flag = line[i+1].flag;
	 line[i].alp = line[i+1].alp;
	 line[i].sig = line[i+1].sig;
         line[i].gammar = line[i+1].gammar;
         line[i].gammas = line[i+1].gammas;
         line[i].gammaw = line[i+1].gammaw;
	 line[i].ai = line[i+1].ai;
         strcpy(line[i].T,line[i+1].T);
	 for(j=0;j<Ntau;j++) {
	   line[i].xnum[j] = line[i+1].xnum[j];
	   line[i].a[j] = line[i+1].a[j];
	   line[i].dopp[j] = line[i+1].dopp[j];
	   line[i].capnu[j] = line[i+1].capnu[j];
	 }
      }
      (*nline)--;
    } else n++;
  }

}
