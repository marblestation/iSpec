#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "spectrum.h"
#include <string.h>
int approx();
void getisotope();
int maxcharge();

int lline(wave,line,atom,qf,start,end,isotope)
double *wave;
linedata *line;
atominfo *atom;
isodata *isotope;
FILE *qf;
double start,end;
{
  int i,icode,j,l;
  double lambda,code,code1,El,Eu,loggf,df,alp,sig,SA;
  double gammar,gammas,gammaw,gam; 
  char T[5];
  char tmp[10];
  char buffer[120];
  int iso = 0;
  double relabund = 1.0;
  extern int flagI;
  extern int flagSq;
  extern char buf2[120];

  if(fgets(buffer,100,qf) == NULL) return(0);
  if(buffer[0] == '#') return(3);

  strcpy(buf2,buffer);
  lambda = atof(strtok(buffer," "));
  if(lambda < start) return(3);
  if(lambda > end) return(0);
  code = atof(strtok(NULL," "));
  if(code < 0.0) flagSq = 1;
  else flagSq = 0;
  code = fabs(code);  
  i = NATOM+1;
  for(l=0;l<NATOM;l++) {
    if((int)floor(code) == atom[l].code) {
       i = l;
       break;
    }
  } 
  if(i >= NATOM) return(2);
  if(maxcharge(atom,code) == -1) return(2);
  iso = 0;
  if(flagI == 1) iso = atoi(strtok(NULL," "));
  El = atof(strtok(NULL," "));
  Eu = atof(strtok(NULL," "));
  loggf = atof(strtok(NULL," "));
  df = atof(strtok(NULL," "));
  strcpy(T,strtok(NULL," "));
  if(strcmp(T,"AI") == 0) line[0].ai = 1;
  else line[0].ai = 0;
  /* If transition type is AO, read in alp and sig parameters */
  if(strcmp(T,"AO") == 0) {
    SA = atof(strtok(NULL," "));
    sig = floor(SA);
    alp = SA - floor(SA);
  } else alp = sig = 0.0;
  /* If transition type is GA, read in individual gammas, where
     Gamma stark is per electron number, Gamma van der Waals per
     neutral hydrogen number.  Logarithms of these Gammas should
     appear in the linelist */
  /* If the transition type is AI (autoionizing), gammar = Gamma_shore, 
     gammas = Ashore and gammaw = Bshore.  However, if Ashore, the asymmetry
     parameter, is negative, -log(|Ashore|), which will be a positive
     quantity, should be entered into the linelist file.  If Ashore is 
     positive, the entry is log(Ashore), which is negative.  The
     code below translates the various possibilities for GA and AI
     transition types */   
  if(strcmp(T,"GA") == 0 || strcmp(T,"AI") == 0) {
    gammar = pow(10.0,atof(strtok(NULL," ")));
    gam = atof(strtok(NULL," "));
    if(strcmp(T,"AI") == 0) {
      if(gam >= 0.0) gammas = -pow(10.0,-gam);
      else gammas = pow(10.0,gam);
    } else gammas = pow(10.0,gam);
    gammaw = pow(10.0,atof(strtok(NULL," ")));
  } else gammar = gammas = gammaw = 0.0;  

  *wave = lambda;
  line[0].wave = lambda;
  line[0].code = code;
  line[0].iso = iso;
  line[0].El = El*1.23981e-04;
  line[0].Eu = Eu*1.23981e-04;
  line[0].gf = pow(10.0,loggf);
  line[0].fac = df;
  strcpy(line[0].T,T);
  line[0].alp = alp;
  line[0].sig = sig;
  line[0].gammar = gammar;
  line[0].gammas = gammas;
  line[0].gammaw = gammaw;
  if(approx(code,floor(code),0.001) == 1) icode = 0;
  else if(approx(code,floor(code)+0.1,0.001) == 1) icode = 1;
  else if(approx(code,floor(code)+0.2,0.001) == 1) icode = 2;
  else if(approx(code,floor(code)+0.3,0.001) == 1) icode = 3;
  else icode = 0;
  i = 0;
  while((int)code != atom[i].code) i++;
  line[0].atomass = atom[i].amass;
  if(flagI == 1 && code < 100.0) {
    relabund = 1.0;
    getisotope(isotope,code,line[0].iso,&line[0].atomass,
		   &relabund);
    line[0].gf *= relabund;
  }
  line[0].chi1 = atom[i].I1;
  line[0].chi2 = atom[i].I2;
  line[0].chi3 = atom[i].I3;
  line[0].chi4 = atom[i].I4;
  if(icode == 0) line[0].chi = atom[i].I1;
  else if(icode == 1) line[0].chi = atom[i].I2;
  else if(icode == 2) line[0].chi = atom[i].I3;
  else if(icode == 3) line[0].chi = atom[i].I4;
  else line[0].chi = atom[i].I2;
  line[0].abund = atom[i].abund;
  line[0].flag = 0;
  for(j=0;j<NTAU;j++) line[0].xnum[j] = line[0].a[j] =
		      line[0].dopp[j] = 0.0;
  return(1);
}
