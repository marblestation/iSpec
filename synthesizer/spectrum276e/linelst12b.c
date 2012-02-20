#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "spectrum.h"
double abund();
double dmax();
double dmin();
void qround();
int maxcharge();
int approx();
void getisotope();
void pop();
double initrad(linelist *list, int n, atominfo *atom, double loggf, double El,
double Eu, double radmax, double teff, double logg, double code, double lambda);
double polishrad(linelist *list, int n, atominfo *atom, double rad, 
		 isodata *isotope, atmosphere *model, double Flux, pfunc *V,
                 population *POP);
void broad(atmosphere *model,linedata *line,int N,double sig,
           double alp,double fac);
void capnu();
void eqtaukap();
double depth();

void linelst(wave,list,nlist,atom,teff,logg,qf,reset,isotope,model,Flux,V,POP,
dwave)
double wave,teff,logg;
linelist *list;
atominfo *atom;
isodata *isotope;
atmosphere *model;
pfunc *V;
population *POP;
int *nlist,reset;
double Flux;
double dwave;
FILE *qf;
{
  int i,l,n;
  extern double inc;
  extern int flagr;
  extern int flagI;
  double lambda,code,El,Eu,loggf,ab,dampfac,SA,alp,sig;
  double gammar,gammas,gammaw,gam; 
  /* Gamma stark per electron number */
  /* Gamma van der waals per neutral hydrogen number */
  char tr[5];
  char tmp[10],buffer[100],src[20];
  static double lastwave = 0.0;
  double radmax = 20.0;
  double wave1,wave2,gfA,rad,rad2,lograd,neff,chi;
  double charg;
  double err;
  double hfsfac = 1.0;
  static int flag = 0;
  static int qfeof = 0;
  double fac = 1.0;
  int iso = 0;
  int flagai = 0;

  if(flag == 0 || reset == 1) lastwave = wave - radmax;

  if(wave < 3646.0) radmax = 25.0;
  else if(wave >= 8000.0) radmax = 25.0;
  else radmax = 20.0;
  wave1 = wave + radmax;
  wave2 = wave - radmax;

  if(qfeof == 1) return;
  if(lastwave > wave1) return;
  while(lastwave <= wave1) {
    if(fgets(buffer,100,qf) == NULL) {
      qfeof = 1;
      break;
    }
    /* printf("%s",buffer); */
    if(buffer[0] == '#') continue;
    lambda = atof(strtok(buffer," "));
    if(lambda < wave2) {
      lastwave = lambda;
      continue;
    }
    hfsfac = 1.0;
    code = atof(strtok(NULL," "));
    /* a provision for hyperfine structure components */
    if(code < 0.0) {
      hfsfac = 3.0;
      code = fabs(code);
    }
    iso = 0;
    if(flagI == 1) iso = atoi(strtok(NULL," "));
    El = atof(strtok(NULL," "));
    Eu = atof(strtok(NULL," "));
    loggf = atof(strtok(NULL," "));
    dampfac = atof(strtok(NULL," "));
    strcpy(tr,strtok(NULL," "));
    if(strcmp(tr,"AI") == 0) flagai = 1;
    else flagai = 0;
    /* If transition type is AO, read in alp and sig parameters */
    if(strcmp(tr,"AO") == 0) {
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
    if(strcmp(tr,"GA") == 0 || strcmp(tr,"AI") == 0) {
      gammar = pow(10.0,atof(strtok(NULL," ")));
      gam = atof(strtok(NULL," "));
      if(strcmp(tr,"AI") == 0) {
	if(gam >= 0.0) gammas = -pow(10.0,-gam);
	else gammas = pow(10.0,gam);
      } else gammas = pow(10.0,gam);
      gammaw = pow(10.0,atof(strtok(NULL," ")));
    } else gammar = gammas = gammaw = 0.0;
    err = atoi(strtok(NULL," "));
    if(err == 99) printf("You may need to use the i switch\n");
    if(lambda < lastwave) continue;
    lastwave = lambda;
    if(flagr == 1) qround(&lambda,inc);
/* Skip over molecules if effective temperature is too high */
    if(teff >= 8500.0 && code > 101.0) continue;
/* If atom is not supported by SPECTRUM, skip over the line */
    if(abund(atom,(int)floor(code)) < 0.0) continue;
/* If an unsupported ion of a supported atom is in the list, skip */
    if(maxcharge(atom,code) == -1) continue;
    list[*nlist].wave = lambda;
    list[*nlist].code = code;
    list[*nlist].iso = iso;
    list[*nlist].El = El;
    El = El*1.23981e-04;
    list[*nlist].Eu = Eu;
    list[*nlist].loggf = loggf;
    list[*nlist].fac = dampfac;
    strcpy(list[*nlist].T,tr);
    list[*nlist].alp = alp;
    list[*nlist].sig = sig;
    list[*nlist].gammar = gammar;
    list[*nlist].gammas = gammas;
    list[*nlist].gammaw = gammaw;
    if(flagai == 1) list[*nlist].ai = 1;
    else list[*nlist].ai = 0;
    rad = 0.0;
    if(flagai == 1) rad = radmax;
    else {
      rad = initrad(list,*nlist,atom,loggf,El,Eu,radmax,teff,logg,code,
                    lambda);
      rad = polishrad(list,*nlist,atom,rad,isotope,model,Flux,V,POP);

      /* Increase computation radius for hyperfine components */
      rad *= hfsfac;
      /* Limit computation radius to radmax */
      rad = (double)dmin(rad,radmax);     
    }
    if(rad == 0.0) continue;
    /* rad = (double)dmax(rad,1.1*dwave); */
    if((lambda+rad) < wave) continue;
    list[*nlist].wavel = lambda - rad;
    list[*nlist].waveh = lambda + rad;
    list[*nlist].flag = 0;

    (*nlist)++;
    if(*nlist >= NLIST-2) break;
  }

  if(flag == 0) {
    flag = 1;
    return;
  }

  if(*nlist <= 1) return;
  while(list[0].wave <= wave) {
    for(n=0;n<(*nlist)-1;n++) {
       list[n].wave = list[n+1].wave;
       list[n].code = list[n+1].code;
       list[n].iso = list[n+1].iso;
       list[n].loggf = list[n+1].loggf;
       list[n].El = list[n+1].El;
       list[n].Eu = list[n+1].Eu;
       list[n].wavel = list[n+1].wavel;
       list[n].waveh = list[n+1].waveh;
       list[n].fac = list[n+1].fac;
       list[n].flag = list[n+1].flag;
       list[n].alp = list[n+1].alp;
       list[n].sig = list[n+1].sig;
       list[n].gammar = list[n+1].gammar;
       list[n].gammas = list[n+1].gammas;
       list[n].gammaw = list[n+1].gammaw;
       list[n].ai = list[n+1].ai;
       strcpy(list[n].T,list[n+1].T);
    }
    (*nlist)--;
    if(*nlist <= 1) break;
  }

}

double initrad(linelist *list, int n, atominfo *atom, double loggf, double El,
double Eu, double radmax, double teff, double logg, double code, double lambda)
{
  double ab,rad,gfA,charge,lograd,fac,neff,charg,chi;
  int i,l;

  ab = abund(atom,(int)floor(list[n].code));
  gfA = loggf + log10(ab);
  charg = list[n].code - floor(list[n].code);
  if(list[n].code - floor(list[n].code) >= 0.01) lograd = 3.96 +
      (0.75861 + 0.0258773*gfA)*gfA -0.242223*El;
  else lograd = 4.46 + (0.741125 + 0.0215996*gfA)*gfA +
      ((0.0639498*El - 0.311767)*El + 0.0891479)*El;
  if(lograd > log10(radmax)) rad = radmax;
  else rad = pow(10.0,lograd);
  /* Make a rough temperature adjustment in the radius */
  if(teff > 5700.0) fac = 1.0;
  else fac = 1.0+(5700.0-teff)/400.0;
  rad *= fac;
  if(logg > 4.0) fac = 1.0+ logg-4.0;
  rad *= fac;
  /* Make adjustment for high excitation lines like Mg II 4481 */
  if(code < 100.0 && Eu > 7.0e+04) rad *= 6.0;
  /* Make adjustment for neutral lines near ionization limit */
  if(code - floor(code) < 0.05 && code < 100.0) {
    for(l=0;l<NATOM;l++) {
       if((int)floor(code) == atom[l].code) {
          i = l;
          break;
       }
    } 
    chi = atom[i].I1;
    if(chi > Eu*1.23981e-04) {
       neff = sqrt(13.585/(chi-Eu*1.23981e-04));
       if(neff > 10.0) rad *= pow(neff/10.0,2.0);
    }
  } 
  /* Make adjustment for CNO atomic lines */
  if(code >= 6.0 && code < 9.0) rad *= 4.0;
  rad = (double)dmax(rad,0.20);
  /* if(code > 100.0) rad = (double)dmin(rad,0.5); */
  /* Account for increasing doppler widths */
  if(lambda > 5000.0) rad *= lambda/5000.0;
  /* Make adjustment for minimum in H- opacity */
  if(lambda > 8000.0 && lambda <= 13000.0) rad *= 2.7;
  if(lambda > 13000.0 && lambda < 19000.0) rad *= 3.8;
  if(lambda >= 19000.0 && lambda < 50000.0) rad *= 2.7;
  /* Make adjustment for non-neutral species */
  if(code < 100 && charg > 0.0) 
     rad *= 3.5;
  if(code < 100 && charg > 0.1)
     rad *= 1.5;
  /* Increase computation radius for low-excitation lines */
  if(El < 0.124) rad *= 3.0;

  return(rad);
}

double polishrad(linelist *list, int n, atominfo *atom, double rad, 
isodata *isotope, atmosphere *model, double Flux, pfunc *V, population *POP)
{
  extern linedata *oneline;
  extern int flagI;
  extern int Ntau;
  double code,relabund,wave;
  int i,j,k,l,icode;
  double Depth,Depth1,Depth2,radnew,A;

  oneline[0].wave = list[n].wave;
  code = oneline[0].code = list[n].code;
  oneline[0].iso = list[n].iso;
  oneline[0].El = list[n].El*1.23981e-04;
  oneline[0].Eu = list[n].Eu*1.23981e-04;
  oneline[0].gf = pow(10.0,list[n].loggf);
  wave = list[n].wave;
  oneline[0].wavel = list[n].wave - rad;
  oneline[0].waveh = list[n].wave + rad;
  oneline[0].fac = list[n].fac;
  oneline[0].alp = list[n].alp;
  oneline[0].sig = list[n].sig;
  oneline[0].gammar = list[n].gammar;
  oneline[0].gammas = list[n].gammas;
  oneline[0].gammaw = list[n].gammaw;
  oneline[0].ai = list[n].ai;
  strcpy(oneline[0].T,list[n].T);
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
  }
  oneline[0].atomass = atom[i].amass;
  if(flagI == 1 && code < 100.0) {
     relabund = 1.0;
     getisotope(isotope,code,oneline[0].iso,&oneline[0].atomass,
		&relabund);
     oneline[0].gf *= relabund;
   }
   oneline[0].chi1 = atom[i].I1;
   oneline[0].chi2 = atom[i].I2;
   oneline[0].chi3 = atom[i].I3;
   oneline[0].chi4 = atom[i].I4;
   if(icode == 0) oneline[0].chi = atom[i].I1;
   else if(icode == 1) oneline[0].chi = atom[i].I2;
   else if(icode == 2) oneline[0].chi = atom[i].I3;
   else if(icode == 3) oneline[0].chi = atom[i].I4;
   else oneline[0].chi = atom[i].I2;
   oneline[0].abund = atom[i].abund;
   oneline[0].flag = 0;
   for(j=0;j<Ntau;j++) oneline[0].xnum[j] = oneline[0].a[j] =
		       oneline[0].dopp[j] = 0.0;

   pop(oneline,0,model,V,POP);
   broad(model,oneline,0,oneline[0].sig,oneline[0].alp,oneline[0].fac);
   capnu(oneline,0,model);
   eqtaukap(wave,model,oneline);
   Depth = depth(model,wave,Flux);
   if(Depth < 0.0001) {
     rad = 0.0;
     return(rad);
   }

   /* Polish rad by finding where in the wing depth ~ 0.0001 */
   eqtaukap(wave+rad,model,oneline);
   Depth = depth(model,wave+rad,Flux);
   k = 0;
   while(fabs(10000.0*(Depth-0.0001)) > 0.20) {
     if(k > 5) break;
     A = rad*rad*Depth;
     radnew = sqrt(A/0.0001);
     rad = (rad + radnew)/2.0; /* Some lines will not converge without this */
     eqtaukap(wave+rad,model,oneline);
     Depth = depth(model,wave+rad,Flux);
     k++;
   }

   return(rad);
  
}
