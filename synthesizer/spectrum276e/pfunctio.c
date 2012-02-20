#include "spectrum.h"
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
int approx();
int imin();

double pfunction(V,code,N)
pfunc *V;
double code;
int N;
{
  int i;
  int k,nk;

  nk = 0;
  k = (int)imin((int)(3*code),TNATOM-1);
  while(approx(code,V->species[k],0.001) == 0) {
    nk++;
    if(nk >= TNATOM) {
      printf("No information on ion of code %4.1f\n",code);
      exit(1);
    }
    if(code > V->species[k]) k++;
    else k--;
  }
  return((double)V->pf[k][N]);
}


