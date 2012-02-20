#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "spectrum.h"
int approx();

int maxcharge(atom,code)
atominfo *atom;
double code;
{
  int i,k;
  int n;
  int Code;
  int charge;
  double fcharge;

  charge = 0;
  Code = (int)floor(code);
  fcharge = 10.0*(code - floor(code));
  n = (int)fcharge + 2;
  for(i=0;i<n;i++) {
     if(approx(fcharge,(double)i,0.05) == 1) charge = i;
  }
  k = 0;
  while(atom[k].code != Code) {
    k++;
    if(k >= NATOM) return(-1);
  }

  if(charge > atom[k].maxcharge) return(-1);
  else return(0);
}
