#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "spectrum.h"

void getisotope(isotope,code,iso,atmass,relabund)
     isodata *isotope;
     double code;
     int iso;
     double *atmass,*relabund;
{
  extern int NI;
  int flagiso = 0;

  int i;

  if(NI == 0) return;
  if(iso == 0) return;

  for(i=0;i<NI;i++) {
    if(isotope[i].code == floor(code) && isotope[i].iso == iso) {
      *atmass = isotope[i].atomass;
      *relabund = isotope[i].relabund;
      flagiso = 1;
      break;
    }
  }

  if(flagiso == 0) {
    printf("\nThis isotope not recognized by SPECTRUM %5.1f %d\n",code,iso);
    exit(1);
  }

  return;
} 
