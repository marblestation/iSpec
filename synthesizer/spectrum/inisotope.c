#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "spectrum.h"
#define TINY 1.0e-09

void inisotope(isofile,isotope)
char isofile[];
isodata *isotope;
{

  extern int NI;
  FILE *iso;
  int i;

  if((iso = fopen(isofile,"r")) == NULL) {
    printf("Cannot open isotope data file %s\n",isofile);
    exit(1);
  }

  i = 0;

  while(fscanf(iso,"%lf %d %lf %lf",&isotope[i].code,&isotope[i].iso,
	       &isotope[i].atomass,&isotope[i].relabund) != EOF) {
    if(isotope[i].relabund <= 0.0) isotope[i].relabund = TINY;
    i++;
  }

  NI = i;

  fclose(iso);

  return;
}
