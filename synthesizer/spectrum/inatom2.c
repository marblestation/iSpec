#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "spectrum.h"

void inatom(atmdat,atom,MH,ah,ahe)
char atmdat[];
atominfo *atom;
double MH;
double *ah,*ahe;
{
  int i,code;
  FILE *ap;
  double amass,I1,I2,I3,I4,labund;
  int maxcharge;
  char buffer[80];
  int ni;

  if((ap = fopen(atmdat,"r")) == NULL) {
    printf("Cannot open atom data file %s\n",atmdat);
    exit(1);
  }

  if(fgets(buffer,80,ap) == NULL) {
    printf("Error in file access in inatom\n");
    exit(1);
  }
  for(i=0;i<NATOM;i++) {
    ni = fscanf(ap,"%d %lf %lf %lf %lf %lf %lf %d",&code,&labund,&amass,&I1,&I2,&I3,&I4,&maxcharge);
    atom[i].code = code;
    if(i < 2) atom[i].abund = pow(10.0,labund);
    else atom[i].abund = pow(10.0,labund+MH);
    atom[i].amass = amass;
    atom[i].I1 = I1;
    atom[i].I2 = I2;
    atom[i].I3 = I3;
    atom[i].I4 = I4;
    atom[i].maxcharge = maxcharge;
  }
  *ah = atom[0].abund;
  *ahe = atom[1].abund;
  fclose(ap);
}


