#include <stdio.h>
#include "spectrum.h"

void invelgrad(vgrad)
     char vgrad[80];
{
  int i,ni;
  FILE *invel;
  extern float *velgrad;
  extern int Ntau;

  invel = fopen(vgrad,"r");

  for(i=0;i<Ntau;i++) ni = fscanf(invel,"%f",&velgrad[i]);
}
