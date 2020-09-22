#include <math.h>
#include "spectrum.h"
#include <errno.h>
#include <stdio.h>

void capnu(line,N,model)
linedata *line;
int N;
atmosphere *model;
{
  int i;
  double deltnu;
  extern int Ntau;

  for(i=0;i<Ntau;i++) {
    deltnu = line[N].dopp[i]/line[N].wave;
    line[N].capnu[i] = 2.65386e-10*line[N].xnum[i]*line[N].gf*
		model->stim[i]/deltnu;
  }
  line[N].flag = 1;
}

