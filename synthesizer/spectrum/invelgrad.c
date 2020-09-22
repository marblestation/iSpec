#include <stdio.h>
#include "spectrum.h"

void invelgrad(vgrad)
char vgrad[80];
{
    int i;
    FILE *invel;
    extern float *velgrad;
    extern int Ntau;

    invel = fopen(vgrad, "r");

    for (i = 0; i < Ntau; i++)
        fscanf(invel, "%f", &velgrad[i]);
}
