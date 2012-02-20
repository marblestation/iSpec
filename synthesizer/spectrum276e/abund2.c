#include <stdio.h>
#include <stdlib.h>
#include "spectrum.h"
#include <errno.h>

double abund(atom,code)
atominfo *atom;
int code;
{
   int k;

   k = 0;
   while(atom[k].code != code) {
     k++;
     if(k >= NATOM) return(-1.0);
   }
   return(atom[k].abund);
}


