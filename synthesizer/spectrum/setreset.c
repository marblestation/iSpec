#include "spectrum.h"

void setreset(k)
int k;
{
  extern memo reset;

  reset.lyman = reset.balmer = reset.paschen = reset.brackett = reset.pfund = 
     reset.humphreys =  reset.hprofl = reset.helium = reset.strong = 
       reset.interval = k;
}

