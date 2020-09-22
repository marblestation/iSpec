#include <math.h>
#include <errno.h>

/* This programs transfers num to the nearest incremental point given
by inc */

void qround(num,inc)
double *num,inc;

{
  double frac,nnum;
  double fnum,cnum;

  nnum = (*num)/inc;
  fnum = floor(nnum);
  frac = nnum - fnum;
  if(frac < 0.5) {
    *num = fnum*inc;
    return;
  } else {
    cnum = ceil(nnum);
    *num = cnum*inc;
    return;
  }
}
