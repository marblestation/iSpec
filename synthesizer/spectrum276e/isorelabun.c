#include <stdio.h>
#include <math.h>
#include "spectrum.h"
void getisotope();

void isorelabun(isotope)
isodata *isotope;
{
  extern double ra1H,ra2H,ra12C,ra13C,ra14N,ra15N,ra16O,ra17O,ra18O;
  extern double ra24Mg,ra25Mg,ra26Mg,ra28Si,ra29Si,ra30Si,ra40Ca,ra42Ca;
  extern double ra43Ca,ra44Ca,ra46Ca,ra48Ca,ra46Ti,ra47Ti,ra48Ti,ra49Ti;
  extern double ra50Ti;
  double atmass;

  getisotope(isotope,1.0,1,&atmass,&ra1H);
  getisotope(isotope,1.0,2,&atmass,&ra2H);
  getisotope(isotope,6.0,12,&atmass,&ra12C);
  getisotope(isotope,6.0,13,&atmass,&ra13C);
  getisotope(isotope,7.0,14,&atmass,&ra14N);
  getisotope(isotope,7.0,15,&atmass,&ra15N);
  getisotope(isotope,8.0,16,&atmass,&ra16O);
  getisotope(isotope,8.0,17,&atmass,&ra17O);
  getisotope(isotope,8.0,18,&atmass,&ra18O);
  getisotope(isotope,12.0,24,&atmass,&ra24Mg);
  getisotope(isotope,12.0,25,&atmass,&ra25Mg);
  getisotope(isotope,12.0,26,&atmass,&ra26Mg);
  getisotope(isotope,14.0,28,&atmass,&ra28Si);
  getisotope(isotope,14.0,29,&atmass,&ra29Si);
  getisotope(isotope,14.0,30,&atmass,&ra30Si);
  getisotope(isotope,20.0,40,&atmass,&ra40Ca);
  getisotope(isotope,20.0,42,&atmass,&ra42Ca);
  getisotope(isotope,20.0,43,&atmass,&ra43Ca);
  getisotope(isotope,20.0,44,&atmass,&ra44Ca);
  getisotope(isotope,20.0,46,&atmass,&ra46Ca);
  getisotope(isotope,20.0,48,&atmass,&ra48Ca);
  getisotope(isotope,22.0,46,&atmass,&ra46Ti);
  getisotope(isotope,22.0,47,&atmass,&ra47Ti);
  getisotope(isotope,22.0,48,&atmass,&ra48Ti);
  getisotope(isotope,22.0,49,&atmass,&ra49Ti);
  getisotope(isotope,22.0,50,&atmass,&ra50Ti);

  return;
}  
  

