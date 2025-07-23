*
* Include file for TurboSpectrum
*                                         27_octobre_1996
      INTEGER NDP,NRAYS,lpoint
*
*   NDP should be larger than the number of layers of the model atmosphere
*   NRAYS is the maximum number of rays used to solve the radiative transfer
*         it should be at least equal to NDP + the number of core rays.
*         The number of core rays is the numebr immediately following the 
*         keyword SPHERICAL in the script to run bsyn.f 
*         (i.e. '30' in the case below)
*             'SPHERICAL:'  'T'
*             30
*             300.00
*             15
*             1.30
*   lpoint is the maximum number of wavelength over which the spectrum
*          will be computed
* As there is no dynamic allocation of memory (old fortran!), you may want
* to decrease these numbers, as arrays will be of size 
* real*4: NDP*lpoint
*
      PARAMETER (NDP=150,NRAYS=180,lpoint=500000)
*
*  equation of state parameters   (for tabs and absko)
*
*  mkomp  = max # of absorption coefficient contributors
*  mkompr = max # of temp. indep. absorption coefficient contributors
*
      integer ifadim,kfadim,mkomp,mkompr,nabdim,maxnlatb
      parameter (mkomp=52,mkompr=16,nabdim=15000,
     &  ifadim=(mkomp-mkompr)*ndp,kfadim=ifadim*3+mkompr*ndp,
     &  maxnlatb=100)
*
* number of sets for continuous opacity in common cxlset
* this 400 is a security margin. What is needed is actually
* at maximum lpoint/20 +1 +1.
*
      integer numbset
      parameter (numbset=(lpoint+400)/20)

