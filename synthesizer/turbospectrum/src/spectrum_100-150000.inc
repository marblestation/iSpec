*
* Include file for TurboSpectrum
*                                         27_octobre_1996
      INTEGER NDP,NRAYS,lpoint
*      PARAMETER (NDP=60,NRAYS=80,lpoint=360000)
*      PARAMETER (NDP=150,NRAYS=150,lpoint=150000)
      PARAMETER (NDP=100,NRAYS=120,lpoint=150000)
*      PARAMETER (NDP=750,NRAYS=10,lpoint=150000)
*      PARAMETER (NDP=645,NRAYS=10,lpoint=150000)
*      PARAMETER (NDP=275,NRAYS=275,lpoint=150000)
      
*
*
*  equation of state parameters   (for tabs and absko)
*
*  mkomp  = max # of absorption coefficient contributors
*  mkompr = max # of temp. indep. absorption coefficient contributors
*
      integer ifadim,kfadim,mkomp,mkompr,nabdim
      parameter (mkomp=52,mkompr=16,nabdim=15000,
     &  ifadim=(mkomp-mkompr)*ndp,kfadim=ifadim*3+mkompr*ndp)
* number of sets for continuous opacity in common cxlset
* this 400 is a security margin. What is needed is actually
* at maximum lpoint/20 +1 +1.

      integer numbset
      parameter (numbset=(lpoint+400)/20)

