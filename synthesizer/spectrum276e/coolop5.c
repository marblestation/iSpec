#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "spectrum.h"
double seaton();
double c1op_av(double nu,atmosphere *model,int k);
double al1op();
double mg1op_av(double nu,atmosphere *model,int k);
double si1op();
double ca1op_av(double nu,atmosphere *model,int k);
double fe1op();
double chop();
double ohop();
double mghop();
double coolop();
double partfn(double code,double T,double Ne);

double seaton(nu,nu0,xsect,power,a)
double nu,nu0,xsect,power,a;
{
  double seaton;

  seaton = xsect*(a+(1.0-a)*(nu0/nu))*
	   sqrt(pow(nu0/nu,floor(2.0*power+0.01)));
  return(seaton);
}

double coolop(nu,model,j)
double nu;
atmosphere *model;
int j;
{

  return(c1op_av(nu,model,j)*model->NCI[j] + mg1op_av(nu,model,j)*model->NMgI[j] +
	 al1op(nu,model,j)*model->NAlI[j] + si1op(nu,model,j)*model->NSiI[j] +
	 ca1op_av(nu,model,j)*model->NCaI[j] +
	 fe1op(nu,model,j)*model->NFeI[j] + chop(nu,model,j)*model->NCH[j] +
         ohop(nu,model,j)*model->NOH[j] + 
         mghop(nu,model,j)*model->NMgH[j]);
}





