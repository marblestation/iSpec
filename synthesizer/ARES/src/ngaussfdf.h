/* 
 * File:   ngaussdfd.h
 * Author: sousasag
 *
 * Created on April 23, 2011, 12:28 PM
 */
/* expfit.c -- compute residual for exponential
   +background model */

struct data {
  size_t n;
  int para;
  double * t;
  double * y;
  double * sigma;
};

int
expb_f (const gsl_vector * x, void *data,
        gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *) data)->sigma;
  int para = ((struct data *)data)->para;

  double ac[para];
  int ijk;
  for (ijk=0;ijk<para;ijk++)
	ac[ijk]=gsl_vector_get (x, ijk);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = sum(i) : Ai * exp(-lambdai * (i-ci)*(i-ci))*/

	double Yi=0.0;
	int j;
        for(j=0;j<para;j+=3)
		{
		 Yi+= ac[j] * exp (-ac[j+1] * (t[i]-ac[j+2])*(t[i]-ac[j+2]) )  ;
		}
	gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);

    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *data,
         gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *) data)->sigma;
  double *t = ((struct data *)data)->t;
  int para = ((struct data *)data)->para;

  double ac[para];
  int ijk;
  for (ijk=0;ijk<para;ijk++)
	ac[ijk]=gsl_vector_get (x, ijk);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /* and the xj are the parameters (Ai,lambdai,ci,...) */

      int j=0;
      double ti = t[i];
      double s = sigma[i];

      for (j=0; j<para; j+=3)
	{
	      double e1 = exp(-ac[j+1]*(ti-ac[j+2])*(ti-ac[j+2]) ) ;

	      gsl_matrix_set (J, i, j, e1/s);
	      gsl_matrix_set (J, i, j+1, - ac[j] * (ti-ac[j+2])*(ti-ac[j+2]) * e1/s);
	      gsl_matrix_set (J, i, j+2,  ac[j] * ac[j+1] * 2. * (ti-ac[j+2]) * e1/s);
	}

    }
  return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *data,
          gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, data, f);
  expb_df (x, data, J);

  return GSL_SUCCESS;
}