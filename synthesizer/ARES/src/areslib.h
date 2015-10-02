/* 
 * File:   areslib.h
 * Author: sousasag
 *
 * Created on April 25, 2011, 2:58 PM
 */

#ifndef _ARESLIB_H
#define	_ARESLIB_H

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_interp.h>
#include "ngaussfdf.h"
#include "filesio.h"
#include "aresplot.h"
#define max(a,b) (((a)>(b))?(a):(b))
#ifdef	__cplusplus
extern "C" {
#endif
    
void clean_zero_gaps(double* flux, long np);
void arraysubcp(double *, double *, long, long);
void poly_fitn(double *, double *, double *, long, long, double *);
int continuum_det5 (double *, double *, double *, long, double *, double, int);
void deriv(double *, double *, double *, long);
void smooth(double *, long, int, double *);
void zeroscenterfind(double *, double *, double *, double *, long, long *, long *);
double maxele_vec(double *, long);
void fitngauss(double *, double *, double *, long,  double *, double *, int, int *);

long find_pixel_line(double * xpixels, double linha);

void getMedida(double * xpixels, double * pixels, float linha, double space, double tree, int* plots_flag2, double smoothder, double distlinha, FILE * pFile3, int ilinha, double *aponta, double lambdai, double lambdaf);



void clean_zero_gaps(double* flux, long np){
    int i;
    for (i=0; i<np; i++)
        if (flux[i] <= 0 || flux[i] != flux[i]) flux[i]=1.;
}


long find_pixel_line(double * xpixels, double linha){
    long nctest;
    double restest[2];
    double cdelta1=xpixels[1]-xpixels[0];
    double crval1=xpixels[0];
    restest[0]=1./cdelta1;
    restest[1]=-crval1/cdelta1;
    nctest=(long) (restest[0]*linha+restest[1]);
    nctest++;
// implementar verificação de proximidade da linha. Para o caso de cdeltas nao equidistantes. Neste caso podemos implementar um if para ver se está suficientemente perto.
// Se nao estiver perto um while até se encontrar perto depois de verificar se tem de somar ou subtrair (cuidado com os limites)
    return nctest;
}

void getMedida(double * xpixels, double * pixels, float linha, double space, double rejt, int* plots_flag2, double smoothder, double distlinha, FILE * pFile3 , int ilinha, double *aponta, double lambdai, double lambdaf){

    //definicao dos pontos do intervalo local para normalizar o espectro a volta da linha
            int i, status2;
            int plots_flag=*plots_flag2;
            long nctest=find_pixel_line(xpixels, linha);
            long nx1test=find_pixel_line(xpixels, linha-space);
            long nx2test=find_pixel_line(xpixels, linha+space);

            char strLinhaInicial[100];
            strcpy(strLinhaInicial,"  ");
            sprintf(strLinhaInicial,"\n\nline nº %i searching for line %.2f in the interval [%.2f,%.2f]. Using rejt %f \n\n",ilinha+1, linha, lambdai, lambdaf, rejt);
            printf("%s",strLinhaInicial);

            double xltest[nx2test-nx1test], atest[nx2test-nx1test];
            arraysubcp(xltest, xpixels,nx1test,nx2test );
            arraysubcp(atest ,  pixels,nx1test,nx2test );

    // encontrar o continuum
            // the -1 in the (nx1test-1) and the +1 in the nx is to keep the same points as in ARES v1.
            long nx=nx2test-nx1test+1;

            double x[nx],y[nx], ynorm[nx];
            arraysubcp(x, xpixels,nx1test-1,nx2test );
            arraysubcp(y,  pixels,nx1test-1,nx2test );
            double res[4];
            int testflag = continuum_det5(x,y,ynorm,nx,res,rejt,plots_flag);
            if (testflag == -1) {
            	printf("Problem with the normalization\n Ignoring this line\n");
            	aponta[ilinha*9+4]=-1;
                //Escrever no ficheiro de Log:
                pFile3 = fopen ("logARES.txt","a");
                fprintf(pFile3,"%s%s",strLinhaInicial,"Problem with the normalization\n Ignoring this line\n");
                fclose (pFile3);
                //Nothing more to do here
            	return;
            }

            for (i=0; i<nx; i++)
                    y[i]=y[i]/(res[0]+res[1]*x[i]+res[2]*x[i]*x[i]+res[3]*x[i]*x[i]*x[i]);

            //encontro dos pontos extremos(xind1,xind2) para o calculo das derivadas...  Encontrar os extremos para o fit. 
            //Nao se usa o space todo para o fit. O space todo apenas e usado para a determinacao local do continuum

            int xind1=0,xind2=nx-1,hjk;
            float klo=0.1;
            for (hjk=0; hjk < nx; hjk++){
                if ( (y[hjk] > rejt) && (x[hjk]-(linha-klo) > x[xind1] - (linha-klo)) && (x[hjk] - (linha-klo) < 0) )
                    xind1=hjk;
                if ( (y[hjk] > rejt) && (x[hjk]-(linha+klo) < x[xind2] - (linha+klo)) && (x[hjk] - (linha+klo) > 0) )
                    xind2=hjk;
            }

            int nlin=xind2-xind1;
            double xlin[nlin], iylin[nlin], ylin[nlin], dylin[nlin], ddylin[nlin];
            double ylincaga[nx], dylincaga[nx], ddylincaga[nx], tmpcaga[nx];
            arraysubcp(xlin, x,xind1,xind2 );
            arraysubcp(iylin, y,xind1,xind2 );


            // Calculo das derivadas e respectivo smooth para a zona para o fit...
            deriv(x,y,tmpcaga,nx);
            smooth(tmpcaga, nx, (int)smoothder, ylincaga);
            arraysubcp(ylin, ylincaga,xind1,xind2 );

            deriv(x,ylincaga,tmpcaga,nx);
            smooth(tmpcaga, nx, (int)smoothder, dylincaga);
            arraysubcp(dylin, dylincaga,xind1,xind2 );

            deriv(x,dylincaga,tmpcaga,nx);
            smooth(tmpcaga, nx, (int)smoothder, ddylincaga);
            arraysubcp(ddylin, ddylincaga,xind1,xind2 );


//	procura das riscas que ha a volta da risca que queremos

            double cont[nlin], zeros[nlin];
            long ncont=nlin, nzeros=nlin, ncenter=nlin, center[nlin];

            zeroscenterfind(ylin, iylin, dylin, ddylin, nlin, center, &ncenter);

//	calculo, interpolacao da posicao das riscas no espectro

            if (center[0] != -1 & ncenter != 0) {
                double xlinhas[ncenter], ylinhas[ncenter];
                int i1, i2;
                for (i=0; i<ncenter; i++) {
                    i1=center[i];
                    xlinhas[i]= ( -ddylin[center[i]-1] + ( ddylin[center[i]] - ddylin[center[i]-1] )/( xlin[center[i]]-xlin[center[i]-1] ) * xlin[center[i]] )  / ( ( ddylin[center[i]] - ddylin[center[i]-1] )/(xlin[center[i]]-xlin[center[i]-1]) );
                    ylinhas[i]= ( iylin[center[i]] - iylin[center[i]-1] )/( xlin[center[i]] -xlin[center[i]-1] ) * xlinhas[i] + iylin[center[i]-1] - ( iylin[center[i]]- iylin[center[i]-1] )/( xlin[center[i]] -xlin[center[i]-1]) * xlin[center[i]] ;
                }

                char strLinhaFound[ncenter*8+30];
                strcpy(strLinhaFound,"\n LINES FOUND TO FIT \n");
                
		for (i=0; i<ncenter; i++) {
                    char strtmp[9];
                    sprintf(strtmp,"%.2f ", xlinhas[i]);
                    strcat(strLinhaFound,strtmp);
                }
                strcat(strLinhaFound,"\n");
                printf("%s",strLinhaFound);
                
                //RESAMPLING, Eliminacao das riscas que estao muito juntas...
		double xvec2[ncenter], yvec2[ncenter];

                int nvec2,j;
		xvec2[0]=xlinhas[0];
		yvec2[0]=ylinhas[0];
		j=0;
		for(i=1;i<ncenter;i++) {
                    if (fabs(xvec2[j]-xlinhas[i]) < distlinha ) {
                        xvec2[j]=(xvec2[j]+xlinhas[i])/2.;
                        yvec2[j]=(yvec2[j]+ylinhas[i])/2.;
                    } else {
                        j++;
                        xvec2[j]=xlinhas[i];
                        yvec2[j]=ylinhas[i];
                    }
		}
		nvec2=j+1;
                
                
                char strLinhaResample[nvec2*8+30];
                strcpy(strLinhaResample,"\n RESAMPLING \n");

		for (i=0; i<nvec2; i++)	{
                    char strtmp[9];
                    sprintf(strtmp,"%.2f ", xvec2[i]);
                    strcat(strLinhaResample,strtmp);
                }
                strcat(strLinhaResample,"\n");
                printf("%s",strLinhaResample);

		ncenter=nvec2;
		int para=3*ncenter;
		int npara=0;
		double acoef[para], acoef_er[para];     //initial guesses
		for (i=0;i<ncenter;i++) {
                    acoef[3*npara]=yvec2[i]-1.;
                    acoef[3*npara+1]=400.;
                    acoef[3*npara+2]=xvec2[i];
                    npara++;
		}

		double xfit[nlin], yfit[nlin], sigma[nlin];
		for (i=0;i<nlin;i++) {
                    xfit[i]=xlin[i];
                    yfit[i]=iylin[i]-1.0;
//                    sigma[i]=0.1;   //NEED to DEFINE a better sigma (dependent on the S/N)
                    sigma[i]=1.-rejt;   //NEED to DEFINE a better sigma (dependent on the S/N)
		}
                
                char strLinhaGuess[para*65+30];
                strcpy(strLinhaGuess,"\n GUESS COEFS :\n");
		for (i=0;i<para;i+=3){
                    char strtmp[65];
                    sprintf(strtmp,"acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);
                    strcat(strLinhaGuess,strtmp);
                }
                printf("%s",strLinhaGuess);
                
		fitngauss(xfit,yfit,sigma,nlin,acoef,acoef_er,para,&status2);

		char strLinhaFitted[para*200+30];
                strcpy(strLinhaFitted,"\n FITTED COEFS :\n");
		for (i=0;i<para;i+=3){
                    char strtmp[2000];
                    sprintf(strtmp,"::acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);
                    strcat(strLinhaFitted,strtmp);
                    sprintf(strtmp,"+-ac_er[%2i]:  %.5f ac_er[%2i]:  %9.5f ac_er[%2i]:  %9.6f \n", i, acoef_er[i], i+1, acoef_er[i+1], i+2, acoef_er[i+2]);
                    strcat(strLinhaFitted,strtmp);
                }
                printf("%s",strLinhaFitted);
                
		double yfit2[nx];
		for (i=0;i<nx;i++) {
                    yfit2[i]=1.0;
                    for (j=0;j<ncenter;j++)
                        yfit2[i]+=acoef[j*3]* exp (- acoef[j*3+1] * (x[i]-acoef[j*3+2]) * (x[i]-acoef[j*3+2]) );
		}

		double medida=0, medida_er_square=0, medida_er=0;
		int nmed=0, hj=0, hjl=0;

		for (hj=0; hj<ncenter;hj++) {
                    if ( fabs(linha-acoef[3*hj+2]) < distlinha ) {
                        medida+=acoef[3*hj]*sqrt(3.1415927/acoef[3*hj+1]);
                        medida_er_square+=medida*medida * ( acoef_er[3*hj]*acoef_er[3*hj]/acoef[3*hj]/acoef[3*hj] + (0.5*0.5*acoef_er[3*hj+1]*acoef_er[3*hj+1]/acoef[3*hj+1]/acoef[3*hj+1]));
                        nmed++;
                        hjl=hj;
                    }
                }

		medida=medida*(-1000.);
                medida_er=sqrt(medida_er_square)*1000.;
                char strLinhaResult[300];
                char strtmp[100];
                sprintf(strtmp,"\n---------------------------\nline result: %.5f \n", linha);
                strcpy(strLinhaResult,strtmp);
                sprintf(strtmp,"ew (mA)      :  %.5f \n", medida);
                strcat(strLinhaResult,strtmp);
                sprintf(strtmp,"ew error (mA):  %.5f \n", medida_er);
                strcat(strLinhaResult,strtmp);
                sprintf(strtmp,"nfit : %ld \n", ncenter);
                strcat(strLinhaResult,strtmp);


		if (nmed == 1) {
                    sprintf(strtmp,"line depth : %.5f \n", -acoef[3*hjl]);
                    strcat(strLinhaResult,strtmp);
                    // FWHM para a gaussiana defenida: F(X)=Aexp(-Lambda(x-c)^2) => FWHM=2*sqrt(ln(2)/lambda)
                    sprintf(strtmp,"FWHM : %.5f \n-------------------------\n", 2.*sqrt(log(2)/acoef[3*hjl+1]));
                    strcat(strLinhaResult,strtmp);
                    // FWHM para a gaussiana defenida: F(X)=Aexp(-Lambda(x-c)^2) => FWHM=2*sqrt(ln(2)/lambda)
		}

		sprintf(strtmp,"int 2 status: %i", status2);
                strcat(strLinhaResult,strtmp);
                printf("%s",strLinhaResult);

		if (plots_flag == 1) {
			double xcvec[ncenter], ycvec[ncenter];
			for (i=0;i<ncenter;i++) {
				xcvec[i]=acoef[i*3+2];
				ycvec[i]=acoef[i*3]+1.;
			}
                        plotxyover2(x,y,nx,x,yfit2,nx,linha-space,linha+space);
                        int pausav;
			printf ("\n\nTo Close the plots, click on it.\n 1-continue to show plots, 0-stop plots\n Make your choise:");
			scanf("%i", &pausav);
			plots_flag=pausav;
                        int nprocs=1;
                        if (plots_flag==0){
                                nprocs=omp_get_num_procs();
                                omp_set_num_threads( nprocs );
                        }
                }
                
                
                if (status2 == 0) {
                    aponta[ilinha*9+0]=linha;
                    aponta[ilinha*9+1]=ncenter;
                    aponta[ilinha*9+2]=-acoef[3*hjl];
                    aponta[ilinha*9+3]=2.*sqrt(log(2)/acoef[3*hjl+1]);
                    aponta[ilinha*9+4]=medida;
                    aponta[ilinha*9+5]=acoef[3*hjl];
                    aponta[ilinha*9+6]=acoef[3*hjl+1];
                    aponta[ilinha*9+7]=acoef[3*hjl+2];
                    aponta[ilinha*9+8]=medida_er;
                } else aponta[ilinha*9+4]=-1;
                
                //Escrever no ficheiro de Log:
                pFile3 = fopen ("logARES.txt","a");
                    fprintf(pFile3,"%s%s%s%s%s%s",strLinhaInicial,strLinhaFound,strLinhaResample,strLinhaGuess,strLinhaFitted,strLinhaResult);
                fclose (pFile3);
                
            } else {
                printf("\n line not found\n");
                pFile3 = fopen ("logARES.txt","a");
                fprintf(pFile3,"%s%s",strLinhaInicial,"\n line not found\n");
                fclose (pFile3);
                aponta[ilinha*9+4]=-1;
            }
        *plots_flag2=plots_flag;

}

void arraysubcp(double retarr[], double arr[],long a, long b){
    long i;
    for (i=0; i < b-a; i++)
        retarr[i]=arr[a+i];
}



void poly_fitn(double xvec[], double yvec[], double err[], long n, long ord, double coefs[]) {
    int i, j, k;
    double xi, yi, ei, chisq,xi2;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;
    ord++;
    X = gsl_matrix_alloc (n, ord);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    c = gsl_vector_alloc (ord);
    cov = gsl_matrix_alloc (ord, ord);
    for (i = 0; i < n; i++) {
      xi=xvec[i];
      yi=yvec[i];
      ei=err[i];
      for (j = 0; j < ord; j++) {
        xi2=1.0;
        for (k=0; k<j; k++) xi2*=xi;
        gsl_matrix_set (X, i, j, xi2);
        }
      gsl_vector_set (y, i, yi);
      gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, ord);
    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    #define C(i) (gsl_vector_get(c,(i)))
    #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    for (j = 0; j < ord; j++)
        coefs[j]=C(j);
    
}


void deriv(double x[], double y[], double dy[], long n) {
    int i;
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp *interp = gsl_interp_alloc (gsl_interp_cspline, n);
    //gsl_interp *interp = gsl_interp_alloc (gsl_interp_akima, n);
    gsl_interp_init (interp, x, y, n);
    for (i=0; i<n; i++)
        dy[i]=gsl_interp_eval_deriv (interp, x, y,x[i],acc);
    gsl_interp_free (interp);
    gsl_interp_accel_free (acc);
}


int continuum_det5 (double x[], double y[], double ynorm[], long nxele, double res[], double tree, int plots_flag){
//    clean_zero_gaps(y,nxele);
    int order,i,j;
    order=2;
    double err[nxele], coefs[order+1];
    long nvec;

    for(i=0;i<nxele;i++)
        err[i]=1.;
//    plotxy(x,y,nxele,x[0],x[nxele-1]);
    poly_fitn(x,y,err,nxele,order,coefs);

    double xi=1.;
    for(i=0;i<nxele;i++) {
        ynorm[i]=0.;
        xi=1.;
        for (j=0;j<order+1;j++) {
            ynorm[i]+=coefs[j]*xi;
            xi*=x[i];
        }
    }
    double vecx[nxele],vecy[nxele];
    int jk;
    for (jk=0; jk<5; jk++) {
        nvec=0;
        for (i=0; i<nxele-1; i++) {
    // testes foram feitos com 0.01, nao deve causar problemas por maior. Puz 0.1
            if (y[i] > ynorm[i]*tree && fabs(y[i]-y[i+1]) < 0.1*y[i]){
                vecx[nvec]=x[i];
                vecy[nvec]=y[i];
                nvec++;
            }
        }
	    if (nvec <= 2)
	    	return -1;
        poly_fitn(vecx,vecy,err,nvec,order,coefs);

        for(i=0;i<nxele;i++) {
            ynorm[i]=0.;
            xi=1.;
            for (j=0;j<order+1;j++) {
                ynorm[i]+=coefs[j]*xi;
                xi*=x[i];
            }
        }
    }

    for (i=0; i < order+1; i++)
        res[i]=coefs[i];

    res[3]=0.;
    if (plots_flag == 1)
        plotxyover3(x,y,nxele,x,ynorm,nxele,vecx,vecy,nvec,x[0],x[nxele-1]);
    
    for (i=0; i<nxele; i++)
	ynorm[i]=y[i]/(res[0]+res[1]*x[i]+res[2]*x[i]*x[i]+res[3]*x[i]*x[i]*x[i]);


return 1;

}


void smooth(double vec[], long n, int w, double svec[]) {
    int i,j;
    double soma;

    if (w%2 != 1)
            w++;
    for (i=0; (i < (w-1)/2);i++)
        svec[i]=vec[i];
    for (i=(w-1)/2;i<n-((w-1)/2);i++) {
        soma=0.;
        for (j=i-((w-1)/2); j<=i+((w-1)/2);j++)
            soma+=vec[j];
        svec[i]=soma/w;
    }
    for (i=n-((w-1)/2); i<n;i++)
        svec[i]=vec[i];
}


void zeroscenterfind(double y[], double iy[], double dy[], double ddy[], long n, long center[], long *ncenter) {
    double zerostot[n], contot[n], tutezerostot[n][2], maxdy;
    long ntot=0, nctot=0, ctot=0, i, centertot[n];
    int signal=0, signalc=0, signal_ant, signalc_ant;
    if (y[0] == abs(y[0]))
        signal=1;
    if (ddy[0] == abs(ddy[0]))
        signalc=1;
    signal_ant=signal;
    signalc_ant=signalc;
    maxdy=maxele_vec(dy,n);

    for (i=0; i<n; i++) {
        signalc=0;
        if ( (float) ddy[i] == fabs( (float) ddy[i]) )
            signalc=1;
        // quando muda de sinal, é um maximo local na 2a derivada esta abaixo do ruido e a 3 derivada já negativa o suficiente (devido a oscilacao do ruido)
        // no 0.98 a ideia era ter o tree, mas a coisa nao funcionava bem. Identicaria muitas riscas para o caso de termos bom S/N
        // Assim so aceitamos riscas identificadas que tenham uma dept de pelo menos 0.98
        if ( (signalc != signalc_ant) && (dy[i] > 0.01*maxdy) && (iy[i] < 0.98) && (ddy[i] < -0.1) ) {
            centertot[ctot]=i;
            ctot++;
        }
//        signal=0;
//        if ( (float) y[i] == (float) fabs(y[i]))
//            signal=1;
//        if ( signal != signal_ant) {
//            tutezerostot[ntot+nctot][0]=i;
//            if (iy[i] < 0.98) {
//                zerostot[ntot]=i;
//                if (dy[i] <= 0)		tutezerostot[ntot+nctot][1]=0;
//                else			tutezerostot[ntot+nctot][1]=0.5;
//                ntot++;
//            } else {
//                contot[nctot]=i;
//                tutezerostot[ntot+nctot][1]=1.;
//                nctot++;
//            }
//        }
//        signal_ant=signal;
        signalc_ant=signalc;

    }


    if (ctot != 0) {
        *ncenter=ctot;
        for (i=0;i<ctot;i++) 	center[i]=centertot[i];
    } else {
        center[0]=-1;
        *ncenter=0;
    }
}


double maxele_vec(double vec[], long nvec) {
    long i;
    double maxi=vec[1+nvec/20];
    for (i=1+nvec/20; i<nvec-nvec/20; i++)
        maxi = max(maxi,vec[i]);
    return maxi;
}


void fitngauss(double t[], double y[], double sigma[], long nvec, double acoef[], double acoef_er[], int para, int *status2)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  int status;
  size_t i, iter = 0;
  long N=nvec;
  const size_t n = N;
  const size_t p = para;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double dy[N];
  struct data d = { n, para, t, y, sigma};
  gsl_multifit_function_fdf f;

  double x_init[para];
  for (i=0; i<para; i++)
	x_init[i]=acoef[i];

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  T = gsl_multifit_fdfsolver_lmder;

  s = gsl_multifit_fdfsolver_alloc (T, n, p);

  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      int i=0;
      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-6, 1e-6);
    }
  while (status == GSL_CONTINUE && iter < 5000);

  gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  {

    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));

    for (i=0; i<para; i++)
	{
	acoef[i]=FIT(i);
	acoef_er[i]=c*ERR(i);
        }

  printf ("status = %s\n", gsl_strerror (status));
  printf ("int status= %i", status);
  *status2=status;
  *status2=0;   //sometimes we have bad fits but the result is perfectably acceptable
  }
  gsl_multifit_fdfsolver_free (s);
}

#ifdef	__cplusplus
}
#endif

#endif	/* _ARESLIB_H */

