/* 
 * File:   rvcor.h
 * Author: sousasag
 *
 * Created on April 20, 2012, 11:45 AM
 */

#ifndef RVCOR_H
#define	RVCOR_H

#ifdef	__cplusplus
extern "C" {
#endif


#include "areslib.h"    
#include <gsl/gsl_interp.h>

struct maskstruct{
    double center;
    double space;
    double* lines;
    int nl;
};    
    
void correct_lambda(double * lambda , int np, double vrad);
int get_local_rvo(double *lambda, double* flux, long np, double * mask, int nmask, double* rvel);
int get_local_rv(double *lambda, double* flux, long np, struct maskstruct mask, double* rvel);
double get_rv(double *lambda, double* flux, long np, char* rvmask);
void interpollin(double* x,double* y,double* xi, double* yi, int n, int ni);


double get_rv(double *lambda, double* flux, long np, char* rvmask){
    double radvel=0;
    char *pch;
    int rvflag;
    pch = strtok (rvmask,",");
    rvflag=atoi(pch);
    
    if (rvflag == 0) {
        pch = strtok (NULL, " ");
        radvel=atof(pch);
        return radvel;
    } else {
        struct maskstruct mask;        
        double masklocal[rvflag];
        int i;
        for (i=0; i<rvflag; i++){
                pch = strtok (NULL,",");
                masklocal[i]=atof(pch);
        }
        pch = strtok (NULL,",");
        mask.center=atof(pch);
        pch = strtok (NULL, " ");
        mask.space=atof(pch);;
        mask.nl=rvflag;
        mask.lines=masklocal;
        int flag=0;
        flag=get_local_rv(lambda, flux, np, mask, &radvel);        
        if (flag == 1) 
            return radvel;
        else {
            printf("Radial Velocity not determined.\n");
            exit(1);
        }
    }
}



void correct_lambda(double * lambda, int np, double vrad){
    //vrad in Km/s
    int i;
    for (i=0; i<np; i++){
        lambda[i]=lambda[i]/(1+vrad/3.e5);
    }
}
    
//mask
//0-> center
//1-> space
//2--n+2-> center of the lines in the mask
// return {-1,0} if no radial velocity derived 
int get_local_rvo(double *lambda, double* flux, long np, double * mask, int nmask, double* rvel){
    int radvel[2];
    radvel[0]=-1;
    *rvel=0;
    // we can change the space accordingly with the mask
    double linecenter=mask[0];
    double spacerv=mask[1];
    
    double lambdai=lambda[0];
    double lambdaf=lambda[np-1];
    
    // get center line for mask
    
    if ( (linecenter > lambdai+2*spacerv) && (linecenter < lambdaf-2*spacerv) ){
            long indexi=find_pixel_line(lambda, linecenter-spacerv);
            long indexf=find_pixel_line(lambda, linecenter+spacerv);
            long nl=indexf-indexi;
            double lambda_loc[nl], flux_loc[nl];
            arraysubcp(lambda_loc, lambda,indexi,indexf);
            arraysubcp(flux_loc , flux,indexi,indexf);
            // continuum determination
            double flux_loc_norm[nl];
            double res[4];
            continuum_det5(lambda_loc,flux_loc,flux_loc_norm,nl,res,0.985,1);

            int nlm=nmask-2;
            double lambdaccfi[nlm];
            arraysubcp(lambdaccfi,mask,2, nmask);
            double lambdaccfv[nlm];
            double lambdaflxv[nlm];
            int nccf=8000;  //number of points for the local ccf
            double ccfvel[nccf];
            double ccfflx[nccf];
            double vmin=0; double cflxmin=200000;
            double stepv=0.1; //0.1 Km/s
            double inivel=-400;
            int i;
            for (i=0;i<nccf;i++){
                    ccfvel[i]=i*stepv+inivel;
                    int oo;
                    for (oo=0; oo < nlm; oo++)
                        lambdaccfv[oo]=lambdaccfi[oo]*(1+ccfvel[i]/3.e5);

                    interpollin(lambda_loc,flux_loc_norm,lambdaccfv,lambdaflxv,nl,nlm);
                    ccfflx[i]=lambdaflxv[0]+lambdaflxv[1]+lambdaflxv[2];

//                    printf("%d - %f, %f - %f -- %f - %f ::: %f\n",i,ccfvel[i], lambdaccfv[0], lambdaccfv[1], lambdaflxv[0], lambdaflxv[1],ccfflx[i] );
                    
                    if (ccfflx[i] < cflxmin) { 
                        vmin=ccfvel[i]; 
                        cflxmin=ccfflx[i];
                    }
            }
            
//            plotxy(ccfvel, ccfflx, nccf, -400., 400.);
            printf("Velocidade radial: %f\n", vmin);
            radvel[0]=1;
            *rvel=vmin;
            
    }
    
    return radvel[0];
}
    

int get_local_rv(double *lambda, double* flux, long np, struct maskstruct mask, double* rvel){
    int flag=-1;
    *rvel=0;
    // we can change the space accordingly with the mask
    double linecenter=mask.center;
    double spacerv=mask.space;
    
    double lambdai=lambda[0];
    double lambdaf=lambda[np-1];
    
    // get center line for mask
    
    if ( (linecenter > lambdai+2*spacerv) && (linecenter < lambdaf-2*spacerv) ){
            long indexi=find_pixel_line(lambda, linecenter-spacerv);
            long indexf=find_pixel_line(lambda, linecenter+spacerv);
            long nl=indexf-indexi;
            double lambda_loc[nl], flux_loc[nl];
            arraysubcp(lambda_loc, lambda,indexi,indexf);
            arraysubcp(flux_loc , flux,indexi,indexf);
            // continuum determination
            double flux_loc_norm[nl];
            double res[4];
            continuum_det5(lambda_loc,flux_loc,flux_loc_norm,nl,res,0.985,0);

            int nlm=mask.nl;
            double lambdaccfi[nlm];
            arraysubcp(lambdaccfi,mask.lines,0, mask.nl);
            double lambdaccfv[nlm];
            double lambdaflxv[nlm];
            int nccf=8000;  //number of points for the local ccf
            double ccfvel[nccf];
            double ccfflx[nccf];
            double vmin=0; double cflxmin=200000;
            double stepv=0.1; //0.1 Km/s
            double inivel=-400;
            int i;
//            printf("i  ccfvel[i] lambdaccfv[0] lambdaccfv[1] lambdaflxv[0] lambdaflxv[1] ccfflx[i]\n" );
            for (i=0;i<nccf;i++){
                    ccfvel[i]=i*stepv+inivel;
                    int oo;
                    for (oo=0; oo < nlm; oo++)
                        lambdaccfv[oo]=lambdaccfi[oo]*(1+ccfvel[i]/3.e5);

                    interpollin(lambda_loc,flux_loc_norm,lambdaccfv,lambdaflxv,nl,nlm);
                    ccfflx[i]=lambdaflxv[0]+lambdaflxv[1]+lambdaflxv[2];
                    ccfflx[i]=0;
                    for (oo=0; oo < nlm; oo++)
                        ccfflx[i]+=lambdaflxv[oo];
//                    printf("%d   %f   %f   %f   %f   %f   %f\n",i,ccfvel[i], lambdaccfv[0], lambdaccfv[1], lambdaflxv[0], lambdaflxv[1],ccfflx[i] );
                    
                    if (ccfflx[i] < cflxmin) { 
                        vmin=ccfvel[i]; 
                        cflxmin=ccfflx[i];
                    }
            }
            
//            plotxy(ccfvel, ccfflx, nccf, -400., 400.);
//            int pausav;
//            scanf("%i", &pausav);
//            printf("Velocidade radial: %f\n", vmin);
            flag=1;
            *rvel=vmin;
    }
    
    return flag;
}


void interpollin(double* x,double* y,double* xi, double* yi, int n, int ni){
      gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, n);

      gsl_interp_init(interpolation, x, y, n);
      gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

      //get interpolation for x = 1981
	int i;
	for (i=0; i<ni; i++){
		yi[i] = gsl_interp_eval(interpolation, x, y, xi[i], accelerator); 
//		printf("%f\n",yi[i]);
	}
      gsl_interp_free (interpolation);
      gsl_interp_accel_free (accelerator);
}
  
    
    
    
    


#ifdef	__cplusplus
}
#endif

#endif	/* RVCOR_H */

