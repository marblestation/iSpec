/* 
 * File:   sn_rejt_estimator.h
 * Author: sousasag
 *
 * Created on December 20, 2012, 3:50 PM
 */

#ifndef SN_REJT_ESTIMATOR_H
#define	SN_REJT_ESTIMATOR_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_interp.h>
#include "areslib.h"
#include "filesio.h"

#define SMOOTH_NOISE 8

double avg(double* array, long nelem);
double median(double* array, long nelem);
double stddev(double* array, long nelem);
void remove_outlier(double* arrayin, long nin,double* arrayout,long* nout, double out_sigma);
double compute_sn_range(double li,double lf,double* lambda, double* flux, long np);
double get_rejt(char* tree, double* lambda, double* flux, long np);
double get_rejt_from_sn(double sn);
double get_rejt_lambda_file(double lambda);

double get_rejt_lambda_file(double lambda){
    double *lambdavec, *rejtvec;   
    double cdelta1mean=0;
    double crval1=0;
    long npoints=0;
    char filetest[50] = "lambda_rejt.opt";
    read_ascii_file(filetest, &npoints, &rejtvec, &lambdavec, &cdelta1mean, &crval1);
/*    int i;
    for (i=0; i<npoints;i++){
        printf("%d - %f - %f\n", i,lambdavec [i], rejtvec[i]);
    }
*/
    
    gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,npoints);
    gsl_interp_init(interpolation, lambdavec, rejtvec, npoints);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    double rejt_lambda = gsl_interp_eval(interpolation,lambdavec, rejtvec, lambda, accelerator);

//    printf("%f  -  %f  \n",lambda,rejt_lambda);
   
    gsl_interp_free (interpolation);
    gsl_interp_accel_free (accelerator);
    
    free(lambdavec);
    free(rejtvec);
    return rejt_lambda;

}

double get_rejt_from_sn(double sn){
    return 1.-(1./sn);
}





double get_rejt(char* tree, double* lambda, double* flux, long np){
    double rejt=0;
    char* pch;
    int i;
    if (strchr(tree, ';') == NULL){  // do stuff
        rejt=atof(tree);
        if (rejt >= 1) rejt=get_rejt_from_sn((double)rejt);
        printf("Computed rejt: %f \n", rejt);
    } else {
        printf("Computing rejt for spectrum...");
        char temp[200];
        strcpy(temp, tree);
        printf("Tree: %s\n",temp);
        pch = strtok (temp,";");
        int nranges=atoi(pch);
        printf("Nranges: %d\n", nranges);
        double snvec[nranges]; double snmax=0;
        for (i=0;i<nranges;i++){
            pch = strtok (NULL,",");
            double li=atof(pch);
            pch = strtok (NULL,",");
            double lf=atof(pch);
//            printf("Range %d: [%f - %f]\n",i,li,lf);
            snvec[i]=compute_sn_range(li,lf,lambda,flux,np);
            printf("Range %d: [%f - %f] : SN: %f\n",i,li,lf,snvec[i]);
            if (snvec[i] > snmax) snmax=snvec[i];
        }
        int ranges_valid=0;
        for (i=0; i<nranges; i++){
            snvec[ranges_valid]=snvec[i];
            if (snvec[ranges_valid] > 0) ranges_valid++;
        }
//        printf("Ranges Valid: %d\n",ranges_valid);
        double mediansn=median(snvec,ranges_valid);
        printf("Median SN: %f\n", mediansn);
//        printf("Max SN:    %f\n", snmax);
        rejt=get_rejt_from_sn(mediansn);
    }
    return rejt;
}


double compute_sn_range(double li,double lf,double* lambda, double* flux, long np){
    double sn=-1.0;
    if ( lambda[0] < li && lf < lambda[np-1]){
        long pi=find_pixel_line(lambda, li);
        long pf=find_pixel_line(lambda, lf);
        long npr=pf-pi;
        double fluxrange[npr],lambdarange[npr],fluxrangesmooth[npr],noise[npr],temp[npr],noise_clean[npr];
        arraysubcp(fluxrange,flux,pi,pf);
        arraysubcp(lambdarange,lambda,pi,pf);
        smooth(fluxrange,npr,SMOOTH_NOISE,fluxrangesmooth);
        double average=avg(fluxrange,npr);
        int i;
        for (i=0;i<npr;i++) 
            noise[i]=fluxrange[i]-fluxrangesmooth[i]+average;
        long n_clean=0;
        average=avg(noise,npr);
        double sigma=stddev(noise,npr);
//        printf("Average: %f , Stddev: %f , SN: %f \n",average, sigma, average/sigma);
//        remove_outlier(noise,npr,noise_clean,&n_clean,2.0);
//        average=avg(noise_clean,n_clean);
//        sigma=stddev(noise_clean,n_clean);
//        printf("Average: %f , Stddev: %f , SN: %f \n",average, sigma, average/sigma);
        sn=average/sigma;
    } else {
        sn=-1.0;
    }
    return sn;
}


double avg(double* array, long nelem) {
    return gsl_stats_mean(array, 1, nelem);
}

double median(double* array, long nelem){
    gsl_sort(array, 1, nelem);
    return gsl_stats_median_from_sorted_data(array,1,nelem);
}

double stddev(double* array, long nelem){
    return sqrt(gsl_stats_variance(array, 1,nelem));
}


void remove_outlier(double* arrayin, long nin,double* arrayout,long* nout, double out_sigma){
    double meanvec=avg(arrayin,nin);
    double sigmavec=stddev(arrayin,nin);
    long i;
    long ind=0;
    double out=out_sigma*sigmavec;
    for (i=0; i<nin; i++){
        if ( fabs(arrayin[i]-meanvec) <= out) {
            arrayout[ind]=arrayin[i];
            ind++;
        }
    }
//    printf("Total Points: %ld Pontos in: %ld Pontos out: %ld\n",nin, ind, nin-ind);
    *nout=ind;
}



#ifdef	__cplusplus
}
#endif

#endif	/* SN_REJT_ESTIMATOR_H */

