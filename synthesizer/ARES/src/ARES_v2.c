/* 
 * File:   ARES_v2.c
 * Author: Sérgio Sousa
 *
 * Created on Julho 31, 2014
 */
// ARES v2.0 Automatic Routine for Equivalent widths of Spectra
// method and implementation :Sergio G. Sousa (IA - CAUP) sousasag@astro.up.pt
// 
// Based on the original code ARES from Sergio Sousa.
// supervised by Nuno C. Santos and Mario J. P. F. G. Monteiro
//
// v.2.0 - fully compatible with ares v.1.0 
// New features:
//   - parallelization with OpenMP 
//   - reads spectra in ascii file (still is recommended to be evenly spaced in wavelength
//   - cleaning of zeros in the spectra - to avoid crashing when looking for lines in spectral gaps. 
//   - Radial velocity Correction is possible within ARES - either estimated or using RV value
//   - new option to estimate rejt value automatically
//   - new option to use rejt depending on wavelenght (need to set up file: 'lambda_rejt.opt')
//   - estimation of the error on the EWs. based on the S/N of the spectra and the errors on the fit of the Gaussians
// 
//
// Dependences: Gnu Scientific Library (gsl); CFITSIO; OpenMP; PlotUtils to show simple plots
//
// Compilations instructions:
// If the libraries are well set up in your system, simply:
//     gcc -o ARES ARES_v2.c -lcfitsio -lgsl -lgslcblas -lm -lgomp -fopenmp
// 
// Otherwise need to define the path of the libraries: Example:
// gcc -o ARES ARES_v2.c -L/usr/local/cfitsio/lib/ -I/usr/local/cfitsio/include/ -lcfitsio -lgsl -lgslcblas -lm -lgomp -fopenmp
// 
// See the paper and the webpage for more information on the code.
// You can also contact the author for any questions/suggestions related with ARES: sousasag@astro.up.pt
// 
// Ps. Most of the code is not fully commented, and still mostly in Portuguese (sorry for that)


#include <omp.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <gsl/gsl_rng.h> 

#include "areslib.h"
#include "rvcor.h"
#include "sn_rejt_estimator.h"


int main() {

/*Declaraçao de variaveis*/


    FILE * pFile3;
    int status2; 

    long npoints;
    double *pixels, *xpixels, cdelta1, crval1;
/* As minhas variaveis*/
    char filetest[200], fileleitura[200], fileout[200], rvmask[200], tree[200];
    double lambdai, lambdaf, smoothder, rejt, distlinha, miniline;
    double space;
    int plots_flag, plots_stop,nchar;
    int spectrum_type=-1;

/* FIM DE Declaracao de variaveis*/

    
/* leitura das opcoes , mine.opt file read*/

    read_mine("mine.opt", filetest, fileleitura, fileout, &lambdai, &lambdaf, &smoothder, &space, tree, &distlinha, &miniline, &plots_flag, rvmask);

    
    
    plots_stop=0;
    if (plots_flag == 2)  {
            plots_flag=1;
            plots_stop=1;
    }

    printf("File: %s\n", filetest);
    char* ext=getFileExtension(filetest);
    spectrum_type=set_tipo_espectro(ext);

    printf("Extension: %s - %d\n",ext, spectrum_type );

    read_spectrum_file(filetest, &npoints, &pixels, &xpixels, &cdelta1, &crval1, spectrum_type);

    printf("Cleaning zero gaps\n");
    clean_zero_gaps(pixels,npoints);
    
    double radvel=get_rv(xpixels, pixels, npoints, rvmask);        
    printf("Velocidade radial: %f\n", radvel);
    printf("ponto %d, lambda: %f \n",1000,xpixels[1000]);
    correct_lambda(xpixels,npoints,radvel);
    crval1=xpixels[0];
    printf("ponto %d, lambda: %f \n",1000,xpixels[1000]);
    
    rejt=get_rejt(tree, xpixels,pixels, npoints);
    if (rejt == -2)
        printf("USING lambda_rejt.opt file to rejt dependence on lambda\n");
    else {
        if (rejt <= 0)      exit(2);
    }
    printf("Rejt used for computations: %f\n",rejt);
    
    
// verificação dos limites do espectro
    if (lambdai < xpixels[0]) lambdai=xpixels[0];
    if (lambdaf > xpixels[npoints-1]) lambdaf=xpixels[npoints-1];
    
/* leitura do ficheiro das riscas a medir laboratory.dat*/

    long nl;
    double* linhas;
    read_lines_list(fileleitura, &nl, &linhas);
    pFile3 = fopen ("logARES.txt","a");
    fprintf(pFile3,"\nLog Result ARES...\n");
    fprintf(pFile3,"Velocidade radial: %f\n", radvel);
    if (rejt == -2)
        fprintf(pFile3,"Rejt used for computations: USING lambda_rejt.opt file to rejt dependence on lambda\n");
    else
        fprintf(pFile3,"Rejt used for computations: %f\n",rejt);
    fprintf(pFile3,"Updated limit spectra: [%f - %f ]\n",lambdai,lambdaf);    
    fclose (pFile3);
    

    double aponta[nl*9];
    
    int nprocs=1;
    omp_set_num_threads( nprocs );
    
    if (plots_flag==0){
        nprocs=omp_get_num_procs();
        omp_set_num_threads( nprocs );            
    }

    double rejtin=rejt;
    int i=0,fgh;

/*
// Tried to shufle the lines to see if performence improved.
// May improve in some senses. Output is missing 2-3 lines don't know the reason
    const gsl_rng_type * T;
    T = gsl_rng_default;
    gsl_rng * r;
    r = gsl_rng_alloc (T);
    int ir[nl];
    for (i=0;i<nl;i++){
        ir[i]=i;
    }
    gsl_ran_shuffle (r, ir, nl, sizeof(int));
    for (i=0;i<nl;i++){
        printf("%d\n", ir[i]);
    }
    gsl_rng_free (r);
*/
    
    #pragma omp parallel for 
//    for (i=0;i<nl;i++) {
//        fgh=ir[i];        
    for (fgh=0;fgh<nl;fgh++) {
        float linha=linhas[fgh];
        
        if ( (linha > lambdai+2*space) && (linha < lambdaf-2*space) )  {        //não mede riscas nos limites do espectro
            if (rejt == -2) {
                rejtin = get_rejt_lambda_file(linha);
                printf("\nUsing rejt = %f for line: %f\n", rejtin, linha);
            }
            getMedida(xpixels, pixels, linha, space, rejtin, &plots_flag, smoothder, distlinha, pFile3, fgh, aponta, lambdai, lambdaf);
                
        } else aponta[fgh*9+4]=-1;
    }

//  Write of the result in the output file
    write_outfile(fileout, aponta, nl, miniline);


//  Some free memory
    system("rm tmp tmp2 tmp3 tmp20 tmp22 tmp23");
    free(pixels);
    free(xpixels);
    
    printf("\n FINISH. HAVE A GOOD DAY...\n");
    return 0;
}

