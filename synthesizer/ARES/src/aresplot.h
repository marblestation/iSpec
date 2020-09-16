/* 
 * File:   areslib.h
 * Author: sousasag
 *
 * Created on April 25, 2011, 2:58 PM
 */

#ifndef _ARESPLOT_H
#define	_ARESPLOT_H

// 1- for plot_utils
// 2- for gnuplot
// 3- for gnuplot saving png plots in plotdir
#define PLOT_TYPE 2


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

#ifdef	__cplusplus
extern "C" {
#endif
    
void plotxy(double *, double *, long, double, double);
void plotxyover(double *, double *, long, double *, double *, long, double, double);
void plotxyover2(double *, double *, long, double *, double *, long, double, double);
void plotxyover3(double *, double *, long, double *, double *, long, double *, double *, long, double, double);
void create_dir_plot();


void create_dir_plot(){

    if (PLOT_TYPE == 3){
		struct stat st = {0};

		if (stat("plotdir", &st) == -1) {
    		mkdir("plotdir", 0700);
		} else{
        system("rm -rf plotdir");
        mkdir("plotdir", 0700);
		}
    }
}

void plotxy(double xvec[], double yvec[], long np, double xi, double xf){
    FILE * pFile2;
    long t;
    char str[200];
    pFile2 = fopen ("tmp","wt");
    for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
    fclose (pFile2);
    char *buffer;
    int decimal, sign;
    buffer = fcvt (xi, 0, &decimal, &sign);
//    strcpy (str,"graph -T X -x -");
    strcpy (str,"graph -T X -x ");
    strcat (str,buffer);
    strcat (str," ");
    buffer = fcvt (xf, 0, &decimal, &sign);
    strcat (str,buffer);
    strcat (str," < tmp");
    system(str); 
/*
    Py_Initialize();
      PyRun_SimpleString("import pylab");
      PyRun_SimpleString("import numpy");
      PyRun_SimpleString("x,y = numpy.loadtxt('tmp')");
      PyRun_SimpleString("pylab.plot(x,y)");
      PyRun_SimpleString("pylab.show()");
    Py_Exit(0);
*/    
    printf("%s\n",str);
//    system("rm tmp");
}

void plotxyover(double xvec[], double yvec[],long np, double xvec2[], double yvec2[],long np2,double xi, double xf){
    FILE * pFile;
    FILE * pFile2;
    long t;
    char str[200];
    pFile2 = fopen ("tmp","wt");
    for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
    fclose (pFile2);
    pFile = fopen ("tmp2","wt");
    for(t=0;t<np2;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec2[t],yvec2[t]);
    fclose (pFile);
    char *buffer;
    int decimal, sign;
    buffer = fcvt (xi, 0, &decimal, &sign);
    strcpy (str,"graph -T X -x ");
    strcat (str,buffer);
    strcat (str," ");
    buffer = fcvt (xf, 0, &decimal, &sign);
    strcat (str,buffer);
    //	strcat (str," tmp -S 2 -C -m 42 tmp2");
    strcat (str," tmp -S 1 -C -m 43 tmp2");
    system(str);
    system("rm tmp tmp2");
}

void plotxyover2(double xvec[], double yvec[],long np, double xvec2[], double yvec2[],long np2,double xi, double xf){
    FILE * pFile;
    FILE * pFile2;
    FILE * pFile3;
    long t;
    char str[200];
    pFile2 = fopen ("tmp","wt");
    for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
    fclose (pFile2);
    pFile = fopen ("tmp2","wt");
    for(t=0;t<np2;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec2[t],yvec2[t]);
    fclose (pFile);
    pFile3 = fopen ("tmp3","wt");
    fprintf(pFile3," %15.6f  %15.6f \n",xi,1.);
    fprintf(pFile3," %15.6f  %15.6f \n",xf,1.);
    fclose (pFile3);


    if (PLOT_TYPE == 1) {
    	char *buffer;
    	int decimal, sign;
    	buffer = fcvt (xi, 0, &decimal, &sign);
    	strcpy (str,"graph -T X -x ");
    	strcat (str,buffer);
    	strcat (str," ");
    	buffer = fcvt (xf, 0, &decimal, &sign);
    	strcat (str,buffer);
    	strcat (str," tmp tmp3 -S 1 -C -m 43 tmp2 ");
    	system(str);
    } else {
//  GNUPLOT:
		if (PLOT_TYPE == 3) {
            double line = xi + (xf-xi)/2.;
			FILE *pipe = popen("gnuplot -persist","w");
			fprintf(pipe, "set term png \n");	
			fprintf(pipe, "set output \"plotdir/plot_%08.2f.png\" \n", line);	
			fprintf(pipe, "plot 'tmp' with lines, 'tmp2' with lines, 'tmp3' with lines \n");	
			fprintf(pipe, "set term x11 \n");	
			pclose(pipe);
            char str_norm[300];
            sprintf(str_norm, "cp tmp plotdir/spec_%08.2f.dat", line);
            system(str_norm);
		} else { 
			FILE *pipe = popen("gnuplot -persist","w");
			fprintf(pipe, "plot 'tmp' with lines, 'tmp2' with lines, 'tmp3' with lines \n");	
			pclose(pipe);
		}
	}

//    system("rm tmp tmp2 tmp3");


}

void plotxyover3(double xvec[], double yvec[],long np, double xvec2[], double yvec2[], long np2, double xvec3[], double yvec3[], long np3, double xi, double xf){
    FILE * pFile;
    FILE * pFile2;
    FILE * pFile3;
    long t;
    char str[200];
    pFile2 = fopen ("tmp20","wt");
    for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
    fclose (pFile2);

    pFile = fopen ("tmp22","wt");
    for(t=0;t<np2;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec2[t],yvec2[t]);
    fclose (pFile);

    pFile3 = fopen ("tmp23","wt");
    for(t=0;t<np3;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec3[t],yvec3[t]);
    fclose (pFile3);


    if (PLOT_TYPE == 1) {
    	char *buffer;
    	int decimal, sign;
    	buffer = fcvt (xi, 0, &decimal, &sign);
    	strcpy (str,"graph -T X -x ");
    	strcat (str,buffer);
    	strcat (str," ");
    	buffer = fcvt (xf, 0, &decimal, &sign);
    	strcat (str,buffer);
    	//	strcat (str," tmp -S 2 -C -m 42 tmp2 -S 4 -C -m 0 tmp3");
    	strcat (str," tmp20 -S 1 -C -m 43 tmp22 -S 1 -C -m 0 tmp23");
	//    printf("%s\n",str);
	    system(str);
	} else {
	//  GNUPLOT:
		FILE *pipe = popen("gnuplot -persist","w");
		fprintf(pipe, "plot 'tmp20' with lines, 'tmp22' with lines, 'tmp23' with points\n");	
		pclose(pipe);  
	//    system("rm tmp20 tmp22 tmp23");
		
	}

}




#ifdef	__cplusplus
}
#endif

#endif	/* _ARESPLOT_H */

