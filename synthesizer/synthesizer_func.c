/**
    This file is part of Spectra.
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
    
    Spectra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spectra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Spectra.  If not, see <http://www.gnu.org/licenses/>.
**/
/*
    This file has been adapted from spectrum.c:
     
        SPECTRUM, a Stellar Spectral Synthesis Program
        (C) Richard O. Gray 1992 - 2010 Version 2.76e
        May 3, 2010
     
    It works always on the Integrated Disk mode (normalized Intensity) 
*/
#include <stdio.h>
#include "spectrum276e/spectrum.h"
#include "synthesizer_func.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

void inmodel(atmosphere *model, char *file, int flagw);
void Density(atmosphere *model, atominfo *atom, double ah, double ahe, int flagw);
void hotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw);
void veryhotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw);
void tauref(atmosphere *model, double wavref);
void tauwave(atmosphere *model, double wave);
void tauflx(atmosphere *model, double wstart, double wend);
double flux(atmosphere *model, double wave);
void fluxflx(atmosphere *model, double wstart, double wend, double *flxstart, double *flxend);
double intensity(atmosphere *model, double wave, double mu);
void intenint(atmosphere *model, double wstart, double wend, double *intenstart, double *intenend, double mu);
void inlin(double wave, linedata *line, int *nline, linelist *list, int nlist, atominfo *atom, isodata *isotope);
void pop(linedata *line, int N, atmosphere *model, pfunc *V, population *POP);
void broad(atmosphere *model,linedata *line,int N,double sig,double alp,double fac);
void capnu(linedata *line, int N, atmosphere *model);
void taukap(double wave, atmosphere *model, atominfo *atom, linedata *line, int N, linedata *strgln, pfunc *V, Helium *He, population *POP);
void inatom(char atmdat[], atominfo *atom, double MH, double *ah, double *ahe);
void linelst(double wave, linelist *list, int *nlist, atominfo *atom, double teff, double logg, FILE *qf, int reset, isodata *isotope, atmosphere *model, double Flux, pfunc *V, population *POP, double dwave);
void pfinit(pfunc *V, atominfo *atom, atmosphere *model, int flagw);
void popinit(population *POP, atominfo *atom, atmosphere *model, pfunc *V, int flagw);
double depth(atmosphere *model, double wave, double Flux);
double depthflx(atmosphere *model, double wave, double wstart, double wend);
double depthmu(atmosphere *model, double wave, double Intensity, double mu);
double idepthmu(atmosphere *model, double wave, double mu, double wstart, double wend);
double opacity(atmosphere *model, double lambda, int ntau);
float **cmatrix(int nrl,int nrh,int ncl,int nch);
void nrerror(char error_text[]);
char *ggets(char *s);
void inisotope(char isofile[], isodata *isotope);
void isorelabun(isodata *isotope);
void invelgrad(char vgrad[80]);
void infix(char fixfile[], atominfo *atom, double ah);
double interval(double wave, double dwave, int flagf);
void setreset(int k);
double dmax(double x, double y),dmin(double x, double y);
int imax(int x, int y),imin(int x, int y);

int synthesize_spectrum(char *atmosphere_model_file, char *linelist_file, char *abundances_file, double microturbulence_vel, int verbose, int num_measures, const double waveobs[], double fluxes[], progressfunc user_func, void *user_data) {
    int i;
    int nline = 0;
    int nlist = 0;
    
    
    int flagw = verbose; // Verbose mode (1: True, 0: False)
    int flaga = 0;
    int flagG = 0;
    int flagN = 0;
    int flagx = 0;
    double ah,ahe,waveref,wave,dwave=0,Flux,Depth;
    double vturb = 2.0e+5;
    atmosphere *model;
    atominfo *atom;
    linelist *list;
    linedata *line;
    linedata *strgln;
    isodata *isotope;
    pfunc *V;
    population *POP;
    Helium *He;
    char *file,name[80],lines[80],atmdat[80];
    char fixfile[80];
    char vgrad[80],isofile[80];
    FILE *qf;
    double opwave;

    setreset(0);
    ///////////////////////// Arguments
    strcpy(name, atmosphere_model_file); // stellar atmosphere data file
    file = name;
    strcpy(lines, linelist_file);
    strcpy(atmdat, abundances_file);        // atomic data file
    strcpy(isofile, "isotope.iso"); // only used if global flagI == 1 but then it produces segmentation fault (original SPECTRUM problem)
    //strcpy(oname, output_file); // output file
    //ofile = oname;
    vturb = microturbulence_vel; // km/s, if flagu == 0
    //start = wave_base; // Begin wavelength
    //end = wave_top; // End wavelength
    //dwave = wave_step; // step
    //inc = dwave;
    strcpy(vgrad,"velgrad.dat"); // velocity gradient file only used if flagg == 1
    /////////////////////////


    if(flagG == 1) flagg = 1;
    if(flagP == 1) flagp = 1;
    if(flagN == 1) flagw = 0;
    if((model = (atmosphere *) calloc(1,sizeof(atmosphere))) == NULL)
        nrerror("Allocation of memory for atmosphere failed");
    if((V = (pfunc *) calloc(1,sizeof(pfunc))) == NULL)
        nrerror("Allocation of memory for pfunc failed");
    if((atom = (atominfo *) calloc(NATOM,sizeof(atominfo))) == NULL)
        nrerror("Allocation of memory for atom failed");
    if((POP = (population *) calloc(NATOM-NMOL,sizeof(population))) == NULL)
        nrerror("Allocation of memory for population failed");
    if((line = (linedata *)
    calloc((unsigned long)NLINE,(unsigned long)sizeof(linedata))) == NULL)
        nrerror("Allocation of memory for line failed");
    if((oneline = (linedata *)
    calloc((unsigned long)1,(unsigned long)sizeof(linedata))) == NULL)
        nrerror("Allocation of memory for oneline failed");
    if((strgln = (linedata *)
        calloc(NSTRG,(unsigned long)sizeof(linedata))) == NULL)
        nrerror("Allocation of memory for caii failed");
    if((list = (linelist *) calloc(NLIST,sizeof(linelist))) == NULL)
        nrerror("Allocation of memory for list failed");
    if((He = (Helium *) calloc(NHE,sizeof(Helium))) == NULL)
        nrerror("Allocation of memory for He failed");
    if((velgrad = (float *) calloc(NTAU,sizeof(float))) == NULL)
        nrerror("Allocation of memory for velgrad failed");
    if((isotope = (isodata *) calloc(500,sizeof(isodata))) == NULL)
        nrerror("Allocation of memory for isotope failed");
    bkap = cmatrix(0,3,0,NTAU);
    bkap2 = cmatrix(0,3,0,NTAU);
    bkap3 = cmatrix(0,3,0,NTAU);
    bkap4 = cmatrix(0,3,0,NTAU);
    if(flagN != 1) {
        printf("\nSPECTRUM, a Stellar Spectral Synthesis Program"); 
        printf("\n(C) Richard O. Gray 1992 - 2010 Version 2.76e");
        printf("\nMay 3, 2010");
        printf("\n* Linked to Python by Sergi Blanco Cuaresma - February 2012\n");
        printf("\nIntegrated Disk mode (normalized Intensity)\n\n");
        if(flagc == 1) printf("Output will be continuum only (no line absorption)\n");
        if(flagw == 0) printf("Silent mode\n");
        if(flagg == 1) printf("Velocity gradient mode\n");
        if(flagI == 1) printf("Isotopes enabled\n");
        if(flagu == 1) printf("Reading microturbulent velocity from atmosphere model\n");
    }

    // stellar atmosphere data file
    inmodel(model,file,flagw);
    // Velocity gradient    
    if(flagg == 1) invelgrad(vgrad);
    if(flagg == 0) for(i=0;i<Ntau;i++) velgrad[i] = 0.0;
    
    // Enter name of line list file: (default = luke.lst)
    if((qf = fopen(lines,"r")) == NULL) {
        printf("Cannot find line data file\n");
        return(1);
    }
    if(flagI == 1) {
        // Enter name of isotope data file (default = isotope.iso)
        inisotope(isofile,isotope);
        isorelabun(isotope);
    }
    
    if(flagu == 0) {
        // microturbulence from km/s to cm/s
        vturb *= 1.0e+05;
        for(i=0;i<Ntau;i++) model->mtv[i] = vturb;
    }

    /* Normal abundances for hydrogen and helium; inatom may change these */
    ah = 0.911;
    ahe = 0.089;
    inatom(atmdat,atom,model->MH,&ah,&ahe);
    if(flaga == 1) printf("\nHydrogen abundance = %5.3f     Helium = %5.3f\n",ah,ahe);

    if(flagx == 1) infix(fixfile,atom,ah);

    if(flagw == 1) printf("\n");
    pfinit(V,atom,model,flagw);
    nline = 0;
    nlist = 0;
    if(flagw == 1) printf("Calculating Number Densities\n");
    if(model->teff >= 23500.0) veryhotDensity(model,atom,ah,ahe,flagw);
    else if(model->teff <= 8500.0) Density(model,atom,ah,ahe,flagw);
    else hotDensity(model,atom,ah,ahe,flagw);
    
    popinit(POP, atom, model, V, flagw);
    

    if(flagO == 1) {
        // Enter wavelength for opacity output
        opwave = 5000.0;
        opout = fopen("opacity.out","w");
        for(i=0;i<Ntau;i++) opacity(model,opwave,i);
        fclose(opout);
    }

    waveref = 5000.0;
    if(flagw == 1) printf("Calculating Reference Opacities\n");
    tauref(model,waveref);
    
    
    if(flagw == 1) printf("Entering Main Loop\n");

    int pos = 0;
    int dwave1 = dwave;
    int dwave2 = dwave;
    while(pos < num_measures) {
        wave = waveobs[pos];
        
        if ((pos > 0) && (pos < num_measures-1)) {
            dwave1 = waveobs[pos] - waveobs[pos-1];
            dwave2 = waveobs[pos+1] - waveobs[pos];
            if (dwave1 < dwave2) {
                dwave = dwave1;
            } else {
                dwave = dwave2;
            }
            inc = dwave;
        } else if (pos == 0) {
            dwave = 0.02;
        } // else: pos == num_measures-1 => use the last dwave calculated value
        
        Depth = 1.0;

        tauwave(model,wave);
        
        Flux = flux(model,wave);
        
        if (pos == 0) {
            // Reset static vars such as last wave (argument 8)
            linelst(wave,list,&nlist,atom,model->teff,model->logg,qf,1,isotope,model,Flux,V,POP,dwave);
        } else {
            linelst(wave,list,&nlist,atom,model->teff,model->logg,qf,0,isotope,model,Flux,V,POP,dwave);
        }
        inlin(wave,line,&nline,list,nlist,atom,isotope);
        for(i=0;i<nline;i++) {
            if(line[i].flag == 0) {
                pop(line,i,model,V,POP);
                broad(model,line,i,line[i].sig,line[i].alp,line[i].fac);
                capnu(line,i,model);
            }
        }
        taukap(wave,model,atom,line,nline,strgln,V,He,POP);
        Depth = depth(model,wave,Flux);
        if(flagw == 1) printf("%9.3f %d %d\n",wave,nline,nlist);
        
        fluxes[pos] = 1.0 - Depth;
        
        if (pos % 100 == 0) {
            user_func(((1.0*pos)/num_measures)*100.0, user_data);
        }
    
        pos++;        
    }
    fclose(qf);
    
    free(model);
    free(V);
    free(atom);
    free(POP);
    free(line);
    free(oneline);
    free(strgln);
    free(list);
    free(He);
    free(velgrad);
    free(isotope);
    
    return(0);
}


