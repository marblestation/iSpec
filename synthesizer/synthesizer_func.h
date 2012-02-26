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
int Ntau;
float **bkap;
float **bkap2;
float **bkap3;
float **bkap4;
double inc;

int flagr;
int flagc;
int flagk;
int flagg;
int flagmgh;
int flagI;
int flagt;
int flagp;
int flagP;
int flagu;
int flagO;
int flagC;
int mghla;
int mghlb;
float *velgrad;
double mu;
int NI;
// variables for isotopes
double ra1H,ra2H,ra12C,ra13C,ra14N,ra15N,ra16O,ra17O,ra18O;
double ra24Mg,ra25Mg,ra26Mg,ra28Si,ra29Si,ra30Si,ra40Ca,ra42Ca;
double ra43Ca,ra44Ca,ra46Ca,ra48Ca,ra46Ti,ra47Ti,ra48Ti,ra49Ti;
double ra50Ti;
//memo reset;
FILE *opout;
//linedata *oneline;
memo reset;
FILE *opout;
linedata *oneline;

typedef void (*progressfunc)(double num, void *user_data);

int synthesize_spectrum(char *atmosphere_model_file, char *linelist_file, char *abundances_file, double microturbulence_vel, int verbose, int num_measures, const double waveobs[], double fluxes[], progressfunc user_func, void *user_data);

