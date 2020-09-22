#include <stdio.h>
#include "spectrum.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
int Ntau = 72;
void inmodel();
void Density();
void hotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw);
void veryhotDensity(atmosphere *model,atominfo *atom,double ah,double ahe,int flagw);
void tauref();
void tauwave();
void tauflx();
double flux();
void fluxflx();
double intensity();
void intenint();
void inlin();
void pop();
void broad(atmosphere *model,linedata *line,int N,double sig,
           double alp,double fac);
void capnu();
void taukap();
void inatom();
void linelst();
void pfinit();
void popinit();
double depth();
double depthflx();
double depthmu();
double idepthmu();
double opacity();
float **cmatrix();
void nrerror();
char *ggets(char *s);
void inisotope();
void isorelabun();
void invelgrad();
void infix(char fixfile[], atominfo *atom, double ah);
float **bkap;
float **bkap2;
float **bkap3;
float **bkap4;
double interval();
memo reset;
void setreset();
double dmax(),dmin();
int imax(),imin();
double inc;
int flagr = 0;
int flagc = 0;
int flagk = 0;
int flagg = 0;
int flagmgh = 0;
int flagI = 0;
int flagt = 0;
int flagp = 0;
int flagP = 0;
int flagu = 0;
int flagO = 0;
int flagC = 0;
/* flagt = 1 and SPECTRUM parses ATLAS9 headers in atmosphere models */
int mghla = 0;
int mghlb = 0;
float *velgrad;
double mu = 1.0;
int NI = 0;
/* variables for isotopes */
double ra1H,ra2H,ra12C,ra13C,ra14N,ra15N,ra16O,ra17O,ra18O;
double ra24Mg,ra25Mg,ra26Mg,ra28Si,ra29Si,ra30Si,ra40Ca,ra42Ca;
double ra43Ca,ra44Ca,ra46Ca,ra48Ca,ra46Ti,ra47Ti,ra48Ti,ra49Ti;
double ra50Ti;
FILE *opout;
linedata *oneline;

int main(int argc, char *argv[])
{
  int i,k;
  int nq;
  int nline = 0;
  int nlist = 0;
  int code;
  int flagf = 0;
  int flagm = 0;
  int flagi = 1;
  int flagb = 0;
  int flagw = 1;
  int flaga = 0;
  int flagd = 0;
  int flagM = 0;
  int flagG = 0;
  int flagN = 0;
  int flagC = 0;
  int flagfM = 0;
  int flagx = 0;
  int ninit = 0;
  int nbytes;
  double ah,ahe,waveref,start,end,wave,dwave,dw,Flux,Depth,wstart,wend;
  double vturb = 2.0e+5;
  double nabund,DW=20.0,Intensity,flxstart,flxend,q,y,intenstart,intenend;
  float qd;
  unsigned len_written;
  atmosphere *model;
  atominfo *atom;
  linelist *list;
  linedata *line;
  linedata *strgln;
  isodata *isotope;
  pfunc *V;
  population *POP;
  Helium *He;
  char *file,name[80],oname[80],lines[80],c,tmp[20],atmdat[80];
  char fixfile[80];
  char *vg,vgrad[80],isofile[80];
  char *ofile;
  FILE *qf,*fp;
  int fd;
  char junk[10];
  double opwave,kap;

  setreset(0);
  strcpy(atmdat,"stdatom.dat");

  if(argc > 1) {
    if(strcspn(argv[1],"f") != strlen(argv[1])) flagf = 1;
    if(strcspn(argv[1],"m") != strlen(argv[1])) flagm = 1;
    if(strcspn(argv[1],"M") != strlen(argv[1])) flagM = 1;
    if(strcspn(argv[1],"b") != strlen(argv[1])) flagb = 1;
    if(strcspn(argv[1],"n") != strlen(argv[1])) flagw = 0;
    if(strcspn(argv[1],"a") != strlen(argv[1])) flaga = 1;
    if(strcspn(argv[1],"d") != strlen(argv[1])) flagd = 1;
    if(strcspn(argv[1],"r") != strlen(argv[1])) flagr = 1;
    if(strcspn(argv[1],"c") != strlen(argv[1])) flagc = 1;
    if(strcspn(argv[1],"k") != strlen(argv[1])) flagk = 1;
    if(strcspn(argv[1],"g") != strlen(argv[1])) flagg = 1;
    if(strcspn(argv[1],"t") != strlen(argv[1])) flagt = 1;
    if(strcspn(argv[1],"G") != strlen(argv[1])) flagG = 1;
    if(strcspn(argv[1],"i") != strlen(argv[1])) flagI = 1;
    if(strcspn(argv[1],"N") != strlen(argv[1])) flagN = 1;
    if(strcspn(argv[1],"p") != strlen(argv[1])) flagp = 1;
    if(strcspn(argv[1],"P") != strlen(argv[1])) flagP = 1;
    if(strcspn(argv[1],"u") != strlen(argv[1])) flagu = 1;
    if(strcspn(argv[1],"O") != strlen(argv[1])) flagO = 1;
    if(strcspn(argv[1],"x") != strlen(argv[1])) flagx = 1;
  }
  if(flagG == 1) flagg = 1;
  if(flagP == 1) flagp = 1;
  if(flagN == 1) flagw = 0;
  if(flagm == 1 || flagM == 1) flagf = 0;
  if(flagm == 1 || flagf == 1 || flagM == 1) flagi = 0;
  if(flagf == 1 || flagM == 1) flagfM = 1;
  else flagfM = 0;
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
    printf("\n(C) Richard O. Gray 1992 - 2018 Version 2.77");
    printf("\nFeb 25, 2018");
    if(flagi == 1) printf("\nIntegrated Disk mode (normalized Intensity)\n\n");
    if(flagf == 1) printf("\nAbsolute Flux Mode (Integrated Disk)\n\n");
    if(flagM == 1) printf("\nSpecific Intensity Mode\n\n");
    if(flagm == 1) printf("\nSpecific Intensity Mode (normalized Intensity)\n\n");
    if(flagb == 1) printf("Output file will be binary\n");
    if(flagc == 1) printf("Output will be continuum only (no line absorption)\n");
    if(flagw == 0) printf("Silent mode\n");
    if(flagg == 1) printf("Velocity gradient mode\n");
    if(flagI == 1) printf("Isotopes enabled\n");
    if(flagu == 1) printf("Reading microturbulent velocity from atmosphere model\n");
  }

  if(flagw == 1) printf("\nEnter name of stellar atmosphere data file > ");
  file = ggets(name);
  inmodel(model,file,flagw);
  if(flagG == 1) {
    printf("\nEnter name of velocity gradient file > ");
    ggets(vgrad);
  } else if(flagg == 1) strcpy(vgrad,"velgrad.dat");
  if(flagg == 1) invelgrad(vgrad);
  if(flagg == 0) for(i=0;i<Ntau;i++) velgrad[i] = 0.0;
  if(flagw == 1) 
    printf("\nEnter name of line list file: (default = luke.lst) > ");
  ggets(lines);
  if(strcmp(lines,"") == 0) strcpy(lines,"luke.lst");
  if((qf = fopen(lines,"r")) == NULL) {
   printf("Cannot find line data file\n");
   exit(1);
  }
  if(flaga == 1) {
    if(flagw == 1) 
      printf("\nEnter name of atom data file (default = stdatom.dat) > ");
    ggets(atmdat);
    if(strcmp(atmdat,"") == 0) strcpy(atmdat,"stdatom.dat");
  }
  if(flagI == 1) {
    if(flagw == 1) 
       printf("\nEnter name of isotope data file (default = isotope.iso) > ");
    ggets(isofile);
    if(strcmp(isofile,"") == 0) strcpy(isofile,"isotope.iso");
    inisotope(isofile,isotope);
    isorelabun(isotope);
  }
  if(flagx == 1) {
    printf("\nEnter name of fixed abundance file > ");
    ggets(fixfile);
  }
  if(flagw == 1) printf("\nEnter name of output file > ");
  ofile = ggets(oname);
  if(flagw == 1 && flagu == 0) printf("\nEnter microturbulence (km/s) > ");
  if(flagu == 0) {
    nq = scanf("%lf",&vturb);
    vturb *= 1.0e+05;
    for(i=0;i<Ntau;i++) model->mtv[i] = vturb;
  }

  if(flagm == 1 || flagM == 1) {
    if(flagw == 1) printf("\nEnter mu = cos(theta) > ");
    nq = scanf("%lf",&mu);
  }
  if(flagw == 1) printf("\nEnter beginning and ending wavelengths (A)> ");
  nq = scanf("%lf,%lf",&start,&end);
  if(flagw == 1) printf("\nEnter wavelength step (A)> ");
  nq = scanf("%lf",&dwave);
  inc = dwave;

  if(flagb == 0) fp = fopen(ofile,"w");
  if(flagb == 1) {
    fd = open(oname,O_CREAT|O_TRUNC|O_RDWR|S_IREAD|S_IWRITE,0666);
    if(fd == -1) {
     perror("Error:");
     exit(1);
    }
    nbytes = write(fd,&model->teff,sizeof(double));
    nbytes = write(fd,&model->logg,sizeof(double));
    nbytes = write(fd,&model->MH,sizeof(double));
    q = vturb/1.0e+05;
    nbytes = write(fd,&q,sizeof(double));
    nbytes = write(fd,&start,sizeof(double));
    nbytes = write(fd,&dwave,sizeof(double));
  }
  /* Normal abundances for hydrogen and helium; inatom may change these */
  ah = 0.911;
  ahe = 0.089;
  inatom(atmdat,atom,model->MH,&ah,&ahe);
  if(flaga == 1) printf("\nHydrogen abundance = %5.3f   Helium = %5.3f\n",ah,ahe);

  if(flagx == 1) infix(fixfile,atom,ah);

  if(flagw == 1) printf("\n");
  pfinit(V,atom,model,flagw);
  nline = 0;
  nlist = 0;
  if(flagw == 1) printf("Calculating Number Densities\n");
  if(model->teff >= 23500.0) veryhotDensity(model,atom,ah,ahe,flagw);
  else if(model->teff <= 8500.0) Density(model,atom,ah,ahe,flagw);
  else hotDensity(model,atom,ah,ahe,flagw);
  popinit(POP,atom,model,V,flagw);

  if(flagO == 1) {
    printf("Enter wavelength for opacity output > ");
    nq = scanf("%lf",&opwave);
    opout = fopen("opacity.out","w");
    for(i=0;i<Ntau;i++) kap = opacity(model,opwave,i);
    fclose(opout);
  }

  waveref = 5000.0;
  if(flagw == 1) printf("Calculating Reference Opacities\n");
  tauref(model,waveref);
  wstart = start;
  DW = interval(start,dwave,flagfM);
  wend = start + DW;
  wend = (double)dmin(wend,end);
  if(flagw == 1) printf("Entering Main Loop\n");

  while(wstart < end - 0.0001) {
    Intensity = Depth = 1.0;
     wave = (wstart + wend)/2.0;
     if(flagi == 1 || flagm == 1 || ninit == 0) tauwave(model,wave);
     if(flagf == 1 || flagM == 1) tauflx(model,wstart,wend);
     Flux = flux(model,wave);
     /* if(flagi == 1) Flux = flux(model,wave); */
     if(flagf == 1) fluxflx(model,wstart,wend,&flxstart,&flxend);
     if(flagm == 1) Intensity = intensity(model,wave,mu);
     if(flagM == 1) intenint(model,wstart,wend,&intenstart,&intenend,mu);

     wave = wstart;
     ninit = 1;

     while(wend-wave > dwave/2.0 + 0.0001) {
	linelst(wave,list,&nlist,atom,model->teff,model->logg,qf,0,isotope,
                model,Flux,V,POP,dwave);
	inlin(wave,line,&nline,list,nlist,atom,isotope);
	for(i=0;i<nline;i++) {
	  if(line[i].flag == 0) {
	  pop(line,i,model,V,POP);
	  broad(model,line,i,line[i].sig,line[i].alp,line[i].fac);
	  capnu(line,i,model);
	  }
	}
	taukap(wave,model,atom,line,nline,strgln,V,He,POP);
	if(flagi == 1) Depth = depth(model,wave,Flux);
	if(flagf == 1) Depth = depthflx(model,wave,wstart,wend);
	if(flagm == 1) Depth = depthmu(model,wave,Intensity,mu);
	if(flagM == 1) Depth = idepthmu(model,wave,wstart,wend,mu);
	if(flagw == 1) printf("%9.3f %d %d\n",wave,nline,nlist);
	if(flagb == 0 && (flagi == 1 || flagm == 1))
	       fprintf(fp,"%9.3f   %f\n",wave,1.0 - Depth);
	if(flagb == 1 && (flagi == 1 || flagm == 1)) {
	       qd = Depth;
	       nbytes = write(fd,&qd,sizeof(float));
	}
	if(flagf == 1) {
	  Flux = flxstart + (flxend-flxstart)*(wave-wstart)/(wend-wstart);
          if(Depth > Flux) Depth = Flux;
	  if(flagb == 0) {
	    if(Depth == 0.0 && Flux == 0.0) fprintf(fp,"%9.3f   %e\n",wave,0.0);
	    else fprintf(fp,"%9.3f   %e\n",wave,Flux*(1.0-Depth/Flux));
	  }
	  if(flagb == 1) {
	    if(Depth == 0.0 && Flux == 0.0) qd = 0.0;
	    else qd = Flux*(1.0-Depth/Flux);
	    nbytes = write(fd,&qd,sizeof(float));
	  }
	}
        if(flagM == 1) {
	  Intensity = intenstart + 
                      (intenend-intenstart)*(wave-wstart)/(wend-wstart);
          if(Depth > Intensity) Depth = Intensity;
	  if(flagb == 0) {
	    if(Depth == 0.0 && Intensity == 0.0) fprintf(fp,"%9.3f   %e\n",wave,0.0);
	    else fprintf(fp,"%9.3f   %e\n",wave,Intensity*(1.0-Depth/Intensity));
	  }
	  if(flagb == 1) {
	    if(Depth == 0.0 && Intensity == 0.0) qd = 0.0;
	    else qd = Intensity*(1.0-Depth/Intensity);
	    nbytes = write(fd,&qd,sizeof(float));
	  }
	}
	wave += dwave;
     }

     wstart = wave;
     DW = interval(wave,dwave,flagf);
     DW = (double)dmax(dwave,DW);
     wend = wstart + DW;
     wend = (double) dmin(wend,end);
  /* if(flagb == 0) {
       fclose(fp);
       if((fp = fopen(ofile,"a")) == NULL) {
	 printf("Cannot open %s for appending\n");
	 exit(1);
       }
     }  */
  }
  fclose(qf);
  if(flagb == 0) fclose(fp);
  if(flagb == 1) close(fd);
  return(0);
}



