/* 
 * File:   filesio.h
 * Author: sousasag
 *
 * Created on April 23, 2011, 2:21 PM
 */

#ifndef _FILESIO_H
#define	_FILESIO_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <string.h>

    /* leitura das opcoes , mine.opt file read*/
//void read_mine(char * mine, char* filetest, char* fileleitura, char* fileout, double* lambdai, double* lambdaf, double* smoothder, double* tree, double* space, double* tree, double* distlinha, double* miniline, int* plots_flag);

long file_lines(char * file);
char* getFileExtension(char * name);
int set_tipo_espectro(char * extensao);


void read_mine(char * mine, char* filetest, char* fileleitura, char* fileout, double* lambdai, double* lambdaf, double* smoothder, double* space, char* tree, double* distlinha, double* miniline, int* plots_flag, char* rvmask);
void read_spectrum_file(char * filetest, long* npoints, double** pixels2, double** xpixels2, double* cdelta1, double* crval1, int spectrum_type);
void read_fits_file(char * filetest, long* npoints, double** pixels2, double** xpixels2, double* cdelta1, double* crval1);
void read_ascii_file(char * filetest, long* npoints, double** pixels2, double** xpixels2, double* cdelta1mean, double* crval1);
void read_lines_list(char * fileleitura, long* nl2, double** linhas);
void write_outfile(char *fileout, double* aponta, int nl, int miniline);

void write_outfile(char *fileout, double* aponta, int nl, int miniline){
   
    
    FILE *pFile2;
    pFile2 = fopen (fileout,"wt");
    int ilinha;
    for (ilinha=0;ilinha<nl;ilinha++) {
        if (aponta[ilinha*9+4] > miniline && aponta[ilinha*9+4] < 500. && aponta[ilinha*9+1] > 0)
            fprintf(pFile2," %10.2f  %ld  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.2f\n", aponta[ilinha*9+0], (long) aponta[ilinha*9+1], aponta[ilinha*9+2], aponta[ilinha*9+3], aponta[ilinha*9+4], aponta[ilinha*9+8], aponta[ilinha*9+5], aponta[ilinha*9+6], aponta[ilinha*9+7]);

    }
    fclose (pFile2);
}

long file_lines(char * file) {
    FILE * pFile;
    long n = 0;
    char str[200];
    pFile = fopen (file,"rt");
    if (pFile==NULL) perror ("Error opening file");
    else
    {
    while (!feof(pFile)) {
            fgets (str , 200, pFile);
            n++;
            }
    fclose (pFile);
    }
    return (n);
}



void read_lines_list(char * fileleitura, long* nl2, double** linhas){
    long nl;
    char str[200], * pch;
    nl=file_lines(fileleitura);
    double* lab=(double *) malloc((nl-2) * sizeof(double));
    if (lab == NULL) {
            printf("Memory allocation error\n");
            exit(1);
    }

    FILE * pFile = fopen (fileleitura,"rt");
    fgets (str , 200, pFile);
    fgets (str , 200, pFile);
    int il;
    for (il=0; il < nl-3; il++) {
            fgets (str , 200, pFile);
            pch = strtok (str," ");
            char *stopstring;
            lab[il]=strtod(pch,&stopstring);
    }
    fclose (pFile);
    *linhas=lab;
    *nl2=nl-2;
}

int set_tipo_espectro(char * extensao){
    if (!strcmp(extensao, ".fits") || !strcmp(extensao, ".FITS"))
        return 0; // FITS file input
    else return 1; // ASCII file input
}


char* getFileExtension(char * name){
    printf("FILE: %s\n",name);
    char * extension= (char *) malloc(10*sizeof(char));
    char* peek = name + strlen (name) - 1;
    while (peek >= name) {
        if (*peek == '.') {
            strcpy (extension, peek);
            break;
        }
        peek--;
    }
    return extension;
}



void read_spectrum_file(char * filetest, long* npoints, double** pixels2, double** xpixels2, double* cdelta1, double* crval1, int spectrum_type){
    if (spectrum_type == 0){
        read_fits_file(filetest, npoints, pixels2, xpixels2, cdelta1, crval1);
    } else {
        read_ascii_file(filetest, npoints, pixels2, xpixels2, cdelta1, crval1);
        printf("\n ASCII READ\n");
    }

}


void read_ascii_file(char * filetest, long* npoints, double** pixels2, double** xpixels2, double* cdelta1mean, double* crval1){    

    long nl;
    char str[200], * pch, *stopstring;
    double *pixels, *xpixels;
    nl=file_lines(filetest)-1;
    

    pixels = (double *) malloc((nl-2) * sizeof(double));
    xpixels= (double *) malloc((nl-2) * sizeof(double));
    if (pixels == NULL || xpixels == NULL) {
            printf("Memory allocation error\n");
            exit(1);
    }
    *pixels2=pixels;
    *xpixels2=xpixels;
    FILE * pFile = fopen (filetest,"rt");

    double meanx=0.;
    double xtmp;
    double ftmp;

    fgets (str , 200, pFile);
    fgets (str , 200, pFile);
    long il=0;
            fgets (str , 200, pFile);
            sscanf (str,"%lf %lf",&xtmp,&ftmp);
            xpixels[il]=xtmp;
            pixels[il]=ftmp;
    for (il=1; il < nl-2; il++) {
            fgets (str , 200, pFile);
//            printf("%ld: %s\n",il,str);
            sscanf (str,"%lf %lf",&xtmp,&ftmp);
            xpixels[il]=xtmp;
            pixels[il]=ftmp;
            meanx+=xpixels[il]-xpixels[il-1];
    }
    fclose (pFile);
    *crval1=xpixels[0];
    *npoints=nl-2;
    *cdelta1mean=meanx/(nl-2);
    
}

void read_fits_file(char * filetest, long* npoints, double** pixels2, double** xpixels2, double* cdelta1, double* crval1){
    printf("READING FITS FILE: %s \n", filetest);
    long naxes[2] = {1,1}, fpixel[2] = {1,1};
    int bitpix, naxis, ii;
    char format[20], hdformat[20];
    double *pixels, *xpixels;
    double caga;
    int status = 0;/* CFITSIO status value MUST be initialized to zero! */
    fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    if (!fits_open_file(&fptr, filetest, READONLY, &status)) {
        printf("Ficheiro Aberto...");
            if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) ) {
                if (naxis > 2 || naxis == 0)
                    printf("Error: only 1D or 2D images are supported\n");
                else {
                    /* get memory for 1 row */
                    pixels = (double *) malloc(naxes[0] * sizeof(double));
                    xpixels = (double *) malloc(naxes[0] * sizeof(double));
                    *pixels2=pixels;
                    *xpixels2=xpixels;
                    if (pixels == NULL || xpixels == NULL) {
                            printf("Memory allocation error\n");
                            exit(1);
                    }

                    if (bitpix > 0) {  /* set the default output format string */
                           strcpy(hdformat, " %7d");
                           strcpy(format,   " %7.0f");
                    } else {
                        strcpy(hdformat, " %15d");
                        strcpy(format,   " %15.5f");
                    }

                    printf("Lendo ficheiro...\n");          /* print column header */
                    status=fits_read_key(fptr, TDOUBLE, "CDELT1", &caga, card, &status);
                    *cdelta1=caga;
                    status=fits_read_key(fptr, TDOUBLE, "CRVAL1", &caga, card, &status);
                    *crval1=caga;
                    *npoints=naxes[0];

                    for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--) {
                        if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status) )  /* read row of pixels */
                            break;  /* jump out of loop on error */
                        for (ii = 0; ii < naxes[0]; ii++)
                            xpixels[ii]= (double)(ii * (*cdelta1) + (*crval1));
                    }
                }
            }
          fits_close_file(fptr, &status);
    }
    if (status) fits_report_error(stderr, status); /* print any error message */


}

void read_mine(char * mine, char* filetest, char* fileleitura, char* fileout, double* lambdai, double* lambdaf, double* smoothder, double* space, char* tree, double* distlinha, double* miniline, int* plots_flag, char* rvmask) {

    
    long nl;
    
    nl = file_lines("mine.opt")-1;
    printf("number of options: %ld\n", nl);
    
    FILE * pFile3;
    pFile3= fopen ("logARES.txt","wt");
    fprintf(pFile3,"Input File ARES...\n");


    
    FILE * fopt = fopen("mine.opt","rt");
    char str[200];
    system("clear");
    printf("Input Parameters:\n\n");
    fprintf(pFile3,"Input Parameters:\n\n");

    fgets (str , 200 , fopt);
    char *pch;
    pch = strtok (str,"' ");
    pch = strtok (NULL, "' ");
    strcpy (filetest,pch);
    printf("specfits: %s\n",filetest);
    fprintf(pFile3,"specfits: %s\n",filetest);


    fgets (str , 200 , fopt);
    pch = strtok (str,"' ");
    pch = strtok (NULL, "' ");
    strcpy (fileleitura,pch);
    printf("readlinedat: %s\n",fileleitura);
    fprintf(pFile3,"readlinedat: %s\n",fileleitura);


    fgets (str , 200 , fopt);
    pch = strtok (str,"' ");
    pch = strtok (NULL, "' ");
    strcpy (fileout,pch);
    printf("fileout: %s\n",fileout);
    fprintf(pFile3,"fileout: %s\n",fileout);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *lambdai=atof(pch);
    printf("lambdai: %6.1f\n",*lambdai);
    fprintf(pFile3,"lambdai: %6.1f\n",*lambdai);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *lambdaf=atof(pch);
    printf("lambdaf: %6.1f\n",*lambdaf);
    fprintf(pFile3,"lambdaf: %6.1f\n",*lambdaf);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *smoothder=atof(pch);
    printf("smoothder: %3.1f\n",*smoothder);
    fprintf(pFile3,"smoothder: %3.1f\n",*smoothder);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *space=atol(pch);
    printf("space: %.2f \n",(float) *space);
    fprintf(pFile3,"space: %.2f \n",(float) *space);


//    fgets (str , 200 , fopt);
//    pch = strtok (str,"=");
//    pch = strtok (NULL, " ");
//    *tree=atof(pch);
//    printf("tree: %5.3f\n",*tree);
//    fprintf(pFile3,"tree: %5.3f\n",*tree);
    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    strcpy (tree,pch);
//    printf("tree: %s\n",tree);
//    fprintf(pFile3,"tree: %s\n",tree);
    printf("rejt: %s\n",tree);
    fprintf(pFile3,"rejt: %s\n",tree);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *distlinha=atof(pch);
    printf("lineresol: %4.3f\n",*distlinha);
    fprintf(pFile3,"lineresol: %4.3f\n",*distlinha);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *miniline=atof(pch);
    printf("miniline: %3.1f\n",*miniline);
    fprintf(pFile3,"miniline: %3.1f\n",*miniline);


    fgets (str , 200 , fopt);
    pch = strtok (str,"=");
    pch = strtok (NULL, " ");
    *plots_flag=atoi(pch);
    printf("plots_flag: %i\n",*plots_flag);
    fprintf(pFile3,"plots_flag: %i\n",*plots_flag);
    
    if (nl > 11){
        fgets (str , 200 , fopt);
        pch = strtok (str,"' ");
        pch = strtok (NULL, "' ");
        strcpy (rvmask,pch);
        printf("rvmask: %s\n",rvmask);
        fprintf(pFile3,"rvmask: %s\n",rvmask);
    } else {
        strcpy (rvmask,"0,0");
    }


    fclose (fopt);
    fclose (pFile3);


}

#ifdef	__cplusplus
}
#endif

#endif	/* _FILESIO_H */

