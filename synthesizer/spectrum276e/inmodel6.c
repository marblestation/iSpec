#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include "spectrum.h"
double partfn(double code, double T, double Ne);

void inmodel(model, file, flagw)
atmosphere *model;
char *file;
int flagw;
{
    int i, N;
    FILE *fp;
    double teff, logg, MH;
    double k = 8.617084e-05;
    char buffer[150], *tmp;
    char buf2[150];
    extern int Ntau;
    extern int flagt;
    extern int flagu;

    /* if flagt == 1, inmodel parses Kurucz (ATLAS9) headers */

    if ((fp = fopen(file, "r")) == NULL) {
        printf("Cannot open atmosphere data file\n");
        exit(1);
    }

    if (flagt == 0) {
        /* This section deals with models with the traditional SPECTRUM simple 
           header but can deal with Castelli models with extra columns */
        fgets(buffer, 80, fp);
        model->teff = atof(strtok(buffer, " "));
        if (model->teff > 60000.0 || model->teff < 3000.0)
            printf("It looks like this model does not have the traditional SPECTRUM\nheader.  You may need to use the t switch, especially if you\nget a segmentation fault.\n");
        model->logg = atof(strtok(NULL, " "));
        model->MH = atof(strtok(NULL, " "));
        Ntau = atoi(strtok(NULL, " "));
        if (flagw == 1)
            printf("Teff = %6.0f log(g) = %5.2f [M/H] = %5.2f\n", model->teff, model->logg, model->MH);

        for (i = 0; i < Ntau; i++) {
            fgets(buffer, 120, fp);
            model->mass[i] = atof(strtok(buffer, " "));
            model->T[i] = atof(strtok(NULL, " "));
            model->kT[i] = k * model->T[i];
            model->P[i] = atof(strtok(NULL, " "));
            model->Ne[i] = atof(strtok(NULL, " "));
            model->U[i] = partfn(1.0, model->T[i], model->Ne[i]);
            tmp = strtok(NULL, " ");
            tmp = strtok(NULL, " ");
            if (flagu == 1)
                model->mtv[i] = atof(strtok(NULL, " "));
        }
    } else if (flagt == 1) {
        /* This section attempts to parse a Kurucz ATLAS9 header. Success not
           guaranteed! Can deal with Castelli models with extra columns */
        do {
            fgets(buffer, 120, fp);
        }
        while (strstr(buffer, "TEFF") == NULL);
        strcpy(buf2, strstr(buffer, "TEFF"));
        tmp = strtok(buf2, " ");
        model->teff = atof(strtok(NULL, " "));
        strcpy(buf2, strstr(buffer, "GRAVITY"));
        tmp = strtok(buf2, " ");
        model->logg = atof(strtok(NULL, " "));
        do {
            fgets(buffer, 120, fp);
        }
        while (strstr(buffer, "SCALE") == NULL);
        strcpy(buf2, strstr(buffer, "SCALE"));
        tmp = strtok(buf2, " ");
        model->MH = log10(atof(strtok(NULL, " ")));
        if (flagw == 1)
            printf("Teff = %6.0f log(g) = %5.2f [M/H] = %5.2f\n", model->teff, model->logg, model->MH);
        do {
            fgets(buffer, 120, fp);
        }
        while (strstr(buffer, "READ") == NULL);
        i = 0;
        while (fgets(buffer, 120, fp) != NULL) {
            if (strstr(buffer, "PRADK") != NULL)
                break;
            model->mass[i] = atof(strtok(buffer, " "));
            model->T[i] = atof(strtok(NULL, " "));
            model->kT[i] = k * model->T[i];
            model->P[i] = atof(strtok(NULL, " "));
            model->Ne[i] = atof(strtok(NULL, " "));
            model->U[i] = partfn(1.0, model->T[i], model->Ne[i]);
            tmp = strtok(NULL, " ");
            tmp = strtok(NULL, " ");
            if (flagu == 1)
                model->mtv[i] = atof(strtok(NULL, " "));
            i++;
        }
        Ntau = i;
    } else {
        /* This should never happen! */
        printf("\nModel atmosphere format ambiguous, exiting\n");
        exit(1);
    }

    fclose(fp);
}
