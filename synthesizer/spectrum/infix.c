#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "spectrum.h"

void infix(char fixfile[], atominfo *atom, double ah)
{
  char key[20];
  int i,k,flag=0;
  int code;
  FILE *fx;
  double labund,ltot;

  if((fx = fopen(fixfile,"r")) == NULL) {
    printf("Cannot open fixed abundance file %s\n",fixfile);
    exit(1);
  }

  if(fgets(key,20,fx) == NULL) {
    printf("File access error in infix\n");
    exit(1);
  }
  if(strcmp(key,"HYDROGEN\n") == 0) flag = 1;
  else if(strcmp(key,"TOTAL\n") == 0) flag = 2;
  else {
    printf("Keyword in fixed abundance file not recognized\n");
    exit(1);
  }

  if(flag == 1) {
    while(fscanf(fx,"%d %lf",&code,&labund) != EOF) {
      ltot = labund - 12.0 + log10(ah);
      for(i=0;i<NATOM;i++) {
	if(atom[i].code == code) {
	  k = i;
	  break;
	}
      }
      if(i >= NATOM) {
	printf("\nError 1 in INFIX\n");
	exit(1);
      }
      atom[k].abund = pow(10.0,ltot);
    }
  } else if(flag == 2) {
    while(fscanf(fx,"%d %lf",&code,&ltot) != EOF) {
      for(i=0;i<NATOM;i++) {
	if(atom[i].code == code) {
	  k = i;
	  break;
	}
      }
      if(i >= NATOM) {
	printf("\nError 1 in INFIX\n");
	exit(1);
      }
      atom[k].abund = pow(10.0,ltot);
    }
  } else {
    printf("\nError 2 in INFIX\n");
    exit(1);
  }
  return;
} 
