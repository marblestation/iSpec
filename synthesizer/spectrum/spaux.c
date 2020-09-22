#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXIT 100
#define EULER 0.5772156649015328606
#define FPMIN 1.0e-200
#define EPS 1.0e-07
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

double dmax(),dmin();
int imax(),imin();
double poly();
char *ggets(char *s);
double expint();
void gaussj();
void nrerror();
void newlocate();
void newhunt();
int imin();

double dmax(x,y)
double x,y;
{
  if(x > y) return(x);
  else return(y);
}

int imax(x,y)
int x,y;
{
  if(x > y) return(x);
  else return(y);
}

double dmin(x,y)
double x,y;
{
  if(x < y) return(x);
  else return(y);
}

int imin(x,y)
int x,y;
{
  if(x < y) return(x);
  else return(y);
}

double poly(x,n,coeffs)
double x,coeffs[];
int n;
{
  double result;
  int i;

  result = coeffs[n];

  for(i=n-1;i>=0;i--) {
    result *= x;
    result += coeffs[i];
  }

  return(result);
}

char *ggets(char *s)
{
  size_t n;
  char tmp[80];

  if(fgets(tmp,80,stdin) == NULL) {
    printf("Access error in ggets\n");
    exit(1);
  }
  n = strlen(tmp);
  strncpy(s,tmp,n-1);
  s[n-1] = '\0';
  return(s);
}

double expint(x,n)
     double x;
     int n;
{
  void nrerror();
  int i,ii,nm1;
  double a,b,c,d,del,fact,h,psi,ans;

  nm1 = n-1;
  if(x < 0.0) x = 0.0;
  if(n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
    nrerror("bad arguments in expint");
  else {
    if(n == 0) ans = exp(-x)/x;
    else {
      if(x == 0.0) ans = 1.0/nm1;
      else {
	if(x > 1.0) {
	  b = x+n;
	  c = 1.0/FPMIN;
	  d=1.0/b;
	  h=d;
	  for(i=1;i<=MAXIT;i++) {
	    a = -i*(nm1+i);
	    b += 2.0;
	    d = 1.0/(a*d+b);
	    c = b+a/c;
	    del = c*d;
	    h *= del;
	    if(fabs(del-1.0) < EPS) {
	      ans = h*exp(-x);
	      return ans;
	    }
	  }
      printf("expint: %f %i\n", x, n); // SBC
	  nrerror("continued fraction failed in expint");
	} else {
	  ans = (nm1 != 0 ? 1.0/nm1 : -log(x)-EULER);
	  fact = 1.0;
	  for(i=1;i<=MAXIT;i++) {
	    fact *= -x/i;
	    if(i != nm1) del = -fact/(i-nm1);
	    else {
	      psi = -EULER;
	      for(ii=1;ii<=nm1;ii++) psi += 1.0/ii;
	      del = fact*(-log(x)+psi);
	    }
	    ans += del;
	    if(fabs(del) < fabs(ans)*EPS) return ans;
	  }
	  nrerror("series failed in expint");
	}
      }
    }
  } 
  return ans;
}

void gaussj(a,n,b,m)
double **a,**b;
int n,m;
{
   int *indxc,*indxr,*ipiv;
   int i,icol,irow,j,k,l,ll,*ivector();
   double big,dum,pivinv;
   void nrerror(),free_ivector();

   indxc = ivector(1,n);
   indxr = ivector(1,n);
   ipiv = ivector(1,n);
   for(j=1;j<=n;j++) ipiv[j] = 0;
   for(i=1;i<=n;i++) {
      big = 0.0;
      for(j=1;j<=n;j++)
	 if(ipiv[j] != 1)
	    for(k=1;k<=n;k++) {
	       if(ipiv[k] == 0) {
		 if(fabs(a[j][k]) >= big) {
		   big = fabs(a[j][k]);
		   irow = j;
		   icol = k;
		 }
	       } else if(ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
	    }
      ++(ipiv[icol]);
      if(irow != icol) {
	for(l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
	for(l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
      }
      indxr[i] = irow;
      indxc[i] = icol;
      if(a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
      pivinv = 1.0/a[icol][icol];
      a[icol][icol] = 1.0;
      for(l=1;l<=n;l++) a[icol][l] *= pivinv;
      for(l=1;l<=m;l++) b[icol][l] *= pivinv;
      for(ll=1;ll<=n;ll++)
	 if(ll != icol) {
	   dum = a[ll][icol];
	   a[ll][icol] = 0.0;
	   for(l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	   for(l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
	 }
   }
   for(l=n;l>=1;l--) {
      if(indxr[l] != indxc[l])
	for(k=1;k<=n;k++)
	   SWAP(a[k][indxr[l]],a[k][indxc[l]]);
   }
   free_ivector(ipiv,1,n);
   free_ivector(indxr,1,n);
   free_ivector(indxc,1,n);
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
   int i;
   float **m;
   m = (float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i] = (float*) malloc((unsigned) (nch-ncl+1)*sizeof(float));
      if(!m[i]) nrerror("allocation failure 2 in matrix()");
      m[i] -= ncl;
   }
   return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
   int i;
   double **m;

   m = (double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if(!m) nrerror("allocation failure 1 in dmatrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
       m[i] = (double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
       if(!m[i]) nrerror("allocation failure 2 in dmatrix()");
       m[i] -= ncl;
   }
   return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
   int i,j;
   float **m;

   m = (float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
   if (!m) nrerror("allocation failure in submatrix()");
   m -= newrl;

   for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j] = a[i]+oldcl-newcl;
   return m;
}

void free_submatrix(b,nrl,nrh,ncl,nch)
int **b;
int nrl,nrh,ncl,nch;
{
   free((char*) (b+nrl));
}


void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

float **cmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
   int i;
   float **cm;
   if(nrl != 0 || ncl != 0) {
     nrerror("Matrices declared by cmatrix must be zero-offset");
     exit(1);
   }
   cm = (float **) calloc(nrh+1,sizeof(float*));
   if (!cm) nrerror("allocation failure 1 in cmatrix()");

   for(i=0;i<=nrh;i++) {
      cm[i] = (float*) calloc(nch+1,sizeof(float));
      if(!cm[i]) nrerror("allocation failure 2 in cmatrix()");
   }
   return cm;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
  int i,**m;

  m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
  if(!m) nrerror("Allocation failure 1 in imatrix");
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i] = (int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
    if (!m[i]) nrerror("Allocation failure 2 in imatrix");
    m[i] -= ncl;
  }
  return m;
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
  int i;
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}

void nrerror(error_text)
char error_text[];
{
   void exit();
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,". . . now exiting to system . . .\n");
   exit(1);
}

float *vector(nl,nh)
int nl,nh;
{
    float *v;

    v = (float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
    if(!v) nrerror("allocation failure in vector()");
    return v-nl;
}

void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
    free((char*) (v+nl));
}


double *dvector(nl,nh)
int nl,nh;
{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if(!v) nrerror("allocation failure in dvector()");
    return v-nl;
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
    free((char*) (v+nl));
}

int *ivector(nl,nh)
int nl,nh;
{
    int *v;

    v = (int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if(!v) nrerror("allocation failure in ivector()");
    return v-nl;
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
    free((char*) (v+nl));
}

double *cvector(nl,nh)
int nl,nh;
{
    double *v;

    v = (double *)calloc((nh-nl+1),sizeof(double));
    if(!v) nrerror("allocation failure in cvector()");
    return v-nl;
}

int approx(x1,x2,prec)
double x1,x2,prec;
{
  double diff;

  diff = fabs(x1 - x2);
  if(diff > -1.0*prec && diff < prec) return(1);
  else return(0);
}


void polynn(xa,ya,n,x,y)
float xa[],ya[],x,y[];
int n;
{
  float polyn = 0.0;
  int npoly = 4;
  int nm,j,k,l,nm1,nup,lll,m;
  float term;

  nm = (npoly+1)/2;
  nup = n + nm - npoly + 1;

  newlocate(xa,n,x,&j,nm,nup);

  l = (int) imin(j - nm,n-npoly);
  lll = l + npoly - 1;

  for(k=l;k<=lll;k++) {
    term = 1.0;
    for(m=l;m<=lll;m++) {
      if(k == m) continue;
      term *= (x - xa[m])/(xa[k] - xa[m]);
    }
    term *= ya[k];
    polyn += term;
  }
  *y = polyn;
}


void polynn2(x1a,x2a,ya,m,n,x1,x2,y)
float x1a[],x2a[],**ya,x1,x2,*y;
int m,n;
{
  int k,j,jm,l,lll,nm,nup;
  int npoly = 4;
  float *ymtmp,*vector();
  void polynn(),free_vector();

  nm = (npoly+1)/2;
  nup = m + nm - npoly + 1;
  ymtmp = vector(0,m-1);

  newlocate(x1a,m,x1,&jm,nm,nup);
  l = (int)imin(jm-nm,m-npoly);
  lll = l + npoly - 1;

  for(j=0;j<m;j++) {
    if(j>=l && j <=lll) polynn(x2a,ya[j],n,x2,&ymtmp[j]);
    else ymtmp[j] = 0.0;
  }
  polynn(x1a,ymtmp,m,x1,y);
  free_vector(ymtmp,0,m-1);
}

void polynn3(x0a,x1a,x2a,ya,l,m,n,x0,x1,x2,y)
float x0a[],x1a[],x2a[],***ya,x0,x1,x2,*y;
int l,m,n;
{
  int k,j,jl,l1,lll,nm,nup;
  float *ymtmp,*vector();
  void polynn2(),polynn(),free_vector();
  int npoly = 4;

  nm = (npoly+1)/2;
  nup = l + nm - npoly + 1;

  ymtmp = vector(0,l-1);

  newlocate(x0a,l,x0,&jl,nm,nup);
  l1 = (int)imin(jl-nm,l-npoly);
  lll = l1 + npoly - 1;

  for(j=0;j<l;j++) {
    if(j >= l1 && j <= lll) polynn2(x1a,x2a,ya[j],m,n,x1,x2,&ymtmp[j]);
    else ymtmp[j] = 0.0;
  }
  polynn(x0a,ymtmp,l,x0,y);
  free_vector(ymtmp,0,l-1);
}

void newlocate(xx,n,x,j,nlow,nhigh)
float xx[],x;
int n,*j,nlow,nhigh;
/* Given an array xx[0..n-1], and a given value of x, returns a value j */
/* such that x is between xx[j-1] and xx[j] and nlow <= j <= nhigh  */
{
   int ascnd,ju,jm,jl;

   jl = nlow-1;
   ju = nhigh + 1;
   ascnd = xx[nhigh] > xx[0];
   while(ju-jl > 1) {
      jm = (ju + jl) >> 1;
      if(x > xx[jm] == ascnd)
	 jl=jm;
      else
	 ju=jm;
   }
   *j=jl+1;
}

void dlocate(xx,n,x,j,nlow,nhigh)
double xx[],x;
int n,*j,nlow,nhigh;
/* Given an array xx[0..n-1], and a given value of x, returns a value j */
/* such that x is between xx[j-1] and xx[j] and nlow <= j <= nhigh  */
{
   int ascnd,ju,jm,jl;

   jl = nlow-1;
   ju = nhigh + 1;
   ascnd = xx[nhigh] > xx[0];
   while(ju-jl > 1) {
      jm = (ju + jl) >> 1;
      if(x > xx[jm] == ascnd)
	 jl=jm;
      else
	 ju=jm;
   }
   *j=jl+1;
}

/* Given xx[0...n-1] and x, returns a value jlo such that x lies between */
/* xx[jlo], xx[jlo+1]   */



void newhunt(xx,n,x,jlo,nlow,nhigh)
int n,*jlo,nlow,nhigh;
float xx[],x;
{
   int jm,jhi,inc,ascnd;

   ascnd = (xx[nhigh] > xx[0]);
   if(*jlo < 0 || *jlo > nhigh) {
      *jlo = nlow-1;
      jhi = nhigh+1;
   } else {
      inc = 1;
      if(x >= xx[*jlo] == ascnd) {
	 if(*jlo == nhigh) return;
	 jhi = (*jlo) + 1;
	 while(x >= xx[jhi] == ascnd) {
	    *jlo = jhi;
	    inc += inc;
	    jhi = (*jlo) + inc;
	    if(jhi > nhigh) {
	       jhi = nhigh+1;
	       break;
	    }
	 }
      } else {
	 if(*jlo == 0) {
	    return;
	 }
	 jhi = (*jlo);
	 *jlo -= 1;
	 while(x < xx[*jlo] == ascnd) {
	    jhi = (*jlo);
	    inc += inc;
	    *jlo = jhi - inc;
	    if(*jlo < 0) {
	       *jlo = -1;
	       break;
	    }
	 }
      }
   }
   while(jhi-(*jlo) != 1) {
      jm = (jhi+(*jlo)) >> 1;
      if(x > xx[jm] == ascnd)
	 *jlo = jm;
      else
	 jhi = jm;
   }
   *jlo += 1;
}


void opolynn(xa,ya,n,x,y,h)
float xa[],ya[],x,y[];
int n,*h;
{
  float polyn = 0.0;
  int npoly = 4;
  int nm,j,k,l,nm1,nup,lll,m;
  float term;

  nm = (npoly+1)/2;
  nup = n + nm - npoly + 1;

  newhunt(xa,n,x,h,nm,nup);
  j = *h;

  l = (int) imin(j - nm,n-npoly);
  lll = l + npoly - 1;

  for(k=l;k<=lll;k++) {
    term = 1.0;
    for(m=l;m<=lll;m++) {
      if(k == m) continue;
      term *= (x - xa[m])/(xa[k] - xa[m]);
    }
    term *= ya[k];
    polyn += term;
  }
  *y = polyn;
}


void opolynn2(x1a,x2a,ya,m,n,x1,x2,y,jh,kh)
float x1a[],x2a[],**ya,x1,x2,*y;
int m,n,*jh,*kh;
{
  int k,j,jm,l,lll,nm,nup;
  int npoly = 4;
  float *ymtmp,*vector();
  void polynn(),free_vector();

  nm = (npoly+1)/2;
  nup = m + nm - npoly + 1;
  ymtmp = vector(0,m-1);

  newhunt(x1a,m,x1,jh,nm,nup);
  jm = *jh;
  l = (int)imin(jm-nm,m-npoly);
  lll = l + npoly - 1;

  for(j=0;j<m;j++) {
    if(j>=l && j <=lll) opolynn(x2a,ya[j],n,x2,&ymtmp[j],kh);
    else ymtmp[j] = 0.0;
  }
  opolynn(x1a,ymtmp,m,x1,y,jh);
  free_vector(ymtmp,0,m-1);
}

void opolynn3(x0a,x1a,x2a,ya,l,m,n,x0,x1,x2,y,ih,jh,kh)
float x0a[],x1a[],x2a[],***ya,x0,x1,x2,*y;
int l,m,n,*ih,*jh,*kh;
{
  int k,j,jl,l1,lll,nm,nup;
  float *ymtmp,*vector();
  void polynn2(),polynn(),free_vector();
  int npoly = 4;

  nm = (npoly+1)/2;
  nup = l + nm - npoly + 1;

  ymtmp = vector(0,l-1);

  newhunt(x0a,l,x0,ih,nm,nup);
  jl = *ih;
  l1 = (int)imin(jl-nm,l-npoly);
  lll = l1 + npoly - 1;

  for(j=0;j<l;j++) {
    if(j >= l1 && j <= lll) opolynn2(x1a,x2a,ya[j],m,n,x1,x2,&ymtmp[j],jh,kh);
    else ymtmp[j] = 0.0;
  }
  opolynn(x0a,ymtmp,l,x0,y,ih);
  free_vector(ymtmp,0,l-1);
}


void search(float xx[], int n, float x, int *jlo)
{
  int jm,jhi,inc;
  int ascnd;

  ascnd=(xx[n] >= xx[1]);
  if (*jlo <= 0 || *jlo > n) {
     *jlo=0;
     jhi=n+1;
  } else {
     inc=1;
     if (x >= xx[*jlo] == ascnd) {
	if (*jlo == n) return;
	   jhi=(*jlo)+1;
	   while (x >= xx[jhi] == ascnd) {
	     *jlo=jhi;
	     inc += inc;
	     jhi=(*jlo)+inc;
	     if (jhi > n) {
		jhi=n+1;
		break;
	     }
	    }
      } else {
	 if (*jlo == 1) {
	    *jlo=0;
	    return;
	 }
	    jhi=(*jlo)--;
	    while (x < xx[*jlo] == ascnd) {
	      jhi=(*jlo);
	      inc <<= 1;
	      if (inc >= jhi) {
		 *jlo=0;
		 break;
	      }
	      else *jlo=jhi-inc;
	  }
       }
    }
    while (jhi-(*jlo) != 1) {
      jm=(jhi+(*jlo)) >> 1;
      if (x >= xx[jm] == ascnd)
	 *jlo=jm;
      else
	 jhi=jm;
      }
      if (x == xx[n]) *jlo=n-1;
      if (x == xx[1]) *jlo=1;
}
