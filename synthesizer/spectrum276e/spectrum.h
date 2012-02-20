#define NTAU 100
#define NATOM 107
#define NMOL 15
#define TNATOM 308
#define NLINE 5000
#define NLIST 50000
#define NHE 5
#define NSTRG 19

typedef struct {
  double wave;
  double code;
  int iso;
  double atomass;
  double abund;
  double chi1;
  double chi2;
  double chi3;
  double chi4;
  double chi;
  double Eu;
  double El;
  double gf;
  double wavel;
  double waveh;
  float  xnum[NTAU];
  float  a[NTAU];
  float  dopp[NTAU];
  float  capnu[NTAU];
  float  dlg[NTAU];
  char   T[5];
  double alp;
  double sig;
  double gammar;
  double gammas;
  double gammaw;
  double fac;
  int ai;
  int flag;
} linedata;

typedef struct {
  int code;
  double amass;
  double abund;
  double I1;
  double I2;
  double I3;
  double I4;
  int maxcharge;
} atominfo;

typedef struct {
  double atom;
  float Q4[NTAU];
  float R1[NTAU];
  float R2[NTAU];
  float R3[NTAU];
} population;

typedef struct {
  double teff;
  double logg;
  double MH;
  double mass[NTAU];
  double tauref[NTAU];
  double x[NTAU];
  double tauwave[NTAU];
  double tauwstart[NTAU];
  double tauwend[NTAU];
  double T[NTAU];
  float  kT[NTAU];
  float  U[NTAU];
  float  P[NTAU];
  float  rho[NTAU];
  float  mtv[NTAU];
  double kapparef[NTAU];
  double kappawave[NTAU];
  double kappawstart[NTAU];
  double kappawend[NTAU];
  double Ne[NTAU];
  double NA[NTAU];
  double NHI[NTAU];
  double Np[NTAU];
  double NH2[NTAU];
  double N1[NTAU];
  double N2[NTAU];
  double N3[NTAU];
  double NHminus[NTAU];
  double NHeI[NTAU];
  double NHeII[NTAU];
  double NHeIII[NTAU];
  double NCI[NTAU];
  double NNI[NTAU];
  double NOI[NTAU];
  double NCH[NTAU];
  double NOH[NTAU];
  double NMgH[NTAU];
  double NH2O[NTAU];
  float  NMgI[NTAU];
  float  NAlI[NTAU];
  float  NSiI[NTAU];
  float  NMgII[NTAU];
  float  NCaII[NTAU];
  float  NSiII[NTAU];
  float  NFeI[NTAU];
  float  NCaI[NTAU];
  float  NTiI[NTAU];
  float  NZrI[NTAU];
  double kapnu[NTAU];
  double taunu[NTAU];
  double stim[NTAU];
  int nmax[NTAU];
} atmosphere;

typedef struct {
  double lambda[4];
  double kappa[4][NTAU];
  int flag;
} inkappa;

typedef struct {
  double wave;
  double code;
  int iso;
  double El;
  double Eu;
  double loggf;
  double wavel;
  double waveh;
  double fac;
  double alp;
  double sig;
  double gammar;
  double gammas;
  double gammaw;
  char T[5];
  int flag;
  int ai;
} linelist;

typedef struct {
  double linecenter,end;
  double w[NTAU],dw[NTAU],a[NTAU],sig[NTAU],kap[NTAU];
  int flag;
} Helium;

typedef struct {
  double species[TNATOM];
  float pf[TNATOM][NTAU];
} pfunc;

typedef struct {
  int lyman;
  int balmer;
  int paschen;
  int brackett;
  int pfund;
  int humphreys;
  int hprofl;
  int helium;
  int strong;
  int interval;
} memo;

typedef struct {
  double code;
  int iso;
  double atomass;
  double relabund;
} isodata;
