#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "spectrum.h"
double Hfnm();
double Vcse1f();
double Sofbet();
double Wave0();
double corr();
double air();
double wair();
double dmax(),dmin();
int imin(),imax();
double partfn(double code,double T,double Ne);

/* Fine Structure Version of Hprofl.  Includes approximate treatment of */
/* fine structure in the doppler cores of the hydrogen lines.  The exact  */
/* pattern is used for the alpha lines, and the m = infinity pattern for */
/* all other lines.  Function translated into C from the original */
/* subroutine by Deane Peterson */

/* This has been further modified to include the Lyman alpha quasi H2 
   satellites.  Translated from new version of Synthe by Robert Kurucz */ 

/* The new version in Synthe uses a normalization different from the
original (which involves multiplying the profile by the factor
1.77245*dop.  This translation preserves the original normalization */

float propbm[7][15][5] = 
  {{{-.980,-.967,-.948,-.918,-.873},
    {-.968,-.949,-.921,-.879,-.821},
    {-.950,-.922,-.883,-.830,-.764},
    {-.922,-.881,-.830,-.770,-.706},
    {-.877,-.823,-.763,-.706,-.660},
    {-.806,-.741,-.682,-.640,-.625},
    {-.691,-.628,-.588,-.577,-.599},
    {-.511,-.482,-.484,-.514,-.568},
    {-.265,-.318,-.382,-.455,-.531},
    {-.013,-.167,-.292,-.394,-.478},
    { .166,-.056,-.216,-.332,-.415},
    { .251, .035,-.122,-.237,-.320},
    { .221, .059,-.068,-.168,-.247},
    { .160, .055,-.037,-.118,-.189},
    { .110, .043,-.022,-.085,-.147}},
   {{-.242, .060, .379, .671, .894},
    { .022, .314, .569, .746, .818},
    { .273, .473, .605, .651, .607},
    { .432, .484, .489, .442, .343},
    { .434, .366, .294, .204, .091},
    { .304, .184, .079,-.025,-.135},
    { .167, .035,-.082,-.189,-.290},
    { .085,-.061,-.183,-.287,-.374},
    { .032,-.127,-.249,-.344,-.418},
    {-.024,-.167,-.275,-.357,-.420},
    {-.061,-.170,-.257,-.327,-.384},
    {-.047,-.124,-.192,-.252,-.306},
    {-.043,-.092,-.142,-.190,-.238},
    {-.038,-.070,-.107,-.146,-.187},
    {-.030,-.049,-.075,-.106,-.140}},
   {{-.484,-.336,-.206,-.111,-.058},
    {-.364,-.264,-.192,-.154,-.144},
    {-.299,-.268,-.250,-.244,-.246},
    {-.319,-.333,-.337,-.336,-.337},
    {-.397,-.414,-.415,-.413,-.420},
    {-.456,-.455,-.451,-.456,-.478},
    {-.446,-.441,-.446,-.469,-.512},
    {-.358,-.381,-.415,-.463,-.522},
    {-.214,-.288,-.360,-.432,-.503},
    {-.063,-.196,-.304,-.394,-.468},
    { .063,-.108,-.237,-.334,-.409},
    { .151,-.019,-.148,-.245,-.319},
    { .149, .016,-.091,-.177,-.246},
    { .115, .023,-.056,-.126,-.189},
    { .078, .021,-.036,-.091,-.145}},
   {{-.082, .163, .417, .649, .829},
    { .096, .316, .515, .660, .729},
    { .242, .393, .505, .556, .534},
    { .320, .373, .394, .369, .290},
    { .308, .274, .226, .152, .048},
    { .232, .141, .052,-.046,-.154},
    { .148, .020,-.094,-.200,-.299},
    { .083,-.070,-.195,-.299,-.385},
    { .031,-.130,-.253,-.348,-.422},
    {-.023,-.167,-.276,-.359,-.423},
    {-.053,-.165,-.254,-.326,-.384},
    {-.038,-.119,-.190,-.251,-.306},
    {-.034,-.088,-.140,-.190,-.239},
    {-.032,-.066,-.103,-.144,-.186},
    {-.027,-.048,-.075,-.106,-.142}},
   {{-.819,-.759,-.689,-.612,-.529},
    {-.770,-.707,-.638,-.567,-.498},
    {-.721,-.659,-.595,-.537,-.488},
    {-.671,-.617,-.566,-.524,-.497},
    {-.622,-.582,-.547,-.523,-.516},
    {-.570,-.545,-.526,-.521,-.537},
    {-.503,-.495,-.496,-.514,-.551},
    {-.397,-.418,-.448,-.492,-.547},
    {-.246,-.315,-.384,-.453,-.522},
    {-.080,-.210,-.316,-.406,-.481},
    { .068,-.107,-.239,-.340,-.418},
    { .177,-.006,-.143,-.246,-.324},
    { .184, .035,-.082,-.174,-.249},
    { .146, .042,-.046,-.123,-.190},
    { .103, .036,-.027,-.088,-.146}},
   {{-.073, .169, .415, .636, .809},
    { .102, .311, .499, .639, .710},
    { .232, .372, .479, .531, .514},
    { .294, .349, .374, .354, .279},
    { .278, .253, .212, .142, .040},
    { .215, .130, .044,-.051,-.158},
    { .141, .015,-.097,-.202,-.300},
    { .080,-.072,-.196,-.299,-.385},
    { .029,-.130,-.252,-.347,-.421},
    {-.022,-.166,-.275,-.359,-.423},
    {-.050,-.164,-.253,-.325,-.384},
    {-.035,-.118,-.189,-.252,-.306},
    {-.032,-.087,-.139,-.190,-.240},
    {-.029,-.064,-.102,-.143,-.185},
    {-.025,-.046,-.074,-.106,-.142}},
   {{ .005, .128, .260, .389, .504},
    { .004, .109, .220, .318, .389},
    {-.007, .079, .162, .222, .244},
    {-.018, .041, .089, .106, .080},
    {-.026,-.003, .003,-.023,-.086},
    {-.025,-.048,-.087,-.148,-.234},
    {-.008,-.085,-.165,-.251,-.343},
    { .018,-.111,-.223,-.321,-.407},
    { .032,-.130,-.255,-.354,-.431},
    { .014,-.148,-.269,-.359,-.427},
    {-.005,-.140,-.243,-.323,-.386},
    { .005,-.095,-.178,-.248,-.307},
    {-.002,-.068,-.129,-.187,-.241},
    {-.007,-.049,-.094,-.139,-.186},
    {-.010,-.036,-.067,-.103,-.143}}};


double Hprofl(wave,n,m,j,model)
atmosphere *model;
int n,m,j;
double wave;
{
  
  static double FO[NTAU],pp[NTAU],y1b[NTAU],y1s[NTAU],t3nhe[NTAU];
  static double t3nh2[NTAU];
  static double dopph[NTAU],c1d[NTAU],c2d[NTAU],gcon1[NTAU],gcon2[NTAU];
  static double xknmtb[3][4] = {{0.0001716,0.009019,0.1001,0.5820},
			       {0.0005235,0.01772,0.171,0.866},
			       {0.0008912,0.02507,0.223,1.02}};
  static double y1wtm[2][2] = {{1.0e+18,1.0e+17},{1.0e+16,1.0e+14}};
  /* Fine Structure Components */
  static double stalph[34] = {-730.,370.,188.,515.,327.,619.,-772.,-473.,
  -369.,120.,256.,162.,285.,-161.,-38.3,6.82,-174.,-147.,-101.,-77.5,
  55.,126.,75.,139.,-60.,3.7,27.,-69.,-42.,-18.,-5.5,-9.1,-33.,-24.};
  /* Alpha component Weights */
  static double stwtal[34] = {1.,2.,1.,2.,1.,2.,1.,2.,3.,1.,2.,1.,2.,1.,4.,
  6.,1.,2.,3.,4.,1.,2.,1.,2.,1.,4.,6.,1.,7.,6.,4.,4.,4.,5.};
  static int istal[4] = {1,3,10,21};
  static int lnghal[4] = {2,7,11,14};
  /* Fine structure for m = infinity */
  static double stcomp[4][5] = {{0.,0.,0.,0.,0.},{468.,576.,-522.,0.,0.,},
  {260.,290.,-33.,-140.,0.0},{140.,150.,18.,-27.,-51.}};
  static double stcpwt[4][5] = {{1.,0.,0.,0.,0.},{1.,1.,2.,0.,0.},
  {1.,1.,4.,3.,0.},{1.,1.,4.,6.,4.}};
  static int lncomp[4] = {1,3,4,5};
  static double cutoffh2plus[111] = 
    {-15.14,-15.06,-14.97,-14.88,-14.80,-14.71,-14.62,-14.53,
     -14.44,-14.36,-14.27,-14.18,-14.09,-14.01,-13.92,-13.83,
     -13.74,-13.65,-13.57,-13.48,-13.39,-13.30,-13.21,-13.13,
     -13.04,-12.95,-12.86,-12.77,-12.69,-12.60,-12.51,-12.40,
     -12.29,-12.15,-12.02,-11.90,-11.76,-11.63,-11.53,-11.41,
     -11.30,-11.22,-11.15,-11.09,-11.07,-11.06,-11.07,-11.09,
     -11.12,-11.16,-11.19,-11.21,-11.24,-11.27,-11.30,-11.33,
     -11.36,-11.39,-11.42,-11.45,-11.48,-11.48,-11.48,-11.48,
     -11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,
     -11.48,-11.48,-11.48,-11.48,-11.41,-11.40,-11.39,-11.38,
     -11.37,-11.36,-11.35,-11.34,-11.33,-11.32,-11.30,-11.29,
     -11.28,-11.27,-11.27,-11.27,-11.26,-11.25,-11.24,-11.23,
     -11.22,-11.21,-11.20,-11.19,-11.18,-11.17,-11.15,-11.14,
     -11.13,-11.12,-11.11,-11.10,-11.09,-11.08,-11.07};
  static double cutoffh2[91] = 
    {-13.64,-13.52,-13.39,-13.27,-13.14,-13.01,-12.87,-12.74,
     -12.63,-12.56,-12.51,-12.48,-12.47,-12.49,-12.52,-12.55,
     -12.57,-12.61,-12.65,-12.69,-12.72,-12.76,-12.79,-12.82,
     -12.84,-12.85,-12.87,-12.90,-12.93,-12.94,-12.93,-12.95,
     -12.95,-12.96,-12.97,-12.96,-12.96,-12.95,-12.95,-12.96,
     -12.98,-12.99,-12.95,-12.96,-13.00,-13.00,-12.98,-12.97,
     -13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,
     -13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,
     -13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-12.89,-12.88,
     -12.87,-12.86,-12.85,-12.84,-12.83,-12.81,-12.80,-12.79,
     -12.78,-12.76,-12.74,-12.72,-12.70,-12.68,-12.65,-12.62,
     -12.59,-12.56,-12.53};
  static double asumlyman[100] = 
    {0.000E+00, 6.265E+08, 1.897E+08, 8.126E+07, 4.203E+07, 2.450E+07,
     1.236E+07, 8.249E+06, 5.782E+06, 4.208E+06, 3.158E+06, 2.430E+06,
     1.910E+06, 1.567E+06, 1.274E+06, 1.050E+06, 8.752E+05, 7.373E+05,
     6.269E+05, 5.375E+05, 4.643E+05, 4.038E+05, 3.534E+05, 3.111E+05,
     2.752E+05, 2.447E+05, 2.185E+05, 1.959E+05, 1.763E+05, 1.593E+05,
     1.443E+05, 1.312E+05, 1.197E+05, 1.094E+05, 1.003E+05, 9.216E+04,
     8.489E+04, 7.836E+04, 7.249E+04, 6.719E+04, 6.239E+04, 5.804E+04,
     5.408E+04, 5.048E+04, 4.719E+04, 4.418E+04, 4.142E+04, 3.888E+04,
     3.655E+04, 3.440E+04, 3.242E+04, 3.058E+04, 2.888E+04, 2.731E+04,
     2.585E+04, 2.449E+04, 2.322E+04, 2.204E+04, 2.094E+04, 1.991E+04,
     1.894E+04, 1.804E+04, 1.720E+04, 1.640E+04, 1.566E+04, 1.496E+04,
     1.430E+04, 1.368E+04, 1.309E+04, 1.254E+04, 1.201E+04, 1.152E+04,
     1.105E+04, 1.061E+04, 1.019E+04, 9.796E+03, 9.419E+03, 9.061E+03,
     8.721E+03, 8.398E+03, 8.091E+03, 7.799E+03, 7.520E+03, 7.255E+03,
     7.002E+03, 6.760E+03, 6.530E+03, 6.310E+03, 6.100E+03, 5.898E+03,
     5.706E+03, 5.522E+03, 5.346E+03, 5.177E+03, 5.015E+03, 4.860E+03,
     4.711E+03, 4.569E+03, 4.432E+03, 4.300E+03};
  static double asum[100] =
    {0.000E+00, 4.696E+08, 9.980E+07, 3.017E+07, 1.155E+07, 5.189E+06,
     2.616E+06, 1.437E+06, 8.444E+05, 5.234E+05, 3.389E+05, 2.275E+05,
     1.575E+05, 1.120E+05, 8.142E+04, 6.040E+04, 4.560E+04, 3.496E+04,
     2.719E+04, 2.141E+04, 1.711E+04, 1.377E+04, 1.119E+04, 9.166E+03,
     7.572E+03, 6.341E+03, 5.338E+03, 4.523E+03, 3.854E+03, 3.302E+03,
     2.844E+03, 2.460E+03, 2.138E+03, 1.866E+03, 1.635E+03, 1.438E+03,
     1.269E+03, 1.124E+03, 9.983E+02, 8.894E+02, 7.947E+02, 7.120E+02,
     6.396E+02, 5.759E+02, 5.198E+02, 4.703E+02, 4.263E+02, 3.873E+02,
     3.526E+02, 3.215E+02, 2.938E+02, 2.689E+02, 2.465E+02, 2.264E+02,
     2.082E+02, 1.918E+02, 1.769E+02, 1.634E+02, 1.512E+02, 1.400E+02,
     1.298E+02, 1.206E+02, 1.121E+02, 1.043E+02, 9.720E+01, 9.066E+01,
     8.465E+01, 7.912E+01, 7.403E+01, 6.933E+01, 6.498E+01, 6.097E+01,
     5.725E+01, 5.381E+01, 5.061E+01, 4.765E+01, 4.489E+01, 4.232E+01,
     3.994E+01, 3.771E+01, 3.563E+01, 3.369E+01, 3.188E+01, 3.019E+01,
     2.860E+01, 2.712E+01, 2.572E+01, 2.442E+01, 2.319E+01, 2.204E+01,
     2.096E+01, 1.994E+01, 1.898E+01, 1.808E+01, 1.722E+01, 1.642E+01,
     1.566E+01, 1.495E+01, 1.427E+01, 1.363E+01};
  static int flag1 = 0;
  static int n1 = 0;
  static int m1 = 0;
/*  static float rydh = 3.288052195e+15; */
/* The Rydberg constant below gives consistency between the Rydberg
   formula and the correction to vacuum equation; even though it is wrong,
   it gives excellent wavelengths for Lyman and Balmer series on that
   basis */
  static double rydh = 3.2880928e+15;
  static double fnm,y1num,y1wht,freqnm,dbeta,wavenm,c1con,c2con,radamp;
  static double resont,vdw,hwvdw,hwrad,stark,wave0;
  int i,k,nwid,ifcore;
  static int mmn,ifins;
  static double finest[14],finswt[14];
  double t4,t43,xn,xn2,xm,xm2,xmn2,xm2mn2,gnm,xknm;
  double freq,xne16,del,delw,wl;
  double hwstk,hwlor,hfwid,hprof4,dop,d,ff,hhw,wty1,y1scal,c1,c2,g1,gnot;
  double beta,y1,y2,gam,prqs,f,p1,fns,hwdop,hwres,hprofres,cutoff,spacing;
  double freq22000,cutfreq,hprofrad,hprofvdw,hproflor;
  double top,beta4000,prqsp4000,cutoff4000,freq15000;
  int ipos,icut;
  double K = 1.38054e-16;
  extern int Ntau;
  extern memo reset;

  if(n <= 3) rydh = 3.2880928e+15;
  else rydh = 3.288065e+15;

  /* Convert wavelength to vacuum wavelength */
  wave = wair(wave);
  if(flag1 == 0 || reset.hprofl == 1) {
    /* Set up Depth Vectors */
    flag1 = 1;
    reset.hprofl = 0;
    for(k=0;k<Ntau;k++) {
      xne16 = pow(model->Ne[k],0.1666667);
      pp[k] = 0.08989*xne16/sqrt(model->T[k]);
      FO[k] = 1.25e-09*pow(xne16,4.0);
      y1b[k] = 2.0/(1.0 + 0.012/model->T[k]*sqrt(model->Ne[k]/model->T[k]));
      t4 =  model->T[k]/10000.0;
      t43 = pow(t4,0.3);
      y1s[k] = t43/xne16;
      t3nhe[k] = t43*model->NHeI[k];
      t3nh2[k] = t43*model->NH2[k];
      dopph[k] = sqrt(2.0*K*model->T[k]/1.008/1.660e-24 + 
		      pow(model->mtv[k],2.0))/2.99792458e+10;
      c1d[k] = 78940.0*FO[k]/model->T[k];
      c2d[k] = FO[k]*FO[k]/5.96e-23/model->Ne[k];
      gcon1[k] = 0.2 + 0.09*sqrt(t4)/(1.0 + model->Ne[k]/1.0e+13);
      gcon2[k] = 0.2/(1.0 + model->Ne[k]/1.0e+15);
    }
  }
  /* Set Up for this Line */
  if(n != n1 || m != m1) {
    n1 = n;
    m1 = m;
    fnm = 0.0265384*Hfnm(n,m);
    mmn = m-n;
    xn = (double)n;
    xn2 = xn*xn;
    xm = (double)m;
    xm2 = xm*xm;
    xmn2 = xm2*xn2;
    xm2mn2 = xm2-xn2;
    gnm = xm2mn2/xmn2;
    if(mmn <= 3 && n <= 4) xknm = xknmtb[mmn-1][n-1];
    if(mmn > 3 || n > 4) xknm = 5.5e-05/gnm*xmn2/(1.0 + 0.13/(double)mmn);
    y1num = 320.0;
    if(m == 2) y1num = 550.0;
    if(m == 3) y1num = 380.0;
    y1wht = 1.0e+13;
    if(mmn <= 3) y1wht = 1.0e+14;
    if(mmn <= 2 && n <= 2) y1wht = y1wtm[mmn-1][n-1];
    freqnm = rydh*gnm;
    /*    wave0 = Wave0(n,m); */
    dbeta = 2.997924562e+18/(freqnm*freqnm)/xknm;
    wavenm = 2.997924562e+18/freqnm;
    c1con = xknm/wavenm*gnm*xm2mn2;
    c2con = pow(xknm/wavenm,2.0);
    radamp = 1.389e+09/pow(xm,4.53)/(1.0 + 5.0/xm2/xm);
    if(n != 1) radamp = radamp+1.389e+09/pow(xn,4.53)/(1.0 + 5.0/xn2/xn);
    /* radamp = radamp/freqnm; */
    /*   radamp = asum[n-1]+asum[m-1];
	 if(n == 1) radamp = asumlyman[m-1]; */
    /* radamp /= 12.5664; */
    radamp /= freqnm;
    resont = Hfnm(1,m)/xm/(1.0 - 1.0/xm2); 
    if(n != 1) resont = resont + Hfnm(1,n)/xn/(1.0 - 1.0/xn2);
    resont = resont*3.579e-24/gnm;  /* modified in new hprofl */
    vdw = 4.45e-26/gnm*pow((xm2*(7.0*xm2+5.0)),0.4);
    hwvdw = vdw*t3nhe[j]+2.0*vdw*t3nh2[j];
    hwrad = radamp;
    stark = 1.6678e-18*freqnm*xknm;
/* Fine Structure Components */
    if(n > 4 || m > 10) {
      ifins = 1;
      finest[0] = 0.0;
      finswt[0] = 1.0;
    } else {
      if(mmn != 1) {
	ifins = lncomp[n-1];
	for(i=0;i<ifins;i++) {
	  finest[i] = stcomp[n-1][i]*1.0e+07;
	  finswt[i] = stcpwt[n-1][i]/xn2;
	}
      } else {
	/* alpha lines */
	ifins = lnghal[n-1];
	ipos = istal[n-1];
	for(i=0;i<ifins;i++) {
	  k = ipos-1+i;
	  finest[i] = stalph[k]*1.0e+07;
	  finswt[i] = stwtal[k]/xn2/3.0;
	}
      }
    }
  }
/* Now Do This Depth */
  delw = wave - wavenm;
  /* delw = wave - wave0; */
  wl = wavenm+delw;
  freq = 2.99792458e+18/wl;
  del = fabs(freq-freqnm);
/* Convert wl to nm */
  wl /= 10.0;
/* These half-widths are really dnu/nu */
  hwstk = stark*FO[j];
  hwvdw = vdw*t3nhe[j]+2.0*vdw*t3nh2[j];
  hwrad = radamp;
  hwres = resont*model->N1[j];
  hwlor = hwres + hwvdw + hwrad;
  hwdop = dopph[j];
  /* Specify the largest half width in case of core calculation */
  /* nwid =1, doppler; =2, lorentz; =3, stark  */
  if(hwdop >= hwstk && hwdop >= hwlor) nwid = 1;
  else if(hwlor >= hwstk) nwid = 2;
  else nwid = 3;

  hfwid = freqnm*(double)dmax((double)dmax(hwdop,hwlor),hwstk);
  /* Sets flag if in a line core */
  hprof4 = 0.0;
  ifcore = 0;
  if(fabs(del) <= hfwid) ifcore = 1;
  /* Do Doppler */
  dop = freqnm*hwdop;
  /* Next is beginning of computed goto */
  if(ifcore == 0 || (ifcore == 1 && nwid == 1 )) {
/* Put fine structure in Doppler core */
    for(i=0;i<ifins;i++) {
      d = fabs(freq-freqnm-finest[i])/dop;
      if(d <= 7.0) hprof4 = hprof4 + exp(-d*d)/1.77245/dop*finswt[i];
    }
    if(ifcore == 1) return(hprof4*fnm);
  }
  /* -------------------------------------------------------------- */
  /* Do Lorentz */
  if(ifcore == 0 || (ifcore == 1 && nwid == 2 )) {
    if(n == 1 && m == 2) {   
      /* Lyman Alpha; modify old resonance broadening to match at 4000 cm-1 */
      hwres = hwres*4.0;
      hwlor = hwres+hwvdw+hwrad;
      hhw = freqnm*hwlor;
      hprofres = 0.0;
      if(freq > (82259.105-4000.)*2.99792458e+10) 
        hprofres = hwres*freqnm/3.14159/(del*del+hhw*hhw);
      else {
	/* only in far wing; data from Allard 1997 */
	cutoff = 0.0;
	if(freq < 50000.0*2.99792458e+10) hprofres = 0.0;
	else {
	  spacing = 200.0*2.99792458e+10;
	  freq22000 = (82259.105-22000.)*2.99792458e+10;
	  if(freq < freq22000) {
	    cutoff = (cutoffh2[1]-cutoffh2[0])/spacing*(freq-freq22000)
                     +cutoffh2[0];
	  } else {
	    icut = (freq-freq22000)/spacing;
	    cutfreq = icut*spacing+freq22000;
	    cutoff = (cutoffh2[icut+1]-cutoffh2[icut])/spacing*
                     (freq-cutfreq) + cutoffh2[icut];
	  }
	  cutoff = pow(10.0,cutoff-14.0)*model->N1[j]/2.99792458e+10;
	}
	hprofres = cutoff;
      }
      hprofrad = 0.0;
      /* Rayleigh scattering */
      if(freq >= 2.4190611e+15 && freq < 0.77*3.288051e+15) 
	hprofrad = hwrad*freqnm/3.14159/(del*del + hhw*hhw);
      /* Note: above calculation of hprofrad gives some problems in 
         Lyman alpha core near 10000K.  Seems okay at 20000K*/
      hprofvdw = hwvdw*freqnm/3.14159/(del*del+hhw*hhw);
	/* hprofrad = hwrad*freqnm/3.14159/(del*del + hhw*hhw)*1.77245*dop;
	   hprofvdw = hwvdw*freqnm/3.14159/(del*del+hhw*hhw)*1.77245*dop; */
      if(freq < 1.8e+15) hprofvdw=0.0;
      hproflor = hprofres+hprofrad+hprofvdw;
      hprof4 = hprof4+hproflor;
      if(ifcore == 1) return(hprof4*fnm);
    } else {
    /* Not Lyman alpha */
      hhw = freqnm*hwlor;
      top = hhw;
      if(n == 1 && m <= 5) {
        /* Lyman beta */
        if(m == 3 && (freq > 0.885*3.288051e+15 && freq < 0.890*3.288051e+15)) 
	  top = hhw - freqnm*hwrad;
        /* Lyman gamma */
        if(m == 4 && (freq > 0.936*3.288051e+15 && freq < 0.938*3.288051e+15))
	  top = hhw - freqnm*hwrad;
        /* Lyman delta */
        if(m == 5 && (freq > 0.959*3.288051e+15 && freq < 0.961*3.288051e+15))
	  top = hhw - freqnm*hwrad;
      }
    /*    hproflor = top/3.14159/(del*del+hhw*hhw)*1.77245*dop;  */
      hproflor = top/3.14159/(del*del+hhw*hhw);
      hprof4=hprof4+hproflor;
      if(ifcore == 1) return(hprof4*fnm);
    }
  }

  /* Do Stark */
  if(ifcore == 0 || (ifcore == 1 && nwid == 3)) {
    wty1 = 1.0/(1.0+model->Ne[j]/y1wht);
    y1scal = y1num*y1s[j]*wty1 + y1b[j]*(1.0 - wty1);
    c1 = c1d[j]*c1con*y1scal;
    c2 = c2d[j]*c2con;
    g1 = 6.77*sqrt(c1);
    gnot = g1*(double)dmax(0.0,0.2114+log(sqrt(c2)/c1))*
						 (1.0-gcon1[j]-gcon2[j]);
    beta = fabs(del)/FO[j]*dbeta;
    y1 = c1*beta;
    y2 = c2*beta*beta;
    gam = gnot;
    if(y2 > 0.0001 || y1 > 0.00001) {
      gam = g1*(0.5*exp(-(double)dmin(80.0,y1))+Vcse1f(y1) - 0.5*Vcse1f(y2))*
	      (1.0-gcon1[j]/(1.0+pow(90.0*y1,3.0))-gcon2[j]/(1.0+2000.0*y1));
      if(gam <= 1.0e-20) gam = 0.0;
    }
    prqs = Sofbet(beta,pp[j],n,m);
    if(m <= 2) {
      prqs = 0.5*prqs;
      cutoff = 0.0;
      /*Lyman Alpha quasi H2+ cutoff.  Data from Allard 1997 */
      if(freq >= (82259.105-20000.0)*2.99792458e+10) { 
	 if(freq <= (82259.105-4000.0)*2.99792458e+10) {
	   freq15000 = (82259.105-15000.0)*2.99792458e+10;
	   spacing = 100.0*2.99792458e+10;
	   if(freq < freq15000) cutoff = cutoffh2plus[0] + (freq-freq15000)*
                  (cutoffh2plus[1]-cutoffh2plus[0])/spacing;
	   else {
             icut = (freq-freq15000)/spacing;
             cutfreq = icut*spacing+freq15000;
	     cutoff = (freq-cutfreq)*(cutoffh2plus[icut+1]-cutoffh2plus[icut])/
	       spacing + cutoffh2plus[icut];
	   }
	   cutoff = model->Np[j]*pow(10.0,cutoff-14.0)/2.99792458e+10;
           hprof4=hprof4+cutoff;
         } else {
	   beta4000 = 4000.0*2.99792458e+10/FO[j]*dbeta;
	   prqsp4000 = Sofbet(beta4000,pp[j],n,m)*0.5*dbeta/FO[j]; 
	   cutoff4000 = pow(10.0,-11.07-14.0)*model->Np[j]/2.99792458e+10;
           hprof4 = hprof4 + cutoff4000/prqsp4000*prqs/FO[j]*dbeta;
         }
      }
    }
    f = 0.0;
    if(gam > 0.0) f = gam/3.14159/(gam*gam + beta*beta);
    p1 = pow(0.9*y1,2.0);
    fns = (p1+0.03*sqrt(y1))/(p1+1.0);
    hprof4 = hprof4 + (prqs*(1.0+fns)+f)/FO[j]*dbeta;
    return(hprof4*fnm);
  }
}


double Hfnm(n,m)
int n,m;
{
  double xn,ginf,gca,fkn,wtc,xm,xmn,fk,xmn12,wt,fnm;

  if(m <= n) return(0.0);
  xn = (double)n;
  ginf = 0.2027/pow(xn,0.71);
  gca = 0.124/xn;
  fkn = xn*1.9603;
  wtc = 0.45 - (2.4/pow(xn,3.0))*(xn-1.0);
  xm = (double)m;
  xmn = (double)(m-n);
  fk = fkn*pow((xm/(xmn*(xm+xn))),3.0);
  xmn12 = pow(xmn,1.2);
  wt = (xmn12 - 1.0)/(xmn12+wtc);
  fnm = fk*(1.0-wt*ginf-(0.222+gca/xm)*(1.0-wt));
  return(fnm);
}

double Vcse1f(x)
double x;
{
  double vcse1f;

  vcse1f = 0.0;
  if(x <= 0.0) return(vcse1f);
  if(x <= 0.01) {
    vcse1f = -log(x) - 0.577215 + x;
    return(vcse1f);
  }
  if(x <= 1.0) {
    vcse1f = -log(x) - 0.57721566+x*(0.99999193+x*(-0.24991055+x*(0.05519968+
	     x*(-0.00976004+x*0.00107857))));
    return(vcse1f);
  }
  if(x <= 30.0) {
    vcse1f = (x*(x+2.334733)+0.25062)/(x*(x+3.330657)+1.681534)/x*exp(-x);
    return(vcse1f);
  }
  return(vcse1f);
}



double Sofbet(b,p,n,m)
double b,p;
int n,m;
{
  double corr,b2,sb,wtpp,wtpm,wtbp,wtbm,cbp,cbm,pr1,pr2,wt,sofbet,cc,dd;
  int indx,mmn,im,ip,jm,jp,j;

  static float c[7][5] = {{-18.396, 84.674,-96.273,  3.927, 55.191},
			 { 95.740, 18.489, 14.902, 24.466, 42.456},
			 {-25.088,145.882,-50.165,  7.902, 51.003},
			 { 93.783, 10.066,  9.224, 20.685, 40.136},
			 {-19.819, 94.981,-79.606,  3.159, 52.106},
			 {111.107, 11.910,  9.857, 21.371, 41.006},
			 {511.318,  1.532,  4.044, 19.266, 41.812}};
  static float d[7][5] = {{ 11.801,  9.079, -0.651,-11.071,-26.545},
			 { -6.665, -7.136,-10.605,-15.882,-23.632},
			 {  7.872,  5.592, -2.716,-12.180,-25.661},
			 { -5.918, -6.501,-10.130,-15.588,-23.570},
			 { 10.938,  8.028, -1.267,-11.375,-26.047},
			 { -5.899, -6.381,-10.044,-15.574,-23.644},
			 { -6.070, -4.528, -8.759,-14.984,-23.956}};
  static float pp[] = {0.0,0.2,0.4,0.6,0.8};
  static float beta[] = {1.000,1.259,1.585,1.995,2.512,3.162,3.981,5.012,
			6.310,7.943,10.0,12.59,15.85,19.95,25.12};

  corr = 1.0;
  b2 = b*b;
  sb = sqrt(b);
  if(b <= 500) {
    indx = 7;
    mmn = m-n;
    if(n <= 3 && mmn <= 2) indx = 2*(n-1) + mmn;
    /* Determine Relevant Debye Range */
    im = (int)imin((int)(5.0*p)+1,4);
    ip = im+1;
    wtpp = 5.0*(p - pp[im-1]);
    wtpm = 1.0 - wtpp;
    if(b <= 25.12) {
      for(j=1;j<15;j++) {
	if(b <= beta[j]) break;
      }
      jm = j-1;
      jp = j;
      wtbp = (b - beta[jm])/(beta[jp] - beta[jm]);
      wtbm = 1.0 - wtbp;
      cbp = propbm[indx-1][jp][ip-1]*wtpp +
	    propbm[indx-1][jp][im-1]*wtpm;
      cbm = propbm[indx-1][jm][ip-1]*wtpp +
	    propbm[indx-1][jm][im-1]*wtpm;
      corr = 1.0 + cbp*wtbp + cbm*wtbm;
      /* Get inner approximate profile */
      pr1 = 0.0;
      pr2 = 0.0;
      wt = (double)dmax((double)dmin(0.5*(10-b),1.0),0.0);
      if(b <= 10) pr1 = 8.0/(83.0 + (2.0 + 0.95*b2)*b);
      if(b >= 8) pr2 = (1.5/sb + 27.0/b2)/b2;
      sofbet = (pr1*wt + pr2*(1.0-wt))*corr;
      return(sofbet);
    }
    /* Asymptotic Parts */
    cc = c[indx-1][ip-1]*wtpp + c[indx-1][im-1]*wtpm;
    dd = d[indx-1][ip-1]*wtpp + d[indx-1][im-1]*wtpm;
    corr = 1.0 + dd/(cc + b*sb);
  }
  sofbet = (1.5/sb + 27.0/b2)/b2*corr;
  return(sofbet);
}

double corr(freq)
double freq;
{
  return(3.04798e-04 - 5.50637798e+09/freq);
}

double air(freq)
double freq;
{
  double n,wn,wn2,wnair;
  double c = 2.997924562e+10;

  wn = freq/c;
  wn2 = wn*wn;
  n = 1.0 + 6432.8e-08 + 2949810.0/(146.0e+08 - wn2) +
      25540.0/(41.0e+08 - wn2);
  wnair = wn*n;
  return(wnair*c);
}

double wair(wave)
double wave;
{
  double n,wn,wn2,wvac,freq;
  double cA = 2.997924562e+18;
  double ccm = 2.997924562e+10;

  if(wave < 2800.0) return(wave);
  /* find approximate wavenumber and index of refraction */

  freq = cA/wave;
  wn = freq/ccm;
  wn2 = wn*wn;
  n = 1.0 + 6432.8e-08 + 2949810.0/(146.0e+08 - wn2) +
      25540.0/(41.0e+08 - wn2);
  wvac = wave + wave*(n-1.0);
  /* Do it again to get better approximation */
  freq = cA/wvac;
  wn = freq/ccm;
  wn2 = wn*wn;
  n = 1.0 + 6432.8e-08 + 2949810.0/(146.0e+08 - wn2) +
      25540.0/(41.0e+08 - wn2);
  wvac = wave + wave*(n-1.0);
  return(wvac);
}

double Wave0(n,m)
int n,m;
{
  double R = 3.288052195e+15;
  double N,M;
  double freqnm;

  N = (double)n;
  M = (double)m;
  R = R/2.99792458e+10;
  freqnm = R*(1.0/(n*n) - 1.0/(m*m));
  return(1.0e+08/freqnm);
}
