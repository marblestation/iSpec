#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
int approx();
double ply();
double partfn(double code,double T,double Ne);
double molpartfn();
double qrot(double T);
double h2opartfn(double T);

double Xi(code,T,Ne)
double code,T,Ne;
{
  static double k = 8.617084e-05;
  double xi,logt,logK,t;
  double Uab,Ua,Ub,mu;
  double Uabc,Uc,mu1,mu2;
  static double ksi = 1.380658e-23;

  t = 5040.0/T;
  logt = log10(t);

  if(approx(code,101.0,0.001) == 1) {
      Uab = molpartfn(101.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(1.0,T,Ne); 
      mu = 0.50397;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(4.4781/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,106.0,0.001) == 1) {
      Uab = molpartfn(106.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(6.0,T,Ne); 
      mu = 0.9299;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(3.470/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,107.0,0.001) == 1) {
      Uab = molpartfn(107.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(7.0,T,Ne); 
      mu = 0.94028;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(3.47/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,108.0,0.001) == 1) {
      Uab = molpartfn(108.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(8.0,T,Ne); 
      mu = 0.94820;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(4.392/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,112.0,0.001) == 1) {
      Uab = molpartfn(112.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(12.0,T,Ne); 
      mu = 0.96780;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(1.34/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,113.0,0.001) == 1) {
      Uab = molpartfn(113.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(13.0,T,Ne); 
      mu = 0.97170;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(3.06/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,114.0,0.001) == 1) {
      Uab = molpartfn(114.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(14.0,T,Ne); 
      mu = 0.97302;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(3.06/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,120.0,0.001) == 1) {
      Uab = molpartfn(120.0,T);
      Ua = partfn(1.0,T,Ne);
      Ub = partfn(20.0,T,Ne); 
      mu = 0.98321;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(1.70/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,606.0,0.001) == 1) {
      Uab = molpartfn(606.0,T);
      Ua = partfn(6.0,T,Ne);
      Ub = partfn(6.0,T,Ne); 
      mu = 6.0055;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(6.15/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,607.0,0.001) == 1) {
      Uab = molpartfn(607.0,T);
      Ua = partfn(6.0,T,Ne);
      Ub = partfn(7.0,T,Ne); 
      mu = 6.4627;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(7.66/(k*T))*pow(T,-1.5);
  }
      
  else if(approx(code,608.0,0.001) == 1) {
      Uab = molpartfn(608.0,T);
      Ua = partfn(6.0,T,Ne);
      Ub = partfn(8.0,T,Ne); 
      mu = 6.8604;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(11.108/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,614.0,0.001) == 1) {
      logK = ply(10.8445,-0.9184,0.1532,-0.3771,0.0000,0.0000,logt)-4.6400*t;
      xi = 1.0e+06*ksi*T/pow(10.0,logK);
  }

  else if(approx(code,707.0,0.001) == 1)
      xi = exp(9.763/(k*T) - 4.8696e+01 + T*(1.9224e-03 + T*(-4.9143e-07 +
	   T*(7.4692e-11 - T*4.5399e-15))) - 1.5*log(T));

  else if(approx(code,708.0,0.001) == 1)
      xi = exp(6.508/(k*T) - 4.7365e+01 + T*(1.9840e-03 + T*(-4.9280e-07 +
	   T*(7.2933e-11 - T*4.3354e-15))) - 1.5*log(T));

  else if(approx(code,722.0,0.001) == 1) {
      logK = ply(11.7421,-1.7267,1.5855,-1.3639,0.0000,0.0000,logt)-4.9000*t;
      xi = 1.0e+06*ksi*T/pow(10.0,logK);
  }

  else if(approx(code,740.0,0.001) == 1) {
      logK = ply(11.6496,-1.8280,0.8730,-0.7237,0.0000,0.0000,logt)-5.8100*t;
      xi = 1.0e+06*ksi*T/pow(10.0,logK);
  }

  else if(approx(code,808.0,0.001) == 1)
      xi = exp(5.116/(k*T) - 4.8625e+01 + T*(1.6321e-03 + T*(-3.3915e-07 +
	   T*(4.6025e-11 - T*2.5839e-15))) - 1.5*log(T));

  else if(approx(code,812.0,0.001) == 1) {
      logK = ply(9.7780,-0.0290,1.1725,-1.1724,0.0000,0.0000,logt)-3.5300*t;
      xi = 1.0e+06*ksi*T/pow(10.0,logK);
  }

  else if(approx(code,813.0,0.001) == 1) {
      Uab = molpartfn(813.0,T);
      Ua = partfn(8.0,T,Ne);
      Ub = partfn(13.0,T,Ne); 
      mu = 10.0436;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(5.27/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,814.0,0.001) == 1) {
      Uab = molpartfn(814.0,T);
      Ua = partfn(8.0,T,Ne);
      Ub = partfn(14.0,T,Ne); 
      mu = 10.1897;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(8.26/(k*T))*pow(T,-1.5);
  }

  else if(approx(code,820.0,0.001) == 1) {
      logK = ply(10.2335,0.1684,1.8170,-6.3417,2.6646,1.9124,logt)-4.7600*t;
      xi = 1.0e+06*ksi*T/pow(10.0,logK);
  }

  else if(approx(code,822.0,0.001) == 1) {
      Uab = molpartfn(822.0,T);
      Ua = partfn(8.0,T,Ne);
      Ub = partfn(22.0,T,Ne); 
      mu = 11.9921;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(6.87/(k*T))*pow(T,-1.5);
  }

  /* VO   12.17542 */
  else if(approx(code,823.0,0.001) == 1) {
      Uab = molpartfn(823.0,T);
      Ua = partfn(8.0,T,Ne);
      Ub = partfn(23.0,T,Ne); 
      mu = 12.17542;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(6.41/(k*T))*pow(T,-1.5);
  }

  /* YO   13.55929 */
  else if(approx(code,839.0,0.001) == 1) {
      Uab = molpartfn(839.0,T);
      Ua = partfn(8.0,T,Ne);
      Ub = partfn(39.0,T,Ne); 
      mu = 13.55929;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(7.29/(k*T))*pow(T,-1.5);
  }

  /* ZrO  13.61204 */
  else if(approx(code,840.0,0.001) == 1) {
      Uab = molpartfn(840.0,T);
      Ua = partfn(8.0,T,Ne);
      Ub = partfn(40.0,T,Ne); 
      mu = 13.61204;
      xi = 5.322e-21*(Uab/(Ua*Ub))*pow(mu,-1.5)*exp(7.85/(k*T))*pow(T,-1.5);
  }

  /*  else if(approx(code,840.0,0.001) == 1) {
      logK = ply(11.5031,-1.1916,0.5692,-0.1795,0.0000,0.0000,logt)-7.8500*t;
      xi = 1.0e+06*ksi*T/pow(10.0,logK);
      } */

  else if(approx(code,10108.,0.001) == 1) {
    Uabc = molpartfn(10108.,T);
    Ua = Uc = partfn(1.0,T,Ne);
    Ub = partfn(8.0,T,Ne);
    mu1 = 0.94820;   /* Mean molecular weight of OH */
    mu2 = 0.95220;   /* Mean molecular weight of H-OH = H2O */
    /* Below, bond dissociation energies of water and OH are taken from
       Ruscic et al 2002 J Chem Phys 106, page ? */
    xi = 2.83237e-41*(Uabc/(Ua*Ub*Uc))*pow(mu1,-1.5)*pow(mu2,-1.5)*
      exp(4.41/(k*T))*exp(5.096/(k*T))*pow(T,-3.0);
  }

  else {
    printf("The molecule with code = %8.1f is unknown to this program\n",code);
    exit(1);
  }

  return(xi);
}

double ply(b0,b1,b2,b3,b4,b5,x)
double b0,b1,b2,b3,b4,b5,x;
{
  double K;

  K = b0 + x*(b1 + x*(b2 + x*(b3 + x*(b4 + x*b5))));
  return(K);
}

double molpartfn(code,T)
double code,T;
{
  double a0,a1,a2,a3,a4,a5;
  double lgQ,logt,Q,t;
 
  Q = lgQ = 1.0; 
  logt = log10(5040/T);
  
  if(approx(code,101.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(1.6498,-1.6265,0.7472,-0.2751,0.0,0.0,logt);
  }

  if(approx(code,106.0,0.001) == 1) {
    /* Partition Function from Jorgensen et al (1996) */
    lgQ = ply(3.29689,-1.68566,0.117806,0.699041,-0.429731,0.0,logt);
  }

  if(approx(code,107.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(3.0735,-1.8501,0.9607,-0.3935,0.0,0.0,logt);
  }

  if(approx(code,108.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(3.0929,-1.6778,0.6743,-0.1874,0.0,0.0,logt);
  }

  if(approx(code,112.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(3.6704,-2.2682,0.9354,-0.2597,0.0,0.0,logt);
  }

  if(approx(code,113.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(3.3209,-2.5909,1.7415,-0.7636,0.0,0.0,logt);
  }

  if(approx(code,114.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(3.6908,-1.9801,0.7704,-0.2247,0.0,0.0,logt);
  }

  if(approx(code,120.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(3.8411,-2.3891,1.3578,-0.6893,0.0,0.0,logt);
  }

  if(approx(code,606.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(4.3091,-2.2406,0.4865,-0.2049,0.0,0.0,logt);
  }

  if(approx(code,607.0,0.001) == 1) {
    /* Partition Function from Jorgensen & Larsson (1990) */
    lgQ = ply(3.99679,-2.12747,0.996815,-0.294187,0.0,0.0,logt);
  }

  /*
  if(approx(code,608.0,0.001) == 1) {
    Partition Function from Sauval & Tatum (1984)
    lgQ = ply(3.6076,-1.7608,0.4172,0.0,0.0,0.0,logt);
  } */

  if(approx(code,608.0,0.001) == 1) {
    /* Partition Function for 12C16O from Goorvitch (1994) */
    lgQ = ply(3.6142,-1.7736,0.4108,0.0176,0.0,0.0,logt);
  } 

  if(approx(code,813.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(4.9191,-2.6291,0.5831,0.3163,0.0,0.0,logt);
  }

  if(approx(code,814.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) */
    lgQ = ply(4.2275,-1.9144,0.7201,-1.3099,1.1657,0.0,logt);
  }

  if(approx(code,822.0,0.001) == 1) {
    /* Partition Function from Jorgensen (1994) */
    lgQ = ply(5.27769,-2.14689,0.277945,0.215438,-0.122058,0.0,logt);
  }

  if(approx(code,823.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) VO */
    lgQ = ply(5.0687,-2.2186,0.9545,-0.4592,0.0,0.0,logt);
  }  

  if(approx(code,839.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) YO */
    lgQ = ply(4.9515,-2.0866,0.6565,-0.3082,0.0,0.0,logt);
  }  

  if(approx(code,840.0,0.001) == 1) {
    /* Partition Function from Sauval & Tatum (1984) ZrO */
    lgQ = ply(5.3279,-2.4694,0.2164,-0.2313,0.0,0.0,logt);
  }  
  
  Q = pow(10.0,lgQ);

  if(approx(code,10108.,0.001) == 1) {
    /* H2O Partition Function from Jorgensen et al. (2001) */
    Q = h2opartfn(T);
  }
  return(Q);
}



/* Partition function for water from Jorgensen et al. (2001) */
double h2opartfn(double t)
{
  static double qvib[40] = 
     { 1.0000, 1.0032, 1.0227, 1.0633, 1.1239, 1.2036, 1.3028, 1.4186, 1.5541, 
       1.7085, 1.8825, 2.0767, 2.2919, 2.5289, 2.7886, 3.0721, 3.3804, 3.7144,
       4.0752, 4.4637, 4.8807, 5.3272, 5.8038, 6.3112, 6.8499, 7.4203, 8.0227,
       8.6572, 9.3238,10.0225,10.7531,11.5153,12.3085,13.1325,13.9864,14.8698,
      15.7819,16.7218,17.6889,18.6821};
  static double T[40] = 
                 { 200.0, 400.0, 600.0, 800.0,1000.0,1200.0,1400.0,1600.0,
		  1800.0,2000.0,2200.0,2400.0,2600.0,2800.0,3000.0,3200.0,
		  3400.0,3600.0,3800.0,4000.0,4200.0,4400.0,4600.0,4800.0,
		  5000.0,5200.0,5400.0,5600.0,5800.0,6000.0,6200.0,6400.0,
		  6600.0,6800.0,7000.0,7200.0,7400.0,7600.0,7800.0,8000.0};
  double Qrot,Qvib,partfn;
  int i,jlo;

  if(t < 200) t = 200.0;
  if(t > 7800) t = 7800.0;

  jlo = (int)floor((t-200.0)/200.0);

  if(jlo < 0) jlo = 0;
  if(jlo > 38) jlo = 38;

  Qvib = qvib[jlo] + (qvib[jlo+1]-qvib[jlo])*(t-T[jlo])/200.0;
  Qrot = qrot(t);
  partfn = Qrot*Qvib;
  return(partfn);

}
 

double qrot(double T)
{
  static double pi = 3.141592;
  static double rc01 = 27.88063;
  static double rc02 = 14.52177;
  static double rc03 = 9.27771;
  static double hck = 1.43877;
  double hckt,Qra,A,B,Qrb;


  hckt = hck/T;
  Qra = rc01*rc02*rc03;
  Qra = sqrt(pi/Qra*pow(1.0/hckt,3.0));
  A = rc01;
  B = sqrt(rc02*rc03);
  Qrb = B*hckt;
  B = (1.0 - B/A)*Qrb;
  Qrb = Qra*exp(0.25*Qrb)*(1.0+B/12.0 + B*B*7.0/480.0);
  Qrb = Qrb + Qrb;
  return(Qrb);
}





