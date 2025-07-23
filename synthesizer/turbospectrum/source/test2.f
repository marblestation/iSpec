1     print*,'a, dppler?'
      read(5,*) a,dop

      x0=0.
      call voigt(a,0.,h)
      sum=h*0.005
      do i=1,1000
        x=x0+i*0.01
        v=x/dop
        call voigt(a,v,h)
        sum=sum+h*0.01
         print*,x,sum
      enddo
      print*,2.*sum
      goto 1
      end
    
      subroutine voigt(a,v,h)
*
* Calculates the quickest reliable voigt function for the current
* values of a and v, assuming the precicion of the Apollo dn10k.
* The result is normalized to an area of 1/sqrt(pi).
*
      real a,v,h,hh,hvoigt
      common/count/icount1,icount2,icount3,icount4
      data icount1/0/,icount2/0/,icount3/0/,icount4/0/
*
* on peut penser a mettre a=0 if a<0.001 par ex.
      if(a.gt.100..or.v.gt.100.) then
        h=hh(a,v)
        icount1=icount1+1
* Reliable, but not quick. (strange for a<=0.001)
      else if(a.lt.0.5.and.v.lt.4.0) then
        call qvoigt(a,v,h)
        icount4=icount4+1
      else if (a.gt.0.1) then
        call nvoigt(a,v,h)
        icount2=icount2+1
* Reasonably quick and reliable.
* Crashes on dn10k for large values of a or v
* because of numerical imprecision.
* strange for small a.
      else
       h=hvoigt(a,v)
       icount3=icount3+1
* Quick and reliable for small values of a.
      endif
      return
      end
C
C
C
      FUNCTION HVOIGT (AA,XX)
C
C---- A QUICK AND DIRTY VOIGT FUNCTION GENERATOR, GOOD TO 3 PLACES OR
C BETTER FOR AA .LT. 0.1
C USES HARRIS EXPANSION UP TO H2
C
      integer i
      real aa,a2,aq,c1,expt,h1,hvoigt,x,xx,x2,xq,yq
      real A(7),B(17),C(12)
*
      DATA A/-.636619772,-.575045176,-.517467538,-.463663043,
     &       -.413393140,-.366415179,-.322659217 /
      DATA B/-.366415179,-.273927401,-.181624005,-.094625424,
     &       -.016995346,.048484334,.100366097,.138437039,.163505313,
     &        .177120522,.181275706,.178131090,.169787631,.158124873,
     &        .144705210,.130737227,.1168 /
      DATA C/ .318309886,.331066371,.345857391,.363751789,.386756825,
     &        .415509920,.446971438,.476384118,.499908503,.515558594,
     &        .522948903,.53 /
      DATA C1/.5641895835/
*
      if(aa.gt.1.) then
    5   HVOIGT = AA/( 3.1415965*(AA*AA+XX*XX) )*1.77245385
        RETURN
      else
   10   CONTINUE
        AQ = AA
        X = ABS(XX)
        X2 = X**2
        IF (X.GT.0.5) GO TO  15
        XQ = X2/.05
        I = INT(XQ)
        XQ = XQ-I
        I = I+1
        H1 = A(I)+XQ*( A(I+1)-A(I) )
        GO TO  25
   15   IF (X.GT.2.0) GO TO  20
        XQ = (X-.5)/.1
        I = INT(XQ)
        XQ = XQ-I
        I = I+1
        H1 = B(I)+XQ*( B(I+1)-B(I) )
        GO TO  25
   20   YQ = 1./X2
        XQ = YQ/.025
        I = INT(XQ)
        XQ = XQ-I
        I = I+1
        H1 = C(I)+XQ*( C(I+1)-C(I) )
        H1 = H1*YQ
        HVOIGT = 0.0
        IF (AQ.LT.1.E-6) GO TO  25
        IF (X.GT.5.) GO TO  30
   25   EXPT = EXP(-X2)
        A2 = AQ**2
        HVOIGT = C1*EXPT*( 1.+A2*(1.-2.*X2) )
   30   HVOIGT = (HVOIGT+H1*AQ)*1.77245385
        RETURN
      endif
      END
*
*
*
      subroutine nvoigt(a,v,h)
*
*  this vectorizable voigt function is based on the paper by
*  hui, armstrong and wray, jqsrt 19, 509 (1977).  it has been
*  checked against the old landolt & boernstein standard voigt.
*  errors become significant (at the < 1 % level) around the
*  "knee" between the doppler core and the damping wings for a
*  smaller than 1.e-3.  note that, as written here, the routine
*  has no provision for the trivial case a=0. the normalization
*  is such that the integral is sqrt(pi); i.e., a plain exp(-x**2)
*  should be used at a=0.  the new landolt & boernstein gives a
*  useful expression for the small but finite a case.
*
*  coded by: a. nordlund/16-dec-82.
*
*     dimension a(n),v(n),h(n)
*
      complex z
      real a,v,h
      real a0,a1,a2,a3,a4,a5,a6
      real b0,b1,b2,b3,b4,b5,b6
*
      save a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6
*
      data a0/122.607931777104326/
      data a1/214.382388694706425/
      data a2/181.928533092181549/
      data a3/93.155580458138441/
      data a4/30.180142196210589/
      data a5/5.912626209773153/
      data a6/0.564189583562615/
      data b0/122.607931773875350/
      data b1/352.730625110963558/
      data b2/457.334478783897737/
      data b3/348.703917719495792/
      data b4/170.354001821091472/
      data b5/53.992906912940207/
      data b6/10.479857114260399/
*
*     do 100 i=1,n
        z=cmplx(a,-abs(v))
        h=real(
     &   ((((((a6*z+a5)*z+a4)*z+a3)*z+a2)*z+a1)*z+a0)
     &   /(((((((z+b6)*z+b5)*z+b4)*z+b3)*z+b2)*z+b1)*z+b0)
     &   )
* 100 continue
*
      return
      end
      FUNCTION HH(AA,UU)
*
*
*
*        THIS FUNCTION COMPUTES THE VOIGT FUNCTION USING THE ALORITHM
*        OF MATTA AND REICHEL, MATH. COMPUT. 25,339
*
*         CODED BY K.ERIKSSON (JULY-81)
*
      integer i
*
      real a,aa,a2,a1,b1,c1,d1,p2,u,uu,u2,a2pu2,a2mu2,twoau,fa2u2,hh,h1
      real pi,twopi,fourpi,sqrtpi,fpia,fpiu
      real H2N2(12), EXH2N2(12)
*
      DATA PI/3.14159265/,TWOPI/6.28318531/,FOURPI/12.5663706/,
     &     SQRTPI/1.7724539/,
     &     H2N2/0.25, 1., 2.25, 4., 6.25, 9., 12.25, 16., 20.25, 25.,
     &          30.25, 36./,
     & EXH2N2/7.7880078E-01, 3.6787944E-01, 1.0539922E-01, 1.8315638E-02
     &     , 1.9304541E-03, 1.2340980E-04, 4.7851173E-06, 1.1253517E-07,
     &       1.6052280E-09, 1.3887943E-11, 7.2877240E-14, 2.3195228E-16/
*
*
      A=AA
      U=UU
      A2=A*A
      U2=U*U
      A2PU2=A2+U2
      A2MU2=A2-U2
      FA2U2=4.*A2*U2
      HH=A/(TWOPI*A2PU2)
      H1=0.
      DO 10 I=1,12
      H1=H1 + EXH2N2(I)*(A2PU2 + H2N2(I))/((A2MU2 + H2N2(I))**2 + FA2U2)
 10   CONTINUE
      HH = HH + A*H1/PI
*
      IF(A.GT.TWOPI) RETURN
*
      TWOAU = 2.*A*U
      FPIA=FOURPI*A
      FPIU=FOURPI*U
      A1=COS(TWOAU)
      B1=SIN(TWOAU)
      C1=EXP(-FPIA) - COS(FPIU)
      D1=SIN(FPIU)
      P2=2.*EXP(A2MU2-FPIA)*((A1*C1-B1*D1)/(C1*C1+D1*D1))
      HH= HH + P2
*
      RETURN
      END
*
      subroutine qvoigt(aa,vv,hh)
*
* beware of the inversion in notation
* h(i,j) is the voigt function for v(i) and a(j).
* first group in table is in fact h(a=1.e-6,v)
      dimension h(26,11),a(11),v(26),h0(26)
      equivalence (h0(1),h(1,1))
      data a/0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5/
      data v/0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,
     & 2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0/
      data (h(i1,1),i1=1,26)/
     & 0.100000E+01,0.990050E+00,0.960790E+00,0.913931E+00,0.852144E+00,
     & 0.778801E+00,0.697676E+00,0.612626E+00,0.527292E+00,0.444858E+00,
     & 0.367879E+00,0.236928E+00,0.140858E+00,0.773047E-01,0.391639E-01,
     & 0.183156E-01,0.790705E-02,0.315111E-02,0.115923E-02,0.393669E-03,
     & 0.123410E-03,0.357128E-04,0.954016E-05,0.235258E-05,0.535535E-06,
     & 0.112535E-06/
      data (h(i1,2),i1=1,26)/
     & 0.945990E+00,0.937090E+00,0.910898E+00,0.868892E+00,0.813378E+00,
     & 0.747281E+00,0.673889E+00,0.596580E+00,0.518571E+00,0.442706E+00,
     & 0.371305E+00,0.248130E+00,0.155561E+00,0.923139E-01,0.526474E-01,
     & 0.295805E-01,0.169754E-01,0.103700E-01,0.694476E-02,0.510894E-02,
     & 0.404434E-02,0.335883E-02,0.287204E-02,0.250089E-02,0.220489E-02,
     & 0.196217E-02/
      data (h(i1,3),i1=1,26)/
     & 0.896457E+00,0.888479E+00,0.864983E+00,0.827246E+00,0.777267E+00,
     & 0.717588E+00,0.651076E+00,0.580698E+00,0.509299E+00,0.439421E+00,
     & 0.373170E+00,0.257374E+00,0.168407E+00,0.105843E+00,0.650991E-01,
     & 0.402014E-01,0.256782E-01,0.173969E-01,0.126351E-01,0.977827E-02,
     & 0.794268E-02,0.667002E-02,0.572753E-02,0.499482E-02,0.440595E-02,
     & 0.392175E-02/
      data (h(i1,4),i1=1,26)/
     & 0.850936E+00,0.843769E+00,0.822646E+00,0.788675E+00,0.743590E+00,
     & 0.689605E+00,0.629227E+00,0.565064E+00,0.499636E+00,0.435217E+00,
     & 0.373717E+00,0.264903E+00,0.179580E+00,0.118002E+00,0.765640E-01,
     & 0.501819E-01,0.339989E-01,0.242101E-01,0.182104E-01,0.143860E-01,
     & 0.118070E-01,0.996116E-02,0.857019E-02,0.747990E-02,0.660055E-02,
     & 0.587641E-02/
      data (h(i1,5),i1=1,26)/
     & 0.809020E+00,0.802567E+00,0.783538E+00,0.752895E+00,0.712146E+00,
     & 0.663223E+00,0.608322E+00,0.549739E+00,0.489710E+00,0.430271E+00,
     & 0.373153E+00,0.270928E+00,0.189247E+00,0.128895E+00,0.870899E-01,
     & 0.595313E-01,0.419269E-01,0.307919E-01,0.236535E-01,0.189182E-01,
     & 0.156268E-01,0.132245E-01,0.113944E-01,0.995197E-02,0.878554E-02,
     & 0.782374E-02/
      data (h(i1,6),i1=1,26)/
     & 0.770347E+00,0.764525E+00,0.747348E+00,0.719652E+00,0.682752E+00,
     & 0.638337E+00,0.588335E+00,0.534770E+00,0.479629E+00,0.424735E+00,
     & 0.371658E+00,0.275638E+00,0.197562E+00,0.138620E+00,0.967260E-01,
     & 0.682635E-01,0.494563E-01,0.371289E-01,0.289499E-01,0.233623E-01,
     & 0.193922E-01,0.164528E-01,0.141947E-01,0.124070E-01,0.109579E-01,
     & 0.976135E-02/
      data (h(i1,7),i1=1,26)/
     & 0.734599E+00,0.729337E+00,0.713801E+00,0.688720E+00,0.655244E+00,
     & 0.614852E+00,0.569238E+00,0.520192E+00,0.469480E+00,0.418736E+00,
     & 0.369386E+00,0.279199E+00,0.204662E+00,0.147272E+00,0.105522E+00,
     & 0.763959E-01,0.565857E-01,0.432110E-01,0.340875E-01,0.277075E-01,
     & 0.230945E-01,0.196393E-01,0.169661E-01,0.148411E-01,0.131146E-01,
     & 0.116869E-01/
      data (h(i1,8),i1=1,26)/
     & 0.701496E+00,0.696730E+00,0.682651E+00,0.659895E+00,0.629471E+00,
     & 0.592673E+00,0.550998E+00,0.506027E+00,0.459332E+00,0.412382E+00,
     & 0.366469E+00,0.281758E+00,0.210674E+00,0.154939E+00,0.113528E+00,
     & 0.839489E-01,0.633170E-01,0.490312E-01,0.390566E-01,0.319442E-01,
     & 0.267257E-01,0.227776E-01,0.197037E-01,0.172506E-01,0.152527E-01,
     & 0.135982E-01/
      data (h(i1,9),i1=1,26)/
     & 0.670788E+00,0.666463E+00,0.653680E+00,0.632996E+00,0.605295E+00,
     & 0.571717E+00,0.533581E+00,0.492289E+00,0.449244E+00,0.405763E+00,
     & 0.363020E+00,0.283443E+00,0.215711E+00,0.161702E+00,0.120793E+00,
     & 0.909445E-01,0.696549E-01,0.545853E-01,0.438494E-01,0.360644E-01,
     & 0.302788E-01,0.258620E-01,0.224029E-01,0.196319E-01,0.173696E-01,
     & 0.154930E-01/
      data (h(i1,10),i1=1,26)/
     & 0.642252E+00,0.638320E+00,0.626692E+00,0.607858E+00,0.582593E+00,
     & 0.551903E+00,0.516953E+00,0.478988E+00,0.439260E+00,0.398955E+00,
     & 0.359136E+00,0.284370E+00,0.219876E+00,0.167638E+00,0.127364E+00,
     & 0.974063E-01,0.756067E-01,0.598715E-01,0.484603E-01,0.400614E-01,
     & 0.337474E-01,0.288874E-01,0.250596E-01,0.219818E-01,0.194626E-01,
     & 0.173693E-01/
      data (h(i1,11),i1=1,26)/
     & 0.615690E+00,0.612109E+00,0.601513E+00,0.584333E+00,0.561252E+00,
     & 0.533157E+00,0.501079E+00,0.466127E+00,0.429418E+00,0.392021E+00,
     & 0.354900E+00,0.284638E+00,0.223262E+00,0.172820E+00,0.133288E+00,
     & 0.103359E+00,0.811817E-01,0.648898E-01,0.528850E-01,0.439298E-01,
     & 0.371264E-01,0.318490E-01,0.276699E-01,0.242970E-01,0.215291E-01,
     & 0.192249E-01/
*
      if (vv.lt.1.0) then
        iv=int(vv*10.)+1
        tt=(vv-v(iv))*10.
      else
        iv=int(vv*5.)+6
        tt=(vv-v(iv))*5.
      endif
cc      tt=(vv-v(iv))/(v(iv+1)-v(iv))
      if (aa.le.0.000001) then
       hh=tt*(h0(iv+1)-h0(iv))+h0(iv)
      else
       ia=int(aa*20.)+1
cc       uu=(aa-a(ia))/(a(ia+1)-a(ia))
       uu=(aa-a(ia))*20.
       hh=(1.-tt)*(1.-uu)*h(iv,ia)+tt*(1.-uu)*h(iv+1,ia)+
     &  tt*uu*h(iv+1,ia+1)+(1.-tt)*uu*h(iv,ia+1)
      endif
      return
      end


