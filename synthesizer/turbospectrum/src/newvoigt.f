*
*     program voigts
**
**
*      real    rnewvoigt
*   10 print *,'a,v?'
*        read(*,*), a,v
*        h=rnewvoigt(a,v)
*        print *, a,v,h
*      goto 10
*      end
*
*
*
*
      REAL FUNCTION newVOIGT(A,V)
* Champion Voigt routine in speed and accuracy test by BE+NP 97-11-18
C 
C  VOIGT FUNCTION CALCULATION: HUMLICEK'S APPROXIMATION 
* The results are normalized to an area of 1/sqrt(pi).
C 
      IMPLICIT NONE
      COMPLEX*8 TAV,UAV,W4,V4 
      REAL*8    SAV,AA,VV
      REAL      A,V
C 
cc      if (A.lt.1.e-8) then
cc        newvoigt=exp(-v*v)
cc      else
      AA=A
      VV=abs(V)
      TAV=CMPLX(AA,-VV) 
      SAV=VV+AA 
      UAV=TAV*TAV 
      IF(SAV.GE.15.D0) THEN
        W4=TAV*0.5641896D0/(0.5D0+UAV)
      ELSE IF(SAV.GE.5.5D0) THEN
        W4=TAV*(1.410474D0+UAV*0.5641896D0)/(0.75D0+UAV*(3.D0+UAV))
      ELSE IF(AA.GE.0.195D0*VV-0.176D0) THEN
        W4=(16.4955D0+TAV*(20.20933D0+TAV*(11.96482D0+
     &     TAV*(3.778987D0+TAV*0.5642236D0))))/(16.4955D0+
     &     TAV*(38.82363D0+TAV*(39.27121D0+TAV*(21.69274D0+
     &     TAV*(6.699398D0+TAV)))))
      ELSE
        W4=TAV*(36183.31D0-UAV*(3321.9905D0-UAV*(1540.787D0-
     &     UAV*(219.0313D0-UAV*(35.76683D0-UAV*(1.320522D0-
     &     UAV*0.56419D0))))))
        V4=(32066.6D0-UAV*(24322.84D0-UAV*(9022.228D0-
     &     UAV*(2186.181D0-UAV*(364.2191D0-
     &     UAV*(61.57037D0-UAV*(1.841439D0-UAV)))))))
        W4=EXP(UAV)-W4/V4
      END IF
      newVOIGT=REAL(W4) 
cc      endif
C 
      RETURN 
      END 
*
