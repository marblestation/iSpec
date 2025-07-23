!---------*************---------------------------------------------+-------
   MODULE cubint_module
!
!SUBROUTINES: (contained)
!   cubintg     : Cubic integration   routine
!   cubintp     : Cubic interpolation routine
!
   INTERFACE cubintp
     module procedure cubintp04, cubintp14, cubintp08, cubintp18
   end interface ! cubintp
!
   CONTAINS

!-------------*******-----------------------------------------------+-------
   SUBROUTINE cubintg(yin, xin, uin, qin, nin, intmod, ibound)
!
! This subroutine computes the integral for a given general 
! table of argument and function values
!
! Author: M. Steffen, July 2005
!
!   input  parameters:
!      xin      - vector of argument values
!      yin      - vector of function values
!      uin(1)   - input dy/dx at inner boundary (ibound=1)
!      uin(nin) - input dy/dx at outer boundary (ibound=1)
!      qin(1)   - integration constant
!      nin      - length of all vectors
!      intmod   - interpolation mode
!                 INTMOD = 1: Derivatives at inner points are given by
!                             slope of the secant through neighbour points.
!                 INTMOD = 2: Derivatives at inner points are given by
!                             slope of a parabola through the point
!                             and its neighbours. If the xin-values are
!                             equidistant, INTMOD=1 and INTMOD=2 have the
!                             same effect.
!                 INTMOD = 3: The derivatives of INTMOD = 2 are modified to
!                             assure monotonic behaviour of the 
!                             interpolating cubic polynomial in all 
!                             intervals.
!      ibound   - if /= 0, slope at the boundary must be given
!                    as uin(1), uin(nin) at subroutine entry;
!                    otherwise boundary slopes are computed internally,
!                    depending on intmod.
!   output parameters:
!      uin      - vector with derivatives dy/dx
!      qin      - vector with integrals of yin over xin
!                 from xin(1) .. xin(i)
!
      IMPLICIT NONE
!
!--- Dummy arguments:
!
      INTEGER,                   INTENT(IN)    :: ibound, intmod, nin
      REAL,      DIMENSION(nin), INTENT(IN)    :: xin, yin
      REAL,      DIMENSION(nin), INTENT(INOUT) :: uin, qin
!
!--- Local scalar variables:
      INTEGER         :: i, n, n1, n2
      REAL            :: dxm, dxn, sxn, sxm, s3, s4, s5, sx1, sx2
      REAL            :: dxi, dsi
!
      IF (nin < 2) THEN
         WRITE(*,*) 'CUBINTG error: NIN must be greater than 1'
         RETURN
      ENDIF
      n1 = 2
      n2 = nin - 1
      dxm =  xin(n1)-xin(n1-1)
      sxm = (yin(n1)-yin(n1-1))/dxm
      sx1 = sxm
!
      IF (nin == 2) THEN
! --- Derivatives at boundary points:
         IF (ibound == 0) THEN
            uin(n1-1) = sx1
            uin(n2+1) = sx1
         ENDIF
      ELSEIF(INTMOD == 1) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           uin(n) = (yin(n+1)-yin(n-1))/(xin(n+1)-xin(n-1))
         ENDDO
! --- Derivatives at boundary points, slope of boundary interval:
         IF (ibound == 0) THEN
            uin(n1-1) = (yin(n1)-yin(n1-1))/(xin(n1)-xin(n1-1))
            uin(n2+1) = (yin(n2)-yin(n2+1))/(xin(n2)-xin(n2+1))
         ENDIF
      ELSEIF(INTMOD == 2) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           uin(n) = (sxn*dxm+sxm*dxn)/(dxm+dxn)
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, slope of parabola:
         IF (ibound == 0) THEN
            uin(n1-1) = 2.0*sx1 - 1.0*uin(n1)
            uin(n2+1) = 2.0*sx2 - 1.0*uin(n2)
         ENDIF
      ELSEIF(INTMOD == 3) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           s3  = (sxn*dxm+sxm*dxn)/(2.D0*(dxm+dxn))
           s4  = MIN(s3,sxn,sxm)
           s5  = MAX(s3,sxn,sxm)
           uin(n) = 2.D0*(MAX(s4,0.D0)+MIN(s5,0.D0))
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, vanishing second derivative:
         IF (ibound == 0) THEN
            uin(n1-1) = 1.5*sx1 - 0.5*uin(n1)
            uin(n2+1) = 1.5*sx2 - 0.5*uin(n2)
         ENDIF
      ENDIF
!
! Main integration loop
!
      DO  i=n1-1,n2
        dxi = xin(i+1)-xin(i)
        dsi = dxi*(0.5D0*(yin(i+1)+yin(i)) - dxi*(uin(i+1)-uin(i))/12.D0)
        qin(i+1)=qin(i) + dsi
      ENDDO
!      
      END SUBROUTINE cubintg
!--------------------*******----------------------------------------+-------

!-------------*********---------------------------------------------+-------
   SUBROUTINE cubintp04(yin, xin, uin, xout, yout, iout, &
                        nin, intmod, ibound, iquant, iflag, extrap)
!
! This subroutine computes interpolated y values (yout) for a given 
! x-value, xout (scalar version).
! It also returns the corresponding interval index (iout) if iflag==0
!
! Author: M. Steffen, July 2005
!
!   input  parameters:
!      xin      - vector of argument values
!      yin      - vector of function values
!      uin(1)   - input dy/dx at inner boundary (ibound=1)
!      uin(nin) - input dy/dx at outer boundary (ibound=1)
!      iout     - vector with interval indices for loactions xout (iflag /=0)
!      nin      - length of all  input vectors
!      intmod   - interpolation mode
!                 INTMOD = 1: Derivatives at inner points are given by
!                             slope of the secant through neighbour points.
!                 INTMOD = 2: Derivatives at inner points are given by
!                             slope of a parabola through the point
!                             and its neighbours. If the xin-values are
!                             equidistant, INTMOD=1 and INTMOD=2 have the
!                             same effect.
!                 INTMOD = 3: The derivatives of INTMOD = 2 are modified to
!                             assure monotonic behaviour of the 
!                             interpolating cubic polynomial in all 
!                             intervals.
!      ibound   - if /= 0, slope at the boundary must be given
!                 as uin(1), uin(nin) at subroutine entry;
!                 otherwise boundary slopes are computed internally,
!                 depending on intmod.
!      iquant   - 0 -> interpolate y
!                 1 -> interpolate dy/dx
!      iflag    - if /= 0, the interval index must be given
!                 as iout  at subroutine entry;
!                 otherwise interval indices are computed internally,
!   output parameters:
!      iout     - vector with interval indices for loactions xout
!      yout     - vector with interpolated  y   at loactions xout
!      extrap   - logical; true if extrapolation was performed
!
      IMPLICIT NONE
!
!--- Dummy arguments:
!
      INTEGER,                    INTENT(IN)    :: ibound, iflag, intmod, &
                                                   iquant, nin
      INTEGER,                    INTENT(INOUT) :: iout
      REAL,      DIMENSION(nin),  INTENT(IN)    :: xin, yin
      REAL,      DIMENSION(nin),  INTENT(INOUT) :: uin
      REAL,                       INTENT(IN)    :: xout
      REAL,                       INTENT(OUT)   :: yout
      LOGICAL                                   :: extrap

!--- Local scalar variables:
      INTEGER         :: i, n, n1, n2
      REAL            :: dx, dxm, dxn, sxn, sxm, s3, s4, s5, sx1, sx2
      REAL            :: delta, dxrel, dxi, dsi, u1, u2, y1, y2
!
      IF (nin < 2) THEN
         WRITE(*,*) 'CUBINTP error: NIN must be greater than 1'
         RETURN
      ENDIF
!
      extrap=.false.
!
      n1 = 2
      n2 = nin - 1
      dxm =  xin(n1)-xin(n1-1)
      sxm = (yin(n1)-yin(n1-1))/dxm
      sx1 = sxm
!
      IF (nin == 2) THEN
! --- Derivatives at boundary points:
         IF (ibound == 0) THEN
            uin(n1-1) = sx1
            uin(n2+1) = sx1
         ENDIF
      ELSEIF(INTMOD == 1) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           uin(n) = (yin(n+1)-yin(n-1))/(xin(n+1)-xin(n-1))
         ENDDO
! --- Derivatives at boundary points, slope of boundary interval:
         IF (ibound == 0) THEN
            uin(n1-1) = (yin(n1)-yin(n1-1))/(xin(n1)-xin(n1-1))
            uin(n2+1) = (yin(n2)-yin(n2+1))/(xin(n2)-xin(n2+1))
         ENDIF
      ELSEIF(INTMOD == 2) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           uin(n) = (sxn*dxm+sxm*dxn)/(dxm+dxn)
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, slope of parabola:
         IF (ibound == 0) THEN
            uin(n1-1) = 2.0*sx1 - 1.0*uin(n1)
            uin(n2+1) = 2.0*sx2 - 1.0*uin(n2)
         ENDIF
      ELSEIF(INTMOD == 3) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           s3  = (sxn*dxm+sxm*dxn)/(2.D0*(dxm+dxn))
           s4  = MIN(s3,sxn,sxm)
           s5  = MAX(s3,sxn,sxm)
           uin(n) = 2.D0*(MAX(s4,0.D0)+MIN(s5,0.D0))
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, vanishing second derivative:
         IF (ibound == 0) THEN
            uin(n1-1) = 1.5*sx1 - 0.5*uin(n1)
            uin(n2+1) = 1.5*sx2 - 0.5*uin(n2)
         ENDIF
      ENDIF
!
! Find interval indices if not given as input
!
      IF (iflag == 0 .AND. xin(1) < xin(nin)) THEN
         iout=0
         DO i=1,nin
            IF (xin(i) <= xout) iout=i
         ENDDO
      ENDIF
!      
      IF (iflag == 0 .AND. xin(1) > xin(nin)) THEN
         iout=0
         DO i=1,nin
            IF (xin(i) >= xout) iout=i
         ENDDO
      ENDIF
!      
! Main interpolation loop
!
         i = iout
         IF     (i==0  ) THEN
! --- Linear extrapolation ---
           dx = xout-xin(1)
           u1 = uin(1)
           IF (iquant==0) THEN
             yout = yin(1  ) + u1*dx
           ELSE
             yout = u1
           ENDIF
           extrap=.true.
         ELSEIF (i==nin) THEN
! --- Linear extrapolation ---
           dx = xout-xin(nin)
           u2 = uin(nin)
           IF (iquant==0) THEN
             yout = yin(nin) + u2*dx
           ELSE
             yout = u2
           ENDIF
           extrap=.true.
         ELSE
! --- Cubic interpolation ---
           y1 = yin(i)
           y2 = yin(i+1)
           u1 = uin(i)
           u2 = uin(i+1)
           delta =  xin(i+1)-xin(i)
           dxrel = (xout-xin(i))/delta
           IF (iquant==0) THEN
! --- Intepolate y ---
             yout = (((     (u1+u2)*delta - 2.0*(y2-y1))*dxrel  + &
                      ( -2.0*u1-u2)*delta + 3.0*(y2-y1))*dxrel  + &
                             u1    *delta              )*dxrel  + y1
           ELSE
! --- Intepolate dy/dx ---
             yout =  ( 3.0*u1+3.0*u2 - 6.0*(y2-y1)/delta)*dxrel + &
                     (-4.0*u1-2.0*u2 + 6.0*(y2-y1)/delta)*dxrel + u1
           ENDIF
         ENDIF
!
      END SUBROUTINE cubintp04
!--------------------*********--------------------------------------+-------

!-------------*********---------------------------------------------+-------
   SUBROUTINE cubintp08(yin, xin, uin, xout, yout, iout, &
                        nin, intmod, ibound, iquant, iflag, extrap)
!
! This subroutine computes interpolated y values (yout) for a given 
! x-value, xout (scalar version).
! It also returns the corresponding interval index (iout) if iflag==0
!
! Author: M. Steffen, July 2005
!
!   input  parameters:
!      xin      - vector of argument values
!      yin      - vector of function values
!      uin(1)   - input dy/dx at inner boundary (ibound=1)
!      uin(nin) - input dy/dx at outer boundary (ibound=1)
!      iout     - vector with interval indices for loactions xout (iflag /=0)
!      nin      - length of all  input vectors
!      intmod   - interpolation mode
!                 INTMOD = 1: Derivatives at inner points are given by
!                             slope of the secant through neighbour points.
!                 INTMOD = 2: Derivatives at inner points are given by
!                             slope of a parabola through the point
!                             and its neighbours. If the xin-values are
!                             equidistant, INTMOD=1 and INTMOD=2 have the
!                             same effect.
!                 INTMOD = 3: The derivatives of INTMOD = 2 are modified to
!                             assure monotonic behaviour of the 
!                             interpolating cubic polynomial in all 
!                             intervals.
!      ibound   - if /= 0, slope at the boundary must be given
!                 as uin(1), uin(nin) at subroutine entry;
!                 otherwise boundary slopes are computed internally,
!                 depending on intmod.
!      iquant   - 0 -> interpolate y
!                 1 -> interpolate dy/dx
!      iflag    - if /= 0, the interval index must be given
!                 as iout  at subroutine entry;
!                 otherwise interval indices are computed internally,
!   output parameters:
!      iout     - vector with interval indices for loactions xout
!      yout     - vector with interpolated  y   at loactions xout
!      extrap   - logical; true if extrapolation was performed
!
      IMPLICIT NONE
!
!--- Dummy arguments:
!
      INTEGER, PARAMETER :: k8=SELECTED_REAL_KIND(12)   ! "double precision"
!
      INTEGER,                         INTENT(IN)    :: ibound, iflag, intmod, &
                                                        iquant, nin
      INTEGER,                         INTENT(INOUT) :: iout
      REAL(KIND=k8),  DIMENSION(nin),  INTENT(IN)    :: xin, yin
      REAL(KIND=k8),  DIMENSION(nin),  INTENT(INOUT) :: uin
      REAL(KIND=k8),                   INTENT(IN)    :: xout
      REAL(KIND=k8),                   INTENT(OUT)   :: yout
      LOGICAL                                   :: extrap

!--- Local scalar variables:
      INTEGER         :: i, n, n1, n2
      REAL(KIND=k8)   :: dx, dxm, dxn, sxn, sxm, s3, s4, s5, sx1, sx2
      REAL(KIND=k8)   :: delta, dxrel, dxi, dsi, u1, u2, y1, y2
!
      IF (nin < 2) THEN
         WRITE(*,*) 'CUBINTP error: NIN must be greater than 1'
         RETURN
      ENDIF
!
      extrap=.false.
!
      n1 = 2
      n2 = nin - 1
      dxm =  xin(n1)-xin(n1-1)
      sxm = (yin(n1)-yin(n1-1))/dxm
      sx1 = sxm
!
      IF (nin == 2) THEN
! --- Derivatives at boundary points:
         IF (ibound == 0) THEN
            uin(n1-1) = sx1
            uin(n2+1) = sx1
         ENDIF
      ELSEIF(INTMOD == 1) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           uin(n) = (yin(n+1)-yin(n-1))/(xin(n+1)-xin(n-1))
         ENDDO
! --- Derivatives at boundary points, slope of boundary interval:
         IF (ibound == 0) THEN
            uin(n1-1) = (yin(n1)-yin(n1-1))/(xin(n1)-xin(n1-1))
            uin(n2+1) = (yin(n2)-yin(n2+1))/(xin(n2)-xin(n2+1))
         ENDIF
      ELSEIF(INTMOD == 2) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           uin(n) = (sxn*dxm+sxm*dxn)/(dxm+dxn)
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, slope of parabola:
         IF (ibound == 0) THEN
            uin(n1-1) = 2.0*sx1 - 1.0*uin(n1)
            uin(n2+1) = 2.0*sx2 - 1.0*uin(n2)
         ENDIF
      ELSEIF(INTMOD == 3) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           s3  = (sxn*dxm+sxm*dxn)/(2.D0*(dxm+dxn))
           s4  = MIN(s3,sxn,sxm)
           s5  = MAX(s3,sxn,sxm)
           uin(n) = 2.D0*(MAX(s4,0.D0)+MIN(s5,0.D0))
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, vanishing second derivative:
         IF (ibound == 0) THEN
            uin(n1-1) = 1.5*sx1 - 0.5*uin(n1)
            uin(n2+1) = 1.5*sx2 - 0.5*uin(n2)
         ENDIF
      ENDIF
!
! Find interval indices if not given as input
!
      IF (iflag == 0 .AND. xin(1) < xin(nin)) THEN
         iout=0
         DO i=1,nin
            IF (xin(i) <= xout) iout=i
         ENDDO
      ENDIF
!      
      IF (iflag == 0 .AND. xin(1) > xin(nin)) THEN
         iout=0
         DO i=1,nin
            IF (xin(i) >= xout) iout=i
         ENDDO
      ENDIF
!      
! Main interpolation loop
!
         i = iout
         IF     (i==0  ) THEN
! --- Linear extrapolation ---
           dx = xout-xin(1)
           u1 = uin(1)
           IF (iquant==0) THEN
             yout = yin(1  ) + u1*dx
           ELSE
             yout = u1
           ENDIF
           extrap=.true.
         ELSEIF (i==nin) THEN
! --- Linear extrapolation ---
           dx = xout-xin(nin)
           u2 = uin(nin)
           IF (iquant==0) THEN
             yout = yin(nin) + u2*dx
           ELSE
             yout = u2
           ENDIF
           extrap=.true.
         ELSE
! --- Cubic interpolation ---
           y1 = yin(i)
           y2 = yin(i+1)
           u1 = uin(i)
           u2 = uin(i+1)
           delta =  xin(i+1)-xin(i)
           dxrel = (xout-xin(i))/delta
           IF (iquant==0) THEN
! --- Intepolate y ---
             yout = (((     (u1+u2)*delta - 2.0*(y2-y1))*dxrel  + &
                      ( -2.0*u1-u2)*delta + 3.0*(y2-y1))*dxrel  + &
                             u1    *delta              )*dxrel  + y1
           ELSE
! --- Intepolate dy/dx ---
             yout =  ( 3.0*u1+3.0*u2 - 6.0*(y2-y1)/delta)*dxrel + &
                     (-4.0*u1-2.0*u2 + 6.0*(y2-y1)/delta)*dxrel + u1
           ENDIF
         ENDIF
!
      END SUBROUTINE cubintp08
!--------------------*********--------------------------------------+-------

!-------------*********---------------------------------------------+-------
   SUBROUTINE cubintp14(yin, xin, uin, xout, yout, iout, &
                        nin, nout, intmod, ibound, iquant, iflag, extrap)
!
! This subroutine computes interpolated y values (yout) for a given 
! x-vector (xout) (vector version).
! It also returns the corresponding interval indices (iout) if iflag==0
!
! Author: M. Steffen, July 2005
!
!   input  parameters:
!      xin      - vector of argument values
!      yin      - vector of function values
!      uin(1)   - input dy/dx at inner boundary (ibound=1)
!      uin(nin) - input dy/dx at outer boundary (ibound=1)
!      iout     - vector with interval indices for loactions xout (iflag /=0)
!      nin      - length of all  input vectors
!      nout     - length of all output vectors
!      intmod   - interpolation mode
!                 INTMOD = 1: Derivatives at inner points are given by
!                             slope of the secant through neighbour points.
!                 INTMOD = 2: Derivatives at inner points are given by
!                             slope of a parabola through the point
!                             and its neighbours. If the xin-values are
!                             equidistant, INTMOD=1 and INTMOD=2 have the
!                             same effect.
!                 INTMOD = 3: The derivatives of INTMOD = 2 are modified to
!                             assure monotonic behaviour of the 
!                             interpolating cubic polynomial in all 
!                             intervals.
!      ibound   - if /= 0, slope at the boundary must be given
!                 as uin(1), uin(nin) at subroutine entry;
!                 otherwise boundary slopes are computed internally,
!                 depending on intmod.
!      iquant   - 0 -> interpolate y
!                 1 -> interpolate dy/dx
!      iflag    - if /= 0, the vector of interval indices must be given
!                 as iout(1:nout)  at subroutine entry;
!                 otherwise interval indices are computed internally,
!   output parameters:
!      iout     - vector with interval indices for loactions xout
!      yout     - vector with interpolated  y   at loactions xout
!      extrap   - logical; true if extrapolation was performed
!
      IMPLICIT NONE
!
!--- Dummy arguments:
!
      INTEGER,                    INTENT(IN)    :: ibound, iflag, intmod, &
                                                   iquant, nin, nout
      INTEGER,   DIMENSION(nout), INTENT(INOUT) :: iout
      REAL,      DIMENSION(nin),  INTENT(IN)    :: xin, yin
      REAL,      DIMENSION(nin),  INTENT(INOUT) :: uin
      REAL,      DIMENSION(nout), INTENT(IN)    :: xout
      REAL,      DIMENSION(nout), INTENT(OUT)   :: yout

!--- Local scalar variables:
      INTEGER         :: i, n, n1, n2
      REAL            :: dx, dxm, dxn, sxn, sxm, s3, s4, s5, sx1, sx2
      REAL            :: delta, dxrel, dxi, dsi, u1, u2, y1, y2
      LOGICAL                                   :: extrap
!
      IF (nin < 2) THEN
         WRITE(*,*) 'CUBINTP error: NIN must be greater than 1'
         RETURN
      ENDIF
!
      extrap=.false.
!
      n1 = 2
      n2 = nin - 1
      dxm =  xin(n1)-xin(n1-1)
      sxm = (yin(n1)-yin(n1-1))/dxm
      sx1 = sxm
!
      IF (nin == 2) THEN
! --- Derivatives at boundary points:
         IF (ibound == 0) THEN
            uin(n1-1) = sx1
            uin(n2+1) = sx1
         ENDIF
      ELSEIF(INTMOD == 1) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           uin(n) = (yin(n+1)-yin(n-1))/(xin(n+1)-xin(n-1))
         ENDDO
! --- Derivatives at boundary points, slope of boundary interval:
         IF (ibound == 0) THEN
            uin(n1-1) = (yin(n1)-yin(n1-1))/(xin(n1)-xin(n1-1))
            uin(n2+1) = (yin(n2)-yin(n2+1))/(xin(n2)-xin(n2+1))
         ENDIF
      ELSEIF(INTMOD == 2) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           uin(n) = (sxn*dxm+sxm*dxn)/(dxm+dxn)
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, slope of parabola:
         IF (ibound == 0) THEN
            uin(n1-1) = 2.0*sx1 - 1.0*uin(n1)
            uin(n2+1) = 2.0*sx2 - 1.0*uin(n2)
         ENDIF
      ELSEIF(INTMOD == 3) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           s3  = (sxn*dxm+sxm*dxn)/(2.D0*(dxm+dxn))
           s4  = MIN(s3,sxn,sxm)
           s5  = MAX(s3,sxn,sxm)
           uin(n) = 2.D0*(MAX(s4,0.D0)+MIN(s5,0.D0))
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, vanishing second derivative:
         IF (ibound == 0) THEN
            uin(n1-1) = 1.5*sx1 - 0.5*uin(n1)
            uin(n2+1) = 1.5*sx2 - 0.5*uin(n2)
         ENDIF
      ENDIF
!
! Find interval indices if not given as input
!
      IF (iflag == 0 .AND. xin(1) < xin(nin)) THEN
         DO  n=1,nout
         iout(n)=0
            DO i=1,nin
               IF (xin(i) <= xout(n)) iout(n)=i
            ENDDO
         ENDDO
      ENDIF
!      
      IF (iflag == 0 .AND. xin(1) > xin(nin)) THEN
         DO  n=1,nout
         iout(n)=0
            DO i=1,nin
               IF (xin(i) >= xout(n)) iout(n)=i
            ENDDO
         ENDDO
      ENDIF
!      
! Main interpolation loop
!
      DO n=1,nout
         i = iout(n)
         IF     (i==0  ) THEN
! --- Linear extrapolation ---
           dx = xout(n)-xin(1)
           u1 = uin(1)
           IF (iquant==0) THEN
             yout(n) = yin(1  ) + u1*dx
           ELSE
             yout(n) = u1
           ENDIF
           extrap=.true.
         ELSEIF (i==nin) THEN
! --- Linear extrapolation ---
           dx = xout(n)-xin(nin)
           u2 = uin(nin)
           IF (iquant==0) THEN
             yout(n) = yin(nin) + u2*dx
           ELSE
             yout(n) = u2
           ENDIF
           extrap=.true.
         ELSE
! --- Cubic interpolation ---
           y1 = yin(i)
           y2 = yin(i+1)
           u1 = uin(i)
           u2 = uin(i+1)
           delta =  xin(i+1)-xin(i)
           dxrel = (xout(n)-xin(i))/delta
           IF (iquant==0) THEN
! --- Intepolate y ---
             yout(n) = (((     (u1+u2)*delta - 2.0*(y2-y1))*dxrel  + &
                         ( -2.0*u1-u2)*delta + 3.0*(y2-y1))*dxrel  + &
                                u1    *delta              )*dxrel  + y1
           ELSE
! --- Intepolate dy/dx ---
             yout(n) =  ( 3.0*u1+3.0*u2 - 6.0*(y2-y1)/delta)*dxrel + &
                        (-4.0*u1-2.0*u2 + 6.0*(y2-y1)/delta)*dxrel + u1
           ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE cubintp14
!--------------------*********--------------------------------------+-------

!-------------*********---------------------------------------------+-------
   SUBROUTINE cubintp18(yin, xin, uin, xout, yout, iout, &
                        nin, nout, intmod, ibound, iquant, iflag, extrap)
!
! This subroutine computes interpolated y values (yout) for a given 
! x-vector (xout) (vector version).
! It also returns the corresponding interval indices (iout) if iflag==0
!
! Author: M. Steffen, July 2005
!
!   input  parameters:
!      xin      - vector of argument values
!      yin      - vector of function values
!      uin(1)   - input dy/dx at inner boundary (ibound=1)
!      uin(nin) - input dy/dx at outer boundary (ibound=1)
!      iout     - vector with interval indices for loactions xout (iflag /=0)
!      nin      - length of all  input vectors
!      nout     - length of all output vectors
!      intmod   - interpolation mode
!                 INTMOD = 1: Derivatives at inner points are given by
!                             slope of the secant through neighbour points.
!                 INTMOD = 2: Derivatives at inner points are given by
!                             slope of a parabola through the point
!                             and its neighbours. If the xin-values are
!                             equidistant, INTMOD=1 and INTMOD=2 have the
!                             same effect.
!                 INTMOD = 3: The derivatives of INTMOD = 2 are modified to
!                             assure monotonic behaviour of the 
!                             interpolating cubic polynomial in all 
!                             intervals.
!      ibound   - if /= 0, slope at the boundary must be given
!                 as uin(1), uin(nin) at subroutine entry;
!                 otherwise boundary slopes are computed internally,
!                 depending on intmod.
!      iquant   - 0 -> interpolate y
!                 1 -> interpolate dy/dx
!      iflag    - if /= 0, the vector of interval indices must be given
!                 as iout(1:nout)  at subroutine entry;
!                 otherwise interval indices are computed internally,
!   output parameters:
!      iout     - vector with interval indices for loactions xout
!      yout     - vector with interpolated  y   at loactions xout
!      extrap   - logical; true if extrapolation was performed
!
      IMPLICIT NONE
!
!--- Dummy arguments:
!
      INTEGER, PARAMETER :: k8=SELECTED_REAL_KIND(12)   ! "double precision"
!
      INTEGER,                    INTENT(IN)    :: ibound, iflag, intmod, &
                                                   iquant, nin, nout
      INTEGER,        DIMENSION(nout), INTENT(INOUT) :: iout
      REAL(KIND=k8),  DIMENSION(nin),  INTENT(IN)    :: xin, yin
      REAL(KIND=k8),  DIMENSION(nin),  INTENT(INOUT) :: uin
      REAL(KIND=k8),  DIMENSION(nout), INTENT(IN)    :: xout
      REAL(KIND=k8),  DIMENSION(nout), INTENT(OUT)   :: yout
      LOGICAL                                   :: extrap

!--- Local scalar variables:
      INTEGER         :: i, n, n1, n2
      REAL(KIND=k8)   :: dx, dxm, dxn, sxn, sxm, s3, s4, s5, sx1, sx2
      REAL(KIND=k8)   :: delta, dxrel, dxi, dsi, u1, u2, y1, y2
!
      IF (nin < 2) THEN
         WRITE(*,*) 'CUBINTP error: NIN must be greater than 1'
         RETURN
      ENDIF
!
      extrap=.false.
!
      n1 = 2
      n2 = nin - 1
      dxm =  xin(n1)-xin(n1-1)
      sxm = (yin(n1)-yin(n1-1))/dxm
      sx1 = sxm
!
      IF (nin == 2) THEN
! --- Derivatives at boundary points:
         IF (ibound == 0) THEN
            uin(n1-1) = sx1
            uin(n2+1) = sx1
         ENDIF
      ELSEIF(INTMOD == 1) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           uin(n) = (yin(n+1)-yin(n-1))/(xin(n+1)-xin(n-1))
         ENDDO
! --- Derivatives at boundary points, slope of boundary interval:
         IF (ibound == 0) THEN
            uin(n1-1) = (yin(n1)-yin(n1-1))/(xin(n1)-xin(n1-1))
            uin(n2+1) = (yin(n2)-yin(n2+1))/(xin(n2)-xin(n2+1))
         ENDIF
      ELSEIF(INTMOD == 2) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           uin(n) = (sxn*dxm+sxm*dxn)/(dxm+dxn)
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, slope of parabola:
         IF (ibound == 0) THEN
            uin(n1-1) = 2.0*sx1 - 1.0*uin(n1)
            uin(n2+1) = 2.0*sx2 - 1.0*uin(n2)
         ENDIF
      ELSEIF(INTMOD == 3) THEN
! --- Derivatives at interior points:
         DO n=n1, n2
           dxn =  xin(n+1)-xin(n)
           sxn = (yin(n+1)-yin(n))/dxn
           s3  = (sxn*dxm+sxm*dxn)/(2.D0*(dxm+dxn))
           s4  = MIN(s3,sxn,sxm)
           s5  = MAX(s3,sxn,sxm)
           uin(n) = 2.D0*(MAX(s4,0.D0)+MIN(s5,0.D0))
           dxm = dxn
           sxm = sxn
         ENDDO
         sx2 = sxn
! --- Derivatives at boundary points, vanishing second derivative:
         IF (ibound == 0) THEN
            uin(n1-1) = 1.5*sx1 - 0.5*uin(n1)
            uin(n2+1) = 1.5*sx2 - 0.5*uin(n2)
         ENDIF
      ENDIF
!
! Find interval indices if not given as input
!
      IF (iflag == 0 .AND. xin(1) < xin(nin)) THEN
         DO  n=1,nout
         iout(n)=0
            DO i=1,nin
               IF (xin(i) <= xout(n)) iout(n)=i
            ENDDO
         ENDDO
      ENDIF
!      
      IF (iflag == 0 .AND. xin(1) > xin(nin)) THEN
         DO  n=1,nout
         iout(n)=0
            DO i=1,nin
               IF (xin(i) >= xout(n)) iout(n)=i
            ENDDO
         ENDDO
      ENDIF
!      
! Main interpolation loop
!
      DO n=1,nout
         i = iout(n)
         IF     (i==0  ) THEN
! --- Linear extrapolation ---
           dx = xout(n)-xin(1)
           u1 = uin(1)
           IF (iquant==0) THEN
             yout(n) = yin(1  ) + u1*dx
           ELSE
             yout(n) = u1
           ENDIF
           extrap=.true.
         ELSEIF (i==nin) THEN
! --- Linear extrapolation ---
           dx = xout(n)-xin(nin)
           u2 = uin(nin)
           IF (iquant==0) THEN
             yout(n) = yin(nin) + u2*dx
           ELSE
             yout(n) = u2
           ENDIF
           extrap=.true.
         ELSE
! --- Cubic interpolation ---
           y1 = yin(i)
           y2 = yin(i+1)
           u1 = uin(i)
           u2 = uin(i+1)
           delta =  xin(i+1)-xin(i)
           dxrel = (xout(n)-xin(i))/delta
           IF (iquant==0) THEN
! --- Intepolate y ---
             yout(n) = (((     (u1+u2)*delta - 2.0*(y2-y1))*dxrel  + &
                         ( -2.0*u1-u2)*delta + 3.0*(y2-y1))*dxrel  + &
                                u1    *delta              )*dxrel  + y1
           ELSE
! --- Intepolate dy/dx ---
             yout(n) =  ( 3.0*u1+3.0*u2 - 6.0*(y2-y1)/delta)*dxrel + &
                        (-4.0*u1-2.0*u2 + 6.0*(y2-y1)/delta)*dxrel + u1
           ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE cubintp18
!--------------------*********--------------------------------------+-------

   END MODULE cubint_module
!-------------*************-----------------------------------------+-------
