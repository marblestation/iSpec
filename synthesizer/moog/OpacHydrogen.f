c******************************************************************************
c  The subroutines needed to calculate the H I b-f, H I f-f, H- b-f, and
c  H- f-f opacities are in this file.  These are from ATLAS9.
c******************************************************************************





      subroutine opacH1
c******************************************************************************
c  This routine computes the H- bound-free and free-free opacities.
c  It returns the cross-section times the number density of H- particles.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com' 
      include 'Kappa.com'
      include 'Mol.com'
      real*8 cont(8), bolt(100,8), exlim(100), freet(100), boltex(100)
      save
      data modcount/0/

c  set up some data upon first entrance with a new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            do n=1,8
               xn2 = dfloat(n*n)
               bolt(i,n) = dexp(-13.595*(1.-1./xn2)/tkev(i))*2.*
     .                     xn2*numdens(1,1,i)/u(1,1,i)
            enddo
            freet(i) = ne(i)*numdens(1,2,i)/u(1,2,i)/dsqrt(t(i))
            xr = numdens(1,1,i)/u(1,1,i)*(2./2./13.595)*tkev(i)
            boltex(i) = dexp(-13.427/tkev(i))*xr
            exlim(i) = dexp(-13.595/tkev(i))*xr
         enddo
      endif

      do n=1,8
         cont(n) = coulx(n,freq,1.d0)
      enddo

      freq3 = freq**3
      cfree = 3.6919d8/freq3
      c = 2.815d29/freq3
      do i=1,ntau
         ex = boltex(i)
         if (freq .lt. 4.05933d13) ex = exlim(i)/evhkt(i)
         h = (cont(7)*bolt(i,7) + cont8*bolt(i,8) + 
     .       (ex - exlim(i))*c + 
     .       coulff(1,tlog(i),freq)*freet(i)*cfree)*(1.-evhkt(i))
         do n=1,6
            h = h + cont(n)*bolt(i,n)*(1.-evhkt(i))
         enddo
         aH1(i) = h
      enddo

      return
      end






      subroutine opacHminus
c******************************************************************************
c  This routine computes the H- bound-free and free-free opacities.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com' 
      include 'Kappa.com'
      real*8 wbf(85), bf(85), fflog(22,11), ff(11,22)
      real*8 ffbeg(11,11), ffend(11,11), fftt(11), wfflog(22)
      real*8 fftheta(100), thetaff(11), wavek(22)
      real*8 xhmin(100)
      equivalence (ff(1,1),ffbeg(1,1)), (ff(1,12),ffend(1,1))
      save
c  From Mathisen (1984), after Wishart (1979) and Broad & Reinhardt (1976) 
      data wbf/  18.00,  19.60,  21.40,  23.60,  26.40,  29.80,  34.30,
     .   40.40,  49.10,  62.60, 111.30, 112.10, 112.67, 112.95, 113.05,
     .  113.10, 113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00,
     .  175.00, 200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00,
     .  375.00, 400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00,
     .  575.00, 600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00,
     .  775.00, 800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00,
     .  975.00,1000.00,1025.00,1050.00,1075.00,1100.00,1125.00,1150.00,
     . 1175.00,1200.00,1225.00,1250.00,1275.00,1300.00,1325.00,1350.00,
     . 1375.00,1400.00,1425.00,1450.00,1475.00,1500.00,1525.00,1550.00,
     . 1575.00,1600.00,1610.00,1620.00,1630.00,1643.91/
      data bf/   0.067,  0.088,  0.117,  0.155,  0.206,  0.283,  0.414,
     .   0.703,   1.24,   2.33,  11.60,  13.90,  24.30,  66.70,  95.00,
     .   56.60,  20.00,  14.60,   8.50,   7.10,   5.43,   5.91,   7.29,
     .   7.918,  9.453,  11.08,  12.75,  14.46,  16.19,  17.92,  19.65,
     .   21.35,  23.02,  24.65,  26.24,  27.77,  29.23,  30.62,  31.94,
     .   33.17,  34.32,  35.37,  36.32,  37.17,  37.91,  38.54,  39.07,
     .   39.48,  39.77,  39.95,  40.01,  39.95,  39.77,  39.48,  39.06,
     .   38.53,  37.89,  37.13,  36.25,  35.28,  34.19,  33.01,  31.72,
     .   30.34,  28.87,  27.33,  25.71,  24.02,  22.26,  20.46,  18.62,
     .   16.74,  14.85,  12.95,  11.07,  9.211,  7.407,  5.677,  4.052,
     .   2.575,  1.302, 0.8697, 0.4974, 0.1989,    0. /
C  Bell and Berrington (1987, J.Phys.B, 20, 801-806)
      data wavek/ .50,.40,.35,.30,.25,.20,.18,.16,.14,.12,.10,.09,.08,
     .            .07,.06,.05,.04,.03,.02,.01,.008,.006/
      data thetaff/
     .  0.5,  0.6, 0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.8,  3.6/
      data ffbeg/
     ..0178,.0222,.0308,.0402,.0498,.0596,.0695,.0795,.0896, .131, .172,   1823
     ..0228,.0280,.0388,.0499,.0614,.0732,.0851,.0972, .110, .160, .211,   2278
     ..0277,.0342,.0476,.0615,.0760,.0908, .105, .121, .136, .199, .262,   2604
     ..0364,.0447,.0616,.0789,.0966, .114, .132, .150, .169, .243, .318,   3038
     ..0520,.0633,.0859, .108, .131, .154, .178, .201, .225, .321, .418,   3645
     ..0791,.0959, .129, .161, .194, .227, .260, .293, .327, .463, .602,   4557
     ..0965, .117, .157, .195, .234, .272, .311, .351, .390, .549, .711,   5063
     . .121, .146, .195, .241, .288, .334, .381, .428, .475, .667, .861,   5696
     . .154, .188, .249, .309, .367, .424, .482, .539, .597, .830, 1.07,   6510
     . .208, .250, .332, .409, .484, .557, .630, .702, .774, 1.06, 1.36,   7595
     . .293, .354, .468, .576, .677, .777, .874, .969, 1.06, 1.45, 1.83/   9113
      data ffend/
     . .358, .432, .572, .702, .825, .943, 1.06, 1.17, 1.28, 1.73, 2.17,  10126
     . .448, .539, .711, .871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.60,  11392
     . .579, .699, .924, 1.13, 1.33, 1.51, 1.69, 1.86, 2.02, 2.67, 3.31,  13019
     . .781, .940, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48, 2.69, 3.52, 4.31,  15189
     . 1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.20, 3.51, 3.80, 4.92, 5.97,  18227
     . 1.73, 2.08, 2.74, 3.37, 3.90, 4.50, 5.01, 5.50, 5.95, 7.59, 9.06,  22784
     . 3.04, 3.65, 4.80, 5.86, 6.86, 7.79, 8.67, 9.50, 10.3, 13.2, 15.6,  30378
     . 6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2, 23.0, 29.5, 35.0,  45567
     . 27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117., 140.,  91134
     . 42.3, 50.6, 66.4, 80.8, 94.5, 107., 120., 131., 142., 183., 219., 113918
     . 75.1, 90.0, 118., 144., 168., 191., 212., 234., 253., 325., 388./ 151890
      data (xhmin(i),i=1,100)/100*0.0/
      data modcount,istart/ 0,0/
         
c  fill some arrays once and for all
      if (istart .eq. 0) then
         istart = 1
c  91.134 number taken from Bell & Berrington
         do iwave=1,22
            wfflog(iwave) = dlog(91.134d0/wavek(iwave))
            do itheta=1,11
               fflog(iwave,itheta) = dlog(ff(itheta,iwave)*1.d-26)
            enddo
         enddo
      endif

c  initialize some quantities for each new model atmosphere
c  .754209 Hotop & Lineberger (1985, J. Phys. Chem. Ref. Data, 14,731-752)
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            xhmin(i) = dexp(.754209/tkev(i))/(2.*2.4148d15*t(i)**1.5)*
     .                 numdens(1,1,i)/u(1,1,i)*ne(i)
         enddo
      endif

c  main opacity computation yielding "aHminus"
      wave = 2.99792458d17/freq
      wavelog = dlog(wave)
      do itheta=1,11
         nnnn = 22
         call linter (wfflog,fflog(1,itheta),nnnn,wavelog,fftlog,1)
         fftt(itheta) = dexp(fftlog)/thetaff(itheta)*5040.*1.380658E-16
      enddo
      hminbf = 0.
      nnnn = 85
      if (freq  .gt.  1.82365d14) 
     .              maxwave = map1(wbf,bf,nnnn,wave,hminbf,1) 
      do i=1,ntau
         nnnn = 11
         call linter (thetaff,fftt,nnnn,theta(i),fftheta(i),1)
         hminff = fftheta(i)*numdens(1,1,i)/u(1,1,i)*2.*ne(i)
         h = hminbf*1.d-18*(1.-evhkt(i))*xhmin(i)
         aHminus(i) = h + hminff
      enddo

      return
      end





      subroutine linter (xold,yold,nold,xnew,ynew,nnew)
c******************************************************************************
c  this is a linear interpolation scheme.  The arrays "xold" and "xnew" 
c  should be increasing order.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      real*8 xold(*), yold(*), xnew(*), ynew(*)
        
      iold = 2
      do inew=1,nnew
1       if (xnew(inew).lt.xold(iold) .or. iold.eq.nold) then
           ynew(inew) = yold(iold-1) + (yold(iold)-yold(iold-1))/
     .              (xold(iold)-xold(iold-1))*(xnew(inew)-xold(iold-1))
           return
        else
           iold = iold + 1
           go to 1
        endif
      enddo
      
      end






      integer function map1 (xold,fold,nold,xnew,fnew,nnew)
c******************************************************************************
c  This is an interpolation scheme that is copied blindly from ATLAS. 
c  I decided that it would be too dangerous to re-write this one.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      real*8 xold(*), fold(*), xnew(*), fnew(*)

      l = 2
      ll = 0
      do 50 k=1,nnew
10       if(xnew(k) .lt. xold(l)) go to 20
         l = l + 1
         if (l .gt. nold) go to 30
         go to 10
20       if (l .eq. ll) go to 50
         if (l .eq. 2) go to 30
         if (l .eq. 3) go to 30
         l1 = l - 1
         if (l.gt.ll+1 .or. l.eq.3) go to 21
         if (l.gt.ll+1 .or. l.eq.4) go to 21
         cbac = cfor
         bbac = bfor
         abac = afor
         if (l .eq. nold) go to 22
         go to 25
21       l2 = l - 2
         d = (fold(l1)-fold(l2))/(xold(l1)-xold(l2))
         cbac = fold(l)/((xold(l)-xold(l1))*(xold(l)-xold(l2)))+
     .      (fold(l2)/(xold(l)-xold(l2))-fold(l1)/(xold(l)-xold(l1)))/
     .      (xold(l1)-xold(l2))
         bbac = d - (xold(l1)+xold(l2))*cbac
         abac = fold(l2) - xold(l2)*d + xold(l1)*xold(l2)*cbac
         if (l .lt. nold) go to 25
22       c = cbac
         b = bbac
         a = abac
         ll = l
         go to 50
25       d = (fold(l)-fold(l1))/(xold(l)-xold(l1))
         cfor =  fold(l+1)/((xold(l+1)-xold(l))*(xold(l+1)-xold(l1)))+
     .      (fold(l1)/(xold(l+1)-xold(l1))-fold(l)/(xold(l+1)-xold(l)))/
     .      (xold(L)-xold(l1))
         bfor = d - (xold(l)+xold(l1))*cfor
         afor = fold(l1) - xold(l1)*d + xold(l)*xold(l1)*cfor
         wt = 0.
         if (dabs(cfor) .ne. 0.) wt = dabs(cfor)/(dabs(cfor)+dabs(cbac))
         a = afor + wt*(abac-afor)
         b = bfor + wt*(bbac-bfor)
         c = cfor + wt*(cbac-cfor)
         ll = l
         go to 50
30       if (l .eq. ll) go to 50
         l = min0(nold,l)
         c = 0.
         b = (fold(l)-fold(l-1))/(xold(l)-xold(l-1))
         a = fold(l) - xold(l)*b
         ll = l
50       fnew(k)= a + (b + c*xnew(k))*xnew(k)

         map1 = ll - 1 
         return
         end




