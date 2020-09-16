
c******************************************************************************
c  The subroutines needed to calculate the Mg I, Mg II, Al I, Si I, Si II,
c  and Fe I b-f opacities are in this file.  These are from ATLAS9.
c******************************************************************************





      subroutine opacC1
c******************************************************************************
c     This routine computes the bound-free absorption due to C I.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      real*8 c1240(100), c1444(100)
      save
      data freq1, modcount/0. ,0/

c  initialize some quantities for each new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            c1240(i) = 5.*dexp(-1.264/tkev(i))
            c1444(i) = dexp(-2.683/tkev(i))
         enddo
      endif

c  initialize some quantities for each new model atmosphere or new frequency;
c  Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
c  Burke, P.G. and Taylor, K.T. 1979, J. Phys. B, 12, 2971-2984.
      if (modelnum.ne.modcount .or. freq.ne.freq1) then
         freq1 = freq
         ryd = 109732.298
         waveno = freq/2.99792458d10
         xs0 = 0.
         xs1 = 0.
         xd0 = 0.
         xd1 = 0.
         xd2 = 0.
         x1444 = 0.
         x1240 = 0.
         x1100 = 0.
c        P2 3P   1
c  I AM NOT SURE WHETHER THE CALL TO SEATON IN THE NEXT STATEMENT IS
c  CORRECT, BUT IT ONLY AFFECTS THINGS BELOW 1100A
         if (freq .ge. 2.7254d15) x1100 =
     .      10.**(-16.80-(waveno-90777.000)/3./ryd)*
     .      seaton (2.7254d15,1.219d-17,2.d0,3.317d0)
c        P2 1D   2
         if (freq .ge. 2.4196d15) then
            xd0 = 10.**(-16.80-(waveno-80627.760)/3./ryd)
            eeps = (waveno-93917.)*2./9230.
            aa = 22.d-18
            bb = 26.d-18
            xd1 = (aa*eeps+bb)/(eeps**2+1.)
            eeps = (waveno-111130.)*2./2743.
            aa = -10.5d-18
            bb = 46.d-18
            xd2 = (aa*eeps+bb)/(eeps**2+1.)
            x1240 = xd0 + xd1 + xd2
         endif
c        P2 1S   3
         if (freq .ge. 2.0761d15) then
            xs0 = 10.**(-16.80-(waveno-69172.400)/3./ryd)
            eeps = (waveno-97700.)*2./2743.
            aa = 68.d-18
            bb = 118.d-18
            xs1 = (aa*eeps+bb)/(eeps**2+1.)
            x1444 = xs0 + xs1
         endif
      endif

      do i=1,ntau
         if (freq .ge. 2.0761d15) then
            aC1(i) = (x1100*9. + x1240*c1240(i) + x1444*c1444(i))*
     .               numdens(3,1,i)/u(6,1,i)
         endif
      enddo

      return
      end








      subroutine opacMg1
c******************************************************************************
c     This routine computes the bound-free absorption due to Mg I.  
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      real*8 flog(9), freqMg(7), peach(7,15), xx(7), tlg(7)
      real*8 dt(100)
      integer nt(100)
      save
      data peach/
c         4000     5000     6000     7000     8000     9000    10000 
     . -42.474, -42.350, -42.109, -41.795, -41.467, -41.159, -40.883,
     . -41.808, -41.735, -41.582, -41.363, -41.115, -40.866, -40.631,
     . -41.273, -41.223, -41.114, -40.951, -40.755, -40.549, -40.347,
     . -45.583, -44.008, -42.957, -42.205, -41.639, -41.198, -40.841,
     . -44.324, -42.747, -41.694, -40.939, -40.370, -39.925, -39.566,
     . -50.969, -48.388, -46.630, -45.344, -44.355, -43.568, -42.924,
     . -50.633, -48.026, -46.220, -44.859, -43.803, -42.957, -42.264,
     . -53.028, -49.643, -47.367, -45.729, -44.491, -43.520, -42.736,
     . -51.785, -48.352, -46.050, -44.393, -43.140, -42.157, -41.363,
     . -52.285, -48.797, -46.453, -44.765, -43.486, -42.480, -41.668,
     . -52.028, -48.540, -46.196, -44.507, -43.227, -42.222, -41.408,
     . -52.384, -48.876, -46.513, -44.806, -43.509, -42.488, -41.660,
     . -52.363, -48.856, -46.493, -44.786, -43.489, -42.467, -41.639,
     . -54.704, -50.772, -48.107, -46.176, -44.707, -43.549, -42.611,
     . -54.359, -50.349, -47.643, -45.685, -44.198, -43.027, -42.418/
      data freqMg/ 1.9341452d15, 1.8488510d15, 1.1925797d15,
     .             7.9804046d14, 4.5772110d14, 4.1440977d14,
     .             4.1113514d14/
      data flog/ 35.23123, 35.19844, 35.15334, 34.71490, 34.31318,
     .           33.75728, 33.65788, 33.64994, 33.43947/
      data tlg/ 8.29405, 8.51719, 8.69951, 8.85367, 8.98720, 9.10498, 
     .          9.21034/
      data freq1, modcount/0. ,0/

c  initialize some quantities for each new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            n = max0(min0(6,nint(t(i)/1000.)-3),1)
            nt(i) = n
            dt(i) = (tlog(i)-tlg(n))/(tlg(n+1)-tlg(n))
         enddo
      endif
   
c  initialize some quantities for each new model atmosphere or new frequency;
      if (modelnum.ne.modcount .or. freq.ne.freq1) then
         freq1 = freq
         do n=1,7
            if (freq .gt. freqMg(n)) go to 23
         enddo
         n = 8
23       dd = (freqlg-flog(n))/(flog(n+1)-flog(n))
         if (n .gt. 2) n = 2*n -2
         dd1 = 1.0 - dd
         do it=1,7
            xx(it) = peach(it,n+1)*dd + peach(it,n)*dd1
         enddo
      endif

      do i=1,ntau
         if (freq .ge. 2.997925d+14) then
            n = nt(i)
            aMg1(i) = dexp(xx(n)*(1.d0-dt(i))+xx(n+1)*dt(i))*
     .                numdens(4,1,i)/u(12,1,i)
         endif
      enddo

      return
      end







      subroutine opacMg2
c******************************************************************************
c     This routine computes the bound-free absorption due to Mg II.  
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      real*8 c1169(100)
      save
      data freq1, modcount/0. ,0/

c  initialize some quantities for each new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            c1169(i) = 6.*dexp(-4.43d+0/tkev(i))
         enddo
      endif

c  initialize some quantities for each new model atmosphere or new frequency;
c  there are two edges, one at 824 A and the other at 1169 A
      if (modelnum.ne.modcount .or. freq.ne.freq1) then
         freq1 = freq
         if (freq .ge. 3.635492d15) then
            x824 = seaton (3.635492d15,1.40d-19,4.d0,6.7d0)
         else
            x824 = 1.d-99
         endif
         if (freq .ge. 2.564306d15) then
            x1169 = 5.11d-19*(2.564306d15/freq)**3
         else
            x1169 = 1.d-99
         endif
      endif
      
      do i=1,ntau
         if (x1169 .ge. 1.d-90) then
            aMg2(i) = (x824*2. + x1169*c1169(i))*
     .                numdens(4,2,i)/u(12,2,i)
         endif
      enddo

      return
      end








      subroutine opacAl1
c******************************************************************************
c     This routine computes the bound-free absorption due to Al I.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
   
      do i=1,ntau
         if (freq .ge. 1.443d15) then
            aAl1(i) = 6.5d-17*(1.443d15/freq)**5*6.*
     .                numdens(5,1,i)/u(13,1,i)
         endif
      enddo

      return
      end









      subroutine opacSi1
c******************************************************************************
c     This routine computes the bound-free absorption due to Si I.  
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      real*8 flog(11), freqSi(9), xx(9), tlg(9), dt(100) 
      integer nt(100)
      real*8 peach(9,19)
      save
      data peach/                                                     
c       4000   5000   6000   7000   8000   9000   10000  11000  12000 
     . 38.136,38.138,38.140,38.141,38.143,38.144,38.144,38.145,38.145,
     . 37.834,37.839,37.843,37.847,37.850,37.853,37.855,37.857,37.858,
     . 37.898,37.898,37.897,37.897,37.897,37.896,37.895,37.895,37.894,
     . 40.737,40.319,40.047,39.855,39.714,39.604,39.517,39.445,39.385,
     . 40.581,40.164,39.893,39.702,39.561,39.452,39.366,39.295,39.235,
     . 45.521,44.456,43.753,43.254,42.878,42.580,42.332,42.119,41.930,
     . 45.520,44.455,43.752,43.251,42.871,42.569,42.315,42.094,41.896,
     . 55.068,51.783,49.553,47.942,46.723,45.768,44.997,44.360,43.823,
     . 53.868,50.369,48.031,46.355,45.092,44.104,43.308,42.652,42.100,
     . 54.133,50.597,48.233,46.539,45.261,44.262,43.456,42.790,42.230,
     . 54.051,50.514,48.150,46.454,45.176,44.175,43.368,42.702,42.141,
     . 54.442,50.854,48.455,46.733,45.433,44.415,43.592,42.912,42.340,
     . 54.320,50.722,48.313,46.583,45.277,44.251,43.423,42.738,42.160,
     . 55.691,51.965,49.444,47.615,46.221,45.119,44.223,43.478,42.848,
     . 55.661,51.933,49.412,47.582,46.188,45.085,44.189,43.445,42.813,
     . 55.973,52.193,49.630,47.769,46.349,45.226,44.314,43.555,42.913,
     . 55.922,52.141,49.577,47.715,46.295,45.172,44.259,43.500,42.858,
     . 56.828,52.821,50.110,48.146,46.654,45.477,44.522,43.730,43.061,
     . 56.657,52.653,49.944,47.983,46.491,45.315,44.360,43.569,42.901/
c     3P, 1D, 1S, 1D, 3D, 3F, 1D, 3P
      data freqSi/ 2.1413750d15, 1.9723165d15, 1.7879689d15,
     .             1.5152920d15, 5.5723927d14, 5.3295914d14,
     .             4.7886458d14, 4.7216422d14, 4.6185133d14/
      data flog/ 35.45438, 35.30022, 35.21799, 35.11986, 34.95438,
     .           33.95402, 33.90947, 33.80244, 33.78835, 33.76626,
     .           33.70518/
      data tlg/ 8.29405, 8.51719, 8.69951, 8.85367, 8.98720, 9.10498,
     .          9.21034, 9.30565, 9.39266/
      data freq1, modcount/0. ,0/

c  initialize some quantities for each new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            n = max0(min0(8,nint(t(i)/1000.)-3),1)
            nt(i) = n
            dt(i) = (tlog(i)-tlg(n))/(tlg(n+1)-tlg(n))
         enddo
      endif

c  initialize some quantities for each new model atmosphere or new frequency
      if (modelnum.ne.modcount .or. freq.ne.freq1) then
         freq1 = freq
         do n=1,9
            if (freq .gt. freqSi(n)) go to 23
         enddo
         n = 10
23       dd = (freqlg-flog(n))/(flog(n+1)-flog(n))
         if (n .gt. 2) n = 2*n - 2
         dd1 = 1.0 - dd
         do it=1,9
            xx(it) = peach(it,n+1)*dd + peach(it,n)*dd1
         enddo
      endif

      do i=1,ntau
         if (freq .ge. 2.997925d+14) then
            n = nt(i)
            aSi1(i) = (dexp(-(xx(n)*(1.-dt(i)) + xx(n+1)*dt(i)))*9.)*
     .                numdens(6,1,i)/u(14,1,i)
         endif
      enddo

      return
      end






      subroutine opacSi2
c******************************************************************************
c     This routine computes the bound-free absorption due to Si II.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      real*8 flog(9), freqSi(7), peach(6,14), xx(6), tlg(6), dt(100)
      integer nt(100)
      save
      data peach/
c         10000     12000     14000     16000     18000     20000  
     .  -43.8941, -43.8941, -43.8941, -43.8941, -43.8941, -43.8941,
     .  -42.2444, -42.2444, -42.2444, -42.2444, -42.2444, -42.2444,
     .  -40.6054, -40.6054, -40.6054, -40.6054, -40.6054, -40.6054,
     .  -54.2389, -52.2906, -50.8799, -49.8033, -48.9485, -48.2490,
     .  -50.4108, -48.4892, -47.1090, -46.0672, -45.2510, -44.5933,
     .  -52.0936, -50.0741, -48.5999, -47.4676, -46.5649, -45.8246,
     .  -51.9548, -49.9371, -48.4647, -47.3340, -46.4333, -45.6947,
     .  -54.2407, -51.7319, -49.9178, -48.5395, -47.4529, -46.5709,
     .  -52.7355, -50.2218, -48.4059, -47.0267, -45.9402, -45.0592,
     .  -53.5387, -50.9189, -49.0200, -47.5750, -46.4341, -45.5082,
     .  -53.2417, -50.6234, -48.7252, -47.2810, -46.1410, -45.2153,
     .  -53.5097, -50.8535, -48.9263, -47.4586, -46.2994, -45.3581,
     .  -54.0561, -51.2365, -49.1980, -47.6497, -46.4302, -45.4414,
     .  -53.8469, -51.0256, -48.9860, -47.4368, -46.2162, -45.2266/
      data freqSi/ 4.9965417d15, 3.9466738d15, 1.5736321d15,
     .             1.5171539d15, 9.2378947d14, 8.3825004d14,
     .             7.6869872d14/
c     2P, 2D, 2P, 2D, 2P
      data flog/ 36.32984, 36.14752, 35.91165, 34.99216, 34.95561,
     .           34.45951, 34.36234, 34.27572, 34.20161/
      data tlg/ 9.21034, 9.39266, 9.54681, 9.68034, 9.79813, 9.90349/
      data freq1, modcount/0., 0/

c  set up some data upon first entrance with a new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            n = max0(min0(5,nint(t(i)/2000.)-4),1)
            nt(i) = n
            dt(i) = (tlog(i)-tlg(n))/(tlg(n+1)-tlg(n))
         enddo
      endif

c  initialize some quantities for each new model atmosphere or new frequency
      if (modelnum.ne.modcount .or. freq.ne.freq1) then
         freq1 = freq
         do n=1,7
            if (freq .gt. freqSi(n)) go to 23
         enddo
         n = 8
23       dd = (freqlg-flog(n))/(flog(n+1)-flog(n))
         if (n .gt. 2) n = 2*n - 2
         if (n .eq. 14) n = 13
         dd1 = 1.0 - dd
         do it=1,6
            xx(it) = peach(it,n+1)*dd + peach(it,n)*dd1
         enddo
      endif

      do i=1,ntau
         if (freq .ge. 7.6869872d14) then
            n = nt(i)
            aSi2(i) = (dexp(xx(n)*(1.-dt(i)) + xx(n+1)*dt(i))*6.)*
     .             numdens(6,2,i)/u(14,2,i)
         endif
      enddo

      return
      end







      subroutine opacFe1
c******************************************************************************
c     This routine computes the bound-free absorption due to Fe I.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Kappa.com'
      real*8 bolt(48,100), gg(48), ee(48), wno(48), xsect(48)
      save
      data gg/25.,35.,21.,15.,9.,35.,33.,21.,27.,49.,9.,21.,27.,9.,9.,
     . 25.,33.,15.,35.,3.,5.,11.,15.,13.,15.,9.,21.,15.,21.,25.,35.,
     . 9.,5.,45.,27.,21.,15.,21.,15.,25.,21.,35.,5.,15.,45.,35.,55.,25./
      data ee/500.,7500.,12500.,17500.,19000.,19500.,19500.,21000.,
     . 22000.,23000.,23000.,24000.,24000.,24500.,24500.,26000.,26500.,
     . 26500.,27000.,27500.,28500.,29000.,29500.,29500.,29500.,30000.,
     . 31500.,31500.,33500.,33500.,34000.,34500.,34500.,35000.,35500.,
     . 37000.,37000.,37000.,38500.,40000.,40000.,41000.,41000.,43000.,
     . 43000.,43000.,43000.,44000./
      data wno/63500.,58500.,53500.,59500.,45000.,44500.,44500.,43000.,
     . 58000.,41000.,54000.,40000.,40000.,57500.,55500.,38000.,57500.,
     . 57500.,37000.,54500.,53500.,55000.,34500.,34500.,34500.,34000.,
     . 32500.,32500.,32500.,32500.,32000.,29500.,29500.,31000.,30500.,
     . 29000.,27000.,54000.,27500.,24000.,47000.,23000.,44000.,42000.,
     . 42000.,21000.,42000.,42000./
      data freq1, modcount/0., 0/

c  set up some data upon first entrance with a new model atmosphere
      if (modelnum .ne. modcount) then
         modcount = modelnum
         do i=1,ntau
            hkt = 6.6256d-27/(1.38065d-16*t(i))
            do k=1,48
               bolt(k,i) = gg(k)*dexp(-ee(k)*2.99792458d10*hkt)
            enddo
         enddo
      endif

c  initialize some quantities for each new model atmosphere or new frequency;
c  the absorption begins at 4762 A.
      if (modelnum.ne.modcount .or. freq.ne.freq1) then
         freq1 = freq
         waveno = freq/2.99792458d10
         if (waveno .ge. 21000.) then
            do k=1,48
               xsect(k) = 0.
               if (wno(k) .lt. waveno) xsect(k)= 3.d-18/
     .                      (1.+((wno(k)+3000.-waveno)/wno(k)/.1)**4)
            enddo
         endif
      endif

      do i=1,ntau
         if (waveno .ge. 21000.) then
            do k=1,48
               aFe1(i) = aFe1(i) + xsect(k)*bolt(k,i)*
     .             numdens(7,1,i)/u(26,1,i)
            enddo
         endif
            
      enddo

      return
      end







      real*8 function seaton (freq0,xsect,power,a)
c******************************************************************************
c     This function is a general representation for b-f absorption above
c     a given ionization limit freq0, using cross-section xsect, 
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Kappa.com'
      freqratio = freq0/freq
      seaton = xsect*(a + freqratio*(1.0-a))*
     .         dsqrt(freqratio**(nint(2.*power+0.01)))

      return
      end


