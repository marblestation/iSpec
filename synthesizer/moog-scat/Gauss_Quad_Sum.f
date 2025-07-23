      subroutine gauss_quad_sum (k,a,b,ai,xmyi)
c******************************************************************************
c     For integration between the upper bound B and the lower bound A, this 
c     Routine gives Gaussian weights and integration points. 
 
c     Source for data: Lowan, Davids, Levenson,  Bull Amer Math Soc  48 page 739  (1942)
c     Source for original code: Olof Morell, Uppsala (03-1988)

c     ai/wtmu = weights, xmyi/mu = integration points

c     Until code below is altered, AngWeight works only with odd numbers of rays.
c*****************************************************************************

      implicit real*8 (a-h,o-z)
      real*8 ai(5), xmyi(5), ap(29), xmyp(29)
      integer kud, ioev, ined, indov(9)
      data ap/1.0,0.55555555555555,.88888888888888,.347854845137,
     &        0.65214515486254,0.23692688505618,0.47862867049936,
     &        0.56888888888888,0.17132449237917,0.36076157304813,
     &        0.46791393457269,0.12948496616887,0.27970539148927,
     &        0.38183005050511,0.41795918367346,0.10122853629037,
     &        0.22238103445337,0.31370664587788,0.36268378337836,
     &        0.08127438836157,0.18064816069485,0.26061069640293,
     &        0.31234707704000,0.33023935500126,0.06667134430868,
     &        0.14945134915058,0.21908636251598,0.26926671930999,
     &        0.29552422471475/
      data xmyp/
     &     0.57735026918962,.77459666924148,.0,0.86113631159405,
     &     0.33998104358485,.90617984593866,.53846931010568,.0,
     &     0.93246951420315,.66120938646626,.23861918608319,
     &     0.94910791234275,.74153118559939,.40584515137739,.0,
     &     0.96028985649753,.79666647741362,.52553240991632,
     &     0.18343464249565,.96816023950762,.83603110732663,
     &     0.61337143270059,.32425342340380,.0,0.97390652851717,
     &     0.86506336668898,.67940956829902,.43339539412924,
     &     0.14887433898163/
      data indov/1,3,5,8,11,15,19,24,29/
 

      if (k .eq. 1) then
         xmyi(1) = (b + a)/2.
         ai(1) = b - a
         return
      endif

      kud = 0
      flk = dfloat(k)/2.
      k2 = k/2
      fk = dfloat(k2)
      if (dabs(flk-fk) .ge. -1.0d-7) then
         k2 = k2 + 1
         kud = 1
      endif
      ioev = indov(k-1)
      ined = ioev - k2
      do i=1,k2
         ip = ined + i
         xmyi(i) = (-xmyp(ip)*(b-a) + (b+a))/2.
         ai(i) = (b-a)*ap(ip)/2.
      enddo
      k2 = k2 + 1
      do i=k2,k
        ip = ioev + k2 - i
        if (kud .gt. 0) ip = ip - 1
        xmyi(i)= (xmyp(ip)*(b-a) + (b+a))/2.
        ai(i) = (b-a)*ap(ip)/2.
      enddo


      return
      end


