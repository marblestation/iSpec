
      subroutine eqlib                                                  
c******************************************************************************
c     This routine solves the equations of molecular dissociative           
c     equilibrium. the equations can include neutral and ionized molecules  
c     and atoms 
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Mol.com'
      include 'Quants.com'
      include 'Dummy.com'
      real*8 xfic(30), xcorr(30), deltax(30), c(30,30), ans(30,30)      
      real*8 uu(2)
      real*8 lth
      integer ident(30,110)


c*****clear the arrays 
      iorder(1) = 0                                                
      neq = iorder(1) 
      tdel = (t(ntau)-t(1))/4.0 + 1.0                                   
      do k=1,30                                                       
         iorder(k) = 0                                                  
         xcorr(k) = 0.                                                  
         do kk=1,110                       
            ident(k,kk) = 0                                             
         enddo
      enddo


c*****the number of species to be considered has been defined in "Inmodel";
c*****either read in the dissociation data for a molecular species
      do jmol=1,nmol                                                   
         if (amol(jmol) .ge. 100.) then
            do k=1,110                                                  
               if (datmol(1,k) .eq. amol(jmol)) go to 11                
            enddo
            write (nf1out,1001) amol(jmol)                   
            stop                                                        
11          do kk=1,6                                                   
               const(kk,jmol) = datmol(kk+1,k)                          
            enddo 


c*****or read the ionization data for an atomic species
         else
            iatom1 = amol(jmol) + 0.0001
            atom = iatom1                                            
            const(1,jmol) = xchi1(iatom1)                               
            do kk=1,5               
               ti = t(1) + (kk-1)*tdel 
               it = ti
               do jj=1,2
                  att = atom+0.1*(jj-1)
                  if (partflag(iatom1,jj) .gt. 0) then
                     uu(jj) = partnew(att,jj,it)
                  else
                     uu(jj) = ucalc(att,it)
                  endif
               enddo
               const(kk+1,jmol) = uu(2)/uu(1)
            enddo 
         endif
      enddo
      if (molopt .ge. 2) 
     .   write (nf1out,1002) (amol(i),(const(j,i),j=1,6),i=1,nmol)


c*****set up the molecular equilibrium array. each element of array         
c*****'ident' contains the identifier (subscript of an element of           
c*****array 'amol') of an ion or a molecule. the neutral atom is always     
c*****understood, but is not explicitly contained in 'ident'.               
      nmax = 1                                                          
      do jmol=1,nmol                                                  
         atom = amol(jmol)                                              
2        call sunder(atom,iatom1,iatom2)                                
         do k=1,30                                                      
            if (iatom1 .eq. iorder(k)) go to 4                          
         enddo
         neq = neq + 1                                                  
         iorder(neq) = iatom1                                           
         ident(neq,1) = jmol                                            
         go to 6                                                        
4        do kk=1,nmax                                                   
            if (ident(k,kk).eq.0 .or. ident(k,kk).eq. jmol) go to 7     
         enddo 
         nmax = nmax + 1                                                
         kk = nmax                                                
7        ident(k,kk) = jmol                                             
6        if (iatom2 .ne. 0) then
            atom = iatom2                                               
            go to 2                                                     
         endif
      enddo                                                           
      if (molopt .ge. 2) then
         do i=1,neq
            dummy1(i) = iorder(i)
         enddo
         write (nf1out,1003) (dummy1(i),i=1,neq),(amol(i),i=1,nmol)
      endif


c*****now begin the loop that goes through all the atmosphere tau layers
      do 21 kev=1,ntau                                                  


c*****calculate *xfic* and make a first guess at *xatom*                    
      i = ntau + 1 - kev                                           
      lev = i
      tk = 1.38065d-16*t(i)                                             
      do k=1,neq                                                      
         korder = iorder(k)                                             
         xfic(k) = xabund(korder)*nhtot(i)                              
      enddo
      if (i .lt. ntau) then
         do k=1,neq                                                     
            xatom(k) = xatom(k)*nhtot(i)/nhtot(i+1)                     
         enddo
      else
         do k=1,neq                                                     
            xatom(k) = xfic(k)                                          
         enddo
      endif


c*****compute the number of molecules:
c*****Here is some information about the equilibrium constants.
c*****The equilibrium constants, Kp, are defined as follows:
c
c
c                 P(A)*P(B)     N(A)*N(B)
c       Kp(AB)  = --------- = kT--------- =
c                  P(AB)         N(AB)
c
c             2*pi*kT 3/2    M(A)*M(B) 3/2   Q(A)*Q(B)
c       = kT*(-------)    * (---------)    * --------- * exp(-D(AB)/kT)
c               h^2            M(AB)           Q(AB)
c
c 
c*****MOOG uses: n(AB)/n(A)n(B) = kt/Kp
c 
c        Kp - dissociation constant, Q - partition functions, M - masses
c        P - partial pressures, N - number densities, T - temperature,
c        D - dissociation energy, h - plank constant. Remember to use
c        masses in grams (1 amu = 1.660540E-24 g) and energy in ergs
c        (1 eV = 1.60219E-12 ergs). Also, k = 1.38065E-16 erg/K,
c        h = 6.626076E-27 erg s, and pi = 3.1415926536.
27    do jmol=1,nmol                                                  
         atom = amol(jmol)                                              
         if (atom .ge. 100.) then
            if (t(i) .gt. 12000.) then
               xmol(jmol,i) = 1.0d-20
            else   
               xmol(jmol,i) = 1.
               count = 0.      
37             call sunder(atom,iatom1,iatom2)
               count = count + 1.            
               do k=1,neq                   
                  if (iorder(k) .eq. iatom1) 
     .               xmol(jmol,i) = xmol(jmol,i)*xatom(k)
               enddo                                    
               if (iatom2 .ne. 0) then
                  atom = iatom2                        
                  go to 37                                         
               endif
               hion = 10.0*(amol(jmol) - dint(amol(jmol)))  
               th = 5040./t(i)
               lth = log10(th)
               xmol(jmol,i) = xmol(jmol,i)*(((1.38065d-16 * t(i))**
     .             (count-1.0))/(10.0**(const(2,jmol)+(const(3,jmol)*
     .             lth)+(const(4,jmol)*(lth**2))+(const(5,jmol)*
     .             (lth**3))+(const(6,jmol)*(lth**4))-
     .             (const(1,jmol)*th))))/(ne(i)**hion)
            endif


c*****compute the number of ions:
         else
            delt = (t(i)-t(1))/tdel                               
            m = min0(idint(delt)+2,5)                          
            delt = delt - idint(delt)                           
             u1 = const(m,jmol) + 
     .           (const(m+1,jmol)-const(m,jmol))*delt          
            iatom1 = atom                                           
            do k=1,neq         
               if (iorder(k) .eq. iatom1) xmol(jmol,i) = 
     .            4.825d15*u1*t(i)**1.5/ne(i)*dexp(-1.1605d4* 
     .            const(1,jmol)/t(i))*xatom(k)               
            enddo
         endif
      enddo


c*****compute matrix *c*, which is the derivative of each equation with     
c*****respect to each atom.                                                 
      do k=1,neq                                                      
         deltax(k) = -xfic(k) + xatom(k)                                
         do kk=1,neq                                                    
            c(k,kk) = 0.                                                
         enddo
         c(k,k) = 1.                                                    
         korder = iorder(k)                                             
         do kk=1,neq                                                    
            kderiv = iorder(kk)                                         
            do 28 j=1,nmax                                              
               jmol = ident(k,j)                                        
               if (jmol .eq. 0) go to 28                                
               call discov(amol(jmol),kderiv,num2)                      
               if (num2 .eq. 0) go to 28                                
               call discov(amol(jmol),korder,num1)                      
               c(k,kk) = c(k,kk) + xmol(jmol,i)*num1*num2/xatom(kk)     
               if (k .eq. kk) deltax(k) = deltax(k) + num1*xmol(jmol,i) 
28          continue                                                    
         enddo
      enddo


c*****calculate array 'xcorr', the change in 'xatom'. array 'xcorr' is      
c*****'deltax' multiplied by the inverse of 'c'                            
      call invert (neq,c,ans,30)                                        
      do k=1,neq                                                      
         x1 = xcorr(k)                                                  
         xcorr(k) = 0.                                                  
         do kk=1,neq                                                    
            xcorr(k) = xcorr(k) + ans(k,kk)*deltax(kk)                  
         enddo
      enddo


c*****decide if another iteration is needed                                 
      iflag = 0                                                         
      do k=1,neq                                                      
c*****here oscillations are damped out                                      
         if (x1*xcorr(k) .lt. -0.5*x1**2) xcorr(k) = 0.5*xcorr(k)       
         x1 = xatom(k)                                                  
         if (dabs(xcorr(k)/xatom(k)) .gt. 0.005) iflag = 1              
         xatom(k) = xatom(k) - xcorr(k)                                 
c*****fix element number densities which are ridiculous                     
         if (xatom(k).le.0.0 .or. xatom(k).ge.1.001*xfic(k)) then
            iflag = 1                                                   
            xatom(k) = 1.0d-2*dabs(xatom(k)+xcorr(k))                   
         endif
      enddo                                                           
      if (iflag .ne. 0) go to 27                                        


c*****print out atomic and molecular partial pressures. *xamol* is the      
c*****number density for each neutral atom                                  
      do k=1,neq                                                     
         xamol(k,i) = xatom(k)                                          
         patom(k) = dlog10(xatom(k)*tk)                                 
      enddo
      do jmol=1,nmol                                                 
         pmol(jmol) = dlog10(xmol(jmol,i)*tk)                           
      enddo
      if (molopt .ge. 2) then
         pglog = dlog10(pgas(lev))
         write (nf1out,1004) lev,int(t(lev)+0.001),pglog,
     .                       (patom(i),i=1,neq), (pmol(i),i=1,nmol)
      endif


c*****here the big loop in tau ends
21    continue                                                          


c*****finally, transfer some of this number density output to other arrays
      call setmols
      return                                                            


c*****format statements
1001  format ('I do not know this molecule: ',f10.0)              
1002  format (/'INPUT EQUILIBRIUM DATA:'/1x, 'species', 2x, 
     .        'D0/Chi ', 3x, 'const1', 6x, 'const2', 6x,
     .        'const3', 6x, 'const4', 6x, 'const5'/
     .        (0pf8.1, f8.3, 1p5d12.4))
1003  format (/'MOLECULAR EQUILIBRIUM SOLUTIONS:  (log partial',
     .        ' pressures listed under names)'/
     .        2x,'i',5x,'T',2x,'Pgas',8f8.1/(15x,8f8.1))
1004  format (i3,i6,f6.2,8f8.2/(15x,8f8.2))

      end                                                               



