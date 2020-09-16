
      subroutine nearly (numpass)
c******************************************************************************
c     This routine calculates the quantity 'kapnu0' for each line and       
c     each tau; this is the line opacity at line center.  When kapnu0 is 
c     multiplied by the voigt function, one gets the line opacity at an 
c     arbitrary wavelength.  Then, the Doppler factor 'dopp' and damping
c     parameter 'a' are also computed ('a' is done in a call to another
c     routine
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Quants.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Dummy.com'
      real*8 xnum(100)
      integer numpass
      equivalence (xnum,dummy1)


c*****load in data for damping factors to be computed from the data
c     of Barklem, if desired
      if (numpass.eq.1 .and. dampingopt.eq.1) call gammabark


c*****Locate the atmosphere level where taulam is near 0.5
      if (numpass.eq.1 .or. numpass.eq.3) then
         call opacit (2,wave1(1))
         do i=1,ntau
            if (taulam(i) .ge. 0.5) go to 180
         enddo
180      jtau5 = i
      endif


c*****if reading in a new line list, do all lines (numpass=1); 
c     if iterating "abfind" on a species affected by mol. eq.,
c     do just that species (numpass=2); if doing a fake line, or
c     maybe one line for a special purpose, do just the first line 
c     in the "list" (numpass=3)
      if     (numpass .eq. 1) then
         j1 = 1
         j2 = nlines + nstrong
      elseif (numpass .eq. 2) then
         j1 = lim1line
         j2 = lim2line
      elseif (numpass .eq. 3) then
         j1 = 1
         j2 = 1
      endif


c*****now make the calculations:  set up some parameters
      do j=j1,j2
         ich = idint(charge(j) + 0.1)                    
         iatom = atom1(j) + 0.0001
         factoriso = 1.0
    
                  
c*****compute the Doppler factors
         if (numpass.eq.1 .or. numpass.eq.3) then
            do i=1,ntau
               dopp(j,i) = dsqrt(1.6631d8*t(i)/amass(j)+vturb(i)**2)
            enddo


c*****compute damping parameters in a separate routine
            call damping (j)
         endif


c*****either: compute lower state number densities for atomic lines;
c     q21 is the ion/neutral ratio, etc., and q is the ratio of the total
c     to the species of interest; do the Saha equation first, then
c     the Boltzmann equation
         if     (iatom .lt. 100) then
            do i=1,ntau
               q21 = 4.825d15*u(iatom,2,i)/(u(iatom,1,i)*ne(i))*
     .               t(i)**1.5*dexp(-chi(j,1)/tkev(i)) 
               q32 = 4.825d15*u(iatom,3,i)/(u(iatom,2,i)*ne(i))*
     .               t(i)**1.5*dexp(-chi(j,2)/tkev(i)) 
               q43 = 4.825d15*u(iatom,4,i)/(u(iatom,3,i)*ne(i))*
     .               t(i)**1.5*dexp(-chi(j,3)/tkev(i)) 
               if (neq .ne. 0) then
                  do n=1,neq
                     if (iatom .eq. iorder(n)) then
                        xxnum = xamol(n,i)
                        if     (ich .eq. 1) then
                           q = 1.0
                        elseif (ich .eq. 2) then
                           q = 1.0/q21
                        elseif (ich .eq. 3) then
                           q = 1.0/(q21*q32) + 1.0/q32 + 1.0 + q43
                        endif
                        xnum(i) = xxnum/q*dexp(-e(j,1)/tkev(i))/
     .                            u(iatom,ich,i)   
                     endif
                  enddo
               endif
               if     (ich .eq. 1) then
                  q = 1.0 + q21 + q32*q21 + q43*q32*q21
               elseif (ich .eq. 2) then
                  q = 1.0/q21 + 1.0 + q32 + q43*q32
               elseif (ich .eq. 3) then
                  q = 1.0/(q21*q32) + 1.0/q32 + 1.0 + q43
               endif
               if (control .eq. 'abandy') then
                  xxab = xabund(iatom)*10**deltaabund
               else
                  xxab = xabund(iatom)
               endif
               xnum(i) = xxab*nhtot(i)/q*
     .                 dexp(-e(j,1)/tkev(i))/u(iatom,ich,i)
            enddo


c*****or: compute lower state number densities for molecular lines
         elseif (iatom .lt. 10000) then
            call sunder(atom1(j),iaa,ibb)
            do n=1,neq
               if(iorder(n) .eq. iaa) ia = n
               if(iorder(n) .eq. ibb) ib = n
            enddo
            do i=1,ntau
               psipri = 
     .             1.38065d-16*t(i)*10.0**(d0(j)*theta(i)-13.670)*
     .             theta(i)**2.5/
     .             (rdmass(j)**1.5*u(iaa,1,i)*u(ibb,1,i)) 
               xnum(i) = dexp(-e(j,1)/tkev(i))*psipri*
     .                   xamol(ia,i)*xamol(ib,i)
            enddo
         else
            if     (iatom .eq. 10108) then
               do i=1,ntau
                  xnum(i) = xnh2o(i)/uh2o(i)*dexp(-e(j,1)/tkev(i))
               enddo
            elseif (iatom .eq. 60808) then
               do i=1,ntau
                  xnum(i) = xnco2(i)/uco2(i)*dexp(-e(j,1)/tkev(i))
               enddo
            else
               write (*,1001) iatom
               stop
            endif
         endif


c*****finally, compute line opacities at line centers
         if (atom1(j)-float(iatom) .ge. 0.0) then
            do n=1,numiso
               if (atom1(j) .eq. isotope(n)) then
                  factoriso = isoabund(n,isorun)
               endif
            enddo
         endif
         do i=1,ntau
            kapnu0(j,i) = 2.65386d-2*xnum(i)*gf(j)*wave1(j)*1.0d-8/
     .                    dopp(j,i)*(1.0-dexp(-1.43879d8/
     .                    (wave1(j)*t(i))))/factoriso
         enddo
         strength(j) = kapnu0(j,jtau5)
      enddo


c*****output regular line information, and strong line information 
c     if appropriate; exit routine normally
      if (numpass.eq.1 .or. numpass.eq.3) then
         if (linprintopt .ge. 0) then
            if (dostrong .gt. 0) call lineinfo (2)
            call lineinfo (1)
         endif
      endif
      return


c*****format statements
1001  format (i6, 'IS NOT AN ALLOWED TRIATOMIC (H_2O AND CO_2);',
     .        ' I QUIT!')


      end


