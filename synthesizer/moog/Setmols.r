
      subroutine inmodel 
c******************************************************************************
c     This subroutine reads in the model 
c****************************************************************************** 

      implicit real*8 (a-h,o-z) 
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Quants.com'
      include 'Factor.com'
      include 'Dummy.com'
      include 'Pstuff.com'
      real*8 element(95), logepsilon(95)
      real*8 kaprefmass(100)
      real*8 bmol(110)
      integer nel(7)
      character list*80, list2*70
      integer append
      data (nel(i),i=1,7)/  1,  2,  6, 12, 13, 14, 26/


c*****Read in the key word to define the model type
      modelnum = modelnum + 1
      rewind nfmodel
      read (nfmodel,2001) modtype
      write (nf1out,1010) modtype
      if (modtype .eq. 'begn      ' .or.  modtype .eq. 'BEGN      ') 
     .   write (nf1out,1011)


c*****Read a comment line (usually describing the model)
      read (nfmodel,2002) moditle


c*****Read the number of depth points
      read (nfmodel,2002) list
      list2 = list(11:)
      read (list2,*) ntau
      if (ntau .gt. 100) then
         write (array,1012)
         call prinfo (10)
         stop
      endif


c*****EITHER: Read in a model from the output of the experimental new
c     MARCS code.  This modtype is called "NEWMARCS".  On each line 
c     the numbers are:
c     tau(5000), t, pe, pgas, rho,  model microtrubulent velocity,
c     and mean opacity (cm^2/gm) at the reference wavelength (5000A).
      if (modtype .eq. 'NEWMARCS  ') then
         read (nfmodel,*) wavref    
         do i=1,ntau
            read (nfmodel,*) tauref(i),t(i),ne(i),pgas(i),rho(i),
     .                         vturb(1),kaprefmass(i)
         enddo
c*****OR: Read in a model from the output of the on-line new
c     MARCS code.  This modtype is called "WEBMARCS".  On each line
c     the numbers are:
c     layer number (not needed), log{tau(Rosseland)} (not needed),
c     log{tau(5000)}, depth, t, pe, pgas, prad (not read in) and
c     pturb (not read in)
      elseif (modtype .eq. 'WEBMARCS') then
         read (nfmodel,*) wavref
         do i=1,ntau
            read (nfmodel,*) k, dummy1(k), tauref(i), dummy2(k), t(i),
     .                       ne(i), pgas(i)
         enddo
c*****OR: Read in a model from an alternative form of on-line new
c     MARCS code.  This modtype is called "WEB2MARC".  On each line
c     the numbers are:
c     atmospheric layer number (not needed), log{tau(5000)}, t, 
c     log(Pe), log(Pgas), rhox
      elseif (modtype .eq. 'WEB2MARC') then
         read (nfmodel,*) wavref
         do i=1,ntau
            read (nfmodel,*) k,tauref(i),t(i),ne(i),pgas(i),rhox(i)
         enddo
c     OR: Read in a model from the output of the ATLAS code.  This
c     modtype is called "KURUCZ".  On each line the numbers are:
c     rhox, t, pgas, ne, and Rosseland mean opacity (cm^2/gm), and
c     two numbers not used by MOOG.  
      elseif (modtype .eq. 'KURUCZ    ') then
         do i=1,ntau
            read (nfmodel,*) rhox(i),t(i),pgas(i),ne(i),kaprefmass(i)
         enddo
c     OR: Read in a model from the output of the NEXTGEN code.  This
c     modtype is called "NEXTGEN".  These models have tau scaled at a 
c     specific wavelength that is read in before the model. MOOG will 
c     need to generate the opacities internally.On each line the numbers 
c     are: tau, t, pgas, pe, density, mean molecular weight, two numbers
c     not used by MOOG, and Rosseland mean opacity (cm^2/gm).
      elseif (modtype .eq. 'NEXTGEN   ') then
         read (nfmodel,*) wavref
         do i=1,ntau
            read (nfmodel,*) tauref(i),t(i),pgas(i),ne(i), rho(i),  
     .                       molweight(i), x2, x3, kaprefmass(i)
         enddo
c     OR: Read in a model from the output of the MARCS code.  This modtype
c     type is called "BEGN".  On each line the numbers are:
c     tauross, t, log(pg), log(pe), mol weight, and kappaross.
      elseif (modtype .eq. 'BEGN      ') then
         do i=1,ntau
            read (nfmodel,*) tauref(i),t(i),pgas(i),ne(i),
     .                          molweight(i),  kaprefmass(i)
         enddo
c     OR: Read in a model generated from ATLAS, but without accompanying
c     opacities.  MOOG will need to generate the opacities internally,
c     using a reference wavelength that it reads in before the model.
      elseif (modtype .eq. 'KURTYPE') then
         read (nfmodel,*) wavref    
         do i=1,ntau
            read (nfmodel,*) rhox(i),t(i),pgas(i),ne(i)
         enddo
c     OR: Read in a model generated from ATLAS, with output generated
c     in Padova.  The columns are in somewhat different order than normal
      elseif (modtype .eq. 'KUR-PADOVA') then
         read (nfmodel,*) wavref
         do i=1,ntau
            read (nfmodel,*) tauref(i),t(i),kaprefmass(i),
     .                         ne(i),pgas(i),rho(i)
         enddo
c     OR: Read in a generic model that has a tau scale at a specific 
c     wavelength that is read in before the model.  
c     MOOG will need to generate the opacities internally.
      elseif (modtype .eq. 'GENERIC   ') then
         read (nfmodel,*) wavref    
         do i=1,ntau
            read (nfmodel,*) tauref(i),t(i),pgas(i),ne(i)
         enddo
c     OR: quit in utter confusion if those model types are not specified
      else
         write (*,1001)
         stop
      endif


c*****Compute other convenient forms of the temperatures
      do i=1,ntau
          theta(i) = 5040./t(i)
          tkev(i) = 8.6171d-5*t(i)
          tlog(i) = dlog(t(i))
      enddo


c*****Convert from logarithmic Pgas scales, if needed
      if (pgas(ntau)/pgas(1) .lt. 10.) then
         do i=1,ntau                                                    
            pgas(i) = 10.0**pgas(i)
         enddo
      endif


c*****Convert from logarithmic Ne scales, if needed
      if(ne(ntau)/ne(1) .lt. 20.) then
         do i=1,ntau                                                    
            ne(i) = 10.0**ne(i)
         enddo
      endif


c*****Convert from Pe to Ne, if needed
      if(ne(ntau) .lt. 1.0e7) then
         do i=1,ntau                                                    
            ne(i) = ne(i)/1.38054d-16/t(i)
         enddo
      endif


c*****compute the atomic partition functions
      do j=1,95
         elem(j) = dble(j)
         call partfn (elem(j),j)
      enddo


c*****Read the microturbulence (either a single value to apply to 
c     all layers, or a value for each of the ntau layers). 
c     Conversion to cm/sec from km/sec is done if needed
      read (nfmodel,2003) (vturb(i),i=1,6)
      if (vturb(2) .ne. 0.) then
         read (nfmodel,2003) (vturb(i),i=7,ntau) 
      else
         do i=2,ntau                                                    
            vturb(i) = vturb(1)
         enddo
      endif
      if (vturb(1) .lt. 100.) then
         write (moditle(55:62),1008) vturb(1)
         do i=1,ntau
            vturb(i) = 1.0e5*vturb(i)
         enddo
      else
         write (moditle(55:62),1008) vturb(1)/1.0e5
      endif


c*****Read in the abundance data, storing the original abundances in xabu
c*****The abundances not read in explicity are taken from the default
c*****solar ones contained in array xsolar.
      read (nfmodel,2002) list
      list2 = list(11:)
      read (list2,*) natoms,abscale
      write (moditle(63:73),1009) abscale
      if(natoms .ne. 0) 
     .         read (nfmodel,*) (element(i),logepsilon(i),i=1,natoms) 
      xhyd = 10.0**xsolar(1)
      xabund(1) = 1.0
      xabund(2) = 10.0**xsolar(2)/xhyd
      do i=3,95                                                      
         xabund(i) = 10.0**(xsolar(i)+abscale)/xhyd
         xabu(i) = xabund(i)
      enddo
      if (natoms .ne. 0) then
         do i=1,natoms                                                  
            xabund(idint(element(i))) = 10.0**logepsilon(i)/xhyd
            xabu(idint(element(i))) = 10.0**logepsilon(i)/xhyd
         enddo
      endif

c*****Compute the mean molecular weight, ignoring molecule formation
c     in this approximation (maybe make more general some day?)
      wtnum = 0.
      wtden = 0.
      do i=1,95
         wtnum = wtnum + xabund(i)*xam(i)
         wtden = wtden + xabund(i)
      enddo
      wtmol = wtnum/(xam(1)*wtden)
      nomolweight = 0
      if (modtype .eq. 'BEGN      ' .or. modtype .eq. 'NEXTGEN   ') then
         nomolweight = 1
      endif
      if (nomolweight .ne. 1) then
         do i=1,ntau
             molweight(i) = wtmol
         enddo
      endif

c*****Compute the density 
      if (modtype .ne. 'NEXTGEN   ') then
         do i=1,ntau                                                    
            rho(i) = pgas(i)*molweight(i)*1.6606d-24/(1.38054d-16*t(i))
         enddo
      endif


c*****Calculate the fictitious number density of hydrogen
c     Note:  ph = (-b1 + dsqrt(b1*b1 - 4.0*a1*c1))/(2.0*a1)
      iatom = 1
      call partfn (dble(iatom),iatom)
      do i=1,ntau    
         th = 5040.0/t(i)         
         ah2 = 10.0**(-(12.7422+(-5.1137+(0.1145-0.0091*th)*th)*th))
         a1 = (1.0+2.0*xabund(2))*ah2
         b1 = 1.0 + xabund(2)
         c1 = -pgas(i)         
         ph = (-b1/2.0/a1)+dsqrt((b1**2/(4.0*a1*a1))-(c1/a1))
         nhtot(i) = (ph+2.0*ph*ph*ah2)/(1.38054d-16*t(i))
      enddo


c*****Molecular equilibrium called here.
c     First, a default list of ions and molecules is considered. Then the 
c     user's list is read in. A check is performed to see if any of these 
c     species need to be added. If so, then they are appended to the 
c     default list. The molecular equilibrium routine is then called.
c     Certain species are important for continuous opacities and damping
c     calculations - these are read from the equilibrium solution and 
c     saved.


c*****Set up the default molecule list
      if (molset .eq. 0) then
         nmol = 30
      else
         nmol = 59
      endif
      if     (molset .eq. 0) then
         do i=1,110
            amol(i) = smallmollist(i)
         enddo
      elseif (molset .eq. 1) then
         do i=1,110
            amol(i) = largemollist(i)
         enddo
      else
         array = 'molset = 0 or 1 only; I quit!'
         call putasci (29,6)
         stop
      endif


c*****Read in the names of additional molecules to be used in 
c     molecular equilibrium if needed.
      read (nfmodel,2002,end=101) list
      list2 = list(11:)
      read (list2,*) moremol
      if (moremol .ne. 0) then
         read (nfmodel,*) (bmol(i),i=1,moremol)
         append = 1
         do k=1,moremol
            do l=1,nmol
               if (nint(bmol(k)) .eq. nint(amol(l))) 
     .         append = 0  
            enddo
            if (append .eq. 1) then 
               nmol = nmol + 1
               amol(nmol) = bmol(k)
            endif
            append = 1
         enddo  
      endif


c*****do the general molecular equilibrium
101   call eqlib
      call resetmol


c     In the number density array "numdens", the elements denoted by
c     the first subscripts are named in at the ends of the assignment
c     lines; at present these are the only ones needed for continuous 
c     opacities
c     the desired elements are H, He, C, Mg, Al, Si, Fe, along
c     with the H_2 molecule
c     first do the neutrals
      do j=1,7
         do k=1,neq
            if (nel(j) .eq. iorder(k)) then
               do i=1,ntau
                  numdens(j,1,i) = xamol(k,i)
               enddo
               exit
            endif
         enddo
      enddo
c     then do the ions
      do j=1,7
         ispec10 = nint(dble(10*nel(j)+1))
         do k=1,nmol
            kmol10 = nint(10.*amol(k))
            if (ispec10 .eq. kmol10) then
               do i=1,ntau
                  numdens(j,2,i) = xmol(k,i)
               enddo
               exit
            endif
         enddo
      enddo
c     finally add in H_2
      do k=1,nmol
         ispec = 101
         if (ispec .eq. nint(amol(k))) then
            do i=1,ntau
               numdens(8,1,i) = xmol(k,i)
            enddo
            exit
         endif
      enddo


c*****compute partitiion functions for H_2O and CO_2; 
      do i=1,ntau
         h2olog = 0.
         co2log = 0.
         if (t(i) .gt. 5000.) then
            uh2o(i) = 1.d8
            uco2(i) = 1.d8
         else
            do j=1,5
               h2olog = h2olog + h2ocoeff(j)*t(i)**(j-1)
               co2log = co2log + co2coeff(j)*t(i)**(j-1)
            enddo
            uh2o(i) = 10**h2olog
            uco2(i) = 10**co2log
         endif
      enddo


c*****transfer H_2O and CO_2 number densities from the molecular 
c     equilibrium output
c     note:  HITRAN partition functions are given only for T < 5000K
      do j=1,nmol
         if (nint(amol(j)) .eq. 10108) ih2o = j
         if (nint(amol(j)) .eq. 60808) ico2 = j
      enddo
      do i=1,ntau
         xnh2o(i) = xmol(ih2o,i)
         xnco2(i) = xmol(ico2,i)
      enddo
         
      
c*****SPECIAL NEEDS: for NEWMARCS models, to convert kaprefs to our units
      if (modtype .eq. 'NEWMARCS  ') then
         do i=1,ntau
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
c     SPECIAL NEEDS: for KURUCZ models, to create the optical depth array,
c     and to convert kaprefs to our units
      elseif (modtype .eq. 'KURUCZ    ') then
         first = rhox(1)*kaprefmass(1)
         tottau = rinteg(rhox,kaprefmass,tauref,ntau,first) 
         tauref(1) = first
         do i=2,ntau
            tauref(i) = tauref(i-1) + tauref(i)
         enddo
         do i=1,ntau
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
c     SPECIAL NEEDS: for NEXTGEN models, to convert kaprefs to our units
      elseif (modtype .eq. 'NEXTGEN   ') then
         do i=1,ntau                                                    
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
c     SPECIAL NEEDS: for BEGN models, to convert kaprefs to our units
      elseif (modtype .eq. 'BEGN      ') then
         do i=1,ntau                                                    
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
c     SPECIAL NEEDS: for KURTYPE models, to create internal kaprefs,
c     and to compute taurefs from the kaprefs converted to mass units
      elseif (modtype .eq. 'KURTYPE   ') then
         call opacit (1,wavref)
         do i=1,ntau                                                    
            kaprefmass(i) = kapref(i)/rho(i)
         enddo
         first = rhox(1)*kaprefmass(1)
         tottau = rinteg(rhox,kaprefmass,tauref,ntau,first) 
         tauref(1) = first
         do i=2,ntau
            tauref(i) = tauref(i-1) + tauref(i)
         enddo
c     SPECIAL NEEDS: for NEWMARCS models, to convert kaprefs to our units
      elseif (modtype .eq. 'KUR-PADOVA') then
         do i=1,ntau
            kapref(i) = kaprefmass(i)*rho(i)
         enddo
c     SPECIAL NEEDS: for generic models, to create internal kaprefs,
      elseif (modtype .eq. 'GENERIC   ' .or.
     .        modtype .eq. 'WEBMARCS  ' .or.
     .        modtype .eq. 'WEB2MARC  ') then
         call opacit (1,wavref)
      endif


c*****Convert from logarithmic optical depth scales, or vice versa.
c     xref will contain the log of the tauref
      if(tauref(1) .lt. 0.) then
         do i=1,ntau                                                    
            xref(i) = tauref(i)
            tauref(i) = 10.0**xref(i)
         enddo
      else
         do i=1,ntau                                                    
            xref(i) = dlog10(tauref(i))
         enddo
      endif


c*****Write information to output files
      if (modprintopt .lt. 1) return
      write (nf1out,1002) moditle
      do i=1,ntau
         dummy1(i) = dlog10(pgas(i))
         dummy2(i) = dlog10(ne(i)*1.38054d-16*t(i))
      enddo
      write (nf1out,1003) wavref,(i,xref(i),tauref(i),t(i),dummy1(i),
     .                    pgas(i),dummy2(i),ne(i),vturb(i),i=1,ntau)
      write (nf1out,1004)
      do i=1,95
         dummy1(i) = dlog10(xabund(i)) + 12.0
      enddo
      write (nf1out,1005) (names(i),i,dummy1(i),i=1,95)
      write (nf1out,1006) modprintopt, molopt, linprintopt, fluxintopt
      write (nf1out,1007) (kapref(i),i=1,ntau)
      return


c*****format statements
2001  format (a10)
2002  format (a80)
2003  format (6d13.0)
1001  format('permitted model types are:'/'KURUCZ, BEGN, ',
     .       'KURTYPE, KUR-PADOVA, NEWMARCS, WEBMARCS, NEXTGEN, ',
     .       'WEB2MARC, or GENERIC'/ 'MOOG quits!')
1002  format (/'MODEL ATMOSPHERE HEADER:'/a80/)
1003  format ('INPUT ATMOSPHERE QUANTITIES', 10x,
     .        '(reference wavelength =',f10.2,')'/3x, 'i', 2x, 'xref',
     .        3x, 'tauref', 7x, 'T', 6x, 'logPg', 4x, 'Pgas',
     .        6x, 'logPe', 5x, 'Ne', 9x, 'Vturb'/
     .        (i4, 0pf6.2, 1pd11.4, 0pf9.1, f8.3, 1pd11.4, 0pf8.3,
     .        1pd11.4, d11.2))
1004  format (/'INPUT ABUNDANCES: (log10 number densities, log H=12)'/
     .       '      Default solar abundances: Asplund et al. 2009')
1005  format (5(3x,a2, '(',i2,')=', f5.2))
1006  format (/'OPTIONS: atmosphere = ', i1, 5x, 'molecules  = ', i1/
     .        '         lines      = ', i1, 5x, 'flux/int   = ', i1)
1007  format (/'KAPREF ARRAY:'/(6(1pd12.4)))
1008  format ('vt=', f5.2)
1009  format (' M/H=', f5.2)
1010  format (13('-'),'MOOG OUTPUT FILE', 10('-'),
     .        '(MOOG version from 01/28/09)', 13('-')//
     .        'THE MODEL TYPE: ', a10)
1011  format ('   The Rosseland opacities and optical depths have ',
     .        'been read in')
1012  format ('HOUSTON, WE HAVE MORE THAN 100 DEPTH POINTS! I QUIT!')


      end

