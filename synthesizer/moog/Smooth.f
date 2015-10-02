
      subroutine smooth (line,ncall)
c******************************************************************************
c     This subroutine prepares synthesis data for plotting by smoothing
c     the data in various ways.
c******************************************************************************
 
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      real*4 shortnum
      character*5 abchars
      character*4 isochars


c*****initialize parameters
      write (abitle(1:400),1081)
      write (isoitle(1:240),1082)
      nsyn = 1
      if (ncall .eq. 1) then
         gaussflag = 'f'
         rotateflag = 'f'
         lorenflag = 'f'
         macroflag = 'f'
         if (choice .ne. 'l') addflux = 0.
      endif


c*****on entering, figure out what kind of smoothing is desired, unless
c     the default smoothing options have been set for first pass 
      if (line .gt. 0) then
2        write (array,1007)
         istat = ivwrite (line,3,array,67)
         write (array,1004)
         istat = ivwrite (line+1,3,array,50)
         array = 'What is your choice? '
         nchars = 21
         call getasci (nchars,line+2)
         smtype = chinfo(1:1)
         if (smtype.ne.'n' .and. smtype.ne.'g' .and.
     .       smtype.ne.'l' .and. smtype.ne.'v' .and.
     .       smtype.ne.'m' .and. smtype.ne.'c' .and.
     .       smtype.ne.'d' .and. smtype.ne.'r' .and.
     .       smtype.ne.'p') go to 2
      endif


c     if a user-specified variable Gaussian smoothing over the spectral range
c     is called for, option 'p', branch to a different routine
      if (smtype .eq. 'p') then
         call vargauss (line+1)
         return
      endif
 

c*****rewind the unsmoothed and smoothed spectrum output 
c     files, and get the synthesis range parameters from the 'dump' file
      rewind nf2out
      rewind nf3out
      do i=1,20
         read (nf2out,1002) array
         if     (array(1:7).eq.'element' .or.  
     .       array(1:7).eq.'Changin' .or.
     .       array(1:7).eq.'ALL abu' .or.
     .       array(1:7).eq.'Isotopi') then
            cycle
         elseif (array(1:7).eq.'MODEL: ') then
             moditle(1:73) = array(8:80)
            read (nf2out,*) start, sstop, step
            kount = int((sstop - start + (step/4.0) )/step) + 1
            exit
         endif
      enddo


c*****branch to the desired smoothing function
      write (smitle,1010) smtype
      ism = 11
      if     (smtype .eq. 'l') then
         lorenflag = 't'
      elseif (smtype .eq. 'g') then
         gaussflag = 't'
      elseif (smtype .eq. 'v') then
         rotateflag = 't'
      elseif (smtype .eq. 'c') then
         rotateflag = 't'
         gaussflag = 't'
      elseif (smtype .eq. 'm') then
         macroflag = 't'
      elseif (smtype .eq. 'd') then
         macroflag = 't'
         gaussflag = 't'
      elseif (smtype .eq. 'r') then
         macroflag = 't'
         rotateflag = 't'
         gaussflag = 't'
      endif


c*****compute a stellar rotational broadening function; this follows 
c     D. F. Gray, 1976, "The Obs. & Anal. of Stell. Phot", p394-9
      if (rotateflag .eq. 't') then
32       array = 'GIVE THE STELLAR vsini [0.0]: '
         nchars = 30
         if (line .gt. 0) then
            call getnum (nchars,line+2,vsini,shortnum)
            if (vsini .eq. -9999.) vsini = 0.
         endif
         if (vsini .lt. 0.0) go to 32
         write (smitle(ism+1:ism+13),1011) vsini
         ism = ism + 13
31       array = 'GIVE THE LIMB DARKENING COEFFICIENT [0.0]: '
         nchars = 43
         if (line .gt. 0) then
            call getnum (nchars,line+2,limbdark,shortnum)
            if (limbdark .eq. -9999.) limbdark = 0.
         endif
         if (limbdark .lt. 0.0) go to 31
         write (smitle(ism+1:ism+13),1012) limbdark
         ism = ism + 13
         dlamlim = (start+sstop)/2.*vsini/3.0e5
         if (step .ge. dlamlim) then
            rotateflag = 'f'
         else
            pi = 3.141527
            bottom = dlamlim*pi*(1.-limbdark/3.)
            c1 = 2.*(1.-limbdark)/bottom
            c2 = 0.5*limbdark*pi/bottom
            prot0 = c1 + c2
            powerrot = prot0
            jdelrot = idint(dlamlim/step)
            if (jdelrot .gt. 1000) then
               write (*,1026) 
               smtype = 'e'
               return
            elseif (jdelrot .gt. kount/4) then
               write (*,1028)
               smtype = 'e'
               return
            endif
            do i=1,jdelrot
               dlam = (i-1)*step
               term = 1. - (dlam/dlamlim)**2
               prot(i) = c1*dsqrt(term) + c2*term
               powerrot = powerrot + 2.0*prot(i)
            enddo
         endif
      endif

         
c*****compute a macroturbulent smoothing function (uses subroutine vmacro)
      if (macroflag .eq. 't') then
51       array = 'GIVE THE MACROTURBULENT VELOCITY [0.0]: '
         nchars = 39
         if (line .gt. 0) then
            call getnum (nchars,line+2,vmac,shortnum)
            if (vmac .eq. -9999.) vmac = 0.
         endif
         if (vmac .lt. 0.0) go to 51
         write (smitle(ism+1:ism+13),1013) vmac
         ism = ism + 13
         if (vmac .eq. 0.) then
            macroflag = 'f'
         else
            wavemac = (start+sstop)/2.*vmac/3.0e5
            powermac = 1.0
            do i=1,1000
               wavei = step*i/wavemac
               pmac(i) = vmacro(wavei)
               powermac = powermac + 2.0 *pmac(i)
               if (pmac(i) .lt. 0.002) then
                  jdelmac = i
                  exit
               endif 
            enddo
            if (jdelmac .gt. 1000) then
               write (*,1025) wavemac
               smtype = 'e'
               return
            elseif (jdelmac .gt. kount/4) then
               write (*,1022)
               smtype = 'e'
               return
            endif   
         endif
      endif


c*****compute a Gaussian smoothing function
      if (gaussflag .eq. 't') then
11       array = 'GIVE THE FWHM OF THE GAUSSIAN FUNCTION [0.0]: ' 
         nchars = 46
         if (line .gt. 0) then
            call getnum (nchars,line+2,fwhmgauss,shortnum)
            if (fwhmgauss .eq. -9999.) fwhmgauss = 0.
         endif
         if (fwhmgauss .lt. 0.0) go to 11
         write (smitle(ism+1:ism+18),1014) fwhmgauss
         ism = ism + 18
         if (fwhmgauss .eq. 0.) then
            gaussflag = 'f'
         else
            sigma = fwhmgauss/2.
            aa = 0.6932/sigma**2
            power = 1.0
            do i=1,1000
               p(i) = dexp(-aa*(step*i)**2 )
               power = power + 2*p(i)
               if (p(i) .lt. 0.02) then
                  jdel = i
                  exit
               endif
            enddo
            if (jdel .gt. 1000) then
               write (*,1029) sigma
               smtype = 'e'
               return   
            elseif (jdel .gt. kount/4) then
               write (*,1021)
               smtype = 'e'
               return
            endif
         endif
      endif


c*****compute a Lorenzian smoothing function
      if (lorenflag .eq. 't') then
21       array = 'GIVE THE FWHM OF THE LORENTZIAN FUNCTION [0.0]: '
         nchars = 48
         if (line .gt. 0) then
            call getnum (nchars,line+2,fwhmloren,shortnum)
            if (fwhmloren .eq. -9999.) fwhmloren = 0.
         endif
         if (fwhmloren .lt. 0.0) go to 21
         write (smitle(ism+1:ism+20),1015) fwhmloren
         ism = ism + 20
         if (fwhmloren .eq. 0.) then
            lorenflag = 'f'
         else
            sigma = fwhmloren/2.
            power = 1.0
            do i=1,1000
               p(i) = ((sigma**2)/((sigma**2)+((step*i)**2)))
               power = power + 2.0 *p(i)
               if (p(i) .lt. 0.02) then
                  jdel = i
                  exit
               endif
            enddo    
            if (jdel .gt. 1000) then
               write (*,1030) sigma
               smtype = 'e'
               return  
            elseif (jdel .gt. kount/4) then
               write (*,1031)
               smtype = 'e'
               return
            endif
         endif
      endif


c*****finally there is a big loop that will read in the raw synthetic spectra,
c*****beginning with grabbing the information that appears before the array
c*****spectrum depth array,  Note that after the depth array is input, the
c*****code will flip to a flux scale
      rewind nf2out
      do nsyn=1,100


c*****here is the reading/grabbing of stuff preceding the depth array:
         naboff = 80*(nsyn-1)
         nabunds = 0
         nisos = 0
         do i=1,20
            read (nf2out,1002,end=2000) array
            if     (array(1:7).eq.'ALL abu') then
               cycle
            elseif (array(1:7).eq.'Changin') then
               abitle (naboff+1:naboff+23) = '[M/H] FOR ALL ELEMENTS:'
               abitle (naboff+24:naboff+29) = array(32:37)
               cycle
            elseif (array(1:7).eq.'element') then
               nabunds = nabunds + 1
               if (control .eq. 'binary ') then
                  if (nabunds .le. 5) then
                     ioff = naboff + 16*(nabunds-1)
                     abitle(ioff+1:ioff+2) = array(9:10)
                     abitle(ioff+3:ioff+14) = array(26:37)
                     abitle(ioff+15:ioff+16) = '  '
                  endif
               else
                  if (nabunds .le. 8) then
                     ioff = naboff + 9*(nabunds-1)
                     abitle(ioff+1:ioff+2) = array(9:10)
                     read (array(26:32),*) abnum
                     if     (abnum .gt. 0) then
                        write (abchars,1061) abnum
                     elseif (abnum .le. -10.) then
                        write (abchars,1062) abnum
                     else
                        write (abchars,1063) abnum
                     endif
                     abitle(ioff+3:ioff+7) = abchars
                     abitle(ioff+8:ioff+9) = '  '
                  endif
               endif
               cycle
            elseif (array(1:7).eq.'Isotopi') then
               nisos = nisos + 1
               if (nisos .le. 6) then
                  read (array(37:46),1050) ratio
                  if     (ratio .ge. 1000.) then 
                     write (isochars,1054) int(ratio)
                  elseif (ratio .ge.  100.) then
                     write (isochars,1051) ratio
                  elseif (ratio .ge.   10.) then 
                     write (isochars,1052) ratio
                  else
                     write (isochars,1053) ratio
                  endif
                  if (nsyn .eq. 1) then
                     ioff = 40*(nisos-1) + 5*(nsyn-1)
                     isoitle(ioff+1:ioff+10) = array(23:32)
                     isoitle(ioff+11:ioff+12) = ': '
                     isoitle(ioff+13:ioff+16) = isochars(1:4)
                     isoitle(ioff+17:ioff+17) = '/'
                  else
                     ioff = 12 + 40*(nisos-1) + 5*(nsyn-1)
                     isoitle(ioff+1:ioff+4) = isochars(1:4)
                     isoitle(ioff+5:ioff+5) = '/'
                  endif
               endif
            elseif (array(1:7).eq.'MODEL: ') then
               read (nf2out,1002) array
               exit
            endif
         enddo


c*****here is the actual reading of the depth array
         read (nf2out,1003,end=2000) (y(i),i=1,kount)
         do i=1,kount
            y(i) = 1.0 - y(i)
         enddo


c*****here a veiling addition can be added in
         if (addflux .gt. 0.0) then
             do i=1,kount
                 y(i) = (y(i) + addflux)/(1.0+addflux)
             enddo
         endif


c*****apply the rotational broadening if desired
         if (rotateflag .eq. 't') then
            min = jdelrot + 1
            max = kount - jdelrot
            do i=1,jdelrot
               z(i) = 1.
               z(kount-i+1) = 1.
            enddo
            do i=min,max
               z(i) = prot0*y(i)
               do j=1,jdelrot
                  z(i) = z(i) + prot(j)*(y(i-j) + y(i+j))
               enddo
               z(i) = z(i)/powerrot
            enddo
            do i=1,kount
               y(i) = z(i)
            enddo
         endif


c*****apply the macroturbulent broadening if desired
         if (macroflag .eq. 't') then
            min = jdelmac + 1
            max = kount - jdelmac
            do i=1,jdelmac
               z(i) = 1.
               z(kount-i+1) = 1.
            enddo
            do i=min,max
               z(i) = y(i)
               do j=1,jdelmac
                  z(i) = z(i) + pmac(j)*(y(i-j) + y(i+j))
               enddo
               z(i) = z(i)/powermac
            enddo
            do i=1,kount
               y(i) = z(i)
            enddo
         endif


c*****apply the Gaussian or Lorenzian smoothing if desired (this
c     is an either/or situation; only one of these can apply.
         if (gaussflag .eq. 't' .or. lorenflag .eq. 't') then
            min = jdel + 1
            max = kount - jdel
            do i=1,jdel
               z(i) = 1.
               z(kount-i+1) = 1.
            enddo
            do i=min,max
               z(i) = y(i)
               do j=1,jdel
                  z(i) = z(i) + p(j)*(y(i-j) + y(i+j))
               enddo
               z(i) = z(i)/power
            enddo
            do i=1,kount
               y(i) = z(i)
            enddo
         endif


c*****move the smoothed spectrum (or unsmoothed, if nothing has
c     been done to the y-array) into the appropriate array
         do i = 1,kount
            chunk(i,nsyn) = y(i)
         enddo
   

c*****compute the wavelength array; must be done for each synthetic
c     spectrum because of the way the equivalences were set up
         if     (iunits .eq. 1) then
            do i=1,kount
               xsyn(i) = 1.d-4*(start + (i-1)*step)
            enddo
         else
            do i=1,kount
               xsyn(i) = start + (i-1)*step
            enddo
         endif


c*****dump the smoothed spectrum in a MONGO-style set of 
c     (wavelength,flux) point pairs
         write (nf3out,1005) kount,start,sstop,step
         if (xsyn(1) .le. 100.0) then
            write (nf3out,1009) (xsyn(i),chunk(i,nsyn),i=1,kount)
         else 
            write (nf3out,1008) (xsyn(i),chunk(i,nsyn),i=1,kount)
         endif
      enddo


c*****exit the routine normally
2000  nsyn = nsyn - 1
      return


c*****a problem has developed in user-specified parameter (like a Gaussian
c     FWHM being too big); print out a warning and let user decide what
c     to do next


c*****format statements
1002  format (a80)
1003  format (10f7.4)
1004  format ('           c=v+g, d=m+g, r=m+v+g, p=VARIABLE GAUSS')
1005  format ('the number of points per synthesis = ',i5/
     .        'start = ',f10.3,5x,'stop = ',f10.3,5x,'step = ',f10.3)
1007  format ('SMOOTHING: n=NONE, g=GAUSS, l=LORENZ, ',
     .        'v=ROTATION, m=MACROTURBULENCE')
1008  format (f10.3,'  ',f10.5)
1009  format (f10.6,'  ',f10.5)
1010  format ('smoothing=',a1,69x)
1011  format ('  Vsini=',f5.1)
1012  format ('  L.D.C.=',f4.2)
1013  format ('  Vmacro=',f4.1)
1014  format ('  FWHMgauss=',f6.3)
1015  format ('  FWHMlorentz=',f6.3)
1021  format ('ERROR: GAUSSIAN PROFILE COVERS WHOLE ',
     .        'SPECTRUM!'/
     .        '       INCREASE SYNTHESIS LENGTH OR DECREASE ',
     .        'STEP SIZE?')
1022  format ('ERROR: MACROTURBULENT PROFILE COVERS WHOLE ',
     .        'SPECTRUM!'/
     .        '       INCREASE SYNTHESIS LENGTH OR DECREASE ',
     .        'STEP SIZE?')
1025  format ('ERROR: MACROTURBULENT PROFILE TOO BIG! ',
     .        '(WAVELENGTH WIDTH=',f6.2,')'/
     .        '       INCREASE SPECTRUM STEP SIZE?')
1026  format ('ERROR: ROTATIONAL PROFILE TOO BIG!'/  
     .        '       INCREASE SPECTRUM STEP SIZE?')
1028  format ('ERROR: ROTATIONAL PROFILE COVERS WHOLE ',
     .        'SPECTRUM!'/
     .        '       INCREASE SYNTHESIS LENGTH OR DECREASE ',
     .        'STEP SIZE?')
1029  format ('ERROR: GAUSSIAN PROFILE TOO BIG! ',
     .        '(HALF WIDTH=',f6.2,')'/
     .        '       INCREASE SPECTRUM STEP SIZE?')
1030  format ('ERROR: LORENZIAN PROFILE TOO BIG! ',
     .        '(HALF WIDTH=',f6.2,')'/
     .        '       INCREASE SPECTRUM STEP SIZE?')
1031  format ('ERROR: LORENZIAN PROFILE COVERS WHOLE ',
     .        'SPECTRUM!'/
     .        '       INCREASE SYNTHESIS LENGTH OR DECREASE ',
     .        'STEP SIZE?')
1050  format (f10.3)
1051  format (f4.0)
1052  format (f4.1)
1053  format (f4.2)
1054  format (i4)
1061  format ('+', f4.2)
1062  format (f5.1)
1063  format (f5.2)
1081  format (400(' '))
1082  format (240(' '))

      end





