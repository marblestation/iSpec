
      subroutine trudamp (linnum)
c*************************************************************************
c     This routine calculates damping parameters for lines that have
c     accurately known laboratory damping parameters
c     THE USER IS WARNED THAT THESE FORMULAE ARE OLD AND NOT EASILY 
c     DEFENDED; THE BARKLEM NUMBERS ARE TO BE PREFERRED.
c*************************************************************************
      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'


c*****begin with some calculations leading to a c6 value ("unsold")
      j = linnum
      iwave = int(wave1(j))
      iatom10 = nint(10.*atom1(j))
      ich = nint(charge(j))
      unsold = dabs(1.61d-33*(13.5*charge(j)/(chi(j,ich) -
     .         e(j,1)))**2 - 1.61d-33*(13.5*charge(j)/
     .         (chi(j,ich)-e(j,2)))**2)


c*****Ca II "IR triplet" lines at the Ca II K line at 3934 A
      if (iatom.eq.201 .and. iwave.eq.3933) then
         do i=1,ntau
            gammaa = 1.45d+8 
            gnature = gammaa + 0.5*1.5d-9*t(i)**(1/3)*numdens(1,1,i)
            gvander = 1.6d-8 * (t(i)/5000.)**0.3 * numdens(1,1,i)
            gstark = 3.0d-6 * ne(i)
            gammadamp = anature + gvander + gstark
            a(j,i) = gammadamp*wave1(j)*1.0e-8/(12.56636*dopp(j,i))
            write (nf1out,1001) gnature, gstark, gvander,
     .                          gammadamp, a(j,i)             
         enddo                                                
         return


c*****Ca II "IR triplet" lines at 8498, 8542, and 8662 A
      elseif (iatom10.eq.201 .and.
     .        (iwave.eq.8498.or.
     .         iwave.eq.8542.or.
     .         iwave.eq.8662)) then
         write (nf1out,1000)
         do i=1,ntau
            nhe = xabund(2)*nhtot(i)
            gnature = 1.5e+08
            gstark  = 1.5e-06*ne(i)*(t(i)/5000.)**0.1666
            ghelium = 3.0e-09*nhe*(t(i)/5000.)**0.4 
            ghydro  = 1.0e-08*nhtot(i)*(t(i)/5000.)**0.4 
            gammadamp = gnature/2. + gstark + ghelium + ghydro 
            a(j,i) = gammadamp*wave1(j)*1.0e-8/(12.56636*dopp(j,i))
            v1 = dsqrt(2.1175d8*t(i)*(1.0/amass(j)+1.008))  
            gvander = 17.0*unsold**0.4*v1**0.6*nhtot(i)       
            avander = gvander*wave1(j)*1.0d-8/(12.56636*dopp(j,i))  
            write (nf1out,1001) gnature, gstark, ghelium, ghydro,
     .                          gammadamp, a(j,i), avander
         enddo
         return


c*****Ca I 6717 A
      elseif (iatom10.eq.200 .and. iwave.eq.6717) then
         write (nf1out,1002) iwave
         do i=1,ntau
            nhe = xabund(2)*nhtot(i)
            gnature = 0.4e-08
            ghelium = 1.0e-09*nhe*(t(i)/5000.)**0.4 
            ghydro  = 2.0e-08*nhtot(i)*(t(i)/5000.)**0.4 
            gammadamp = gnature/2. + ghelium + ghydro 
            gammadamp = gammadamp*2.
            a(j,i) = gammadamp*wave1(j)*1.0e-8/(12.56636*dopp(j,i))
            v1 = dsqrt(2.1175d8*t(i)*(1.0/amass(j)+1.008)) 
            gvander = 17.0*unsold**0.4*v1**0.6*nhtot(i)       
            avander = gvander*wave1(j)*1.0d-8/(12.56636*dopp(j,i))  
            write (nf1out,1001) gnature, ghelium, ghydro,
     .                          gammadamp, a(j,i), avander
         enddo
         return


c*****Ca I 6318, 6343, 6361 A autoionization lines
      elseif (iatom10.eq.200 .and. 
     .        (iwave.eq.6318 .or.
     .         iwave.eq.6343 .or.
     .         iwave.eq.6361)) then
         write (nf1out,1005) iwave
         do i=1,ntau
            if (dampnum(j) .eq. 0) then
               gnature = 1.5d12
            else
               gnature = dampnum(j)*1.5d12
            endif
            gammadamp = gnature
            a(j,i) = gammadamp*wave1(j)*1.0e-8/(12.56636*dopp(j,i))
            write (nf1out,1001) gnature, gammadamp, a(j,i)
         enddo
         return


c*****Na I lines
      elseif (iatom .eq. 110) then
         write (nf1out,1003) iwave
         do i=1,ntau
            gnature = 2.21e+15/wave1(j)**2
            v1 = dsqrt(2.1175d8*t(i)*(1.0/amass(j)+1.008))
            gvander = 17.0*unsold**0.4*v1**0.6*nhtot(i)       
            gcoll   = gvander*2.1
            gammadamp = gnature/2. + gcoll
            a(j,i) = gammadamp*wave1(j)*1.0e-8/(12.56636*dopp(j,i))
            avander = gvander*wave1(j)*1.0d-8/(12.56636*dopp(j,i))  
            write (nf1out,1001) gnature, gcoll, gammadamp, 
     .                          a(j,i), avander
         enddo
         return


c*****CH autoionization line at 3693 A
      elseif (iatom10.eq.1060 .and.
     .        iwave.eq.3693) then
         write (nf1out,1006) iwave
         do i=1,ntau
            if (dampnum(j) .eq. 0) then
               gnature = 4.0d11
            else
               gnature = dampnum(j)*4.0d11
            endif
            gammadamp = gnature
            a(j,i) = gammadamp*wave1(j)*1.0e-8/(12.56636*dopp(j,i))
            write (nf1out,1001) gnature, gammadamp, a(j,i)
         enddo
         return
      endif


c*****format statements
1000  format(//' LINE BROADENING PARAMETERS FOR ',
     .       'CaII INFRARED TRIPLET'/
     .       4x,'natural',6x,'stark',5x,'helium',3x,'hydrogen',
     .       6x,'gammadamp',5x,'a(j,i)',5x,'a(VdW)')
1001  format (1p7e11.3)
1002  format(//' LINE BROADENING PARAMETERS FOR ',
     .       'CaI LINE AT',i6/
     .       4x,'natural',5x,'helium',3x,'hydrogen',
     .       6x,'gammadamp',5x,'a(j,i)',5x,'a(VdW)')
1003  format(//' LINE BROADENING PARAMETERS FOR ',
     .       'NaI LINE AT',i6/
     .       4x,'natural',4x,'VdW*1.2',6x,'gammadamp',
     .       5x,'a(j,i)',5x,'a(VdW)')
1004  format(//' LINE BROADENING PARAMETERS FOR ',
     .       'CaI LINE AT 6318A'/
     .       4x,'natural',6x,'Stark',4x,'vdwaals',
     .       6x,'gammadamp',5x,'a(j,i)')
1005  format(//' LINE BROADENING PARAMETERS FOR ',
     .       'THE CaI AUTOIONIZATION LINE AT',i6/
     .       4x,'natural', 2x,'gammadamp',5x,'a(j,i)')
1006  format(//' LINE BROADENING PARAMETERS FOR ',
     .       'THE CH AUTOIONIZATION LINE AT',i6/
     .       4x,'natural', 2x,'gammadamp',5x,'a(j,i)')


      end




