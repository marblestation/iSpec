
      subroutine lineinfo (number)
c******************************************************************************
c     This routine is an output subroutine for things having to do with
c     lines.  
c            number = 1 gives input line data information
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Mol.com'
      include 'Dummy.com'
      include 'Factor.com'
      include 'Pstuff.com'
      include 'Quants.com'
      real*8 loggf, logstrength
      integer ifresh
      character*8 molname
      character*4 ion(3)
      character*1 name*1
      data ifresh /0/
      data ion/' I  ', ' II ', ' III'/


      go to (1,2,3), number


c*****here the line data are output to "standard_out"; all relevant 
c     drivers use this
c     if you don't want any line output, linprintopt=0 will exit the routine

1     if (linprintopt .lt. 1) return

c     if you want standard output, linprintopt=1 is chosen
c     linprintopt>=2 outputs ionization potentials, charges, masses,
c                            reduced masses for molecules, 
c     linprintopt>=3 outputs partition functions
c     lineprintop =4 outputs line-center opacities
      write (nf1out,1001) nlines
      if (linprintopt .ge. 2) write (nf1out,1002)
      do j=1,nlines
         ich = idint(charge(j) + 0.1)
         iatom = idint(atom1(j))
         loggf = dlog10(gf(j))
         logstrength = dlog10(strength(j))
         if     (iatom .lt. 100) then
            if (iunits .eq. 1) then
               write (nf1out,1003) j, 1.d-4*wave1(j), names(iatom),
     .           ion(ich), atom1(j), e(j,1), loggf, damptype(j), 
     .           logstrength, 1000.*width(j)
            else
               write (nf1out,1004) j, wave1(j), names(iatom),
     .           ion(ich), atom1(j), e(j,1), loggf, damptype(j),
     .           logstrength, 1000.*width(j)
            endif
            if (linprintopt .ge. 2) write (nf1out,1005) 
     .                 (chi(j,k),k=1,3), charge(j), amass(j), rdmass(j)
         elseif (iatom .lt. 10000) then
            call sunder (atom1(j),i1,i2)
            if (i1 .eq. 1) then
               l = i1
               i1 = i2
               i2 = l
            endif
            leftovr = idint(10000.*(atom1(j)-iatom)+0.1)
            if (i1 .lt. 10) then
               read (names(i1),1006) name
               write (molname,1007) name,names(i2),leftovr
            else
               write (molname,1008) names(i1),names(i2),leftovr
            endif
            if (iunits .eq. 1) then
               write (nf1out,1009) j, 1.d-4*wave1(j), molname, 
     .               atom1(j), e(j,1), loggf, damptype(j), 
     .               logstrength, 1000.*width(j)
            else
               write (nf1out,1010) j, wave1(j), molname, 
     .               atom1(j), e(j,1), loggf, damptype(j), 
     .               logstrength, 1000.*width(j)
            endif
            if (linprintopt .ge. 2) 
     .         write (nf1out,1005) 
     .               d0(j), (chi(j,k),k=1,2), charge(j), amass(j), 
     .               rdmass(j)
         elseif (iatom .lt. 1000000) then
            call sunder (atom1(j),i1,i2)
            xia = dble(i2)
            call sunder (xia,i2,i3)
            if (iatom .eq. 10108) then
               molname = 'H_2O    '
            else
               molname = 'CO_2    '
            endif
            if (iunits .eq. 1) then
               write (nf1out,1009) j, 1.d-4*wave1(j), molname,
     .               atom1(j), e(j,1), loggf, damptype(j),
     .               logstrength, 1000.*width(j)
            else
               write (nf1out,1010) j, wave1(j), molname,
     .               atom1(j), e(j,1), loggf, damptype(j),
     .               logstrength, 1000.*width(j)
            endif
            if (linprintopt .ge. 2)
     .         write (nf1out,1005)
     .               d0(j), (chi(j,k),k=1,2), charge(j), amass(j),
     .               rdmass(j)
         endif
      enddo    
      if (start.ne.0.0 .or. sstop.ne.0.0) then
         if (iunits .eq. 1) then
            write (nf1out,1011) oldstart,oldstop,oldstep,olddelta
         else 
            write (nf1out,1012) start,sstop,step,delta  
         endif
      if (rwlow .ne. 0.) write (nf1out,1013) rwlow, rwhigh, rwstep
      endif  
      if (linprintopt .ge. 3) then
         write (nf1out,1014)
         do j=1,95
            if (elem(j) .ne. 0.) then
               iatom = int(elem(j))
               write (nf1out,1015) iatom, names(iatom), xam(j),
     .                            xchi1(j), xchi2(j), xchi3(j)
               do k=1,4
                  write (nf1out,1016) k-1,(u(j,k,i),i=1,ntau)   
               enddo
            endif
         enddo
      endif
      if (linprintopt .ge. 4) then
         write (nf1out,1001)
         do j=1,nlines
            write (nf1out,1002) j,(kapnu0(j,i),i=1,ntau)
         enddo
      endif
      return


c*****here the STRONG line data are output; MOOG assumes that no
c     molecular line can possibly be in this category
2     write (nf1out,2001) nstrong
      do j=nlines+1,nlines+nstrong
         ich = idint(charge(j) + 0.1)
         iatom = idint(atom1(j))
         loggf = dlog10(gf(j))
         logstrength = dlog10(strength(j))
         if (iatom .lt. 100) then
            if (iunits .eq. 1) then
               write (nf1out,1003) j-nlines,1.d-4*wave1(j),names(iatom),
     .                             ion(ich), atom1(j), e(j,1), loggf,
     .                             damptype(j), logstrength
            else
               write (nf1out,1004) j-nlines, wave1(j),names(iatom),
     .                             ion(ich), atom1(j), e(j,1), loggf,
     .                             damptype(j), logstrength
            endif
         else
            write (*,2004) iatom
            stop
         endif
      enddo
      printstrong = 1
      return


c*****results of force-fitting EW to yield abundances are output here
c     look here also for the calls to the trend line calculations
3     if (ifresh .eq.0) then
         write (nf2out,3001) linitle,moditle
         ifresh = 1
      endif
      if (cogatom .eq. 0.) then
         iatom = iabatom
      else
         iatom = idint(cogatom)
      endif
      xab = dlog10(xabund(iatom)) + 12.
      ich = idint(charge(lim1obs) + 0.1)
      if (atom1(lim1obs) .lt. 100.) then
         write (array,3002) names(iatom), ion(ich) ,xab
         line = 1
         call prinfo (line)
         write (nf2out,*)
         write (nf2out,3002) names(iatom), ion(ich), xab
         write (array,3003)
         line = 2
         call prinfo (line)
         write (nf2out,3003)
      else
         call sunder (atom1(lim1obs),ia,ib)
         if (ia .eq. 1) then
            l = ia
            ia = ib
            ib = l
         endif
         leftovr = idint(10000.*(atom1(lim1obs)-iatom)+0.1)
         if (ia .lt. 10) then
            read (names(ia),1006) name
            write (molname,1007) name,names(ib)
         else
            write (molname,1008) names(ia),names(ib)
         endif
         write (array,3004) molname,xab
         line = 1
         call prinfo (line)
         write (nf2out,*)
         write (nf2out,3004) molname,xab
         write (array,3005) names(iabatom)
         line = 2
         call prinfo (line)
         write (nf2out,3005) names(iabatom)
         write (array,3006)
         line = 3
         call prinfo (line)
         write (nf2out,3006)
      endif
      do l=lim1obs,lim2obs
         if (abundout(l) .ne. 999.99) then
            diff = abundout(l) - average
         else
            diff = 999.99
         endif
         ew = 1000.*width(l)
         rw = dlog10(width(l)/wave1(l))
         loggf = dlog10(gf(l))
         write (array,3007) wave1(l), atom1(l), e(l,1), loggf,
     .         ew, rw, abundout(l), diff
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,3007) wave1(l), atom1(l), e(l,1), loggf,
     .         ew, rw, abundout(l), diff
      enddo
      write (array,3008) average, deviate, kount
      line = line + 1
      if (errmess(1:9) .ne. 'stopinfo!') call prinfo (line)
      write (nf2out,3008) average, deviate, kount
      if (kount .gt. 2 .and. deltaep .gt. 1.5) then
         write (array,3009) xxm1, xxb1, xxr1 
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,3009) xxm1, xxb1, xxr1
      else
         write (array,*) 'No statistics done for E.P. trends' 
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,*) 'No statistics done for E.P. trends' 
      endif
      if (kount .gt. 2 .and. deltarw .gt. 0.5) then
         write (array,3010) xxm2, xxb2, xxr2   
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,3010) xxm2, xxb2, xxr2   
      else
         write (array,*) 'No statistics done for R.W. trends' 
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,*) 'No statistics done for R.W. trends' 
      endif
      if (kount .gt. 2 .and. deltawv .gt. 500.) then
         write (array,3011) xxm3, xxb3, xxr3
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,3011) xxm3, xxb3, xxr3
      else
         write (array,*) 'No statistics done for wavelength trends'
         if (errmess(1:9) .ne. 'stopinfo!') then
            line = line + 1
            call prinfo (line)
         endif
         write (nf2out,*) 'No statistics done for wavelength trends'
      endif
      return


c*****format statements
1001  format (/'INPUT LINES DATA FOR ' ,i5, ' LINES'/ 
     .        '   #', 5x, 'wave1', 3x, 'spec', 9x, 'spec#', 
     .        3x, 'E.P.', 3x, 'loggf', 5x, 'damp', 4x, 'logSTR', 
     .        5x, 'E.W.')
1002  format (20x, 6x, 'chi1', 4x, 'chi2', 6x, 'chi3', 4x, 'charge',
     .        6x, 'mass', 4x, 'rdmass')
1003  format (i4, f10.6, 2x, a2, a4, f13.5, f7.3, f8.3, 2x, a7,
     .        f9.1, f9.2)
1004  format (i4, f10.3, 2x, a2, a4, f13.5, f7.3, f8.3, 2x, a7,
     .        f10.2, f8.2)
1005  format (20x, f10.3, f8.3, f10.3, f10.1, f10.2, f10.4)
1006  format (a1)
1007  format (1x,a1,a2,i4)
1008  format (2a2,i4)
1009  format (i4, f10.6, 3x, a4, 1x, f13.5, f7.3, f8.3, 2x, a7,
     .        f10.2, f8.2)
1010  format (i4, f10.3, 3x, a4, 1x, f13.5, f7.3, f8.3, 2x, a7,
     .        f10.2, f8.2)
1011  format (/'SYNTHETIC SPECTRUM PARAMETERS (units=1/cm)'/
     .        10x,'start =',f11.3,' ',5x,'stop =',f11.3,' '/  
     .        'step size in the spectrum =',f11.4,' '/
     .        'at each point, opacity will include lines' ,
     .        ' within',f11.4,' of the point')          
1012  format (/'SYNTHETIC SPECTRUM PARAMETERS (units=A)'/
     .        10x,'start =',f11.3,' ',5x,'stop =',f11.3,' '/  
     .        'step size in the spectrum =',f11.3,' '/
     .        'at each point, opacity will include lines' ,
     .        ' within',f11.3,' of the point')          
1013  format (/'CURVE-OF-GROWTH PARAMETERS'/
     .        10x,'log(R.W) lower bound =',f7.3, 
     .        10x,'upper bound =',f7.3/
     .        10x,'step size in the curve =',f7.3)
1014  format (/'PARTITION FUNCTIONS')
1015  format (/'Z =',i2,' (',a2,'),  mass=',f8.3,'  I.P.s=', 3f7.3)
1016  format ('    ionization state = ',i1/(10f8.3))  
2001  format (/'INPUT LINES DATA FOR',i4,' STRONG LINES'/
     .        '   #', 5x, 'wave1', 3x, 'spec', 9x, 'spec#',
     .        3x, 'E.P.', 3x, 'loggf', 5x, 'damp', 4x, 'logSTR')
2004  format ('SPECIES = ', i5, ' IS A MOLECULE, NOT ALLOWED AS A ',
     .        'STRONG LINE; I QUIT!')
3001  format (a80)
3002  format ('Abundance Results for Species ',a2,a4,
     .        '       (input abundance = ',f7.3,')')
3003  format ('wavelength', 9x, 'ID', 6x, 'EP', 3x, 'logGF', 5x, 'EWin',
     .        3x, 'logRWin', 5x, 'abund', 3x, 'delavg')
3004  format ('Abundance Results for Species ',a8,
     .        '       (input abundance = ',f6.2,')')
3005  format ('From these data, the abundance of ',a2, 
     .        ' will be altered')
3006  format ('wavelength', 9x, 'ID', 6x, 'EP', 3x, 'logGF', 5x, 'EWin',
     .        3x, 'logRWin', 5x, 'abund', 3x, 'delavg')
3007  format (f10.3, f11.5, f8.3, f8.3, f9.2, f10.3, f10.3, f9.3)
3008  format ('average abundance = ',f6.3,'     std. ',
     .        'deviation = ',f6.3,'     #lines = ',i3)
3009  format ('E.P. correlation:  slope = ',f7.3,'  intercept = ',
     .        f7.3,'  corr. coeff. = ',f7.3)
3010  format ('R.W. correlation:  slope = ',f7.3,'  intercept = ',
     .        f7.3,'  corr. coeff. = ',f7.3)
3011  format ('wav. correl.:  slope = ',1pd11.3,'  intercept = ',
     .        0pf7.3,'  corr. coeff. = ',f7.3)



      return
      end



