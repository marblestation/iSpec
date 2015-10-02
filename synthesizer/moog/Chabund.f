
      subroutine chabund
c*****************************************************************************
c  this routine changes the input abundances according to the users request
c  and the re-runs the synthesis with the new abundances OVERWRITING the 
c  original output files, but using the same model, linelist, etc
c*****************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Pstuff.com'
      include 'Multistar.com'
      real*8 xnum
      real*4 shortnum
      character choice2


c*****set the local temporary parameters to what they were for the last synth
      newnumiso = numiso
      newnumisosyn = numisosyn
      newnumpecatom = numpecatom
      newnumatomsyn = numatomsyn
      if (newnumatomsyn .eq. 0) newnumatomsyn = 1
      do i=3,95
         newpec(i) = pec(i)
         do j=1,newnumatomsyn
            newpecabund(i,j) = pecabund(i,j)
         enddo
      enddo
      do i=1,newnumiso
         if(newnumiso .ge. 1) then
            newisotope(i)=isotope(i)
            do j=1,newnumisosyn
               newisoabund(i,j)=isoabund(i,j)
            enddo
         endif
      enddo
      if (newnumatomsyn .le. 0) newnumatomsyn = 1


c*****SPECIAL CASE
c     for the "binary" driver, present a limited set of options; only
c     the abundances of the elements already chosen for variation in
c     the parameter file can be altered
11    if (control .eq. 'binary ') then
         istat = ivcleof(4,1)
         if (syncount .eq. 1) then
            array = 'Abundance Alterations for the PRIMARY STAR:   '
         else
            array = 'Abundance Alterations for the SECONDARY STAR: '
         endif
         istat = ivwrite(5,5,array,46)
         array = 'Atom   Abundance (relative to initial model)'
         istat = ivwrite(6,1,array,49)
         line = 6
         do i=3,95
            if(newpec(i) .eq. 1) then
               line = line + 1
               write (array,1001) i,(newpecabund(i,j),j=1,newnumatomsyn)
               istat = ivwrite(line,1,array,57)
            endif
         enddo
         array = 'Options: c = change abundances '
         istat = ivwrite(line+1,1,array,41)
         array = '         x = exit no changes    q = rerun synth'
         istat = ivwrite(line+2,1,array,47)
         array = 'What is your choice? '
         nchars = 21
         call getasci (nchars,line+4)
         choice = chinfo(1:1)
         if (choice.ne.'c' .and. choice.ne.'x' .and.
     .       choice.ne.'q') go to 11
      endif
         
        
c*****for the "synth" driver, present the user with the current abundance 
c     alterations, and the options; find out what is desired
1     if (control .eq. 'synth'  .or. 
     .    control .eq. 'isotop' .or.
     .    control .eq. 'synpop') then
         line = 4
         istat = ivcleof(line,1)
         write (array,1006)
         line = line + 1
         istat = ivwrite(6,1,array,72)
         line = line + 1
         do i=3,95
            if(newpec(i) .eq. 1) then
               line = line + 1
               write (array,1001) i,(newpecabund(i,j),j=1,newnumatomsyn)
               istat = ivwrite(line,1,array,57)
            endif
         enddo
         if (newnumiso .gt. 0) then
            do i=1,newnumiso
               line = line + 1
               write (array,1002) i, newisotope(i), 
     .                            (newisoabund(i,j),j=1,newnumisosyn)
               istat = ivwrite(line,1,array,68)
            enddo
         endif
         write (array,1007)
         istat = ivwrite(line+1,1,array,61)
         write (array,1008)
         istat = ivwrite(line+2,1,array,72)
         array = 'What is your choice? '
         nchars = 21
         call getasci (nchars,line+3)
         choice = chinfo(1:1)
         if (choice.ne.'c' .and. choice.ne.'i' .and.
     .       choice.ne.'n' .and. choice.ne.'x' .and.
     .       choice.ne.'q') go to 1
      endif

         
c*****for "synth", change elemental abundances
      if     (choice .eq. 'c') then
20       istat = ivcleof(line+1,1)
         array = 'Which element to change? '
         nchars = 25
         call getnum (nchars,line+1,xnum,shortnum)
         if (xnum.lt.2.0 .or. xnum.gt.95.) go to 20
         j = nint(xnum)
         if (newpec(j) .eq. 1) then
25          array = 'n = new abundances, or z = zero offsets? '
            nchars = 41
            call getasci (nchars,line+2)
            choice2 = chinfo(1:1)
            if     (choice2 .eq. 'z') then 
               newpec(j) = 0
               do i=1,5
                  newpecabund(j,i) = 0.0
               enddo
               newnumpecatom = newnumpecatom - 1
            elseif (choice2 .eq. 'n') then
               write (array,1004)
               istat = ivwrite(line+3,1,array,40)
               read (*,*) (newpecabund(j,i),i=1,newnumatomsyn)
            else
               go to 25
            endif
         else
            newpec(j) = 1
            write (array,1004)
            istat = ivwrite(line+2,1,array,40)
            read (*,*) (newpecabund(j,i),i=1,newnumatomsyn) 
            newnumpecatom = newnumpecatom + 1
         endif


c*****for "synth", change isotopic factors
      elseif (choice .eq. 'i') then
         istat = ivcleof(line+1,1)
         array = 'Options: c = change an isotopic factor'
         istat = ivwrite(line+1,1,array,49)
         array = '         n = enter a new isotope'
         istat = ivwrite(line+2,1,array,48)
30       array = 'What is your choice? '
         nchars = 22
         call getasci (nchars,line+3)
         choice2 = chinfo(1:1)
         if (choice2 .eq. 'c') then
35          istat = ivcleof(line+1,1)
            array =  'Which isotope number from the list? '
            nchars = 36
            call getnum (nchars,line+1,xnum,shortnum)
            j = nint(xnum)
            if (j.lt.1 .or. j.gt.newnumiso) go to 35
            istat = ivcleof(line+2,1)
            array = 'What are the new division factors? '
            istat = ivwrite(line+2,1,array,35)
            read (*,*) (newisoabund(j,i),i=1,newnumisosyn)
         elseif (choice2 .eq. 'n') then
            newnumiso = newnumiso + 1
            istat = ivcleof(line+1,1)
            array = 'What is the new isotope designation? '
            nchars = 37
            call getnum (nchars,line+4,xnum,shortnum)
            newisotope(newnumiso) = xnum
            array = 'What are its division factors? '
            istat = ivwrite(line+2,1,array,31)
            read (*,*) (newisoabund(newnumiso,i),i=1,newnumisosyn)
         else
            go to 30
         endif
            

c*****for "synth", change the number of syntheses
      elseif (choice .eq. 'n') then
55       array = 'How many synths? '
         nchars = 17
         call getnum (nchars,line+5,xnum,shortnum)
         if (xnum .gt. 5.) go to 55
         newnumatomsyn = nint(xnum)
         go to 1


c*****for "synth", exit the routine without changing anything
      elseif (choice .eq. 'x') then
         return


c*****for "synth", make the proposed alterations permanent; then
c     return to the calling routine, which allegedly will
c     redo the syntheses.
      elseif (choice .eq. 'q') then
         numiso = newnumiso
         numisosyn = newnumisosyn
         numpecatom = newnumpecatom
         numatomsyn = newnumatomsyn
         do i=3,95
            pec(i) = newpec(i)
            do j=1,numatomsyn
               pecabund(i,j) = newpecabund(i,j)
            enddo
         enddo
         do i=1,numiso
            if(numiso .ge. 1) then
               isotope(i) = newisotope(i)
               do j=1,numisosyn
                  isoabund(i,j)=newisoabund(i,j)
               enddo
            endif
         enddo
         return
      endif


c*****loop back and print out the main menu again
      if (control .eq. 'synth  ') then
         go to 1
      else
         go to 11
      endif


c*****format statements
1001  format (i3, 4x, 5(f8.3, 2x))
1002  format (i2, 2x, f10.5, 5(3x, f8.3))
1004  format ('Enter the new offsets on the line below:')
1006     format ('element, abundance offsets   OR   ',
     .            'isotope number, isotope name, factors')
1007     format ('Options: c = change abundance       ',
     .           'i = change isotopic ratio')
1008     format('         n = change # syntheses     q = rerun ',
     .          'syntheses        x = exit')

      end



