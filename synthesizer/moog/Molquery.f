
      subroutine molquery
c******************************************************************************
c     This routine discovers whether a species being analyzed in EW 
c     predicting of force-fitting mode is involved in molecular equilibrium 
c     (=M.E.) calculations as well; the variables is set appropriately
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      include 'Mol.com'
      include 'Pstuff.com'

      molflag = 0
      iabatom = int(atom1(lim1obs)+0.0001)


c*****the species is an atom:
c*****search for it in the list of elements done in M.E.
      if (atom1(lim1line) .lt. 100.) then
         if (neq .eq. 0) then
            return
         else
            do n=1,neq
               if (iabatom .eq. iorder(n)) molflag = 1
               return
            enddo
         endif


c*****the species is a molecule:
c*****halt if M.E. wasn't done or didn't include this species
      else
         if (neq .eq. 0) then
            lscreen = lscreen + 2
            write (array,1001) iabatom
            call prinfo (lscreen)
            stop
         endif
         call sunder(atom1(lim1obs),ia,ib)
         iaa = ia
         ibb = ib
         do n=1,neq
            if (ia.eq.iorder(n) .or. ib.eq.iorder(n)) molflag = ia
         enddo
         if (molflag .eq. 0) then
            lscreen = lscreen + 2
            write (array,1002) iabatom
            call prinfo (lscreen)
            stop
         endif
         molflag = 1

c*****if molecule is a hydride, the non-H element will be varied
         if (ia.eq.1 .or. ib.eq.1) then
            if (ia .eq. 1) then
               iabatom = ib
            else
               iabatom = ia
            endif


c*****for other molecules, the user specifies which element will be varied
         else
            write (array,1003) iabatom
            nchars = 56
            call getnum (nchars,ikount+1,xnum,shortnum)
            iabatom = int(xnum+0.0001)
         
            if (iabatom.ne.ia .and. iabatom.ne.ib) then
               write (array,1003)
               stop
            endif
         endif
         return
      endif


c*****format statements
1001  format ('YOU FORGOT TO DO MOLECULAR EQUILIBRIUM FOR ',
     .        'SPECIES ', i3, '; I QUIT!')
1002  format ('YOUR MOLECULAR EQUILIBRIUM DOES NOT INCLUDE ',
     .        'THE ATOMS FOR SPECIES ', i3, '; I QUIT!')
1003  format ('LINES OF SPECIES ', i3, ' ARE NEXT: ',
     .        '   CHANGE ATOMIC NUMBER? ')


      end

