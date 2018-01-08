
      subroutine weedout                   
c******************************************************************************
c     This routine goes through an initial line list and culls from it
c     absorption lines that are not substantial contributors.  This is
c     done in a simple fashion by eliminating lines weaker than X, where
c     X = kapnu/kaplam at the approximate line wavelength, calculated
c     at a continuue optical depth of 0.5.  The user will be prompted for
c     the desired value of X.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      real*8 xratio
      real*4 shortnum


c*****examine the parameter file
      call params


c*****open the files for: standard output, lines to be kept, and lines
c     to be heaved out
      nf1out = 20     
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf8out = 21               
      lscreen = lscreen + 2
      array = 'KEPT LINES'
      nchars = 10
      call infile ('output ',nf8out,'formatted  ',0,nchars,
     .             f8out,lscreen)
      nf9out = 22               
      lscreen = lscreen + 2
      array = 'DISCARDED LINES'
      nchars = 15
      call infile ('output ',nf9out,'formatted  ',0,nchars,
     .             f9out,lscreen)


c*****open and read the model atmosphere file
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel
      call eqlib


c*****open the initial line list file
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE INITIAL LINE LIST'
      nchars = 21
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)

c*****ask the user about the value of X
      lscreen = lscreen + 2
      array = 'GIVE THE MINIMUM LINE/CONTINUUM OPACITY RATIO TO KEEP: '
      nchars = 55
      call getnum (nchars,lscreen,xratio,shortnum)
      write (nf1out,1001) xratio
      write (nf8out,1002)
      write (nf9out,1003)
      

c*****compute the line opacities
      call inlines (3)
1     call nearly (1)


c*****calculate continuum quantities at the line list wavelength middle
      wave = (wave1(1)+wave1(nlines))/2.
      call opacit (2,wave)
      if (modprintopt .ge. 2) 
     .   write(nf1out,1006) wave,(kaplam(i),i=1,ntau)


c*****divide the lines into keepers and discards
      do j=1,nlines
         residual = 10.*atom1(j) - dble(nint(10.*(atom1(j))))
         if (strength(j)/kaplam(jtau5) .ge. xratio) then
            if (atom1(j) .lt. 100.) then
               if (residual .gt. 0. .and. dampnum(j) .gt. 0.) then
                  write (nf8out,1007) wave1(j), atom1(j), e(j,1), 
     .            dlog10(gf(j)), dlog10(dampnum(j)), dlog10(strength(j))
               else if (residual .gt. 0. .and. dampnum(j) .le. 0.) then
                  write (nf8out,1007) wave1(j), atom1(j), e(j,1), 
     .            dlog10(gf(j)), 0.0, dlog10(strength(j))
               else if (residual .le. 0. .and. dampnum(j) .gt. 0.) then
                  write (nf8out,1004) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), dlog10(dampnum(j)), dlog10(strength(j))
               else 
                  write (nf8out,1004) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), 0.0, dlog10(strength(j))
               endif
            else 
               if (residual .gt. 0.) then
                  write (nf8out, 1008) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), d0(j), dlog10(strength(j))
               else
                  write (nf8out, 1005) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), d0(j), dlog10(strength(j))
               endif
            endif
         else
            if (atom1(j) .lt. 100.) then
               if (residual .gt. 0. .and. dampnum(j) .gt. 0.) then
                  write (nf9out,1007) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), dlog10(dampnum(j)), dlog10(strength(j))
               else if (residual .gt. 0. .and. dampnum(j) .le. 0.) then
                  write (nf9out,1007) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), 0.0, dlog10(strength(j))
               else if (residual .le. 0. .and. dampnum(j) .gt. 0.) then
                  write (nf9out,1004) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), dlog10(dampnum(j)), dlog10(strength(j))
               else
                  write (nf9out,1004) wave1(j), atom1(j), e(j,1),
     .            dlog10(gf(j)), 0.0, dlog10(strength(j))
               endif
            else
               if (residual .gt. 0.) then
                  write (nf9out, 1008) wave1(j), atom1(j), e(j,1),
     .                     dlog10(gf(j)), d0(j), dlog10(strength(j))
               else
                  write (nf9out, 1005) wave1(j), atom1(j), e(j,1),
     .                     dlog10(gf(j)), d0(j), dlog10(strength(j))
               endif
            endif
         endif
      enddo
      if (nlines +nstrong .eq. 2500) then
         call inlines (6)
         go to 1
      endif
         

c*****finish
      call finish (0)


c*****format statements
1001  format (/'DESIRED LINE-TO-CONTINUUM MINIMUM OPACITY RATIO: ', 
     .        1pe10.2)
1002  format ('THIS IS THE KEEPER LINE LIST')
1003  format ('THIS IS THE DISCARDED LINE LIST')
1004  format (f10.4, f10.1, f10.3, f10.3, f10.3, 20x, f9.1)
1005  format (f10.4, f10.1, f10.3, f10.3, 10x, f10.3, 10x, f9.1)
1006  format ('  kaplam from 1 to ntau at wavelength',f10.2/
     .        (6(1pd12.4)))
1007  format (f10.4, f10.4, f10.3, f10.3, f10.3, 20x, f9.1)
1008  format (f10.4, f10.5, f10.3, f10.3, 10x, f10.3, 10x, f9.1)


      end


