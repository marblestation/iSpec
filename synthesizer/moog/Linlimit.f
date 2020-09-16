
      subroutine linlimit 
c******************************************************************************
c     This routine marks the range of lines to be considered in a 
c     particular line calculations, depending on the type of calculation
c     (e.g. synthetic spectrum, single line curve-of-growth, etc.)
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Linex.com'
      include 'Pstuff.com'


      lineflag = 0
c*****for single-line computations, the line rage is the whole line set;
c     this will be called from "ewfind"
      if     (mode .eq. 1) then
         lim1line = 1
         lim2line = nlines
          return


c*****for deriving abundances from single of lines of one species, delimit
c     the lines of that species as the range; called from "abfind"
      elseif (mode .eq. 2) then
         if (lim2line .eq. nlines) then
            if (nlines .eq. 1) then
               lim1line = 1
               lim2line = 1
               mode = -1
            else
               lim1line = -1
            endif
            return
         endif
         if (lim1line .eq. 0) then
            lim1line = 1 
         else
            lim1line = lim2line + 1
         endif
         if (lim1line .eq. nlines) then
            lim2line = lim1line
            return
         else
            oldatom = atom1(lim1line)
            do j=lim1line+1,nlines
               if (atom1(j) .ne. oldatom) then
                  lim2line = j - 1
                  return
               endif
            enddo
         endif
         lim2line = nlines
         return


c*****for spectrum synthesis, find the range of lines to include at each
c     wavelength step; called from "synspec"; if requested synthesis
c     begins more than 10A blueward or goes on more than 105A redward of the
c     linelist limits, the synthesis aborts with a message;
      elseif (mode .eq. 3) then
         wavelo = wave - delta
         wavehi = wave + delta
c        requested synthesis too far from linelist limits
         if     (wavehi .lt. wave1(1)-10.0) then
            write (*,1004)
            stop
         elseif (wavelo .gt. wave1(nlines)+10.0 .and.
     .           nlines+nstrong .lt. 2500) then
            write (*,1005)
            stop
         endif
c        blank synthesis at start or end of requested wavelength range
         if     (wavehi .lt. wave1(1)) then
            lim1line = 1
            lim2line = 1
            lineflag = -1
            return
         elseif (wavelo .gt. wave1(nlines)) then
            lim1line = nlines
            lim2line = nlines
            lineflag = -1
            return
         endif
c        requested synthesis region is at least partially within the line list
c        first set the lower synthesis limit; at the beginning of a synthesis
c        of with a new chunk of lines, lim1line = 0
         if (lim1line.eq.0 .or. wavelo .lt. wave1(1)) lim1line = 1
         do j=lim1line,nlines
            if (wavelo .lt. wave1(j)) then
               lim1line = j
               exit
            endif
         enddo
c        now set the upper synthesis limit
         do j=lim1line,nlines
            if (wavehi .lt. wave1(j)) then
               lim2line = j - 1
               if (lim1line.eq.lim2line .and.
     .             wavelo.gt. wave1(lim1line)) then
                  lineflag = -1
               endif
               return
            endif
         enddo   
c        here the end of the line list has been reached; decide whether
c        to read in more lines
         if (nlines+nstrong .eq. 2500) then
            lim2line = -1
         else
            lim2line = nlines
         endif
         return


c*****for blended line force fits to the EWs, the range is a set of
c     lines in a particular blend
      elseif (mode .eq. 4) then
         if (lim1line .eq. 0) lim1line = 1
         if (group(lim1line) .ne. 0.) then
            write (array,1001)
            call prinfo (10)
            write (array,1002)
            call prinfo (11)
            stop
         endif
         if (lim1line .eq. nlines) then
            lim2line = lim1line
         else
            do j=lim1line+1,nlines
               if (group(j) .ne. 1) then
                  lim2line = j - 1
                  return
               endif
            enddo
            lim2line = nlines
         endif
         return
      else
         write (array,1003) mode
         call prinfo (10)
         stop
      endif


c*****format statements
1001  format ('TROUBLE! THE FIRST LINE IN THE GROUP')
1002  format ('DOES NOT DEFINE A NEW GROUP OF LINES!  I QUIT!')
1003  format ('mode = ', i3, ' BUT IT MUST BE 1-4 ONLY; I QUIT!')
1004  format ('synthesis begins more than 10A blueward of the',
     .        ' linelist start; I QUIT!')
1005  format ('synthesis ends more than 10A redward of the',
     .        ' linelist end; I QUIT!')


      end

