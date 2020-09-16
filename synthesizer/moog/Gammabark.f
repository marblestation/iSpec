
      subroutine gammabark
c******************************************************************************
c     This subroutine pulls in damping factors from Barklem data
c     So far, these data have been compiled from:
c                  1. Barklem, P. S., Piskunov, N., & O'Mara, B. J. 2000, 
c                     A&ApS, 142, 467 for (mostly) neutral species
c                  2. Barklem, P. S., & Aspelund-Johansson, J. 2005, A&Ap, 
c                     435, 373 for Fe II lines with E_lower < 70000 cm-1
c  Added column to Barklem.dat for radiative damping, Gamma_rad.  -AMcW 12/2013
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Pstuff.com'
      include 'Atmos.com'
      include 'Linex.com'
      include 'Dampdat.com'
      data firstread/0/
      character*80 line


c*****on first entry to this routine, read damping data from either 
c    'Barklem.dat' or 'BarklemUV.dat', depending on the wavelength region 
c    of the linelist; read Van der Waals params, if radiative damping 
c    data are present
      if (firstread .eq. 0) then
         if (wave1(nlines) .gt. 3000.) then
            nwant = 35
         else
            nwant = 36
         endif
         do k=1,30000
            call blankstring (line)
            read (nwant,1001,end=10) line 
            read (line,*) wavebk(k), idbk(k), gammabk(k), alphabk(k)
            if (line(34:) .ne. '    ') then
               read(line(34:),*) gammarad(k)
            else
               gammarad(k) = 0.0
            endif
         enddo
10       numbark = k -1
         firstread = 1
      endif
     

c*****identify the Barklem list positions of the wavelength limits of
c     the input line list
      wavemin = 10000000.
      do j=1,nlines+nstrong
         if (wave1(j) .lt. wavemin) wavemin = wave1(j)
      enddo
      wavemax = 0.
      do j=1,nlines+nstrong
         if (wave1(j) .gt. wavemax) wavemax = wave1(j)
      enddo
      do k=1,numbark
         if (wavemin-wavebk(k) .lt. 1.0) then
            nummin = k
            exit
         endif
      enddo
      do k=nummin,numbark
         if (wavebk(k)-wavemax .gt. 1.0) then
            nummax = k
            exit
         endif
      enddo


c*****search for Barklem data
       do j=1,nlines+nstrong
         gambark(j) = -1.
         alpbark(j) = -1.
         gamrad(j)  = -1.
         if (atom1(j) .gt. 100.) cycle
         iatom10 = nint(10.*atom1(j))
         do k=nummin,nummax
            waveerror = -(wave1(j) - wavebk(k))/wavebk(k)
            iii = nint(10.*idbk(k))
            if (dabs(waveerror).lt.5.0d-06 .and.
     .          iii .eq. iatom10) then
               gamrad(j)  = gammarad(k)
               gambark(j) = 10.**gammabk(k)
               alpbark(j) = (1.-alphabk(k))/2.
               exit
            endif
            if (waveerror .gt. 5.0d-06) exit
         enddo
      enddo


c*****exit normally
      return


c*****format statements
1001  format (a80)
      end




