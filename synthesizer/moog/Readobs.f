
      subroutine readobs (line)
c******************************************************************************
c     this routine reads in observed spectra.  the format is indicated by
c     the value of "specfileopt"
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      include 'Obspars.com'
      byte      int1(2880), onebyte
      integer*2 int2(1440)
      integer*4 int4(720)
      real*4 real4(720)
      real*4 xtemp, ytemp
      equivalence (int1,int2,int4,real4)
      character head*2880


c*****first set some initial parameters
      write (obsitle,1021)
      iend = 0
      naxis = 1
      naxis1 = 0
      bscale = 0.
      bzero = 0.
      do j=1,9
         disp(j) = 0.
      enddo


c*****branch to the appropriate data file type
      go to (100,200,300,400,500), specfileopt


c*****here is a pure FITS file type. Reading is done with ordinary
c     FORTRAN read statements
c     first get the header records and search for key parameters
100   irec = 1
101   read (unit=nfobs,rec=irec,err=1002,iostat=ierr) head
      call obshead (head,iend,line)
      if (lount .eq. -1) return
      if (iend .eq. 0) then
         irec = irec + 1
         go to 101
      endif


c     next read the flux array from the file
      nrec = lount/nblock
      ipt = 0
      if (mod(lount,nblock) .ne. 0) nrec = nrec + 1
      do j=1,nrec
         irec = irec + 1
         jpt = min0(nblock,lount-ipt)
         if (ibits .eq. 16) then
            read (unit=nfobs,rec=irec,err=1006,iostat=ierr) int1
            if (byteswap .eq. 1) then
               do k=2,2880,2
                  onebyte = int1(k)
                  int1(k) = int1(k-1)
                  int1(k-1) = onebyte
               enddo
            endif
            do k=1,jpt
               yobs(ipt+k) = bzero + bscale*real(int2(k))
            enddo
         elseif (ibits .eq. 32) then
            read (unit=nfobs,rec=irec,err=1006,iostat=ierr) int1
            if (byteswap .eq. 1) then
               do k=4,2880,4
                  onebyte = int1(k)
                  int1(k) = int1(k-3)
                  int1(k-3) = onebyte
                  onebyte = int1(k-1)
                  int1(k-1) = int1(k-2)
                  int1(k-2) = onebyte
               enddo
            endif
            do k=1,jpt
               yobs(ipt+k) = bzero + bscale*real(int4(k))
            enddo
         elseif(ibits .eq. -32) then
            read (unit=nfobs,rec=irec,err=1006,iostat=ierr) int1
            if (byteswap .eq. 1) then
               do k=4,2880,4
                  onebyte = int1(k)
                  int1(k) = int1(k-3)
                  int1(k-3) = onebyte
                  onebyte = int1(k-1)
                  int1(k-1) = int1(k-2)
                  int1(k-2) = onebyte
               enddo
            endif
            do k=1,jpt
               yobs(ipt+k) = real4(k)
            enddo
         endif
         ipt = ipt + jpt
      enddo     


c     now adjust the continuum, if desired
      do j=1,lount
         yobs(j) = yobs(j)*contnorm
      enddo


c     now fill in the wavelength array
      do j=1,lount
         xobs(j) = sngl(wavecalc(sngl(dfloat(j)),lount,disp))
      enddo
      lim1obs = 1
      lim2obs = lount
      if(xobs(2) .lt. xobs(1)) then
         do j=1,lount/2
            xtemp = xobs(j)
            ytemp = yobs(j)
            xobs(j) = xobs(lount+1-j)
            yobs(j) = yobs(lount+1-j)
            xobs(lount+1-j) = xtemp
            yobs(lount+1-j) = ytemp
         enddo
      endif
      return
         

c     or if trouble, give error messages, and abort
1002  write(array,1003) ierr
      istat = ivwrite(line+2,3,array,39)
      go to 1007
1006  write(array,1004) ierr
      istat = ivwrite(line+2,3,array,39)
1007  lount = -1
      return


c*****this is a dummy for future file use
200   return


c*****this is a dummy for future file use
300   return


c*****this is a dummy for future file use
400   return


c*****here is a MONGO-style input array
500   rewind nfobs
      read (nfobs,5001) obsitle
      i = 1
501   read (nfobs,*,end=525) xobs(i),yobs(i)
      i = i + 1
      go to 501
525   lount = i - 1
      lim1obs = 1
      lim2obs = lount
      if(xobs(2) .lt. xobs(1)) then
         do j=1,lount/2
            xtemp = xobs(j)
            ytemp = yobs(j)
            xobs(j) = xobs(lount+1-j)
            yobs(j) = yobs(lount+1-j)
            xobs(lount+1-j) = xtemp
            yobs(lount+1-j) = ytemp
         enddo
      endif


c     now adjust the continuum, if desired
      do j=1,lount
         yobs(j) = yobs(j)*contnorm
      enddo
      return


c*****format statements
1003  format('ERROR IN FITS HEADER READ:   ERROR=',i4)
1004  format('ERROR IN DATA READ:  ERROR=',i4,'/')
1021  format (80(' '))
5001  format (a80)


      end



