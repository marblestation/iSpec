
      subroutine obshead (head,iend,line)
c******************************************************************************
c     this routine decodes the header records of the observed FITS
c     spectrum file
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Pstuff.com'
      include 'Obspars.com'
      include 'Atmos.com'
      character head*2880

      do j=1,36
         k = 80*(j-1)
         if     (head(k+1:k+8) .eq. 'SIMPLE  ') then
            if (head(k+30:k+30) .ne. 'T') then
               write(array,1029) head(k+1:k+58)
               istat = ivwrite (line+2,3,array,79)
               go to 1007
            endif
         elseif (head(k+1:k+8) .eq. 'BITPIX  ') then
               read (head(k+1:k+80),1025) ibits
               if (ibits .eq. 16) then
                  nblock = 1440
               elseif (ibits .eq. 32) then
                  nblock = 720
               elseif (ibits .eq. -32) then
                  nblock = 720
               else
                  write(array,1026) ibits
                  istat = ivwrite (line+2,3,array,32)
                  go to 1007
               endif
         elseif (head(k+1:k+8) .eq. 'NAXIS   ') then
               read (head(k+1:k+80),1025) naxis
               if (naxis .ne. 1) then
                  write(array,1028) head(k+1:k+58) 
                  go to 1007
               endif
         elseif (head(k+1:k+8) .eq. 'NAXIS1  ') then
               read (head(k+1:k+80),1025) lount
         elseif (head(k+1:k+8) .eq. 'OBJECT  ') then
               write (obsitle,1027) head(k+12:k+80)
         elseif (head(k+1:k+8) .eq. 'BZERO   ') then
               read (head(k+1:k+80),1024) bzero
         elseif (head(k+1:k+8) .eq. 'BSCALE  ') then
               read (head(k+1:k+80),1024) bscale
         elseif ((head(k+1:k+8) .eq. 'W0      ') .or. 
     .      (head(k+1:k+8) .eq. 'CRVAL1  ')) then
               read (head(k+1:k+80),1024) disp(1)
         elseif ((head(k+1:k+8) .eq. 'WPC     ') .or. 
     .      (head(k+1:k+8) .eq. 'CDELT1  ')) then
               read (head(k+1:k+80),1024) dval
               if (dval .ne. 1.) disp(2) = dval
         elseif (head(k+1:k+8) .eq. 'CD1_1    ') then
               read (head(k+1:k+80),1024) disp(2)
         elseif (head(k+1:k+8) .eq. 'FILENAME') then
               write (obsitle(39:80),1023) head(k+12:k+53)
         elseif (head(k+1:k+8) .eq. 'HISTORY ') then
               if (head(k+24:k+28) .eq. 'DISP=') then
                  read (head(k+1:k+80),1022) (disp(i),i=1,4)
               elseif (head(k+20:k+26) .eq. 'D1,2,3:') then
                  read (head(k+1:k+80),1042) (disp(i),i=1,3)
               elseif (head(k+20:k+26) .eq. 'D4,5,6:') then
                  read (head(k+1:k+80),1042) (disp(i),i=4,6)
               elseif (head(k+20:k+26) .eq. 'D7,8,9:') then
                  read (head(k+1:k+80),1042) (disp(i),i=7,9)
                  if (disp(7).ne.0.0 .and. disp(8).eq.0.0 .and.
     .                disp(9).eq.0.0) then
                     disp(8) = 1.0
                     disp(9) = lount
                  endif
               endif
            elseif (head(k+1:k+8) .eq. 'END     ') then
               iend = 1
               return
         endif
      enddo
      return
         

1007  lount = -1
      return


c*****format statements
1022  format(28x,1p4d13.5)
1023  format (a41)
1024  format (10x,d20.10)
1025  format (10x,i20)
1026  format('SORRY: I CANT HANDLE BITPIX=',i4)
1027  format (a68)
1028  format('ILLEGAL NAXIS ENTRY: ',a58)
1029  format('ILLEGAL FILE FORMAT: ',a58)
1042  format(26x,1p3d18.11)


      end



