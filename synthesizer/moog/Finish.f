
      subroutine finish (number)
c******************************************************************************
c     This routine simply wraps up MOOG
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'


c  close the files
      if (nfmodel .ne. 0)      close (unit=nfmodel)
      if (nflines .ne. 0)      close (unit=nflines)
      if (nfslines .ne. 0)     close (unit=nfslines)
      if (nftable .ne. 0)      close (unit=nftable)
      if (nfobs .ne. 0)        close (unit=nfobs)
      if (nf1out .ne. 0)       close (unit=nf1out)
      if (nf2out .ne. 0)       close (unit=nf2out)
      if (nf3out .ne. 0)       close (unit=nf3out)
      if (nf4out .ne. 0)       close (unit=nf4out)
      if (nf5out .ne. 0)       close (unit=nf5out)
      if (control .ne. 'gridsyn' .and. control .ne. 'gridplo') then
         if (nf6out .ne. 0)    close (unit=nf6out)
         if (nf7out .ne. 0)    close (unit=nf7out)
         if (nf8out .ne. 0)    close (unit=nf8out)
         if (nf9out .ne. 0)    close (unit=nf9out)
         if (nf10out .ne. 0)   close (unit=nf10out)
      endif

c  write the closing message
      if (number .eq. 0) then
         istat = ivcleof (4,1)
         write (array,1001) 
         istat = ivwrite (5,1,array,79)
      endif
      return


c*****format statements
1001  format (22('<'),10x,'MOOG HAS ENDED!',10x,22('>'))


      end
      

