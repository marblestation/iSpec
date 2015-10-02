
      subroutine infile (type,iunit,mode,irec,charcount,fname,line)
c******************************************************************************
c     this routine serves to open files for reading, and does minimal
c     error checking.
c******************************************************************************

      include 'Atmos.com'
      include 'Pstuff.com'
      integer charcount
      character type*7,kstat*7,yesno*1,mode*11
      character fname*80

c  decide on the file status desired
      jstat = 0
      if     (type .eq. 'input  ') then
         kstat = 'old    '
      elseif (type .eq. 'output ') then
         kstat = 'new    '
      elseif (type .eq. 'overout') then
         kstat = 'unknown'
      endif

c  write out the appropriate message about this file
5     nchars = charcount
      if (fname .eq. 'optional_output_file') then
         return
      elseif (fname .eq. 'no_filename_given') then
         array(charcount+1:charcount+24) ='; what is the filename? '
         charcount = charcount + 24
         call getasci (charcount,line)
         fname = chinfo
      else
         array(charcount+1:charcount+24) ='; here is the filename: '
         array(charcount+25:79) = fname
         charcount = 79
         call putasci (charcount,line)
         if (type .ne. 'input  ') kstat = 'unknown'
      endif
     
c  open the file specified by the user, earlier or now
6     if (mode .eq. 'formatted  ') then
         open (unit=iunit,file=fname,access='sequential',
     .         form=mode,blank='null',status=kstat,
     .         iostat=jstat,err=10)
      else
         open (unit=iunit,file=fname,access='direct',
     .         form=mode,status=kstat,recl=irec,
     .         iostat=jstat,err=10)
      endif
      istat = ivmove (line+1,1)
      istat = ivcleol ()
      return

c  here are the file reading error messages;
c  if an expected file is not found, 118 is the error code for SunOS, 1018
c  is for Solaris, and 2 is for Redhat Linux operating systems.
10    if (jstat .eq. 118 .or. jstat .eq. 1018 .or.
     .    jstat .eq. 2) then
         write (errmess,1001) jstat
         istat = ivwrite (line+2,3,errmess,44)
         fname = 'no_filename_given'
         charcount = nchars
         go to 5
c  if a file is in danger of being over-written, 117 is the error code for 
c  SunOS, 1017 is for Solaris, and 128 is for Redhat Linux operating systems.
      elseif (jstat .eq. 117 .or. jstat .eq. 1017 .or.
     .        jstat .eq. 128) then
         write (errmess,1002) jstat
         istat = ivwrite (line+2,3,errmess,41)
         read (*,*) yesno
         if (yesno .eq. 'y') then
            kstat = 'unknown'
            go to 6
         else
            write (errmess,1003) 
            istat = ivwrite (line+2,3,errmess,32)
            fname = 'no_filename_given'
            charcount = nchars
            go to 5
         endif
      else
         write (errmess,1004) jstat
         istat = ivwrite (line+2,3,errmess,46)
      endif


c*****format statements
1001  format ('ERROR ',i4,': FILE DOES NOT EXIST! TRY AGAIN.') 
1002  format ('ERROR ',i4,': FILE EXISTS! OVERWRITE (y/n)?')
1003  format ('PLEASE GIVE A DIFFERENT FILENAME')
1004  format ('ERROR ',i4,': UNKNOWN FILE READING ERROR; ABORT!')


      end




