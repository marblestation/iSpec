
      subroutine begin
c***************************************************************************
c     This routine simply starts up MOOG
c     THIS VERSION IS FOR LINUX REDHAT MACHINES
c***************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'
      character*80 line, systemcall
      integer num


c*****define the number of text screen lines for silent mode;
c     this number is hardwired, since it is not really needed at run time.
      if (silent .eq. 'y') then
         maxline = 24
         go to 10
      endif


c*****define the number of lines available on the text screen for
c     interactive mode; this number is discovered from the "stty -a"
c     command, for which the output format is unique to the operating
c     system.
      write (systemcall,*) 'stty -a > tmpsize'
      call system (systemcall)
      open (99,file='tmpsize')
5     read (99,1010,end=15) line
      do i=1,77
         if (line(i:i+3) .eq. 'rows') then
            if     (machine .eq. 'pcl') then
            read (line(i+4:i+6),1011) maxline
            elseif (machine .eq. 'mac') then
               read (line(i-4:i-2),1011) maxline
            elseif (machine .eq. 'uni') then
               read (line(i+6:i+8),1011) maxline
            endif
            go to 10
         endif
      enddo
      go to 5
15    array = 'SCREEN ROW COUNT UNKNOWN; USE 24 ([y]/n)? '
      nchars = 42
      ikount = 2
      call getasci (nchars,ikount)
      choice = chinfo(1:1)
      if (choice.eq.'y' .or. nchars.le.0) then
         go to 10
      else
         call finish (0)
      endif
10    close (99,status='delete')
      write (systemcall,*) '\\rm -f tmpsize'
      call system (systemcall)
      if (maxline .lt. 10) then
         maxline = 24
      else
         maxline = maxline - 2
      endif


c*****clear the text screen
      write (systemcall,*) 'clear'
      call system (systemcall)


c*****open data files carried with the source code: Barklem damping
      nfbarklem = 35
      num = 60
      call getcount (num,moogpath)
      if (moogpath(num:num) .ne. '/') then
         num = num + 1
         moogpath(num:num) = '/'
      endif
      fbarklem(1:num) = moogpath(1:num)
      fbarklem(num+1:num+11) = 'Barklem.dat'
      open (nfbarklem,file=fbarklem)


c*****open data files carried with the source code: Barklem UV damping
      nfbarklemUV = 36
      num = 60
      call getcount (num,moogpath)
      if (moogpath(num:num) .ne. '/') then
         num = num + 1
         moogpath(num:num) = '/'
      endif
      fbarklemUV(1:num) = moogpath(1:num)
      fbarklemUV(num+1:num+13) = 'BarklemUV.dat'
      open (nfbarklemUV,file=fbarklemUV)
 

c  write a header and find the appropriate parameter file, and exit normally
      write (array,1001)
      istat = ivwrite (1,1,array,79)
      write (array,1004)
      istat = ivwrite (2,1,array,79)
      array = 'MOOG PARAMETERS? ' 
      nchars = 15
      nfparam = 50     
      lscreen = 4
      if (silent .eq. 'y') then
         fparam = 'batch.par'
      else
         fparam = 'no_filename_given'     
      endif
      call infile ('input  ',nfparam,'formatted  ',0,nchars,
     .             fparam,lscreen)
      read (nfparam,1002) control
      write (array,1003) control
      istat = ivwrite (2,1,array,58)
      write (array,1001)
      istat = ivwrite (3,1,array,79)
      return


c*****format statements
1001  format (79('*'))
1002  format (a7)
1003  format (22x,'MOOG IS CONTROLLED BY DRIVER ',a7)
1004  format (25(' '),'MOOG LTE VERSION (FEB 2017)',26(' '))   
1010  format (a80)
1011  format (i3)


      end
      







