      subroutine read_segments(iunit,segmentsfile,
     &             nsegmax,nsegment,xlsegmin,xlsegmax)
!
! read wavelengths segments to be calculated
! BPz 03/07-2020
!
      implicit none
      integer nsegment,nsegmax,i,iunit
      doubleprecision xlsegmin(*),xlsegmax(*)
      character segmentsfile*256,header*256
      logical head
!
      open(iunit,file=segmentsfile,status='old')
      head = .true.
      do while (head)
        read(iunit,10) header
        if (header(1:1).eq.';') then
          head=.true.
        else
          head=.false.
        endif
      enddo
      backspace(iunit)
      do i=1,nsegmax+1
        read(iunit,*,end=99) xlsegmin(i),xlsegmax(i)
        nsegment=i
      enddo
      print*,'number of segments too high!'
      print*,'increase nsegmax!'
      stop ' stop in read_segments'
99    continue
      print*,'Read ',nsegment,' wavelengths segments'
      return
10    format(a)
      end
