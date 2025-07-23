!
      subroutine read_departure(iunit,departurefile,maxlevel,modnlevel,
     &                   ndepth,ndepth_read,taumod,b_departure)
!
! read NLTE departure coefficients for one model atom, and one atmosphere model
! created by BPlez 17/04-2020
!
! unit number to open, file name, maximum number of levels, number of levels to read, 
! maximum number of depth points, number of depth points found in file,
! optical depth, departure coefficients
!
      implicit none
      character departurefile*256,oneline*256
      integer iunit,modnlevel,i,ndepth_read,j,ndepth,maxlevel
      real taumod(ndepth),b_departure(maxlevel,ndepth)
      logical header

      open(iunit,file=departurefile,status='old')
      header=.true.
      do while (header)
        read(iunit,10,end=99) oneline
        if (oneline(1:1).ne.'!') then
          header=.false.
          backspace(iunit)
        endif
      enddo
      i=1
      do while (.true.)
        if (i.gt.ndepth) then
          print*,' more depth points in departure file than can ',
     &           'be accomodated in read_departure.f !'
          stop 'Stopping here!'
        endif
        read(iunit,*,end=99) taumod(i),(b_departure(j,i),j=1,modnlevel)
        taumod(i)=10.**taumod(i)
        i=i+1
      enddo

10    format(a)
99    ndepth_read=i-1
      print*,'read departure coefficients for ',i,
     &       ' depths in read_departure'

      close(iunit)

      return
      end
