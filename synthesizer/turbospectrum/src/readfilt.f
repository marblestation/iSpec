      subroutine readfilt(iunit)
*  BPz 10/02/93
      logical chk,limbdark
      character*80 filttitle
      character*256 filterfil
      common/filter/limbdark,ifilt,filtlam(1000),filttrans(1000),
     &              filterfil,filttitle
      data chk/.true./

      rewind(iunit)
      read(iunit,*) filttitle
      i=1
      do while (.true.)
        read(iunit,*,end=9) filtlam(i),filttrans(i)
        i=i+1
      enddo
9     continue
      ifilt=i-1
      if (ifilt.eq.1) stop 'readfilt: problem  with the filter'
      if (ifilt.gt.1000) stop 'readfilt: dimension too small'
      if (chk) then
*test:
        print*,' Filter transmission: '
        print*,filttitle
        do i=1,ifilt
          print*,i,filtlam(i),filttrans(i)
        enddo
        print*,' End of filter'
*endtest
      endif
      return
      end
