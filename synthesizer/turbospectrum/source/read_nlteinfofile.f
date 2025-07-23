!
      subroutine read_nlteinfofile(iunit,nlteinfofile,iel,
     &       modelatomfile,departurefile,departbin,nlte_species)
!
! read NLTE information for each species
! created by BPlez 23/02-2021
!
! unit number to open, file name, 
!
      implicit none
      integer iel,iunit,at_number
      logical departbin,nlte_species
      character departurefile*256,nlteinfofile*256,bin_flag*16,
     &          oneline*512,pathmodel*256,pathdepart*256,
     &          element*2,nlte_flag*4,dep_file*256,at_file*256,
     &          modelatomfile*256

      open(iunit, file=nlteinfofile,status='old')
      read(iunit,10) oneline
      backspace(iunit)
      do while (oneline(1:27).ne.'# path for model atom files') 
        read(iunit,10,end=98) oneline
      enddo
      read(iunit,10,end=98) pathmodel
      do while (oneline(1:26).ne.'# path for departure files') 
        read(iunit,10) oneline
      enddo
      read(iunit,10) pathdepart
      do while (oneline(1:1).eq.'#')
        read(iunit,10,end=99) oneline
      enddo
      backspace(iunit)
      do while (.true.)
        read(iunit,*,end=99) at_number,element,nlte_flag,at_file,
     &                dep_file,bin_flag
        if (at_number.eq.iel) then
! found the species
          if (nlte_flag(1:4).eq.'nlte') then
            nlte_species=.true.
          else if (nlte_flag(1:3).eq.'lte') then
            nlte_species=.false.
          else
            stop 'unrecognised keyword nlte_species 
     & in read_nlteinfofile'
          endif
          if (bin_flag(1:3).eq.'bin') then
            departbin=.true.
          else if (bin_flag(1:3).eq.'asc') then
            departbin=.false.
          else
            stop 'unrecognised keyword departbin
     &  in read_nlteinfofile'
          endif
          modelatomfile=''
          if (at_file.ne.'') then
            modelatomfile=
     &       trim(adjustl(pathmodel))//trim(adjustl(at_file))
          endif
          departurefile=''
          if (dep_file.ne.'') then
            departurefile=
     &       trim(adjustl(pathdepart))//trim(adjustl(dep_file))
          endif
          close(iunit)
          return
        endif
      enddo
99    continue
!
! this species was not found in file. LTE by default
      nlte_species=.false.
      departbin=.false.
      departurefile=''
      modelatomfile=''

      close(iunit)
      return

98    print*,'The NLTE info file seems totally empty!'
      close(iunit)
      stop 'stop in read_nlteinfofile.f'

10    format(a)
      end
