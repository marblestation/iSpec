!
      subroutine read_modelatom(iunit,modelatomfile,maxn,modnlevel,
     &                          modenergy,modg,modion,modid,
     &                          species)
!
! read NLTE model atom
! created by BPlez 15/04-2020
!
! unit number to open, file name, max number of levels that can be accomodated,
! number of level read, energy, statistical weight, ionisation stage, and 
! identification of the levels
!
      implicit none
      integer maxn
      character modelatomfile*256,oneline*256,species*20
      character modid(maxn)*40,bla*20
      integer iunit,modnlevel,modion(maxn),i
      real modenergy(maxn),abundance,mass,modg(maxn)
      logical header,header2,header3

      open(iunit,file=modelatomfile,status='old')
      header=.true.
      do while (header)
        read(iunit,10,end=99) oneline
        if (oneline(1:1).ne.'*') then
          read(oneline,10) species
          if (species(1:1).eq.' ') then
! move leading blank, if any, to end, in order to conform 
! to aname standard (see makeabund.f)
            write(bla,20) species(2:20)
            species=bla
          endif
          header2=.true.
          do while (header2)
            read(iunit,10,end=99) oneline
            if (oneline(1:1).ne.'*') then
              read(oneline,*) abundance,mass
              header2=.false.
            endif
          enddo
          header3=.true.
          do while (header3)
            read(iunit,10,end=99) oneline
            if (oneline(1:1).ne.'*') then
              read(oneline,*) modnlevel
              header3=.false.
            endif
          enddo
          header=.false.
        endif
      enddo
      i=0
      do while (i.lt.modnlevel)
        read(iunit,10,end=99) oneline
        if (oneline(1:1).ne.'*') then
          i=i+1
          read(oneline,*) modenergy(i),modg(i),modid(i),modion(i)
        endif
      enddo

      print*,'read_modelatom ',species,modnlevel,' levels'
      do i=1,modnlevel
        print*,i,modion(i),modid(i),modenergy(i),modg(i)
      enddo
10    format(a)
20    format(a19,' ')

      close(iunit)

      return

99    print*,'unexpected end of file in read_modelatom'
      stop 'stop in read_modelatom'

      end
