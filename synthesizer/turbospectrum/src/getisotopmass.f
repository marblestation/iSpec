cc      program test
cc
cc      implicit none
cc      real amass(92,0:250)
cc
cc      call getisotopmass(amass)
cc
cc      end
****************
      subroutine getisotopmass(amass)

      implicit none

      integer i,j,natom,debut,fin
      logical first
      character firstchar*1,oneline*72,cleaned*72
      parameter (natom=92)
      real amass(natom,0:250),frac(natom,0:250),massmix(natom)
      real masscheck

      do j=1,250
        do i=1,natom
          amass(i,j)=0.0
        enddo
      enddo

      open(72,file='DATA/atomicweights.dat')

* read header 
      read(72,10) firstchar
10    format(a)
      do while (firstchar.ne.'_')
        read(72,10) firstchar
      enddo
      backspace(72)

      i=0
* read data
      do while (i.le.natom) 
        read(72,10) oneline
        read(oneline,10) firstchar
        if (firstchar.eq.'_') then
* new element
          first=.true.
        else
* read data for new element; first line is different
          if (first) then
            read(oneline,20) i,j
20          format(i3,5x,i3)
            first=.false.
          else
            read(oneline,30) j
30          format(8x,i3)
          endif
          if (i.le.natom) then
            fin=index(oneline,'(')
            cleaned=oneline(12:fin-1)
            read(cleaned,*) amass(i,j)
            oneline=oneline(fin+1:72)
            debut=index(oneline,')')
            oneline=oneline(debut+1:72)
            if (oneline(1:1).eq.'#') oneline=oneline(2:72)
            fin=0
            fin=index(oneline,'(')
* First make sure there is data at end of line
            if (fin.ne.0) then
              cleaned=oneline(1:fin-1)
              read(cleaned,*) frac(i,j)
              debut=index(oneline,')')
              oneline=oneline(debut+1:72)
              if (oneline(1:1).eq.'#') oneline=oneline(2:72)
              fin=0
              fin=index(oneline,'(')
              if (fin.ne.0) then
                cleaned=oneline(1:fin-1)
                read(cleaned,*) massmix(i)
              endif
            endif
          endif
        endif
      enddo

c      print*,'TEST mass isotopes'
c      do i=1,92
c        masscheck=0.0
c        do j=1,250
c          if (amass(i,j).ne.0.0) then
c            print*,i,j,amass(i,j),frac(i,j)
c            masscheck=masscheck+amass(i,j)*frac(i,j)
c          endif
c        enddo
c        print*,i,massmix(i),masscheck
c        if (1.-massmix(i)/masscheck.gt.0.01) print*,'PROBLEM!'
c      enddo

      return
      end
