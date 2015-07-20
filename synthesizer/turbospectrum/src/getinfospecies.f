*
      subroutine getinfospecies(species,iel,natom,atom,isotope)
*
* BEWARE: species = character*26   ! 26=natom*5+1
* the largest molecules are 5-atomic (dim of atom()).
* species names are e.g.:   
*     '608.012016'    = 12C16O
* or '0608.012016'
*    '0822.000048'    = 48TiO (ALL ISOTOPES OF O !!)
*    '0822.000000'    = TiO   (ALL ISOTOPES OF Ti AND O!)
*     '101.001001'    = 1H1H
*       '3.006'       = 6Li
*   '10607.001012014' = 1H12C14N
*
* output: iel  = 608, 822, 822, 101, 3, 10607 in above cases
*         atom(1:natom)    = 6,8  1,1  3  1,6,7 in above cases
*         isotope(1:natom) = 12,16  1,1  6  1,12,14 in above cases
*
* 12/10-1995 BPz

      implicit none

      integer            iel,isotope(5),atom(5),natom,i,pos,shift
      character*26       species,spec

      do i=1,5
        atom(i)=0
        isotope(i)=0
      enddo
* clean up
      shift=0
      do while(species(1:1).eq.' ')
        spec=species
        write(species,'(a25,'' '')') spec(2:26)
        shift=shift+1
      enddo
      do while(species(1:1).eq.'0'.and.species(2:2).eq.'0')
        spec=species
        write(species,'(a25,'' '')') spec(2:26)
        shift=shift+1
      enddo
        
      pos=index(species,'.')
      if (pos/2.-pos/2.eq.0.) then
        if (2*pos+4+shift.gt.26) then
          stop 'getinfospecies: species string truncated'
        endif
        spec=species
        write(species,10) spec(1:25)
10      format('0',25a)
      else
        if (2*pos+3+shift.gt.26) then
          stop 'getinfospecies: species string truncated'
        endif
      endif
      pos=index(species,'.')
      natom=pos/2
      if (natom.gt.5) stop 'getinfospecies: natom >5'
      read(species,'(5i2)') (atom(i),i=1,natom)
      write(spec,'(a)') species((natom+1)*2:26)
      read(spec,'(5i3)') (isotope(i),i=1,natom)
      iel=0
      do i=1,natom
        iel=iel+100**(natom-i)*atom(i)
      enddo
      
      return
      end
