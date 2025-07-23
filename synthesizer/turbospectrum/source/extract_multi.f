      program extract_multi

* read multi output from bsyn
* BPz 15/07-04
*
      implicit none
      character*128 filein
      logical star
      integer ieliel,ntau,maxlam,k,j
      real absave(83),absos(150,10000),abso(150,10000),
     &       absocont(150,10000),tau(150),xi(150),xls
      doubleprecision xl1,xl2,del,lambda(10000)

      star=.true.

      print*,' Input file with multi output from bsyn?'
      read(5,10) filein
10    format (A)

      open(10,form='unformatted", file=filein)

      do while (star)
        read(10) bla
        if (bla(1:1).ne.'*') then
          star=.false.
          backspace(10)
        endif
      enddo

      read(10) (ieliel,absave(ieliel),ieliel=1,60),
     &          (ieliel,absave(ieliel),ieliel=62,83)
      read(10) ntau,maxlam,xl1,xl2,del
      read(10) xls
      read(10) (tau(k),k=1,ntau),(xi(k),k=1,ntau)
      read(10) (absocont(k,j),k=1,ntau),j=1,maxlam)
      read(10) (abso(k,j),k=1,ntau),j=1,maxlam)
      read(10) ((absos(k,j)*ross(k),k=1,ntau),j=1,maxlam)

      do j=1,maxlam
        lambda(j)=xl1+(j-1)*del
      enddo

      stop
      end

