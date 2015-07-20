      program xy2moog
*
      character file*128,title*80
      doubleprecision x
      dimension x(10),y(10)
*
      print *,' input file to convert?'
      read(5,10) file
      open(10,file=file,status='old')
      print*,' output file (will be in format suitable for moog)'
      read(5,10) file
      open(20,file=file,status='new')
      print*,'title to give to the output spectrum?'
      read(5,10) title
      write(20,11) title
      rewind(10)
      read(10,*) x(1),y(1)
      read(10,*) x(2),y(2)
      icount=2
      dlam=x(2)-x(1)
1     read(10,*,end=99) x(2),y(2)
      icount=icount+1
      goto 1
99    continue
      write(20,12) x(1),x(2),dlam,dummy
      rewind(10)
      do 2 i=1,int(icount/10.)
        do 3 j=1,10
          read(10,*) x(j),y(j)
3       continue
        write(20,13) (1.-y(j),j=1,10)
2     continue
      jcount=0
      do 4 j=1,10
        read(10,*,end=999) x(j),y(j)
        jcount=jcount+1
4     continue
999   continue
      if (jcount.gt.0) write(20,13) (1.-y(j),j=1,jcount)
*
10    format (a)
11    format (a80)
12    format (4f10.2)
13    format (10f7.4)
*
      end

