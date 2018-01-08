
      subroutine tablepop (option)
c******************************************************************************
c     this routine opens the table file containing information for a stellar
c     population, and reads the data in that file; the information is a
c     bit different for "abpop" and "synpop"
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Multimod.com'
      include 'Pstuff.com'
      integer option
      character*80 line

  
c*****open the model table input file and the summary table output file
      nftable = 18
      lscreen = 4
      array = 'MODEL TABLE INPUT FILE'
      nchars = 22
      call infile ('input ',nftable,'formatted  ',0,nchars,
     .             ftable,lscreen)      
      nf7out = 24
      lscreen = 6
      array = 'MODEL TABLE OUTPUT FILE'
      nchars = 18
      call infile ('output ',nf7out,'formatted  ',0,nchars,
     .             f7out,lscreen)


c*****read the table input for integrated light EW matching
      if (option .eq. 1) then
         read (nftable,1001) line
         if (line(1:5) .ne. 'abpop' ) then
            write(*,*) 'OOPS!  WRONG TABLE FOR ABPOP; I QUIT!'
            stop
         endif
         do i=1,1000
            read (nftable,1001) line
            if     (line(1:5) .eq. 'modpr') then
               call blankstring (modpre)
               modpre(1:70) = line(11:80)
            elseif (line(1:5) .eq. 'synpr') then
               call blankstring (synpre)
               synpre(1:70) = line(11:80)
            elseif (line(1:5) .eq. 'title') then
               call blankstring (abitle)
               popitle(1:74) = line(7:80)
               write (nf7out,1002) popitle(1:73)
            elseif (line(1:5) .eq. 'model') then
               do mmod=1,100
                  if (mmod .eq. 100) then
                     write(*,*) 'MORE THAN 99 MODELS; I QUIT!'
                     stop
                  endif
                  read (nftable,*,end=10) j, radius(mmod), 
     .                                    relcount(mmod)
               enddo
            endif
         enddo


c*****read the table input for integrated light spectrum syntheses
      else
         nabs = 0
         nisos = 0
         read (nftable,1001) line(1:80)
         if (line(1:6) .ne. 'synpop') then
            write (*,*) 'OOPS!  WRONG TABLE FOR SYNPOP; I QUIT!'
            stop
         endif
         do k=1,1000
            call blankstring (line)
            read (nftable,1001) line(1:80)
            if     (line(1:5) .eq. 'modpr') then
               call blankstring (modpre)
               modpre(1:70) = line(11:80)
            elseif (line(1:5) .eq. 'synpr') then
               call blankstring (synpre)
               synpre(1:70) = line(11:80)
            elseif (line(1:5) .eq. 'title') then
               call blankstring (popitle)
               popitle(1:74) = line(7:80)
               write (nf7out,1003) popitle(1:74)
            elseif (line(1:5) .eq. 'abund') then
               read (line(12:80),*) nabs
               if (nabs .gt. 0) then
                  read (nftable,1001) line(1:80)
                  read (line(1:80),*) (elspecial(i),i=1,nabs)
                  write (nf7out,1009) (nint(elspecial(i)),i=1,nabs)
               else
                  write (nf7out,1010)
               endif
            elseif (line(1:5) .eq. 'isoto') then
               read (line(10:80),*) nisos
               if (nisos .gt. 0) then
                  read (nftable,1001) line(1:80)
                  read (line(1:80),*) (isospecial(i),i=1,nisos)
                  write (nf7out,1011) (isospecial(i),i=1,nisos)
               else
                  write (nf7out,1012)
               endif
            elseif (line(1:5) .eq. 'model') then
               do mmod=1,100
                  if (mmod .eq. 100) then
                     write (*,*) 'MORE THAN 99 MODELS; I QUIT!'
                     stop
                  endif
                  if     (nabs.le.0 .and. nisos.le.0) then
                     read (nftable,*,end=10) j, radius(mmod),
     .                                       relcount(mmod)
                     write (nf7out,1006) j, radius(mmod), 
     .                                   relcount(mmod)     
                  elseif (nabs.gt.0 .and. nisos.le.0) then       
                     read (nftable,*,end=10) j, radius(mmod),    
     .                                       relcount(mmod),     
     .                                      (abspecial(mmod,i),i=1,nabs)
                     write (nf7out,1006) j, radius(mmod), 
     .                                   relcount(mmod),       
     .                                   (abspecial(mmod,i),i=1,nabs)  
                  elseif (nabs.le.0 .and. nisos.gt.0) then       
                     read (nftable,*,end=10) j, radius(mmod),    
     .                                       relcount(mmod),     
     .                                   (fracspecial(mmod,i),i=1,nisos)
                     write (nf7out,1006) j, radius(mmod), relcount(mmod)
                     write (nf7out,1015) (fracspecial(mmod,i),i=1,nisos)
                  else                                           
                     read (nftable,*,end=10) j, radius(mmod),    
     .                                       relcount(mmod),     
     .                                   (abspecial(mmod,i),i=1,nabs),
     .                                   (fracspecial(mmod,i),i=1,nisos)
                     write (nf7out,1006) j, radius(mmod), 
     .                                   relcount(mmod),       
     .                                   (abspecial(mmod,i),i=1,nabs)
                     write (nf7out,1015) (fracspecial(mmod,i),i=1,nisos)
                  endif                                          
               enddo                                             
            endif                                                
         enddo                                                   
      endif


c*****close the model input file, exit normally
 10   close (unit=nftable)
      modtot = mmod-1
      return


c*****format statements
1001  format (a80)
1002  format ('EQUIVALENT WIDTH ANALYSIS FOR INTEGRATED-LIGHT SPECTRA'//
     .        a80)
1003  format ('POPULATION SYNTHESIS FOR INTEGRATED-LIGHT SPECTRA'/a80)
1006  format (i3, 1x, f8.0, 2f8.2, 10f6.2)
1009  format ('SPECIAL ELEMENTS OR NAMES OF ISOTOPES'/(10i8))
1010  format ('NO DECLARED SPECIAL ELEMENTS')
1011  format ('SPECIAL ISOTOPE NAMES'/((5f10.5)))
1012  format ('NO DECLARED SPECIAL ISOTOPES')
1015  format (6f13.5)


      end








