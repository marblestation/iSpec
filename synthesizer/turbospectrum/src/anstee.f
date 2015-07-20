*     program testanstee
*     character recipe*1
*   1 continue
*     print *,'give iel ion wavelength elower & eionization'
*     read(*,*) iel,ion,xlb,elo,eion
*     call anstee(iel,ion,xlb,elo,eion,sigma,velexp,levlo,levup,recipe)
*     write(*,*) 'iel ion xlb elo eion sigma velexp levlo levup recipe'
*     write(*,*) iel,ion,xlb,elo,eion,sigma,velexp,levlo,levup,recipe
*     goto 1
*     end
      subroutine anstee(iel,ion,xlb,elo,eion,sigma,velexp,levlo,levup,
     &                                                          recipe)
* Get Hydrogen damping data for s-p, p-s, p-d, d-p, d-f and f-d transitions in neutral species
* and for a large number of ionized lines
* ref Anstee S.D., O'Mara B.J. 1995, MNRAS 276, 859-866         'A'
* Barklem P.S., O'Mara B.J. 1997, MNRAS 290, 102-106            'B'
* Barklem P.S., O'Mara B.J., Ross J.E. 1998, MNRAS 296, 1057    'C'
* Barklem P.S., O'Mara B.J. 1998, MNRAS 300, 863                'S'
* Barklem 1997, 1999, private comm (approx. for some lines)     'S'
* Barklem P.S., Piskunov N., O'Mara B.J. 2000, A&AS 142, 467    'S'
* Barklem P.S. 2003, private communication (IR lines)           'S'
* Barklem, Aspelund-Johansson 2005, Fe II, astro-ph/0502098     'S'
* Barklem 2008, private communication Cr II as for Fe II above  'S'
* 
* table1 and table2 gives sigma(s,p) and velexp(s,p)
* table3 and table4 gives sigma(p,d) and velexp(p,d)
* table5 and table6 gives sigma(d,f) and velexp(d,f)
*
* input:
* iel    = atomic number
* ion    = species (neutral=1)
* xlb    = wavelength in Aangstroem
* elo    = excitation energy of the lower level in eV
* eion   = ionisation energy in eV
* levlo  = state designation of the lower level, s p d f g etc
* levup  = state designation of the upper level, s p d f g etc
* output:
* sigma  = collision crossection in atomic units at 10^4 m/s
* velexp = velocity exponent "alpha"
* recipe = A if Anstee & O'Mara (s-p), B if Barklem & O'Mara (p-d) or 1998,
*          C if Barklem, O'Mara & Ross (d-f), S Barklem & O'Mara (ions),
*          Barklem, private comm. or Barklem & Aspelund-Johansson 2005.
*          Else U (Unsoeld)
*
* Updates 2008-02-27 /B. Edvardsson:
* Added 13165 Cr II lines from Paul + 226 long-wavelength Fe II lines
* Much improved identification of 'special' Fe II and Cr II lines by
*      use of also the excitation energy/8067.
* Homogenized different values of the H I ionization energy to 13.595 eV
*
      implicit none
      integer nspecial,nfe2,ncr2
      parameter (nspecial=5024)
      parameter (nfe2=24188)
      parameter (ncr2=13165)
      real table1(10:50,10:50),table2(10:50,10:50)
      real table3(10:50,10:50),table4(10:50,10:50)
      real table5(10:50,10:50),table6(10:50,10:50)
      real elo,eup,eion,xlb,nstars,nstarp,nstard,nstarf,sigma,velexp
      real xlbs(nspecial),sigmas(nspecial),alphas(nspecial)
      real xlbfe2(nfe2),sigmafe2(nfe2),alphafe2(nfe2)
      real xlbcr2(ncr2),sigmacr2(ncr2),alphacr2(ncr2)
      real elos(nspecial),elofe2(nfe2),elocr2(ncr2)
      integer i,j,iel,ion,ielfe2,ionfe2,ielcr2,ioncr2
      integer nsp,iels(nspecial),ions(nspecial),nfe,ncr
      character levlo*1,levup*1,string*132,check*1,recipe*1
      logical first
*
      data first / .true. /
*
      recipe='U'
      nstars=0.0
      nstarp=0.0
      nstard=0.0
      nstarf=0.0
      sigma=0.0
      velexp=0.0
*     eup=12398.426/xlb+elo (not consistent with eqwidth.f)
      eup=12399.8/xlb+elo
      if(eup.ge.eion) then
        print *,'subr. anstee: Upper excit. energy .ge. ioniz. energy;'
        print *,'el,ion,wavelength,chi,I=',iel,ion,xlb,eup,eion
        return
      endif
      if(ion.gt.2.or.ion.lt.1) then
        return
* no such data, use Unsoeld if anything
      endif
*
*
      if(first) then
*
* Read s-p and p-s data from Anstee.dat for sigma and velexp
*
        open(612,
     &  file='DATA/Anstee-1802-2014.dat',
     &           status='old',form='formatted')
*
* First read data from table 1, not commentary lines
* Horizontally: p-states 1.3 - 3.0, step 0.1  (18 values)
* Vertically:   s-states 1.0 - 3.0, step 0.1  (21 values)
*
    1   continue
        read(612,1000) string
 1000   format(a132)
        read(string,1010) check
 1010   format(a1)
        if (check.ne.'*') then
          read(string,*) (table1(j,10),j=13,30)
          do i=11,30
            read(612,*) (table1(j,i),j=13,30)
          enddo
        else
          goto 1
        endif
*
* Read data from table 2, not commentary lines
*
    2   continue
        read(612,1000) string
        read(string,1010) check
        if (check.ne.'*') then
          read(string,*) (table2(j,10),j=13,30)
          do i=11,30
            read(612,*) (table2(j,i),j=13,30)
          enddo
        else
          goto 2
        endif
*
* Then read data from table 3, not commentary lines
* Horizontally: d-states 2.3 - 4.0, step 0.1  (18 values)
* Vertically:   p-states 1.3 - 3.0, step 0.1  (18 values)
*
    3   continue
        read(612,1000) string
        read(string,1010) check
        if (check.ne.'*') then
          read(string,*) (table3(j,13),j=23,40)
          do i=14,30
            read(612,*) (table3(j,i),j=23,40)
          enddo
        else
          goto 3
        endif
*
* Read data from table 4, not commentary lines
*
    4   continue
        read(612,1000) string
        read(string,1010) check
        if (check.ne.'*') then
          read(string,*) (table4(j,13),j=23,40)
          do i=14,30
            read(612,*) (table4(j,i),j=23,40)
          enddo
        else
          goto 4
        endif
*
* Then read data from table 5, not commentary lines
* Horizontally: f-states 3.3 - 5.0, step 0.1  (18 values)
* Vertically:   d-states 2.3 - 4.0, step 0.1  (18 values)
*
    5   continue
        read(612,1000) string
        read(string,1010) check
        if (check.ne.'*') then
          read(string,*) (table5(j,23),j=33,50)
          do i=24,40
            read(612,*) (table5(j,i),j=33,50)
          enddo
        else
          goto 5
        endif
*
* Read data from table 6, not commentary lines
*
    6   continue
        read(612,1000) string
        read(string,1010) check
        if (check.ne.'*') then
          read(string,*) (table6(j,23),j=33,50)
          do i=24,40
            read(612,*) (table6(j,i),j=33,50)
          enddo
        else
          goto 6
        endif
*
* read data for special lines and store
*
        do i=1,14
          read(612,*)
        enddo
        do i=1,nspecial
          read(612,1020) iels(i),ions(i),xlbs(i),elos(i),
     &                  sigmas(i),alphas(i)
 1020     format(i2,3x,i2,f10.3,f11.3,59x,f8.2,f8.3)
          nsp=i
        enddo
        if(nsp.ne.nspecial) then
          print *,'nsp .ne. nspecial:',nsp,nspecial
          stop 'Why do nsp and nspecial differ?'
        endif
* Now read and store the large number of Fe II lines:
        do i=1,5
          read(612,*)
        enddo
        do i=1,nfe2
          read(612,1030) ielfe2,ionfe2,xlbfe2(i),elofe2(i),
     &                  sigmafe2(i),alphafe2(i)
 1030     format(i2,3x,i2,f10.3,f11.3,59x,f8.2,f8.3)
          nfe=i
        enddo
        if(nfe.ne.nfe2) then
          print *,'nfe .ne. nfe2:',nfe,nfe2
          stop 'Why do nfe and nfe2 differ?'
        endif
* Finally read and store the large number of Cr II lines:
        do i=1,6
          read(612,*)
        enddo
        do i=1,ncr2
          read(612,1030) ielcr2,ioncr2,xlbcr2(i),elocr2(i),
     &                  sigmacr2(i),alphacr2(i)
          ncr=i
        enddo
        if(ncr.ne.ncr2) then
          print *,'ncr .ne. ncr2:',ncr,ncr2
          stop 'Why do ncr and ncr2 differ?'
        endif
*
        close(612)
        first=.false.
      endif
* 
      if(ion.eq.1) then
*
        if ((levlo.eq.'s'.or.levlo.eq.'S').and.
     &        (levup.eq.'p'.or.levup.eq.'P')) then
          nstars=1./(sqrt((eion-elo)/13.595))
          nstarp=1./(sqrt((eion-eup)/13.595))
          if ((nstarp.ge.1.3.and.nstarp.le.3.0).and.
     &         (nstars.ge.1.0.and.nstars.le.3.0))then
            call interpolation(table1,table2,nstars,nstarp,
     &                                          sigma,velexp)
            recipe='A'
          endif
        else if ((levlo.eq.'p'.or.levlo.eq.'P').and.
     &             (levup.eq.'s'.or.levup.eq.'S')) then
          nstarp=1./(sqrt((eion-elo)/13.595))
          nstars=1./(sqrt((eion-eup)/13.595))
          if ((nstarp.ge.1.3.and.nstarp.le.3.0).and.
     &          (nstars.ge.1.0.and.nstars.le.3.0)) then
            call interpolation(table1,table2,nstars,nstarp,
     &                                         sigma,velexp)
            recipe='A'
          endif
*
        else if ((levlo.eq.'p'.or.levlo.eq.'P').and.
     &             (levup.eq.'d'.or.levup.eq.'D')) then
          nstarp=1./(sqrt((eion-elo)/13.595))
          nstard=1./(sqrt((eion-eup)/13.595))
          if ((nstarp.ge.1.3.and.nstarp.le.3.0).and.
     &          (nstard.ge.2.3.and.nstard.le.4.0)) then
            call interpolation(table3,table4,nstarp,nstard,
     &                                         sigma,velexp)
            recipe='B'
          endif
        else if ((levlo.eq.'d'.or.levlo.eq.'D').and.
     &             (levup.eq.'p'.or.levup.eq.'P')) then
          nstard=1./(sqrt((eion-elo)/13.595))
          nstarp=1./(sqrt((eion-eup)/13.595))
          if ((nstarp.ge.1.3.and.nstarp.le.3.0).and.
     &          (nstard.ge.2.3.and.nstard.le.4.0)) then
            call interpolation(table3,table4,nstarp,nstard,
     &                                         sigma,velexp)
            recipe='B'
          endif
*
        else if ((levlo.eq.'d'.or.levlo.eq.'D').and.
     &             (levup.eq.'f'.or.levup.eq.'F')) then
          nstard=1./(sqrt((eion-elo)/13.595))
          nstarf=1./(sqrt((eion-eup)/13.595))
          if ((nstard.ge.2.3.and.nstard.le.4.0).and.
     &          (nstarf.ge.3.3.and.nstarf.le.5.0)) then
            call interpolation(table5,table6,nstard,nstarf,
     &                                         sigma,velexp)
            recipe='C'
          endif
        else if ((levlo.eq.'f'.or.levlo.eq.'F').and.
     &             (levup.eq.'d'.or.levup.eq.'D')) then
          nstarf=1./(sqrt((eion-elo)/13.595))
          nstard=1./(sqrt((eion-eup)/13.595))
          if ((nstard.ge.2.3.and.nstard.le.4.0).and.
     &          (nstarf.ge.3.3.and.nstarf.le.5.0)) then
            call interpolation(table5,table6,nstard,nstarf,
     &                                         sigma,velexp)
            recipe='C'
          endif
        endif
*
      endif
*
* Check also whether the line has been individually treated:
* Then this data is preferred.
*
      if(iel.eq.24 .and. ion.eq.2) then
* Cr II
        do i=1,ncr
          if(iel.eq.ielcr2 .and. ion.eq.ioncr2) then
            if(abs(xlb-xlbcr2(i)).le.0.05) then
* wavelengths within 0.05 Aangstroem
              if(abs(elo-elocr2(i)/8067.).lt.0.02) then
* energies within 0.02 eV
                recipe='S'
                velexp=alphacr2(i)
                sigma= sigmacr2(i)
                print *,'Special Cr 2 Barklem priv.',iel,ion,xlb
                goto 8
              endif
            endif
          else
            stop 'This was supposedly a Cr II line!'
          endif
        enddo
    8   continue
      else if(iel.eq.26 .and. ion.eq.2) then
* Fe II
        do i=1,nfe
          if(iel.eq.ielfe2 .and. ion.eq.ionfe2) then
            if(abs(xlb-xlbfe2(i)).le.0.05) then
              if(abs(elo-elofe2(i)/8067.).lt.0.02) then
                recipe='S'
                velexp=alphafe2(i)
                sigma= sigmafe2(i)
                print *,'Special Fe 2 Barklem & A-J',iel,ion,xlb
                goto 9
              endif
            endif
          else
            stop 'This was supposedly an Fe II line!'
          endif
        enddo
    9   continue
      else
* all other 'special' lines
        do i=1,nsp
          if(iel.eq.iels(i) .and. ion.eq.ions(i)) then
            if(abs(xlb-xlbs(i)).le.0.05) then
              if(abs(elo-elos(i)/8067.).lt.0.02) then
                recipe='S'
                velexp=alphas(i)
                sigma= sigmas(i)
                print *,'Special line P. Barklem',iel,ion,xlb
                goto 10
              endif
            endif
          endif
        enddo
   10   continue
      endif
*
      if(ion.eq.1 .and. recipe.eq.'U') then
* the last resort
        if((levlo.eq.'d'.or.levlo.eq.'D').and.
     &         (levup.eq.'p'.or.levup.eq.'P')) then
          nstard=1./(sqrt((eion-elo)/13.595))
          if(nstard.lt.2.0) then
* The upper p-state will dominate, use that one alone
* since there is (as yet) no data for the lower level:
* ref: Paul Barklem, email 1997 
            nstarp=1./(sqrt((eion-eup)/13.595))
            if(nstarp.ge.2.6.and.nstarp.le.2.7) then
              sigma=762.+960.0*(nstarp-2.6)
              velexp=0.270-0.04*(nstarp-2.6)
              recipe='S'
              print *,'Special line, only upper p level',iel,ion,xlb
            else if(nstarp.gt.2.7.and.nstarp.le.2.8) then
              sigma=858.+1110.*(nstarp-2.7)
              velexp=0.266+0.02*(nstarp-2.7)
              recipe='S'
              print *,'Special line, only upper p level',iel,ion,xlb
            else if(nstarp.gt.2.8.and.nstarp.le.2.9) then
              sigma=969.+610.0*(nstarp-2.8)
              velexp=0.268-0.10*(nstarp-2.8)
              recipe='S'
              print *,'Special line, only upper p level',iel,ion,xlb
            else if(nstarp.gt.2.9.and.nstarp.le.3.0) then
              sigma=1030.+440.*(nstarp-2.9)
              velexp=0.258-0.04*(nstarp-2.9)
              recipe='S'
              print *,'Special line, only upper p level',iel,ion,xlb
            endif
          endif
        endif
      endif
*
      if(recipe.ne.'A' .and. recipe.ne.'B' .and. recipe.ne.'C'
     &    .and. recipe.ne.'S') then
*
* OK, we have to stick to Unsolds recipe
*
        recipe='U'
      endif
*
      return
      end
*
************************************************************************
*
      subroutine interpolation(tablea,tableb,nstar1,nstar2,sigma,velexp)
* interpolate in broadening data tables to get cross section and
* velocity exponent for s-p, p-s, p-d, d-p, d-f or f-d transitions
*
      implicit none
      real tablea(10:50,10:50),tableb(10:50,10:50),sigma,velexp
      real nstar1,nstar2,fract1,fract2,a,b,c,d,x,y
      integer i1,i2
*
      i2=int(10.*nstar2)
      i1=int(10.*nstar1)
      fract2=(10.*nstar2-float(i2))
      fract1=(10.*nstar1-float(i1))
* interpolate in cross section table, first in the 2 direction
      a=tablea(i2,i1)
      b=tablea(i2+1,i1)
      c=tablea(i2,i1+1)
      d=tablea(i2+1,i1+1)
      x=a*(1.-fract2)+b*fract2
      y=c*(1.-fract2)+d*fract2
* interpolate in the 1 direction
      sigma=x*(1.-fract1)+y*fract1
* interpolate in velocity exponent table, first in the 2 direction
      a=tableb(i2,i1)
      b=tableb(i2+1,i1)
      c=tableb(i2,i1+1)
      d=tableb(i2+1,i1+1)
      x=a*(1.-fract2)+b*fract2
      y=c*(1.-fract2)+d*fract2
* interpolate in the 1 direction
      velexp=x*(1.-fract1)+y*fract1
*
      return
      end
*
************************************************************************
*
