      subroutine input
*
* Ecrire un mode d'emploi.
*
* reads input from file 5
*
      parameter(maxfil=100)
      character*2 atominclude(100)
      character*30 keyword,answer
      character*256 charvalue
      character*128 filmet,filmol,filwavel
      character*256 linefil,detout,inatom,inmod,inabun,inspec,
     &              outfil,mongofil,filterfil,continopac,inpmod
      logical tsuji,spherical,limbdark,abfind,multidump,xifix,mrxf,
     &        hydrovelo,pureLTE
      integer iint,k
      real    isoch(1000),isochfact(1000),xic,xmyc,scattfrac
      doubleprecision xl1,xl2,del,xlmarg,xlboff
      integer nsegments_lambda_min, nsegments_lambda_max !SBC
      doubleprecision segments_lambda_min(400), segments_lambda_max(400)
      common/inputdata/mmaxfil,tsuji,filmet,filmol,noffil,
     &                 linefil(maxfil),spherical,mihal,taum,ncore,
     &                 diflog,detout,inatom,
     &                 inmod,inabun,inspec,outfil,
     &                 mongofil,abch(100),limbdark,filterfil,
     &                 overall,abfind,multidump,isoch,isochfact,
     &                 helium,alpha,rabund,sabund,xifix,xic,mrxf,
     &                 inpmod,continopac,filwavel,hydrovelo,
     &                 xl1,xl2,del,xlmarg,xlboff,iint,xmyc,scattfrac,
     &                 pureLTE
      common/species/atominclude
      data atominclude 
     &   /'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &    'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     &    'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     &    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &    'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     &    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &    'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     &    'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &    'Pa','U ','  ','  ','  ','  ','  ','  ','  ','  ' /
*
      scattfrac=0.0
      pureLTE=.false.
      do k=1,100
        abch(k)=-99.9
      enddo
      isoch=-1.0
      isochfact=-1.0
      abfind=.false.
*
      if (maxfil.ne.mmaxfil) stop 'parameter maxfil wrong'
      iread=5
      iwrite=6
* default data files and setup:

      continopac='DATA/jonabs_vac_v07.3.dat'
      filmol=
     &  'DATA/LISTE_molecules_all_v12.1.dat'
ccc     & '/b1/plez/BIGGRID/DATA/LISTE_molecules_all.dat'
      filmet=' '
cc     &'/b1/plez/BIGGRID/DATA/tsuji.atoms_AndersGrev'
      detout='dummy-output.dat'
ccc      inatom='/b1/plez/BIGGRID/DATA/atomdata'
ccc      inatom='DATA/atomdata-v12.1'
        inatom=' '
      filwavel='DATA/wavelengths.UVIS'
      outfil='syntspec/synout'
      hydrovelo=.false.
      limbdark=.false.
      tsuji=.true.
      spherical=.false.
      multidump=.false.
      mrxf=.true.
      xifix=.true.
      xic=0.0
      overall=0.0
      helium=0.0
      alpha =0.0
      rabund=0.0
      sabund=0.0
      nchabund=0
      iint=0
      xmyc=1.0
      xlboff=0.d0
      xlmarg=2.d0
      xl1=-10.d0
      xl2=-10.d0
      del=-10.d0
      nsegments_lambda_min = 0 ! SBC
      nsegments_lambda_max = 0 ! SBC
*
1     read(iread,*,end=99) keyword,charvalue

      print*, keyword(1:lenstr(keyword)),charvalue(1:lenstr(charvalue))

      if (keyword(1:23).eq.'INCLUDE FOLLOWING ATOMS') then
        read(charvalue,*) nlines
        if (nlines.gt.0.and.nlines.le.10) then
* if atoms to be included given in input, we use them. Otherwise we 
* use the default list provided in the DATA block of this routine
          do jj=1,100
            atominclude(jj)='  '
          enddo
          do i=1,nlines
            read(5,*) (atominclude((i-1)*10+jj),jj=1,10)
          enddo
        else if (nlines.gt.10) then
          print*,'nlines = ',nlines,' too large. Limit = 100'
          stop 
        endif
      else if (keyword(1:11).eq.'METALLICITY') then
        read(charvalue,*) overall
      else if (keyword(1:8).eq.'ALPHA/Fe') then
        read(charvalue,*) alpha
      else if (keyword(1:6).eq.'HELIUM') then
        read(charvalue,*) helium
      else if (keyword(1:9).eq.'R-PROCESS') then
        read(charvalue,*) rabund
      else if (keyword(1:9).eq.'S-PROCESS') then
        read(charvalue,*) sabund
      else if (keyword(1:21).eq.'INDIVIDUAL ABUNDANCES') then
        read(charvalue,*) nchabund
        do i=1,nchabund
          read(5,*) iii,fifi
          abch(iii)=fifi
          print*,iii,fifi
        enddo
      else if (keyword(1:5).eq.'TSUJI') then
        read(charvalue,*) tsuji
      else if (keyword(1:9).eq.'SPHERICAL') then
        read(charvalue,*) spherical
        read(iread,*) mihal
        read(iread,*) taum
        read(iread,*) ncore
        read(iread,*) diflog
      else if (keyword(1:6).eq.'NFILES') then
        read(charvalue,*) noffil
        do i=1,noffil
          read(iread,'(a)') linefil(i)
        enddo
      else if (keyword(1:7).eq.'LOGFILE') then
        read(charvalue,10) detout
      else if (keyword(1:8).eq.'ATOMDATA') then
        read(charvalue,10) inatom
      else if (keyword(1:9).eq.'MODELOPAC') then
        read(charvalue,10) inmod
      else if (keyword(1:9).eq.'PARAMETER') then
        read(charvalue,10) inabun
      else if (keyword(1:8).eq.'SPECDATA') then
        read(charvalue,10) inspec
      else if (keyword(1:10).eq.'RESULTFILE') then
        read(charvalue,10) outfil
      else if (keyword(1:9).eq.'MOLECULES') then
        read(charvalue,10) filmol
      else if (keyword(1:8).eq.'LIMBDARK') then
        read(charvalue,*) limbdark
c We now save limb-darkening for each wavelength in the interval.
ccc        read(iread,*) filterfil
      else if (keyword(1:9).eq.'MULTIDUMP') then
        read(charvalue,*) multidump
      else if (keyword(1:8).eq.'ISOTOPES') then
        read(charvalue,*) niso
        if (niso.gt.1000) then
          stop 'Error! too many isotopic changes requested'
        endif
        do i=1,niso
          read(iread,*) isoch(i),isochfact(i)
        enddo
      else if (keyword(1:10).eq.'MARCS-FILE') then
        read(charvalue,*) mrxf
      else if (keyword(1:5).eq.'XIFIX') then
        read(charvalue,*) xifix
        read(iread,*) xic
      else if (keyword(1:20).eq.'CONTINUOUS-OPACITIES') then
        read(charvalue,10) continopac
      else if (keyword(1:10).eq.'MODELINPUT') then
        read(charvalue,10) inpmod
      else if (keyword(1:6).eq.'ABFIND') then
        read(charvalue,*) abfind
      else if (keyword(1:13).eq.'C_WAVELENGTHS') then
        read(charvalue,10) filwavel
      else if (keyword(1:14).eq.'HYDRODYN_DEPTH') then
        read(charvalue,*) hydrovelo
      else if (keyword(1:10).eq.'LAMBDA_MIN') then
        if (xl1.lt.0) then
            read(charvalue,*) xl1
        else
            nsegments_lambda_min = nsegments_lambda_min + 1 !SBC
            read(charvalue,*) segments_lambda_min(nsegments_lambda_min)
        endif
      else if (keyword(1:10).eq.'LAMBDA_MAX') then
        if (xl2.lt.0) then
            read(charvalue,*) xl2
        else
            nsegments_lambda_max = nsegments_lambda_max + 1 !SBC
            read(charvalue,*) segments_lambda_max(nsegments_lambda_max)
        endif
      else if (keyword(1:11).eq.'LAMBDA_STEP') then
        read(charvalue,*) del
      else if (keyword(1:14).eq.'INTENSITY/FLUX') then
        read(charvalue,*) answer
        if (answer(1:1).eq.'I') then
          iint=1
        else if (answer(1:1).eq.'F') then
          iint=0
        else
          print*,'INPUT: bad answer: ',answer,' for INTENSITY/FLUX'
          print*,'INPUT: answer: should be I or F'
          stop
        endif
      else if (keyword(1:10).eq.'COS(THETA)') then
        read(charvalue,*) xmyc
      else if (keyword(1:9).eq.'SCATTFRAC') then
! fraction of line opacity to include in scattering
        read(charvalue,*) scattfrac
      else if (keyword(1:8).eq.'PURE-LTE') then
! pureLTE=.true. means take all opacity into kappa. No scattering.
! This ensures that S=B
        read(charvalue,*) pureLTE
        if (pureLTE) then
          scattfrac=0.0
        endif
      else
        print*,'input.f; bad keyword: ',keyword
        stop
      endif

      goto 1
99    continue

      if (xl1.lt.0.d0.or.xl2.lt.0.d0.or.del.lt.0.d0) then
        print*,'ERROR! Provide lambda min max and deltalambda!'
        print*,'xl1 xl2 del =',xl1,xl2,del
        stop
      endif

      ! SBC:
        if (nsegments_lambda_min.ne.nsegments_lambda_max) then
          print*,'ERROR! Provide provide the same number of'
          print*,' lambda min max for segments'
          stop
        endif

cc* store atomic species to use in molecular equilibrium into file 26.
cc      open(26,status='scratch')
cc      do jj=1,100
cc        if (atominclude(jj).ne.'  ') then
cc          write(26,111) atominclude(jj)
cc111       format('''',a2,'''')
cc        endif
cc      enddo
      print*,
     & 'Following atomic species included in molecular equilibrium:'
      do i=1,10
        print*, (atominclude((i-1)*10+jj),' ',jj=1,10)
      enddo


10    format (a)

      return
 
      end
