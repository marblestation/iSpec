      REAL FUNCTION NHOP(WAVE,T,PATH,PATHLEN,NHtbl_unit)
C
C Cross-sections of NH photoionization in A^2/molecule in LTE
C This routine is based on data provided by Phillip Stancil
C BPz 30/07-2019
C
C Input parameters:
C  WAVE    - (R*8) wavelength in Angstroems
C  T       - (R*4) temperature in K
C  PATH    - (CHAR) path to data file
C  PATHLEN - (I*4) length of the PATH string (not PATH variable!)
C  NHtbl_unit - (I*4) unit number on which the file file be open at the first call.
C               It just has to have no conflict with other open units. No need to
C               keep it for the next calls or otherwise save it.
C
      implicit none
      real t,t_tbl(15),gauss_fwhm
      real*8 wave,wl0,wlstep,gcross(4701,15),wl(4701),wlstepold
      real*8 factor_wl,factor_temp,f1,f2
      integer n_WL,n_Temp,NHtbl_unit,pathlen,i
      integer i_wl,i_temp
      character*(*) path
      character*(512) filename
      character*(2048) head
      logical first
C
      save first,wl0,wlstep,n_wl,n_temp,t_tbl,gcross
      data first/.true./
C
      nhop=0.
C
C Read data file
C
      if (first) then
        if (NHtbl_unit < 1.or.NHtbl_unit == 5) then
          write(*,*) 'NHOP: illegal unit number for cross-section data'
          return
        endif
        filename=path(1:pathlen)//'NHcont_Stancil_2019.dat'//' '
        write(*,'(A)') trim(filename)
        open(unit=NHtbl_unit,file=trim(filename),status='OLD',
     *       action='READ')

        read(NHtbl_unit,*)
        read(NHtbl_unit,*)
        read(NHtbl_unit,*) n_WL
        if (n_wl > 4701) stop 'NHOP.f: n_wl too large in file!'
        read(NHtbl_unit,*) n_Temp
        if (n_temp > 15) stop 'NHOP.f: n_temp too large in file!'
        read(NHtbl_unit,*) (T_TBL(i_temp),i_temp=1,n_Temp)
        do i_wl=1,n_wl
          read(NHtbl_unit,*) wl(i_wl),(gcross(i_wl,i_temp),
     &         i_temp=1,n_temp)
        enddo
        wlstepold=wl(2)-wl(1)
        do i_wl=3,n_wl
          wlstep=wl(i_wl)-wl(i_wl-1)
          if ((wlstep-wlstepold)/wlstep > 1.e-4) then
            stop 'problem with wavelength step in NHOP.f!'
          endif
          wlstepold=wlstep
        enddo
        wl0=wl(1)
        close(NHtbl_unit)
        first=.false.
      endif
C
C No extrapolation
C
      if (wave.lt.wl0.or.wave.gt.wl0+wlstep*(n_wl-1)) return
      if (t.lt.t_tbl(1).or.t.gt.t_tbl(n_temp)) return

      i_wl=(wave-wl0)/wlstep
      factor_wl=(wave-wl0-i_wl*wlstep)/wlstep

      do i=1,n_temp-1
        i_temp=i
        if(t_tbl(i_temp+1).gt.t) go to 1
      enddo

  1   factor_temp=(t-t_tbl(i_temp))/(t_tbl(i_temp+1)-t_tbl(i_temp))

      f1=(gcross(i_wl  ,i_temp+1)-gcross(i_wl  ,i_temp))*factor_temp
     +   +gcross(i_wl  ,i_temp)
      f2=(gcross(i_wl+1,i_temp+1)-gcross(i_wl+1,i_temp))*factor_temp
     +   +gcross(i_wl+1,i_temp)
      NHop=(f2-f1)*factor_wl+f1
C
      RETURN
      END
