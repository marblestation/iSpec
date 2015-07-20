*
************************************************************************
*
cc      SUBROUTINE CLOCK(message)
      SUBROUTINE CLOCK
*
*  TIME SINCE LAST CALL,ACCUMULATED EXECUTION TIME
*
*
* statement for portability of etime and dtime to intel Fortran
      use ifport

c      real dtime,etime
      real tt(2),ett(2),tempspart,tempsacc
      character*80 message
      logical info
      data info/.true./

      message='no message'
*
      if (info) then
        tempspart=dtime(tt)
        tempsacc=etime(ett)
*
        WRITE(7,50) tt(1),ett(1),tt(2),ett(2),tempspart,tempsacc
        write(7,'(a)') message
50      FORMAT('0','user_TIME ',F9.3,2x,'ACCMLTD ',F9.3,2x,
     &    'syst_TIME ',F9.3,2x,'ACCMLTD ',F9.3,2x,
     &    'totl_TIME ',f9.3,2x,'ACCMLTD',f9.3)
      endif

      return

      END
