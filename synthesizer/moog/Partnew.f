
      real*8 function partnew (atom,k,level)
c******************************************************************************
c     This routine computes partition functions for those species that have
c     been updated from those that are in ATLAS9.  The polynomial 
c     representations are in the form used by Irwin (1981, ApJS, 45, 621):
c            log10(U) = SUM{C_j*log10(T)**(j-1)}    for j = 1-->6
c     where T = temperature and C_j are the six polynomial coefficients.
c     Information on the updated partition functions are given by
c     Lawler & Sneden (2002, in preparation).  This routine will be 
c     increasingly invoked as new new data are added to the *newpartdata*
c     array,  until a full replacement for the older partition functions
c     is implemented.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Quants.com'


      iatom = nint(atom)
      iarray = partflag(iatom,k)

      if (level .gt. 500) then
         temp = dlog(dble(level))
      else
         temp = tlog(level)
      endif

      ulog = 0.
      do j=1,6
         ulog = ulog + newpartdata(iarray,j)*temp**(j-1)
      enddo
      partnew = dexp(ulog)

      return
      end                                  



