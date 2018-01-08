
c******************************************************************************
c     these common block holds parameters and arrays need in the 
c     plotting routines; the second common block puts character variables
c     in their own separate common statement.
c******************************************************************************

      real*8        xxm1, xxb1, xxr1, xxm2, xxb2, xxr2, 
     .              xxm3, xxb3, xxr3,
     .              deltaep, deltarw, deltawv,
     .              p(1000), prot(1000), pmac(1000),
     .              vsini, limbdark, vmac, 
     .              fwhmgauss, fwhmloren,
     .              oldvsini, oldlimbdark, oldvmac, 
     .              oldfwhmgauss, oldfwhmloren,
     .              origvsini, origlimbdark, origvmac, 
     .              origfwhmgauss, origfwhmloren,
     .              half, power, powerrot, powermac,
     .              average, deviate, oldhalf
      real*4        xobs(500000), yobs(500000)
      real*4        xlo, xhi, ylo, yhi, 
     .              origxlo, origxhi, origylo, origyhi,
     .              oldxlo, oldxhi, oldylo, oldyhi,
     .              xadd, yadd, ymult, veladd, deltavel,
     .              origxadd, origyadd, origymult, origveladd,
     .              oldxadd, oldyadd, oldymult, oldveladd,
     .              bigxtic, bigytic, smlxtic, smlytic,
     .              xplotpos, yplotpos, addflux, ydelta
      integer       lount, kount, nsyn, lim1obs, lim2obs,
     .              iterm, iscale, maxline, histoyes, deviations,
     .              syncount, kounter, maxshift

      character*400 abitle
      character*240 isoitle
      character*80  smterm, smtotalterm, smt1, smt2
      character*80  moditle, obsitle, linitle, smitle, popitle,
     .              array, chinfo, errmess, plotroutine
      character*4   whichwin, origwhichwin, oldwhichwin
      character*1   choice, smtype, oldsmtype, origsmtype,
     .              gaussflag, rotateflag, lorenflag, macroflag

      common/pstuff/xxm1, xxb1, xxr1, xxm2, xxb2, xxr2,
     .              xxm3, xxb3, xxr3,
     .              deltaep, deltarw, deltawv,
     .              p, prot, pmac,
     .              vsini, limbdark, vmac,
     .              fwhmgauss, fwhmloren,
     .              oldvsini, oldlimbdark, oldvmac,
     .              oldfwhmgauss, oldfwhmloren,
     .              origvsini, origlimbdark, origvmac,
     .              origfwhmgauss, origfwhmloren,
     .              half, power, powerrot, powermac,
     .              average, deviate, oldhalf,
     .              xobs, yobs,
     .              xlo, xhi, ylo, yhi,
     .              origxlo, origxhi, origylo, origyhi,
     .              oldxlo, oldxhi, oldylo, oldyhi,
     .              xadd, yadd, ymult, veladd, deltavel,
     .              origxadd, origyadd, origymult, origveladd,
     .              oldxadd, oldyadd, oldymult, oldveladd,
     .              bigxtic, bigytic, smlxtic, smlytic,
     .              xplotpos, yplotpos, addflux, ydelta,
     .              lount, kount, nsyn, lim1obs, lim2obs,
     .              iterm, iscale, maxline, histoyes, deviations,
     .              syncount, kounter, maxshift

      common/cstuff/abitle, isoitle,
     .              smterm, smtotalterm, smt1, smt2,
     .              moditle, obsitle, linitle, smitle, popitle,
     .              array, chinfo, errmess, plotroutine,
     .              whichwin, origwhichwin, oldwhichwin,
     .              choice, smtype, oldsmtype, origsmtype,
     .              gaussflag, rotateflag, lorenflag, macroflag



