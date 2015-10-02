
c******************************************************************************
c     this common block has variables related to the lines.  Most
c     input quantities typically have single dimensions, while the 
c     things that are computed for each line at each atmosphere level
c     have double dimensions.  The variables "a", "dopp", and 
c     "kapnu0" are often over-written with plotting data,
c     so leave them alone or suffer unspeakable programming tortures.
c******************************************************************************

      real*8       a(2500,100), dopp(2500,100), kapnu0(2500,100)
      real*8       gf(2500), wave1(2500), atom1(2500), e(2500,2),
     .             chi(2500,3), amass(2500), charge(2500), d0(2500),
     .             dampnum(2500), gf1(2500), width(2500), 
     .             abundout(2500), widout(2500), strength(2500), 
     .             rdmass(2500), gambark(2500), alpbark(2500),
     .             gamrad(2500), wid1comp(2500)
      real*8       kapnu(100), taunu(100), cd(100), sline(100)
      real*8       d(5000), dellam(400), w(100),
     .             rwtab(3000), gftab(3000), gfhold
      real*8       delta, start, sstop, step, contnorm,
     .             oldstart, oldstop, oldstep, olddelta
      real*8       rwlow, rwhigh, rwstep, wavestep, cogatom,
     .             delwave, wave, waveold, st1
      real*8       gammatot, gammav, gammas, gammar
      integer      group(2500), dostrong, gfstyle, lineflag, molflag,
     .             lim1, lim2, mode, nlines, nstrong, ndepths, ncurve,
     .             lim1line, lim2line, n1marker, ntabtot, 
     .             iabatom, iaa, ibb
      character*7  damptype(2500)

      common/linex/a, dopp, kapnu0,   
     .             gf, wave1, atom1, e,
     .             chi, amass, charge, d0,
     .             dampnum, gf1, width, 
     .             abundout, widout, strength,
     .             rdmass, gambark, alpbark, 
     .             gamrad, wid1comp,
     .             kapnu, taunu, cd, sline,
     .             d, dellam, w,
     .             rwtab, gftab, gfhold,
     .             delta, start, sstop, step, contnorm,
     .             oldstart, oldstop, oldstep, olddelta,
     .             rwlow, rwhigh, rwstep, wavestep, cogatom,
     .             delwave, wave, waveold, st1,
     .             gammatot, gammav, gammas, gammar,
     .             group, dostrong, gfstyle, lineflag, molflag,
     .             lim1, lim2, mode, nlines, nstrong, ndepths, ncurve,
     .             lim1line, lim2line, n1marker, ntabtot,
     .             iabatom, iaa, ibb
      common/lindamp/damptype

