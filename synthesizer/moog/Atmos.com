
c******************************************************************************
c     these common blocks mainly have quantities related to the model
c     atmosphere and input/output file names/numbers, and program control
c******************************************************************************

      real*8       t(100), theta(100), tkev(100), tlog(100),
     .             pgas(100), ne(100), nhtot(100), numdens(8,2,100),
     .             molweight(100), vturb(100), scont(100), 
     .             kapref(100), kaplam(100), tauref(100), taulam(100), 
     .             kaplamabs(100), kaplamsca(100),
     .             rho(100), rhox(100), xref(100), xdepth(100)
      real*8       elem(95), xabund(95), xabu(95), u(95,4,100)
      real*8       flux, fudge, wavref, abscale, deltaabund
      integer      ntau, jtau5, iunits, itru, iraf, modelnum
      integer      nfparam, nfmodel, nflines, nfslines, nfobs, nftable
      integer      nf1out, nf2out, nf3out, nf4out, nf5out, nf6out,
     .             nf7out, nf8out, nf9out, nf10out,
     .             nfbarklem, nfbarklemUV
      integer      modprintopt, linprintopt, linprintalt,
     .             fluxintopt, plotopt, dampingopt, specfileopt, 
     .             linfileopt, printstrong, linecount, oldcount,
     .             scatopt
      character*80 f1out, f2out, f3out, f4out, f5out, f6out,
     .             f7out, f8out, f9out, f10out,
     .             fparam, fmodel, flines, fslines, fobs, ftable,
     .             fbarklem, fbarklemUV
      character*60 moogpath
      character*2  names(95)
      character*10 modtype
      character*7  control
      character*3  machine
      character*1  silent

      common/atmos/t, theta, tkev, tlog, 
     .             pgas, ne, nhtot, numdens,
     .             molweight, vturb, scont,
     .             kapref, kaplam, tauref, taulam,
     .             kaplamabs, kaplamsca,
     .             rho, rhox, xref, xdepth,
     .             elem, xabund, xabu, u,
     .             flux, fudge, wavref, abscale, deltaabund,
     .             ntau, jtau5, iunits, itru, iraf, modelnum,
     .             nfparam, nfmodel, nflines, nfslines, nfobs, nftable,
     .             nf1out, nf2out, nf3out, nf4out, nf5out, nf6out,
     .             nf7out, nf8out, nf9out, nf10out,
     .             nfbarklem, nfbarklemUV,
     .             modprintopt, linprintopt, linprintalt,
     .             fluxintopt, plotopt, dampingopt, specfileopt, 
     .             linfileopt, printstrong, linecount, oldcount, scatopt

      common/charstuff/ f1out, f2out, f3out, f4out, f5out, f6out,
     .                  f7out, f8out, f9out, f10out,
     .                  fparam, fmodel, flines, fslines, fobs, ftable,
     .                  fbarklem, fbarklemUV,
     .                  moogpath,
     .                  names,
     .                  modtype,
     .                  control,
     .                  machine,
     .                  silent

