

c******************************************************************************
c     this common block has data relating (so far) to Barklem damping
c     quantities; in the future it may be expanded
c******************************************************************************

      real*8 wavebk(30000), idbk(30000), gammabk(30000), alphabk(30000)
      real*8 gammarad(30000)
      real*8 wavemin, wavemax
      integer firstread, numbark, nummin, nummax

      common/dampdat/ wavebk, idbk, gammabk, alphabk,
     .                gammarad,
     .                wavemin, wavemax,
     .                firstread, numbark, nummin, nummax

