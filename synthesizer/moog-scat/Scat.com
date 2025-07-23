c*****this common block has quantities utilized in the computation
c     of the source function and the MOOG SCAT version of the code 
c     (with associated subroutines Sourcefunc_*, Cdcalc, etc.).

      real*8          dtau1(100),  
     .                wtmu(5), mu(5), Flux_line, 
     .                Flux_cont, adepth, 
     .                Flux_cont_moog, Flux_line_moog,                 
     .                S_cont(100), S_line(100),
     .                J_cont(100), J_line(100),
     .                B_Planck(100)
      integer         mmu

      common/intense/ dtau1,
     .                wtmu, mu, Flux_line, 
     .                Flux_cont, adepth,
     .                Flux_cont_moog, Flux_line_moog,
     .                S_cont, S_line,
     .                J_cont, J_line,
     .                B_Planck,
     .                mmu

