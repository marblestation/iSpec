
      subroutine calmod                   
c******************************************************************************
c     this program re-calculates BEGN models on a tau5000 scale.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Pstuff.com'


c*****open the files for standard output and the re-calculated model
      nf1out = 20
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21
      lscreen = lscreen + 2
      array = 'RE-CALCULATED MODEL OUTPUT'
      nchars = 26
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)


c*****open and read the model atmosphere
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****do the tau scale conversion
      call opacit(2,wavref)   
      write (nf2out,1001) wavref,(taulam(i),t(i),pgas(i),
     .                ne(i),kaplam(i),i=1,ntau)
      return


c*****format statements
1001  format ('output model for wavelength = ',f5.0/
     .        (1pe11.4,0pf7.0,1p3e11.4))


      end              
