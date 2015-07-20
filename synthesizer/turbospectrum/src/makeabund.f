* Date: Mon, 4 May 1998 13:22:30 +0200 (MET DST)
* From: Bengt Edvardsson <be@astro.uu.se>
* To: Plez@Ferrum.fysik.lu.se
************************************************************************
*
      subroutine makeabund(overall,alpha,helium,rabund,sabund,fixabund,
     &                  abund,amass,aname,isotopfrac)
*
* Create a formatted abundance file
* The programme scales groups of elements automatically
* The input solar abundances are Grevesse et al. logarithmic
*     number abundances on a scale where N(Hydrogen) == 12.00
* The output abundances are on the same scale (there has to be some H)
* A null abundance is signalled as -29.0
* STD?=1 indicates than no abundances have been fiddeled and
* that the first line describes the whole list
* STD?=0 mean that some abundance(s) has been modified by hand
*
* small changes by BPz 04/05-1998 for inclusion in MARCS35
* more changes (isotopfrac inclusion), for SPECTRUM version (BPz 010798)
*
      implicit none
*
      integer natom,nelalpha,i,j
      parameter (natom = 92)
      parameter (nelalpha = 8)
      real amass(natom,0:250),abund(natom),sunabund(natom),
     &     rfract(31:92),tempmass,
     &     fixabund(natom),damass(natom),isotopfrac(natom,0:250),
     &     sunabund_1998(natom),sunabund_2007(natom),checkiso
      real overall,alpha,helium,rabund,sabund,abundr,abunds,hfactor
      integer iel,ialpha,ielalpha(nelalpha)
      character*2 aname(natom),daname(natom)
      logical abund_2007
*
* isotopfrac(22,48) is isotopic fraction of 48Ti 
* sum_over_i>0{ isotopfrac(x,i) } = 1.000 
* isotopfrac(x,0)==1.0; this is used in the case of no isotopes wanted
* in calculation. It will ensure that the total abundance of the species
* will be used.(if for example isotopic factors were included in gf-values)
*

*
CCCC before 24/03-2004*                     C  O  Ne  Mg  Si   S  Ar  Ca  Ti
CCCC before 24/03-2004      data ielalpha / 6, 8, 10, 12, 14, 16, 18, 20, 22 /
*                      O  Ne  Mg  Si   S  Ar  Ca  Ti
      data ielalpha /  8, 10, 12, 14, 16, 18, 20, 22 /

*
      data daname /                                                      | atomic nr
     &           'H ','He','Li','Be','B ','C ','N ','O ','F ',          |  1 -  9
     &           'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar',          | 10 - 18
     &           'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co',          | 19 - 27
     &           'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',          | 28 - 36
     &           'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh',          | 37 - 45
     &           'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',          | 46 - 54
     &           'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',          | 55 - 63
     &           'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',          | 64 - 72
     &           'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl',          | 73 - 81
     &           'Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',          | 82 - 90
     &           'Pa','U '/                                             | 91 - 92
*
      data damass /                                                      | atomic mass
     &   1.008, 4.003, 6.940, 9.010, 10.81, 12.01, 14.01, 16.00, 19.00, |  1 -  9
     &   20.18, 22.99, 24.31, 26.98, 28.08, 30.97, 32.06, 35.45, 39.95, | 10 - 18
     &   39.10, 40.08, 44.96, 47.90, 50.94, 52.00, 54.94, 55.85, 58.93, | 19 - 27
     &   58.70, 63.54, 65.37, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80, | 28 - 36
     &   85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 98.91, 101.0, 102.9, | 37 - 45
     &   106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9, 131.3, | 46 - 54
     &   132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 146.0, 150.4, 152.0, | 55 - 63
     &   157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 170.0, 175.0, 178.5, | 64 - 72
     &   180.9, 183.9, 186.0, 190.0, 192.2, 195.1, 197.0, 200.6, 204.4, | 73 - 81
     &   207.2, 209.0, 210.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0, | 82 - 90
     &   230.4, 238.0 /                                                 | 91 - 92
*
* Solar abundances ref: Grevesse Sauval 1998, Space Science Rev, in press
      data sunabund_1998 /
CCCC before 19/01-2004 changed to agree with Asplund 2003, and Uppsala MARCS version and GRID
CCCC     & 12.00, 10.93,  1.10,  1.40,  2.55,  8.52,  7.92,  8.83,  4.56,   |  1 -  9
CCCC     &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.50,  6.40,   | 10 - 18
CCCC
CCCC
CCCC
CCCC     &  1.69,  0.94,  1.77,  1.66,  2.00,  1.00,  2.24,  1.51,  2.17,   | 46 - 54
CCCC
CCCC     &  1.12, -0.10,  1.14,  0.26,  0.93,  0.00,  1.08,  0.06,  0.88,   | 64 - 72
CCCC     & -0.13,  1.11,  0.28,  1.45,  1.35,  1.80,  1.01,  1.13,  0.90,   | 73 - 81
CCCC
CCCC     & -99.0, -0.47 /                                                   | 91 - 92
     & 12.00, 10.93,  1.10,  1.40,  2.79,  8.41,  7.80,  8.66,  4.48,   |  1 -  9
     &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.28,  6.40,   | 10 - 18
     &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,   | 19 - 27
     &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.41,  2.63,  3.31,   | 28 - 36
     &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -99.0,  1.84,  1.12,   | 37 - 45
     &  1.69,  1.24,  1.77,  0.82,  2.00,  1.00,  2.24,  1.51,  2.17,   | 46 - 54
     &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -99.0,  1.01,  0.51,   | 55 - 63
     &  1.12,  0.35,  1.14,  0.51,  0.93,  0.15,  1.08,  0.06,  0.88,   | 64 - 72
     & -0.13,  0.69,  0.28,  1.45,  1.35,  1.80,  0.85,  1.13,  0.83,   | 73 - 81
     &  1.95,  0.71, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.09,   | 82 - 90
     & -99.0, -0.50 /                                                   | 91 - 92
*
      data sunabund_2007 /
     & 12.00, 10.93,  1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,   |  1 -  9
     &  7.84,  6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,   | 10 - 18
     &  5.08,  6.31,  3.17,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,   | 19 - 27
     &  6.23,  4.21,  4.60,  2.88,  3.58,  2.29,  3.33,  2.56,  3.25,   | 28 - 36
     &  2.60,  2.92,  2.21,  2.58,  1.42,  1.92, -99.0,  1.84,  1.12,   | 37 - 45
     &  1.66,  0.94,  1.77,  1.60,  2.00,  1.00,  2.19,  1.51,  2.24,   | 46 - 54
     &  1.07,  2.17,  1.13,  1.70,  0.58,  1.45, -99.0,  1.00,  0.52,   | 55 - 63
     &  1.11,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08,  0.06,  0.88,   | 64 - 72
     & -0.17,  1.11,  0.23,  1.25,  1.38,  1.64,  1.01,  1.13,  0.90,   | 73 - 81
     &  2.00,  0.65, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.06,   | 82 - 90
     & -99.0, -0.52 /                                                   | 91 - 92
*
*
* Fractional r-process abundance for Ga-Bi (r+s simply assumed == 100%) | Date 2000-01-18
* (Note: Ga-Sr (31-38) was just copied from Kaeppeler et al. 1989, below)
* s-process from Stellar models: Arlandini C., Kaeppeler F., Wisshak K.,
* Gallino R., Busso M., Straniero O., 1999, Astrophys J. 525, 886-900
*   Fractions corrected to the revised meteoritic abundances
*   of Grevesse N., Sauval A.J. 1998, Space Science Review 85, 161-174
      data rfract /
*    & H  - F                                                           |  1 -  9
*    & Ne - Ar                                                          | 10 - 18
*    & K  - Co                                                          | 19 - 27
     &                       .43, .47, .81, .85, .39, .47,              | 28 - 36
     &        .41, .11, .08, .17, .15, .50,-.99, .68, .86,              | 37 - 45
     &        .54, .80, .48, .65, .35, .75, .83, .80, .80,              | 46 - 54
     &        .85, .19, .38, .23, .51, .44,-.99, .71, .93,              | 55 - 63
     &        .85, .93, .85, .92, .83, .87, .67, .80, .44,              | 64 - 72
     &        .59, .44, .91, .91, .99, .95, .94, .41, .24,              | 73 - 81
     &        .54, .95,-.99,-.99,-.99,-.99,-.99,-.99, 1.0,              | 82 - 90
     &       -.99, 1.0 /                                                | 91 - 92

c* fractional r-process abundance for Cu-U (r+s simply assumed == 100%)
c* r-process from Kaeppeler F., Beer H., Wisshak K. 1989, Reports on
c*   Progress in Physics 52, 945. Fraction found by division by the
c*   total meteoritic abundances given by Anders Grevesse 1989, Geochim.
c*   Cosmochim. Acta 53, 197. Th+U guess 1.0
c* (Perhaps some more data in Beer et al 1997, ApJ 474,843)
c* (Sneden et al 1998 ApJ 496,235: Sm 76%, Eu 97%)
c      data rfract /
c*    & H  - F                                                           |  1 -  9
c*    & Ne - Ar                                                          | 10 - 18
c*    & K  - Co                                                          | 19 - 27
c     &             .21, .35, .43, .47, .81, .85, .39, .47,              | 28 - 36
c     &        .41, .11, .28, .18, .16, .25,-.99, .50, .83,              | 37 - 45
c     &        .56, .88, .57, .64, .36, .99, .82, .93, .84,              | 46 - 54
c     &        .85, .12, .25, .23, .56, .51,-.99, .68, .93,              | 55 - 63
c     &        .81, .91, .90, .92, .84, .84, .65, .83, .53,              | 64 - 72
c     &        .35, .46, .90, .96, .98, .96, .92, .42, .28,              | 73 - 81
c     &        .20, .64,-.99,-.99,-.99,-.99,-.99,-.99, 1.0,              | 82 - 90
c     &       -.99, 1.0 /                                                | 91 - 92
*
*
*
cc* Solar abundances ref: Grevesse Noels Sauval 1996 ASP Conf 99, 117:
cc      data sunabund /
cc     & 12.00, 10.93,  1.16,  1.15,  2.60,  8.55,  7.97,  8.87,  4.56,   |  1 -  9
cc     &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.50,  6.52,   | 10 - 18
cc     &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,   | 19 - 27
cc     &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.38,  2.63,  3.23,   | 28 - 36
cc     &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -29.0,  1.84,  1.12,   | 37 - 45
cc     &  1.69,  0.94,  1.77,  1.66,  2.00,  1.00,  2.24,  1.51,  2.23,   | 46 - 54
cc     &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -29.0,  1.01,  0.51,   | 55 - 63
cc     &  1.12, -0.10,  1.14,  0.26,  0.93,  0.00,  1.08,  0.76,  0.88,   | 64 - 72
cc     & -0.13,  1.11,  0.28,  1.45,  1.35,  1.80,  1.01,  1.17,  0.90,   | 73 - 81
cc     &  1.95,  0.71, -29.0, -29.0, -29.0, -29.0, -29.0, -29.0,  0.09,   | 82 - 90
cc     & -29.0, -0.47 /                                                   | 91 - 92
cc*
cc* fractional r-process abundance for Cu-U (r+s assumed == 100%)
cc* r-process from Kaeppeler F., Beer H., Wisshak K. 1989, Reports on
cc*   Progress in Physics 52, 945. Fraction found by division by the
cc*   total meteoritic abundances given by Anders Grevesse 1989, Geochim.
cc*   Cosmochim. Acta 53, 197 (Xe=1.34 Dy=1.08 Lu=1.69 !) Th+U guess 1.0
cc*   Isotopic fractions from nuclide chart U.S. Dept of Energy 1988
cc*   Ba r-process=1.46 from McWilliam et al. 1998, AJ 115, 1640
cc* (Perhaps some more data in Beer et al 1997, ApJ 474,843)
cc      data rfract /
cc*    & H  - F                                                           |  1- 9
cc*    & Ne - Ar                                                          | 10-18
cc*    & K  - Co                                                          | 19-27
cc     &             .21, .25, .43, .50, .81, .78, .39, .59,              | 28-36
cc     &        .41, .11, .28, .18, .16, .34,-.99, .55, .83,              | 37-45
cc     &        .73, .88, .54, .64, .37, .50, .75, .93, .50,              | 46-54
cc     &        .85, .21, .25, .11, .56, .61,-.99, .85, .93,              | 55-63
cc     &        .82, .91, .50, .92, .84, .84, .59, .50, .53,              | 64-72
cc     &        .35, .46, .90, .96, .98, .90, .92, .67, .28,              | 73-81
cc     &        .20, .64,-.99,-.99,-.99,-.99,-.99,-.99, 1.0,              | 82-90
cc     &       -.99, 1.0 /                                                | 91-92
*
cc      print *,'Overall metallicity? (-1.0 = Solar/10)'
cc      read(*,*) overall
cc      print *,'Helium abundance? (0.0 = Solar)'
cc      read(*,*) helium
cc      print *,'Alpha element scaling [alpha/Fe]? (+0.3 = Solar*2)'
cc      read(*,*) alpha
cc      print *,'r-element abundance [r/Fe]? (0.0 = Solar)'
cc      read(*,*) rabund
cc      print *,'s-element abundance [s/Fe]? (0.0 = Solar)'
cc      read(*,*) sabund
*
* chose set of abundances
      abund_2007=.true.
      if (abund_2007) then
        do iel=1,natom
          sunabund(iel)=sunabund_2007(iel)
        enddo
      else
        do iel=1,natom
          sunabund(iel)=sunabund_1998(iel)
        enddo
      endif
*
      do iel=1,92
        aname(iel)=daname(iel)
* mass of the default isotopic mix
* This is superseded by the getisotopmass.f routine below.
* BPz 18/06-2012
        amass(iel,0)=damass(iel)
      enddo
      abund(1)=sunabund(1)
      abund(2)=sunabund(2)
* overall metallicity
      do iel=3,natom
        if(sunabund(iel).gt.-29.0) then
          abund(iel)=sunabund(iel)+overall
        else
          abund(iel)=-29.0
        endif
      enddo
*
      if(alpha.ne.0.0) then
* alpha element abundance [alpha/Fe]
        do ialpha=1,nelalpha
          abund(ielalpha(ialpha))=abund(ielalpha(ialpha))+alpha
        enddo
      endif
*
* helium deficient?
      if(helium.ne.0.0) then
        abund(2)=abund(2)+helium
      endif
*
* r-element or s-element anomaly?
      if(rabund.ne.0.0 .or. sabund.ne.0.0) then
      print*,aname(33),sunabund(33),abund(33),rfract(33)
        do iel=31,92  
cc        do iel=29,92
          if(sunabund(iel).gt.-29.0) then
            abundr=10.**rabund*rfract(iel)*10.**abund(iel)
            abunds=10.**sabund*(1.-rfract(iel))*10.**abund(iel)
      if(iel.eq.33) print*,rabund,sabund,abundr,abunds,abund(iel),iel
            abund(iel)=alog10(abundr+abunds)
          endif
        enddo
      endif
*
* If applicable, overrides any prescription by specific input value
      do iel=1,92
cc        abund(iel)=max(abund(iel),fixabund(iel))
        if (fixabund(iel).gt.-98.) then
          abund(iel)=fixabund(iel)
        endif
      enddo
*
* normalize to H=12.00
      hfactor=12.00-abund(1)
      if(hfactor.ne.0.0) then
        do iel=1,natom
          abund(iel)=abund(iel)+hfactor
        enddo
      endif
*
* isotopic ratios
      do i=1,250
        do j=1,92
          isotopfrac(j,i)=0.0
        enddo
      enddo
      do j=1,92
        isotopfrac(j,0)=1.0
      enddo
*
      isotopfrac(1,1)=0.999966
      isotopfrac(1,2)=0.000034
      isotopfrac(2,3)=0.000142
      isotopfrac(2,4)=0.999858
      isotopfrac(3,6)=0.075
      isotopfrac(3,7)=0.925
      isotopfrac(4,9)=1.0
      isotopfrac(5,10)=0.199
      isotopfrac(5,11)=0.801
      isotopfrac(6,12)=0.9890
      isotopfrac(6,13)=0.0110
      isotopfrac(7,14)=0.99634
      isotopfrac(7,15)=0.00366
      isotopfrac(8,16)=0.99762
      isotopfrac(8,17)=0.00038
      isotopfrac(8,18)=0.00200
      isotopfrac(9,19)=1.00
      isotopfrac(10,20)=0.9299
      isotopfrac(10,21)=0.00226
      isotopfrac(10,22)=0.0679
      isotopfrac(11,23)=1.00
      isotopfrac(12,24)=0.7899
      isotopfrac(12,25)=0.1000
      isotopfrac(12,26)=0.1101
      isotopfrac(13,27)=1.00
      isotopfrac(14,28)=0.9223
      isotopfrac(14,29)=0.0467
      isotopfrac(14,30)=0.0310
      isotopfrac(15,31)=1.00
      isotopfrac(16,32)=0.9502
      isotopfrac(16,33)=0.0075
      isotopfrac(16,34)=0.0421
      isotopfrac(16,36)=0.0002
      isotopfrac(17,35)=0.7577
      isotopfrac(17,37)=0.2423
      isotopfrac(18,36)=0.842
      isotopfrac(18,38)=0.158
      isotopfrac(19,39)=0.932581
      isotopfrac(19,40)=0.0001167
      isotopfrac(19,41)=0.067302
      isotopfrac(20,40)=0.96941
      isotopfrac(20,42)=0.00647
      isotopfrac(20,43)=0.00135
      isotopfrac(20,44)=0.02086
      isotopfrac(20,46)=0.00004
      isotopfrac(20,48)=0.00187
      isotopfrac(21,45)=1.00
      isotopfrac(22,46)=0.080
      isotopfrac(22,47)=0.073
      isotopfrac(22,48)=0.738
      isotopfrac(22,49)=0.055
      isotopfrac(22,50)=0.054
      isotopfrac(23,50)=0.00250
      isotopfrac(23,51)=0.99750
      isotopfrac(24,50)=0.04345
      isotopfrac(24,52)=0.83789
      isotopfrac(24,53)=0.09501
      isotopfrac(24,54)=0.02365
      isotopfrac(25,55)=1.00
      isotopfrac(26,54)=0.058
      isotopfrac(26,56)=0.9172
      isotopfrac(26,57)=0.022
      isotopfrac(26,58)=0.0028
      isotopfrac(27,59)=1.00
      isotopfrac(28,58)=0.6827
      isotopfrac(28,60)=0.2610
      isotopfrac(28,61)=0.0113
      isotopfrac(28,62)=0.0359
      isotopfrac(28,64)=0.0091
      isotopfrac(29,63)=0.6917
      isotopfrac(29,65)=0.3083
      isotopfrac(30,64)=0.4863
      isotopfrac(30,66)=0.2790
      isotopfrac(30,67)=0.0410
      isotopfrac(30,68)=0.1875
      isotopfrac(30,70)=0.0062
      isotopfrac(31,69)=0.60108
      isotopfrac(31,71)=0.39892
      isotopfrac(32,70)=0.205
      isotopfrac(32,72)=0.274
      isotopfrac(32,73)=0.078
      isotopfrac(32,74)=0.365
      isotopfrac(32,76)=0.078
      isotopfrac(33,75)=1.00
      isotopfrac(34,74)=0.0088
      isotopfrac(34,76)=0.090
      isotopfrac(34,77)=0.076
      isotopfrac(34,78)=0.236
      isotopfrac(34,80)=0.497
      isotopfrac(34,82)=0.092
      isotopfrac(35,79)=0.5069
      isotopfrac(35,81)=0.4931
      isotopfrac(36,78)=0.00339
      isotopfrac(36,80)=0.0222
      isotopfrac(36,82)=0.1145
      isotopfrac(36,83)=0.1147
      isotopfrac(36,84)=0.5711
      isotopfrac(36,86)=0.1742
      isotopfrac(37,85)=0.72165
      isotopfrac(37,87)=0.27835
      isotopfrac(38,84)=0.0056
      isotopfrac(38,86)=0.0986
      isotopfrac(38,87)=0.0700
      isotopfrac(38,88)=0.8258
      isotopfrac(39,89)=1.00
      isotopfrac(40,90)=0.5145
      isotopfrac(40,91)=0.1122
      isotopfrac(40,92)=0.1715
      isotopfrac(40,94)=0.1738
      isotopfrac(40,96)=0.0280
      isotopfrac(41,93)=1.00
      isotopfrac(42,92)=0.1484
      isotopfrac(42,94)=0.0925
      isotopfrac(42,95)=0.1592
      isotopfrac(42,96)=0.1668
      isotopfrac(42,97)=0.0955
      isotopfrac(42,98)=0.2413
      isotopfrac(42,100)=0.0963
      isotopfrac(44,96)=0.0552
      isotopfrac(44,98)=0.0188
      isotopfrac(44,99)=0.127
      isotopfrac(44,100)=0.126
      isotopfrac(44,101)=0.170
      isotopfrac(44,102)=0.316
      isotopfrac(44,104)=0.187
      isotopfrac(45,103)=1.00
      isotopfrac(46,102)=0.01020
      isotopfrac(46,104)=0.1114
      isotopfrac(46,105)=0.2233
      isotopfrac(46,106)=0.2733
      isotopfrac(46,108)=0.2646
      isotopfrac(46,110)=0.1172
      isotopfrac(47,107)=0.51839
      isotopfrac(47,109)=0.48161
      isotopfrac(48,106)=0.0125
      isotopfrac(48,108)=0.0089
      isotopfrac(48,110)=0.1249
      isotopfrac(48,111)=0.1280
      isotopfrac(48,112)=0.2413
      isotopfrac(48,113)=0.1222
      isotopfrac(48,114)=0.2873
      isotopfrac(48,116)=0.0749
      isotopfrac(49,113)=0.043
      isotopfrac(49,115)=0.957
      isotopfrac(50,112)=0.00973
      isotopfrac(50,114)=0.00659
      isotopfrac(50,115)=0.00339
      isotopfrac(50,116)=0.14538
      isotopfrac(50,117)=0.07672
      isotopfrac(50,118)=0.24217
      isotopfrac(50,119)=0.08587
      isotopfrac(50,120)=0.32596
      isotopfrac(50,122)=0.04632
      isotopfrac(50,124)=0.05787
      isotopfrac(51,121)=0.57362
      isotopfrac(51,123)=0.42638
      isotopfrac(52,120)=0.0009
      isotopfrac(52,122)=0.0257
      isotopfrac(52,123)=0.0089
      isotopfrac(52,124)=0.0476
      isotopfrac(52,125)=0.0710
      isotopfrac(52,126)=0.1889
      isotopfrac(52,128)=0.3173
      isotopfrac(52,130)=0.3397
      isotopfrac(53,127)=1.00
      isotopfrac(54,124)=0.00121
      isotopfrac(54,126)=0.00108
      isotopfrac(54,128)=0.0219
      isotopfrac(54,129)=0.2734
      isotopfrac(54,130)=0.0435
      isotopfrac(54,131)=0.2169
      isotopfrac(54,132)=0.2650
      isotopfrac(54,134)=0.0976
      isotopfrac(54,136)=0.0794
      isotopfrac(55,133)=1.00
      isotopfrac(56,130)=0.00106
      isotopfrac(56,132)=0.00101
      isotopfrac(56,134)=0.02417
      isotopfrac(56,135)=0.06592
      isotopfrac(56,136)=0.07854
      isotopfrac(56,137)=0.1123
      isotopfrac(56,138)=0.7170
      isotopfrac(57,138)=0.00089
      isotopfrac(57,139)=0.99911
      isotopfrac(58,136)=0.0019
      isotopfrac(58,138)=0.0025
      isotopfrac(58,140)=0.8848
      isotopfrac(58,142)=0.1108
      isotopfrac(59,141)=1.00
      isotopfrac(60,142)=0.2713
      isotopfrac(60,143)=0.1218
      isotopfrac(60,144)=0.2380
      isotopfrac(60,145)=0.0830
      isotopfrac(60,146)=0.1719
      isotopfrac(60,148)=0.0576
      isotopfrac(60,150)=0.0564
      isotopfrac(62,144)=0.031
      isotopfrac(62,147)=0.150
      isotopfrac(62,148)=0.113
      isotopfrac(62,149)=0.138
      isotopfrac(62,150)=0.074
      isotopfrac(62,152)=0.267
      isotopfrac(62,154)=0.227
      isotopfrac(63,151)=0.478
      isotopfrac(63,153)=0.522
      isotopfrac(64,152)=0.0020
      isotopfrac(64,154)=0.0218
      isotopfrac(64,155)=0.1480
      isotopfrac(64,156)=0.2047
      isotopfrac(64,157)=0.1565
      isotopfrac(64,158)=0.2484
      isotopfrac(64,160)=0.2186
      isotopfrac(65,159)=1.00
      isotopfrac(66,156)=0.00056
      isotopfrac(66,158)=0.00096
      isotopfrac(66,160)=0.0234
      isotopfrac(66,161)=0.1891
      isotopfrac(66,162)=0.2551
      isotopfrac(66,163)=0.2490
      isotopfrac(66,164)=0.2819
      isotopfrac(67,165)=1.00
      isotopfrac(68,162)=0.0014
      isotopfrac(68,164)=0.0161
      isotopfrac(68,166)=0.336
      isotopfrac(68,167)=0.2295
      isotopfrac(68,168)=0.268
      isotopfrac(68,170)=0.149
      isotopfrac(69,169)=1.00
      isotopfrac(70,168)=0.0013
      isotopfrac(70,170)=0.0305
      isotopfrac(70,171)=0.143
      isotopfrac(70,172)=0.219
      isotopfrac(70,173)=0.1612
      isotopfrac(70,174)=0.318
      isotopfrac(70,176)=0.127
      isotopfrac(71,175)=0.9741
      isotopfrac(71,176)=0.0259
      isotopfrac(72,174)=0.00162
      isotopfrac(72,176)=0.05206
      isotopfrac(72,177)=0.18606
      isotopfrac(72,178)=0.27297
      isotopfrac(72,179)=0.13629
      isotopfrac(72,180)=0.35100
      isotopfrac(73,180)=0.00012
      isotopfrac(73,181)=0.99988
      isotopfrac(74,180)=0.0013
      isotopfrac(74,182)=0.263
      isotopfrac(74,183)=0.143
      isotopfrac(74,184)=0.3067
      isotopfrac(74,186)=0.286
      isotopfrac(75,185)=0.3740
      isotopfrac(75,187)=0.6260
      isotopfrac(76,184)=0.00018
      isotopfrac(76,186)=0.0158
      isotopfrac(76,187)=0.016
      isotopfrac(76,188)=0.133
      isotopfrac(76,189)=0.161
      isotopfrac(76,190)=0.264
      isotopfrac(76,192)=0.410
      isotopfrac(77,191)=0.373
      isotopfrac(77,193)=0.627
      isotopfrac(78,190)=0.000127
      isotopfrac(78,192)=0.0078
      isotopfrac(78,194)=0.329
      isotopfrac(78,195)=0.338
      isotopfrac(78,196)=0.252
      isotopfrac(78,198)=0.0719
      isotopfrac(79,197)=1.00
      isotopfrac(80,196)=0.001534
      isotopfrac(80,198)=0.09968
      isotopfrac(80,199)=0.16873
      isotopfrac(80,200)=0.23096
      isotopfrac(80,201)=0.13181
      isotopfrac(80,202)=0.29863
      isotopfrac(80,204)=0.06865
      isotopfrac(81,203)=0.29524
      isotopfrac(81,205)=0.70476
      isotopfrac(82,204)=0.0194
      isotopfrac(82,206)=0.1912
      isotopfrac(82,207)=0.2062
      isotopfrac(82,208)=0.5831
      isotopfrac(83,209)=1.00
      isotopfrac(90,232)=1.00
      isotopfrac(92,235)=0.007200
      isotopfrac(92,238)=0.992745

* get the mass of each isotope
* and recompute the default mix mass
* BPz 12/06-2012
      call getisotopmass(amass)
      do i=1,92
        tempmass=amass(i,0)
        amass(i,0)=0.
        do j=1,250
          amass(i,0)=amass(i,0)+isotopfrac(i,j)*amass(i,j)
        enddo
        if (amass(i,0).eq.0.) then
          amass(i,0)=tempmass
        endif
      enddo

*      print*,' Atomic masses are masses of the default terrestrial',
*     & ' isotopic mixture. Cf makeabund.f.'
*
*
* print out

c      open(90,file='abundances.tmp',status='unknown',form='formatted')
c      write(90,1020) overall,alpha,helium,rabund,sabund
c 1020 format('* STD? [Fe/H] [He/Fe] [alpha/Fe] [r/Fe] [s/Fe]',/,
c     &      '   1 ', 6f8.2)
c      write(90,1000) (aname(iel),iel= 1,10)
c      write(90,1010) (abund(iel),iel= 1,10)
c      write(90,1000) (aname(iel),iel=11,20)
c      write(90,1010) (abund(iel),iel=11,20)
c      write(90,1000) (aname(iel),iel=21,30)
c      write(90,1010) (abund(iel),iel=21,30)
c      write(90,1000) (aname(iel),iel=31,40)
c      write(90,1010) (abund(iel),iel=31,40)
c      write(90,1000) (aname(iel),iel=41,50)
c      write(90,1010) (abund(iel),iel=41,50)
c      write(90,1000) (aname(iel),iel=51,60)
c      write(90,1010) (abund(iel),iel=51,60)
c      write(90,1000) (aname(iel),iel=61,70)
c      write(90,1010) (abund(iel),iel=61,70)
c      write(90,1000) (aname(iel),iel=71,80)
c      write(90,1010) (abund(iel),iel=71,80)
c      write(90,1000) (aname(iel),iel=81,90)
c      write(90,1010) (abund(iel),iel=81,90)
c      write(90,1000) (aname(iel),iel=91,92)
c      write(90,1010) (abund(iel),iel=91,92)
c 1000 format('*   ',10(3x,a2,2x))
c 1010 format(2x,10f7.2)
c*
c      print *,'makeabund; check Output in file abundances.tmp'
      return
      end
