C
      SUBROUTINE MOLMAT(PK,GGH,GGC,GGN,GGO,GGS,GGK,AAC,AAN,AAO,
     &                  AAS,AAK,FH,FC,FN,FO,FS,FK,FE,XXNEN,F,A)
C
C THIS ROUTINE COMPUTES THE ELEMENTS OF MATRIX A AND VECTOR F IN THE
C NEWTON-RAPHSON PROCEDURE FOR DETERMINING THE MOLECULAR EQUILIBRIUM.
C IT IS CALLED BY SUBR. MOL.
C
CCC
C complete double precision version, after some problems encountered in
C TiO pressure computations (P(0I) was wrong).             BPz 890207
C
      DOUBLE PRECISION PK,FH,FC,FN,FO,FS,FK,FE,F,A
      DOUBLE PRECISION FHE,FCE,FNE,FOE,FSE,FKE,H,HH,C,CC,CCC,XN,XNN,
     &                 O,OO,S,SS,XK,GH,GC,GN,GO,GS,GK,XNEN,AC,AN,AO,
     &                 AS,AK
C
      DIMENSION PK(30),F(7),A(7,7)
CCC
C CONVERT TO DOUBLE PRECISION
C
      GH=DBLE(GGH)
      GC=DBLE(GGC)
      GN=DBLE(GGN)
      GO=DBLE(GGO)
      GS=DBLE(GGS)
      GK=DBLE(GGK)
      XNEN=DBLE(XXNEN)
      AC=DBLE(AAC)
      AN=DBLE(AAN)
      AO=DBLE(AAO)
      AS=DBLE(AAS)
      AK=DBLE(AAK)
CCC
      FHE=FH/FE
      FCE=FC/FE
      FNE=FN/FE
      FOE=FO/FE
      FSE=FS/FE
      FKE=FK/FE
      H=1.0+GH+PK(1)+FCE*PK(6)+FNE*PK(13)+FOE*PK(5)+FSE*PK(18)
     *  +FCE*FNE*PK(15)+FCE*FCE*(PK(16)+FCE*PK(20))+FKE*PK(19)
      HH=2.0*FHE*(PK(2)+GH*PK(3)+FCE*FCE*PK(14)+FOE*PK(4))
      C=1.0+GC+FHE*PK(6)+FNE*PK(8)+FOE*PK(7)+FSE*PK(22)+FKE*PK(23)+FHE*
     *  FNE*PK(15)
      CC=2.0*FCE*(PK(9)+FHE*FHE*PK(14)+FHE*PK(16)+FKE*PK(24))
      CCC=3.0*FCE*FCE*(PK(21)+FHE*PK(20))
      XN=1.0+GN+FHE*PK(13)+FCE*PK(8)+FOE*PK(12)+FSE*PK(25)+FKE*PK(26)+
     *   FHE*FCE*PK(15)
      XNN=2.0*FNE*PK(10)
      O=1.0+GO+FHE*PK(5)+FCE*PK(7)+FNE*PK(12)+FSE*PK(28)+FKE*PK(27)+
     *  FHE*FHE*PK(4)
      OO=2.0*FOE*PK(11)
      S=1.0+GS+FHE*PK(18)+FCE*PK(22)+FNE*PK(25)+FOE*PK(28)+FKE*PK(30)
      SS=2.0*FSE*PK(29)
      XK=1.0+GK+FHE*PK(19)+FCE*PK(23)+FNE*PK(26)+FOE*PK(27)+FSE*
     *   PK(30)+FCE*FCE*PK(24)
      F(1)=FH*(H+HH)-1.0
      F(2)=FC*(C+CC+CCC)-AC
      F(3)=FO*(O+OO)-AO
      F(4)=FN*(XN+XNN)-AN
      F(5)=FH*(GH-PK(1)+FHE*GH*PK(3))+FC*GC+FN*GN+FO*GO+FS*GS+FK*GK
     *     +XNEN-FE
      F(6)=FS*(S+SS)-AS
      F(7)=FK*XK-AK
      A(1,1)=H+2.0*HH
      A(1,2)=FHE*(PK(6)+FNE*PK(15)+2.0*FCE*PK(16)+3.0*FCE*FCE*PK(
     *       20)+4.0*FCE*FHE*PK(14))
      A(1,3)=FHE*(PK(5)+2.0*FHE*PK(4))
      A(1,4)=FHE*(PK(13)+FCE*PK(15))
      A(1,5)=-FHE*(2.0*FHE*(PK(2)+GH*PK(3))+FCE*PK(6)+FNE*PK(13)+
     *       FOE*PK(5)+FSE*PK(18)+FKE*PK(19))-2.0*FHE*(FCE*(FNE*
     *       PK(15)+FCE*PK(16))+2.0*FHE*FOE*PK(4))-3.0*FCE*FCE*
     *       FHE*(FCE*PK(20)+2.0*FHE*PK(14))
      A(1,6)=FHE*PK(18)
      A(1,7)=FHE*PK(19)
      A(2,1)=FCE*(PK(6)+FNE*PK(15)+FCE*(4.0*FHE*PK(14)+2.0*PK(16)+
     *       3.0*FCE*PK(20)))
      A(2,2)=C+2.0*CC+3.0*CCC
      A(2,3)=FCE*PK(7)
      A(2,4)=FCE*(PK(8)+FHE*PK(15))
      A(2,5)=-FCE*(2.0*FCE*PK(9)+FHE*PK(6)+FNE*PK(8)+FOE*PK(7)+
     *       FSE*PK(22)+FKE*PK(23))-2.0*FCE*(FHE*FNE*PK(15)+FCE*(3.0
     *       *FCE*PK(21)+2.0*FHE*PK(16)+2.0*FKE*PK(24)))-3.0*FCE
     *       *FCE*FHE*(2.0*FHE*PK(14)+3.0*FCE*PK(20))
      A(2,6)=FCE*PK(22)
      A(2,7)=FCE*(PK(23)+2.0*FCE*PK(24))
      A(3,1)=FOE*(PK(5)+FHE*PK(4)*2.0)
      A(3,2)=FOE*PK(7)
      A(3,3)=O+2.0*OO
      A(3,4)=FOE*PK(12)
      A(3,5)=-FOE*(2.0*FOE*PK(11)+FHE*PK(5)+FCE*PK(7)+FNE*PK(12)+FSE
     *       *PK(28)+FKE*PK(27)+2.0*FHE*FHE*PK(4))
      A(3,6)=FOE*PK(28)
      A(3,7)=FOE*PK(27)
      A(4,1)=FNE*(PK(13)+FCE*PK(15))
      A(4,2)=FNE*(PK(8)+FHE*PK(15))
      A(4,3)=FNE*PK(12)
      A(4,4)=XN+2.0*XNN
      A(4,5)=-FNE*(2.0*FNE*PK(10)+FHE*PK(13)+FCE*PK(8)+FOE*PK(12)+
     *       FSE*PK(25)+FKE*PK(26)+2.0*FHE*FCE*PK(15))
      A(4,6)=FNE*PK(25)
      A(4,7)=FNE*PK(26)
      A(5,1)=GH*(1.0+2.0*FHE*PK(3))-PK(1)
      A(5,2)=GC
      A(5,3)=GO
      A(5,4)=GN
      A(5,5)=-GH*FHE*FHE*PK(3)-1.0
      A(5,6)=GS
      A(5,7)=GK
      A(6,1)=FSE*PK(18)
      A(6,2)=FSE*PK(22)
      A(6,3)=FSE*PK(28)
      A(6,4)=FSE*PK(25)
      A(6,5)=-FSE*(2.0*FSE*PK(29)+FHE*PK(18)+FCE*PK(22)+FNE*PK(25)+
     *       FOE*PK(28)+FKE*PK(30))
      A(6,6)=S+2.0*SS
      A(6,7)=FSE*PK(30)
      A(7,1)=FKE*PK(19)
      A(7,2)=FKE*(PK(23)+2.0*FCE*PK(24))
      A(7,3)=FKE*PK(27)
      A(7,4)=FKE*PK(26)
      A(7,5)=-FKE*(FHE*PK(19)+FCE*PK(23)+FNE*PK(26)+FOE*PK(27)+FSE*
     *       PK(30)+2.0*FCE*FCE*PK(24))
      A(7,6)=FKE*PK(30)
      A(7,7)=XK
CCC
      RETURN
      END
