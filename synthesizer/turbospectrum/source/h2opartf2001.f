cc      program testpart
cc
cc      real temp(10),fpartition(10)
cc      data temp /100., 204., 1345.,
cc     &     2400., 3425., 4567., 7300., 15000., 17000.,20000./
cc
cc      call h2opartf2001(10,temp,fpartition)
cc
cc      do i=1,10
cc        print*,temp(i),fpartition(i)
cc      enddo
cc
cc      end

      subroutine h2opartf2001(ntau,temp,fpartition)
*
* Calculates partition function for Joergensen et al. 2001
* line list. 
* 
      implicit none
      integer i,ntau,index
      real Qvib(40),Tvib(40),temp(ntau),fpartition(ntau),Qv,Qrot,a,
     &     Qrb,Qra,B
      data Qvib / 1.0000, 1.0032, 1.0227, 1.0633, 1.1239, 1.2036,
     &            1.3018, 1.4186, 1.5541, 1.7085, 1.8825, 2.0767, 
     &            2.2919, 2.5289, 2.7886, 3.0721, 3.3804, 3.7144, 
     &            4.0752, 4.4637, 4.8807, 5.3272, 5.8038, 6.3112,
     &            6.8499, 7.4203, 8.0227, 8.6572, 9.3238, 10.0225,
     &            10.7531, 11.5153, 12.3085, 13.1325, 13.9864,
     &            14.8698, 15.7819, 16.7218, 17.6889, 18.6821/

      do i=1,40
        Tvib(i)=200.*float(i)
      enddo
      do i=1,ntau
        if (temp(i).le.200.) then
          Qv=1.0
        else if (temp(i).ge.8000.) then
          a=(temp(i)-Tvib(40))/(Tvib(40)-Tvib(30))
          Qv=a*(sqrt(Qvib(40))-sqrt(Qvib(30)))+sqrt(Qvib(40))
          Qv=Qv**2
        else
          index=int(temp(i)/200.)
          print*,index,Qvib(index),Qvib(index+1)
          a=(temp(i)-Tvib(index))/(Tvib(index+1)-Tvib(index))
          Qv=a*(sqrt(Qvib(index+1))-sqrt(Qvib(index)))+sqrt(Qvib(index))
          Qv=Qv**2
        endif
        Qra=sqrt(3.14159/(27.88063*14.52177*9.27771)*
     &                            (temp(i)/1.43877)**3)
        print*,'T, Qra',temp(i),Qra
        B=(1.-sqrt(14.52177*9.27771)/27.88063)*
     &         sqrt(14.52177*9.27771)*1.43877/temp(i)
        print*,'T, B',temp(i),B
        Qrb=(Qra*exp(.25*sqrt(14.52177*9.27771)*1.43877/temp(i))*
     &             (1.+B/12.+B**2*7./480.))*2.
        print*,'T, Qrb',temp(i),Qrb
        fpartition(i)=Qrb*Qv
        print*,'T,partf ',temp(i),fpartition(i)
      enddo
      return
      end
