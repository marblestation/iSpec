
      real*8 function ucalc (atom,level)                             
c******************************************************************************
c     This routine decodes the partition function data. the source       
c     is the atlas program of kurucz                                     
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Quants.com'
      dimension scale(4) 
      data scale/0.001,0.01,0.1,1.0/ 

      if (level .gt. 500) then
         temp = level               
      else
         temp = t(level)           
      endif

      iatom = nint(atom)
      ion = nint(10.*(atom - float(iatom)))
      j = 4*(iatom-1) + ion + 1
    
      if (ion .eq. 0) then
         chix = xchi1(iatom)      
      elseif (ion .eq. 1) then
         chix = xchi2(iatom)     
      else
         chix = xchi3(iatom)
      endif

      t2000 = chix*2000./11.    
      it = max0(1,min0(9,idint(temp/t2000-0.5))) 
      dt = temp/t2000 - float(it) - 0.5         
      pmin = 1.                                
      i = (it+1)/2                            
      k1 = nudata(i,j)/100000                
      k2 = nudata(i,j) - k1*100000          
      k3 = k2/10                           
      kscale = k2 - k3*10                 
      if (mod(it,2) .eq. 0) then
         p1 = float(k3)*scale(kscale)    
         k1 = nudata(i+1,j)/100000      
         kscale = mod(nudata(i+1,j),10)
         p2 = float(k1)*scale(kscale) 
      else
         p1 = float(k1)*scale(kscale)
         p2 = float(k3)*scale(kscale)            
         if (dt .ge. 0.) go to 13                
         if (kscale .gt. 1.) go to 13           
         kp1 = p1                              
         if (kp1 .ne. idint(p2+0.5)) go to 13 
         pmin = kp1                          
      endif
13    ucalc = dmax1(pmin,p1+(p2-p1)*dt)

      return                                
      end 


