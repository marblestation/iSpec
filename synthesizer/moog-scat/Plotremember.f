
      subroutine plotremember (option)
c******************************************************************************
c     This subroutine declares or stores or recovers the plotting parameters 
c     that were default, or set in the parameter file or were the ones used 
c     in the most recent call to a plotting routine.
c****************************************************************************** 

      include 'Linex.com'
      include 'Pstuff.com'
      integer option


c*****initialize the plot parameters; iscale=0 is simple default;
c     iscale=1 is when these parameters have been read from the
c     parameter file
      if     (option .eq. 0) then
         if (iscale .eq. 0) then
            xlo       = oldstart
            xhi       = oldstop
            ylo       = 0.
            yhi       = 1.1
            veladd    = 0.
            xadd      = 0.
            yadd      = 0.
            ymult     = 1.
            smtype    = 'g'
            fwhmgauss = 0.1
            vsini     = 0.
            limbdark  = 0.
            vmac      = 0.
            fwhmloren = 0.
            whichwin  = '1of1'
         endif


c*****store the original plot parameters
      elseif (option .eq. 1) then
         origxlo       = xlo
         origxhi       = xhi
         origylo       = ylo
         origyhi       = yhi
         origveladd    = veladd
         origxadd      = xadd
         origyadd      = yadd
         origymult     = ymult
         origsmtype    = smtype
         origfwhmgauss = fwhmgauss
         origvsini     = vsini
         origlimbdark  = limbdark
         origvmac      = vmac
         origfwhmloren = fwhmloren
         origwhichwin  = whichwin


c*****re-set the plot parameters to their original values
      elseif (option .eq. 2) then
         xlo       = origxlo
         xhi       = origxhi
         ylo       = origylo
         yhi       = origyhi
         veladd    = origveladd
         xadd      = origxadd
         yadd      = origyadd
         ymult     = origymult
         smtype    = origsmtype
         fwhmgauss = origfwhmgauss
         vsini     = origvsini
         limbdark  = origlimbdark
         vmac      = origvmac
         fwhmloren = origfwhmloren
         whichwin  = origwhichwin


c*****store the plot parameters from the last entry into pltspec
      elseif (option .eq. 3) then
         oldxlo       = xlo
         oldxhi       = xhi
         oldylo       = ylo
         oldyhi       = yhi
         oldveladd    = veladd
         oldxadd      = xadd
         oldyadd      = yadd
         oldymult     = ymult
         oldsmtype    = smtype
         oldfwhmgauss = fwhmgauss
         oldvsini     = vsini
         oldlimbdark  = limbdark
         oldvmac      = vmac
         oldfwhmloren = fwhmloren
         oldwhichwin  = whichwin


c*****re-set plot parameters to values from the last entry into pltspec
      elseif (option .eq. 4) then
         xlo       = oldxlo
         xhi       = oldxhi
         ylo       = oldylo
         yhi       = oldyhi
         veladd    = oldveladd
         xadd      = oldxadd
         yadd      = oldyadd
         ymult     = oldymult
         smtype    = oldsmtype
         fwhmgauss = oldfwhmgauss
         vsini     = oldvsini
         limbdark  = oldlimbdark
         vmac      = oldvmac
         fwhmloren = oldfwhmloren
         whichwin  = oldwhichwin
      else
         stop 'bad option in plotremember; I quit!'
      endif


      return
      end


      

      

