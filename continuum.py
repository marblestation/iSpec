#!/usr/bin/env python
#import ipdb
import asciitable
import numpy as np
from plotting import *

def find_continuum(spectra, filename=None, wave_increment=0.05, allowed_mean_variation=0.01, max_std_continuum = 0.01, min_mean_continuum = 0.95):
    ##### Parameters
    ## Size of the initial blocks
    #wave_increment = 0.1 # 1 Angstrom (0.1 nm)
    ## Percent of allowed mean variation for increasing the block size
    #allowed_mean_variation = 0.01
    ## Conditions for accepting a final block
    #max_std_continuum = 0.01 # Max. standard deviation for continuum
    #min_mean_continuum = 0.95 # Min. mean for continuum

    wave_base = np.min(spectra['waveobs'])
    wave_top = wave_base + wave_increment
    dirty_continuum_blocks = []

    if filename != None:
        log = open(filename + ".discarded.txt", "w")
        log.write("wave_base\twave_top\tnum_measures\tmean_flux\tstd_flux\n")
    num_obs = 0 # Number of treated observations
    while (wave_base < np.max(spectra['waveobs'])):
        print wave_base, ">",
        # Filter values that belong to the wavelength block
        wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
        # Stats for current block
        mean_flux = np.mean(spectra['flux'][wave_filter])
        std_flux = spectra['flux'][wave_filter].std()
        num_measures = len(spectra['flux'][wave_filter])
        
        if (num_measures > 0 and std_flux < max_std_continuum and mean_flux > min_mean_continuum):
            # Search for a continuous region with similar mean and stdZ
            while True:
                print "*",
                wave_top += wave_increment
                # Stats for the next block
                wave_filter = (spectra['waveobs'] >= wave_top - wave_increment) & (spectra['waveobs'] < wave_top)
                new_mean_flux = np.mean(spectra['flux'][wave_filter])
                new_std_flux = spectra['flux'][wave_filter].std()
                new_num_measures = len(spectra['flux'][wave_filter])
                
                # Compare to the current block, if it is too different do not group it
                mean_increase = np.abs(mean_flux - new_mean_flux) / mean_flux
                std_increase = np.abs(std_flux - new_std_flux) / std_flux
                if (mean_increase > allowed_mean_variation or new_num_measures == 0 or wave_top > spectra['waveobs'].max()):
                    print "!"
                    # Discard the last treated block
                    wave_top -= wave_increment
                    # Final block data stats
                    wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
                    mean_flux = np.mean(spectra['flux'][wave_filter])
                    std_flux = spectra['flux'][wave_filter].std()
                    num_measures = len(spectra['flux'][wave_filter])
                    break
        
            dirty_continuum_blocks.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        else:
            print "Discarded (std = " + str(std_flux) + ", mean = " + str(mean_flux) + ")"
            if filename != None:
                log.write(str(wave_base) + "\t" + str(wave_top) + "\t" + str(num_measures) + "\t" + str(mean_flux) + "\t" + str(std_flux) + "\n")
        
        # Go to next block
        wave_base = wave_top
        wave_top = wave_base + wave_increment
        num_obs += num_measures
    
    if filename != None:
        log.close()        

    ### It can happend that consecutives blocks with different mean_increase are 
    ### selected to belong to the continuum. We can merge them for coherence:
    cleaned_continuum_blocks = []
    i = 0
    if filename != None:
        log = open(filename + ".continuum.txt", "w")
        log.write("wave_base\twave_top\tnum_measures\tmean_flux\tstd_flux\n")
    # For all blocks (except the last one), check the next one is consecutive in wavelengths
    while i < len(dirty_continuum_blocks) - 2:
        j = 0
        # While wave_top of the current is equal to wave_base of the next...
        while ((dirty_continuum_blocks[j+i][1] == dirty_continuum_blocks[j+i+1][0]) and (j < len(dirty_continuum_blocks) - 2 - i)):
            j += 1
        
        wave_base = dirty_continuum_blocks[i][0]
        wave_top = dirty_continuum_blocks[j+i][1]
        
        if j == 1: # No merge needed
            num_measures = dirty_continuum_blocks[i][2]
            mean_flux = dirty_continuum_blocks[i][3]
            std_flux = dirty_continuum_blocks[i][4]
        else:      # Merge and calculate new stats
            wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
            mean_flux = np.mean(spectra['flux'][wave_filter])
            std_flux = spectra['flux'][wave_filter].std()
            num_measures = len(spectra['flux'][wave_filter])
        
        cleaned_continuum_blocks.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        if filename != None:
            log.write(str(wave_base) + "\t" + str(wave_top) + "\t" + str(num_measures) + "\t" + str(mean_flux) + "\t" + str(std_flux) + "\n")
        i += j + 1 # Skip the blocks that have been merged

    if filename != None:
        log.close()

    # Convert result array to numpy array
    continuum_blocks = np.array(cleaned_continuum_blocks,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_blocks


######### TODO: Use spline function to find continuum
## k = 3rd order
#fspline = UnivariateSpline(spectra['waveobs'], spectra['flux'], k=3)

#xs = linspace(np.min(spectra['waveobs']), np.max(spectra['waveobs']), 1000)
#ys = fspline(xs)

#grid = True
#title = None
#ylabel = 'Flux'
#xlabel = 'Wavelength (nm)'
#figure = plt.figure() # 1 inch = 100 pixels

#ax1 = plt.subplot(1, 1, 1)

#if grid:
#    ax1.grid(True, which="both")

#if title != None:
#    ax1.set_title(title, fontsize="10")

#if xlabel != None:
#    ax1.set_xlabel(xlabel, fontsize="10")

#if ylabel != None:
#    ax1.set_ylabel(ylabel, fontsize="10")

##major_tick = np.round((xdata.max() - xdata.min()) / 100, decimals=2)
##minor_tick = major_tick / 2
##ax1.xaxis.set_major_locator(MultipleLocator(major_tick))
##ax1.xaxis.set_minor_locator(MultipleLocator(minor_tick))
##ax1.xaxis.set_major_formatter(FormatStrFormatter("%1.2f"))


#colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
#linestyles = ['', 'steps', 'steps:', '-', '--', ':']

#ax1.plot(spectra['waveobs'], spectra['flux'], lw=1, color=colors[0], linestyle=linestyles[3], marker='', markersize=1, markeredgewidth=0, markerfacecolor=colors[0])
#ax1.plot(xs, ys, lw=1, color=colors[2], linestyle=linestyles[3], marker='', markersize=1, markeredgewidth=0, markerfacecolor=colors[0])

#plt.show()
