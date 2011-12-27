#!/usr/bin/env python
#import ipdb
import asciitable
import numpy as np
from plotting import *

# Find blocks of wavelengths where the fluxes seem to belong to the continuum.
# - It analyses the spectra in blocks
# - The block size is variable in function of 4*fwhm which is derived 
#   from the current wavelength and the resolution
# - For each block, if...
#     a) the median flux is above the continuum model (but not more than 0.08) or below but not more than 0.01
#         * The continuum model can be a fixed flux value or a fitted model (preferable)
#     b) and the standard deviation is less than a given maximum
#   the block is selected
def find_continuum(spectra, resolution, log_filename=None, max_std_continuum = 0.002, continuum_model = 0.95):
    wave_base = np.min(spectra['waveobs'])
    wave_increment = (wave_base / resolution) * 4
    wave_top = wave_base + wave_increment
    dirty_continuum_blocks = []

    if log_filename != None:
        log = open(log_filename, "w")
        log.write("wave_base\twave_top\tnum_measures\tmean_flux\tstd_flux\n")
    
    num_obs = 0 # Number of treated observations
    i = 0
    max_limit = np.max(spectra['waveobs'])
    while (wave_base < max_limit):
        #~ print wave_base, ">"
        # Filter values that belong to the wavelength block
        wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
        # Stats for current block
        mean_flux = np.mean(spectra['flux'][wave_filter])
        std_flux = spectra['flux'][wave_filter].std()
        num_measures = len(spectra['flux'][wave_filter])
        
        # Continuum_model can be a fitted model or a fixed number
        if isinstance(continuum_model, float) or isinstance(continuum_model, int):
            cont_diff = mean_flux - continuum_model
        else:
            cont_diff = mean_flux - continuum_model(wave_base)
        # Flux should be above the continuum model but no more than a given limit
        near_continuum = (cont_diff < 0.08) & (cont_diff > -0.01)
        
        if (num_measures > 0 and std_flux < max_std_continuum and near_continuum):
            dirty_continuum_blocks.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        else:
            #~ print "Discarded (std = " + str(std_flux) + ", mean = " + str(mean_flux) + ")"
            if log_filename != None:
                log.write(str(wave_base) + "\t" + str(wave_top) + "\t" + str(num_measures) + "\t" + str(mean_flux) + "\t" + str(std_flux) + "\n")
        
        # Go to next block
        wave_base = wave_top
        wave_increment = (wave_base / resolution) * 4
        wave_top = wave_base + wave_increment
        num_obs += num_measures
        
        if (i % 200 == 0):
            print "%.2f" % wave_base
        i += 1
    
    if log_filename != None:
        log.close()
    
    continuum_blocks = np.array(dirty_continuum_blocks,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_blocks


# Find blocks of wavelengths where the fluxes seem to belong to the continuum
# but LIMITED to some regions (also called blocks):
# - It analyses the spectra in blocks
# - The block size is variable in function of 4*fwhm which is derived 
#   from the current wavelength and the resolution
# - For each block, if...
#     a) the median flux is above the continuum model (but not more than 0.08) or below but not more than 0.01
#         * The continuum model can be a fixed flux value or a fitted model (preferable)
#     b) and the standard deviation is less than a given maximum
#   the block is selected
def find_continuum_block_limited(spectra, resolution, blocks, log_filename=None, max_std_continuum = 0.002, continuum_model = 0.95):
    
    if log_filename != None:
        log = open(log_filename, "w")
        log.write("wave_base\twave_top\tnum_measures\tmean_flux\tstd_flux\n")
    
    dirty_continuum_blocks = []
    
    for block in blocks:
        wave_base = block['wave_base']
        wave_increment = (wave_base / resolution) * 4
        wave_top = wave_base + wave_increment
        
        num_obs = 0 # Number of treated observations
        i = 0
        max_limit = block['wave_top']
        while (wave_top < max_limit):
            #~ print wave_base, ">"
            # Filter values that belong to the wavelength block
            wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
            # Stats for current block
            mean_flux = np.mean(spectra['flux'][wave_filter])
            std_flux = spectra['flux'][wave_filter].std()
            num_measures = len(spectra['flux'][wave_filter])
            
            # Continuum_model can be a fitted model or a fixed number
            if isinstance(continuum_model, float) or isinstance(continuum_model, int):
                cont_diff = mean_flux - continuum_model
            else:
                cont_diff = mean_flux - continuum_model(wave_base)
            # Flux should be above the continuum model but no more than a given limit
            near_continuum = (cont_diff < 0.08) & (cont_diff > -0.01)
            
            if (num_measures > 0 and std_flux < max_std_continuum and near_continuum):
                dirty_continuum_blocks.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
            else:
                #~ print "Discarded (std = " + str(std_flux) + ", mean = " + str(mean_flux) + ", cont_diff=" + str(cont_diff) + ")"
                if log_filename != None:
                    log.write(str(wave_base) + "\t" + str(wave_top) + "\t" + str(num_measures) + "\t" + str(mean_flux) + "\t" + str(std_flux) + "\n")
            
            # Go to next block
            wave_base = wave_top
            wave_increment = (wave_base / resolution) * 4
            wave_top = wave_base + wave_increment
            num_obs += num_measures
            
            if (i % 200 == 0):
                print "%.2f" % wave_base
            i += 1
        
    if log_filename != None:
        log.close()
    
    continuum_blocks = np.array(dirty_continuum_blocks,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_blocks



# Given a group of continuum blocks of a spectra, merge those that are
# consecutive
def merge_blocks(spectra, dirty_continuum_blocks):
    ### It can happend that consecutives blocks with different mean_increase are 
    ### selected to belong to the continuum. We can merge them for coherence:
    cleaned_continuum_blocks = []
    i = 0
    # For all blocks (except the last one), check the next one is consecutive in wavelengths
    while i < len(dirty_continuum_blocks) - 2:
        j = 0
        # While wave_top of the current is equal to wave_base of the next...
        while ((dirty_continuum_blocks[j+i]['wave_top'] == dirty_continuum_blocks[j+i+1]['wave_base']) and (j < len(dirty_continuum_blocks) - 2 - i)):
            j += 1
        
        wave_base = dirty_continuum_blocks[i]['wave_base']
        wave_top = dirty_continuum_blocks[j+i]['wave_top']
        
        if j == 1: # No merge needed
            num_measures = dirty_continuum_blocks[i]['num_measures']
            mean_flux = dirty_continuum_blocks[i]['mean_flux']
            std_flux = dirty_continuum_blocks[i]['std_flux']
        else:      # Merge and calculate new stats
            wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
            mean_flux = np.mean(spectra['flux'][wave_filter])
            std_flux = spectra['flux'][wave_filter].std()
            num_measures = len(spectra['flux'][wave_filter])
        
        cleaned_continuum_blocks.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        i += j + 1 # Skip the blocks that have been merged

    # Convert result array to numpy array
    continuum_blocks = np.array(cleaned_continuum_blocks,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_blocks


# Considering a cumulative spectra where the 'err' field is the standard 
# deviation of the flux for a group of stars at a given wavelength, 
# identify blocks of wavelength with standard deviation lower than the median
def find_stable_blocks(cumulative_spectra):
    blocks = []
    # Discard blocks that have at least one point with std higher than the median std
    err_limit = np.median(cumulative_spectra['err'])
    total_points = len(cumulative_spectra['err'])
    base = 0
    current = 0

    while current < total_points:
        add_block = False
        while not (current >= total_points or cumulative_spectra['err'][current] > err_limit):
            if not add_block:
                add_block = True
            current += 1
        
        if add_block:
            wave_base = cumulative_spectra['waveobs'][base]
            wave_top = cumulative_spectra['waveobs'][current - 1]
            wave_filter = (cumulative_spectra['waveobs'] >= wave_base) & (cumulative_spectra['waveobs'] < wave_top)
            mean_flux = np.mean(cumulative_spectra['flux'][wave_filter])
            std_flux = cumulative_spectra['flux'][wave_filter].std()
            num_measures = len(cumulative_spectra['flux'][wave_filter])
            
            blocks.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        
        base = current
        current += 1

    # Convert result array to numpy array
    blocks = np.array(blocks,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])

    return blocks

# Read blocks
def read_continuum_blocks(continuum_blocks_filename):
    continuum_blocks = asciitable.read(table=continuum_blocks_filename, delimiter='\t')
    return continuum_blocks

# Write blocks
def write_continuum_blocks(continuum_blocks, continuum_blocks_filename):
    asciitable.write(continuum_blocks, output=continuum_blocks_filename, delimiter='\t')


