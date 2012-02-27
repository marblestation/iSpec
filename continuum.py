"""
    This file is part of Spectra.
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
    
    Spectra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spectra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Spectra.  If not, see <http://www.gnu.org/licenses/>.
"""
#!/usr/bin/env python
#import ipdb
import asciitable
import numpy as np
from plotting import *

# Find regions of wavelengths where the fluxes seem to belong to the continuum.
# - It analyses the spectra in regions
# - The region size is variable in function of 4*fwhm which is derived 
#   from the current wavelength and the resolution (unless a fixed_wave_step is specified)
# - For each region, if...
#     a) the median flux is above the continuum model (but not more than 0.08) or below but not more than 0.01
#         * The continuum model can be a fixed flux value or a fitted model (preferable)
#     b) and the standard deviation is less than a given maximum
#   the region is selected
def find_continuum(spectra, resolution, log_filename=None, max_std_continuum = 0.002, continuum_model = 0.95, fixed_wave_step=None, frame=None):
    
    min_wave = np.min(spectra['waveobs'])
    max_wave = np.max(spectra['waveobs'])
    wave_base = min_wave
    if fixed_wave_step != None:
        wave_increment = fixed_wave_step
    else:
        wave_increment = (wave_base / resolution) * 4
    wave_top = wave_base + wave_increment
    dirty_continuum_regions = []

    if frame != None:
        frame.update_progress(0)
        total_work_progress = max_wave - min_wave
    
    if log_filename != None:
        log = open(log_filename, "w")
        log.write("wave_base\twave_top\tnum_measures\tmean_flux\tstd_flux\n")
    
    i = 0
    max_limit = np.max(spectra['waveobs'])
    while (wave_base < max_limit):
        #~ print wave_base, ">"
        # Filter values that belong to the wavelength region
        wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
        # Stats for current region
        spectra_region = spectra[wave_filter]
        mean_flux = np.mean(spectra_region['flux'])
        std_flux = np.std(spectra_region['flux'])
        num_measures = len(spectra_region['flux'])
        
        # Continuum_model can be a fitted model or a fixed number
        if isinstance(continuum_model, float) or isinstance(continuum_model, int):
            cont_diff = mean_flux - continuum_model
        else:
            cont_diff = mean_flux - continuum_model(wave_base)
        # Flux should be above the continuum model but no more than a given limit
        near_continuum = (cont_diff < 0.08) & (cont_diff > -0.01)
        
        if (num_measures > 0 and std_flux < max_std_continuum and near_continuum):
            last_accepted = True
            dirty_continuum_regions.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        else:
            last_accepted = False
            #~ print "Discarded (std = " + str(std_flux) + ", mean = " + str(mean_flux) + ")"
            if log_filename != None:
                log.write(str(wave_base) + "\t" + str(wave_top) + "\t" + str(num_measures) + "\t" + str(mean_flux) + "\t" + str(std_flux) + "\n")
        
        # Go to next region
        wave_base = wave_top
        if fixed_wave_step != None:
            wave_increment = fixed_wave_step
        else:
            wave_increment = (wave_base / resolution) * 4
        if not last_accepted:
            wave_increment = wave_increment / 2 # Half increase
        wave_top = wave_base + wave_increment
        
        if (i % 200 == 0):
            print "%.2f" % wave_base
            if frame != None:
                current_work_progress = ((wave_base - min_wave) / total_work_progress) * 100
                frame.update_progress(current_work_progress)
    
        i += 1
    
    if log_filename != None:
        log.close()
    
    continuum_regions = np.array(dirty_continuum_regions,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_regions


# Find regions of wavelengths where the fluxes seem to belong to the continuum
# but LIMITED to some regions (also called regions):
# - It analyses the spectra in regions
# - The region size is variable in function of 4*fwhm which is derived 
#   from the current wavelength and the resolution
# - For each region, if...
#     a) the median flux is above the continuum model (but not more than 0.08) or below but not more than 0.01
#         * The continuum model can be a fixed flux value or a fitted model (preferable)
#     b) and the standard deviation is less than a given maximum
#   the region is selected
def find_continuum_on_regions(spectra, resolution, regions, log_filename=None, max_std_continuum = 0.002, continuum_model = 0.95, fixed_wave_step=None, frame=None):
    
    min_wave = np.min(spectra['waveobs'])
    max_wave = np.max(spectra['waveobs'])
    
    if frame != None:
        frame.update_progress(0)
        total_work_progress = max_wave - min_wave
    
    
    if log_filename != None:
        log = open(log_filename, "w")
        log.write("wave_base\twave_top\tnum_measures\tmean_flux\tstd_flux\n")
    
    dirty_continuum_regions = []
    
    for region in regions:
        wave_base = region['wave_base']
        if fixed_wave_step != None:
            wave_increment = fixed_wave_step
        else:
            wave_increment = (wave_base / resolution) * 4
        wave_top = wave_base + wave_increment
        
        i = 0
        max_limit = region['wave_top']
        while (wave_top < max_limit):
            #~ print wave_base, ">"
            # Filter values that belong to the wavelength region
            wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
            # Stats for current region
            spectra_region = spectra[wave_filter]
            mean_flux = np.mean(spectra_region['flux'])
            std_flux = np.std(spectra_region['flux'])
            num_measures = len(spectra_region['flux'])
            
            # Continuum_model can be a fitted model or a fixed number
            if isinstance(continuum_model, float) or isinstance(continuum_model, int):
                cont_diff = mean_flux - continuum_model
            else:
                cont_diff = mean_flux - continuum_model(wave_base)
            # Flux should be above the continuum model but no more than a given limit
            near_continuum = (cont_diff < 0.08) & (cont_diff > -0.01)
            
            if (num_measures > 0 and std_flux < max_std_continuum and near_continuum):
                last_accepted = True
                dirty_continuum_regions.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
            else:
                last_accepted = False
                #~ print "Discarded (std = " + str(std_flux) + ", mean = " + str(mean_flux) + ", cont_diff=" + str(cont_diff) + ")"
                if log_filename != None:
                    log.write(str(wave_base) + "\t" + str(wave_top) + "\t" + str(num_measures) + "\t" + str(mean_flux) + "\t" + str(std_flux) + "\n")
            
            # Go to next region
            wave_base = wave_top
            if fixed_wave_step != None:
                wave_increment = fixed_wave_step
            else:
                wave_increment = (wave_base / resolution) * 4
            if not last_accepted:
                wave_increment = wave_increment / 2 # Half increase
            wave_top = wave_base + wave_increment
            
            
            if (i % 200 == 0):
                #print "%.2f" % wave_base
                if frame != None:
                    current_work_progress = ((wave_base - min_wave) / total_work_progress) * 100
                    frame.update_progress(current_work_progress)
    
            i += 1
        
    if log_filename != None:
        log.close()
    
    continuum_regions = np.array(dirty_continuum_regions,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_regions



# Given a group of continuum regions of a spectra, merge those that are
# consecutive
def merge_regions(spectra, dirty_continuum_regions):
    ### It can happend that consecutives regions with different mean_increase are 
    ### selected to belong to the continuum. We can merge them for coherence:
    cleaned_continuum_regions = []
    i = 0
    # For all regions (except the last one), check the next one is consecutive in wavelengths
    while i < len(dirty_continuum_regions) - 2:
        j = 0
        # While wave_top of the current is equal to wave_base of the next...
        while ((dirty_continuum_regions[j+i]['wave_top'] == dirty_continuum_regions[j+i+1]['wave_base']) and (j < len(dirty_continuum_regions) - 2 - i)):
            j += 1
        
        wave_base = dirty_continuum_regions[i]['wave_base']
        wave_top = dirty_continuum_regions[j+i]['wave_top']
        
        if j == 1: # No merge needed
            num_measures = dirty_continuum_regions[i]['num_measures']
            mean_flux = dirty_continuum_regions[i]['mean_flux']
            std_flux = dirty_continuum_regions[i]['std_flux']
        else:      # Merge and calculate new stats
            wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
            mean_flux = np.mean(spectra['flux'][wave_filter])
            std_flux = spectra['flux'][wave_filter].std()
            num_measures = len(spectra['flux'][wave_filter])
        
        cleaned_continuum_regions.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        i += j + 1 # Skip the regions that have been merged

    # Convert result array to numpy array
    continuum_regions = np.array(cleaned_continuum_regions,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])
    
    return continuum_regions


# Considering a cumulative spectra where the 'err' field is the standard 
# deviation of the flux for a group of stars at a given wavelength, 
# identify regions of wavelength with standard deviation lower than the median
def find_stable_regions(cumulative_spectra):
    regions = []
    # Discard regions that have at least one point with std higher than the median std
    err_limit = np.median(cumulative_spectra['err'])
    total_points = len(cumulative_spectra['err'])
    base = 0
    current = 0

    while current < total_points:
        add_region = False
        while not (current >= total_points or cumulative_spectra['err'][current] > err_limit):
            if not add_region:
                add_region = True
            current += 1
        
        if add_region:
            wave_base = cumulative_spectra['waveobs'][base]
            wave_top = cumulative_spectra['waveobs'][current - 1]
            wave_filter = (cumulative_spectra['waveobs'] >= wave_base) & (cumulative_spectra['waveobs'] < wave_top)
            mean_flux = np.mean(cumulative_spectra['flux'][wave_filter])
            std_flux = cumulative_spectra['flux'][wave_filter].std()
            num_measures = len(cumulative_spectra['flux'][wave_filter])
            
            regions.append((wave_base, wave_top, num_measures, mean_flux, std_flux))
        
        base = current
        current += 1

    # Convert result array to numpy array
    regions = np.array(regions,  dtype=[('wave_base', float), ('wave_top', float), ('num_measures', int), ('mean_flux', float), ('std_flux', float)])

    return regions


