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
import numpy as np
import pyfits
import scipy.ndimage as ndi
from interpolate import *
from plotting import *
from common import *
from radial_velocity import *

# Calculate the FWHM of the gaussian needed to convert
# a spectra from one resolution to another
def get_fwhm(lambda_peak, from_resolution, to_resolution):
    if from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are equal or bigger than original")
    from_delta_lambda = (1.0*lambda_peak) / from_resolution
    to_delta_lambda = (1.0*lambda_peak) / to_resolution
    fwhm = np.sqrt(to_delta_lambda**2 - from_delta_lambda**2)
    return fwhm

# Convert the FWHM to the sigma for building a gaussian
def get_sigma(fwhm):
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    return sigma


## Spectra resolution smoothness/degradation. Procedure:
#   1. Define a bin per measure which marks the wavelength range that it covers
#   2. For each point, identify the window segment to convolve by using the bin widths and the FWHM
#   2. Build a gaussian using the sigma value and the wavelength values of the spectra window
#   3. Convolve the spectra window with the gaussian and save the convolved value
# If "to_resolution" is not specified or its equal to "from_resolution", then the spectra 
# is convolved with the instrumental gaussian defined by "from_resolution"
def convolve_spectra(spectra, from_resolution, to_resolution=None, frame=None):
    if to_resolution != None and from_resolution < to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")

    total_points = len(spectra)
    convolved_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    convolved_spectra['waveobs'] = spectra['waveobs']
    convolved_spectra['err'] = spectra['err']
    
    if frame != None:
        frame.update_progress(0)
    
    flux = spectra['flux']
    # Consider the wavelength of the measurements as the center of the bins
    waveobs = spectra['waveobs']
    # Calculate the wavelength distance between the center of each bin
    wave_distance = waveobs[1:] - waveobs[:-1]
    # Define the edge of each bin as half the wavelength distance to the bin next to it
    edges_tmp = waveobs[:-1] + 0.5 * (wave_distance)
    # Define the edges for the first and last measure which where out of the previous calculations
    first_edge = waveobs[0] - 0.5*wave_distance[0]
    last_edge = waveobs[-1] + 0.5*wave_distance[-1]
    # Build the final edges array
    edges = np.array([first_edge] + edges_tmp.tolist() + [last_edge])
    
    # Bin width
    bin_width = edges[1:] - edges[:-1]          # width per pixel
    
    # FWHM of the gaussian for the given resolution
    if to_resolution == None or to_resolution == from_resolution:
        # Convolve using instrumental resolution (smooth but not degrade)
        fwhm = waveobs / from_resolution
    else:
        # Degrade resolution
        fwhm = get_fwhm(waveobs, from_resolution, to_resolution)
    sigma = get_sigma(fwhm)
    # Convert from wavelength units to bins
    fwhm_bin = fwhm / bin_width
    
    # Round number of bins per FWHM
    nbins = np.ceil(fwhm_bin) #npixels
    
    # Number of measures
    nwaveobs = len(waveobs)
    
    # In theory, len(nbins) == len(spectra)
    for i in np.arange(len(nbins)):
        current_nbins = 2 * nbins[i] # Each side
        current_center = waveobs[i] # Center
        current_sigma = sigma[i]
        
        # Find lower and uper index for the gaussian, taking care of the current spectra limits
        lower_pos = int(max(0, i - current_nbins))
        upper_pos = int(min(nwaveobs, i + current_nbins + 1))
        
        # Select only the flux values for the segment that we are going to convolve
        flux_segment = flux[lower_pos:upper_pos+1]
        waveobs_segment = waveobs[lower_pos:upper_pos+1]
        
        nsegments = len(flux_segment)
        
        # Build the gaussian corresponding to the instrumental spread function
        gaussian = np.exp(- ((waveobs_segment - current_center)**2) / (2*current_sigma**2)) / np.sqrt(2*np.pi*current_sigma**2)
        gaussian = gaussian / np.sum(gaussian)
        
        # Convolve the current position by using the segment and the gaussian
        weighted_flux = flux_segment * gaussian
        current_convolved_flux = weighted_flux.sum()
        
        convolved_spectra['flux'][i] = current_convolved_flux
        
        if (i % 1000 == 0):
            if frame != None:
                current_work_progress = (i*1.0 / total_points) * 100
                frame.update_progress(current_work_progress)
                #print "%.2f" % convolved_spectra['waveobs'][i]                
    
    return convolved_spectra

## Degradation of the spectra resolution. Procedure for each flux value:
#   1. Define a window based on the fwhm size
#   2. Build a gaussian using the sigma value and the wavelength values of the spectra window
#   3. Convolve the spectra window with the gaussian and save the convolved value
# WARNING: implementation poorly efficient, better use convolve_spectra method
def degrade_spectra_resolution(spectra, from_resolution, to_resolution, frame=None):
    if from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are equal or bigger than original")

    total_points = len(spectra['waveobs'])
    convolved_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    convolved_spectra['waveobs'] = spectra['waveobs']
    
    if frame != None:
        frame.update_progress(0)
    
    for i in np.arange(total_points):
        lambda_peak = spectra['waveobs'][i] # Current lambda (wavelength) to be modified
        fwhm = get_fwhm(lambda_peak, from_resolution, to_resolution) # Calculate the needed fwhm for this wavelength
        
        # Only work with a limited window considering 3 times the fwhm in each side of the current
        # position to be modified and saved in the convolved spectra
        wave_filter = (spectra['waveobs'] >= lambda_peak - 3*fwhm) & (spectra['waveobs'] <= lambda_peak + 3*fwhm)
        spectra_window = spectra[wave_filter]
        min = np.min(spectra_window['waveobs'])
        max = np.max(spectra_window['waveobs'])
        # Find the index position of the current lambda peak in the window
        lambda_peak_win_index = spectra_window['waveobs'].searchsorted(lambda_peak)
        
        # Check that the window is large enough (i.e. in the edges of the original spectra, maybe we do not have enough)
        # At least it should be 3 times FWHM in total
        if np.round(max - min, 2) >= np.round(3*fwhm, 2):
            sigma = get_sigma(fwhm)
            
            # Construct the gaussian
            gaussian = np.exp(- ((spectra_window['waveobs'] - lambda_peak)**2) / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
            gaussian = gaussian / np.sum(gaussian)
            
            # Convolve the current position
            convolved_value = 0
            for j in np.arange(len(spectra_window)):
                convolved_value += spectra_window['flux'][j] * gaussian[j]
            
            convolved_spectra['flux'][i] = convolved_value
            #~ print i, spectra['flux'][i], "\t", convolved_value
        else:
            print "Not enough points for", lambda_peak, "(window = ", np.round(max - min, 2), "< 3*fwhm =", np.round(3*fwhm, 2), ")"
            convolved_spectra['flux'][i] = None
        
        if (i % 1000 == 0):
            if frame != None:
                current_work_progress = (i*1.0 / total_points) * 100
                frame.update_progress(current_work_progress)
            else:
                print "%.2f" % spectra['waveobs'][i]
        
    return convolved_spectra


###################################################################################################
################### Convertion of Caroline spectras of test stars to UVES resolution
if __name__ == '__main__':
    #### MP Giant (Arcturus)
    print "MP Giant (Arcturus)"
    atlas = pyfits.open('/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Caroline/ardata.fits')
    total_points = len(atlas[1].data['WAVELENGTH'])
    arcturus_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    arcturus_spectra['flux'] = atlas[1].data['ARCTURUS']
    arcturus_spectra['waveobs'] = atlas[1].data['WAVELENGTH'] / 10 # nm
    atlas.close()

    # Use only UVES range but wider to ensure convolution
    wave_filter1 = (arcturus_spectra['waveobs'] >= 470) & (arcturus_spectra['waveobs'] <= 690)
    # Detector error are marked with negative fluxes (atlas)
    wave_filter2 = (arcturus_spectra['flux'] > 0)
    spectra = arcturus_spectra[wave_filter1 & wave_filter2]
    spectra.sort(order='waveobs') # Make sure it is ordered by wavelength

    # Convert
    from_resolution = 150000
    to_resolution = 47000
    convolved_spectra = convolve_spectra(spectra, from_resolution, to_resolution)

    # Filter out values not convolved because of not having enough points for the window creation (spectra's edges)
    wave_filter1 = (np.logical_not(np.isnan(convolved_spectra['flux'])))
    # Use only UVES range
    wave_filter2 = (convolved_spectra['waveobs'] >= 480) & (convolved_spectra['waveobs'] <= 680)

    # Save
    write_spectra(convolved_spectra[wave_filter1 & wave_filter2], "output/LUMBA/UVES_arcturus.s.gz", compress=True)

    #### MR Dwarf (Sun)
    print "MR Dwarf (Sun)"
    atlas = pyfits.open('/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Caroline/ardata.fits')
    total_points = len(atlas[1].data['WAVELENGTH'])
    sun_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    sun_spectra['flux'] = atlas[1].data['SOLARFLUX']
    sun_spectra['waveobs'] = atlas[1].data['WAVELENGTH'] / 10 # nm
    atlas.close()

    # Use only UVES range but wider to ensure convolution
    wave_filter1 = (sun_spectra['waveobs'] >= 470) & (sun_spectra['waveobs'] <= 690)
    # Detector error are marked with negative fluxes (atlas)
    wave_filter2 = (sun_spectra['flux'] > 0)
    spectra = sun_spectra[wave_filter1 & wave_filter2]
    spectra.sort(order='waveobs') # Make sure it is ordered by wavelength

    # Convert
    from_resolution = 150000
    to_resolution = 47000
    convolved_spectra = convolve_spectra(spectra, from_resolution, to_resolution)

    # Filter out values not convolved because of not having enough points for the window creation (spectra's edges)
    wave_filter1 = (np.logical_not(np.isnan(convolved_spectra['flux'])))
    # Use only UVES range
    wave_filter2 = (convolved_spectra['waveobs'] >= 480) & (convolved_spectra['waveobs'] <= 680)

    # Save
    write_spectra(convolved_spectra[wave_filter1 & wave_filter2], "output/LUMBA/UVES_sun.s.gz", compress=True)

    #### MR Giant (Mu Leo)
    print "MR Giant (Mu Leo)"
    muleo = pyfits.open('/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Caroline/ESPADONS_hd85503_004.fits')

    total_points = len(muleo[0].data)
    muleo_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    muleo_spectra['flux'] = muleo[0].data

    base_wave = muleo[0].header["CRVAL1"] / 10 # nm
    increment = muleo[0].header["CDELT1"] / 10 # nm
    top_wave = base_wave + total_points * increment
    muleo_spectra['waveobs'] = np.linspace(base_wave, top_wave, total_points)

    muleo.close()

    # Use only UVES range but wider to ensure convolution
    wave_filter = (muleo_spectra['waveobs'] >= 470) & (muleo_spectra['waveobs'] <= 690)
    spectra = muleo_spectra[wave_filter]
    spectra.sort(order='waveobs') # Make sure it is ordered by wavelength

    ## In this case, radial vel corrections should be applied
    #radial_vel = 14.03 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+leo
    #spectra = correct_radial_velocity(spectra, radial_vel)

    # Convert
    from_resolution = 68000
    to_resolution = 47000
    convolved_spectra = convolve_spectra(spectra, from_resolution, to_resolution)

    # Filter out values not convolved because of not having enough points for the window creation (spectra's edges)
    wave_filter1 = (np.logical_not(np.isnan(convolved_spectra['flux'])))
    # Use only UVES range
    wave_filter2 = (convolved_spectra['waveobs'] >= 480) & (convolved_spectra['waveobs'] <= 680)

    # Save
    write_spectra(convolved_spectra[wave_filter1 & wave_filter2], "output/LUMBA/UVES_mu_leo.s.gz", compress=True)


    #### MP Dwarf (Mu Cas A)
    print "MP Dwarf (Mu Cas A)"
    mucas_spectra = read_spectra("/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Caroline/mu.cas_narval_26nov09_sp2_Normal_I_001_Normalized.s")

    # Use only UVES range but wider to ensure convolution
    wave_filter = (mucas_spectra['waveobs'] >= 470) & (mucas_spectra['waveobs'] <= 690)
    spectra = mucas_spectra[wave_filter]
    spectra.sort(order='waveobs') # Make sure it is ordered by wavelength

    ## In this case, radial vel corrections should be applied
    #radial_vel = -97.2 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+cas+a
    #spectra = correct_radial_velocity(spectra, radial_vel)

    # Convert
    from_resolution = 65000
    to_resolution = 47000
    convolved_spectra = convolve_spectra(spectra, from_resolution, to_resolution)

    # Filter out values not convolved because of not having enough points for the window creation (spectra's edges)
    wave_filter1 = (np.logical_not(np.isnan(convolved_spectra['flux'])))
    # Use only UVES range
    wave_filter2 = (convolved_spectra['waveobs'] >= 480) & (convolved_spectra['waveobs'] <= 680)

    # Save
    write_spectra(convolved_spectra[wave_filter1 & wave_filter2], "output/LUMBA/UVES_mu_cas_a.s.gz", compress=True)


    ###### Tests
    arcturus_spectra = read_spectra("output/LUMBA/UVES_arcturus.s.gz")
    sun_spectra = read_spectra("output/LUMBA/UVES_sun.s.gz")
    muleo_spectra = read_spectra("output/LUMBA/UVES_mu_leo.s.gz")
    mucas_spectra = read_spectra("output/LUMBA/UVES_mu_cas_a.s.gz")
    #~ plot_spectra([arcturus_spectra, sun_spectra, muleo_spectra, mucas_spectra])

    filename = "/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Official/spectra/MPG_UVESGauss06.06.dat"
    spectra = asciitable.read(table=filename, delimiter=' ', names=['Wavelength', 'Flux'])
    total_points = len(spectra['Wavelength'])
    arcturus_official = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    arcturus_official['waveobs'] = spectra['Wavelength'] / 10 # nm
    arcturus_official['flux'] = spectra['Flux']
    arcturus_official['err'] = 0.0

    filename = "/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Official/spectra/MRD_UVESGauss06.06.dat"
    spectra = asciitable.read(table=filename, delimiter=' ', names=['Wavelength', 'Flux'])
    total_points = len(spectra['Wavelength'])
    sun_official = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    sun_official['waveobs'] = spectra['Wavelength'] / 10 # nm
    sun_official['flux'] = spectra['Flux']
    sun_official['err'] = 0.0

    filename = "/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Official/spectra/MRG_UVESGauss05.16.dat"
    spectra = asciitable.read(table=filename, delimiter=' ', names=['Wavelength', 'Flux'])
    total_points = len(spectra['Wavelength'])
    muleo_official = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    muleo_official['waveobs'] = spectra['Wavelength'] / 10 # nm
    muleo_official['flux'] = spectra['Flux']
    muleo_official['err'] = 0.0

    filename = "/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Official/spectra/MPD_UVESGauss05.16.dat"
    spectra = asciitable.read(table=filename, delimiter=' ', names=['Wavelength', 'Flux'])
    total_points = len(spectra['Wavelength'])
    mucas_official = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    mucas_official['waveobs'] = spectra['Wavelength'] / 10 # nm
    mucas_official['flux'] = spectra['Flux']
    mucas_official['err'] = 0.0

    write_spectra(arcturus_official, "output/LUMBA/UVES_arcturus_official.s.gz", compress=True)
    write_spectra(sun_official, "output/LUMBA/UVES_sun_official.s.gz", compress=True)
    write_spectra(muleo_official, "output/LUMBA/UVES_mu_leo_official.s.gz", compress=True)
    write_spectra(mucas_official, "output/LUMBA/UVES_mu_cas_a_official.s.gz", compress=True)
    
    #~ plot_spectra([arcturus_spectra, arcturus_official])
    #~ plot_spectra([sun_spectra, sun_official])
    #~ plot_spectra([muleo_spectra, muleo_official])
    #~ plot_spectra([mucas_spectra, mucas_official])
    #~ plot_spectra([arcturus_official, sun_official, muleo_official, mucas_official])
