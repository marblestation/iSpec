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
    from_delta_lambda = lambda_peak / from_resolution
    to_delta_lambda = lambda_peak / to_resolution
    fwhm = np.sqrt(to_delta_lambda**2 - from_delta_lambda**2)
    return fwhm

# Convert the FWHM to the sigma for building a gaussian
def get_sigma(fwhm):
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    return sigma

## Degradation of the spectra resolution. Procedure for each flux value:
#   1. Define a window based on the fwhm size
#   2. Build a gaussian using the sigma value and the wavelength values of the spectra window
#   3. Convolve the spectra window with the gaussian and save the convolved value
def degrade_spectra_resolution(spectra, from_resolution, to_resolution):
    if from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are equal or bigger than original")

    total_points = len(spectra['waveobs'])
    convolved_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    convolved_spectra['waveobs'] = spectra['waveobs']
    
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
    convolved_spectra = degrade_spectra_resolution(spectra, from_resolution, to_resolution)

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
    convolved_spectra = degrade_spectra_resolution(spectra, from_resolution, to_resolution)

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
    convolved_spectra = degrade_spectra_resolution(spectra, from_resolution, to_resolution)

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
    convolved_spectra = degrade_spectra_resolution(spectra, from_resolution, to_resolution)

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
