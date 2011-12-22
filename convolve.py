import numpy as np
import pyfits
import scipy.ndimage as ndi
from interpolate import *
from plotting import *

def get_fwhm(lambda_peak, resolution):
    delta_lambda = {}
    delta_lambda['original'] = lambda_peak / resolution['original']
    delta_lambda['final'] = lambda_peak / resolution['final']
    fwhm = np.sqrt(delta_lambda['final']**2 - delta_lambda['original']**2)
    return fwhm

def get_sigma(fwhm):
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    return sigma

atlas = pyfits.open('/mnt/extra/Data/Universidad/Bordeaux1/LUMBA/Test stars/Caroline/ardata.fits')
total_points = len(atlas[1].data['WAVELENGTH'])
arcturus_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
arcturus_spectra['flux'] = atlas[1].data['ARCTURUS']
arcturus_spectra['waveobs'] = atlas[1].data['WAVELENGTH']
atlas.close()

spectra = arcturus_spectra[:1000]
spectra.sort(order='waveobs') # Make sure it is ordered by wavelength

resolution = {}
resolution['original'] = 150000
resolution['final'] = 47000

if resolution['original'] <= resolution['final']:
    raise Exception("This method degrades resolution, it cannot deal with a different situation where final resolution is equal or bigger than original")

total_points = len(spectra['waveobs'])
convolved_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
convolved_spectra['waveobs'] = spectra['waveobs']

#~ PROCEDURE:
#~  1. Interpolation of the flux on a fine spaced grid (Oversampling factor: 17)
#~  2. Convolution of the flux with a gaussian profile 
#~  3. Interpolation of the smoothed flux on the original grid.

for i in np.arange(total_points):
    lambda_peak = spectra['waveobs'][i] # Current lambda (wavelength) to be modified
    fwhm = get_fwhm(lambda_peak, resolution) # Calculate the needed fwhm for this wavelength
    
    # Only work with a limited window considering 2 times the fwhm in each side of the current
    # position to be modified and saved in the convolved spectra
    wave_filter = (spectra['waveobs'] >= lambda_peak - 2*fwhm) & (spectra['waveobs'] <= lambda_peak + 2*fwhm)
    spectra_window = spectra[wave_filter]
    min = np.min(spectra_window['waveobs'])
    max = np.max(spectra_window['waveobs'])
    # Find the index position of the current lambda peak in the window
    lambda_peak_win_index = spectra_window['waveobs'].searchsorted(lambda_peak)
    
    # Check that the window is large enough (i.e. points in the edges of the original spectra)
    if np.round(max - min, 2) >= 3*fwhm:
        sigma = get_sigma(fwhm)
        
        # Construct the gaussian
        gaussian = np.exp(- ((spectra_window['waveobs'] - lambda_peak)**2) / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
        gaussian = gaussian / np.sum(gaussian)
        convolved_value = 0
        for j in np.arange(len(spectra_window)):
            convolved_value += spectra_window['flux'][j] * gaussian[j]
        
        convolved_spectra['flux'][i] = convolved_value
        print i, spectra['flux'][i], "\t", convolved_value
        #~ break
        #~ # Interpolate on a fine spaced grid
        #~ wavelength_step = fwhm / 17.0
        #~ xaxis = np.arange(min, max, wavelength_step)
        #~ resampled_spectra_window = resample_spectra2axis(spectra_window, xaxis)
        #~ 
        #~ # Gauss filter (convolution) for the window (which is centered in the current lambda_peak)
        #~ convolved_spectra_window = np.recarray((len(xaxis), ), dtype=[('waveobs', float),('flux', float),('err', float)])
        #~ convolved_spectra_window['waveobs'] = resampled_spectra_window['waveobs']
        #~ convolved_spectra_window['flux'] = ndi.gaussian_filter1d(resampled_spectra_window['flux'], sigma)
        #~ 
        #~ # Interpolate on the original wavelength grid (which can be not spaced uniformly)
        #~ resampled_convolved_spectra_window = resample_spectra2axis(convolved_spectra_window, spectra_window['waveobs'])
        #~ 
        #~ # Save only the current lambda convolved value
        #~ convolved_spectra['flux'][i] = resampled_convolved_spectra_window[lambda_peak_win_index]['flux']
        #~ print spectra['flux'][i], "\t", convolved_spectra['flux'][i]
    else:
        print "Not enough points for", lambda_peak, "(fwhm=", fwhm, ")"
        convolved_spectra['flux'][i] = 0


#~ plot_spectra([spectra, convolved_spectra])

