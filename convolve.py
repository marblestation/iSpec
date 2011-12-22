import numpy as np
import pyfits
import scipy.ndimage as ndi
from interpolate import *

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

spectra = arcturus_spectra[:100]
spectra.sort(order='waveobs') # Make sure it is ordered by wavelength


## TODO: Select window margin in function of FWHM
window_margin = 20 # Total window width: window_margin + 1 + window_margin
resolution = {}
resolution['original'] = 150000
resolution['final'] = 47000

lambda_max = np.max(spectra['waveobs'])
fwhm_max = get_fwhm(lambda_max, resolution)

lambda_peak = 400.0
wave_filter = (spectra['waveobs'] >= lambda_peak - 1.5*fwhm_max) & (spectra['waveobs'] <= lambda_peak_peak + 1.5*fwhm_max)    


if resolution['original'] <= resolution['final']:
    # This method degrades resolution, it cannot deal with a different situation
    exit()


total_points = len(spectra['waveobs']) - 2 * window_margin

convolved_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
convolved_spectra['waveobs'] = spectra['waveobs'][window_margin:window_margin+total_points]

#~ PROCEDURE:
#~ ;           1. Interpolation of the flux on a fine spaced grid.
#~ ;              Oversamplingfaktor: 17
#~ ;           2. Convolution of the flux with a gaussian profile 
#~ ;           3. Interpolation of the smoothed flux on the original grid.

for i in np.arange(total_points):
    if i >= window_margin and i < total_points - window_margin:
        j = i - window_margin
        
        lambda_peak = spectra['waveobs'][i]
        fwhm = get_fwhm(lambda_peak, resolution)
        sigma = get_sigma(fwhm)
        
        spectra_window = spectra[i-window_margin:i+1+window_margin]
        min = np.min(spectra_window['waveobs'])
        max = np.max(spectra_window['waveobs'])
        wavelength_step = fwhm / 17.0
        xaxis = np.arange(min, max, wavelength_step)
        resampled_spectra_window = resample_spectra2axis(spectra_window, xaxis)
        
        # Gauss filter (convolution) for the window centered in the current lambda_peak
        convoluted_spectra_window = np.recarray((len(xaxis), ), dtype=[('waveobs', float),('flux', float),('err', float)])
        convoluted_spectra_window['waveobs'] = resampled_spectra_window['waveobs']
        convoluted_spectra_window['flux'] = ndi.gaussian_filter1d(resampled_spectra_window['flux'], sigma)
        
        #~ print convolved_spectra['flux'][j], "\t", convoluted_spectra_window[window_margin], "\t", sigma, 3*fwhm
        print 3*fwhm, np.abs(spectra['waveobs'][i-window_margin] - spectra['waveobs'][i+2+window_margin]), np.abs(spectra['waveobs'][i-window_margin] - spectra['waveobs'][i+2+window_margin])/fwhm
        
        resampled_convoluted_spectra_window = resample_spectra2axis(convoluted_spectra_window, spectra_window['waveobs'])
        
        # Save only the peak convolved value
        convolved_spectra['flux'][j] = resampled_convoluted_spectra_window[window_margin]['flux']
        #~ print spectra['flux'][i], "\t", convolved_spectra['flux'][j]



