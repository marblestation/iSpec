
import asciitable
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt

# Get flux for a given spectra and wavelength, interpolating if needed
# For interpolation it considers:
# - 4 points in general
# - 2 when there are not more (i.e. at the beginning of the array or outside)
def get_flux(spectra, wavelength):
    # Objective wavelength
    objective_wavelength = wavelength
    
    # Find the index position of the first wave length equal or higher than the objective
#    index = np.where(spectra['waveobs'] >= objective_wavelength)[0][0]
    index = spectra['waveobs'].searchsorted(objective_wavelength)
    
    total_points = len(spectra)
    if index == total_points:
        # DISCARD: Linear extrapolation using index-1 and index-2
        # flux = spectra['flux'][index-1] + (objective_wavelength - spectra['waveobs'][index-1]) * ((spectra['flux'][index-1]-spectra['flux'][index-2])/(spectra['waveobs'][index-1]-spectra['waveobs'][index-2]))
        # JUST DUPLICATE:
        flux = spectra['flux'][index-1]
    elif index == 1 or index == total_points-1:
        # Linear interpolation between index and index-1
        # http://en.wikipedia.org/wiki/Linear_interpolation#Linear_interpolation_between_two_known_points
        flux = spectra['flux'][index-1] + (objective_wavelength - spectra['waveobs'][index-1]) * ((spectra['flux'][index]-spectra['flux'][index-1])/(spectra['waveobs'][index]-spectra['waveobs'][index-1]))
    elif index == 0 and spectra['waveobs'][index] != objective_wavelength:
        # DISCARD: Linear extrapolation using index+1 and index
        # flux = spectra['flux'][index] + (objective_wavelength - spectra['waveobs'][index]) * ((spectra['flux'][index+1]-spectra['flux'][index])/(spectra['waveobs'][index+1]-spectra['waveobs'][index]))
        # JUST DUPLICATE:
        flux = spectra['flux'][index]
    elif spectra['waveobs'][index] == objective_wavelength:
        flux = spectra['flux'][index]
    else:
        # Bessel's Central-Difference Interpolation with 4 points
        #   p = [(x - x0) / (x1 - x0)]
        #   f(x) = f(x0) + p ( f(x1) - f(x0) ) + [ p ( p - 1 ) / 4 ] ( f(x2) - f(x1) - f(x0) + f(x-1) )
        # where x-1 < x0 < objective_wavelength = x < x1 < x2 and f() is the flux
        #   http://physics.gmu.edu/~amin/phys251/Topics/NumAnalysis/Approximation/polynomialInterp.html
        
        #  x-1= index - 2
        #  x0 = index - 1
        #  x  = objective_wavelength 
        #  x1 = index
        #  x2 = index + 1
        
        ## Array access optimization
        flux_x_1 = spectra['flux'][index - 2]
        wave_x0 = spectra['waveobs'][index-1]
        flux_x0 = spectra['flux'][index - 1]
        wave_x1 = spectra['waveobs'][index]
        flux_x1 = spectra['flux'][index]
        flux_x2 = spectra['flux'][index + 1]
        
        p = (objective_wavelength - wave_x0) / (wave_x1 - wave_x0)
        flux = flux_x0 + p * (flux_x1 - flux_x0) + (p * (p - 1) / 4) * (flux_x2 - flux_x1 - flux_x0 + flux_x_1)
        
        
#    print flux, spectra['flux'][index], wavelength
    return flux, index


# Returns a grid of wavelength from base to top and with the minimum number of
# points needed to hold a spectra of a given resolution
def generate_wavelength_grid(base_wave, top_wave, resolution, points_per_fwhm = 3):
    xaxis = []
    current_wave = base_wave
    while current_wave <= top_wave:
        fwhm = current_wave/resolution
        wave_step = fwhm / points_per_fwhm
        for j in np.arange(points_per_fwhm):
            xaxis.append(current_wave)
            current_wave += wave_step 
    xaxis = np.array(xaxis)
    
    return xaxis

# Returns a new spectra with measures at the given xaxis wavelength
# Interpolation is linear when there are not enough points and Bessel's 
# Central-Difference Interpolation with 4 points for the rest
def resample_spectra(spectra, xaxis):
    total_real_points = len(spectra)
    total_points = len(xaxis)
    
    resampled_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    resampled_spectra['waveobs'] = xaxis
     
    from_index = 0 # Optimization: discard regions already processed
    for i in np.arange(total_points):
        resampled_spectra['flux'][i], index = get_flux(spectra[from_index:], resampled_spectra['waveobs'][i])
        resampled_spectra['err'][i] = -999.9 # TODO: Calculate errors (propagation or interpolation?)
        if index > 4:
            from_index = index - 4
        if (i % 1000 == 0):
            print "%.2f" % resampled_spectra['waveobs'][i]
    
    return resampled_spectra

