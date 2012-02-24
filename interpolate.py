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
    fluxes = spectra['flux']
    waveobs = spectra['waveobs']
    
    # Find the index position of the first wave length equal or higher than the objective
#    index = np.where(waveobs >= objective_wavelength)[0][0]
    index = waveobs.searchsorted(objective_wavelength)
    
    total_points = len(spectra)
    if index == total_points:
        # DISCARD: Linear extrapolation using index-1 and index-2
        # flux = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index-1]-fluxes[index-2])/(waveobs[index-1]-waveobs[index-2]))
        # JUST DUPLICATE:
        flux = fluxes[index-1]
    elif index == 1 or index == total_points-1:
        # Linear interpolation between index and index-1
        # http://en.wikipedia.org/wiki/Linear_interpolation#Linear_interpolation_between_two_known_points
        flux = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index]-fluxes[index-1])/(waveobs[index]-waveobs[index-1]))
    elif index == 0 and waveobs[index] != objective_wavelength:
        # DISCARD: Linear extrapolation using index+1 and index
        # flux = fluxes[index] + (objective_wavelength - waveobs[index]) * ((fluxes[index+1]-fluxes[index])/(waveobs[index+1]-waveobs[index]))
        # JUST DUPLICATE:
        flux = fluxes[index]
    elif waveobs[index] == objective_wavelength:
        flux = fluxes[index]
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
        flux_x_1 = fluxes[index - 2]
        wave_x0 = waveobs[index-1]
        flux_x0 = fluxes[index - 1]
        wave_x1 = waveobs[index]
        flux_x1 = fluxes[index]
        flux_x2 = fluxes[index + 1]
        
        p = (objective_wavelength - wave_x0) / (wave_x1 - wave_x0)
        flux = flux_x0 + p * (flux_x1 - flux_x0) + (p * (p - 1) / 4) * (flux_x2 - flux_x1 - flux_x0 + flux_x_1)
        
        
#    print flux, fluxes[index], wavelength
    return flux, index


# Returns a grid of wavelength from base to top and with the minimum number of
# points needed to hold a spectra of a given resolution
def generate_wavelength_grid(base_wave, top_wave, resolution, points_per_fwhm = 3):
    xaxis = []
    current_wave = base_wave * 1.0 # Ensure it is a float
    top_wave = top_wave * 1.0
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

