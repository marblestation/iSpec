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
cimport numpy as np

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)
def cross_correlate(np.ndarray[np.double_t,ndim=1] mask_wave, np.ndarray[np.double_t,ndim=1] mask_values, np.ndarray[np.double_t,ndim=1] spectrum_wave, np.ndarray[np.double_t,ndim=1] spectrum_flux, np.ndarray[np.double_t,ndim=1] spectrum_err, double velocity):

    ## Result
    cdef double ccf = 0.0
    cdef double ccf_err = 0.0
    
    ## Speed of light in m/s
    cdef double c = 299792458.0
    
    ## Consider masks of X km/s around the peak
    cdef double velocity_mask_size = 2.0
    cdef int points_per_mask = 10
    
    cdef np.ndarray[np.double_t,ndim=1] shifted_mask_wave
    cdef np.ndarray[np.double_t,ndim=1] shifted_mask_wave_step
    cdef np.ndarray[np.double_t,ndim=1] shifted_mask_wave_base
    cdef np.ndarray[np.double_t,ndim=1] shifted_mask_wave_top

    ## Correct mask with the appropiate velocity
    # Relativistic correction
    shifted_mask_wave = mask_wave / np.sqrt((1.-(velocity*1000.)/c)/(1.+(velocity*1000.)/c))
    
    # Mask limits
    shifted_mask_wave_step = (shifted_mask_wave * (1.-np.sqrt((1.-(velocity_mask_size*1000.)/c)/(1.+(velocity_mask_size*1000.)/c))))/2.0
    # NOTE: I don't understand why the mask that works is the second one (with the first one there is a systematic shift of around 1 km/s):
    #   Lines in mask:            x          x          x
    #   Mask that should work:  |   |      |   |      |   |
    #   Mask that works:          |   |      |   |      |   |
    shifted_mask_wave_base = shifted_mask_wave - 0*shifted_mask_wave_step
    shifted_mask_wave_top = shifted_mask_wave + 2*shifted_mask_wave_step
    
    cdef double wstep = spectrum_wave[1] - spectrum_wave[0]
    cdef double wbase = spectrum_wave[0]

    cdef int i = 0
    cdef int base_index = 0
    cdef int top_index = 0
    for i in xrange(len(mask_wave)):
        # Find indexes for the spectrum
        base_index = int_min(int_max(int((shifted_mask_wave_base[i] - wbase)/wstep), 0), len(spectrum_wave))
        top_index = int_min(int_max(int((shifted_mask_wave_top[i] - wbase)/wstep), 0), len(spectrum_wave))
        if base_index >= len(spectrum_wave) or top_index >= len(spectrum_wave) or base_index >= top_index:
            continue
        # Make sure we consider always 100 points per mask
        # to avoid systematic effects
        wrange = np.arange(points_per_mask)*((spectrum_flux[top_index]-spectrum_flux[base_index])/points_per_mask) + spectrum_flux[base_index]
        fluxes = np.interp(wrange, spectrum_wave[base_index:top_index], spectrum_flux[base_index:top_index])
        errors = np.interp(wrange, spectrum_wave[base_index:top_index], spectrum_err[base_index:top_index])
        ccf += np.sum(mask_values[i] * fluxes)
        ccf_err += np.sum(mask_values[i] * errors)

    return ccf, ccf_err

