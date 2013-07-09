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
def create_mask(np.ndarray[np.double_t,ndim=1] spectrum_wave, np.ndarray[np.double_t,ndim=1] mask_wave, np.ndarray[np.double_t,ndim=1] mask_values, velocity_mask_size=2.0):
    """
    It constructs a zero flux spectrum and assign mask values to the wavelengths
    belonging to that value and its surounds (determined by the velocity_mask_size).
    """

    ## Speed of light in m/s
    cdef double c = 299792458.0
    
    cdef np.ndarray[np.double_t,ndim=1] mask_wave_step
    cdef np.ndarray[np.double_t,ndim=1] mask_wave_base
    cdef np.ndarray[np.double_t,ndim=1] mask_wave_top
    cdef np.ndarray[np.double_t,ndim=1] resampled_mask = np.zeros(len(spectrum_wave))

    # Mask limits
    mask_wave_step = (mask_wave * (1.-np.sqrt((1.-(velocity_mask_size*1000.)/c)/(1.+(velocity_mask_size*1000.)/c))))/2.0
    mask_wave_base = mask_wave - 1*mask_wave_step
    mask_wave_top = mask_wave + 1*mask_wave_step
    
    cdef double wstep = spectrum_wave[1] - spectrum_wave[0]
    cdef double wbase = spectrum_wave[0]

    cdef int i = 0
    cdef int j = 0
    for i from 0 <= i < len(mask_wave) by 1:
        #j = 0
        while j < len(spectrum_wave) and spectrum_wave[j] < mask_wave_base[i]:
            j += 1
        while j < len(spectrum_wave) and spectrum_wave[j] >= mask_wave_base[i] and spectrum_wave[j] <= mask_wave_top[i]:
            resampled_mask[j] = mask_values[i]
            j += 1

    return resampled_mask

