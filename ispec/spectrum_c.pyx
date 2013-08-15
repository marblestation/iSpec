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
import logging
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)
    double log(double x)

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
#cdef inline int float_max(float a, float b): return a if a >= b else b
#cdef inline int float_min(float a, float b): return a if a <= b else b

from cpython cimport bool

cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)
def convolve_spectrum(np.ndarray[np.double_t,ndim=1] waveobs, np.ndarray[np.double_t,ndim=1] flux, np.ndarray[np.double_t,ndim=1] err, double to_resolution, from_resolution=None, frame=None):
    if from_resolution != None and from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")
    cdef int total_points = len(waveobs)
    cdef np.ndarray[np.double_t,ndim=1] convolved_flux = np.zeros(total_points)
    cdef np.ndarray[np.double_t,ndim=1] convolved_err = np.zeros(total_points)
    cdef double fwhm
    cdef double sigma

    cdef int i

    cdef double current_work_progress
    cdef double last_reported_progress = -1
    if frame != None:
        frame.update_progress(0)

    cdef double lambda_peak
    cdef int wave_filter_base
    cdef int wave_filter_top
    cdef np.ndarray[np.double_t,ndim=1] waveobs_window
    cdef np.ndarray[np.double_t,ndim=1] flux_window
    cdef np.ndarray[np.double_t,ndim=1] err_window
    cdef np.ndarray[np.double_t,ndim=1] gaussian
    cdef double total_gaussian
    cdef int last_x = 0
    cdef int x
    cdef int y
    for i in np.arange(total_points):
    #for i from 0 <= i < total_points by 1:
        if flux[i] <= 1e-10:
            continue

        # FWHM of the gaussian for the given resolution
        if from_resolution == None or from_resolution == 0:
            # Convolve using instrumental resolution (smooth but not degrade)
            fwhm = waveobs[i] / to_resolution
            sigma = fwhm / 2.3548200450309493 #(2*sqrt(2*log(2)))
        else:
            # Degrade resolution
            fwhm = sqrt(((1.0*waveobs[i]) / to_resolution)**2 - ((1.0*waveobs[i]) / from_resolution)**2)
            sigma = fwhm / 2.3548200450309493 #(2*sqrt(2*log(2)))

        lambda_peak = waveobs[i] # Current lambda (wavelength) to be modified

        # Only work with a limited window considering 2 times the fwhm in each side of the current
        # position to be modified and saved in the convolved spectra
        wave_filter_base = 0
        for x from last_x <= x < len(waveobs) by 1:
            if waveobs[x] >= lambda_peak - 2*fwhm:
                wave_filter_base = x
                break
        last_x = x
        for y from last_x <= y < len(waveobs) by 1:
            if waveobs[y] >= lambda_peak + 2*fwhm:
                wave_filter_top = y
                break
        waveobs_window = waveobs[wave_filter_base:wave_filter_top]
        flux_window = flux[wave_filter_base:wave_filter_top]
        if err[i] > 0:
            # * Propagate error Only if the current value has a valid error value assigned
            #
            # Error propagation considering that measures are dependent (more conservative approach)
            # because it is common to find spectra with errors calculated from a SNR which
            # at the same time has been estimated from all the measurements in the same spectra
            #
            err_window = err[wave_filter_base:wave_filter_top]

        # Construct the gaussian
        gaussian = np.zeros(len(flux_window))
        total_gaussian = 0.0
        for x from 0 <= x < len(flux_window) by 1:
            gaussian[x] = exp(- ((waveobs_window[x] - lambda_peak)**2) / (2*sigma**2)) / sqrt(2*np.pi*sigma)
            if flux_window[x] > 0.0: # Zero or negative values are considered as gaps in the spectrum
                convolved_flux[i] += flux_window[x] * gaussian[x]
                total_gaussian += gaussian[x]
                if err[i] > 0:
                    convolved_err[i] += err_window[x] * gaussian[x]
                    #convolved_err[i] += (err_window[x] * gaussian[x]) * (err_window[x] * gaussian[x]) # Independent measurements
        convolved_flux[i] /= total_gaussian
        #convolved_err[i] = sqrt(convolved_err[i]) # Independent measurements
        convolved_err[i] /= total_gaussian

        current_work_progress = (i*1.0 / total_points) * 100
        if (int(current_work_progress) % 10 == 0 and current_work_progress - last_reported_progress > 10) or last_reported_progress < 0 or current_work_progress == 100:
            last_reported_progress = current_work_progress
            #logging.info("%.2f%%" % current_work_progress)
            if frame != None:
                frame.update_progress(current_work_progress)

    return waveobs, convolved_flux, convolved_err


@cython.boundscheck(False)
@cython.wraparound(False)
def bessel_interpolation(np.ndarray[np.double_t,ndim=1] waveobs, np.ndarray[np.double_t,ndim=1] fluxes, np.ndarray[np.double_t,ndim=1] err, np.ndarray[np.double_t,ndim=1] resampled_waveobs, frame=None):
    """
    Interpolate flux for a given wavelength by using Bessel's Central-Difference Interpolation.
    It considers:

    - 4 points in general
    - 2 when there are not more (i.e. at the beginning of the array or outside)

    * It does not interpolate if any of the fluxes used for interpolation is zero or negative
      this way it can respect gaps in the spectrum
    """
    cdef double current_work_progress
    cdef double last_reported_progress = -1
    if frame != None:
        frame.update_progress(0)

    cdef int total_points = len(waveobs)
    cdef int new_total_points = len(resampled_waveobs)
    cdef np.ndarray[np.double_t,ndim=1] resampled_flux = np.zeros(new_total_points)
    cdef np.ndarray[np.double_t,ndim=1] resampled_err = np.zeros(new_total_points)
    cdef int index = 0
    cdef int last_x = 0
    cdef int x
    cdef int i
    cdef double flux_x_1
    cdef double wave_x0
    cdef double flux_x0
    cdef double wave_x1
    cdef double flux_x1
    cdef double flux_x2
    cdef double p

    for i from 0 <= i < new_total_points by 1:
        # Target wavelength
        objective_wavelength = resampled_waveobs[i]

        # Find the index position of the first wave length equal or higher than the objective
        # index = np.where(waveobs >= objective_wavelength)[0][0]
        index = 0
        for x from last_x <= x < total_points by 1:
            if waveobs[x] >= objective_wavelength:
                index = x
                break
        last_x = x

        if index == total_points:
            # DISCARD: Linear extrapolation using index-1 and index-2
            # flux = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index-1]-fluxes[index-2])/(waveobs[index-1]-waveobs[index-2]))
            # JUST DUPLICATE:
            #resampled_flux[i] = fluxes[index-1]
            # JUST ZERO:
            resampled_flux[i] = 0.0
            resampled_err[i] = 0.0
        elif index == 1 or index == total_points-1:
            # Do not interpolate if any of the fluxes is zero or negative
            if fluxes[index-1] <= 1e-10 or fluxes[index] <= 1e-10:
                resampled_flux[i] = 0.0
                resampled_err[i] = 0.0
            else:
                # Linear interpolation between index and index-1
                # http://en.wikipedia.org/wiki/Linear_interpolation#Linear_interpolation_between_two_known_points
                resampled_flux[i] = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index]-fluxes[index-1])/(waveobs[index]-waveobs[index-1]))
                resampled_err[i] = err[index-1] + (objective_wavelength - waveobs[index-1]) * ((err[index]-err[index-1])/(waveobs[index]-waveobs[index-1]))
        elif index == 0 and waveobs[index] != objective_wavelength:
            # DISCARD: Linear extrapolation using index+1 and index
            # flux = fluxes[index] + (objective_wavelength - waveobs[index]) * ((fluxes[index+1]-fluxes[index])/(waveobs[index+1]-waveobs[index]))
            # JUST DUPLICATE:
            #resampled_flux[i] = fluxes[index]
            # JUST ZERO:
            resampled_flux[i] = 0.0
            resampled_err[i] = 0.0
        # Do not do this optimization because it can produce a value surounded
        # by zeros because of the condition "Do not interpolate if any of the
        # fluxes is zero or negative" implemented in the rest of the cases
        #elif waveobs[index] == objective_wavelength:
            #resampled_flux[i] = fluxes[index]
            #resampled_err[i] = err[index]
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

            err_x_1 = err[index - 2]
            err_x0 = err[index - 1]
            err_x1 = err[index]
            err_x2 = err[index + 1]

            # Do not interpolate if any of the fluxes is zero or negative
            if flux_x_1 <= 1e-10 or flux_x0 <= 1e-10 or flux_x1 <= 1e-10 or flux_x2 <= 1e-10:
                resampled_flux[i] = 0.0
                resampled_err[i] = 0.0
            else:
                p = (objective_wavelength - wave_x0) / (wave_x1 - wave_x0)
                resampled_flux[i] = flux_x0 + p * (flux_x1 - flux_x0) + (p * (p - 1) / 4) * (flux_x2 - flux_x1 - flux_x0 + flux_x_1)
                resampled_err[i] = err_x0 + p * (err_x1 - err_x0) + (p * (p - 1) / 4) * (err_x2 - err_x1 - err_x0 + err_x_1)

        current_work_progress = (i*1.0 / new_total_points) * 100
        if (int(current_work_progress) % 10 == 0 and current_work_progress - last_reported_progress > 10) or last_reported_progress < 0 or current_work_progress == 100:
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)
            if frame != None:
                frame.update_progress(current_work_progress)


    return resampled_waveobs, resampled_flux, resampled_err


