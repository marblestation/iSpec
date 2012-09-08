"""
    This file is part of Spectra Visual Editor (SVE).
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com

    SVE is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SVE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with SVE. If not, see <http://www.gnu.org/licenses/>.
"""
#!/usr/bin/env python
#import ipdb
import asciitable
import numpy as np
from plotting import *
from common import *
from convolve import *
from pymodelfit import UniformCDFKnotSplineModel

# Group the points in ranges of 1 nm (by default) and select the one with the maximum flux
def find_max_value_per_wavelength_range(spectra, base_points, wave_range=1):
    return find_a_value_per_wavelength_range(spectra, base_points, wave_range=wave_range, median=False)

# Group the points in ranges of 1 nm (by default) and select the one with the median flux
def find_median_value_per_wavelength_range(spectra, base_points, wave_range=1):
    return find_a_value_per_wavelength_range(spectra, base_points, wave_range=wave_range, median=True)

# Group the points in ranges of 1 nm (by default) and select the one with the maximum (by default) or median flux
def find_a_value_per_wavelength_range(spectra, base_points, wave_range=1, median=False):
    waveobs = spectra['waveobs']
    flux = spectra['flux']
    wave_step = wave_range # nm
    wave_base = np.min(waveobs)
    wave_top = np.max(waveobs)

    # Group points in bins and use only the one with the higher flux
    num_candidate_base_points = int((wave_top-wave_base)/wave_step) + 1
    candidate_base_points = -9999 * np.ones(num_candidate_base_points, dtype=int)
    i = 0
    while wave_base < wave_top:
        positions = np.where((waveobs[base_points] >= wave_base) & (waveobs[base_points] <= wave_base + wave_step))[0]

        if len(positions) == 0:
            candidate_base_points[i] = -9999
        else:
            if median:
                # Find the median (instead of the max)
                f = flux[base_points[positions]]
                sortedi = np.argsort(np.abs(f - np.median(f))) # The smallest is the median (this is the best way to avoid floating-point imprecisions)
                candidate_base_points[i] = base_points[positions[sortedi[0]]]
            else:
                # Find max for this bin
                sortedi = np.argsort(flux[base_points[positions]])
                candidate_base_points[i] = base_points[positions[sortedi[-1]]]

        #ipdb.set_trace()
        wave_base += wave_step
        i += 1

    candidate_base_points = candidate_base_points[candidate_base_points != -9999]
    return candidate_base_points


# Considering the diference in flux of the points with the next a previous, discard outliers (iterative process)
# - Median value and 3*sigma is used as criteria for filtering
def discard_outliers_for_continuum_candidates(spectra, candidate_base_points, sig=3):
    # The change between consecutive base points for continuum fitting should not be very strong,
    # identify outliers (first and last base point are excluded in this operation):
    flux_diff1 = (spectra['flux'][candidate_base_points][:-1] - spectra['flux'][candidate_base_points][1:]) / (spectra['waveobs'][candidate_base_points][:-1] - spectra['waveobs'][candidate_base_points][1:])
    flux_diff2 = (spectra['flux'][candidate_base_points][1:] - spectra['flux'][candidate_base_points][:-1]) / (spectra['waveobs'][candidate_base_points][1:] - spectra['waveobs'][candidate_base_points][:-1])
    # Recover first and last
    flux_diff1 = np.asarray(flux_diff1.tolist() + [flux_diff2[-1]])
    flux_diff2 = np.asarray([flux_diff1[0]] + flux_diff2.tolist())
    # Identify outliers
    flux_diff1_selected, not_outliers1 = sigma_clipping(flux_diff1, sig=sig, meanfunc=np.median)
    flux_diff2_selected, not_outliers2 = sigma_clipping(flux_diff2, sig=sig, meanfunc=np.median)
    outliers = np.logical_or(np.logical_not(not_outliers1), np.logical_not(not_outliers2))

    # Ensure that first and last points are not filtered out in order to avoid
    # having strange extrapolations in the edges of the spectrum because lack of points in the fit
    outliers = np.asarray([False] + outliers[1:-1].tolist() + [False])
    # Discard outliers
    continuum_base_points = candidate_base_points[~outliers]

    return continuum_base_points


# 1) Determine max points by using a moving window of 3 elements (also used for line determination)
# 2) Group the points:
#     - In ranges of 0.1 nm and select the one with the median flux (usefull to avoid noisy peaks)
#     - In ranges of 1 nm and select the one with the max flux (
# 3) Considering the diference in flux of the points with the next a previous, discard outliers (iterative process)
def determine_continuum_base_points(spectra, discard_outliers=True, median_wave_range=0.1, max_wave_range=1):
    # Find max points in windows of 3 measures
    candidate_base_points = find_local_max_values(spectra['flux'])
    if median_wave_range > 0:
        candidate_base_points = find_median_value_per_wavelength_range(spectra, candidate_base_points, wave_range=median_wave_range)
    if max_wave_range > 0:
        candidate_base_points = find_max_value_per_wavelength_range(spectra, candidate_base_points, wave_range=max_wave_range)
    if discard_outliers:
        candidate_base_points = discard_outliers_for_continuum_candidates(spectra, candidate_base_points)
    continuum_base_points = candidate_base_points


    return continuum_base_points


# 1) If no continuum points are specified for doing the fit, they are determined
# 2) Fit a knot spline model with n knots, if it is not specified, there will be 1 knot every 10 nm
#    - knots are located depending on the Cumulative Distribution Function, so there will be more
#      where more continuum points exist
# 3) Returns the fitted model
def fit_continuum(spectra, nknots=None, median_wave_range=0.1, max_wave_range=1):
    continuum_base_points = determine_continuum_base_points(spectra, discard_outliers=True, median_wave_range=median_wave_range, max_wave_range=max_wave_range)

    if nknots == None:
        # * 1 knot every 10 nm in average
        nknots = np.max([1, int((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 10)])

    if len(spectra['waveobs'][continuum_base_points]) == 0:
        raise Exception("Not enough points to fit")

    # UniformCDFKnotSplineModel:
    # 1) Builds an histogram: # points in 2*nknots bins
    #    - cdf: Number of points in that bin
    #    - xcdf[n], xcdf[n+1]: Limits of a bin in wavelength
    # 2) Filter bins with no points
    # 3) Calculate the cumulative sum of # points and normalize it
    # 4) Interpolate the wavelength in an homogeneus grid of normalized cumulative sum of points
    #    (from 0 to 1, divided in nknots)
    #    - x=0.02 means wavelength where we reach the 2% of the total number of points
    #    - x=0.98 means wavelength where we reach the 98% of the total number of points
    # 5) Use those interpolated wavelengths for putting the knots, therefore we will have a knot
    #    in those regions were there are an increment on the number of points (avoiding empty regions)
    continuum_model = UniformCDFKnotSplineModel(nknots)
    fitting_error = False
    try:
        continuum_model.fitData(spectra['waveobs'][continuum_base_points], spectra['flux'][continuum_base_points])
    except Exception, e:
        fitting_error = True

    # If there is no fit (because too few points)
    if fitting_error or np.any(np.isnan(continuum_model.residuals())):
        # Try without discarding outliers:
        continuum_base_points = determine_continuum_base_points(spectra, discard_outliers=False, median_wave_range=median_wave_range,
max_wave_range=max_wave_range)
        continuum_model.fitData(spectra['waveobs'][continuum_base_points], spectra['flux'][continuum_base_points])
        if np.any(np.isnan(continuum_model.residuals())):
            raise Exception("Not enough points to fit")

    return continuum_model

# Interpolate continuum points from a given continuum spectra
class MultiLinearInterpolationContinuumModel:
    def __init__(self, continuum_spectra):
        self.continuum_spectra = continuum_spectra

    def __call__(self, waveobs):
        flux = np.interp(waveobs, self.continuum_spectra['waveobs'], self.continuum_spectra['flux'])
        return flux


# 1) If no continuum points are specified for doing the interpolation, they are determined
# 2) Smooth (optional):
#    - If resolution is specified, build a interpolated continuum spectra and smooth it and
#      use those new points for the model
# 3) Build a continuum model that will interpolate values using the specified points
# 4) Returns the model
def interpolate_continuum(spectra, smooth=False, median_wave_range=0.1, max_wave_range=1):
    continuum_base_points = determine_continuum_base_points(spectra, discard_outliers=True, median_wave_range=median_wave_range, max_wave_range=max_wave_range)

    if len(spectra['waveobs'][continuum_base_points]) == 0:
        raise Exception("Not enough points to interpolate")

    if smooth:
        from scipy.ndimage import gaussian_filter1d

        # In order to smooth, first we reproduce a complete continuum spectrum
        # considering a uniform wavelength grid (separated by 1nm)
        waveobs = np.arange(np.min(spectra['waveobs']), np.max(spectra['waveobs']), 1)
        continuum_spectra = np.recarray((len(waveobs), ), dtype=[('waveobs', float),('flux', float),('err', float)])
        continuum_spectra['waveobs'] = waveobs
        # Linear interpolation
        continuum_spectra['flux'] = np.interp(continuum_spectra['waveobs'], spectra[continuum_base_points]['waveobs'], spectra[continuum_base_points]['flux'])
        continuum_spectra['err'] = 9999.9

        # Smooth (gaussian 1 sigma)
        continuum_spectra['flux'] = gaussian_filter1d(continuum_spectra['flux'], 1, order=0)
        return MultiLinearInterpolationContinuumModel(continuum_spectra)
    else:
        return MultiLinearInterpolationContinuumModel(spectra[continuum_base_points])

# Calculate flux for a wavelength grid using a given model
def get_spectra_from_model(model, spectra_wave_grid):
    total_points = len(spectra_wave_grid)
    spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    spectra['waveobs'] = spectra_wave_grid
    spectra['flux'] = model(spectra['waveobs'])
    spectra['err'] = np.zeros(total_points)
    return spectra


# Find regions of wavelengths where the fluxes seem to belong to the continuum.
# - It analyses the spectra in regions
# - The region size is variable in function of 4*fwhm which is derived
#   from the current wavelength and the resolution (unless a fixed_wave_step is specified)
# - For each region, if...
#     a) the median flux is above the continuum model (but not more than 0.08) or below but not more than 0.01
#         * The continuum model can be a fixed flux value or a fitted model (preferable)
#     b) and the standard deviation is less than a given maximum
#   the region is selected
def find_continuum(spectra, resolution, log_filename=None, max_std_continuum = 0.002, continuum_model = 0.95, max_continuum_diff=0.01, fixed_wave_step=None, frame=None):

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
        median_flux = np.median(spectra_region['flux'])
        std_flux = np.std(spectra_region['flux'])
        num_measures = len(spectra_region['flux'])

        # Continuum_model can be a fitted model or a fixed number
        if isinstance(continuum_model, float) or isinstance(continuum_model, int):
            cont_diff = median_flux - continuum_model
            cont_diff_percentage = np.abs(cont_diff) / continuum_model
        else:
            c = continuum_model(wave_base)
            cont_diff = median_flux - c
            cont_diff_percentage = np.abs(cont_diff) / c
        # Flux should be above the continuum model but no more than a given limit (1% by default)
        near_continuum = (cont_diff_percentage <= max_continuum_diff)

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
def find_continuum_on_regions(spectra, resolution, regions, log_filename=None, max_std_continuum = 0.002, continuum_model = 0.95, max_continuum_diff=0.01, fixed_wave_step=None, frame=None):

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
            median_flux = np.median(spectra_region['flux'])
            std_flux = np.std(spectra_region['flux'])
            num_measures = len(spectra_region['flux'])

            # Continuum_model can be a fitted model or a fixed number
            if isinstance(continuum_model, float) or isinstance(continuum_model, int):
                cont_diff = median_flux - continuum_model
                cont_diff_percentage = np.abs(cont_diff) / continuum_model
            else:
                c = continuum_model(wave_base)
                cont_diff = median_flux - c
                cont_diff_percentage = np.abs(cont_diff) / c
            # Flux should be above the continuum model but no more than a given limit (1% by default)
            near_continuum = (cont_diff_percentage <= max_continuum_diff)

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


