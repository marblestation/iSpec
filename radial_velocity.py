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
import matplotlib.pyplot as plt
import scipy.ndimage.filters
from common import *
from interpolate import *
from convolve import *
from lines import *
from mpfitmodels import GaussianModel
from pymodelfit import UniformCDFKnotSplineModel

# If velocity_step is specified, create a mask uniformly spaced in terms of velocity:
# - An increment in position (i => i+1) supposes a constant velocity increment (velocity_step)
# - An increment in position (i => i+1) does not implies a constant wavelength increment
# - On each wavelength position (i), the data_values with data_wave below the next wavelength position (i+1)
#   are added up
# - wave_grid is only used to identify the first and last wavelength
# - Useful for building the cross correlate function used for determining the radial velocity of a star
# If velocity_step is not specified, create a mask spaced in the same way of the wave_grid:
# - If the wave_grid is uniformly spaced in wavelength, an increment in position (i => i+1)
#   supposes a constant wavelength increment
# - Useful for cross correlate/compare one spectra with its radial velocity already corrected against
#   a linelist
def create_cross_correlation_mask(data_wave, data_value, wave_grid, velocity_step=None):
    # Speed of light in m/s
    c = 299792458.0

    if velocity_step != None:
        total_points = len(data_wave)
        wave_general_base = wave_grid[0]
        wave_general_top = wave_grid[-1]
        grid = []
        i = 0
        wave_base = wave_general_base
        while wave_base < wave_general_top and i < total_points:
            wave_top = wave_base + wave_base * ((velocity_step*1000) / c) # nm
            while i < total_points and data_wave[i] < wave_base:
                i += 1
            mask_value = 0
            while i < total_points and data_wave[i] >= wave_base and data_wave[i] < wave_top:
                mask_value += data_value[i]
                i += 1
            grid.append((wave_base, mask_value))
            wave_base = wave_top

        mask = np.array(grid, dtype=[('wave', float), ('value', float)])
    else:
        total_points = len(wave_grid)
        mask = np.recarray((total_points, ), dtype=[('wave', float), ('value', float)])
        mask['wave'] = wave_grid
        mask['value'] = np.zeros(total_points)

        mask_value = mask['value']

        max_j = len(data_wave)
        j = 0
        for i in np.arange(total_points-1):
            wave_base = wave_grid[i]
            wave_top = wave_grid[i+1]
            while j < max_j and data_wave[j] < wave_base:
                j += 1
            while j < max_j and data_wave[j] >= wave_base and data_wave[j] < wave_top:
                mask_value[i] += data_value[j]
                j += 1
            if j >= max_j:
                break
    return mask

# Calculates the cross correlation value between the spectra and the specified mask
# by shifting the mask from lower to upper velocity
# - The spectra and the mask should be uniformly spaced in terms of velocity (which
#   implies non-uniformly distributed in terms of wavelength)
# - The velocity step used for the construction of the mask should be the same
#   as the one specified in this function
# - The lower/upper/step velocity is only used to determine how many shifts
#   should be done (in array positions) and return a velocity grid
def cross_correlation_function(spectra, mask, lower_velocity_limit, upper_velocity_limit, velocity_step, frame=None):

    if frame != None:
        frame.update_progress(0)

    # Speed of light in m/s
    c = 299792458.0

    velocity = np.arange(lower_velocity_limit, upper_velocity_limit, velocity_step)
    # 1 shift = 0.5 km/s (or the specified value)
    shifts = np.int32(np.arange(lower_velocity_limit, upper_velocity_limit, velocity_step) / velocity_step)


    num_shifts = len(shifts)
    # Cross-correlation function
    ccf = np.zeros(num_shifts)
    for shift, i in zip(shifts, np.arange(num_shifts)):
        if shift == 0:
            shifted_mask = mask['value']
        elif shift > 0:
            shifted_mask = np.hstack((shift*[0], mask['value'][:-1*shift]))
        else:
            shifted_mask = np.hstack((mask['value'][-1*shift:], -1*shift*[0]))
        ccf[i] = np.correlate(spectra['flux'], shifted_mask)[0]
        if (i % 100) == 0:
            progress = ((i*1.0)/num_shifts) * 100
            print "%.2f%%" % (progress)
            if frame != None:
                frame.update_progress(progress)

    ccf = ccf/np.max(ccf) # Normalize

    return velocity, ccf


## Determines the velocity profile by cross-correlating the spectra with a mask
## built from a line list.
## Return the velocity coordenates and the normalized fluxes (relative intensities)
def build_velocity_profile(spectra, lines, lower_velocity_limit = -200, upper_velocity_limit = 200, velocity_step=1.0, frame=None):
    # Build a mask with non-uniform wavelength increments but constant in terms of velocity
    linelist_mask = create_cross_correlation_mask(lines['wave_peak'], lines['depth'], spectra['waveobs'], velocity_step)

    # Resampling spectra to match the wavelength grid of the mask
    # Speed of light in m/s
    c = 299792458.0
    interpolated_spectra = resample_spectra(spectra, linelist_mask['wave'])

    # Obtain the cross-correlate function by shifting the mask
    velocity, ccf = cross_correlation_function(interpolated_spectra, linelist_mask, lower_velocity_limit = lower_velocity_limit, upper_velocity_limit =upper_velocity_limit, velocity_step = velocity_step, frame=frame)
    return velocity, ccf, len(lines)


## Fits Gaussians to the deepest peaks in the velocity profile
## Return an array of fitted models
#  * For Radial Velocity profiles, more than 1 outlier peak implies that the star is a spectroscopic binary
def modelize_velocity_profile(xcoord, fluxes, only_one_peak=False):
    models = []
    if len(xcoord) == 0 or len(fluxes) == 0:
        return models

    # Smooth flux
    sig = 1
    smoothed_fluxes = scipy.ndimage.filters.gaussian_filter1d(fluxes, sig)
    # Finding peaks and base points
    peaks, base_points = find_peaks_and_base_points(xcoord, smoothed_fluxes)

    if len(peaks) == 0:
        return models

    # Fit continuum to normalize
    nknots = 2
    continuum_model = UniformCDFKnotSplineModel(nknots)
    continuum_model.fitData(xcoord[base_points], fluxes[base_points])
    fluxes /= continuum_model(xcoord)

    if len(peaks) != 0:
        # Adjusting edges
        base = base_points[:-1]
        top = base_points[1:]
        new_base = np.zeros(len(base), dtype=int)
        new_top = np.zeros(len(base), dtype=int)
        for i in np.arange(len(peaks)):
            new_base[i], new_top[i] = improve_linemask_edges(xcoord, fluxes, base[i], top[i], peaks[i])
        base = new_base
        top = new_top

        if only_one_peak:
            # Just try with the deepest line
            selected_peaks_indices = []
        else:
            # Identify peak outliers
            fluxes_not_outliers, selected_fluxes_not_outliers = sigma_clipping(fluxes[peaks], sig=3, meanfunc=np.median)
            selected_peaks_indices = np.arange(len(peaks))[~selected_fluxes_not_outliers]
            # FILTER: Make sure that the outliers selected are outliers because they are
            # deeper than the rest (we do not want outliers because they are less deep than the rest)
            interesting_outliers = fluxes[peaks[selected_peaks_indices]] < np.median(fluxes[peaks])
            selected_peaks_indices = selected_peaks_indices[interesting_outliers]
            # Sort the interesting peaks from more to less deep
            sorted_peaks_indices = np.argsort(fluxes[peaks[selected_peaks_indices]])
            selected_peaks_indices = selected_peaks_indices[sorted_peaks_indices]
            # Discard too small peaks (less than a third of the main peak)
            if len(selected_peaks_indices) >= 2:
                diff = (fluxes[peaks[selected_peaks_indices]] - fluxes[peaks[selected_peaks_indices]][0])
                wfilter = diff < fluxes[peaks[selected_peaks_indices]][0] / 3
                selected_peaks_indices = selected_peaks_indices[wfilter]

        if len(selected_peaks_indices) == 0:
            # Try with the deepest line
            sorted_peak_indices = np.argsort(fluxes[peaks])
            selected_peaks_indices = [sorted_peak_indices[0]]
    else:
        # If no peaks found, just consider the deepest point and mark the base and top
        # as the limits of the whole data
        sorted_fluxes_indices = np.argsort(fluxes)
        peaks = sorted_fluxes_indices[0]
        base = 0
        top = len(xcoord) - 1
        selected_peaks_indices = [0]

    for i in np.asarray(selected_peaks_indices):
        model = GaussianModel()

        # Parameters estimators
        baseline = np.median(fluxes[base_points])
        A = fluxes[peaks[i]] - baseline
        sig = np.abs(xcoord[top[i]] - xcoord[base[i]])/3.0
        mu = xcoord[peaks[i]]

        parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for j in np.arange(4)]
        parinfo[0]['value'] = 1#fluxes[base[i]] # baseline # Continuum
        parinfo[0]['fixed'] = True
        #parinfo[0]['limited'] = [False, True]
        #parinfo[0]['limits'] = [np.min([fluxes[base[i]], fluxes[top[i]]]), np.max([fluxes[base[i]], fluxes[top[i]]])]
        parinfo[1]['value'] = A # Only negative (absorption lines) and greater than the lowest point + 25%
        parinfo[1]['limited'] = [False, True]
        parinfo[1]['limits'] = [0., 0.]
        parinfo[2]['value'] = sig # Only positives (absorption lines)
        parinfo[2]['limited'] = [True, False]
        parinfo[2]['limits'] = [0., 0.]
        parinfo[3]['value'] = mu # Peak only within the xcoord slice
        parinfo[3]['limited'] = [True, True]
        parinfo[3]['limits'] = [xcoord[base[i]], xcoord[top[i]]]

        f = fluxes[base[i]:top[i]+1]
        min_flux = np.min(f)
        # More weight to the deeper fluxes
        if min_flux < 0:
            weights = f + -1*(min_flux) + 0.01 # Above zero
            weights = np.min(weights) / weights
        else:
            weights = min_flux / f

        try:
            #model.fitData(xcoord[base[i]:top[i]+1], fluxes[base[i]:top[i]+1], parinfo=parinfo, weights=weights)
            model.fitData(xcoord[base[i]:top[i]+1], fluxes[base[i]:top[i]+1], parinfo=parinfo)
            models.append(model)
        except Exception, e:
            print e.message

    #plt.plot(xcoord, fluxes)
    #if len(models) >= 1:
        #plt.plot(xcoord, models[0](xcoord))
    #plt.scatter(xcoord[base_points], fluxes[base_points])
    #plt.show()

    return np.asarray(models)

def find_confident_models(models, xcoord, fluxes):
    accept = []
    if len(models) == 0:
        return accept

    ## We want to calculate the mean and standard deviation of the velocity profile
    ## but discounting the effect of the deepest detected lines:
    # Build the fluxes for the composite models but ONLY for lines deeper than 0.97
    line_fluxes = None
    for model in models:
        if model(model.mu()) >= 0.97:
            continue
        if line_fluxes == None:
            line_fluxes = model(xcoord)
            continue

        current_line_fluxes = model(xcoord)
        wfilter = line_fluxes > current_line_fluxes
        line_fluxes[wfilter] = current_line_fluxes[wfilter]
    ### Substract the line models conserving the base level
    if line_fluxes != None:
        values = 1 + fluxes - line_fluxes
    else:
        values = fluxes
    ## Finally, calculate the mean and standard deviation
    check_mean = np.mean(values)
    check_std = np.std(values)
    for (i, model) in enumerate(models):
        mu = model.mu()
        peak = model(mu)

        # Discard peak if it is not deeper than mean flux + 6*standard deviation
        # unless it is the only detected peak
        limit = check_mean - 6*check_std
        if len(models) > 1 and (limit < 0.05 or peak >= limit):
            accept.append(False)
        else:
            accept.append(True)
    return np.asarray(accept)



## Constructs the velocity profile, fits a Gaussian and returns the mean (km/s)
## - If the lines are atomic lines, the radial velocity (RV) will be determined
## - If the lines are telluric lines, the barycentric velocity [+ RV if already applied] will be determined
def calculate_velocity(spectra, lines, lower_velocity_limit = -200, upper_velocity_limit = 200, velocity_step=1.0, renormalize=False, frame=None):
    xcoord, fluxes, num_lines = build_velocity_profile(spectra, lines, lower_velocity_limit=lower_velocity_limit, upper_velocity_limit=upper_velocity_limit, velocity_step=velocity_step, frame=frame)
    models = modelize_velocity_profile(xcoord, fluxes)
    if len(models) == 0:
        return 0
    else:
        # Velocity
        return np.round(models[0].mu, 2)  # km/s

# Velocity in km/s
def correct_velocity(spectra, velocity):
    # Speed of light in m/s
    c = 299792458.0
    # Radial/barycentric velocity from km/s to m/s
    velocity = velocity * 1000

    # Correct wavelength scale for radial velocity
    spectra['waveobs'] = spectra['waveobs'] / ((velocity / c) + 1)
    return spectra

# Velocity in km/s
def correct_velocity_regions(regions, velocity, with_peak=False):
    # Speed of light in m/s
    c = 299792458.0
    # Radial/barycentric velocity from km/s to m/s
    velocity = velocity * 1000

    regions['wave_base'] = regions['wave_base'] / ((velocity / c) + 1)
    regions['wave_top'] = regions['wave_top'] / ((velocity / c) + 1)
    if with_peak:
        regions['wave_peak'] = regions['wave_peak'] / ((velocity / c) + 1)
    return regions

# Shift regions considering the known radial/barycentric velocity of the star
def correct_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, velocity):
    # Speed of light in m/s
    c = 299792458.0
    # Radial/barycentric velocity from km/s to m/s
    velocity = velocity * 1000
    # Oposite sense because we want to correct regions, not the spectra
    velocity = velocity * -1

    continuum_regions = read_continuum_regions(continuum_regions_filename)
    continuum_regions['wave_base'] = continuum_regions['wave_base'] / ((velocity / c) + 1)
    continuum_regions['wave_top'] = continuum_regions['wave_top'] / ((velocity / c) + 1)
    write_continuum_regions(continuum_regions, continuum_regions_filename_out)


    # Lines
    line_regions = read_line_regions(line_regions_filename)
    line_regions['wave_base'] = line_regions['wave_base'] / ((velocity / c) + 1)
    line_regions['wave_top'] = line_regions['wave_top'] / ((velocity / c) + 1)
    line_regions['wave_peak'] = line_regions['wave_peak'] / ((velocity / c) + 1)
    write_line_regions(line_regions, line_regions_filename_out)

    segment_regions = read_segment_regions(segment_regions_filename)
    segment_regions['wave_base'] = segment_regions['wave_base'] / ((velocity / c) + 1)
    segment_regions['wave_top'] = segment_regions['wave_top'] / ((velocity / c) + 1)
    write_segment_regions(segment_regions, segment_regions_filename_out)


# Return the precession matrix needed to go from EQUINOX1 (i.e. 1950.0) to EQUINOX2 (i.e. 1975.0).
# Source: http://code.google.com/p/astrolibpy/source/browse/trunk/astrolib/
def premat(equinox1, equinox2, fk4=False):
    """
     INPUTS:
             EQUINOX1 - Original equinox of coordinates, numeric scalar.
             EQUINOX2 - Equinox of precessed coordinates.
     OUTPUT:
            matrix - double precision 3 x 3 precession matrix, used to precess
                        equatorial rectangular coordinates
     OPTIONAL INPUT KEYWORDS:
             fk4    - If this keyword is set, the FK4 (B1950.0) system precession
                        angles are used to compute the precession matrix.    The
                        default is to use FK5 (J2000.0) precession angles
     EXAMPLES:
             Return the precession matrix from 1950.0 to 1975.0 in the FK4 system
                matrix = premat( 1950.0, 1975.0, fk4=True)
     PROCEDURE:
             FK4 constants from "Computational Spherical Astronomy" by Taff (1983),
             p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
             Supplement 1992, page 104 Table 3.211.1.
    """

    deg_to_rad = pi / 180.0e0
    sec_to_rad = deg_to_rad / 3600.e0

    t = 0.001e0 * (equinox2 - equinox1)

    if not fk4:
        st = 0.001e0 * (equinox1 - 2000.e0)
        #  Compute 3 rotation angles
        a = sec_to_rad * t * (23062.181e0 + st * (139.656e0 + 0.0139e0 * st) + t * (30.188e0 - 0.344e0 * st + 17.998e0 * t))
        b = sec_to_rad * t * t * (79.280e0 + 0.410e0 * st + 0.205e0 * t) + a
        c = sec_to_rad * t * (20043.109e0 - st * (85.33e0 + 0.217e0 * st) + t * (-42.665e0 - 0.217e0 * st - 41.833e0 * t))
    else:
        st = 0.001e0 * (equinox1 - 1900.e0)
        #  Compute 3 rotation angles
        a = sec_to_rad * t * (23042.53e0 + st * (139.75e0 + 0.06e0 * st) + t * (30.23e0 - 0.27e0 * st + 18.0e0 * t))
        b = sec_to_rad * t * t * (79.27e0 + 0.66e0 * st + 0.32e0 * t) + a
        c = sec_to_rad * t * (20046.85e0 - st * (85.33e0 + 0.37e0 * st) + t * (-42.67e0 - 0.37e0 * st - 41.8e0 * t))

    sina = np.sin(a)
    sinb = np.sin(b)
    sinc = np.sin(c)
    cosa = np.cos(a)
    cosb = np.cos(b)
    cosc = np.cos(c)

    r = np.zeros((3, 3))
    r[0,:] = np.array([cosa * cosb * cosc - sina * sinb, sina * cosb + cosa * sinb * cosc, cosa * sinc])
    r[1,:] = np.array([-cosa * sinb - sina * cosb * cosc, cosa * cosb - sina * sinb * cosc, -sina * sinc])
    r[2,:] = np.array([-cosb * sinc, -sinb * sinc, cosc])

    return r

# Calculates heliocentric and barycentric velocity components of Earth.
# Source: http://code.google.com/p/astrolibpy/source/browse/trunk/astrolib/
def baryvel(datetime, deq=0):
    """
     EXPLANATION:
             BARYVEL takes into account the Earth-Moon motion, and is useful for
             radial velocity work to an accuracy of  ~1 m/s.
     CALLING SEQUENCE:
             dvel_hel, dvel_bary = baryvel(dje, deq)
     INPUTS:
             datetime - (year, month, day, [hour, minute, second])
             DEQ - (scalar) epoch of mean equinox of dvelh and dvelb. If deq=0
                        then deq is assumed to be equal to dje.
     OUTPUTS:
             DVELH: (vector(3)) heliocentric velocity component. in km/s
             DVELB: (vector(3)) barycentric velocity component. in km/s

             The 3-vectors DVELH and DVELB are given in a right-handed coordinate
             system with the +X axis toward the Vernal Equinox, and +Z axis
             toward the celestial pole.
     NOTES:
             Algorithm taken from FORTRAN program of Stumpff (1980, A&A Suppl, 41,1)
             Stumpf claimed an accuracy of 42 cm/s for the velocity.     A
             comparison with the JPL FORTRAN planetary ephemeris program PLEPH
             found agreement to within about 65 cm/s between 1986 and 1994
     EXAMPLE:
             Compute the radial velocity of the Earth toward Altair
                (19 50 46.999 +08 52 05.96) on 15-Feb-2012
                using the original Stumpf algorithm

heliocentric_vel, barycentric_vel = baryvel((2012, 2, 15, 00, 00, 00))

ra = (19.0 + 50.0/60 + 46.999/(60*60)) # hours
ra = ra * 360/24 # degrees
ra = ra * ((2*np.pi) / 360) # radians
dec = (8.0 + 52.0/60 + 5.96/(60*60)) # degrees
dec = dec * ((2*np.pi) / 360) # radians

# Project velocity toward star
vb = barycentric_vel
v = vb[0]*np.cos(dec)*np.cos(ra) + vb[1]*np.cos(dec)*np.sin(ra) + vb[2]*np.sin(dec)

measured_rv = -45.0 # Radial velocity measured from a spectrum
corrected_rv = measured_rv + v
    """
    dje = obstools.calendar_to_jd(datetime) # Julian ephemeris date.

    #Define constants
    dc2pi = 2 * np.pi
    cc2pi = 2 * np.pi
    dc1 = 1.0e0
    dcto = 2415020.0e0
    dcjul = 36525.0e0                            #days in Julian year
    dcbes = 0.313e0
    dctrop = 365.24219572e0                    #days in tropical year (...572 insig)
    dc1900 = 1900.0e0
    au = 1.4959787e8

    #Constants dcfel(i,k) of fast changing elements.
    dcfel = np.array([1.7400353e00, 6.2833195099091e02, 5.2796e-6, 6.2565836e00, 6.2830194572674e02, -2.6180e-6, 4.7199666e00, 8.3997091449254e03, -1.9780e-5, 1.9636505e-1, 8.4334662911720e03, -5.6044e-5, 4.1547339e00, 5.2993466764997e01, 5.8845e-6, 4.6524223e00, 2.1354275911213e01, 5.6797e-6, 4.2620486e00, 7.5025342197656e00, 5.5317e-6, 1.4740694e00, 3.8377331909193e00, 5.6093e-6])
    dcfel = np.reshape(dcfel, (8, 3))

    #constants dceps and ccsel(i,k) of slowly changing elements.
    dceps = np.array([4.093198e-1, -2.271110e-4, -2.860401e-8])
    ccsel = np.array([1.675104e-2, -4.179579e-5, -1.260516e-7, 2.220221e-1, 2.809917e-2, 1.852532e-5, 1.589963e00, 3.418075e-2, 1.430200e-5, 2.994089e00, 2.590824e-2, 4.155840e-6, 8.155457e-1, 2.486352e-2, 6.836840e-6, 1.735614e00, 1.763719e-2, 6.370440e-6, 1.968564e00, 1.524020e-2, -2.517152e-6, 1.282417e00, 8.703393e-3, 2.289292e-5, 2.280820e00, 1.918010e-2, 4.484520e-6, 4.833473e-2, 1.641773e-4, -4.654200e-7, 5.589232e-2, -3.455092e-4, -7.388560e-7, 4.634443e-2, -2.658234e-5, 7.757000e-8, 8.997041e-3, 6.329728e-6, -1.939256e-9, 2.284178e-2, -9.941590e-5, 6.787400e-8, 4.350267e-2, -6.839749e-5, -2.714956e-7, 1.348204e-2, 1.091504e-5, 6.903760e-7, 3.106570e-2, -1.665665e-4, -1.590188e-7])
    ccsel = np.reshape(ccsel, (17, 3))

    #Constants of the arguments of the short-period perturbations.
    dcargs = np.array([5.0974222e0, -7.8604195454652e2, 3.9584962e0, -5.7533848094674e2, 1.6338070e0, -1.1506769618935e3, 2.5487111e0, -3.9302097727326e2, 4.9255514e0, -5.8849265665348e2, 1.3363463e0, -5.5076098609303e2, 1.6072053e0, -5.2237501616674e2, 1.3629480e0, -1.1790629318198e3, 5.5657014e0, -1.0977134971135e3, 5.0708205e0, -1.5774000881978e2, 3.9318944e0, 5.2963464780000e1, 4.8989497e0, 3.9809289073258e1, 1.3097446e0, 7.7540959633708e1, 3.5147141e0, 7.9618578146517e1, 3.5413158e0, -5.4868336758022e2])
    dcargs = np.reshape(dcargs, (15, 2))

    #Amplitudes ccamps(n,k) of the short-period perturbations.
    ccamps = np.array([-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5, -2.490817e-7, -3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5, -1.823138e-7, 6.593466e-7, 1.322572e-5, 9.258695e-6, -4.674248e-7, -3.646275e-7, 1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7, 9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7, 7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7, -2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6, -1.655307e-7, -3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6, -3.736225e-7, 3.442177e-7, 2.671323e-6, 1.832858e-6, -2.394688e-7, -3.478444e-7, 8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8, -1.488378e-6, -1.251789e-5, 5.226868e-7, -2.049301e-7, 0.e0, -8.043059e-6, -2.991300e-6, 1.473654e-7, -3.154542e-7, 0.e0, 3.699128e-6, -3.316126e-6, 2.901257e-7, 3.407826e-7, 0.e0, 2.550120e-6, -1.241123e-6, 9.901116e-8, 2.210482e-7, 0.e0, -6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.e0])
    ccamps = np.reshape(ccamps, (15, 5))

    #Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
    ccsec3 = -7.757020e-8
    ccsec = np.array([1.289600e-6, 5.550147e-1, 2.076942e00, 3.102810e-5, 4.035027e00, 3.525565e-1, 9.124190e-6, 9.990265e-1, 2.622706e00, 9.793240e-7, 5.508259e00, 1.559103e01])
    ccsec = np.reshape(ccsec, (4, 3))

    #Sidereal rates.
    dcsld = 1.990987e-7                         #sidereal rate in longitude
    ccsgd = 1.990969e-7                         #sidereal rate in mean anomaly

    #Constants used in the calculation of the lunar contribution.
    cckm = 3.122140e-5
    ccmld = 2.661699e-6
    ccfdi = 2.399485e-7

    #Constants dcargm(i,k) of the arguments of the perturbations of the motion
    # of the moon.
    dcargm = np.array([5.1679830e0, 8.3286911095275e3, 5.4913150e0, -7.2140632838100e3, 5.9598530e0, 1.5542754389685e4])
    dcargm = np.reshape(dcargm, (3, 2))

    #Amplitudes ccampm(n,k) of the perturbations of the moon.
    ccampm = np.array([1.097594e-1, 2.896773e-7, 5.450474e-2, 1.438491e-7, -2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8, 1.148966e-2, 5.658888e-8, 8.249439e-3, 4.063015e-8])
    ccampm = np.reshape(ccampm, (3, 4))

    #ccpamv(k)=a*m*dl,dt (planets), dc1mme=1-mass(earth+moon)
    ccpamv = np.array([8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12])
    dc1mme = 0.99999696e0

    #Time arguments.
    dt = (dje - dcto) / dcjul
    tvec = np.array([1e0, dt, dt * dt])

    #Values of all elements for the instant(aneous?) dje.
    temp = (np.transpose(np.dot(np.transpose(tvec), np.transpose(dcfel)))) % dc2pi
    dml = temp[0]
    forbel = temp[1:8]
    g = forbel[0]                                 #old fortran equivalence

    deps = (tvec * dceps).sum() % dc2pi
    sorbel = (np.transpose(np.dot(np.transpose(tvec), np.transpose(ccsel)))) % dc2pi
    e = sorbel[0]                                 #old fortran equivalence

    #Secular perturbations in longitude.
    dummy = np.cos(2.0)
    sn = np.sin((np.transpose(np.dot(np.transpose(tvec[0:2]), np.transpose(ccsec[:,1:3])))) % cc2pi)

    #Periodic perturbations of the emb (earth-moon barycenter).
    pertl = (ccsec[:,0] * sn).sum() + dt * ccsec3 * sn[2]
    pertld = 0.0
    pertr = 0.0
    pertrd = 0.0
    for k in range(0, 15):
        a = (dcargs[k,0] + dt * dcargs[k,1]) % dc2pi
        cosa = np.cos(a)
        sina = np.sin(a)
        pertl = pertl + ccamps[k,0] * cosa + ccamps[k,1] * sina
        pertr = pertr + ccamps[k,2] * cosa + ccamps[k,3] * sina
        if k < 11:
            pertld = pertld + (ccamps[k,1] * cosa - ccamps[k,0] * sina) * ccamps[k,4]
            pertrd = pertrd + (ccamps[k,3] * cosa - ccamps[k,2] * sina) * ccamps[k,4]

    #Elliptic part of the motion of the emb.
    phi = (e * e / 4e0) * (((8e0 / e) - e) * np.sin(g) + 5 * np.sin(2 * g) + (13 / 3e0) * e * np.sin(3 * g))
    f = g + phi
    sinf = np.sin(f)
    cosf = np.cos(f)
    dpsi = (dc1 - e * e) / (dc1 + e * cosf)
    phid = 2 * e * ccsgd * ((1 + 1.5 * e * e) * cosf + e * (1.25 - 0.5 * sinf * sinf))
    psid = ccsgd * e * sinf / np.sqrt(dc1 - e * e)

    #Perturbed heliocentric motion of the emb.
    d1pdro = dc1 + pertr
    drd = d1pdro * (psid + dpsi * pertrd)
    drld = d1pdro * dpsi * (dcsld + phid + pertld)
    dtl = (dml + phi + pertl) % dc2pi
    dsinls = np.sin(dtl)
    dcosls = np.cos(dtl)
    dxhd = drd * dcosls - drld * dsinls
    dyhd = drd * dsinls + drld * dcosls

    #Influence of eccentricity, evection and variation on the geocentric
    # motion of the moon.
    pertl = 0.0
    pertld = 0.0
    pertp = 0.0
    pertpd = 0.0
    for k in range(0, 3):
        a = (dcargm[k,0] + dt * dcargm[k,1]) % dc2pi
        sina = np.sin(a)
        cosa = np.cos(a)
        pertl = pertl + ccampm[k,0] * sina
        pertld = pertld + ccampm[k,1] * cosa
        pertp = pertp + ccampm[k,2] * cosa
        pertpd = pertpd - ccampm[k,3] * sina

    #Heliocentric motion of the earth.
    tl = forbel[1] + pertl
    sinlm = np.sin(tl)
    coslm = np.cos(tl)
    sigma = cckm / (1.0 + pertp)
    a = sigma * (ccmld + pertld)
    b = sigma * pertpd
    dxhd = dxhd + a * sinlm + b * coslm
    dyhd = dyhd - a * coslm + b * sinlm
    dzhd = -sigma * ccfdi * np.cos(forbel[2])

    #Barycentric motion of the earth.
    dxbd = dxhd * dc1mme
    dybd = dyhd * dc1mme
    dzbd = dzhd * dc1mme
    for k in range(0, 4):
        plon = forbel[k + 3]
        pomg = sorbel[k + 1]
        pecc = sorbel[k + 9]
        tl = (plon + 2.0 * pecc * np.sin(plon - pomg)) % cc2pi
        dxbd = dxbd + ccpamv[k] * (np.sin(tl) + pecc * np.sin(pomg))
        dybd = dybd - ccpamv[k] * (np.cos(tl) + pecc * np.cos(pomg))
        dzbd = dzbd - ccpamv[k] * sorbel[k + 13] * np.cos(plon - sorbel[k + 5])


    #Transition to mean equator of date.
    dcosep = np.cos(deps)
    dsinep = np.sin(deps)
    dyahd = dcosep * dyhd - dsinep * dzhd
    dzahd = dsinep * dyhd + dcosep * dzhd
    dyabd = dcosep * dybd - dsinep * dzbd
    dzabd = dsinep * dybd + dcosep * dzbd

    #Epoch of mean equinox (deq) of zero implies that we should use
    # Julian ephemeris date (dje) as epoch of mean equinox.
    if deq == 0:
        dvelh = au * (np.array([dxhd, dyahd, dzahd]))
        dvelb = au * (np.array([dxbd, dyabd, dzabd]))
        return (dvelh,dvelb)

    #General precession from epoch dje to deq.
    deqdat = (dje - dcto - dcbes) / dctrop + dc1900
    prema = premat(deqdat, deq, fk4=True)

    dvelh = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxhd, dyahd, dzahd])))))
    dvelb = au * (np.transpose(np.dot(np.transpose(prema), np.transpose(np.array([dxbd, dyabd, dzabd])))))

    return (dvelh, dvelb)

# Find the FWHM in km/s of the telluric lines. Options:
# A) Use the FWHM already fitted and saved in the telluric line list (DEFAULT)
# B) Modelize velocity profile of the Telluric lines
def telluric_fwhm_correction(linelist_telluric, wmin, wmax, use_already_fitted_data=True, lower_velocity_limit = -100, upper_velocity_limit = 100, velocity_step=0.5):
    if use_already_fitted_data:
        # Light speed
        c = 299792458.0 # m/s
        lfilter = (linelist_telluric['wave_peak'] >= wmin) & (linelist_telluric['wave_peak'] <= wmax)
        telluric_fwhm = np.mean((c / (linelist_telluric['wave_peak'][lfilter] / linelist_telluric['fwhm'][lfilter])) / 1000.0) # km/s
    else:
        # USED ONLY HERE: lower_velocity_limit = -100, upper_velocity_limit = 100, velocity_step=0.5
        telluric_spectra_file = "input/telluric/standard_atm_air_norm.s.gz"
        telluric_spectra = read_spectra(telluric_spectra_file)

        wfilter = (telluric_spectra['waveobs'] > wmin) & (telluric_spectra['waveobs'] < wmax)
        telluric_spectra = telluric_spectra[wfilter]

        if len(telluric_spectra) == 0:
            raise Exception("Wavelength limits: Out of range for telluric lines")

        xcoord, fluxes, num_lines = build_velocity_profile(telluric_spectra, linelist_telluric, lower_velocity_limit=lower_velocity_limit, upper_velocity_limit=upper_velocity_limit, velocity_step=velocity_step, frame=None)
        models = modelize_velocity_profile(xcoord, fluxes, only_one_peak=True)
        c = 299792458.0 # m/s
        telluric_fwhm = models[0].sig() * (2*np.sqrt(2*np.log(2)))

    telluric_R = np.int(c/(1000.0*np.round(telluric_fwhm, 2)))
    print "Telluric R:", telluric_R
    print "Telluric FWHM correction:", telluric_fwhm

    return telluric_fwhm



if __name__ == '__main__':
    #~ continuum_regions_filename = "input/LUMBA/UVES_MRD_sun_cmask.s.gz"
    #~ line_regions_filename = "input/LUMBA/UVES_MRD_sun_Fe-linelist.s.gz"
    #~ segment_regions_filename = "input/LUMBA/UVES_MRD_sun_segments.s.gz"
#~
    #~ radial_vel = -97.2 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+cas+a
#~
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_cmask.s.gz"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_Fe-linelist.s.gz"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_segments.s.gz"
    #~
    #~ correct_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~
    #~ radial_vel = 14.03 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+leo
#~
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_cmask.s.gz"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_Fe-linelist.s.gz"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_segments.s.gz"
    #~
    #~ correct_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~
    #~ radial_vel = 0 # The original spectra for Arcturus has been corrected for radial velocity
    #~ #radial_vel = -5.19 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=arcturus
#~
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_cmask.s.gz"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_Fe-linelist.s.gz"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_segments.s.gz"
    #~
    #~ correct_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)

    #star, resolution = "input/spectra/examples/espadons_mu_leo.s.gz", 80000
    #star, resolution = "input/spectra/examples/espadons_mu_leo_norm.s.gz", 80000
    #star, resolution = "input/spectra/examples/harps_procyon.s.gz", 115000
    #star, resolution = "input/spectra/examples/harps_procyon_norm.s.gz", 115000
    #star, resolution = "input/spectra/examples/narval_arcturus.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_arcturus_norm.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_mu_cas.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_mu_cas_norm.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_sun.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_sun_norm.s.gz", 80000
    #star, resolution = "input/spectra/binaries/elodie_hd005516A_spectroscopic_binary.s.gz", 42000
    #star, resolution = "input/spectra/binaries/elodie_hd085503_single_star.s.gz", 42000
    #star, resolution = "input/spectra/instruments/elodie_hd146233_SN237_normalized.s.gz", 42000
    #star, resolution = "input/spectra/instruments/elodie_hd146233_SN237.s.gz", 42000
    #star, resolution = "input/spectra/instruments/giraffe_hd107328_normalized.s.gz", 16000
    #star, resolution = "input/spectra/instruments/giraffe_hd107328.s.gz", 16000
    #star, resolution = "input/spectra/instruments/narval_hd146233_normalized.s.gz", 80000
    #star, resolution = "input/spectra/instruments/narval_hd146233.s.gz", 80000
    #star, resolution = "input/spectra/instruments/uves_hd146233.s.gz", 47000
    #star, resolution = "input/spectra/instruments/uves_hd146233_normalized.s.gz", 47000
    #star, resolution = "input/spectra/instruments/espadons_hd85503.s.gz", 80000
    #star, resolution = "input/spectra/instruments/harps_procyon.s.gz", 115000
    #star, resolution = "input/spectra/telluric_standards/narval_hr1567_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr2845_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr3492_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr3982_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr4828_002.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr708_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr7906_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr838_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr4182_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr5867_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr7235_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr8028_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr8976_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr4828_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr6629_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr7528_001.s.gz", 80000
    #star, resolution = "input/spectra/telluric_standards/narval_hr804_001.s.gz", 80000
    #star, resolution = "input/spectra/synthetic/synth_LUMBA_Gustafsson_SME_arcturus.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_LUMBA_Gustafsson_SME_mu_cas_a.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_LUMBA_Gustafsson_SME_mu_leo.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_LUMBA_Gustafsson_SME_sun.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_kurucz_arcturus.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_kurucz_mu_cas_a.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_kurucz_mu_leo.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_kurucz_sun.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_kurucz_arcturus.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_kurucz_mu_cas_a.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_kurucz_mu_leo.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_kurucz_sun.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_castelli_arcturus.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_castelli_arcturus.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_castelli_mu_cas_a.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_castelli_mu_leo.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_SPEC_castelli_sun.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_castelli_mu_cas_a.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_castelli_mu_leo.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/synth_VALD_castelli_sun.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/telluric_standard_atm_air_model.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/telluric_standard_atm_air_model_norm.s.gz", 200000



    print "Reading spectrum..."
    spectra = read_spectra(star)
    wmin = spectra['waveobs'][0]
    wmax = spectra['waveobs'][-1]

    ## Convolution test
    #from_resolution = 42000
    #to_resolution = 20000
    #spectra = convolve_spectra(spectra, from_resolution, to_resolution)


    lower_velocity_limit = -100
    upper_velocity_limit = 100
    velocity_step=0.5

    print "Reading telluric lines..."
    telluric_lines_file = "input/telluric/standard_atm_air.lst"
    linelist_telluric = read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)
    print "Selecting telluric lines..."
    # Limit to region of interest
    # Light speed
    c = 299792458.0 # m/s
    delta_wmin = wmin * (lower_velocity_limit / (c/1000.0))
    delta_wmax = wmax * (upper_velocity_limit / (c/1000.0))
    wfilter = (linelist_telluric['wave_peak'] <= wmax + delta_wmax) & (linelist_telluric['wave_peak'] >= wmin + delta_wmin)
    linelist_telluric = linelist_telluric[wfilter]
    # Discard not fitted lines
    rfilter = linelist_telluric['rms'] == 9999
    linelist_telluric = linelist_telluric[~rfilter]
    # Discard too deep or too small lines
    rfilter = (linelist_telluric['depth'] <= 0.9) & (linelist_telluric['depth'] >= 0.01)
    linelist_telluric = linelist_telluric[rfilter]
    # Discard outliers FWHM in km/s (which is not wavelength dependent)
    telluric_fwhm = (c / (linelist_telluric['wave_peak'] / linelist_telluric['fwhm'])) / 1000.0 # km/s
    fwhm_selected, fwhm_selected_filter = sigma_clipping(telluric_fwhm, meanfunc=np.median)
    linelist_telluric = linelist_telluric[fwhm_selected_filter]

    #### Telluric fwhm
    telluric_fwhm = telluric_fwhm_correction(linelist_telluric, wmin, wmax)

    print "Measuring global resolution..."
    xcoord, fluxes, num_lines = build_velocity_profile(spectra, linelist_telluric, lower_velocity_limit=lower_velocity_limit, upper_velocity_limit=upper_velocity_limit, velocity_step=velocity_step, frame=None)
    models = modelize_velocity_profile(xcoord, fluxes, only_one_peak=True)
    c = 299792458.0 # m/s
    fwhm = models[0].sig() * (2*np.sqrt(2*np.log(2)))
    fwhm -= telluric_fwhm
    global_R = np.int(c/(1000.0*np.round(fwhm, 2)))
    velocity = models[0].mu()
    print spectra['waveobs'][0], ":", spectra['waveobs'][-1], "=>", global_R
    #plt.plot(xcoord, fluxes)
    #plt.plot(xcoord, models[0](xcoord))
    #plt.show()


    ### Convolution test
    #from_resolution = 5000000
    #to_resolution = 42000
    #telluric_spectra = convolve_spectra(telluric_spectra, from_resolution, to_resolution)

    print "Creating in function of wavelength..."
    wave_range = np.arange(wmin, wmax+1, 1)
    R = np.zeros(len(wave_range))

    ## Chunks
    total_chuncks = len(wave_range)
    max_blocks = 50
    for i in np.arange(total_chuncks-1):
        if i+max_blocks >= total_chuncks: break
        blocks = np.min([max_blocks, (total_chuncks-1)-i])
        wmin_block =  wave_range[i]
        wmax_block = wave_range[i+blocks]
        # Line filter
        lfilter = (linelist_telluric['wave_peak'] >= wmin_block) & (linelist_telluric['wave_peak'] <= wmax_block)
        # Spectra filter
        wfilter = (spectra['waveobs'] >= wmin_block) & (spectra['waveobs'] <= wmax_block)
        if len(spectra[wfilter]) == 0 or len(linelist_telluric[lfilter]) == 0:
            models = []
        else:
            xcoord, fluxes, num_lines = build_velocity_profile(spectra[wfilter], linelist_telluric[lfilter], lower_velocity_limit=lower_velocity_limit, upper_velocity_limit=upper_velocity_limit, velocity_step=velocity_step, frame=None)

            models = modelize_velocity_profile(xcoord, fluxes, only_one_peak=True)

        if len(models) >= 1 and np.abs(models[0].mu() - velocity) <= 2:
            c = 299792458.0 # m/s
            fwhm = models[0].sig() * (2*np.sqrt(2*np.log(2)))
            telluric_fwhm = telluric_fwhm_correction(linelist_telluric, wmin_block, wmax_block)
            fwhm -= telluric_fwhm
            R[i] = np.int(c/(1000.0*np.round(fwhm, 2)))
        #else:
            #plt.plot(xcoord, fluxes)
            #plt.show()
        print wmin_block, ":", wmax_block, "=>", R[i]

    import os
    lfilter = (linelist_telluric['wave_peak'] >= wmin) & (linelist_telluric['wave_peak'] <= wmax)
    ax1 = plt.subplot(211)
    ax1.plot([wmin, wmax], [global_R, global_R], label="Global R")
    ax1.plot(wave_range[R!=0], R[R!=0], label="Local averaged R")
    plt.grid()
    plt.legend()
    plt.ylabel("Resolution")
    plt.title("Resolution in blocks of 100 nm for " + os.path.basename(star))
    plt.subplot(212, sharex=ax1)
    plt.bar(linelist_telluric['wave_peak'][lfilter], linelist_telluric['depth'][lfilter], label="Mask", width=0.01, color="black")
    plt.plot([wmin, wmax], [0.0, 0.0], color="black")
    plt.plot(spectra['waveobs'], spectra['flux']-1, label="Spectra")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Depth")
    plt.legend()
    plt.show()


## create_cross_correlation_mask test
#velocity = 0.23
#velocity_step = 0.5
#lines = linelist_telluric
#w1 = create_cross_correlation_mask(lines['wave_peak'], lines['depth'], np.arange(300,900,100), velocity_step)
#r1 = (c/1000) * ((w1['wave'][1:] - w1['wave'][:-1]) / (w1['wave'][:-1] + ((w1['wave'][1:] - w1['wave'][:-1])/2)))
#valid1 = np.max(r1) - np.min(r1) < 0.01
#w2 = create_cross_correlation_mask(lines['wave_peak'], lines['depth'], np.arange(300,900,100) / ((velocity / c) + 1), velocity_step)
#r2 = (c/1000) * ((w2['wave'][1:] - w2['wave'][:-1]) / (w2['wave'][:-1] + ((w2['wave'][1:] - w2['wave'][:-1])/2)))
#valid2 = np.max(r2) - np.min(r2) < 0.01
#r = c * ((w1['wave'] - w2['wave']) / (w1['wave'] + ((w1['wave'] - w2['wave'])/2)))
#valid = np.max(r) - np.min(r) < 0.01

#spectra = np.recarray((len(np.arange(300,900,100)), ), dtype=[('waveobs', float),('flux', float),('err', float)])
#spectra['waveobs'] = np.arange(300,900,100)
#spectra['flux'] = np.arange(len(spectra))

#s = resample_spectra(spectra, w2['wave'])
#resample_spectra(spectra, w1['wave'])
#plt.plot(s['waveobs'], s['flux'])
#plt.plot(np.arange(300,900,100), np.arange(len(spectra))+1)
#plt.show()





