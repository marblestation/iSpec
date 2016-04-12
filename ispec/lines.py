#
#    This file is part of the Integrated Spectroscopic Framework (iSpec).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#

from astropy.io import ascii
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions
from scipy.fftpack import fft
from scipy.fftpack import ifft
import scipy.ndimage.filters
import scipy.stats as stats
import cPickle as pickle
import gzip
import os
import subprocess
import shutil
from common import *
from continuum import *
from lines import *
from spectrum import *
import matplotlib.pyplot as plt
from mpfitmodels import GaussianModel
from mpfitmodels import VoigtModel
import log
import logging
import copy
import re



########################################################################
## [START] LINE LISTS
########################################################################

def read_telluric_linelist(telluric_lines_file, minimum_depth=0.0):
    """
    Read telluric linelist.
    """
    telluric_lines = ascii.read(telluric_lines_file, delimiter="\t").as_array()

    # Convert string to bool
    telluric_lines = rfn.rename_fields(telluric_lines, {'discarded':'discarded_string',})

    # Sort by wave_peak and descending depth (for that we create a temporary field)
    telluric_lines = rfn.append_fields(telluric_lines, "discarded", dtypes=bool, data=(telluric_lines['discarded_string'] == "True"))

    telluric_lines = rfn.drop_fields(telluric_lines, ['discarded_string'])


    if minimum_depth <= 0.0:
        return telluric_lines
    else:
        # Filter
        telluric_lines_limited = telluric_lines[telluric_lines['depth'] >= minimum_depth]
        return telluric_lines_limited


def __get_element(chemical_elements, molecules, species):
    """
    Convert species code form by the atomic number + "." + ionization state to element names type "Fe 1" or "Fe 2".
    Returns "Discard" if not found.
    """
    atomic_num, ion_state = species.split(".")

    tfilter = (chemical_elements['atomic_num'] == int(atomic_num))
    if len(chemical_elements["symbol"][tfilter]) == 0:
        # Symbol not found, maybe it is a molecule
        mfilter = (molecules['atomic_num'] == int(atomic_num))
        if len(molecules["symbol"][mfilter]) == 0:
            print "WARNING: Discarding lines with atomic number", int(atomic_num)
            return "Discard"
        else:
            symbol = str(molecules["symbol"][mfilter][0])
            #print symbol
    else:
        symbol = str(chemical_elements["symbol"][tfilter][0])

    return symbol + " " + str(int(ion_state) + 1)


def __get_upper_state(lower_state, wavelength):
    """
    Calculate upper excitation level from lower and wavelength.
    Units:

    * eV for lower excitation level
    * nm for wavelength
    """
    # Planck constant
    h = 6.62606957 * 10e-34 # m^2 kg / s
    # Light speed in vacuum
    c = 299792458.0 # m/s

    # Wavelength
    l = wavelength * 10e-9 # m
    # Frequency
    f = c/l # Hz
    # Energy
    E = h * f # Joules
    E = E * 6.24150974e18 # electron Volt (eV)
    return lower_state + E # eV


def __inverse_cm_to_eV(value):
    """
    Units transformation from cm^-1 to eV
    """
    return value / 8065.544 # eV

def __eV_to_inverse_cm(value):
    """
    Units transformation from eV to cm^-1.
    """
    return value * 8065.544 # cm^-1



########################################################################
## [END] LINE LIST
########################################################################


def read_isotope_data(isotopes_file):
    """
    Read isotope information.
    """
    isotope = ascii.read(isotopes_file, names=['atomic_code', 'mass_number', 'molecular_weight', 'relative_abundance_in_the_solar_system']).as_array()
    return isotope

def write_isotope_data(isotope, isotope_filename=None, tmp_dir=None):
    """
    Write isotope information.
    """
    if isotope_filename is not None:
        out = open(isotope_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
    #  1.0   1    1.007825  0.999885
    #  1.0   2    2.0140    0.000115
    #  2.0   3    3.016029  0.00000137
    #  2.0   4    4.002603  0.99999863
    out.write("\n".join(["  ".join(map(str, (iso['atomic_code'], iso['mass_number'], iso['molecular_weight'], iso['relative_abundance_in_the_solar_system']))) for iso in isotope]))
    out.close()
    return out.name


def read_cross_correlation_mask(mask_file):
    """
    Read linelist for building a mask and doing cross-correlation with the
    cross correlation function.

    The linelist should contain at least the columns wave_peak (in nm) and depth.
    """
    ccf_mask = ascii.read(mask_file).as_array()
    return ccf_mask

def read_chemical_elements(chemical_elements_file):
    """
    Read information about the chemical elements (periodic table)
    """
    chemical_elements = ascii.read(chemical_elements_file, delimiter="\t").as_array()
    return chemical_elements

def read_molecular_symbols(molecules_file):
    """
    Read information about the molecules (symbol and "atomic number")
      - For diatomic molecules, the atomic_num specifies the atomic makeup of the molecule.
        Thus, H2 is 101.0, the two ``1''s referring to the two hydrogens, CH is 106.0,
        CO 608.0, MgH 112.0, TiO 822.0, etc.
      - The lightest element always comes first in the code, so that 608.0 cannot be
        confused with NdO, which would be written 860.0.
    """
    molecules = ascii.read(molecules_file, delimiter="\t").as_array()
    return molecules


def read_line_regions(line_regions_filename):
    """
    Read line regions.
    Line region files should be plain text files with **tab** character as column delimiter.
    Four columns should exists: 'wave_peak', 'wave_base', 'wave_top' and 'note'
    (the first line should contain those header names).
    They indicate the peak of the line, beginning and end of each region (one per line)
    and a comment (it can be any string comment). For example:
    ::

        wave_peak       wave_base       wave_top        note
        480.8148        480.7970        480.8330        Fe 1
        496.2572        496.2400        496.2820        Fe 1
        499.2785        499.2610        499.2950
        505.8498        505.8348        505.8660        Fe 1

    The note can be blank but the previous **tab** character should exists anyway.

    The function detects automatically if the file being read was saved with
    write_line_regions(..., extended=True), which contains more fields corresponding
    to atomic cross-match and line fits.
    """
    atomic_dtype = _get_atomic_linelist_definition()
    fitted_dtype = __get_fitted_lines_definition()
    try:
        num_cols = len(open(line_regions_filename, "r").readline().split("\t"))
        if num_cols == len(atomic_dtype) + len(fitted_dtype):
            # Atomic + fitted information
            line_regions = np.genfromtxt(line_regions_filename, dtype=atomic_dtype+fitted_dtype, delimiter="\t", skip_header=1, loose=False)
        elif num_cols == 4:
            line_regions = np.array([tuple(line.rstrip('\r\n').split("\t")) for line in open(line_regions_filename,)][1:], dtype=[('wave_peak', float),('wave_base', float),('wave_top', float), ('note', '|S100')])
        else:
            raise Exception()
    except:
        raise Exception("Wrong line regions file format!")

    if np.any(line_regions['wave_top'] - line_regions['wave_base'] <= 0):
        logging.error("Line regions where wave_base is equal or bigger than wave_top")
        raise Exception("Incompatible format")

    if np.any(line_regions['wave_top'] - line_regions['wave_peak'] < 0) or np.any(line_regions['wave_peak'] - line_regions['wave_base'] < 0):
        logging.error("Line regions where wave_peak is outside wave_base and wave_top")
        raise Exception("Incompatible format")

    return line_regions

def write_line_regions(line_regions, line_regions_filename, extended=False):
    """
    Write line regions file with the following format:
    ::

        wave_peak       wave_base       wave_top        note
        480.8148        480.7970        480.8330        Fe 1
        496.2572        496.2400        496.2820        Fe 1
        499.2785        499.2610        499.2950
        505.8498        505.8348        505.8660        Fe 1

    Except if extended is True, in which case the complete cross-matched
    atomic data and fitted data is also saved.
    """
    if extended:
        out = open(line_regions_filename, "w")
        __generic_write_atomic_linelist_header(out, include_fit=True)
        for line in line_regions:
            __generic_write_atomic_linelist_element(out, line, include_fit=True)
        out.close()
    else:
        out = open(line_regions_filename, "w")
        out.write("wave_peak\twave_base\twave_top\tnote\n")
        out.write("\n".join(["\t".join(map(str, (line['wave_peak'], line['wave_base'], line['wave_top'], line['note']))) for line in line_regions]))
        out.close()
        # 5934.6545  26.0  0  31689  48530  -1.07  1.0  AO  7.78  -5.29  959.247  Fe_1  False  0.04  0.0



def __fit_gaussian(spectrum_slice, continuum_model, mu, sig=None, A=None, baseline_margin=0., free_mu=False):
    """
    Fits a gaussian at a given wavelength location using a fitted continuum model.

    - For absorption lines, it will alway be true: model.A() < 0 and model.sig() > 0.
    - A mu parameter (model.mu) outside the region used for the fitting it is also a symptom of bad fit.
    - If baseline_margin is different from 0, the continuum will be adjusted by letting
      free the baseline in between [(1.-baseline_margin)*baseline, (1+baseline_margin)*baseline]
    """
    model = GaussianModel()
    x = spectrum_slice['waveobs']
    y = spectrum_slice['flux']
    min_flux = np.min(y)

    # Parameters estimators
    baseline = np.median(continuum_model(spectrum_slice['waveobs']))
    if A is None:
        A = np.min((min_flux - baseline, -1e-10))
    if sig is None:
        sig = np.max(((x[-1] - x[0])/3.0, 1e-10))
    if free_mu and (mu < x[0] or mu > x[-1]):
        # Correct mu in case it is outside limits
        mu = spectrum_slice['waveobs'][np.argmin(spectrum_slice['flux'])]


    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(4)]
    parinfo[0]['value'] = baseline # Continuum
    if baseline_margin == 0:
        parinfo[0]['fixed'] = True
    else:
        parinfo[0]['fixed'] = False
        parinfo[0]['limited'] = [True, True]
        #parinfo[0]['limits'] = [0.95*baseline, 1.05*baseline]
        parinfo[0]['limits'] = [(1.-baseline_margin)*baseline, (1+baseline_margin)*baseline]
    parinfo[1]['value'] = A # Only negative (absorption lines) and greater than the lowest point + 25%
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [np.min(((min_flux-baseline) * 1.25, -1e-10)), -1e-10]
    parinfo[2]['value'] = sig # Only positives (absorption lines) and lower than the spectrum slice
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [1e-10, x[-1] - x[0]]
    parinfo[3]['value'] = mu # Peak only within the spectrum slice
    if not free_mu:
        parinfo[3]['fixed'] = True
    else:
        parinfo[3]['fixed'] = False
        parinfo[3]['limited'] = [True, True]
        parinfo[3]['limits'] = [x[0], x[-1]]
        #parinfo[3]['limits'] = [mu - 0.001, mu + 0.001]

    # If there are only 3 data point, fix 'mu'
    # - if not, the fit will fail with an exception because there are not
    #   more data points than parameters
    if len(spectrum_slice) == 3:
        parinfo[3]['fixed'] = True


    f = spectrum_slice['flux']
    min_flux = np.min(f)
    # More weight to the deeper fluxes, the upper part tend to be more contaminated by other lines
    if min_flux < 0:
        weights = f + -1*(min_flux) + 0.01 # Above zero
        weights = np.min(weights) / weights
    else:
        weights = np.zeros(len(f))
        zeros = f == 0
        weights[~zeros] = min_flux / f[~zeros]
        weights[zeros] = np.min(weights[~zeros]) # They will be ignored
    weights -= np.min(weights)
    max_weight = np.max(weights)
    if max_weight != 0:
        weights = weights / max_weight
    model.fitData(x, y, parinfo=parinfo, weights=weights)
    #model.fitData(x, y, parinfo=parinfo)

    # TODO: Remove
    if False:
        #def func(x, a, mu, sig):
            #return a*np.exp(-(x-mu)**2/(2*sig**2))

        #def show():
            #plt.plot(spectrum_slice['waveobs'], spectrum_slice['flux'])
            #plt.plot(x, baseline+func(x, model.A(), model.mu(), model.sig()))
            #plt.show()
        #from scipy.signal import find_peaks_cwt
        #peaks = find_peaks_cwt(y, np.arange(1, 10))
        #print "*", len(peaks)
        #print peaks
        #show()

        from scipy.optimize import curve_fit
        # Let's create a function to model and create data
        def func(x, a, mu, sig):
            return a*np.exp(-(x-mu)**2/(2*sig**2))
            #return ((a*1.)/np.sqrt(2*np.pi*sig**2))*np.exp(-(x-mu)**2/(2*sig**2))

        # Executing curve_fit on noisy data
        func2 = lambda x, a, sigma: func(x, a, mu, sigma)
        popt, pcov = curve_fit(func2, x, y-baseline, p0 = [-1, 0.05])
        #popt returns the best fit values for parameters of the given model (func)
        #print popt

        def show():
            plt.plot(spectrum_slice['waveobs'], spectrum_slice['flux'])
            plt.plot(x, baseline+func2(x, popt[0], popt[1]))
            plt.plot(x[:-1], baseline+func(x[:-1], model.A(), model.mu(), model.sig()))
            plt.show()
        show()


    return model

def __fit_voigt(spectrum_slice, continuum_model, mu, sig=None, A=None, gamma=None, baseline_margin=0, free_mu=False):
    """
    Fits a voigt at a given wavelength location using a fitted continuum model.

    - For absorption lines, it will alway be true: model.A() < 0 and model.sig() > 0.
    - For absorption lines, model.gamma() < 0 indicates strange wings and probably a bad fit.
    - A mu parameter (model.mu) outside the region used for the fitting it is also a symptom of bad fit.
    - If baseline_margin is different from 0, the continuum will be adjusted by letting
      free the baseline in between [(1.-baseline_margin)*baseline, (1+baseline_margin)*baseline]
    """
    model = VoigtModel()
    x = spectrum_slice['waveobs']
    y = spectrum_slice['flux']
    min_flux = np.min(y)

    # Parameters estimators
    baseline = np.median(continuum_model(spectrum_slice['waveobs']))
    if A is None:
        A = np.min((min_flux - baseline, -1e-10))
    if sig is None:
        sig = np.max(((x[-1] - x[0])/3.0, 1e-10))
    if gamma is None:
        gamma = np.max(((x[-1] - x[0])/2.0, 1e-10))
    if free_mu and (mu < x[0] or mu > x[-1]):
        # Correct mu in case it is outside limits
        mu = spectrum_slice['waveobs'][np.argmin(spectrum_slice['flux'])]


    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(5)]
    parinfo[0]['value'] = baseline # Continuum
    if baseline_margin == 0:
        parinfo[0]['fixed'] = True
    else:
        parinfo[0]['fixed'] = False
        parinfo[0]['limited'] = [True, True]
        #parinfo[0]['limits'] = [0.95*baseline, 1.05*baseline]
        parinfo[0]['limits'] = [(1.-baseline_margin)*baseline, (1+baseline_margin)*baseline]
    parinfo[1]['value'] = A # Only negative (absorption lines) and greater than the lowest point + 25%
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [(min_flux-baseline) * 1.25, -1e-10]
    parinfo[2]['value'] = sig # Only positives (absorption lines) and lower than the spectrum slice
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [1e-10, x[-1] - x[0]]
    parinfo[3]['value'] = mu # Peak only within the spectrum slice
    if not free_mu:
        parinfo[3]['fixed'] = True
    else:
        parinfo[3]['fixed'] = False
        parinfo[3]['limited'] = [True, True]
        parinfo[3]['limits'] = [x[0], x[-1]]
        #parinfo[3]['limits'] = [mu - 0.001, mu + 0.001]
    parinfo[4]['value'] = gamma # Only positives (not zero, otherwise its a gaussian) and small (for nm, it should be <= 0.01 aprox but I leave it in relative terms considering the spectrum slice)
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [1e-10, x[-1] - x[0]]

    # If there are only 4 data point, fix 'mu'
    # - if not, the fit will fail with an exception because there are not
    #   more data points than parameters
    if len(spectrum_slice) == 4:
        parinfo[3]['fixed'] = True

    model.fitData(x, y, parinfo=parinfo)

    return model


def __fit_line(spectrum_slice, continuum_model, mu, sig=None, A=None, gamma=None, discard_gaussian = False, discard_voigt = False, baseline_margin=0, free_mu=False):
    """
    Fits a gaussian and a voigt at a given wavelength location using a fitted continuum model.

    - It selects the best fitted model (gaussian or voigt) unless one of them is disabled by the discard_gaussian or discard_voigt argument.
    - For absorption lines, it will alway be true: model.A() < 0 and model.sig() > 0.
    - For absorption lines fitted with voigt, model.gamma() < 0 indicates strange wings and probably a bad fit
    - A mu parameter (model.mu) outside the region used for the fitting it is also a symptom of bad fit.
    - If baseline_margin is different from 0, the continuum will be adjusted by letting
      free the baseline in between [(1.-baseline_margin)*baseline, (1+baseline_margin)*baseline]
    """
    if not discard_gaussian:
        # Default values for failed fit:
        rms_gaussian = 9999.0
        discard_gaussian = True

        # If there are more data points than parameters (if not, the fit will fail with an exception)
        # - 2 free parameters: A, sig
        # - 1 fix parameter: mu (but it will set to free if there is enough data)
        if len(spectrum_slice) > 2:
            try:
                gaussian_model = __fit_gaussian(spectrum_slice, continuum_model, mu, sig=sig, A=A, baseline_margin=baseline_margin, free_mu=free_mu)

                residuals = gaussian_model.residuals()
                rms_gaussian = np.sqrt(np.sum(np.power(residuals, 2)) / len(residuals))
                #minimization_value = gaussian_model.m.fnorm
                #degrees_of_freedom = gaussian_model.m.dof,
                #chisq = np.sum(np.power(residuals, 2) * gaussian_model.weights)
                discard_gaussian = False
            except Exception as e:
                pass
                #if len(e.message) > 0:
                    #print e.message

    if not discard_voigt:
        # Default values for failed fit:
        rms_voigt = 9999.0
        discard_voigt = True

        # If there are more data points than parameters (if not, the fit will fail with an excepteion)
        # - 3 free parameters: A, sig, gamma
        # - 1 fix parameter: mu (but it will set to free if there is enough data)
        if len(spectrum_slice) > 3:
            try:
                voigt_model = __fit_voigt(spectrum_slice, continuum_model, mu, sig=sig, A=A, gamma=gamma, baseline_margin=baseline_margin, free_mu=free_mu)
                residuals = voigt_model.residuals()
                rms_voigt = np.sqrt(np.sum(np.power(residuals, 2)) / len(residuals))
                discard_voigt = False
            except Exception as e:
                pass
                #if len(e.message) > 0:
                    #print e.message

    if (not discard_gaussian and not discard_voigt and rms_gaussian <= rms_voigt) or (not discard_gaussian and discard_voigt):
        return gaussian_model, rms_gaussian
    elif (not discard_gaussian and not discard_voigt and rms_gaussian > rms_voigt) or (discard_gaussian and not discard_voigt):
        return voigt_model, rms_voigt
    else:
        raise Exception("Gaussian or Voigt fit for absorption line not possible.")

def __remove_consecutives_features(features):
    """
    Remove features (i.e. peaks or base points) that are consecutive, it makes
    no sense to have two peaks or two base points together.
    """
    if len(features) >= 2:
        duplicated_features = (np.abs(features[:-1] - features[1:]) == 1)
        duplicated_features = np.array([False] + duplicated_features.tolist())
        cleaned_features = features[~duplicated_features]
        return cleaned_features
    else:
        return features


def __detect_false_and_noisy_features(spectrum, linemasks):
    """
    Detect features that are false positive or noise:

    - Peaks that are less deep than it nearby base points (false positives)
    - Peaks too close (2 or less positions) to the next and previous base point (noise)
    """
    peaks = linemasks['peak']
    base_points = np.asarray([linemasks['base'][0]] +  linemasks['top'].tolist())

    ## First feature found in spectrum: base point
    ## Last feature found in spectrum: base point
    # Left
    peak_base_diff_left = spectrum['flux'][base_points[:-1]] - spectrum['flux'][peaks]
    peak_base_index_diff_left = base_points[:-1] - peaks
    # Right
    peak_base_diff_right = spectrum['flux'][base_points[1:]] - spectrum['flux'][peaks]
    peak_base_index_diff_right = base_points[1:] - peaks

    # Find peaks too close (2 positions) to the next or previous peak (noise)
    # First and last peaks are ignored
    peak_peak_index_diff_left = peaks[:-1] - peaks[1:]
    peak_peak_index_diff_left = np.asarray(peak_peak_index_diff_left.tolist() + [9999.0]) # Last peak
    peak_peak_index_diff_right = peaks[1:] - peaks[:-1]
    peak_peak_index_diff_right = np.asarray([9999.0] + peak_peak_index_diff_right.tolist()) # First peak

    # Filter false positive and noise
    # - Peaks that are less deep than it nearby base points (false positives)
    # - Peaks too close (2 or less positions) to the next and previous base point (noise)
    # - Peaks too close (2 positions) to the next or previous peak (noise)
    false_positives = np.logical_or((peak_base_diff_left < 0), (peak_base_diff_right < 0))
    noise1 = np.logical_and((np.abs(peak_base_index_diff_left) <= 2), (np.abs(peak_base_index_diff_right) <= 2))
    noise2 = np.logical_or((np.abs(peak_peak_index_diff_left) == 2), (np.abs(peak_peak_index_diff_right) == 2))
    noise = np.logical_or(noise1, noise2)

    rejected = np.logical_or(false_positives, noise)
    return rejected

def __assert_structure(xcoord, yvalues, peaks, base_points):
    """
    Given a group of peaks and base_points with the following assumptions:

    - base_points[i] < base_point[i+1]
    - peaks[i] < peaks[i+1]
    - base_points[i] < peaks[i] < base_points[i+1]

    The function returns peaks and base_points where:

    - The first and last feature is a base point: base_points[0] < peaks[0] < base_points[1] < ... < base_points[n-1] < peaks[n-1] < base_points[n] where n = len(base_points)
    - len(base_points) = len(peaks) + 1
    """
    if len(peaks) == 0 or len(base_points) == 0:
        return [], []

    # Limit the base_points array to the ones that are useful, considering that
    # the first and last peak are always removed
    first_wave_peak = xcoord[peaks][0]
    first_wave_base = xcoord[base_points][0]
    if first_wave_peak > first_wave_base:
        if len(base_points) - len(peaks) == 1:
            ## First feature found in spectrum: base point
            ## Last feature found in spectrum: base point
            base_points = base_points
            peaks = peaks
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectrum: base point
            ## Last feature found in spectrum: peak
            # - Remove last peak
            #base_points = base_points
            #peaks = peaks[:-1]
            # - Add a base point (last point in the spectrum)
            base_points = np.hstack((base_points, [len(xcoord)-1]))
            peaks = peaks
        else:
            raise Exception("This should not happen")
    else:
        if len(base_points) - len(peaks) == -1:
            ## First feature found in spectrum: peak
            ## Last feature found in spectrum: peak
            # - Remove first and last peaks
            #base_points = base_points
            #peaks = peaks[1:-1]
            # - Add two base points (first and last point in the spectrum)
            base_points = np.hstack(([0], base_points))
            base_points = np.hstack((base_points, [len(xcoord)-1]))
            peaks = peaks
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectrum: peak
            ## Last feature found in spectrum: base point
            # - Remove first peak
            #base_points = base_points
            #peaks = peaks[1:]
            # - Add a base point (first point in the spectrum)
            base_points = np.hstack(([0], base_points))
            peaks = peaks
        else:
            raise Exception("This should not happen")

    return peaks, base_points

def create_linemasks_structure(num_lines):
    return __create_linemasks_structure(num_lines)

def __create_linemasks_structure(num_peaks):
    """
    Creates a linemasks structure compatible with SPECTRUM, Turbospectrum and
    several iSpec functions
    """
    dtype = _get_atomic_linelist_definition() + __get_fitted_lines_definition()
    linemasks = np.recarray((num_peaks, ), dtype=dtype)
    # Initialization
    linemasks['grouped'] = False # Lines that share the same mask (i.e. HFS and isotopes)
    linemasks['reference_for_group'] = False # It is a representant for its group
    linemasks['discarded'] = False
    # Line mask
    linemasks['wave_peak'] = 0.0
    linemasks['wave_base'] = 0.0
    linemasks['wave_top'] = 0.0
    linemasks['note'] = ""
    # Position in the spectrum, usefull for noise detection
    linemasks['peak'] = 0
    linemasks['base'] = 0
    linemasks['top'] = 0
    linemasks['depth'] = 0.0
    # Relative depth is "peak - mean_base_point" with respect to the total continuum
    linemasks['relative_depth'] = 0
    # Default values for telluric affectation
    linemasks['telluric_wave_peak'] = 0
    linemasks['telluric_fwhm'] = 0
    linemasks['telluric_R'] = 0
    linemasks['telluric_depth'] = 0
    # Default values that indicate that the line has not been fitted
    linemasks['base_fit'] = 0
    linemasks['top_fit'] = 0
    linemasks['wave_base_fit'] = 0.0
    linemasks['wave_top_fit'] = 0.0
    linemasks['mu'] = 0.0
    linemasks['mu_err'] = 0.0
    linemasks['sig'] = 0.0
    linemasks['A'] = 0.0
    linemasks['baseline'] = 0.0
    linemasks['gamma'] = 0.0
    linemasks['fwhm'] = 0.0
    linemasks['fwhm_kms'] = 0.0
    linemasks['R'] = 0.0 # Resolving power
    linemasks['depth_fit'] = 0.0
    linemasks['relative_depth_fit'] = 0.0
    linemasks['integrated_flux'] = 0.0
    linemasks['ewr'] = 0.0 # Reduced equivalent width
    linemasks['ew'] = 0.0 # Equivalent width
    linemasks['ew_err'] = 0.0
    linemasks['snr'] = 0.0 # SNR of the region used for EW calculation
    linemasks['mean_flux'] = 0.0
    linemasks['mean_flux_continuum'] = 0.0
    linemasks['diff_wavelength'] = 0.0
    linemasks['rms'] = 9999.0

    # Atomic
    linemasks['element'] = ""
    #linemasks['turbospectrum_element'] = ""
    #for key in ['wave_A', 'wave_nm', 'loggf', 'lower_state_eV', 'lower_state_cm1', 'lower_j', 'lower_g']:
    for key in ['wave_A', 'wave_nm', 'loggf', 'lower_state_eV', 'lower_state_cm1', 'lower_j']:
        linemasks[key] = 0.
    for key in ['upper_state_eV', 'upper_state_cm1', 'upper_j', 'upper_g']:
        linemasks[key] = 0.
    #for key in ['lande_lower', 'lande_upper', 'lande_mean']:
    for key in ['lande_lower', 'lande_upper']:
        linemasks[key] = 0.
    linemasks['spectrum_transition_type'] = "99"
    linemasks['turbospectrum_rad'] = 1e5
    for key in ['rad', 'stark', 'waals', 'waals_single_gamma_format', 'turbospectrum_fdamp', 'spectrum_fudge_factor']:
        linemasks[key] = 0.
    for key in ['theoretical_depth', 'theoretical_ew']:
        linemasks[key] = 0.
    linemasks['lower_orbital_type'] = "X"
    linemasks['upper_orbital_type'] = "X"
    #for key in ['lower_coupling', 'lower_shrinked_designation', 'upper_coupling', 'upper_shrinked_designation', 'designation']:
    #for key in ['lower_coupling', 'lower_shrinked_designation', 'upper_coupling', 'upper_shrinked_designation']:
        #linemasks[key] = ""
    linemasks['molecule'] = "False"
    #for key in ['lower_coupling', 'lower_shrinked_designation', 'upper_coupling', 'upper_shrinked_designation', 'designation']:
    #for key in ['lower_coupling', 'lower_shrinked_designation', 'upper_coupling', 'upper_shrinked_designation']:
        #linemasks[key] = ""
    #linemasks['atomic_number_001'] = 0
    #linemasks['atomic_number_002'] = 0
    #linemasks['isotope_001'] = 0
    #linemasks['isotope_002'] = 0
    linemasks['spectrum_synthe_isotope'] = 0
    linemasks['ion'] = 0
    linemasks['spectrum_moog_species'] = ""
    linemasks['turbospectrum_species'] = ""
    linemasks['width_species'] = ""
    linemasks['reference_code'] = ""
    #linemasks['reference'] = ""
    for key in ['spectrum_support', 'turbospectrum_support', 'moog_support', 'width_support', 'synthe_support', 'sme_support']:
        linemasks[key] = "False"
    return linemasks


def read_atomic_linelist(linelist_filename, wave_base=None, wave_top=None):
    """
    Read atomic linelist.

    The linelist can be filtered by wave_base and wave_top (in nm) to reduce
    memory consumption.
    """
    atomic_dtype = _get_atomic_linelist_definition()
    fitted_dtype = __get_fitted_lines_definition()
    try:
        num_cols = len(open(linelist_filename, "r").readline().split("\t"))
        if num_cols == len(atomic_dtype):
            # Only atomic linelist
            linelist = np.genfromtxt(linelist_filename, dtype=atomic_dtype, delimiter="\t", skip_header=1, loose=False)
        #elif num_cols == len(atomic_dtype) + len(fitted_dtype):
            ## Atomic + fitted information
            #linelist = np.genfromtxt(linelist_filename, dtype=atomic_dtype+fitted_dtype, delimiter="\t", skip_header=1, loose=False)
        else:
            raise Exception()
    except:
        raise Exception("Wrong atomic linelist file format!")
    wfilter = None
    if wave_base is not None and wave_top is not None:
        wfilter = np.logical_and(linelist['wave_nm']  >= wave_base, linelist['wave_nm']  <= wave_top)
    elif wave_base is not None:
        wfilter = linelist['wave_nm']  >= wave_base
    elif wave_top is not None:
        wfilter = linelist['wave_nm']  <= wave_top

    if wfilter is not None:
        linelist = linelist[wfilter]
    if len (linelist) == 0:
        logging.warn("Linelist does not contain any line between %.4f and %.4f nm" % (wave_base, wave_top))

    return linelist

def __get_fitted_lines_definition():
    return [
            ('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100'), \
            ('peak', int), ('base', int), ('top', int), \
            ('depth', float), ('relative_depth', float), \
            ('wave_base_fit', float), ('wave_top_fit', float), \
            ('base_fit', int), ('top_fit', int), \
            ('mu', float), ('sig', float), ('A', float), ('baseline', float), ('gamma', float), \
            ('mu_err', float), \
            ('fwhm', float), ('fwhm_kms', float), ('R', float), \
            ('depth_fit', float), ('relative_depth_fit', float), \
            ('integrated_flux', float), ('ewr', float), ('ew', float), ('ew_err', float), \
            ('snr', float), ('mean_flux', float), ('mean_flux_continuum', float), ('diff_wavelength', float), \
            ('rms', float), \
            ('telluric_wave_peak', float), ('telluric_fwhm', float), ('telluric_R', float), ('telluric_depth', float), \
            ('grouped', '|S5'), ('reference_for_group', '|S5'),
            ('discarded', '|S5'), \
            ]

def _get_atomic_linelist_definition():
    return [('element', '|S4'),
     #('turbospectrum_element', '|S6'),
     ('wave_A', '<f8'),
     ('wave_nm', '<f8'),
     ('loggf', '<f8'),
     ('lower_state_eV', '<f8'),
     ('lower_state_cm1', '<f8'),
     ('lower_j', '<f8'),
     #('lower_g', '<f8'),
     ('upper_state_eV', '<f8'),
     ('upper_state_cm1', '<f8'),
     ('upper_j', '<f8'),
     ('upper_g', '<f8'),
     ('lande_lower', '<f8'),
     ('lande_upper', '<f8'),
     #('lande_mean', '<f8'),
     ('spectrum_transition_type', '|S2'),
     ('turbospectrum_rad', '<f8'),
     ('rad', '<f8'),
     ('stark', '<f8'),
     ('waals', '<f8'),
     ('waals_single_gamma_format', '<f8'),
     ('turbospectrum_fdamp', '<f8'),
     ('spectrum_fudge_factor', '<f8'),
     ('theoretical_depth', '<f8'),
     ('theoretical_ew', '<f8'),
     ('lower_orbital_type', '|S1'),
     ('upper_orbital_type', '|S1'),
     #('lower_coupling', '|S2'),
     #('lower_shrinked_designation', '|S41'),
     #('upper_coupling', '|S2'),
     #('upper_shrinked_designation', '|S45'),
     #('designation', '|S93'),
     ('molecule', '|S1'),
     #('molecule', '|S5'),
     #('atomic_number_001', '<i4'),
     #('atomic_number_002', '<f8'),
     #('isotope_001', '<f8'),
     #('isotope_002', '<f8'),
     ('spectrum_synthe_isotope', '<i4'),
     ('ion', '<i4'),
     ('spectrum_moog_species', '|S6'),
     ('turbospectrum_species', '|S14'),
     ('width_species', '|S6'),
     ('reference_code', '|S10'),
     #('reference', '|S39'),
     ('spectrum_support', '|S1'),
     ('turbospectrum_support', '|S1'),
     ('moog_support', '|S1'),
     ('width_support', '|S1'),
     ('synthe_support', '|S1'),
     ('sme_support', '|S1')]
     #('spectrum_support', '|S5'),
     #('turbospectrum_support', '|S5'),
     #('moog_support', '|S5'),
     #('width_support', '|S5'),
     #('synthe_support', '|S5')]

def write_atomic_linelist(linelist, linelist_filename=None, code=None, tmp_dir=None):
    """
    Write atomic linelist.
    If code is specified ('spectrum', 'turbospectrum', 'moog'), then it is saved
    in the file format compatible with the indicated code (thus not all the atomic
    information will be stored).
    """
    if code is not None:
        code = code.lower()
        if code not in ['spectrum', 'turbospectrum', 'moog', 'moog_barklem', 'synthe']:
            raise Exception("Unknown radiative transfer code: %s" % (code))

        if code == "moog":
            #logging.info("MOOG file format")
            return __moog_write_atomic_linelist(linelist, linelist_filename=linelist_filename, tmp_dir=tmp_dir)
        elif code == "moog_barklem":
            #logging.info("MOOG Barklem.dat file format")
            return __moog_barklem_write_atomic_linelist(linelist, linelist_filename=linelist_filename, tmp_dir=tmp_dir)
        elif code == "turbospectrum":
            #logging.info("Turbospectrum file format")
            return __turbospectrum_write_atomic_linelist(linelist, linelist_filename=linelist_filename, tmp_dir=tmp_dir)
        elif code == "synthe":
            #logging.info("Synthe file format")
            return __synthe_write_atomic_linelist(linelist, linelist_filename=linelist_filename, tmp_dir=tmp_dir)
        else:
            #logging.info("SPECTRUM file format")
            return __spectrum_write_atomic_linelist(linelist, linelist_filename=linelist_filename, tmp_dir=tmp_dir)
    else:
        if linelist_filename is not None:
            out = open(linelist_filename, "w")
        else:
            # Temporary file
            out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
        __generic_write_atomic_linelist_header(out, include_fit=False)
        for line in linelist:
            __generic_write_atomic_linelist_element(out, line, include_fit=False)
        out.close()
        return out.name

def __generic_write_atomic_linelist_header(out, include_fit):
    for i, (key, dtype) in enumerate(_get_atomic_linelist_definition()):
        if i == 0:
            out.write("%s" % (key))
        else:
            out.write("\t%s" % (key))
    if include_fit:
        for key, dtype in __get_fitted_lines_definition():
            out.write("\t%s" % (key))
    out.write("\n")

def __generic_write_atomic_linelist_element(out, data, include_fit):
    out.write("%s\t" % (data['element']))
    #out.write("%s\t" % (data['turbospectrum_element']))
    out.write("%.3f\t" % (data['wave_A']))
    out.write("%.4f\t" % (data['wave_nm']))
    out.write("%.3f\t" % (data['loggf']))
    out.write("%.4f\t" % (data['lower_state_eV']))
    out.write("%.3f\t" % (data['lower_state_cm1']))
    out.write("%.1f\t" % (data['lower_j']))
    #out.write("%.1f\t" % (data['lower_g']))
    out.write("%.4f\t" % (data['upper_state_eV']))
    out.write("%.3f\t" % (data['upper_state_cm1']))
    out.write("%.1f\t" % (data['upper_j']))
    out.write("%.1f\t" % (data['upper_g']))
    out.write("%.3f\t" % (data['lande_lower']))
    out.write("%.3f\t" % (data['lande_upper']))
    #out.write("%.3f\t" % (data['lande_mean']))
    out.write("%s\t" % (data['spectrum_transition_type']))
    out.write("%.2E\t" % (data['turbospectrum_rad']))
    out.write("%.3f\t" % (data['rad']))
    out.write("%.3f\t" % (data['stark']))
    out.write("%.3f\t" % (data['waals']))
    out.write("%.3f\t" % (data['waals_single_gamma_format']))
    out.write("%.3f\t" % (data['turbospectrum_fdamp']))
    out.write("%.3f\t" % (data['spectrum_fudge_factor']))
    out.write("%.3f\t" % (data['theoretical_depth']))
    out.write("%.3f\t" % (data['theoretical_ew']))
    out.write("%s\t" % (data['lower_orbital_type']))
    out.write("%s\t" % (data['upper_orbital_type']))
    #out.write("%s\t" % (data['lower_coupling']))
    #out.write("%s\t" % (data['lower_shrinked_designation']))
    #out.write("%s\t" % (data['upper_coupling']))
    #out.write("%s\t" % (data['upper_shrinked_designation']))
    #out.write("%s\t" % (data['designation']))
    out.write("%s\t" % (str(data['molecule'])[0]))
    #out.write("%s\t" % (str(data['molecule'])))
    #out.write("%i\t" % (data['atomic_number_001']))
    #out.write("%.0f\t" % (data['atomic_number_002'])) # Zero decimal float instead of integer because it can be a float('nan')
    #out.write("%.0f\t" % (data['isotope_001'])) # Zero decimal float instead of integer because it can be a float('nan')
    #out.write("%.0f\t" % (data['isotope_002'])) # Zero decimal float instead of integer because it can be a float('nan')
    out.write("%i\t" % (data['spectrum_synthe_isotope']))
    out.write("%i\t" % (data['ion']))
    out.write("%s\t" % (data['spectrum_moog_species']))
    out.write("%s\t" % (data['turbospectrum_species']))
    out.write("%s\t" % (data['width_species']))
    out.write("%s\t" % (data['reference_code']))
    #out.write("%s\t" % (data['reference']))
    out.write("%s\t" % (str(data['spectrum_support'])[0])) # Only 'T' or 'F' to save memory
    out.write("%s\t" % (str(data['turbospectrum_support'])[0]))
    out.write("%s\t" % (str(data['moog_support'])[0]))
    out.write("%s\t" % (str(data['width_support'])[0]))
    out.write("%s\t" % (str(data['synthe_support'])[0]))
    out.write("%s" % (str(data['sme_support'])[0]))
    if include_fit:
        out.write("\t%.4f\t" % (data['wave_peak']))
        out.write("%.4f\t" % (data['wave_base']))
        out.write("%.4f\t" % (data['wave_top']))
        out.write("%s\t" % (data['note']))
        out.write("%i\t" % (data['peak']))
        out.write("%i\t" % (data['base']))
        out.write("%i\t" % (data['top']))
        out.write("%.3f\t" % (data['depth']))
        out.write("%.3f\t" % (data['relative_depth']))
        out.write("%.4f\t" % (data['wave_base_fit']))
        out.write("%.4f\t" % (data['wave_top_fit']))
        out.write("%i\t" % (data['base_fit']))
        out.write("%i\t" % (data['top_fit']))
        out.write("%.4f\t" % (data['mu']))
        out.write("%.4f\t" % (data['sig']))
        out.write("%.4f\t" % (data['A']))
        out.write("%.3f\t" % (data['baseline']))
        out.write("%.4f\t" % (data['gamma']))
        out.write("%.4f\t" % (data['mu_err']))
        out.write("%.4f\t" % (data['fwhm']))
        out.write("%.4f\t" % (data['fwhm_kms']))
        out.write("%.0f\t" % (data['R']))
        out.write("%.3f\t" % (data['depth_fit']))
        out.write("%.3f\t" % (data['relative_depth_fit']))
        out.write("%.4f\t" % (data['integrated_flux']))
        out.write("%.2f\t" % (data['ewr']))
        out.write("%.1f\t" % (data['ew']))
        out.write("%.1f\t" % (data['ew_err']))
        out.write("%.1f\t" % (data['snr']))
        out.write("%.3f\t" % (data['mean_flux']))
        out.write("%.3f\t" % (data['mean_flux_continuum']))
        out.write("%.4f\t" % (data['diff_wavelength']))
        out.write("%.3f\t" % (data['rms']))
        out.write("%.4f\t" % (data['telluric_wave_peak']))
        out.write("%.4f\t" % (data['telluric_fwhm']))
        out.write("%.0f\t" % (data['telluric_R']))
        out.write("%.3f\t" % (data['telluric_depth']))
        out.write("%s\t" % (str(data['grouped'])))
        out.write("%s\t" % (str(data['reference_for_group'])))
        out.write("%s" % (str(data['discarded'])))
    out.write("\n")


def __spectrum_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None):
    """
    Saves a SPECTRUM linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.
    """
    supported = linelist['spectrum_support'] == "T"
    linelist = linelist[supported]
    linelist = linelist.copy()
    # http://www.appstate.edu/~grayro/spectrum/spectrum276/node14.html
    #  Since only the energy of the lower state is used in molecular calculations,
    #  this entry (upper state) is sometimes used to encode the molecular band information
    molecules = linelist['molecule'] == "T"
    linelist['upper_state_cm1'][molecules] = 0.

    if linelist_filename is not None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
    #4750.196  26.0 0  36078  57130  -3.662  1.0  GA  8.09  -4.61  -7.32  Fe_1
    out.write("\n".join(["  ".join(map(str, (line['wave_A'], line['spectrum_moog_species'], line['spectrum_synthe_isotope'], line['lower_state_cm1'], line['upper_state_cm1'], line['loggf'], line['spectrum_fudge_factor'], line['spectrum_transition_type'], line['rad'], line['stark'], line['waals'], "_".join(line['element'].split())))) for line in linelist]))
    out.close()
    return out.name


def find_linemasks(spectrum, continuum_model, atomic_linelist=None, max_atomic_wave_diff=0.0005, telluric_linelist=None, vel_telluric=0.0, minimum_depth=None, maximum_depth=None, discard_gaussian = False, discard_voigt = False, check_derivatives=False, smoothed_spectrum=None, accepted_for_fitting=None, consider_omara=False, frame=None):
    """
    Generate a line masks for a spectrum by finding peaks and base points.

    It is recommended that the spectrum is normalized (better results are obtained).

    It is recommended to provide a smoothed version of the spectra, so that the
    line edges will be better determined.

    It tries to fit a gaussian and a voigt model and selects the best unless one of them is disabled by
    the discard_gaussian or discard_voigt argument (if both of them are disabled, there will be
    no fit information).

    If 'check_derivates' is True, the fit will be limited not to the region limits
    but to a refined smaller limits calculated from the second and third derivative
    (convex+concave regions) from 'smoothed_spectrum' if present or 'spectrum' if not.

    To save computation time, min and max depth can be indicated and all the lines out
    of this range will not be considered for fit process (save CPU time) although
    the rest of the information of the line will be conserved in the output.

    Additionally, features that are false positive or noise are not fitted. This
    implies ignoring:

    - Peaks that are less deep than it nearby base points (false positives)
    - Peaks too close (2 or less positions) to the next and previous base point (noise)

    Returns a complete structure with all the necessary information to
    determine if it is a line of interest.
    """
    logging.info("Finding peaks and base points...")
    # Peaks and base points will agree with the following criteria:
    # - The first and last feature is a base point
    #    base_points[0] < peaks[0] < base_points[1] < ... < base_points[n-1] < peaks[n-1] < base_points[n]
    #    where n = len(base_points)
    # - len(base_points) = len(peaks) + 1
    if smoothed_spectrum is None:
        peaks, base_points = __find_peaks_and_base_points(spectrum['waveobs'], spectrum['flux'])
    else:
        peaks, base_points = __find_peaks_and_base_points(smoothed_spectrum['waveobs'], smoothed_spectrum['flux'])
    # If no peaks found, just finnish
    if len(peaks) == 0 or len(base_points) == 0:
        return None

    # Depth of the peak with respect to the total continuum in % over the total continuum
    # - In case that the peak is higher than the continuum, depth < 0
    continuum_at_peak = continuum_model(spectrum['waveobs'][peaks])
    flux_at_peak = spectrum['flux'][peaks]
    depth = 1 - (flux_at_peak /continuum_at_peak)
    dfilter = depth < 0.0
    depth[dfilter] = 0.0

    if minimum_depth is None:
        minimum_depth = np.max((0., np.min(depth)))
    if maximum_depth is None:
        maximum_depth = np.min((1., np.max(depth)))

    # To save computation time, min and max depth can be indicated and all the lines out
    # of this range will not be considered for fit process (save CPU time) although
    # the rest of the information of the line will be conserved in the output
    accepted_for_fitting = np.logical_and(depth >= minimum_depth, depth <= maximum_depth)

    last_reported_progress = -1
    if frame is not None:
        frame.update_progress(0)

    num_peaks = len(peaks)
    linemasks = __create_linemasks_structure(num_peaks)

    # Line mask
    linemasks['wave_peak'] = spectrum['waveobs'][peaks]
    linemasks['wave_base'] = spectrum['waveobs'][base_points[:-1]]
    linemasks['wave_top'] = spectrum['waveobs'][base_points[1:]]
    # Position in the spectrum, usefull for noise detection
    linemasks['peak'] = peaks
    linemasks['base'] = base_points[:-1]
    linemasks['top'] = base_points[1:]
    linemasks['depth'] = depth
    # Relative depth is "peak - mean_base_point" with respect to the total continuum
    # - In case that the mean base point is higher than the continuum, relative_depth < 0
    # - relative_depth < depth is always true
    flux_from_top_base_point_to_continuum = np.abs(continuum_at_peak - np.mean([spectrum['flux'][base_points[:-1]], spectrum['flux'][base_points[1:]]]))
    linemasks['relative_depth'] = ((continuum_at_peak - (flux_at_peak + flux_from_top_base_point_to_continuum)) / continuum_at_peak)

    # To save computation time, exclude false positives and noise from the fitting process
    rejected_by_noise = __detect_false_and_noisy_features(spectrum, linemasks)
    accepted_for_fitting = np.logical_and(accepted_for_fitting, np.logical_not(rejected_by_noise))

    linemasks = fit_lines(linemasks, spectrum, continuum_model, \
                                atomic_linelist = atomic_linelist, \
                                max_atomic_wave_diff = max_atomic_wave_diff, \
                                telluric_linelist = telluric_linelist, \
                                vel_telluric = vel_telluric, \
                                discard_gaussian = discard_gaussian, \
                                discard_voigt = discard_voigt, \
                                smoothed_spectrum = smoothed_spectrum, \
                                check_derivatives = check_derivatives, \
                                free_mu=True, crossmatch_with_mu=True, \
                                accepted_for_fitting=accepted_for_fitting)

    # Identify peaks higher than continuum
    # - Depth is negative if the peak is higher than the continuum
    # - Relative depth is negative if the mean base point is higher than the continuum
    rejected_by_depth_higher_than_continuum = (linemasks['depth'] < 0)

    # Identify peaks with a depth inferior/superior to a given limit (% of the continuum)
    # - Observed depth
    rejected_by_depth_limits1 = np.logical_or((linemasks['depth'] <= minimum_depth), (linemasks['depth'] >= maximum_depth))
    # - Fitted depth
    rejected_by_depth_limits2 = np.logical_or((linemasks['depth_fit'] <= minimum_depth), (linemasks['depth_fit'] >= maximum_depth))
    rejected_by_depth_limits = np.logical_or(rejected_by_depth_limits1, rejected_by_depth_limits2)

    # Identify bad fits (9999.0: Cases where has not been possible to fit a gaussian/voigt)
    rejected_by_bad_fit = (linemasks['rms'] >= 9999.0)

    discarded = rejected_by_depth_higher_than_continuum
    discarded = np.logical_or(discarded, rejected_by_depth_limits)
    discarded = np.logical_or(discarded, rejected_by_bad_fit)

    linemasks = linemasks[~discarded]

    # Correct potential overlapping between line masks
    linemasks.sort(order=['wave_peak'])
    for i, line in enumerate(linemasks[:-1]):
        if line['wave_top'] > linemasks['wave_base'][i+1]:
            mean = (line['wave_top'] + linemasks['wave_base'][i+1]) / 2.0
            line['wave_top'] = mean
            linemasks['wave_base'][i+1] = mean

    return linemasks

def __find_wave_indices(regions, spectrum):
    """
    Find index in spectrum for base, top and peak
    """
    regions['base'] = 0
    regions['top'] = 0
    regions['peak'] = 0
    for i in np.arange(len(regions)):
        where_base = np.where(spectrum['waveobs'] >= regions['wave_base'][i])
        if len(where_base[0]) > 0:
            regions['base'][i] = where_base[0][0]
        where_top = np.where(spectrum['waveobs'] >= regions['wave_top'][i])
        if len(where_top[0]) > 0:
            regions['top'][i] = where_top[0][0]
        where_peak = np.where(spectrum['waveobs'] >= regions['wave_peak'][i])
        if len(where_peak[0]) > 0:
            regions['peak'][i] = where_peak[0][0]
    return regions

def __calculate_depths(regions, spectrum, continuum_model):
    """
    Calculate depth and relative depth for regions where 'peak', 'base' and 'top'
    is already indicated.
    """
    # Depth of the peak with respect to the total continuum in % over the total continuum
    # - In case that the peak is higher than the continuum, depth < 0
    peaks = regions['peak']
    top = regions['top']
    base = regions['base']
    continuum_at_peak = continuum_model(spectrum['waveobs'][peaks])
    flux_at_peak = spectrum['flux'][peaks]
    depth = 1 - (flux_at_peak /continuum_at_peak)
    dfilter = depth < 0.0
    depth[dfilter] = 0.0
    regions['depth'] = depth
    # Relative depth is "peak - mean_base_point" with respect to the total continuum
    # - In case that the mean base point is higher than the continuum, relative_depth < 0
    # - relative_depth < depth is always true
    flux_from_top_base_point_to_continuum = np.abs(continuum_at_peak - np.mean((spectrum['flux'][base], spectrum['flux'][top])))
    regions['relative_depth'] = ((continuum_at_peak - (flux_at_peak + flux_from_top_base_point_to_continuum)) / continuum_at_peak)

    return regions

def reset_fitted_data_fields(regions):
    """
    Set to the default value (e.g. zeros) the columns corresponding to fitted
    lines (e.g. mu, sigma, equivalent width) for regions coming from
    find_linemasks or fit_lines
    """
    regions = regions.copy()
    total_regions = len(regions)
    # Regions coming from a previous find_linemasks or fit_lines
    # Clean fitted lines data
    zeroed_regions = __create_linemasks_structure(total_regions)
    fitted_lines_dtype = __get_fitted_lines_definition()
    for key, vtype in fitted_lines_dtype:
        if key in ('wave_base', 'wave_top', 'wave_peak'):
            continue
        regions[key] = zeroed_regions[key]
    return regions

def fit_lines(regions, spectrum, continuum_model, atomic_linelist, max_atomic_wave_diff=0.0005, telluric_linelist=None, vel_telluric=None, discard_gaussian = False, discard_voigt = False, check_derivatives=False, smoothed_spectrum=None, accepted_for_fitting=None, continuum_adjustment_margin=0.0, free_mu=False, crossmatch_with_mu=False, closest_match=False, frame=None):
    """
    Fits gaussians/voigt models in the specified line regions.
    * 'regions' should be an array with 'wave_base', 'wave_peak' and 'wave_top' columns, but it can be also the result of a previous fit_lines or find_linemasks.
    * If 'check_derivates' is True, the fit will be limited not to the region limits but to a refined smaller limits calculated from the second and third derivative (convex+concave regions) from 'smoothed_spectrum' if present or 'spectrum' if not.
    * If 'accepted_for_fitting' array is present, only those regions that are set to true
    will be fitted

    If atomic_linelist is specified, the
    lines will be cross-matched with the atomic information. In that case, the spectrum
    MUST be already radial velocity corrected.
        * The crossmatch will use 'wave_peak' by default unless crossmatch_with_mu is True
          then the peak of the fitted gaussian will be used (which might be different
          to weak_peak if free_mu was also set to True)
        * When closest_match is False, the atomic line is chosen checking the theoretical
          depth and equivalent width (if present).

          NOTE: If we want to crossmatch with lines that we are completely sure
          that exists in the line regions and linelist with the same wavelength,
          then the ideal is to set crossmatch_with_mu=False and closest_match=True

    If telluric_linelist is specified, then it will try to identify which lines
    could be potentially affected by tellurics. It is mandatory to specify the
    velocity relative to the telluric lines.

    If the spectrum has reported errors, the EW error will be calculated (it needs
    the SNR that will be estimated snr = fluxes/errors)

    NOTE: It is recommended to use an already normalized spectrum

    - If baseline_margin is different from 0, the continuum will be adjusted by letting
      free the baseline in between [(1.-baseline_margin)*baseline, (1+baseline_margin)*baseline]

    :returns:
        Array with additional columns such as 'mu', 'sig', 'A', 'baseline'...
    """
    last_reported_progress = -1
    total_regions = len(regions)
    #normalized_spectrum = ispec.create_spectrum_structure(spectrum['waveobs'], spectrum['flux'], spectrum['err'])
    #normalized_spectrum['flux'] /= continuum_model(spectrum['waveobs'])
    #normalized_spectrum['err'] /= continuum_model(spectrum['waveobs'])
    #spectrum = normalized_spectrum
    #continuum_model = ispec.fit_continuum(self.active_spectrum.data, fixed_value=1.0, model="Fixed value")

    logging.info("Fitting line models...")

    atomic_data_dtype = _get_atomic_linelist_definition()
    fitted_lines_dtype = __get_fitted_lines_definition()
    if regions.dtype == atomic_data_dtype+fitted_lines_dtype:
        regions = reset_fitted_data_fields(regions)
    elif regions.dtype.fields.has_key('wave_base') \
            and regions.dtype.fields.has_key('wave_top') \
            and regions.dtype.fields.has_key('wave_peak'):
        # If it is not a complete regions (it has not been created with
        # __create_linemasks_structure and it only contains wave base, top and peak)
        # we create a complete one
        regions_tmp = regions
        regions = __create_linemasks_structure(total_regions)
        regions['wave_base'] = regions_tmp['wave_base']
        regions['wave_top'] = regions_tmp['wave_top']
        regions['wave_peak'] = regions_tmp['wave_peak']
    else:
        raise Exception("Wrong fields in regions.")
    regions = __find_wave_indices(regions, spectrum)
    regions = __calculate_depths(regions, spectrum, continuum_model)

    i = 0
    wfilter = np.logical_and(spectrum['flux'] > 0., spectrum['err'] > 0.)
    if len(np.where(wfilter)[0]) > 0:
        snr = np.median(spectrum['flux'][wfilter] / spectrum['err'][wfilter])
    else:
        snr = None
    # Model: fit gaussian/voigt
    for i in np.arange(total_regions):
        fitting_not_possible = False
        if accepted_for_fitting is None or accepted_for_fitting[i]:
            if check_derivatives:
                # Adjust edges
                if smoothed_spectrum is not None:
                    new_base, new_top = __improve_linemask_edges(smoothed_spectrum['waveobs'], smoothed_spectrum['flux'], regions['base'][i], regions['top'][i], regions['peak'][i])
                else:
                    new_base, new_top = __improve_linemask_edges(spectrum['waveobs'], spectrum['flux'], regions['base'][i], regions['top'][i], regions['peak'][i])
            else:
                new_base = regions['base'][i]
                new_top = regions['top'][i]
            # Actually, modify line mask:
            regions['base_fit'][i] = new_base
            regions['top_fit'][i] = new_top
            regions['wave_base_fit'][i] = spectrum['waveobs'][new_base]
            regions['wave_top_fit'][i] = spectrum['waveobs'][new_top]

            try:
                if new_base == 0 and new_top == 0:
                    # If the original 'region' structure was a simple one not created
                    # with __create_linemasks_structure, the values of base, peak and top
                    # are set to 0 and we have to find the spectra window by using
                    # wave_base and wave_top
                    wave_filter = (spectrum['waveobs'] >= regions['wave_base'][i]) & (spectrum['waveobs'] <= regions['wave_top'][i])
                    spectrum_window = spectrum[wave_filter]
                else:
                    spectrum_window = spectrum[new_base:new_top+1]

                line_model, rms = __fit_line(spectrum_window, continuum_model, regions['wave_peak'][i], discard_gaussian = discard_gaussian, discard_voigt = discard_voigt, baseline_margin=continuum_adjustment_margin, free_mu=free_mu)

                if free_mu and (line_model.mu() <= spectrum_window['waveobs'][0] or line_model.mu() >= spectrum_window['waveobs'][-1]):
                    raise Exception("Fitted wave peak (mu) outside the limits!")

                # Calculate wave_peak position error derived from the velocity error calculation based on:
                # Zucker 2003, "Cross-correlation and maximum-likelihood analysis: a new approach to combining cross-correlation functions"
                # http://adsabs.harvard.edu/abs/2003MNRAS.342.1291Z
                nbins = len(spectrum_window)
                inverted_fluxes = 1-spectrum_window['flux']
                distance = spectrum_window['waveobs'][1] - spectrum_window['waveobs'][0]
                first_derivative = np.gradient(inverted_fluxes, distance)
                second_derivative = np.gradient(first_derivative, distance)
                ## Using the exact velocity, the resulting error are less coherents (i.e. sometimes you can get lower errors when using bigger steps):
                #second_derivative_peak = np.interp(line_model.mu(), spectrum_window['waveobs'], second_derivative)
                #inverted_fluxes_peak = line_model.mu()
                ## More coherent results:
                peak = spectrum_window['waveobs'].searchsorted(line_model.mu())
                inverted_fluxes_peak = inverted_fluxes[peak]
                second_derivative_peak = second_derivative[peak]
                if inverted_fluxes_peak == 0:
                    inverted_fluxes_peak = 1e-10
                if second_derivative_peak == 0:
                    second_derivative_peak = 1e-10
                sharpness = second_derivative_peak / inverted_fluxes_peak

                denominator = (1 - np.power(inverted_fluxes_peak, 2))
                if denominator != 0:
                    line_snr = np.power(inverted_fluxes_peak, 2) / denominator
                else:
                    line_snr = 0.

                denominator = (nbins * sharpness * line_snr)
                if denominator != 0:
                    # Use abs instead of a simple '-1*' because sometime the result is negative and the sqrt cannot be calculated
                    error = np.sqrt(np.abs(1 / denominator))
                else:
                    error = 0
                #print line_model.mu(), error, "=", nbins, sharpness, line_snr
                line_model.set_emu(error)

                regions['mu'][i] = line_model.mu()
                regions['mu_err'][i] = line_model.emu()
                regions['sig'][i] = line_model.sig()
                regions['A'][i] = line_model.A()
                regions['baseline'][i] = line_model.baseline()
                if type(line_model) == VoigtModel:
                    regions['gamma'][i] = line_model.gamma()
                else:
                    # The model is Gaussian, do not use 'gamma'
                    regions['gamma'][i] = 9999.0

                regions['fwhm'][i], regions['fwhm_kms'][i] = line_model.fwhm()
                regions['R'][i] = regions['mu'][i] / regions['fwhm'][i] # Resolving power

                # Depth of the peak with respect to the total continuum in % over the total continuum
                # - In case that the peak is higher than the continuum, depth < 0
                continuum = line_model.baseline()
                flux = line_model(line_model.mu())
                regions['depth_fit'][i] = 1. - (flux / continuum)
                # Relative depth is "peak - mean_base_point" with respect to the total continuum
                # - In case that the mean base point is higher than the continuum, relative_depth < 0
                # - relative_depth < depth is always true
                flux_from_top_base_point_to_continuum = np.abs(continuum - np.max(spectrum_window['flux']))
                regions['relative_depth_fit'][i] = ((continuum - (flux + flux_from_top_base_point_to_continuum)) / continuum)

                # Equivalent Width
                # - Include 99.9999998% of the gaussian area
                from_x = regions['mu'][i] - 6*regions['sig'][i]
                to_x = regions['mu'][i] + 6*regions['sig'][i]
                if type(line_model) is GaussianModel:
                    # If it is a gaussian we can directly use a formule (but not if it is a voigt!)
                    regions['integrated_flux'][i] = -1.*regions['A'][i]*np.sqrt(2*np.pi*regions['sig'][i]**2) # nm
                    regions['ew'][i] = regions['integrated_flux'][i] / line_model.baseline() # nm
                else:
                    regions['integrated_flux'][i] = -1 * line_model.integrate(from_x, to_x) # nm^2
                    regions['ew'][i] = regions['integrated_flux'][i] / line_model.baseline() # nm
                regions['ewr'][i] = np.log10(regions['ew'][i] / regions['mu'][i])
                regions['ew'][i] *= 10000. # from nm to mA
                # EW error from
                # Vollmann & Eversberg, 2006: http://adsabs.harvard.edu/abs/2006AN....327..862V
                if snr is not None:
                    # Wavelength difference used in the fit
                    diff_wavelength_fit = spectrum['waveobs'][regions['top_fit'][i]] - spectrum['waveobs'][regions['base_fit'][i]]
                    # Wavelength difference covered by the 99.9999998% of the gaussian area
                    diff_wavelength_gaussian = to_x - from_x
                    if diff_wavelength_gaussian > diff_wavelength_fit:
                        # The line mask was too narrow, it is more fair to consider the whole gaussian range
                        diff_wavelength = 6*regions['sig'][i]
                        wfilter = np.logical_and(spectrum['waveobs'] >= from_x, spectrum['waveobs'] <= to_x)
                        local_err = spectrum['err'][wfilter]
                        local_flux = spectrum['flux'][wfilter]
                        local_waveobs = spectrum['waveobs'][wfilter]
                    else:
                        diff_wavelength = diff_wavelength_fit
                        local_err = spectrum['err'][regions['base_fit'][i]:regions['top_fit'][i]]
                        local_flux = spectrum['flux'][regions['base_fit'][i]:regions['top_fit'][i]]
                        local_waveobs = spectrum['waveobs'][regions['base_fit'][i]:regions['top_fit'][i]+1]
                    zeros = local_err == 0
                    if not np.all(zeros):
                        local_snr = np.median(local_flux[~zeros] / local_err[~zeros])
                    else:
                        local_snr = snr
                    regions['snr'][i] = local_snr
                    regions['diff_wavelength'][i] = diff_wavelength

                    mean_flux = np.mean(local_flux)
                    mean_flux_continuum = np.mean(continuum_model(local_waveobs))
                    regions['mean_flux'][i] = mean_flux
                    regions['mean_flux_continuum'][i] = mean_flux_continuum
                    if mean_flux != 0  and local_snr != 0:
                        regions['ew_err'][i] = np.sqrt(1 + mean_flux_continuum / mean_flux) * ((diff_wavelength*10000 - (regions['ew'][i]))/local_snr)
                    else:
                        regions['ew_err'][i] = 0.
                    #print "%.2f\t%.2f\t%.2f" % (regions['ew'][i], regions['ew_err'][i], regions['ew_err'][i] / regions['ew'][i])
                # RMS
                regions['rms'][i] = rms
            except Exception as e:
                #print "WARNING: Bad line fit (", i, ") - ", e.message
                fitting_not_possible = True


        current_work_progress = ((i*1.0)/total_regions) * 100
        if report_progress(current_work_progress, last_reported_progress):
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)
            if frame is not None:
                frame.update_progress(current_work_progress)

    if atomic_linelist is not None:
        logging.info("Cross matching with atomic data...")
        regions = __fill_linemasks_with_atomic_data(regions, atomic_linelist, diff_limit=max_atomic_wave_diff, vel_atomic=0.0, crossmatch_with_mu=crossmatch_with_mu, closest_match=closest_match)
    if telluric_linelist is not None and vel_telluric is not None:
        logging.info("Cross matching with telluric data...")
        regions = __fill_linemasks_with_telluric_info(regions, telluric_linelist, vel_telluric=vel_telluric)

    return regions


def __fill_linemasks_with_telluric_info(linemasks, telluric_linelist, vel_telluric=0.0):
    """
    Adds info to those linemasks that can be affected by tellurics.
    Different nearby linemasks can be affected by the same telluric line.
    """
    # Sort before treating
    linemasks.sort(order=['wave_peak'])

    if vel_telluric != 0:
        # Speed of light in m/s
        c = 299792458.0
        # Radial/barycentric velocity from km/s to m/s
        vel_telluric = vel_telluric * 1000

        # Correct wavelength scale for radial velocity
        original_wave_peak = linemasks['wave_peak'].copy()
        original_mu = linemasks['mu'].copy()
        linemasks['wave_peak'] = linemasks['wave_peak'] / ((vel_telluric / c) + 1)
        linemasks['mu'] = linemasks['mu'] / ((vel_telluric / c) + 1)

    # The discarded flag in the telluric linelist is not a good one because there are clear lines mark as true (i.e. 628.0392 nm)
    #telluric_linelist = telluric_linelist[telluric_linelist['discarded'] != True]
    # It is better to clear the not fitted telluric lines:
    telluric_linelist = telluric_linelist[telluric_linelist['rms'] < 9999]
    telluric_linelist.sort(order=['wave_peak'])

    clean_linemasks = linemasks[linemasks['wave_peak'] != 0]
    max_wave_peak = np.max(clean_linemasks['wave_peak'])
    min_wave_peak = np.min(clean_linemasks['wave_peak'])

    if telluric_linelist['wave_peak'][0] > min_wave_peak or telluric_linelist['wave_peak'][-1] < max_wave_peak:
        logging.warn("Telluric linelist does not cover the whole linemask wavelength range")
        logging.warn("- Telluric range from " + str(telluric_linelist['wave_peak'][0]) + " to " + str(telluric_linelist['wave_peak'][-1]) + " nm")
        logging.warn("- Linemask range from " + str(min_wave_peak) + " to " + str(max_wave_peak) + " nm")

    diff_limit = np.max(telluric_linelist["fwhm"])
    wfilter = (telluric_linelist['wave_peak'] >= min_wave_peak - diff_limit) & (telluric_linelist['wave_peak'] <= max_wave_peak + diff_limit)
    telluric_linelist = telluric_linelist[wfilter]

    if len(telluric_linelist) > 0:
        for j in np.arange(len(linemasks)):
            if linemasks['mu'][j] == 0:
                # This lines has not been fitted correctly, it will be discarded
                continue
            # Find index of the nearest line
            diff = telluric_linelist['wave_peak'] - linemasks['mu'][j]
            diff_limit = np.max((linemasks["fwhm"][j], 0.005)) # At least 0.005
            abs_diff = np.abs(diff)
            i = np.argmin(abs_diff)
            ### Save the information
            if abs_diff[i] <= diff_limit:
                linemasks["telluric_wave_peak"][j] = telluric_linelist['wave_peak'][i]
                linemasks["telluric_depth"][j] = telluric_linelist["depth"][i]
                linemasks["telluric_fwhm"][j] = telluric_linelist["fwhm"][i]
                linemasks['telluric_R'][j] = linemasks['mu'][j] / (linemasks['fwhm'][j] - linemasks['telluric_fwhm'][j]) # Resolving power

    if vel_telluric != 0:
        linemasks['wave_peak'] = original_wave_peak
        linemasks['mu'] = original_mu

    return linemasks


def __fill_linemasks_with_atomic_data(linemasks, atomic_linelist, diff_limit=0.0005, vel_atomic=0.0, crossmatch_with_mu=False, closest_match=False):
    """
    Cross-match linemasks with a atomic linelist in order to find
    the nearest lines and copy the information into the linemasks structure.
    """
    # Sort before treating
    linemasks.sort(order=['wave_peak'])

    if vel_atomic != 0:
        # Speed of light in m/s
        c = 299792458.0
        # Radial/barycentric velocity from km/s to m/s
        vel_atomic = vel_atomic * 1000

        # Correct wavelength scale for radial velocity
        original_wave_peak = linemasks['wave_peak'].copy()
        original_mu = linemasks['mu'].copy()
        linemasks['wave_peak'] = linemasks['wave_peak'] / ((vel_atomic / c) + 1)
        linemasks['mu'] = linemasks['mu'] / ((vel_atomic / c) + 1)


    ## Sort
    #atomic_linelist.sort(order=['wave_nm']) # Ascending

    clean_linemasks = linemasks[linemasks['wave_peak'] != 0]
    max_wave_peak = np.max(clean_linemasks['wave_peak'])
    min_wave_peak = np.min(clean_linemasks['wave_peak'])

    if atomic_linelist['wave_nm'][0] > min_wave_peak or atomic_linelist['wave_nm'][-1] < max_wave_peak:
        print "WARNING: The atomic linelist does not cover the whole linemask wavelength range"
        print "- Atomic line's range from", atomic_linelist['wave_nm'][0], "to", atomic_linelist['wave_nm'][-1], "nm"
        print "- Linemask range from", min_wave_peak, "to", max_wave_peak, "nm"

    wfilter = (atomic_linelist['wave_nm'] >= min_wave_peak - diff_limit) & (atomic_linelist['wave_nm'] <= max_wave_peak + diff_limit)
    atomic_linelist = atomic_linelist[wfilter]

    if len(atomic_linelist) > 0:
        for j in np.arange(len(linemasks)):
            if crossmatch_with_mu:
                diff = atomic_linelist['wave_nm'] - linemasks['mu'][j]
            else:
                diff = atomic_linelist['wave_nm'] - linemasks['wave_peak'][j]

            # Find index of the nearest lines
            abs_diff = np.abs(diff)
            imin_diff = np.where(abs_diff <= diff_limit)[0]
            if len(imin_diff) == 0:
                continue
            else:
                if closest_match or \
                    (np.all(atomic_linelist['theoretical_depth'][imin_diff] == 0) and \
                        np.all(atomic_linelist['theoretical_ew'][imin_diff] == 0)):
                    # Just select the closest line
                    i = np.argmin(abs_diff)
                else:
                    # Select the line that has the biggest theoretical depth (usually calculated with a solar spectrum)
                    max_depth = np.max(atomic_linelist['theoretical_depth'][imin_diff])
                    imax_depth = np.where(atomic_linelist['theoretical_depth'][imin_diff] >= max_depth)[0]
                    if len(imax_depth) > 1:
                        # For the same theoretical depth, select the line that has the biggest EW (usually calculated with a solar spectrum)
                        iii = np.argmax(atomic_linelist['theoretical_ew'][imin_diff[imax_depth]])
                        ii = imax_depth[iii]
                        i = imin_diff[ii]
                    else:
                        ii = np.argmax(atomic_linelist['theoretical_depth'][imin_diff])
                        i = imin_diff[ii]

                # Copy atomic information
                for key, dtype in _get_atomic_linelist_definition():
                    linemasks[key][j] = atomic_linelist[key][i]

    if vel_atomic != 0:
        linemasks['wave_peak'] = original_wave_peak
        linemasks['mu'] = original_mu

    return linemasks

def __find_peaks_and_base_points(xcoord, yvalues):
    """
    Find peaks and base points. It works better with a smoothed spectrum (i.e. convolved using 2*resolution).
    """
    if len(yvalues[~np.isnan(yvalues)]) == 0 or len(yvalues[~np.isnan(xcoord)]) == 0:
        #raise Exception("Not enough data for finding peaks and base points")
        print "WARNING: Not enough data for finding peaks and base points"
        peaks = []
        base_points = []
    else:
        # Determine peaks and base points (also known as continuum points)
        peaks = find_local_min_values(yvalues)
        base_points = find_local_max_values(yvalues)

        # WARNING: Due to three or more consecutive values with exactly the same flux
        # find_local_max_values or find_local_min_values will identify all of them as peaks or bases,
        # where only one of the should be marked as peak or base.
        # These cases break the necessary condition of having the same number of
        # peaks and base_points +/-1
        # It is necessary to find those "duplicates" and remove them:
        peaks = __remove_consecutives_features(peaks)
        base_points = __remove_consecutives_features(base_points)

        if not (len(peaks) - len(base_points)) in [-1, 0, 1]:
            raise Exception("This should not happen")

        # Make sure that
        peaks, base_points = __assert_structure(xcoord, yvalues, peaks, base_points)

    return peaks, base_points


def __improve_linemask_edges(xcoord, yvalues, base, top, peak):
    """
    Given a spectrum, the position of a peak and its limiting region where:

    - Typical shape: concave + convex + concave region.
    - Peak is located within the convex region, although it is just in the border with the concave region (first element).
    """
    # Try to use two additional position, since we are going to lose them
    # by doing the first and second derivative
    original_top = top
    top = np.min([top+2, len(xcoord)])
    y = yvalues[base:top+1]
    x = xcoord[base:top+1]
    # First derivative (positive => flux increases, negative => flux decreases)
    dy_dx = (y[:-1] - y[1:]) / (x[:-1] - x[1:])
    # Second derivative (positive => convex, negative => concave)
    d2y_dx2 = (dy_dx[:-1] - dy_dx[1:])  / (x[:-2] - x[2:])
    # Peak position inside the linemask region
    peak_relative_pos = peak - base
    # The peak should be in a convex region => the second derivative should be positive
    # - It may happen that the peak falls in the beginning/end of a concave region, accept also this cases
    if peak_relative_pos < len(d2y_dx2)-1 and peak_relative_pos > 0 and (d2y_dx2[peak_relative_pos-1] > 0 or d2y_dx2[peak_relative_pos] > 0 or d2y_dx2[peak_relative_pos+1] > 0):
        # Find the concave positions at both sides of the peak
        concave_pos = np.where(d2y_dx2<0)[0]
        if len(concave_pos) == 0:
            # This should not happen, but just in case...
            new_base = base
            new_top = original_top
        else:
            # Concave regions to the left of the peak
            left_concave_pos = concave_pos[concave_pos-peak_relative_pos < 0]
            if len(left_concave_pos) == 0:
                # This should not happen, but just in case...
                new_base = base
            else:
                # Find the edges of the concave regions to the left of the peak
                left_concave_pos_diff = left_concave_pos[:-1] - left_concave_pos[1:]
                left_concave_pos_diff = np.asarray([-1] + left_concave_pos_diff.tolist())
                left_concave_edge_pos = np.where(left_concave_pos_diff != -1)[0]
                if len(left_concave_edge_pos) == 0:
                    # There is only one concave region, we use its left limit
                    new_base = left_concave_pos[0] + base
                else:
                    # There is more than one concave region, select the nearest edge to the peak
                    new_base = np.max([left_concave_pos[np.max(left_concave_edge_pos)] + base, base])

            # Concave regions to the right of the peak
            right_concave_pos = concave_pos[concave_pos-peak_relative_pos > 0]
            if len(right_concave_pos) == 0:
                # This should not happen, but just in case...
                new_top = original_top
            else:
                # Find the edges of the concave regions to the right of the peak
                right_concave_pos_diff = right_concave_pos[1:] - right_concave_pos[:-1]
                right_concave_pos_diff = np.asarray(right_concave_pos_diff.tolist() + [1])
                right_concave_edge_pos = np.where(right_concave_pos_diff != 1)[0]
                if len(right_concave_edge_pos) == 0:
                    # There is only one concave region, we use its right limit
                    new_top = right_concave_pos[-1] + base
                else:
                    # There is more than one concave region, select the one with the nearest edge to the peak
                    new_top = np.min([right_concave_pos[np.min(right_concave_edge_pos)] + base, original_top])

    else:
        # This will happen very rarely (only in peaks detected at the extreme of a spectrum
        # and one of its basepoints has been "artificially" added and it happens to be
        # just next to the peak)
        new_base = base
        new_top = original_top
        #plt.plot(x, y)
        #l = plt.axvline(x = x[peak_relative_pos], linewidth=1, color='red')
        #print d2y_dx2, d2y_dx2[peak_relative_pos]
        #plt.show()

    return new_base, new_top

def adjust_linemasks(spectrum, linemasks, max_margin=0.5, min_margin=0.0):
    """
    Adjust the line masks borders to the shape of the line in the specified spectrum.
    It will consider by default 0.5 nm around the wave peak.

    It tries to respect the min_margin around the peak but in case of overlapping
    regions, it might not be respected.

    It returns a new linemasks structure.
    """
    linemasks = linemasks.copy()
    for line in linemasks:
        wave_peak = line['wave_peak']
        wfilter = np.logical_and(spectrum['waveobs'] >= wave_peak - max_margin, spectrum['waveobs'] <= wave_peak + max_margin)
        spectrum_window = spectrum[wfilter]
        if len(spectrum_window) < 3:
            continue
        peaks, base_points = __find_peaks_and_base_points(spectrum_window['waveobs'], spectrum_window['flux'])
        ipeak = spectrum_window['waveobs'].searchsorted(wave_peak)
        if len(base_points) > 0:
            # wave_peak could not be exactly in a peak, so consider the nearest peak
            # to have only in consideration the maximum values that are after/before it
            nearest_peak = peaks[np.argmin(np.abs(peaks - ipeak))]
            left_base_points = base_points[base_points < nearest_peak]
            if len(left_base_points) > 0:
                wave_base = spectrum_window['waveobs'][left_base_points[-1]] # nearest to the peak
                # Make sure that the wave_base is smaller than the wave_peak
                i = -2
                while wave_base > wave_peak and -1*i <= len(left_base_points):
                    wave_base = spectrum_window['waveobs'][left_base_points[i]] # nearest to the peak
                    i -= 1
                if wave_base < wave_peak:
                    line['wave_base'] =  wave_base

            right_base_points = base_points[base_points > nearest_peak]
            if len(right_base_points) > 0:
                wave_top = spectrum_window['waveobs'][right_base_points[0]] # nearest to the peak
                # Make sure that the wave_top is bigger than the wave_peak
                i = 1
                while wave_top < wave_peak and i < len(right_base_points):
                    wave_top = spectrum_window['waveobs'][right_base_points[i]] # nearest to the peak
                    i += 1
                if wave_top > wave_peak:
                    line['wave_top'] =  wave_top

            if min_margin > 0:
                line['wave_base'] = np.min((line['wave_peak']-min_margin, line['wave_base']))
                line['wave_top'] = np.max((line['wave_peak']+min_margin, line['wave_top']))

    # Correct potential overlapping between line masks
    # except if their base and top wave are exactly the same, those cases could
    # be isotopic lines that belong to the same absorption line
    linemasks.sort(order=['wave_peak'])
    for i, line in enumerate(linemasks[:-1]):
        if line['wave_top'] > linemasks['wave_base'][i+1] and \
                not (np.abs(line['wave_base'] - linemasks['wave_base'][i+1]) < 1e-5 and\
                    np.abs(line['wave_top'] - linemasks['wave_top'][i+1]) < 1e-5):
            mean_from_edges = (line['wave_top'] + linemasks['wave_base'][i+1]) / 2.0
            # Make sure that the new limit (wave_top/wave_base) is bigger/smaller
            # than the wave_peaks of the two lines
            if mean_from_edges > line['wave_peak'] and mean_from_edges < linemasks['wave_peak'][i+1]:
                line['wave_top'] = mean_from_edges
                linemasks['wave_base'][i+1] = mean_from_edges
            elif line['wave_top'] < linemasks['wave_peak'][i+1]:
                linemasks['wave_base'][i+1] = line['wave_top']
            elif linemasks['wave_base'][i+1] > line['wave_peak']:
                line['wave_top'] = linemasks['wave_base'][i+1]
            else:
                mean_from_peaks = (line['wave_peak'] + linemasks['wave_peak'][i+1]) / 2.0
                line['wave_top'] = mean_from_peaks
                linemasks['wave_base'][i+1] = mean_from_peaks

    return linemasks


############## [start] Radial velocity


def __cross_correlation_function_template(spectrum, template, lower_velocity_limit, upper_velocity_limit, velocity_step, frame=None):
    """
    Calculates the cross correlation value between the spectrum and the specified template
    by shifting the template from lower to upper velocity.

    - The spectrum and the template should be uniformly spaced in terms of velocity (which
      implies non-uniformly distributed in terms of wavelength).
    - The velocity step used for the construction of the template should be the same
      as the one specified in this function.
    - The lower/upper/step velocity is only used to determine how many shifts
      should be done (in array positions) and return a velocity grid.
    """

    last_reported_progress = -1
    if frame is not None:
        frame.update_progress(0)

    # Speed of light in m/s
    c = 299792458.0

    velocity = np.arange(lower_velocity_limit, upper_velocity_limit+velocity_step, velocity_step)
    # 1 shift = 0.5 km/s (or the specified value)
    shifts = np.int32(velocity / velocity_step)

    num_shifts = len(shifts)
    # Cross-correlation function
    ccf = np.zeros(num_shifts)
    ccf_err = np.zeros(num_shifts)
    depth = np.abs(np.max(template['flux']) - template['flux'])
    for i, vel in enumerate(velocity):
        factor = np.sqrt((1.-(vel*1000.)/c)/(1.+(vel*1000.)/c))
        shifted_template = np.interp(spectrum['waveobs'], template['waveobs']/factor, depth, left=0.0, right=0.0)
        ccf[i] = np.correlate(spectrum['flux'], shifted_template)[0]
        ccf_err[i] = np.correlate(spectrum['err'], shifted_template)[0] # Propagate errors
        current_work_progress = ((i*1.0)/num_shifts) * 100
        if report_progress(current_work_progress, last_reported_progress):
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)
            if frame is not None:
                frame.update_progress(current_work_progress)

    max_ccf = np.max(ccf)
    ccf = ccf/max_ccf # Normalize
    ccf_err = ccf_err/max_ccf # Propagate errors

    return velocity, ccf, ccf_err


def _sampling_uniform_in_velocity(wave_base, wave_top, velocity_step):
    """
    Create a uniformly spaced grid in terms of velocity:

    - An increment in position (i => i+1) supposes a constant velocity increment (velocity_step).
    - An increment in position (i => i+1) does not implies a constant wavelength increment.
    - It is uniform in log(wave) since:
          Wobs = Wrest * (1 + Vr/c)^[1,2,3..]
          log10(Wobs) = log10(Wrest) + [1,2,3..] * log10(1 + Vr/c)
      The last term is constant when dealing with wavelenght in log10.
    - Useful for building the cross correlate function used for determining the radial velocity of a star.
    """
    # Speed of light
    c = 299792.4580 # km/s
    #c = 299792458.0 # m/s

    ### Numpy optimized:
    # number of elements to go from wave_base to wave_top in increments of velocity_step
    i = int(np.ceil( (c * (wave_top - wave_base)) / (wave_base*velocity_step)))
    grid = wave_base * np.power((1 + (velocity_step / c)), np.arange(i)+1)

    # Ensure wavelength limits since the "number of elements i" tends to be overestimated
    wfilter = grid <= wave_top
    grid = grid[wfilter]

    ### Non optimized:
    #grid = []
    #next_wave = wave_base
    #while next_wave <= wave_top:
        #grid.append(next_wave)
        ### Newtonian version:
        #next_wave = next_wave + next_wave * ((velocity_step) / c) # nm
        ### Relativistic version:
        ##next_wave = next_wave + next_wave * (1.-np.sqrt((1.-(velocity_step*1000.)/c)/(1.+(velocity_step*1000.)/c)))

    return np.asarray(grid)

def __cross_correlation_function_uniform_in_velocity(spectrum, mask, lower_velocity_limit, upper_velocity_limit, velocity_step, mask_size=2.0, mask_depth=0.01, template=False, fourier=False, frame=None):
    """
    Calculates the cross correlation value between the spectrum and the specified mask
    by shifting the mask from lower to upper velocity.

    - The spectrum and the mask should be uniformly spaced in terms of velocity (which
      implies non-uniformly distributed in terms of wavelength).
    - The velocity step used for the construction of the mask should be the same
      as the one specified in this function.
    - The lower/upper/step velocity is only used to determine how many shifts
      should be done (in array positions) and return a velocity grid.

    If fourier is set, the calculation is done in the fourier space. More info:

        VELOCITIES FROM CROSS-CORRELATION: A GUIDE FOR SELF-IMPROVEMENT
        CARLOS ALLENDE PRIETO
        http://iopscience.iop.org/1538-3881/134/5/1843/fulltext/205881.text.html
        http://iopscience.iop.org/1538-3881/134/5/1843/fulltext/sourcecode.tar.gz
    """

    last_reported_progress = -1
    if frame is not None:
        frame.update_progress(0)

    # Speed of light in m/s
    c = 299792458.0

    # 1 shift = 1.0 km/s (or the specified value)
    shifts = np.arange(np.int32(np.floor(lower_velocity_limit)/velocity_step), np.int32(np.ceil(upper_velocity_limit)/velocity_step)+1)
    velocity = shifts * velocity_step

    waveobs = _sampling_uniform_in_velocity(np.min(spectrum['waveobs']), np.max(spectrum['waveobs']), velocity_step)
    flux = np.interp(waveobs, spectrum['waveobs'], spectrum['flux'], left=0.0, right=0.0)
    err = np.interp(waveobs, spectrum['waveobs'], spectrum['err'], left=0.0, right=0.0)


    if template:
        depth = np.abs(np.max(mask['flux']) - mask['flux'])
        resampled_mask = np.interp(waveobs, mask['waveobs'], depth, left=0.0, right=0.0)
    else:
        selected = __select_lines_for_mask(mask, minimum_depth=mask_depth, velocity_mask_size = mask_size, min_velocity_separation = 1.0)
        resampled_mask = __create_mask(waveobs, mask['wave_peak'][selected], mask['depth'][selected], velocity_mask_size=mask_size)

    if fourier:
        # Transformed flux and mask
        tflux = fft(flux)
        tresampled_mask = fft(resampled_mask)
        conj_tresampled_mask = np.conj(tresampled_mask)
        num = len(resampled_mask)/2+1
        tmp = abs(ifft(tflux*conj_tresampled_mask))
        ccf = np.hstack((tmp[num:], tmp[:num]))

        # Transformed flux and mask powered by 2 (second)
        #ccf_err = np.zeros(len(ccf))
        # Conservative error propagation
        terr = fft(err)
        tmp = abs(ifft(terr*conj_tresampled_mask))
        ccf_err = np.hstack((tmp[num:], tmp[:num]))
        ## Error propagation
        #tflux_s = fft(np.power(flux, 2))
        #tresampled_mask_s = fft(np.power(resampled_mask, 2))
        #tflux_err_s = fft(np.power(err, 2))
        #tresampled_mask_err_s = fft(np.ones(len(err))*0.05) # Errors of 5% for masks

        #tmp = abs(ifft(tflux_s*np.conj(tresampled_mask_err_s)))
        #tmp += abs(ifft(tflux_err_s*np.conj(tresampled_mask_s)))
        #ccf_err = np.hstack((tmp[num:], tmp[:num]))
        #ccf_err = np.sqrt(ccf_err)

        # Velocities
        velocities = velocity_step * (np.arange(len(resampled_mask), dtype=float)+1 - num)

        # Filter to area of interest
        xfilter = np.logical_and(velocities >= lower_velocity_limit, velocities <= upper_velocity_limit)
        ccf = ccf[xfilter]
        ccf_err = ccf_err[xfilter]
        velocities = velocities[xfilter]
    else:
        num_shifts = len(shifts)
        # Cross-correlation function
        ccf = np.zeros(num_shifts)
        ccf_err = np.zeros(num_shifts)

        for shift, i in zip(shifts, np.arange(num_shifts)):
            #shifted_mask = resampled_mask
            if shift == 0:
                shifted_mask = resampled_mask
            elif shift > 0:
                #shifted_mask = np.hstack((shift*[0], resampled_mask[:-1*shift]))
                shifted_mask = np.hstack((resampled_mask[-1*shift:], resampled_mask[:-1*shift]))
            else:
                #shifted_mask = np.hstack((resampled_mask[-1*shift:], -1*shift*[0]))
                shifted_mask = np.hstack((resampled_mask[-1*shift:], resampled_mask[:-1*shift]))
            #ccf[i] = np.correlate(flux, shifted_mask)[0]
            #ccf_err[i] = np.correlate(err, shifted_mask)[0] # Propagate errors
            ccf[i] = np.average(flux*shifted_mask)
            ccf_err[i] = np.average(err*shifted_mask) # Propagate errors
            #ccf[i] = np.average(np.tanh(flux*shifted_mask))
            #ccf_err[i] = np.average(np.tanh(err*shifted_mask)) # Propagate errors

            current_work_progress = ((i*1.0)/num_shifts) * 100
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                logging.info("%.2f%%" % current_work_progress)
                if frame is not None:
                    frame.update_progress(current_work_progress)

    max_ccf = np.max(ccf)
    ccf = ccf/max_ccf # Normalize
    ccf_err = ccf_err/max_ccf # Propagate errors

    return velocity, ccf, ccf_err, len(flux)



def create_filter_for_regions_affected_by_tellurics(wavelengths, linelist_telluric, min_velocity=-30.0, max_velocity=30.0, frame=None):
    """
    Returns a boolean array of the same size of wavelengths. True will be assigned
    to those positions thay may be affected by telluric lines in a range from
    min_velocity to max_velocity
    """
    # Light speed in vacuum
    c = 299792458.0 # m/s

    tfilter = wavelengths == np.nan
    tfilter2 = wavelengths == np.nan
    wave_bases = linelist_telluric['wave_peak'] * np.sqrt((1.-(max_velocity*1000.)/c)/(1.+(max_velocity*1000.)/c))
    wave_tops = linelist_telluric['wave_peak'] * np.sqrt((1.-(min_velocity*1000.)/c)/(1.+(min_velocity*1000.)/c))
    last_reported_progress = -1
    total_regions = len(wave_bases)
    last = 0 # Optimization
    for i, (wave_base, wave_top) in enumerate(zip(wave_bases, wave_tops)):
        begin = wavelengths[last:].searchsorted(wave_base)
        end = wavelengths[last:].searchsorted(wave_top)
        tfilter[last+begin:last+end] = True
        #wfilter = np.logical_and(wavelengths >= wave_base, wavelengths <= wave_top)
        #tfilter = np.logical_or(wfilter, tfilter)

        current_work_progress = ((i*1.0)/total_regions) * 100
        if report_progress(current_work_progress, last_reported_progress):
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)
            if frame is not None:
                frame.update_progress(current_work_progress)
        last += end
    return tfilter



try:
    import pyximport
    import numpy as np
    pyximport.install(setup_args={'include_dirs':[np.get_include()]})
    from lines_c import create_mask as __create_mask
except:
    print "*********************************************************************"
    print "Not optimized version loaded!"
    print "*********************************************************************"

    def __create_mask(spectrum_wave, mask_wave, mask_values, velocity_mask_size=2.0):
        """
        It constructs a zero flux spectrum and assign mask values to the wavelengths
        belonging to that value and its surounds (determined by the velocity_mask_size).
        """
        ## Speed of light in m/s
        c = 299792458.0

        resampled_mask = np.zeros(len(spectrum_wave))

        # Mask limits
        mask_wave_step = (mask_wave * (1.-np.sqrt((1.-(velocity_mask_size*1000.)/c)/(1.+(velocity_mask_size*1000.)/c))))/2.0
        mask_wave_base = mask_wave - 1*mask_wave_step
        mask_wave_top = mask_wave + 1*mask_wave_step

        i = 0
        j = 0
        for i in xrange(len(mask_wave)):
            #j = 0
            while j < len(spectrum_wave) and spectrum_wave[j] < mask_wave_base[i]:
                j += 1
            while j < len(spectrum_wave) and spectrum_wave[j] >= mask_wave_base[i] and spectrum_wave[j] <= mask_wave_top[i]:
                resampled_mask[j] = mask_values[i]
                j += 1

        return resampled_mask

def __select_lines_for_mask(linemasks, minimum_depth=0.01, velocity_mask_size = 2.0, min_velocity_separation = 1.0):
    """
    Select the lines that are goint to be used for building a mask for doing
    cross-correlation. It filters by depth and validate that the lines are
    suficiently apart from its neighbors to avoid overlapping.

    For that purpose, 'velocity_mask_size' represents the masks size in km/s
    around the peak and optionally, 'min_velocity_separation' indicates the
    minimum separation needed between two consecutive masks.

    It returns a boolean array indicating what lines have been finally selected.
    """
    total_velocity_separation = velocity_mask_size + min_velocity_separation / 2.0
    selected = linemasks['depth'] >= minimum_depth

    ## Speed of light in m/s
    c = 299792458.0

    # Mask limits
    mask_wave_step = (linemasks['wave_peak'] * (1.-np.sqrt((1.-(total_velocity_separation*1000.)/c)/(1.+(total_velocity_separation*1000.)/c))))/2.0
    mask_wave_base = linemasks['wave_peak'] - 1*mask_wave_step
    mask_wave_top = linemasks['wave_peak'] + 1*mask_wave_step
    #mask_wave_base = linemasks['wave_base'] - 1*mask_wave_step
    #mask_wave_top = linemasks['wave_top'] + 1*mask_wave_step

    i = 0
    while i < len(linemasks):
        if selected[i]:
            # Right
            r = i
            max_r = r
            max_depth_r = linemasks['depth'][i]
            while r < len(linemasks) and mask_wave_base[r] <= mask_wave_top[i]:
                if selected[r] and linemasks['depth'][r] > max_depth_r:
                    max_depth_r = linemasks['depth'][r]
                    max_r = r
                r += 1
            # Left
            l = i
            max_l = l
            max_depth_l = linemasks['depth'][i]
            while l >= 0 and mask_wave_top[l] >= mask_wave_base[i]:
                if selected[l] and linemasks['depth'][l] > max_depth_l:
                    max_depth_l = linemasks['depth'][l]
                    max_l = l
                l -= 1

            if i - 1 == l and i + 1 == r:
                # No conflict
                i += 1
            else:
                #print "*",
                if i + 1 != l and i - 1 != r:
                    #print "both", i, i - l, r - i
                    for x in xrange(r - i):
                        selected[i+x] = False
                    for x in xrange(i - l):
                        selected[i-x] = False
                    if max_depth_l > max_depth_r:
                        selected[max_l] = True
                    else:
                        selected[max_r] = True
                elif i + 1 != l:
                    #print "left"
                    for x in xrange(i - l):
                        selected[i-x] = False
                    selected[max_l] = True
                else:
                    #print "right"
                    for x in xrange(r - i):
                        selected[i+x] = False
                    selected[max_r] = True
                i = r
        else:
            i += 1
    return selected



def cross_correlate_with_mask(spectrum, linelist, lower_velocity_limit=-200, upper_velocity_limit=200, velocity_step=1.0, mask_size=None, mask_depth=0.01, fourier=False, only_one_peak=False, model='2nd order polynomial + gaussian fit', peak_probability=0.75, frame=None):
    """
    Determines the velocity profile by cross-correlating the spectrum with
    a mask built from a line list mask.

    If mask_size is not specified, the double of the velocity step will be taken,
    which generally it is the recommended value.

    :returns:
        - Array with fitted gaussian models sorted by depth (deepest at position 0)
        - CCF structure with 'x' (velocities), 'y' (relative intensities), 'err'

    """
    if mask_size is None:
        mask_size = 2*velocity_step # Recommended
    return __cross_correlate(spectrum, linelist=linelist, template=None, \
                lower_velocity_limit=lower_velocity_limit, upper_velocity_limit = upper_velocity_limit, \
                velocity_step=velocity_step, \
                mask_size=mask_size, mask_depth=mask_depth, fourier=fourier, \
                only_one_peak=only_one_peak, peak_probability=peak_probability, model=model, \
                frame=None)

def cross_correlate_with_template(spectrum, template, lower_velocity_limit=-200, upper_velocity_limit=200, velocity_step=1.0, fourier=False, only_one_peak=False, model='2nd order polynomial + gaussian fit', peak_probability=0.75, frame=None):
    """
    Determines the velocity profile by cross-correlating the spectrum with
    a spectrum template.

    :returns:
        - Array with fitted gaussian models sorted by depth (deepest at position 0)
        - CCF structure with 'x' (velocities), 'y' (relative intensities), 'err'
    """
    return __cross_correlate(spectrum, linelist=None, template=template, \
            lower_velocity_limit=lower_velocity_limit, upper_velocity_limit = upper_velocity_limit, \
            velocity_step=velocity_step, \
            mask_size=None, mask_depth=None, fourier=fourier, \
            only_one_peak=only_one_peak, peak_probability=peak_probability, model=model, \
            frame=None)

def __cross_correlate(spectrum, linelist=None, template=None, lower_velocity_limit = -200, upper_velocity_limit = 200, velocity_step=1.0, mask_size=2.0, mask_depth=0.01, fourier=False, only_one_peak=False, peak_probability=0.75, model='2nd order polynomial + gaussian fit', frame=None):
    ccf, nbins = __build_velocity_profile(spectrum, \
            linelist = linelist, template = template, \
            lower_velocity_limit = lower_velocity_limit, upper_velocity_limit = upper_velocity_limit, \
            velocity_step=velocity_step, \
            mask_size=mask_size, mask_depth=mask_depth, \
            fourier=fourier, frame=frame)

    models = __model_velocity_profile(ccf, nbins, only_one_peak=only_one_peak, \
                                            peak_probability=peak_probability, model=model)
    # We have improved the peak probability detection using RLM, a priori it is not needed
    # this best selection:
    #best = select_good_velocity_profile_models(models, ccf)
    #return models[best], ccf
    return models, ccf

def __build_velocity_profile(spectrum, linelist=None, template=None, lower_velocity_limit = -200, upper_velocity_limit = 200, velocity_step=1.0, mask_size=2.0, mask_depth=0.01, fourier=False, frame=None):
    """
    Determines the velocity profile by cross-correlating the spectrum with:

    * a mask built from a line list if linelist is specified
    * a spectrum template if template is specified

    :returns:
        CCF structure with 'x' (velocities), 'y' (relative intensities), 'err'
        together with the number of spectrum's bins used in the cross correlation.

    """
    if linelist is not None:
        if template is not None:
            logging.warn("Building velocity profile with mask (ignoring template)")

        linelist = linelist[linelist['depth'] > 0.01]
        lfilter = np.logical_and(linelist['wave_peak'] >= np.min(spectrum['waveobs']), linelist['wave_peak'] <= np.max(spectrum['waveobs']))
        linelist = linelist[lfilter]

        velocity, ccf, ccf_err, nbins = __cross_correlation_function_uniform_in_velocity(spectrum, linelist, lower_velocity_limit, upper_velocity_limit, velocity_step, mask_size=mask_size, mask_depth=mask_depth, fourier=fourier, frame=frame)
    elif template is not None:
        ## Obtain the cross-correlate function by shifting the template
        velocity, ccf, ccf_err, nbins = __cross_correlation_function_uniform_in_velocity(spectrum, template, lower_velocity_limit, upper_velocity_limit, velocity_step, template=True, fourier=False, frame=frame)
        #velocity, ccf, ccf_err = __cross_correlation_function_template(spectrum, template, lower_velocity_limit = lower_velocity_limit, upper_velocity_limit=upper_velocity_limit, velocity_step = velocity_step, frame=frame)

    else:
        raise Exception("A linelist or template should be specified")

    ccf_struct = np.recarray((len(velocity), ), dtype=[('x', float),('y', float), ('err', float)])
    ccf_struct['x'] = velocity
    ccf_struct['y'] = ccf
    ccf_struct['err'] = ccf_err
    return ccf_struct, nbins


def __model_velocity_profile(ccf, nbins, only_one_peak=False, peak_probability=0.55, model='2nd order polynomial + gaussian fit'):
    """
    Fits a model ('Gaussian' or 'Voigt') to the deepest peaks in the velocity
    profile. If it is 'Auto', a gaussian and a voigt will be fitted and the best
    one used.

    In all cases, the peak is located by fitting a 2nd degree polynomial. Afterwards,
    the gaussian/voigt fitting is done for obtaining more info (such as sigma, etc.)

    * For Radial Velocity profiles, more than 1 outlier peak implies that the star is a spectroscopic binary.

    WARNING: fluxes and errors are going to be modified by a linear normalization process

    Detected peaks are evaluated to discard noise, a probability is assigned to each one
    in function to a linear model. If more than one peak is found, those with a peak probability
    lower than the specified by the argument will be discarded.

    :returns:
        Array of fitted models and an array with the margin errors for model.mu() to be able to know the interval
        of 99% confiance.

    """
    models = []
    if len(ccf) == 0:
        return models
    xcoord = ccf['x']
    fluxes = ccf['y']
    errors = ccf['err']

    # Smooth flux
    sig = 1
    smoothed_fluxes = scipy.ndimage.filters.gaussian_filter1d(fluxes, sig)
    #smoothed_fluxes = fluxes
    # Finding peaks and base points
    peaks, base_points = __find_peaks_and_base_points(xcoord, smoothed_fluxes)

    if len(peaks) == 0 or len(base_points) == 0:
        return models

    if len(peaks) != 0:
        base = base_points[:-1]
        top = base_points[1:]
        # Adjusting edges
        new_base = np.zeros(len(base), dtype=int)
        new_top = np.zeros(len(base), dtype=int)
        for i in np.arange(len(peaks)):
            new_base[i], new_top[i] = __improve_linemask_edges(xcoord, smoothed_fluxes, base[i], top[i], peaks[i])
            #new_base[i] = base[i]
            #new_top[i] = top[i]
        base = new_base
        top = new_top

        if only_one_peak:
            # Just try with the deepest line
            selected_peaks_indices = []
        else:
            import statsmodels.api as sm
            #x = np.arange(len(peaks))
            #y = fluxes[peaks]
            x = xcoord
            y = fluxes
            # RLM (Robust least squares)
            # Huber's T norm with the (default) median absolute deviation scaling
            # - http://en.wikipedia.org/wiki/Huber_loss_function
            # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
            linear_model = huber_t.fit()
            selected_peaks_indices = np.where(linear_model.weights[peaks] < 1. - peak_probability)[0]

        if len(selected_peaks_indices) == 0:
            # Try with the deepest line
            sorted_peak_indices = np.argsort(fluxes[peaks])
            selected_peaks_indices = [sorted_peak_indices[0]]
        else:
            # Sort the interesting peaks from more to less deep
            sorted_peaks_indices = np.argsort(fluxes[peaks[selected_peaks_indices]])
            selected_peaks_indices = selected_peaks_indices[sorted_peaks_indices]
    else:
        # If no peaks found, just consider the deepest point and mark the base and top
        # as the limits of the whole data
        sorted_fluxes_indices = np.argsort(fluxes)
        peaks = sorted_fluxes_indices[0]
        base = 0
        top = len(xcoord) - 1
        selected_peaks_indices = [0]

    for i in np.asarray(selected_peaks_indices):
        #########################################################
        ####### 2nd degree polinomial fit to determine the peak
        #########################################################
        poly_step = 0.01
        # Use only 9 points for fitting (4 + 1 + 4)
        diff_base = peaks[i] - base[i]
        diff_top = top[i] - peaks[i]
        if diff_base > 4 and diff_top > 4:
            poly_base = peaks[i] - 4
            poly_top = peaks[i] + 4
        else:
            # There are less than 9 points but let's make sure that there are
            # the same number of point in each side to avoid asymetries that may
            # affect the fitting of the center
            if diff_base >= diff_top:
                poly_base = peaks[i] - diff_top
                poly_top = peaks[i] + diff_top
            elif diff_base < diff_top:
                poly_base = peaks[i] - diff_base
                poly_top = peaks[i] + diff_base
        p = np.poly1d(np.polyfit(xcoord[poly_base:poly_top+1], fluxes[poly_base:poly_top+1], 2))
        poly_vel = np.arange(xcoord[poly_base], xcoord[poly_top]+poly_step, poly_step)
        poly_ccf = p(poly_vel)
        mu = poly_vel[np.argmin(poly_ccf)]
        # Sometimes the polynomial fitting can give a point that it is not logical
        # (far away from the detected peak), so we do a validation check
        if mu < xcoord[peaks[i]-1] or mu > xcoord[peaks[i]+1]:
            mu = xcoord[peaks[i]]
            poly_step = xcoord[peaks[i]+1] - xcoord[peaks[i]] # Temporary just to the next iteration

        #########################################################
        ####### Gaussian/Voigt fit to determine other params.
        #########################################################
        # Models to fit
        gaussian_model = GaussianModel()
        voigt_model = VoigtModel()

        # Parameters estimators
        baseline = np.median(fluxes[base_points])
        A = fluxes[peaks[i]] - baseline
        sig = np.abs(xcoord[top[i]] - xcoord[base[i]])/3.0

        parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for j in np.arange(5)]
        parinfo[0]['value'] = 1.0 #fluxes[base[i]] # baseline # Continuum
        parinfo[0]['fixed'] = True
        #parinfo[0]['limited'] = [True, True]
        #parinfo[0]['limits'] = [fluxes[peaks[i]], 1.0]
        parinfo[1]['value'] = A # Only negative (absorption lines) and greater than the lowest point + 25%
        parinfo[1]['limited'] = [False, True]
        parinfo[1]['limits'] = [0., 0.]
        parinfo[2]['value'] = sig # Only positives (absorption lines)
        parinfo[2]['limited'] = [True, False]
        parinfo[2]['limits'] = [1.0e-10, 0.]
        parinfo[3]['value'] = mu # Peak only within the xcoord slice
        #parinfo[3]['fixed'] = True
        parinfo[3]['fixed'] = False
        parinfo[3]['limited'] = [True, True]
        #parinfo[3]['limits'] = [xcoord[base[i]], xcoord[top[i]]]
        #parinfo[3]['limits'] = [xcoord[peaks[i]-1], xcoord[peaks[i]+1]]
        parinfo[3]['limits'] = [mu-poly_step, mu+poly_step]

        # Only used by the voigt model (gamma):
        parinfo[4]['value'] = (xcoord[top[i]] - xcoord[base[i]])/2.0 # Only positives (not zero, otherwise its a gaussian) and small (for nm, it should be <= 0.01 aprox but I leave it in relative terms considering the spectrum slice)
        parinfo[4]['fixed'] = False
        parinfo[4]['limited'] = [True, True]
        parinfo[4]['limits'] = [0.0, xcoord[top[i]] - xcoord[base[i]]]

        f = fluxes[base[i]:top[i]+1]
        min_flux = np.min(f)
        # More weight to the deeper fluxes
        if min_flux < 0:
            weights = f + -1*(min_flux) + 0.01 # Above zero
            weights = np.min(weights) / weights
        else:
            weights = min_flux / f
        weights -= np.min(weights)
        weights = weights /np.max(weights)


        try:
            # Fit a gaussian and a voigt, but choose the one with the best fit
            if model in ['2nd order polynomial + auto fit', '2nd order polynomial + gaussian fit']:
                gaussian_model.fitData(xcoord[base[i]:top[i]+1], fluxes[base[i]:top[i]+1], parinfo=copy.deepcopy(parinfo[:4]), weights=weights)
                #gaussian_model.fitData(xcoord[base[i]:top[i]+1], fluxes[base[i]:top[i]+1], parinfo=copy.deepcopy(parinfo[:4]))
                rms_gaussian = np.sqrt(np.sum(np.power(gaussian_model.residuals(), 2)) / len(gaussian_model.residuals()))
            if model in ['2nd order polynomial + auto fit', '2nd order polynomial + voigt fit']:
                voigt_model.fitData(xcoord[base[i]:top[i]+1], fluxes[base[i]:top[i]+1], parinfo=copy.deepcopy(parinfo), weights=weights)
                #voigt_model.fitData(xcoord[base[i]:top[i]+1], fluxes[base[i]:top[i]+1], parinfo=copy.deepcopy(parinfo))
                rms_voigt = np.sqrt(np.sum(np.power(voigt_model.residuals(), 2)) / len(voigt_model.residuals()))

            if model == '2nd order polynomial + voigt fit' or (model == '2nd order polynomial + auto fit' and rms_gaussian > rms_voigt):
                final_model = voigt_model
                #logging.info("Voigt profile fitted with RMS %.5f" % (rms_voigt))
            else:
                final_model = gaussian_model
                #logging.info("Gaussian profile fitted with RMS %.5f" % (rms_gaussian))
#
            # Calculate velocity error based on:
            # Zucker 2003, "Cross-correlation and maximum-likelihood analysis: a new approach to combining cross-correlation functions"
            # http://adsabs.harvard.edu/abs/2003MNRAS.342.1291Z
            inverted_fluxes = 1-fluxes
            distance = xcoord[1] - xcoord[0]
            first_derivative = np.gradient(inverted_fluxes, distance)
            second_derivative = np.gradient(first_derivative, distance)
            ## Using the exact velocity, the resulting error are less coherents (i.e. sometimes you can get lower errors when using bigger steps):
            #second_derivative_peak = np.interp(final_model.mu(), xcoord, second_derivative)
            #inverted_fluxes_peak = final_model.mu()
            ## More coherent results:
            peak = xcoord.searchsorted(final_model.mu())
            inverted_fluxes_peak = inverted_fluxes[peak]
            second_derivative_peak = second_derivative[peak]
            if inverted_fluxes_peak == 0:
                inverted_fluxes_peak = 1e-10
            if second_derivative_peak == 0:
                second_derivative_peak = 1e-10
            sharpness = second_derivative_peak / inverted_fluxes_peak
            line_snr = np.power(inverted_fluxes_peak, 2) / (1 - np.power(inverted_fluxes_peak, 2))
            # Use abs instead of a simple '-1*' because sometime the result is negative and the sqrt cannot be calculated
            error = np.sqrt(np.abs(1 / (nbins * sharpness * line_snr)))

            final_model.set_emu(error)
            logging.info("Peak found at %.2f km/s (fitted at %.2f +/- %.2f km/s)" % (xcoord[peaks[i]], final_model.mu(), final_model.emu()))
            models.append(final_model)
        except Exception, e:
            print type(e), e.message


    return np.asarray(models)

def select_good_velocity_profile_models(models, ccf):
    """
    Select the modeled peaks that are not deeper than mean flux + 6*standard deviation
    unless it is the only detected peak.
    """
    accept = []
    if len(models) == 0:
        return accept

    xcoord = ccf['x']
    fluxes = ccf['y']

    ## We want to calculate the mean and standard deviation of the velocity profile
    ## but discounting the effect of the deepest detected lines:
    # Build the fluxes for the composite models
    line_fluxes = None
    for model in models:
        if line_fluxes is None:
            # first peak
            line_fluxes = model(xcoord)
            continue

        current_line_fluxes = model(xcoord)
        wfilter = np.where(line_fluxes > current_line_fluxes)[0]
        line_fluxes[wfilter] = current_line_fluxes[wfilter]
    ### Substract the line models conserving the base level
    if line_fluxes is not None:
        values = 1 + fluxes - line_fluxes
    else:
        values = fluxes
    ## Finally, calculate the mean and standard deviation
    check_mean = np.mean(values)
    check_std = np.std(values)
    for (i, model) in enumerate(models):
        # The first peak is always accepted
        if i == 0:
            accept.append(True)
            continue

        mu = model.mu()
        peak = model(mu)

        # Discard peak if it is not deeper than mean flux + 6*standard deviation
        limit = check_mean - 6*check_std
        if limit < 0.05 or peak >= limit:
            accept.append(False)
        else:
            accept.append(True)
    return np.asarray(accept)




############## [end] Radial velocity


#def refine_ew(spectrum, linemasks, continuum_model, resolution=None, margin=0.05):
    ## Normalize spectrum
    #spectrum = create_spectrum_structure(spectrum['waveobs'], spectrum['flux'])
    #spectrum['flux'] /= continuum_model(spectrum['waveobs'])
    #linemasks = linemasks.copy()
    #diff = []
    #for i in xrange(len(linemasks)):
        #line_region = linemasks[i:i+1]
        #wfilter = spectrum['waveobs'] >= line_region['wave_base'] - margin
        #wfilter = np.logical_and(wfilter, spectrum['waveobs'] <= line_region['wave_base'] + margin)
        #spectrum_window = spectrum[wfilter]
        ##--- Local continuum fit -------------------------------------------------------
        #model = "Polynomy" # Linear model (1 degree) for local continuum
        #degree = 1
        #nknots = None

        ## Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
        #order='median+max'
        #median_wave_range=0.1 # Bigger than for non local continuum fit
        #max_wave_range=1.0

        ## Fit locally in each individual segment
        #local_continuum_model = fit_continuum(spectrum_window, from_resolution=resolution, \
                                    #nknots=nknots, degree=degree,\
                                    #median_wave_range=median_wave_range, \
                                    #max_wave_range=max_wave_range, \
                                    #model=model, order=order, \
                                    #automatic_strong_line_detection=False)
        ##--- Fit lines -----------------------------------------------------------------
        ## Spectrum should be already radial velocity corrected
        #original_ew = linemasks['ew'][i]
        #line_region = fit_lines(line_region, spectrum, local_continuum_model, \
                                    #atomic_linelist = None, \
                                    #max_atomic_wave_diff = 0.005, \
                                    #telluric_linelist = None, \
                                    #smoothed_spectrum = None, \
                                    #check_derivatives = False, \
                                    #vel_telluric = 0, discard_gaussian=False, \
                                    #discard_voigt=True, free_mu=False)
        ##print original_ew - line_region['ew'], "=", original_ew, "-", line_region['ew'], "::", linemasks['ew'][i]
        #diff.append(original_ew - line_region['ew'][0])

    #return linemasks


def __moog_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None):
    """
    Saves a MOOG linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.

    linelist['waals'] ares spected to be in single gamma damping coefficient
    ABO theory is not supported and if waals is bigger than 0, it won't be considered
    """
    supported = linelist['moog_support'] == "T"
    linelist = linelist[supported]

    if linelist_filename is not None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)

    d0 = {} # dissociation energies from Turbospectrum's DATA/IRWIN_molecules_v15.1.dat
    # Molecule list from VALD linelist 300 - 1100 nm
    d0['C2 1'] = 6.210 # eV
    d0['CH 1'] = 3.465 # eV
    d0['CN 1'] = 7.724 # eV
    d0['CO 1'] = 11.092 # eV
    d0['MgH '] = 1.285 # eV
    d0['OH 1'] = 6.220 # eV
    d0['SiH '] = 4.337 # eV
    d0['TiO '] = 6.87 # eV

    with_ew = 'ew' in linelist.dtype.names
    line_ew = 0.
    out.write("wavelength species lower_state_eV loggf damping d0 equivalent_width comment\n")
    for line in linelist:
        #6690.261  24.0   3.888    -2.442  2.0    0.0    0.0    vald0.016 CrI may2008
        #6690.269  607.0  0.542    -3.334  0.00  7.73    0.0  ( 7, 2)Q12 12.5
        if with_ew:
            line_ew = line['ew']
        if line['molecule'] == 'T':
            if line['element'] in d0.keys():
                d0_value = d0[line['element']]
                waals = 0.
            else:
                logging.warn("Unknown dissociation energy for molecule '%s'" % (line['element']))
                continue
        else:
            d0_value = 0 # Dissociation energy [eV] (only for molecules)
            waals = line['waals_single_gamma_format']
        out.write("%10.3f%10s%10.3f%10.3f%10.2f%10.2f%10.2f %10s\n" \
                % (line['wave_A'], line['spectrum_moog_species'], line['lower_state_eV'], line['loggf'], waals, d0_value, line_ew, ""))
    out.close()
    return out.name

def __moog_barklem_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None):
    """
    Saves a MOOG linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.

    linelist['waals'] ares spected to be in single gamma damping coefficient
    ABO theory is not supported and if waals is bigger than 0, it won't be considered
    """
    supported = linelist['moog_support'] == "T"
    linelist = linelist[supported]
    # Do not include molecules in Barklem.dat
    # - I tested it and MOOG does not use it for anything
    molecules = linelist['molecule'] == "T"
    linelist = linelist[~molecules]

    if linelist_filename is not None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)

    for line in linelist:
        # wavelength species waals  alpha     rad
        #  4800.649  26.0  -7.73    0.250     8.13E+07
        alpha = 0.00 # Omara theory (ABO) sigma.alpha (sigma was transformed to waals_single_gamma_format) and alpha is provided separately
        if line['spectrum_transition_type'] == "AO":
            alpha = line['waals'] % 1 # Decimal part
        out.write("%10.3f%10s%10.2f%10.3f%10.2E\n" \
                % (line['wave_A'], line['spectrum_moog_species'], line['waals_single_gamma_format'], alpha, line['turbospectrum_rad']))
    out.close()
    return out.name

def __synthe_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None):
    """
    Saves a Synthe linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.
    """
    supported = linelist['synthe_support'] == "T"
    #supported = np.logical_and(supported, linelist['spectrum_support'] == "T")
    #supported = np.logical_and(supported, linelist['width_support'] == "T")
    linelist = linelist[supported]
    #linelist = linelist[linelist['spectrum_support'] == 'True']
    molecules = linelist['molecule'] == "T"

    there_are_molecules = len(np.where(molecules)[0]) > 0

    if linelist_filename is not None:
        atomic_out = open(linelist_filename[0], "w")
    else:
        # Temporary file
        atomic_out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)

    molecules_out = []
    if there_are_molecules:
        iwidth_species = map(int, map(float, linelist[molecules]['width_species']))
        unique_iwidth_species = np.unique(iwidth_species)
        nfiles = len(unique_iwidth_species)
        for i in xrange(nfiles):
            molecules_out.append(tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir))

    for line in linelist[~molecules]:
        #300.0021 -7.261 26.00   60197.937  4.0 s6D3/5g[4]   26874.548  5.0 5Dsp3P z5F  7.19 -5.69 -7.85K94  0 0  0 0.000  0 0.000    0    0           1059 1399
        #300.0047 -0.617 24.01   99677.930  0.5 a3P)5s f4P   66354.830  1.5 a3P)4p y4P  8.87 -5.64 -7.67K88  0 0  0 0.000  0 0.000    0    0           2617 1671
        #300.0061 -1.463  8.01  232536.060  1.5 3d   4P     265859.000  0.5 5f  *4D     9.61  0.00  0.00NBS  0 0  0 0.000  0 0.000    0    0              0    0
        #300.0082 -2.495 23.00    9824.610  2.5 d3s2 a4P     43147.280  3.5 s4F         7.47 -4.87 -7.49K88  0 0 51-1.333 51-0.001  -35    0F2 -2z0    1550 1240     0
        alpha = 0.00 # Omara theory (ABO) sigma.alpha (sigma was transformed to waals_single_gamma_format) and alpha is provided separately
        if line['spectrum_transition_type'] == "AO":
            alpha = line['waals'] % 1 # Decimal part

        if line['element'] == 'H 1':
            stark = 0.0
            waals = 0.0
            nonLTE_first_level_index = line['stark'] # 14 non-LTE level index for first level   I2
            nonLTE_second_level_index = line['waals'] # 15 non-LTE level index for second level   I2
        else:
            stark = line['stark']
            waals = line['waals_single_gamma_format']
            nonLTE_first_level_index = 0.
            nonLTE_second_level_index = 0.

        atomic_out.write("%11.4f%7.3f%6s%12.3f%5.2f %10s%12.3f%5.2f %10s%6.2f%6.2f%6.2f%4s%2.0f%2.0f%3i 0.000%3i 0.000    0    %5i          %5i    0    %.3f\n" % \
                        (line['wave_nm'], line['loggf'], line['width_species'], \
                            line['lower_state_cm1'], line['lower_j'], "X", \
                            line['upper_state_cm1'], line['upper_j'], "X", \
                            line['rad'], stark, waals, \
                            line['reference_code'][:4], \
                            nonLTE_first_level_index, nonLTE_second_level_index, \
                            line['spectrum_synthe_isotope'], line['spectrum_synthe_isotope'], \
                            int(line['lande_lower']*1000), int(line['lande_upper']*1000), \
                            alpha, \
                        )
                    )

    atomic_out.close()
    if not there_are_molecules:
        return (atomic_out.name, None)
    else:
        for i, iwidth in enumerate(unique_iwidth_species):
            sublinelist = linelist[molecules][iwidth_species == iwidth]
            for line in sublinelist:
                #256.2855 -8.611  1.5    17.770  2.5 -39025.059 106X00F1   C03F2   12 706.000
                molecules_out[i].write("%10.4f%7.3f%5.1f%10.3f%5.1f%11.3f%4iXXXXX   XXXXX  %3i%8.3f\n" % \
                                (line['wave_nm'], line['loggf'], \
                                    line['lower_j'], line['lower_state_cm1'], \
                                    line['upper_j'], line['upper_state_cm1'], \
                                    int(float(line['width_species'])), \
                                    #0, \
                                    line['spectrum_synthe_isotope'], \
                                    line['rad']*100
                                )
                            )
            molecules_out[i].close()

        molecules_out_name = []
        for out in molecules_out:
            molecules_out_name.append(out.name)
        return (atomic_out.name, molecules_out_name)

def __turbospectrum_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None):
    """
    Saves a TurboSpectrum linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.
    """
    supported = linelist['turbospectrum_support'] == "T"
    linelist = linelist[supported]

    if linelist_filename is not None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)

    # There can be different species that correspond to different isotopes
    #   Nd II: 60.142 and 60.143
    # And different ionization with the same species:
    #   26.000: Fe I and Fe II
    # thus we have to separate them looking at 'species' and 'element
    with_ew = 'ew' in linelist.dtype.names
    line_ew = 0.
    for species in  np.unique(linelist['turbospectrum_species']):
        sublinelist = linelist[linelist['turbospectrum_species'] == species]
        for element in np.unique(sublinelist['element']):
            #print element
            subsublinelist = sublinelist[sublinelist['element'] == element]
            molecule = subsublinelist['molecule'][0] == "T"
            # '   3.000            '    1        28
            # 'Li I   '
            out.write("'%20s' %4s %9i\n" % (species, subsublinelist['ion'][0], len(subsublinelist)))
            out.write("'%-7s'\n" % (element))
            for line in subsublinelist:
                # Do not use error or turbospectrum will calculate 3 times the abundance for a single line
                #line_ew_err = line['ew_err']
                line_ew_err = 0.0
                if with_ew:
                    line_ew = line['ew']
                if not molecule:
                    # 4602.826  1.848 -0.613 2006.342    4.0  6.61E+07 'p' 'd'   0.0    1.0 'Li I LS:1s2.2p 2P* LS:1s2.4d 2D'
                    out.write("%10.3f %9.5f %6.3f %8.3f %6.1f %9.2E '%s' '%s' %5.1f %6.1f '%s'\n" % (line['wave_A'], line['lower_state_eV'], \
                                                    line['loggf'], line['turbospectrum_fdamp'], line['upper_g'], \
                                                    line['turbospectrum_rad'], line['lower_orbital_type'], \
                                                    line['upper_orbital_type'], line_ew, line_ew_err, ""))
                else:
                    # 5152.971   1.08000 -3.516    0.000   14.0  2.00E+08 'X' 'X'   0.0    1.0 'X' 'X' 0.0 1.0 'pP2(7.5) 1 6.5 e 2 - 3 7.5 e 8 '
                    out.write("%10.3f %9.5f %6.3f %8.3f %6.1f %9.2E '%s' '%s' %5.1f %6.1f  '%s' '%s' %5.1f %6.1f  '%s'\n" % (line['wave_A'], line['lower_state_eV'], \
                                                    line['loggf'], line['turbospectrum_fdamp'], line['upper_g'], \
                                                    line['turbospectrum_rad'], \
                                                    line['lower_orbital_type'], line['upper_orbital_type'], line_ew, line_ew_err, \
                                                    line['lower_orbital_type'], line['upper_orbital_type'], line_ew, line_ew_err, \
                                                    ""))
    out.close()
    return out.name


def update_ew_with_ares(spectrum, linelist, rejt="0.995", tmp_dir=None, verbose=0):
    """
        rejt="0.995"
        rejt="100" # SNR
        rejt="3;5764,5766,6047,6052,6068,6076"
    """
    if not is_ares_support_enabled():
        raise Exception("ARES support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    ares_dir = ispec_dir + "/synthesizer/ARES/"
    ares_executable = ares_dir + "bin/ARES"

    tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
    ares_conf_file = tmp_execution_dir + "/mine.opt"
    linelist_file = tmp_execution_dir + "/linelist.dat"

    linelist = linelist.copy()
    linelist.sort(order=['wave_A'])
    linelist_filename = write_atomic_linelist(linelist, linelist_filename=linelist_file, code="moog", tmp_dir=tmp_dir)

    # ARES2 supports ASCII spectra with non-regular sampling but it seems to be buggy
    # it is better to homogeneize spectra by resampling if needed:
    spectrum_file = tmp_execution_dir + "/spectrum.fits"
    tmp_spectrum = spectrum.copy()
    tmp_spectrum.sort(order=['waveobs'])
    tmp_spectrum['waveobs'] *= 10. # Ares requires Amstrongs and not nm
    wavelengths = np.arange(np.min(tmp_spectrum['waveobs']), np.max(tmp_spectrum['waveobs']), 0.01)
    tmp_spectrum = resample_spectrum(tmp_spectrum, wavelengths, method="linear", zero_edges=True)
    tmp_spectrum['err'] = 0.
    write_spectrum(tmp_spectrum, spectrum_file)

    # Append microturbulence, solar abundances and metallicity
    ares_conf = open(ares_conf_file, "a")
    ares_conf.write("specfits='spectrum.fits'\n")
    #ares_conf.write("specfits='spectrum.txt'\n")
    ares_conf.write("readlinedat='linelist.dat'\n")
    ares_conf.write("fileout='results.ares'\n")
    ares_conf.write("lambdai=%.2f\n" % (np.min(tmp_spectrum['waveobs'])))
    ares_conf.write("lambdaf=%.2f\n" % (np.max(tmp_spectrum['waveobs'])))
    ares_conf.write("smoothder=4\n")
    ares_conf.write("space=3.0\n")
    ares_conf.write("rejt=%s\n" % (rejt))
    ares_conf.write("lineresol=0.1\n")
    ares_conf.write("miniline=0\n")
    ares_conf.write("plots_flag=0\n")
    ares_conf.close()

    previous_cwd = os.getcwd()
    os.chdir(tmp_execution_dir)

    command = ares_executable
    command_input = ""

    if verbose == 1:
        proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
    else:
        logging.info("Executing ARES...")
        proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    # wait for the process to terminate
    out, err = proc.communicate(input=command_input)
    errcode = proc.returncode


    try:
        data = np.loadtxt(tmp_execution_dir+"/results.ares")
        if len(data.shape) == 1:
            # If there is only one line, do a list or this function will fail
            data = (data,)
    except:
        print out
        sys.stdout.flush()
        raise Exception("ARES failed!")

    os.chdir(previous_cwd)
    shutil.rmtree(tmp_execution_dir)

    # ARES does not return EW for some lines sometimes, we should search for them
    # and add zeros for coherence
    ew_tmp = []
    ew_err_tmp = []
    i = 0
    j = 0
    while i < len(linelist) and j < len(data):
        equal_wave_a = np.abs(linelist['wave_A'][i] - data[j][0]) < 0.01
        if equal_wave_a:
            ew_tmp.append(data[j][4])
            ew_err_tmp.append(data[j][5])
            i += 1
            j += 1
        else:
            ew_tmp.append(np.nan)
            ew_err_tmp.append(np.nan)
            i += 1
    while i < len(linelist):
        ew_tmp.append(np.nan)
        ew_err_tmp.append(np.nan)
        i += 1

    ew = np.asarray(ew_tmp)
    ew_err = np.asarray(ew_err_tmp)
    linelist = linelist.copy()
    linelist['ew'] = ew
    linelist['ew_err'] = ew
    linelist['ewr'] = np.log10(linelist['ew'] / (1000.*linelist['wave_A']))

    return linelist

def van_der_Waals_ABO_to_single_gamma_format(gamvdw, atomic_mass, temperature=10000.):
    """
    This function converts the extended van der Waals in Paul Barklem's
    format sigma.alpha to a single gamma_vdw parameter computed for 10000 K.

    Written: 22-Feb-1999 N.Piskunov based on F77 code by Paul Barklem
    UH 2014-10-16: copied from preselect3.f90, svn Revision: 831
                   translated to IDL
    SBC 2015-10-9: Translated to python
    """
    if gamvdw <= 0:
        return gamvdw
    else:
        k = 1.380658e-23      # boltzmanns constant J/K
        m0 = 1.660540e-27     # unit atomic mass kg (Carbon 12 scale)
        a0 = 5.29177249e-11   # bohr radius m

        # Compute the broadening by hydrogen from cross-section data
        sigma = int(gamvdw)*a0*a0   # change to m^2
        alpha = gamvdw - int(gamvdw)

        # Compute the Gamma function of Z, this function valid over the range 1<Z<2
        # ie 0<ALPHA<2
        z =2. - alpha*0.5
        x = z - 1.0
        gammaf = 1.+(-0.5748646+(0.9512363+(-0.6998588+(0.4245549-0.1010678*x)*x)*x)*x)*x

        # Compute the halfwidth per unit perturber number density for vbar
        gvw = (4./np.pi)**(alpha*0.5)*gammaf*1.E4*sigma

        # Compute a single value for 10000K
        vbar = np.sqrt(8.*k*temperature/np.pi/m0*(1./1.008+1./atomic_mass))
        gvw = gvw * ((vbar/1.E4)**(1.-alpha))

        # Full width per perturber per cm^3
        gvw = gvw*1.E6*2.
        return np.log10(gvw)

