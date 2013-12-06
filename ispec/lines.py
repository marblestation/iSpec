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
#import ipdb
from common import *
from continuum import *
from lines import *
from spectrum import *
import matplotlib.pyplot as plt
from pymodelfit import UniformKnotSplineModel
from pymodelfit import UniformCDFKnotSplineModel
from pymodelfit import LinearModel
from mpfitmodels import GaussianModel
from mpfitmodels import VoigtModel
import log
import logging
import copy



########################################################################
## [START] LINE LISTS
########################################################################

def read_telluric_linelist(telluric_lines_file, minimum_depth=0.0):
    """
    Read telluric linelist.
    """
    telluric_lines = ascii.read(telluric_lines_file, delimiter="\t")._data

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
            return "Discard"
        else:
            symbol = str(molecules["symbol"][mfilter][0])
    else:
        symbol = str(chemical_elements["symbol"][tfilter][0])

    return symbol + " " + str(int(ion_state) + 1)

def __get_abund_rank(chemical_elements, species):
    """
    Convert species code form by the atomic number + "." + ionization state to element names type "Fe 1" or "Fe 2".
    Returns "Discard" if not found.
    """
    atomic_num, ion_state = species.split(".")

    tfilter = (chemical_elements['atomic_num'] == int(atomic_num))
    if len(chemical_elements["symbol"][tfilter]) == 0:
        return np.max(chemical_elements["solar_abund_rank"])
    else:
        return chemical_elements["solar_abund_rank"][tfilter][0]



def __get_upper_state(lower_state, wavelength):
    """
    Calculate upper excitation level from lower and wavelength.
    Units:

    * eV for lower excitation level
    * nm for wavelength
    """
    # Planck constant
    h = 6.626068 * 10e-34 # m^2 kg / s
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
    return value / 8065.73 # eV

def __eV_to_inverse_cm(value):
    """
    Units transformation from eV to cm^-1.
    """
    return value * 8065.73 # cm^-1



########################################################################
## [END] LINE LIST
########################################################################




def read_linelist_mask(mask_file):
    """
    Read linelist for building a mask and doing cross-correlation with the
    cross correlation function.

    The linelist should contain at least the columns wave_peak (in nm) and depth.
    """
    ccf_mask = ascii.read(mask_file)._data
    return ccf_mask

def read_chemical_elements(chemical_elements_file):
    """
    Read information about the chemical elements (periodic table)
    """
    chemical_elements = ascii.read(chemical_elements_file, delimiter="\t")._data
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
    molecules = ascii.read(molecules_file, delimiter="\t")._data
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
    """
    line_regions = np.array([tuple(line.rstrip('\r\n').split("\t")) for line in open(line_regions_filename,)][1:], dtype=[('wave_peak', float),('wave_base', float),('wave_top', float), ('note', '|S100')])

    if np.any(line_regions['wave_top'] - line_regions['wave_base'] <= 0):
        logging.error("Line regions where wave_base is equal or bigger than wave_top")
        raise Exception("Incompatible format")

    if np.any(line_regions['wave_top'] - line_regions['wave_peak'] < 0) or np.any(line_regions['wave_peak'] - line_regions['wave_base'] < 0):
        logging.error("Line regions where wave_peak is outside wave_base and wave_top")
        raise Exception("Incompatible format")

    return line_regions

def write_line_regions(line_regions, line_regions_filename):
    """
    Write line regions file with the following format:
    ::

        wave_peak       wave_base       wave_top        note
        480.8148        480.7970        480.8330        Fe 1
        496.2572        496.2400        496.2820        Fe 1
        499.2785        499.2610        499.2950
        505.8498        505.8348        505.8660        Fe 1
    """
    out = open(line_regions_filename, "w")
    out.write("wave_peak\twave_base\twave_top\tnote\n")
    out.write("\n".join(["\t".join(map(str, (line['wave_peak'], line['wave_base'], line['wave_top'], line['note']))) for line in line_regions]))
    out.close()

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
        weights = min_flux / f
    weights -= np.min(weights)
    weights = weights /np.max(weights)
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

def __create_linemasks_structure(num_peaks):
    linemasks = np.recarray((num_peaks, ), dtype=[ \
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
            ('rms', float), \
            ('telluric_wave_peak', float), ('telluric_fwhm', float), ('telluric_R', float), ('telluric_depth', float), \
            ('wave (nm)', '<f8'), ('wave (A)', '<f8'), \
            ('species', '|S10'), ('element', '|S4'), \
            ('lower state (cm^-1)', int), ('upper state (cm^-1)', int), \
            ('lower state (eV)', float), ('upper state (eV)', float), \
            ('log(gf)', '<f8'), ('fudge factor', float), ('transition type', '|S10'), \
            ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), \
            ('valid_theoretical_ew_depth', bool), ('theoretical_ew', float), ('theoretical_depth', float), \
            ('line_strength', float), \
            ('solar_abund_rank', int), \
            ('discarded', bool)])
    # Initialization
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
    linemasks['rms'] = 9999.0
    # Default values for line identification
    linemasks['wave (nm)'] = 0.0
    linemasks['wave (A)'] = 0.0
    linemasks["species"] = ""
    linemasks['element'] = ""
    linemasks['log(gf)'] = 0
    linemasks["lower state (cm^-1)"] = 0
    linemasks["upper state (cm^-1)"] = 0
    linemasks['lower state (eV)'] = 0
    linemasks['upper state (eV)'] = 0
    linemasks["fudge factor"] = 1.0
    linemasks["transition type"] = "99"
    linemasks["rad"] = 0
    linemasks["stark"] = 0
    linemasks["waals"] = 0
    linemasks["valid_theoretical_ew_depth"] = False
    linemasks["theoretical_ew"] = 0
    linemasks["theoretical_depth"] = 0
    linemasks["line_strength"] = 0
    linemasks["solar_abund_rank"] = 0
    return linemasks


def read_atomic_linelist(linelist_filename, chemical_elements, molecules):
    """
    Load a SPECTRUM format linelist and add complementary information.

    If chemical_elements and molecules is indicated, the field 'element' will be filled
    with the name of the element + ionization state (i.e. 'Fe 1')
    """
    try:
        has_depths = True
        raw_linelist = np.array([tuple(line.rstrip('\r\n').split()) for line in open(linelist_filename,)], \
                            dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), \
                            ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'), \
                            ('transition type', '|S10'), ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), \
                            ('note', '|S100'), \
                            ('valid_theoretical_ew_depth', '|S100'), \
                            ('theoretical_ew', '<f8'), \
                            ('theoretical_depth', '<f8')])
    except Exception, e:
        has_depths = False
        # Without valid_theoretical_ew_depth, ew and depth
        raw_linelist = np.array([tuple(line.rstrip('\r\n').split()) for line in open(linelist_filename,)], \
                            dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), \
                            ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'), \
                            ('transition type', '|S10'), ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), \
                            ('note', '|S100')])
    linelist = np.recarray((len(raw_linelist), ), dtype=[ \
            ('wave (nm)', '<f8'), ('wave (A)', '<f8'), \
            ('species', '|S10'), ('element', '|S4'), \
            ('lower state (cm^-1)', int), ('upper state (cm^-1)', int), \
            ('lower state (eV)', float), ('upper state (eV)', float), \
            ('log(gf)', '<f8'), ('fudge factor', float), ('transition type', '|S10'), \
            ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), \
            ('valid_theoretical_ew_depth', bool), ('theoretical_ew', float), ('theoretical_depth', float), \
            ('line_strength', float), \
            ('solar_abund_rank', int), \
            ('note', '|S100')])
    # Default values for line identification
    linelist['wave (nm)'] = raw_linelist['wave (A)'] / 10.
    linelist['wave (A)'] = raw_linelist['wave (A)']
    linelist["species"] = raw_linelist['species']
    if chemical_elements is not None and molecules is not None:
        for i, species in enumerate(raw_linelist['species']):
            linelist['element'][i] = __get_element(chemical_elements, molecules, species)
            linelist['solar_abund_rank'][i] = __get_abund_rank(chemical_elements, species)
    else:
        linelist['element'] = ""
    linelist["lower state (cm^-1)"] = raw_linelist['lower state (cm^-1)']
    linelist["upper state (cm^-1)"] = raw_linelist['upper state (cm^-1)']
    linelist['lower state (eV)'] = __inverse_cm_to_eV(raw_linelist['lower state (cm^-1)'])
    linelist['upper state (eV)'] = __inverse_cm_to_eV(raw_linelist['upper state (cm^-1)'])
    linelist['log(gf)'] = raw_linelist['log(gf)']
    linelist["fudge factor"] = raw_linelist['fudge factor']
    linelist["transition type"] = raw_linelist['transition type']
    linelist["rad"] = raw_linelist['rad']
    linelist["stark"] = raw_linelist['stark']
    linelist["waals"] = raw_linelist['waals']
    # Line strength is calculated by using a theoretical proxy: log(gf) - 5040/Teff
    # - Mucciarelli et al. 2013 - GALA: An automatic tool for the abundance analysis of stellar spectra
    # - Mucciarelli et al. 2010 - Microturbulent velocity from stellar spectra: a comparison between different approaches
    # - Temperature of reference used by Gaia ESO Survey: 3000
    teff_ref = 3000 # K
    linelist["line_strength"] = raw_linelist['log(gf)'] - (5040./teff_ref)*linelist['lower state (eV)']
    #line_strength_proxy = the oscillator strength - the Boltzmann factor (5040./Teff) * excitation potential
    linelist["note"] = raw_linelist['note']
    if has_depths:
        linelist["valid_theoretical_ew_depth"] = raw_linelist['valid_theoretical_ew_depth'] == "True"
        linelist["theoretical_ew"] = raw_linelist['theoretical_ew']
        linelist["theoretical_depth"] = raw_linelist['theoretical_depth']
    else:
        linelist["valid_theoretical_ew_depth"] = False
        linelist["theoretical_ew"] = 0.0
        linelist["theoretical_depth"] = 0.0
    return linelist

def write_atomic_linelist(linelist, linelist_filename=None):
    """
    Saves a SPECTRUM linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.
    """
    if linelist_filename is not None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    #4750.196  26.0  36078  57130  -3.662  1.0  GA  8.09  -4.61  -7.32  Fe_1  True  0.27  0.01
    out.write("\n".join(["  ".join(map(str, (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'], line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['note'], line['valid_theoretical_ew_depth'], line['theoretical_ew'], line['theoretical_depth']))) for line in linelist]))
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
        minimum_depth = np.min((1., np.max(depth)))

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
                                free_mu=True, \
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

def fit_lines(regions, spectrum, continuum_model, atomic_linelist, max_atomic_wave_diff=0.0005, telluric_linelist=None, vel_telluric=0.0, discard_gaussian = False, discard_voigt = False, check_derivatives=False, smoothed_spectrum=None, accepted_for_fitting=None, continuum_adjustment_margin=0.0, free_mu=False, frame=None):
    """
    Fits gaussians models in the specified line regions.
    * 'regions' should be an array with 'wave_base', 'wave_peak' and 'wave_top' columns.
    * If 'check_derivates' is True, the fit will be limited not to the region limits but to a refined smaller limits calculated from the second and third derivative (convex+concave regions) from 'smoothed_spectrum' if present or 'spectrum' if not.
    * If 'accepted_for_fitting' array is present, only those regions that are set to true
    will be fitted

    If atomic_linelist is specified, the
    lines will be cross-matched with the atomic information. In that case, the spectrum
    MUST be already radial velocity corrected.

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
    if not regions.dtype.fields.has_key('wave_base_fit'):
        # If it is not a complete regions (it has not been created with
        # __create_linemasks_structure and it only contains wave base, top and peak)
        # we create a complete one
        regions_tmp = regions
        regions = __create_linemasks_structure(total_regions)
        regions['wave_base'] = regions_tmp['wave_base']
        regions['wave_top'] = regions_tmp['wave_top']
        regions['wave_peak'] = regions_tmp['wave_peak']

        # Find index in spectrum for base, top and peak
        regions['base'] = 0
        regions['top'] = 0
        regions['peak'] = 0
        for i in np.arange(len(regions_tmp)):
            where_base = np.where(spectrum['waveobs'] >= regions_tmp['wave_base'][i])
            if len(where_base[0]) > 0:
                regions['base'][i] = where_base[0][0]
            where_top = np.where(spectrum['waveobs'] >= regions_tmp['wave_top'][i])
            if len(where_top[0]) > 0:
                regions['top'][i] = where_top[0][0]
            where_peak = np.where(spectrum['waveobs'] >= regions_tmp['wave_peak'][i])
            if len(where_peak[0]) > 0:
                regions['peak'][i] = where_peak[0][0]

    i = 0
    wfilter = spectrum['err'] > 0.
    if len(np.where(wfilter)[0]) > 0:
        snr = np.median(spectrum['flux'][wfilter] / spectrum['err'][wfilter])
    else:
        snr = None
    # Model: fit gaussian
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
                line_snr = np.power(inverted_fluxes_peak, 2) / (1 - np.power(inverted_fluxes_peak, 2))
                # Use abs instead of a simple '-1*' because sometime the result is negative and the sqrt cannot be calculated
                error = np.sqrt(np.abs(1 / (nbins * sharpness * line_snr)))
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
                    mean_flux = np.mean(spectrum['flux'][regions['base_fit'][i]:regions['top_fit'][i]+1])
                    mean_flux_continuum = np.mean(continuum_model(spectrum['waveobs'][regions['base_fit'][i]:regions['top_fit'][i]+1]))
                    diff_wavelength = spectrum['waveobs'][regions['top_fit'][i]] - spectrum['waveobs'][regions['base_fit'][i]]
                    regions['ew_err'][i] = np.sqrt(1 + mean_flux_continuum / mean_flux) * ((diff_wavelength - (regions['ew'][i]/10000.))/snr)
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
        regions = __fill_linemasks_with_atomic_data(regions, atomic_linelist, diff_limit=max_atomic_wave_diff, vel_atomic=0.0)
    if telluric_linelist is not None:
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


def __fill_linemasks_with_atomic_data(linemasks, atomic_linelist, diff_limit=0.0005, vel_atomic=0.0):
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
    #atomic_linelist.sort(order=['wave (nm)']) # Ascending

    clean_linemasks = linemasks[linemasks['wave_peak'] != 0]
    max_wave_peak = np.max(clean_linemasks['wave_peak'])
    min_wave_peak = np.min(clean_linemasks['wave_peak'])

    if atomic_linelist['wave (nm)'][0] > min_wave_peak or atomic_linelist['wave (nm)'][-1] < max_wave_peak:
        print "WARNING: The atomic linelist does not cover the whole linemask wavelength range"
        print "- Atomic line's range from", atomic_linelist['wave (nm)'][0], "to", atomic_linelist['wave (nm)'][-1], "nm"
        print "- Linemask range from", min_wave_peak, "to", max_wave_peak, "nm"

    wfilter = (atomic_linelist['wave (nm)'] >= min_wave_peak - diff_limit) & (atomic_linelist['wave (nm)'] <= max_wave_peak + diff_limit)
    atomic_linelist = atomic_linelist[wfilter]

    if len(atomic_linelist) > 0:
        for j in np.arange(len(linemasks)):
            diff = atomic_linelist['wave (nm)'] - linemasks['mu'][j]
            #diff = atomic_linelist['wave (nm)'] - linemasks['wave_peak'][j]

            # Find index of the nearest lines
            abs_diff = np.abs(diff)
            imin_diff = np.where(abs_diff <= diff_limit)[0]
            if len(imin_diff) == 0:
                continue
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
                linemasks['wave (nm)'][j] = atomic_linelist['wave (nm)'][i]
                linemasks['wave (A)'][j]= atomic_linelist['wave (A)'][i]
                linemasks["species"][j] = atomic_linelist['species'][i]
                linemasks['element'][j] = atomic_linelist['element'][i]
                linemasks['log(gf)'][j] = atomic_linelist['log(gf)'][i]
                linemasks["lower state (cm^-1)"][j] = atomic_linelist['lower state (cm^-1)'][i]
                linemasks["upper state (cm^-1)"][j] = atomic_linelist['upper state (cm^-1)'][i]
                linemasks['lower state (eV)'][j] = atomic_linelist['lower state (eV)'][i]
                linemasks['upper state (eV)'][j] = atomic_linelist['upper state (eV)'][i]
                linemasks["fudge factor"][j] = atomic_linelist['fudge factor'][i]
                linemasks["transition type"][j] = atomic_linelist['transition type'][i]
                linemasks["rad"][j] = atomic_linelist['rad'][i]
                linemasks["stark"][j] = atomic_linelist['stark'][i]
                linemasks["waals"][j] = atomic_linelist['waals'][i]
                linemasks["valid_theoretical_ew_depth"][j] = atomic_linelist['valid_theoretical_ew_depth'][i]
                linemasks["theoretical_ew"][j] = atomic_linelist['theoretical_ew'][i]
                linemasks["theoretical_depth"][j] = atomic_linelist['theoretical_depth'][i]
                linemasks["line_strength"][j] = atomic_linelist['line_strength'][i]
                linemasks["solar_abund_rank"][j] = atomic_linelist['solar_abund_rank'][i]

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
    linemasks.sort(order=['wave_peak'])
    for i, line in enumerate(linemasks[:-1]):
        if line['wave_top'] > linemasks['wave_base'][i+1]:
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


def __sampling_uniform_in_velocity(wave_base, wave_top, velocity_step):
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
    # Speed of light in km/s
    c = 299792.4580
    #c = 299792458.0

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

    waveobs = __sampling_uniform_in_velocity(np.min(spectrum['waveobs']), np.max(spectrum['waveobs']), velocity_step)
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



def cross_correlate_with_mask(spectrum, linelist, lower_velocity_limit=-200, upper_velocity_limit=200, velocity_step=1.0, mask_size=None, mask_depth=0.01, fourier=False, only_one_peak=False, model='2nd order polynomial + gaussian fit', frame=None):
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
                only_one_peak=only_one_peak, peak_probability=0.75, model=model, \
                frame=None)

def cross_correlate_with_template(spectrum, template, lower_velocity_limit=-200, upper_velocity_limit=200, velocity_step=1.0, fourier=False, only_one_peak=False, model='2nd order polynomial + gaussian fit', frame=None):
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
            only_one_peak=only_one_peak, peak_probability=0.75, model=model, \
            frame=None)

def __cross_correlate(spectrum, linelist=None, template=None, lower_velocity_limit = -200, upper_velocity_limit = 200, velocity_step=1.0, mask_size=2.0, mask_depth=0.01, fourier=False, only_one_peak=False, peak_probability=0.75, model='2nd order polynomial + gaussian fit', frame=None):
    ccf, nbins = __build_velocity_profile(spectrum, \
            linelist = linelist, template = template, \
            lower_velocity_limit = lower_velocity_limit, upper_velocity_limit = upper_velocity_limit, \
            velocity_step=velocity_step, \
            mask_size=mask_size, mask_depth=mask_depth, \
            fourier=fourier, frame=frame)

    models = __modelize_velocity_profile(ccf, nbins, only_one_peak=only_one_peak, \
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


def __modelize_velocity_profile(ccf, nbins, only_one_peak=False, peak_probability=0.75, model='2nd order polynomial + gaussian fit'):
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

    # Fit continuum to try to remove trends that may affect the peak detection
    #try:
        ##nknots = 2
        ##continuum_model = UniformCDFKnotSplineModel(nknots)
        ##continuum_model.fitData(xcoord[base_points], fluxes[base_points])
        ##continuum = continuum_model(xcoord)

        #continuum_model = LinearModel()
        #continuum_model.fitData(xcoord[base_points], fluxes[base_points])
        #continuum = continuum_model(xcoord)
        #if not np.any(continuum == 0) and not np.any(np.isnan(continuum)):
            #smoothed_fluxes /= continuum
            #fluxes /= continuum
            #errors /= continuum
        #else:
            #raise Exception()
    #except Exception:
        #logging.warn("Velocity profile cannot be normalized")

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
            #import pudb
            #pudb.set_trace()
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
        ## TODO: REMOVE Ivyes code
        #from ccanalyzer import analyzer

        ##rv = xcoord[base[i]:top[i]+1]
        ##cc = fluxes[base[i]:top[i]+1]
        #rv = xcoord
        ##depth = np.abs(np.max(template['flux']) - template['flux'])
        #cc = np.min(fluxes) - fluxes
        #aa = analyzer(rv, cc, rv.shape[0])
        ##print mu, aa.getResult(),aa.getError()
        #mu = aa.getResult()

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

    ##import pudb
    ##pudb.set_trace()
    #return linemasks
