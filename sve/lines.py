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
import asciitable
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions
import cPickle as pickle
import gzip
import os
#import ipdb
from common import *
from continuum import *
from radial_velocity import *
from convolve import *
import matplotlib.pyplot as plt
from pymodelfit import UniformKnotSplineModel
from mpfitmodels import GaussianModel
from mpfitmodels import VoigtModel
import log
import logging

########################################################################
## [START] LINE LISTS
########################################################################
# Load a VALD linelist (short format) and filter it
# - minimum_depth: filter out all the lines with a depth less than this percentage (0.05 = 5%)
# - data_end can be negative in order to ignore the last nth rows (VALD usually
#            adds references to the end of the file that should be ignored)
def read_VALD_linelist(vald_file, minimum_depth=0.0, data_end=None):
    # Original VALD linelist
    if data_end == None:
        vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "rad", "stark", "waals", "factor", "depth", "Reference"], exclude_names=["Vmic (km/s)", "factor", "Reference"], guess=False)
    else:
        vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, data_end=data_end, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "rad", "stark", "waals", "factor", "depth", "Reference"], exclude_names=["Vmic (km/s)", "factor", "Reference"], guess=False)

    ## Convert wavelengths from armstrong to nm
    #vald['wave (A)'] = vald['wave (A)'] / 10.0
    #vald = rfn.rename_fields(vald, {'wave (A)':'wave_peak',})
    vald = rfn.append_fields(vald, "wave_peak", dtypes=float, data=vald['wave (A)'] / 10.0)
    vald.sort(order=['wave_peak'])

    if minimum_depth <= 0.0:
        return vald
    else:
        # Filter
        vald_limited = vald[vald['depth'] >= minimum_depth]
        return vald_limited


def read_telluric_linelist(telluric_lines_file, minimum_depth=0.0):
    # Original VALD linelist
    telluric_lines = asciitable.read(telluric_lines_file, delimiter="\t")

    # Convert string to bool
    telluric_lines = rfn.rename_fields(telluric_lines, {'discarded':'discarded_string',})

    ## Detect duplicates
    # Sort by wave_peak and descending depth (for that we create a temporary field)
    telluric_lines = rfn.append_fields(telluric_lines, "discarded", dtypes=bool, data=(telluric_lines['discarded_string'] == "True"))

    telluric_lines = rfn.drop_fields(telluric_lines, ['discarded_string'])


    if minimum_depth <= 0.0:
        return telluric_lines
    else:
        # Filter
        telluric_lines_limited = telluric_lines[telluric_lines['depth'] >= minimum_depth]
        return telluric_lines_limited

#from lines import *
#l = "input/linelists/original/lumba-gustafsson.lin"
#o = "input/linelists/lumba-gustafsson.lin"
# Original VALD linelist
#vald_linelist = read_VALD_linelist(l, minimum_depth=0.0, data_end=None)
#linelist = VALD_to_SPECTRUM_format(vald_linelist)
# Filter discarded:
#linelist = linelist[linelist['species'] != "Discard"]
#asciitable.write(linelist, output=o, Writer=asciitable.FixedWidthNoHeader, delimiter=None, bookend=False, formats={'wave (A)': '%4.3f', })

# Convert element names type "Fe 1" or "Fe 2" to species code form by the atomic number + "." + ionization state
# Returns "Discard" if not found
def get_specie(chemical_elements, molecules, element_name):
    element = element_name.split() # Do not specify " " to avoid problems with elements with double spaces like "V  1"

    # Element not present or with a bad format, skip
    if element_name == "" or len(element) != 2:
        return "Discard"

    symbol = element[0]
    try:
        element.remove('') # Make sure there are not additional spaces between the symbol and the ionization state
        element.remove('')
        element.remove('')
    except ValueError as e:
        pass
    ionization = str(int(element[1]) - 1)

    tfilter = (chemical_elements['symbol'] == symbol)
    if len(chemical_elements[tfilter]["atomic_num"]) == 0:
        # Symbol not found, maybe it is a molecule
        mfilter = (molecules['symbol'] == symbol)
        if len(molecules[mfilter]["atomic_num"]) == 0:
            specie = "Discard"
        else:
            specie = str(molecules[mfilter]["atomic_num"][0]) + "." + ionization
    else:
        specie = str(chemical_elements[tfilter]["atomic_num"][0]) + "." + ionization
    return specie


## Calculate upper exciation level from lower and wavelength
# Units: eV for lower excitation level
#        nm for wavelength
def get_upper_state(lower_state, wavelength):
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

# Units transformation from eV to cm^-1
def eV_to_inverse_cm(value):
    return value * 8065.73 # cm^-1


# Convert a VALD linelist (short format) to a format that can be used with SPECTRUM
def VALD_to_SPECTRUM_format(vald_linelist):
    # Periodic table
    chemical_elements = asciitable.read("input/abundances/chemical_elements_symbols.dat", delimiter="\t")
    # Some molecular symbols
    # - For diatomic molecules, the atomic_num specifies the atomic makeup of the molecule.
    #   Thus, H2 is 101.0, the two ``1''s referring to the two hydrogens, CH is 106.0,
    #   CO 608.0, MgH 112.0, TiO 822.0, etc.
    # - The lightest element always comes first in the code, so that 608.0 cannot be
    #   confused with NdO, which would be written 860.0.
    molecules = asciitable.read("input/abundances/molecular_symbols.dat", delimiter="\t")

    # Prepare resulting structure
    linelist = np.recarray((len(vald_linelist), ), dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'),('transition type', '|S10'), ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), ('note', '|S100')])
    linelist['species'] = ""
    linelist['fudge factor'] = 1.0
    linelist['transition type'] = "99"
    linelist['note'] = ""

    linelist['wave (A)'] = vald_linelist['wave (A)']
    linelist['upper state (cm^-1)'] = (eV_to_inverse_cm(get_upper_state(vald_linelist['lower state (eV)'], vald_linelist[ "wave (A)"] / 10.))).astype(int)
    linelist['lower state (cm^-1)'] = (eV_to_inverse_cm(vald_linelist['lower state (eV)'])).astype(int)
    linelist['log(gf)'] = vald_linelist['log(gf)']
    #linelist['transition type'] = "99"
    linelist['transition type'] = "GA"
    linelist['rad'] = vald_linelist['rad']
    linelist['stark'] = vald_linelist['stark']
    linelist['waals'] = vald_linelist['waals']
    # Van der Waals should be zero or negative
    fwaals = vald_linelist['waals'] > 0
    linelist['waals'][fwaals] = 0
    # TODO: Consider AO
    #if ():
        ## Anstee-O'Mara theory
        #linelist['transition type'] = "AO"
        #linelist['sig.alpha'] = vald_linelist['Waals']
        #linelist['rad'] = vald_linelist['Waals']
        #linelist['rad'] = 0
        #linelist['stark'] = -999 # Mark to remove manually after
        #linelist['waals'] = -999 # Mark to remove manually after

    i = 0
    for line in vald_linelist:
        linelist[i]['species'] = get_specie(chemical_elements, molecules, line["element"])
        linelist[i]['note'] = "_".join(line["element"].split())
        i += 1

    return linelist



# Convert a VALD linelist (short format) to a format that can be used to measure radial velocities
# and select the top N deepest lines every wave_step armstrongs (1 nm).
# - data_end can be negative in order to ignore the last nth rows (VALD usually
#            adds references to the end of the file that should be ignored)
def VALD_top_N_to_RV_format(vald_file, output_file, top=1, wave_step=10, data_end=None):
    #vald_limited = read_VALD_linelist(vald_file, minimum_depth=0, data_end=data_end)
    vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "Rad", "Stark", "Waals", "factor", "depth", "Reference"], guess=False)
    # Filter
    vald_limited = vald[vald['depth'] >= 0.0]

    wave_base = np.min(vald_limited["wave (A)"])
    wave_top = np.max(vald_limited["wave (A)"])

    # Boolean array initialized to False
    #selected = (np.zeros(top*np.ceil((wave_top - wave_base) / wave_step)) == 1)
    selected = []

    wave_current = wave_base
    i = 0
    # For each segment
    while wave_current < wave_top:
        wfilter = (vald_limited["wave (A)"] >= wave_current) & (vald_limited["wave (A)"] < wave_current + wave_step)
        vald_filtered = vald_limited[wfilter]
        vald_filtered.sort(order="depth")
        # Select the top N deepest lines
        for j in np.arange(top):
            pos = -1*(j+1)
            selected.append(vald_filtered[pos])
        wave_current += wave_step
        i += top

    selected = np.array(selected, dtype=vald_limited.dtype)

    asciitable.write(selected, output=output_file, delimiter=",", quotechar="'")

#VALD_top_N_to_RV_format("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/rv/VALD.300_1100nm.rv.lst2", top=1, wave_step=10, data_end=None)


def filter_telluric_lines(linelist_telluric, spectra, velocity_lower_limit, velocity_upper_limit):
    # Light speed in vacuum
    c = 299792458.0 # m/s

    ## Select telluric lines of interest
    # Limit to region of interest
    wmin = spectra['waveobs'][0]
    wmax = spectra['waveobs'][-1]
    delta_wmin = wmin * (velocity_lower_limit / (c/1000.0))
    delta_wmax = wmax * (velocity_upper_limit / (c/1000.0))
    wfilter = (linelist_telluric['wave_peak'] <= wmax + delta_wmax) & (linelist_telluric['wave_peak'] >= wmin + delta_wmin)
    linelist = linelist_telluric[wfilter]
    # Discard not fitted lines
    rfilter = linelist['rms'] == 9999
    linelist = linelist[~rfilter]
    # Discard too deep or too small lines
    rfilter = (linelist['depth'] <= 0.9) & (linelist['depth'] >= 0.01)
    linelist = linelist[rfilter]
    # Discard outliers FWHM in km/s (which is not wavelength dependent)
    telluric_fwhm = (c / (linelist['wave_peak'] / linelist['fwhm'])) / 1000.0 # km/s
    fwhm_selected, fwhm_selected_filter = sigma_clipping(telluric_fwhm, meanfunc=np.median)
    linelist = linelist[fwhm_selected_filter]
    return linelist


########################################################################
## [END] LINE LIST
########################################################################

# Fits a gaussian at a given wavelength location using a fitted continuum model
# - For absorption lines, it will alway be true:
#      model.A() < 0 and model.sig() > 0
#   The oposite indicates a potential bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_gaussian(spectra_slice, continuum_model, mu, sig=None, A=None):
    model = GaussianModel()
    x = spectra_slice['waveobs']
    y = spectra_slice['flux']
    min_flux = np.min(y)

    # Parameters estimators
    baseline = np.median(continuum_model(spectra_slice['waveobs']))
    if A == None:
        A = min_flux - baseline
    if sig == None:
        sig = (x[-1] - x[0])/3.0


    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(4)]
    parinfo[0]['value'] = baseline # Continuum
    parinfo[0]['fixed'] = True
    parinfo[1]['value'] = A # Only negative (absorption lines) and greater than the lowest point + 25%
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [(min_flux-baseline) * 1.25, 0.]
    parinfo[2]['value'] = sig # Only positives (absorption lines) and lower than the spectra slice
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [0., x[-1] - x[0]]
    parinfo[3]['value'] = mu # Peak only within the spectra slice
    parinfo[3]['limited'] = [True, True]
    #parinfo[3]['limits'] = [x[0], x[-1]]
    parinfo[3]['limits'] = [mu - 0.005, mu + 0.005]

    # If there are only 3 data point, fix 'mu'
    # - if not, the fit will fail with an exception because there are not
    #   more data points than parameters
    if len(spectra_slice) == 3:
        parinfo[3]['fixed'] = True

    model.fitData(x, y, parinfo=parinfo)

    return model

# Fits a voigt at a given wavelength location using a fitted continuum model
# - For absorption lines, it will alway be true:
#      model.A() < 0 and model.sig() > 0
#   The oposite indicates a potential bad fit
# - For absorption lines, model.gamma() < 0 indicates strange wings and probably a bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_voigt(spectra_slice, continuum_model, mu, sig=None, A=None, gamma=None):
    model = VoigtModel()
    x = spectra_slice['waveobs']
    y = spectra_slice['flux']
    min_flux = np.min(y)

    # Parameters estimators
    baseline = np.median(continuum_model(spectra_slice['waveobs']))
    if A == None:
        A = min_flux - baseline
    if sig == None:
        sig = (x[-1] - x[0])/3.0
    if gamma == None:
        gamma = (x[-1] - x[0])/2.0

    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(5)]
    parinfo[0]['value'] = baseline # Continuum
    parinfo[0]['fixed'] = True
    parinfo[1]['value'] = A # Only negative (absorption lines) and greater than the lowest point + 25%
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [(min_flux-baseline) * 1.25, 0.]
    parinfo[2]['value'] = sig # Only positives (absorption lines) and lower than the spectra slice
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [0., x[-1] - x[0]]
    parinfo[3]['value'] = mu # Peak only within the spectra slice
    parinfo[3]['limited'] = [True, True]
    #parinfo[3]['limits'] = [x[0], x[-1]]
    parinfo[3]['limits'] = [mu - 0.005, mu + 0.005]
    parinfo[4]['value'] = gamma # Only positives (not zero, otherwise its a gaussian) and small (for nm, it should be <= 0.01 aprox but I leave it in relative terms considering the spectra slice)
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [0.001, x[-1] - x[0]]

    # If there are only 4 data point, fix 'mu'
    # - if not, the fit will fail with an exception because there are not
    #   more data points than parameters
    if len(spectra_slice) == 4:
        parinfo[3]['fixed'] = True

    model.fitData(x, y, parinfo=parinfo)

    return model


# Fits a gaussian and a voigt at a given wavelength location using a fitted continuum model
# - It selects the best fitted model (gaussian or voigt) unless one of them is disabled by
#   the discard_gaussian or discard_voigt argument
# - For absorption lines, it will alway be true:
#      model.A() < 0 and model.sig() > 0
#   The oposite indicates a potential bad fit
# - For absorption lines fitted with voigt, model.gamma() < 0 indicates strange wings and probably a bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_line(spectra_slice, continuum_model, mu, sig=None, A=None, gamma=None, discard_gaussian = False, discard_voigt = False):
    if not discard_gaussian:
        # Default values for failed fit:
        rms_gaussian = 9999.0
        discard_gaussian = True

        # If there are more data points than parameters (if not, the fit will fail with an exception)
        # - 2 free parameters: A, sig
        # - 1 fix parameter: mu (but it will set to free if there is enough data)
        if len(spectra_slice) > 2:
            try:
                gaussian_model = fit_gaussian(spectra_slice, continuum_model, mu, sig=sig, A=A)

                residuals = gaussian_model.residuals()
                rms_gaussian = np.sqrt(np.sum(np.power(residuals, 2)) / len(residuals))
                discard_gaussian = False
            except Exception as e:
                print e.message

    if not discard_voigt:
        # Default values for failed fit:
        rms_voigt = 9999.0
        discard_voigt = True

        # If there are more data points than parameters (if not, the fit will fail with an excepteion)
        # - 3 free parameters: A, sig, gamma
        # - 1 fix parameter: mu (but it will set to free if there is enough data)
        if len(spectra_slice) > 3:
            try:
                voigt_model = fit_voigt(spectra_slice, continuum_model, mu, sig=sig, A=A, gamma=gamma)
                residuals = voigt_model.residuals()
                rms_voigt = np.sqrt(np.sum(np.power(residuals, 2)) / len(residuals))
                discard_voigt = False
            except Exception as e:
                print e.message

    if (not discard_gaussian and not discard_voigt and rms_gaussian <= rms_voigt) or (not discard_gaussian and discard_voigt):
        return gaussian_model, rms_gaussian
    elif (not discard_gaussian and not discard_voigt and rms_gaussian > rms_voigt) or (discard_gaussian and not discard_voigt):
        return voigt_model, rms_voigt
    else:
        raise Exception("Gaussian or Voigt fit for absorption line not possible.")

# Remove consecutive features (i.e. peaks or base points)
def remove_consecutives_features(features):
    duplicated_features = (np.abs(features[:-1] - features[1:]) == 1)
    duplicated_features = np.array([False] + duplicated_features.tolist())
    cleaned_features = features[~duplicated_features]
    return cleaned_features


def detect_false_positives_and_noise(spectra, linemasks):
    # Find
    # - Peaks that are less deep than it nearby base points (false positives)
    # - Peaks too close (2 or less positions) to the next and previous base point (noise)
    peaks = linemasks['peak']
    base_points = np.asarray([linemasks['base'][0]] +  linemasks['top'].tolist())

    ## First feature found in spectra: base point
    ## Last feature found in spectra: base point
    # Left
    peak_base_diff_left = spectra['flux'][base_points[:-1]] - spectra['flux'][peaks]
    peak_base_index_diff_left = base_points[:-1] - peaks
    # Right
    peak_base_diff_right = spectra['flux'][base_points[1:]] - spectra['flux'][peaks]
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

# Given a group of peaks and base_points with the following assumptions
# - base_points[i] < base_point[i+1]
# - peaks[i] < peaks[i+1]
# - base_points[i] < peaks[i] < base_points[i+1]
# The function returns peaks and base_points where:
# - The first and last feature is a base point
#     base_points[0] < peaks[0] < base_points[1] < ... < base_points[n-1] < peaks[n-1] < base_points[n]
#   where n = len(base_points)
# - len(base_points) = len(peaks) + 1
def assert_structure(xcoord, yvalues, peaks, base_points):
    # Limit the base_points array to the ones that are useful, considering that
    # the first and last peak are always removed
    first_wave_peak = xcoord[peaks][0]
    first_wave_base = xcoord[base_points][0]
    if first_wave_peak > first_wave_base:
        if len(base_points) - len(peaks) == 1:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: base point
            base_points = base_points
            peaks = peaks
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: peak
            # - Remove last peak
            #base_points = base_points
            #peaks = peaks[:-1]
            # - Add a base point (last point in the spectra)
            base_points = np.hstack((base_points, [len(xcoord)-1]))
            peaks = peaks
        else:
            raise Exception("This should not happen")
    else:
        if len(base_points) - len(peaks) == -1:
            ## First feature found in spectra: peak
            ## Last feature found in spectra: peak
            # - Remove first and last peaks
            #base_points = base_points
            #peaks = peaks[1:-1]
            # - Add two base points (first and last point in the spectra)
            base_points = np.hstack(([0], base_points))
            base_points = np.hstack((base_points, [len(xcoord)-1]))
            peaks = peaks
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: peak
            ## Last feature found in spectra: base point
            # - Remove first peak
            #base_points = base_points
            #peaks = peaks[1:]
            # - Add a base point (first point in the spectra)
            base_points = np.hstack(([0], base_points))
            peaks = peaks
        else:
            raise Exception("This should not happen")

    return peaks, base_points

# The function arguments should agree with:
# - The first and last feature is a base point
#     base_points[0] < peaks[0] < base_points[1] < ... < base_points[n-1] < peaks[n-1] < base_points[n]
#   where n = len(base_points)
# - len(base_points) = len(peaks) + 1
# It is recommended that the spectra is normalized (better results are obtained)
# It tries to fit a gaussian and a voigt model and selects the best unless one of them is disabled by
# the discard_gaussian or discard_voigt argument (if both of them are disabled, there will be
# no fit information)
# To save computation time, min and max depth can be indicated and all the lines out
# of this range will not be considered for fit process (save CPU time) although
# the rest of the information of the line will be conserved in the output
# Returns a complete structure with all the necessary information to
# determine if it is a line of interest
def generate_linemasks(spectra, peaks, base_points, continuum_model, minimum_depth=None, maximum_depth=None, smoothed_spectra=None, vald_linelist_file=None, telluric_linelist_file = None, discard_gaussian = False, discard_voigt = False, vel_atomic=0.0, vel_telluric=0.0, frame=None):
    print "NOTICE: This method can generate overflow warnings due to the Least Square Algorithm"
    print "        used for the fitting process, but they can be ignored."
    if smoothed_spectra == None:
        smoothed_spectra = spectra

    # Depth of the peak with respect to the total continuum in % over the total continuum
    # - In case that the peak is higher than the continuum, depth < 0
    continuum_at_peak = continuum_model(spectra['waveobs'][peaks])
    flux_at_peak = spectra['flux'][peaks]
    depth = ((continuum_at_peak - flux_at_peak) / continuum_at_peak)

    if minimum_depth == None:
        minimum_depth = np.min(depth)
    if maximum_depth == None:
        minimum_depth = np.max(depth)

    # To save computation time, min and max depth can be indicated and all the lines out
    # of this range will not be considered for fit process (save CPU time) although
    # the rest of the information of the line will be conserved in the output
    accepted_for_fitting = np.logical_and(depth >= minimum_depth, depth <= maximum_depth)

    last_reported_progress = -1
    if frame != None:
        frame.update_progress(0)

    num_peaks = len(peaks)
    linemasks = np.recarray((num_peaks, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('peak', int), ('base', int), ('top', int), ('depth', float), ('relative_depth', float), ('wave_base_fit', float), ('wave_top_fit', float), ('base_fit', int), ('top_fit', int), ('mu', float), ('sig', float), ('A', float), ('baseline', float), ('gamma', float), ('fwhm', float), ('fwhm_kms', float), ('R', float), ('depth_fit', float), ('relative_depth_fit', float), ('integrated_flux', float), ('ew', float), ('rms', float), ('VALD_wave_peak', float), ('element', '|S4'), ('lower state (eV)', float), ('log(gf)', float), ('telluric_wave_peak', float), ('telluric_fwhm', float), ('telluric_R', float), ('telluric_depth', float), ('solar_depth', float), ('discarded', bool), ('species', '|S10'), ('lower state (cm^-1)', int), ('upper state (cm^-1)', int), ('fudge factor', float), ('transition type', '|S10'), ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8')])

    linemasks['discarded'] = False
    # Line mask
    linemasks['wave_peak'] = spectra['waveobs'][peaks]
    linemasks['wave_base'] = spectra['waveobs'][base_points[:-1]]
    linemasks['wave_top'] = spectra['waveobs'][base_points[1:]]
    # Position in the spectra, usefull for noise detection
    linemasks['peak'] = peaks
    linemasks['base'] = base_points[:-1]
    linemasks['top'] = base_points[1:]
    linemasks['depth'] = depth
    # Relative depth is "peak - mean_base_point" with respect to the total continuum
    # - In case that the mean base point is higher than the continuum, relative_depth < 0
    # - relative_depth < depth is always true
    flux_from_top_base_point_to_continuum = np.abs(continuum_at_peak - np.mean([spectra['flux'][base_points[:-1]], spectra['flux'][base_points[1:]]]))
    linemasks['relative_depth'] = ((continuum_at_peak - (flux_at_peak + flux_from_top_base_point_to_continuum)) / continuum_at_peak)
    # Default values for line identification
    linemasks['VALD_wave_peak'] = 0
    linemasks['element'] = ""
    linemasks['lower state (eV)'] = 0
    linemasks['log(gf)'] = 0
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
    linemasks['ew'] = 0.0
    linemasks['rms'] = 9999.0
    linemasks["species"] = ""
    linemasks["lower state (cm^-1)"] = 0
    linemasks["upper state (cm^-1)"] = 0
    linemasks["fudge factor"] = 0
    linemasks["transition type"] = ""

    # To save computation time, exclude false positives and noise from the fitting process
    rejected_by_noise = detect_false_positives_and_noise(spectra, linemasks)
    accepted_for_fitting = np.logical_and(accepted_for_fitting, np.logical_not(rejected_by_noise))

    # Model: fit gaussian
    for i in np.arange(num_peaks):
        fitting_not_possible = False
        if accepted_for_fitting[i]:
            # Adjust edges
            new_base, new_top = improve_linemask_edges(smoothed_spectra['waveobs'], smoothed_spectra['flux'], linemasks['base'][i], linemasks['top'][i], linemasks['peak'][i])
            linemasks['base_fit'][i] = new_base
            linemasks['top_fit'][i] = new_top
            linemasks['wave_base_fit'][i] = spectra['waveobs'][new_base]
            linemasks['wave_top_fit'][i] = spectra['waveobs'][new_top]

            try:
                line_model, rms = fit_line(spectra[new_base:new_top+1], continuum_model, linemasks['wave_peak'][i], discard_gaussian = discard_gaussian, discard_voigt = discard_voigt)
                linemasks['mu'][i] = line_model.mu()
                linemasks['sig'][i] = line_model.sig()
                linemasks['A'][i] = line_model.A()
                linemasks['baseline'][i] = line_model.baseline()
                if type(line_model) == VoigtModel:
                    linemasks['gamma'][i] = line_model.gamma()
                else:
                    # The model is Gaussian, do not use 'gamma'
                    linemasks['gamma'][i] = 9999.0

                linemasks['fwhm'][i], linemasks['fwhm_kms'][i] = line_model.fwhm()
                linemasks['R'][i] = linemasks['mu'][i] / linemasks['fwhm'][i] # Resolving power

                # Depth of the peak with respect to the total continuum in % over the total continuum
                # - In case that the peak is higher than the continuum, depth < 0
                continuum = line_model.baseline()
                flux = line_model(line_model.mu())
                linemasks['depth_fit'][i] = ((continuum - flux) / continuum)
                # Relative depth is "peak - mean_base_point" with respect to the total continuum
                # - In case that the mean base point is higher than the continuum, relative_depth < 0
                # - relative_depth < depth is always true
                flux_from_top_base_point_to_continuum = np.abs(continuum - np.mean([spectra['flux'][base_points[i]], spectra['flux'][base_points[i+1]]]))
                linemasks['relative_depth_fit'][i] = ((continuum - (flux + flux_from_top_base_point_to_continuum)) / continuum)

                # Equivalent Width
                # - Include 99.97% of the gaussian area
                from_x = linemasks['mu'][i] - 3*linemasks['sig'][i]
                to_x = linemasks['mu'][i] + 3*linemasks['sig'][i]
                linemasks['integrated_flux'][i] = -1 * line_model.integrate(from_x, to_x)
                linemasks['ew'][i] = linemasks['integrated_flux'][i] / np.mean(continuum_model(spectra['waveobs'][new_base:new_top+1]))
                # RMS
                linemasks['rms'][i] = rms
            except Exception as e:
                #print "WARNING: Bad line fit (", i, ") - ", e.message
                fitting_not_possible = True


        current_work_progress = ((i*1.0)/num_peaks) * 100
        if report_progress(current_work_progress, last_reported_progress):
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)
            if frame != None:
                frame.update_progress(current_work_progress)

    if vald_linelist_file != None:
        linemasks = fill_with_VALD_info(linemasks, vald_linelist_file=vald_linelist_file, diff_limit=0.005, vel_atomic=vel_atomic)
    if telluric_linelist_file != None:
        linemasks = fill_with_telluric_info(linemasks, telluric_linelist_file=telluric_linelist_file, vel_telluric=vel_telluric)

    return linemasks

# Adds info to those linemasks that can be affected by tellurics
# Different nearby linemasks can be affected by the same telluric line
def fill_with_telluric_info(linemasks, telluric_linelist_file = "input/linelists/telluric/standard_atm_air_model.lst", vel_telluric=0.0):
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

    # Discard very small lines (lesser than 1% of the continuum)
    telluric_linelist = read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
    # The discarded flag is not a good one because there are clear lines mark as true (i.e. 628.0392 nm)
    #telluric_linelist = telluric_linelist[telluric_linelist['discarded'] != True]
    # It is better to clear the not fitted telluric lines:
    telluric_linelist = telluric_linelist[telluric_linelist['rms'] < 9999]
    telluric_linelist.sort(order=['wave_peak'])

    clean_linemasks = linemasks[linemasks['wave_peak'] != 0]
    max_wave_peak = np.max(clean_linemasks['wave_peak'])
    min_wave_peak = np.min(clean_linemasks['wave_peak'])

    if telluric_linelist['wave_peak'][0] > min_wave_peak or telluric_linelist['wave_peak'][-1] < max_wave_peak:
        print "WARNING: Telluric linelist does not cover the whole linemask wavelength range"
        print "- Telluric range from", telluric_linelist['wave_peak'][0], "to", telluric_linelist['wave_peak'][-1], "nm"
        print "- Linemask range from", min_wave_peak, "to", max_wave_peak, "nm"

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


# Cross-match linemasks with a VALD linelist in order to find
# the nearest lines and copy the information into the linemasks structure
def fill_with_VALD_info(linemasks, vald_linelist_file="input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", diff_limit=0.005, vel_atomic=0.0):
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


    ## Load original VALD linelist
    # Discard very small lines (lesser than 1% of the continuum)
    vald_linelist = read_VALD_linelist(vald_linelist_file, minimum_depth=0.01)

    # Sort by wave_peak and descending depth (for that we create a temporary field)
    vald_linelist = rfn.append_fields(vald_linelist, "reverse_depth", dtypes=float, data=1-vald_linelist['depth'])
    vald_linelist.sort(order=['wave_peak', 'reverse_depth'])
    vald_linelist = rfn.drop_fields(vald_linelist, ['reverse_depth'])

    clean_linemasks = linemasks[linemasks['wave_peak'] != 0]
    max_wave_peak = np.max(clean_linemasks['wave_peak'])
    min_wave_peak = np.min(clean_linemasks['wave_peak'])

    if vald_linelist['wave_peak'][0] > min_wave_peak or vald_linelist['wave_peak'][-1] < max_wave_peak:
        print "WARNING: VALD linelist does not cover the whole linemask wavelength range"
        print "- VALD range from", vald_linelist['wave_peak'][0], "to", vald_linelist['wave_peak'][-1], "nm"
        print "- Linemask range from", min_wave_peak, "to", max_wave_peak, "nm"

    wfilter = (vald_linelist['wave_peak'] >= min_wave_peak - diff_limit) & (vald_linelist['wave_peak'] <= max_wave_peak + diff_limit)
    vald_linelist = vald_linelist[wfilter]

    if len(vald_linelist) > 0:
        for j in np.arange(len(linemasks)):
            # Find index of the nearest line
            diff = vald_linelist['wave_peak'] - linemasks['mu'][j]
            #diff = vald_linelist['wave_peak'] - linemasks['wave_peak'][j]
            abs_diff = np.abs(diff)
            i = np.argmin(abs_diff)
            # Save the information
            if abs_diff[i] <= diff_limit:
                linemasks["VALD_wave_peak"][j] = vald_linelist["wave_peak"][i]
                linemasks["element"][j] = vald_linelist["element"][i]
                linemasks["lower state (eV)"][j] = vald_linelist["lower state (eV)"][i]
                linemasks["log(gf)"][j] = vald_linelist["log(gf)"][i]
                linemasks["solar_depth"][j] = vald_linelist["depth"][i]
                linemasks["rad"][j] = vald_linelist["rad"][i]
                linemasks["stark"][j] = vald_linelist["stark"][i]
                linemasks["waals"][j] = vald_linelist["waals"][i]

    if vel_atomic != 0:
        linemasks['wave_peak'] = original_wave_peak
        linemasks['mu'] = original_mu

    return linemasks

# Works better with a smoothed spectra (i.e. convolved using 2*resolution)
def find_peaks_and_base_points(xcoord, yvalues):
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
        peaks = remove_consecutives_features(peaks)
        base_points = remove_consecutives_features(base_points)

        if not (len(peaks) - len(base_points)) in [-1, 0, 1]:
            raise Exception("This should not happen")

        # Make sure that
        peaks, base_points = assert_structure(xcoord, yvalues, peaks, base_points)

    return peaks, base_points


def build_fitted_spectrum(waveobs, lines, continuum_model=None):
    if continuum_model == None:
        # Build a continuum using the baseline of the fitted lines
        # - It may contain gaps in regions with few lines
        continuum_spectra = np.recarray((len(lines['mu']), ), dtype=[('waveobs', float),('flux', float),('err', float)])
        continuum_spectra['waveobs'] = lines['mu']
        continuum_spectra['flux'] = lines['baseline']
        continuum_model = MultiLinearInterpolationContinuumModel(continuum_spectra)

    # Build a fitted spectrum
    fluxes = np.zeros(len(waveobs))

    num_lines = len(lines)
    i = 0
    twopi = 2*np.pi
    for line in lines:
        if (i % 100) == 0:
            print "%.2f%%" % (((i*1.0)/num_lines) * 100)

        if line['gamma'] == 9999.0:
            mu = line['mu']
            A = line['A']
            sig = line['sig']
            squared_sig = sig**2
            # Build the gaussian corresponding to the line
            line_flux = (A / np.sqrt(twopi*squared_sig)) * np.exp(- (np.power((waveobs - mu), 2)) / (2*squared_sig))
        else:
            mu = line['mu']
            A = line['A']
            sig = line['sig']
            gamma = line['gamma']
            # Build the Voigt corresponding to the line
            if sig == 0:
                # Equivalent to a Lorentzian model
                line_flux = A*gamma/np.pi/(waveobs*waveobs - 2*waveobs*mu+mu*mu+gamma*gamma)
            else:
                # Voigt model (Gaussian and Lorentzian)
                from scipy.special import wofz
                w = wofz(((waveobs - mu) + 1j*gamma)* 2**-0.5/sig)
                line_flux = A * w.real*(twopi)**-0.5/sig
        wfilter = (line_flux < fluxes) & (line_flux < 0)
        #if np.min(line_flux) < -1:
        #    print i, np.max(line_flux), np.min(line_flux), np.std(line_flux), model_rms[i]
        fluxes[wfilter] = line_flux[wfilter]
        i += 1

    fluxes += continuum_model(waveobs)
    total_wavelengths = len(fluxes)
    fitted_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    fitted_spectra['waveobs'] = waveobs
    fitted_spectra['flux'] = fluxes
    return fitted_spectra

# Given a spectra, the position of a peak and its limiting region where:
# - base
# - Typical shape: concave + convex + concave region
# - Peak is located within the convex region, although it is just in the border with
#   the concave region (first element)
def improve_linemask_edges(xcoord, yvalues, base, top, peak):
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
        # This will happen very rarely (only in peaks detected at the extreme of a spectra
        # and one of its basepoints has been "artificially" added and it happens to be
        # just next to the peak)
        new_base = base
        new_top = original_top
        #plt.plot(x, y)
        #l = plt.axvline(x = x[peak_relative_pos], linewidth=1, color='red')
        #print d2y_dx2, d2y_dx2[peak_relative_pos]
        #plt.show()

    return new_base, new_top


def print_linemasks_stats(linemasks, discarded):
    # Profile types
    gaussian = linemasks['gamma'] == 9999.0
    voigt = linemasks['gamma'] != 9999.0

    print "---------------GENERAL----------------"
    print "\t\tMean\tMedian\tStdev"
    print "RMS:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['rms'][~discarded]), np.median(linemasks['rms'][~discarded]), np.std(linemasks['rms'][~discarded]))
    print "A:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['A'][~discarded]), np.median(linemasks['A'][~discarded]), np.std(linemasks['A'][~discarded]))
    print "sig:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['sig'][~discarded]), np.median(linemasks['sig'][~discarded]), np.std(linemasks['sig'][~discarded]))
    if len(linemasks[voigt & ~discarded]) > 0:
        print "gamma:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['gamma'][voigt & ~discarded]), np.median(linemasks['gamma'][voigt & ~discarded]), np.std(linemasks['gamma'][voigt & ~discarded]))
    else:
        print "gamma:\t\t-\t-\t-"
    print "R:\t\t%i\t%i\t%i" % (np.median(linemasks['R'][~discarded]), np.mean(linemasks['R'][~discarded]), np.std(linemasks['R'][~discarded]))
    print "FWHM (nm):\t%.3f\t%.3f\t%.3f" % (np.median(linemasks['fwhm'][~discarded]), np.mean(linemasks['fwhm'][~discarded]), np.std(linemasks['fwhm'][~discarded]))
    print "FWHM (km/s):\t%.3f\t%.3f\t%.3f" % (np.median(linemasks['fwhm_kms'][~discarded]), np.mean(linemasks['fwhm_kms'][~discarded]), np.std(linemasks['fwhm_kms'][~discarded]))
    if len(linemasks[gaussian & ~discarded]) > 0:
        print "---------------GAUSSIAN----------------"
        print "\t\tMean\tMedian\tStdev"
        print "RMS:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['rms'][gaussian & ~discarded]), np.median(linemasks['rms'][gaussian & ~discarded]), np.std(linemasks['rms'][gaussian & ~discarded]))
        print "A:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['A'][gaussian & ~discarded]), np.median(linemasks['A'][gaussian & ~discarded]), np.std(linemasks['A'][gaussian & ~discarded]))
        print "sig:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['sig'][gaussian & ~discarded]), np.median(linemasks['sig'][gaussian & ~discarded]), np.std(linemasks['sig'][gaussian & ~discarded]))
        print "R:\t\t%i\t%i\t%i" % (np.median(linemasks['R'][gaussian & ~discarded]), np.mean(linemasks['R'][gaussian & ~discarded]), np.std(linemasks['R'][gaussian & ~discarded]))
        print "FWHM (nm):\t%.3f\t%.3f\t%.3f" % (np.median(linemasks['fwhm'][gaussian & ~discarded]), np.mean(linemasks['fwhm'][gaussian & ~discarded]), np.std(linemasks['fwhm'][gaussian & ~discarded]))
        print "FWHM (km/s):\t%.3f\t%.3f\t%.3f" % (np.median(linemasks['fwhm_kms'][gaussian & ~discarded]), np.mean(linemasks['fwhm_kms'][gaussian & ~discarded]), np.std(linemasks['fwhm_kms'][gaussian & ~discarded]))
    if len(linemasks[voigt & ~discarded]) > 0:
        print "---------------VOIGT-------------------"
        print "\t\tMean\tMedian\tStdev"
        print "RMS:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['rms'][voigt & ~discarded]), np.median(linemasks['rms'][voigt & ~discarded]), np.std(linemasks['rms'][voigt & ~discarded]))
        print "A:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['A'][voigt & ~discarded]), np.median(linemasks['A'][voigt & ~discarded]), np.std(linemasks['A'][voigt & ~discarded]))
        print "sig:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['sig'][voigt & ~discarded]), np.median(linemasks['sig'][voigt & ~discarded]), np.std(linemasks['sig'][voigt & ~discarded]))
        print "gamma:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['gamma'][voigt & ~discarded]), np.median(linemasks['gamma'][voigt & ~discarded]), np.std(linemasks['gamma'][voigt & ~discarded]))
        print "R:\t\t%i\t%i\t%i" % (np.median(linemasks['R'][voigt & ~discarded]), np.mean(linemasks['R'][voigt & ~discarded]), np.std(linemasks['R'][voigt & ~discarded]))
        print "FWHM (nm):\t%.3f\t%.3f\t%.3f" % (np.median(linemasks['fwhm'][voigt & ~discarded]), np.mean(linemasks['fwhm'][voigt & ~discarded]), np.std(linemasks['fwhm'][voigt & ~discarded]))
        print "FWHM (km/s):\t%.3f\t%.3f\t%.3f" % (np.median(linemasks['fwhm_kms'][voigt & ~discarded]), np.mean(linemasks['fwhm_kms'][voigt & ~discarded]), np.std(linemasks['fwhm_kms'][voigt & ~discarded]))


if __name__ == '__main__':
    #### VALD line list
    ## Original VALD linelist
    #vald_file = "input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst"
    ##vald_file = "input/linelists/uves_linelist_mpa_v3.sme"
    #output_file = "input/linelists/VALD.300_1100nm.lst"
    #minimum_depth = 0.0
    #vald_linelist = read_VALD_linelist(vald_file, minimum_depth=minimum_depth)
    #linelist = VALD_to_SPECTRUM_format(vald_linelist)
    ## Filter discarded:
    #linelist = linelist[linelist['species'] != "Discard"]
    #asciitable.write(linelist, output=output_file, Writer=asciitable.FixedWidthNoHeader, delimiter=None, bookend=False, formats={'wave (A)': '%4.3f', })
    #import ipdb
    #ipdb.set_trace()

    #VALD_top_3_to_RV_format("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/rv/VALD.300_1100nm.rv.lst", top=1, wave_step=10)

    ### Line detection and fitting example:

    #############################
    ## Stars' spectra
    #############################
    #star, resolution = "input/spectra/examples/espadons_mu_leo.s.gz", 80000
    #star, resolution = "input/spectra/examples/espadons_mu_leo_norm.s.gz", 80000
    #star, resolution = "input/spectra/examples/harps_procyon.s.gz", 115000
    #star, resolution = "input/spectra/examples/harps_procyon_norm.s.gz", 115000
    #star, resolution = "input/spectra/examples/narval_arcturus.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_arcturus_norm.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_mu_cas.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_mu_cas_norm.s.gz", 80000
    #star, resolution = "input/spectra/examples/narval_sun.s.gz", 80000
    star, resolution = "input/spectra/examples/narval_sun_norm.s.gz", 80000
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

    # Normalize spectra
    normalize = False #True
    # Smooth spectra using the instrumental resolution
    smooth_spectra = True

    #############################
    ## Synthetic telluric lines' spectra
    #star, resolution = "input/spectra/synthetic/telluric_standard_atm_air_model.s.gz", 200000
    #star, resolution = "input/spectra/synthetic/telluric_standard_atm_air_model_norm.s.gz", 200000
    # Do NOT smooth spectra using the instrumental resolution
    #smooth_spectra = True
    #############################

    ### Fitting parameters
    discard_gaussian = False
    discard_voigt = True
    generate_fitted_spectra = True

    ### Filtering parameters
    # Discard lines that do not have at least a given depth
    minimum_depth = 0.05 # (% of the continuum)
    maximum_depth = 1.00 # (% of the continuum)
    # Discard potential gaps in the spectra that are identified as lines
    discard_too_big_wavelength_range = False
    # Discard outliers (biggest RMS)
    discard_outlier_fit = False
    # Discard outliers in resolving power
    discard_outlier_R = False

    ## Reading spectra
    print "Reading spectrum..."
    spectra = read_spectra(star)

    if normalize:
        print "Continuum normalization..."
        continuum_model_for_normalization = fit_continuum(spectra)
        spectra['flux'] /= continuum_model_for_normalization(spectra['waveobs'])

    #--- Radial velocity ----------------------------------------------
    print "Radial Velocity...",
    #--- Line lists ----------------------------------------------
    # - Atomic
    vald_linelist_file = "input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst"
    linelist_atomic = read_VALD_linelist(vald_linelist_file, minimum_depth=0.0)
    xcoord, fluxes, num_used_lines = build_velocity_profile(spectra, \
            linelist_atomic, lower_velocity_limit=-200.0, upper_velocity_limit=200.0, velocity_step=0.5)
    models = modelize_velocity_profile(xcoord, fluxes)
    rv = np.round(models[0].mu(), 2) # km/s
    spectra = correct_velocity(spectra, rv)
    print rv, "km/s"

    original_spectra = spectra
    if smooth_spectra:
        print "Convolving (smoothing)..."
        # Use 2 times the resolution to smooth the spectra without loosing too much details
        spectra = convolve_spectra(original_spectra, 2*resolution)

    print "Finding peaks and base points..."
    peaks, base_points = find_peaks_and_base_points(spectra['waveobs'], spectra['flux'])

    # Determine continuum
    print "Determining continuum..."
    continuum_model = fit_continuum(spectra)

    print "Generating linemasks, fitting gaussians/voigt and matching VALD lines..."
    linemasks = generate_linemasks(original_spectra, peaks, base_points, continuum_model, minimum_depth=minimum_depth, maximum_depth=maximum_depth, smoothed_spectra=spectra ,vald_linelist_file="input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", telluric_linelist_file = "input/linelists/telluric/standard_atm_air_model.lst", discard_gaussian = discard_gaussian, discard_voigt = discard_voigt)

    ####################################################################
    ##### FILTERS
    ####################################################################
    # - Mandatory filters
    ####################################################################
    print "Applying peak filters..."
    rejected_by_noise = detect_false_positives_and_noise(spectra, linemasks)

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

    rejected_by_atomic_line_not_found = (linemasks['VALD_wave_peak'] == 0)
    rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)

    discarded = rejected_by_noise
    discarded = np.logical_or(discarded, rejected_by_depth_higher_than_continuum)
    discarded = np.logical_or(discarded, rejected_by_depth_limits)
    discarded = np.logical_or(discarded, rejected_by_bad_fit)
    discarded = np.logical_or(discarded, rejected_by_atomic_line_not_found)
    discarded = np.logical_or(discarded, rejected_by_telluric_line)
    mandatory_discarded = discarded.copy()
    ####################################################################
    # - Optional filters for outliers
    #   * Outliers identification should be done avoiding the already discarded points
    #     otherwise, we will have noise contamination
    ####################################################################
    # Identify linemasks with too big wavelength range (outliers)
    # WARNING: If the spectrum has a small wavelength range with some few strong features
    #          they may be erroneously discarded by this filter
    wave_diff = (linemasks['wave_top'] - linemasks['wave_base']) / (linemasks['top'] - linemasks['base'])
    wave_diff_selected, wave_diff_selected_filter = sigma_clipping(wave_diff[~mandatory_discarded], sig=3, meanfunc=np.median) # Discard outliers
    accepted_index = np.arange(len(wave_diff))[~mandatory_discarded][wave_diff_selected_filter]
    rejected_by_wave_gaps = wave_diff < 0 # Create an array of booleans
    rejected_by_wave_gaps[:] = True # Initialize
    rejected_by_wave_gaps[accepted_index] = False
    # If the peak was already discarded previously by the mandatory filters, do not mark it
    rejected_by_wave_gaps[mandatory_discarded] = False

    if discard_too_big_wavelength_range:
        discarded = np.logical_or(discarded, rejected_by_wave_gaps)

    # Identify outliers considering RMS
    rms_selected, rms_selected_filter = sigma_clipping(linemasks['rms'][~mandatory_discarded], sig=3, meanfunc=np.median) # Discard outliers
    accepted_index = np.arange(len(linemasks['rms']))[~mandatory_discarded][rms_selected_filter]
    rejected_by_outlier_fit = linemasks['rms'] < 0 # Create an array of booleans
    rejected_by_outlier_fit[:] = True # Initialize
    rejected_by_outlier_fit[accepted_index] = False
    # If the peak was already discarded previously by the mandatory filters, do not mark it
    rejected_by_outlier_fit[mandatory_discarded] = False

    if discard_outlier_fit:
        discarded = np.logical_or(discarded, rejected_by_outlier_fit)

    # Identify outliers considering the resolution
    r_selected, r_selected_filter = sigma_clipping(linemasks['R'][~mandatory_discarded], sig=3, meanfunc=np.median) # Discard outliers
    accepted_index = np.arange(len(linemasks['R']))[~mandatory_discarded][r_selected_filter]
    rejected_by_outlier_R = linemasks['R'] < 0 # Create an array of booleans
    rejected_by_outlier_R[:] = True # Initialize
    rejected_by_outlier_R[accepted_index] = False
    # If the peak was already discarded previously by the mandatory filters, do not mark it
    rejected_by_outlier_R[mandatory_discarded] = False

    if discard_outlier_R:
        discarded = np.logical_or(discarded, rejected_by_outlier_R)

    # Base points near the continuum in 5%
    c1 = continuum_model(linemasks['wave_base_fit'])
    f1 = spectra['flux'][linemasks['base_fit']]
    d1 = ((c1 - f1) / c1)
    rejected_by_not_ideal1 = d1 > 0.05
    c2 = continuum_model(linemasks['wave_top_fit'])
    f2 = spectra['flux'][linemasks['top_fit']]
    d2 = ((c2 - f2) / c2)
    rejected_by_not_ideal2 = d2 > 0.05
    # Base point with similar fluxes (less than 1% difference)
    x = np.abs(f1 - f2) / f1
    rejected_by_not_ideal = x > 0.05
    #rejected_by_not_ideal = np.logical_or(rejected_by_not_ideal1, rejected_by_not_ideal2)
    #rejected_by_not_ideal = np.logical_or(rejected_by_not_ideal, rejected_by_not_ideal3)
    discarded = np.logical_or(discarded, rejected_by_not_ideal)

    # Identify outliers in the cross-matching with the VALD line list
    #wave_peak_diff = np.abs(linemasks['wave_peak'] - linemasks['VALD_wave_peak'])
    ##wave_diff_selected, wave_diff_selected_filter = sigma_clipping(wave_peak_diff[~discarded], sig=3, meanfunc=np.median) # Discard outliers
    ##accepted_index = np.arange(len(wave_peak_diff))[~discarded][wave_diff_selected_filter]
    ##bad_VALD_crossmatch = wave_peak_diff < 0 # Create an array of booleans
    ##bad_VALD_crossmatch[:] = True # Initialize
    ##bad_VALD_crossmatch[accepted_index] = False
    #bad_VALD_crossmatch = wave_peak_diff > 0.005

    #telluric_R_selected, telluric_R_selected_filter = sigma_clipping(linemasks['telluric_R'][~discarded], sig=3, meanfunc=np.median) # Discard outliers
    #accepted_index = np.arange(len(linemasks['telluric_R']))[~discarded][telluric_R_selected_filter]
    #bad_telluric_R = linemasks['telluric_R'] < 0 # Create an array of booleans
    #bad_telluric_R[:] = True # Initialize
    #bad_telluric_R[accepted_index] = False


    #telluric_wave_peak_diff = np.abs(linemasks['wave_peak'] - linemasks['telluric_wave_peak'])
    ##wave_diff_selected, wave_diff_selected_filter = sigma_clipping(telluric_wave_peak_diff[~discarded], sig=3, meanfunc=np.median) # Discard outliers
    ##accepted_index = np.arange(len(wave_peak_diff))[~discarded][wave_diff_selected_filter]
    ##bad_telluric_crossmatch = telluric_wave_peak_diff < 0 # Create an array of booleans
    ##bad_telluric_crossmatch[:] = True # Initialize
    ##bad_telluric_crossmatch[accepted_index] = False
    #bad_telluric_crossmatch = (telluric_wave_peak_diff > 0.005) & (linemasks['telluric_fwhm'] != 0)

    ####################################################################
    linemasks['discarded'][discarded] = True
    ####################################################################
    ##### END FILTERS
    ####################################################################

    # Profile types
    gaussian = linemasks['gamma'] == 9999.0
    voigt = linemasks['gamma'] != 9999.0

    ## Peaks
    print "--------------------------------------"
    print "Number of preliminary peaks:\t", len(linemasks)
    print "- Noise and false positive:\t", len(linemasks[rejected_by_noise])
    print "- Peaks higher than continuum:\t", len(linemasks[rejected_by_depth_higher_than_continuum])
    print "- Bad fits:\t\t\t", len(linemasks[rejected_by_bad_fit])
    print "- Out of depth limits:\t\t", len(linemasks[rejected_by_depth_limits]), "\t[%.2f - %.2f]" % (minimum_depth, maximum_depth)
    print "- No atomic line data:\t\t", len(linemasks[rejected_by_atomic_line_not_found])
    print "- Telluric line:\t\t", len(linemasks[rejected_by_telluric_line])
    print "Number of line candidates:\t", len(linemasks[~mandatory_discarded])
    if discard_too_big_wavelength_range:
        print "- Too big wavelength range:\t", len(linemasks[rejected_by_wave_gaps]), "\t[Enabled]"
    else:
        print "- Too big wavelength range:\t", len(linemasks[rejected_by_wave_gaps]), "\t[Disabled]"
    if discard_outlier_fit:
        print "- Outlier fits:\t\t\t", len(linemasks[rejected_by_outlier_fit]), "\t[Enabled]"
    else:
        print "- Outlier fits:\t\t\t", len(linemasks[rejected_by_outlier_fit]), "\t[Disabled]"
    if discard_outlier_R:
        print "- Resolving power outliers:\t", len(linemasks[rejected_by_outlier_R]), "\t[Enabled]"
    else:
        print "- Resolving power outliers:\t", len(linemasks[rejected_by_outlier_R]), "\t[Disabled]"
    print "- Not ideal:\t\t\t", len(linemasks[rejected_by_not_ideal]), "\t[Enabled]"
    print "Final number of lines:\t\t", len(linemasks[~discarded])
    print "- Gaussian profile:\t\t", len(linemasks[gaussian & ~discarded])
    print "- Voigt profile:\t\t", len(linemasks[voigt & ~discarded])
    print_linemasks_stats(linemasks, discarded)

    if generate_fitted_spectra:
        print "Building a fitted spectrum..."
        #waveobs = generate_wavelength_grid(np.min(spectra['waveobs']), np.max(spectra['waveobs']), resolution, points_per_fwhm = 3)
        waveobs = spectra['waveobs']
        fitted_spectra = build_fitted_spectrum(waveobs, linemasks[~discarded], continuum_model)

        # Plot
        fig = plt.figure()
        plt.plot(spectra['waveobs'], original_spectra['flux'])
        plt.plot(spectra['waveobs'], continuum_model(spectra['waveobs']))
        plt.plot(fitted_spectra['waveobs'], fitted_spectra['flux'])
        #plt.scatter(spectra['waveobs'][peaks], original_spectra['flux'][peaks], c='red')
        #plt.scatter(spectra['waveobs'][base_points], original_spectra['flux'][base_points], c='green')
        plt.scatter(spectra['waveobs'][peaks[~discarded]], original_spectra['flux'][peaks[~discarded]], c='red')
        plt.scatter(spectra['waveobs'][base_points[~discarded]], original_spectra['flux'][base_points[~discarded]], c='green')
        plt.scatter(spectra['waveobs'][base_points[1:][~discarded[:-1]]], original_spectra['flux'][base_points[1:][~discarded[:-1]]], c='green')
        plt.show()
    else:
        fitted_spectra = None

    ##### Saving...
    dump_filename = "output/" + os.path.basename(star) + ".dump"
    # Complete dump
    version = 20120402
    data = (os.path.abspath(star), resolution, linemasks, fitted_spectra)
    pickle.dump((version, data), gzip.open(dump_filename, "wb", compresslevel=3), protocol=2)

    ##### Restoring...
    #dump_filename = "output/" + os.path.basename(star) + ".dump"
    #version, data = pickle.load(gzip.open(dump_filename, "rb"))
    #star, resolution, linemasks, fitted_spectra = data

    ## Build a continuum using the baseline of the fitted lines
    ## - It may contain gaps in regions with few lines
    #continuum_spectra = np.recarray((len(linemasks['mu'][~discarded]), ), dtype=[('waveobs', float),('flux', float),('err', float)])
    #continuum_spectra['waveobs'] = linemasks['mu'][~discarded]
    #continuum_spectra['flux'] = linemasks['baseline'][~discarded]
    #continuum_model = MultiLinearInterpolationContinuumModel(continuum_spectra)
    ##################

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['wave_peak'], linemasks['fwhm'][~discarded])
    #plt.scatter(linemasks[~discarded]['wave_peak'], linemasks['fwhm'][~discarded], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('FWHM (nm)')
    #plt.show()

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['wave_peak'], linemasks['fwhm_kms'][~discarded])
    #plt.scatter(linemasks[~discarded]['wave_peak'], linemasks['fwhm_kms'][~discarded], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('FWHM (km/s)')
    #plt.show()

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['wave_peak'], linemasks['R'][~discarded])
    #plt.scatter(linemasks[~discarded]['wave_peak'], linemasks['R'][~discarded], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~rejected_by_outlier_R]['wave_peak'], linemasks['R'][~rejected_by_outlier_R])
    #plt.scatter(linemasks[~rejected_by_outlier_R]['wave_peak'], linemasks['R'][~rejected_by_outlier_R], s=4)
    #plt.plot(linemasks[~rejected_by_outlier_R]['wave_peak'], trend(linemasks[~rejected_by_outlier_R]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #r = (linemasks['telluric_R'] > 0) #& (linemasks['telluric_R'] < 100000)
    ##r = r & (np.abs(linemasks['depth'] - linemasks['telluric_depth']) < 0.1)
    ##r = r & (linemasks['fwhm'] - linemasks['telluric_fwhm'] < 0.04)
    #r = r & (linemasks['VALD_wave_peak'] == 0)
    #trend.fitData(linemasks[r]['wave_peak'], linemasks['telluric_R'][r])
    #plt.plot(linemasks[r]['wave_peak'], linemasks['telluric_R'][r])
    ##plt.plot(linemasks[r]['wave_peak'], trend(linemasks[r]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()

    #print np.mean(linemasks['telluric_R'][r])

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[gaussian & ~discarded]['wave_peak'], linemasks['R'][gaussian & ~discarded])
    #plt.scatter(linemasks[gaussian & ~discarded]['wave_peak'], linemasks['R'][gaussian & ~discarded], s=4)
    #plt.plot(linemasks[gaussian & ~discarded]['wave_peak'], trend(linemasks[gaussian & ~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()

    #show_histogram(linemasks['R'][~discarded][linemasks['R'][~discarded] < 100000])

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[gaussian & ~discarded]['depth'], linemasks['R'][gaussian & ~discarded])
    #plt.scatter(linemasks[gaussian & ~discarded]['depth'], linemasks['R'][gaussian & ~discarded], s=4)
    #plt.plot(linemasks[gaussian & ~discarded]['depth'], trend(linemasks[gaussian & ~discarded]['depth']), color="red", linewidth=3)
    #plt.xlabel('Depth')
    #plt.ylabel('Resolving power')
    #plt.show()

    #from pymodelfit import LinearModel
    #plt.scatter(linemasks['R'][gaussian & ~discarded], linemasks[gaussian & ~discarded]['depth'], s=4)
    #plt.ylabel('Depth')
    #plt.xlabel('Resolving power')
    #plt.show()

    ############ STATS
    #wave_peak_diff = np.abs(linemasks['wave_peak'] - linemasks['VALD_wave_peak'])
    #show_histogram(np.abs(wave_peak_diff), nbins=100)
    #show_histogram(np.abs(linemasks["depth"] - linemasks["relative_depth"]), nbins=100)
    #show_histogram(linemasks["rms"], nbins=100)
    #show_histogram(linemasks["A"], nbins=100)
    #show_histogram(linemasks["sig"], nbins=100)
    #
    #show_histogram(linemasks["depth"], nbins=100)
    #show_histogram(linemasks["relative_depth"], nbins=100)
    #show_histogram(np.abs(linemasks["depth"] - linemasks["solar_depth"]), nbins=100)
    #show_histogram(linemasks[~discarded]['rms'], nbins=100)

    #show_histogram(linemasks[~discarded]['sig'], xlabel='sig', nbins=100)
    #show_histogram(linemasks[~discarded]['A'], xlabel='A', nbins=100)
    #show_histogram(linemasks[~discarded]['gamma'], xlabel='gamma', nbins=100)

    #show_histogram(linemasks['R'][~discarded][linemasks['R'][~discarded] < 200000], xlabel='R', nbins=100)


    #elements = np.unique(linemasks['element'][~discarded])
    #x = np.arange(len(elements))
    #y = np.zeros(len(elements))
    #for element, i in zip(elements, x):
        #subset = linemasks[linemasks['element'][~discarded] == element]
        #y[i] = np.median(subset['depth'])
    #plt.bar(x, y, align='center')
    #plt.xlabel("Elements")
    #plt.ylabel("Median depth")
    #t = plt.xticks(x, elements, size='small')
    #plt.show()

    #elements = np.unique(linemasks['element'][~discarded])
    #x = np.arange(len(elements))
    #y = np.zeros(len(elements))
    #for element, i in zip(elements, x):
        #subset = linemasks[linemasks['element'][~discarded] == element]
        #y[i] = len(subset['depth'])
    #plt.bar(x, y, align='center')
    #plt.xlabel("Elements")
    #plt.ylabel("Count")
    #t = plt.xticks(x, elements, size='small')
    #plt.show()


    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['mu'], linemasks[~discarded]['sig'])
    #plt.plot(linemasks[~discarded]['mu'], linemasks[~discarded]['sig'])
    #plt.plot(linemasks[~discarded]['mu'], trend(linemasks[~discarded]['mu']), color="red", linewidth=4)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Sigma')
    #plt.show()

    #plt.plot(linemasks[~discarded]['mu'], linemasks[~discarded]['rms'])
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('RMS')
    #plt.show()

    #wave_peak_diff = np.abs(linemasks['wave_peak'] - linemasks['VALD_wave_peak'])
    #plt.plot(linemasks[~discarded]['mu'], np.abs(wave_peak_diff[~discarded]))
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Wavelength diff. with VALD')
    #plt.show()

    #plt.scatter(linemasks[~discarded]['depth'], linemasks[~discarded]['log(gf)'], s=4)
    #plt.xlabel('Depth')
    #plt.ylabel('log(gf)')
    #plt.show()

    #plt.scatter(linemasks[~discarded]['depth'], linemasks[~discarded]['lower state (eV)'], s=4)
    #plt.xlabel('Depth')
    #plt.ylabel('lower state (eV)')
    #plt.show()

    #plt.scatter(linemasks[~discarded]['log(gf)'], linemasks[~discarded]['lower state (eV)'], s=4)
    #plt.xlabel('log(gf)')
    #plt.ylabel('lower state (eV)')
    #plt.show()

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot_wireframe(, rstride=10, cstride=10)
    #plt.show()

    ##import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization!
    #fig = plt.figure()
    #ax = Axes3D(fig) #<-- Note the difference from your original code...
    #X, Y, Z = linemasks[~discarded]['log(gf)'], linemasks[~discarded]['lower state (eV)'], linemasks[~discarded]['depth']
    #cset = ax.contour(X, Y, Z, 16, extend3d=True)
    #ax.clabel(cset, fontsize=9, inline=1)
    #plt.show()

    #plt.scatter(linemasks[~discarded]['depth'], linemasks[~discarded]['solar_depth'], s=4)
    #plt.xlabel('Depth')
    #plt.ylabel('Solar depth')
    #plt.show()

    #plt.scatter(linemasks[~discarded]['depth'], linemasks[~discarded]['relative_depth'], s=4)
    #plt.xlabel('Depth')
    #plt.ylabel('Relative depth')
    #plt.show()
