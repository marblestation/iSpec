import asciitable
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions
import cPickle as pickle
import gzip
import os
import ipdb
from common import *
from continuum import *
from radial_velocity import *
from convolve import *
import matplotlib.pyplot as plt
from pymodelfit import UniformKnotSplineModel
from pymodelfit import GaussianModel
from pymodelfit import VoigtModel

########################################################################
## [START] LINE LISTS
########################################################################
# Load a VALD linelist (short format) and filter it
# - depth_limit: filter out all the lines with a depth less than this percentage (0.05 = 5%)
# - data_end can be negative in order to ignore the last nth rows (VALD usually 
#            adds references to the end of the file that should be ignored)
def load_and_filter_VALD(vald_file, depth_limit=0.0, data_end=None):
    # Original VALD linelist
    if data_end == None:
        vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "Rad", "Stark", "Waals", "factor", "depth", "Reference"], exclude_names=["Vmic (km/s)", "Rad", "Stark", "Waals", "factor", "Reference"], guess=False)
    else:
        vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, data_end=data_end, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "Rad", "Stark", "Waals", "factor", "depth", "Reference"], exclude_names=["Vmic (km/s)", "Rad", "Stark", "Waals", "factor", "Reference"], guess=False)

    if depth_limit <= 0.0:
        return vald
    else:
        # Filter
        vald_limited = vald[vald['depth'] >= depth_limit]
        return vald_limited


# Convert a VALD linelist (short format) to a format that can be used with SPECTRUM
# - depth_limit: filter out all the lines with a depth less than this percentage (0.05 = 5%)
# - data_end can be negative in order to ignore the last nth rows (VALD usually 
#            adds references to the end of the file that should be ignored)
def VALD_to_SPECTRUM_format(vald_file, output_file, depth_limit=0.0, data_end=None):
    # Planck constant
    h = 6.626068 * 10e-34 # m^2 kg / s
    # Light speed in vacuum
    c = 299792458.0 # m/s

    # Original VALD linelist
    vald_limited = load_and_filter_VALD(vald_file, depth_limit=depth_limit, data_end=data_end)
    
    # Periodic table
    table = asciitable.read("input/abundances/chemical_elements_symbols.dat", delimiter="\t")

    # Prepare resulting structure
    linelist = np.recarray((len(vald_limited), ), dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'),('transition type', '|S10'), ('note', '|S100')])

    i = 0
    for line in vald_limited:
        linelist[i]['wave (A)'] = line['wave (A)']
        
        element = line['element'].split(" ")
        symbol = element[0]
        try:
            element.remove('') # Make sure there are not additional spaces between the symbol and the ionization state
            element.remove('')
            element.remove('')
        except ValueError as e:
            pass
        ionization = str(int(element[1]) - 1)
        
        tfilter = (table['symbol'] == symbol)
        linelist[i]['species'] = str(table[tfilter]["atomic_num"][0]) + "." + ionization
        
        #print linelist[i]['species']
        linelist[i]['lower state (cm^-1)'] = int(line['lower state (eV)'] * 8065.73) #cm-1
        # Wavelength
        l = (line[ "wave (A)"] / 10) * 10e-9 # m
        # Frequency
        f = c/l # Hz
        # Energy
        E = h * f # Joules
        E = E * 6.24150974e18 # electron Volt (eV)
        E = E * 8065.73 #cm-1
        linelist[i]['upper state (cm^-1)'] = int(linelist[i]['lower state (cm^-1)'] + E)
        linelist[i]['log(gf)'] = line['log(gf)'] *1.5
        linelist[i]['fudge factor'] = 1.0
        linelist[i]['transition type'] = "99"
        linelist[i]['note'] = line["element"].replace(" ", "_")
        i += 1

    asciitable.write(linelist, output=output_file, Writer=asciitable.FixedWidthNoHeader, delimiter=None, bookend=False, formats={'wave (A)': '%4.3f', })


# Convert a VALD linelist (short format) to a format that can be used to measure radial velocities
# and select the top N deepest lines every wave_step armstrongs (1 nm).
# - data_end can be negative in order to ignore the last nth rows (VALD usually 
#            adds references to the end of the file that should be ignored)
def VALD_top_N_to_RV_format(vald_file, output_file, top=1, wave_step=10, data_end=None):

    vald_limited = load_and_filter_VALD(vald_file, depth_limit=0, data_end=data_end)
    
    wave_base = np.min(vald_limited["wave (A)"])
    wave_top = np.max(vald_limited["wave (A)"])
    
    # Prepare resulting structure
    linelist = np.recarray((top*np.ceil((wave_top - wave_base) / wave_step), ), dtype=[('wave_peak', '<f8'), ('depth', '<f8'), ('element', '|S100')])
    
    wave_current = wave_base
    i = 0
    # For each segment
    while wave_current < wave_top:
        wfilter = (vald_limited["wave (A)"] >= wave_current) & (vald_limited["wave (A)"] < wave_current + wave_step)
        vald_filtered = vald_limited[wfilter]
        vald_filtered.sort(order="depth")
        # Select the top 3 deepest lines
        for j in np.arange(top):
            pos = -1*(j+1)
            linelist[i+j]["wave_peak"] = vald_filtered[pos]["wave (A)"] / 10 #nm
            linelist[i+j]["depth"] = vald_filtered[pos]["depth"]
            linelist[i+j]["element"] = vald_filtered[pos]["element"]
        wave_current += wave_step
        i += top
    
    linelist.sort(order="wave_peak")
    asciitable.write(linelist, output=output_file, delimiter="\t")

########################################################################
## [END] LINE LIST
########################################################################

########################################################################
## [START] SIGNAL-TO-NOISE = EQUIVALENT WIDTH / ERROR
##   Schneider et al (1993). The Hubble Space Telescope quasar 
##   absorption line key project. II - Data calibration and 
##   absorption-line selection
########################################################################
# Calculate the Equivalent Width (EW) and EW error per pixel bin (measure)
# and convolve by using the instruction resolution.
# - Flux and error should be already normalized (divided by the continuum)
def get_convolved_ew(spectra, resolution):
    # Consider the wavelength of the measurements as the center of the bins
    waveobs = spectra['waveobs']
    # Calculate the wavelength distance between the center of each bin
    wave_distance = waveobs[1:] - waveobs[:-1]
    # Define the edge of each bin as half the wavelength distance to the bin next to it
    edges_tmp = waveobs[:-1] + 0.5 * (wave_distance)
    # Define the edges for the first and last measure which where out of the previous calculations
    first_edge = waveobs[0] - 0.5*wave_distance[0]
    last_edge = waveobs[-1] + 0.5*wave_distance[-1]
    # Build the final edges array
    edges = np.array([first_edge] + edges_tmp.tolist() + [last_edge])
    
    # Bin width
    bin_width = edges[1:] - edges[:-1]          # width per pixel
    
    # Equivalent width (EW) and EW error per bin
    ew = bin_width * (1-spectra['flux']) # EW positive if it is an absorption line
    ewer = bin_width * spectra['err'] 
    
    # FWHM of the gaussian for the given resolution
    fwhm = waveobs / resolution
    sigma = get_sigma(fwhm)
    # Convert from wavelength units to bins
    fwhm_bin = fwhm / bin_width
    
    # Round number of bins per FWHM
    nbins = np.ceil(fwhm_bin) #npixels
    
    # Number of measures
    nwaveobs = len(waveobs)
    
    convolved_ew = []
    convolved_ewer = []
    for i in np.arange(len(nbins)):
        current_nbins = 2 * nbins[i] # Each side
        current_center = waveobs[i] # Center
        current_sigma = sigma[i]
        
        # Find lower and uper index for the gaussian, taking care of the current spectra limits
        lower_pos = int(max(0, i - current_nbins))
        upper_pos = int(min(nwaveobs, i + current_nbins + 1))
        
        # Select only the EW values for the segment that we are going to convolve
        ew_segment = ew[lower_pos:upper_pos+1]
        ewer_segment = ewer[lower_pos:upper_pos+1]
        waveobs_segment = waveobs[lower_pos:upper_pos+1]
        
        nsegments = len(ew_segment)
        
        # Build the gaussian corresponding to the instrumental spread function
        gaussian = np.exp(- ((waveobs_segment - current_center)**2) / (2*current_sigma**2)) / np.sqrt(2*np.pi*current_sigma**2)
        gaussian = gaussian / np.sum(gaussian)
        
        # Convolve the current position by using the segment and the gaussian
        weighted_ew = ew_segment * gaussian
        weighted_ew_err = ewer_segment * gaussian
        current_convolved_ew = weighted_ew.sum() / nsegments # Normalized
        current_convolved_ewer = np.sqrt((weighted_ew_err**2).sum()) / nsegments # Normalized
        
        convolved_ew.append(current_convolved_ew)
        convolved_ewer.append(current_convolved_ewer)
    
    return np.array(convolved_ew), np.array(convolved_ewer)

# Determine the Equivalent Width (EW) and EW error per pixel bin (measure)
# and calculate the Singnal-to-Noise Ratio (SNR) per pixel bin.
# Finally, adjust a uniformly divided spline to smooth the SNR and correct
# small differences derived from original's spectrum inhomogeneites continuum normalization
def get_corrected_convolved_ew_snr(spectra, resolution):
    # Find EW (bigger EW == deeper absorption line) and convolved EW (smoother) using the instrument resolution
    convolved_ew, convolved_ewer = get_convolved_ew(spectra, resolution)
    valid = ~np.isnan(convolved_ew) & ~np.isnan(convolved_ewer) & (convolved_ewer > 0)
    convolved_ew = convolved_ew[valid]
    convolved_ewer  = convolved_ewer[valid]

    # Determine the Signal-to-Noise Ratio (SNR)
    convolved_ew_snr = convolved_ew / convolved_ewer

    # Find candidates to continuum points by looking for peaks in an inversed SNR
    continuum_candidates = find_min_win(convolved_ew_snr, span=3)

    # Discard continuum candidates below the mean
    limit = np.mean(convolved_ew_snr) # 80
    significant = (convolved_ew_snr[continuum_candidates] < limit)
    continuum_points = continuum_candidates[significant]

    ## TODO: Use fit_continuum() which now has a better implementation for this purposes
    # Fit a uniformly divided splines (n knots) in order to increase smoothness and
    # homogeneity (reducing differences between far away wavelength segments due to
    # small differences in the continuum normalization process) that will improve
    # the determination of peaks
    # * 1 knot every 10 nm
    nknots = np.max([1, int((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 10)])
    continuum_model = UniformKnotSplineModel(nknots=nknots)
    continuum_model.fitData(spectra['waveobs'][continuum_points], convolved_ew_snr[continuum_points])

    # Correct SNR
    corrected_convolved_ew_snr = convolved_ew_snr - continuum_model(spectra['waveobs'])

    return corrected_convolved_ew_snr
########################################################################
## [END] SIGNAL-TO-NOISE = EQUIVALENT WIDTH / ERROR
########################################################################

# Fits a gaussian at a given wavelength location using a fitted continuum model
# - For absorption lines, it will alway be true:
#      model.A < 0 and model.sig > 0
#   The oposite indicates a potential bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_gaussian(spectra_slice, continuum_model, mu, sig=0.02, A=-0.025, prioritize_deeper_fluxes=False):
    model = GaussianModel()
    model.mu = mu
    model.sig = sig
    model.A = A
    cont = continuum_model(spectra_slice['waveobs'])
    conterr = 0
    
    if prioritize_deeper_fluxes:
        # More weight to the deeper fluxes
        min_flux = np.min(spectra_slice['flux'])
        if min_flux < 0:
            weights = spectra_slice['flux'] + -1*(min_flux) + 0.01 # Above zero
            weights = np.min(weights) / weights
        else:
            weights = min_flux / spectra_slice['flux']
            
        # Priorizing deeper fluxes:
        model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=weights, fixedpars=[])
    else:
        # Without weigths:
        model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, fixedpars=[])
    
    # A gaussian with a positive A and a negative Sigma is exactly the same as the inverse
    # - We inverse the signs since it is more intuitive for absorption lines (negative amplitudes)
    #   and it will be easier to apply filtering rules after (such as, discard all gaussian/voigt with positive A)
    if model.A > 0 and model.sig < 0:
        model.A = -1 * model.A
        model.sig = -1 * model.sig
    
    return model

# Fits a voigt at a given wavelength location using a fitted continuum model
# - If prioritize_deeper_fluxes is True, values with lower flux values will
#   have more weight in the fitting process
# - For absorption lines, it will alway be true:
#      model.A < 0 and model.sig > 0
#   The oposite indicates a potential bad fit
# - For absorption lines, model.gamma < 0 indicates strange wings and probably a bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_voigt(spectra_slice, continuum_model, mu, sig=0.02, A=-0.025, gamma=0.025, prioritize_deeper_fluxes=False):
    model = VoigtModel()
    model.gamma = gamma
    model.mu = mu
    model.sig = sig
    model.A = A
    cont = continuum_model(spectra_slice['waveobs'])
    conterr = 0
    
    if prioritize_deeper_fluxes:
        # More weight to the deeper fluxes
        min_flux = np.min(spectra_slice['flux'])
        if min_flux < 0:
            weights = spectra_slice['flux'] + -1*(min_flux) + 0.01 # Above zero
            weights = np.min(weights) / weights
        else:
            weights = min_flux / spectra_slice['flux']
        # Priorizing deeper fluxes:
        model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=weights, fixedpars=[])
    else:
        # Without weigths:
        model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, fixedpars=[])
    
    # A voigt with a positive A, a negative Sigma and a negative gamma is exactly the same as the inverse
    # - We inverse the signs since it is more intuitive for absorption lines (negative amplitudes)
    #   and it will be easier to apply filtering rules after (such as, discard all gaussian/voigt with positive A)
    if model.A > 0 and model.sig < 0 and model.gamma < 0:
        model.A = -1 * model.A
        model.sig = -1 * model.sig
        model.gamma = -1 * model.gamma
    
    return model
    

# Fits a gaussian and a voigt at a given wavelength location using a fitted continuum model
# - It selects the best fitted model (gaussian or voigt) unless one of them is disabled by
#   the discard_gaussian or discard_voigt argument
# - If prioritize_deeper_fluxes is True, values with lower flux values will
#   have more weight in the fitting process
# - For absorption lines, it will alway be true:
#      model.A < 0 and model.sig > 0
#   The oposite indicates a potential bad fit
# - For absorption lines fitted with voigt, model.gamma < 0 indicates strange wings and probably a bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_line(spectra_slice, continuum_model, mu, sig=0.02, A=-0.025, gamma=0.025, discard_gaussian = False, discard_voigt = False, prioritize_deeper_fluxes=False):
    if not discard_gaussian:
        try:
            gaussian_model = fit_gaussian(spectra_slice, continuum_model, mu, sig=sig, A=A, prioritize_deeper_fluxes=prioritize_deeper_fluxes)
            rms_gaussian = np.sqrt(np.sum(np.power(gaussian_model.residuals(), 2)) / len(gaussian_model.residuals()))
            mu_residual_gaussian = np.min([9999.0, calculate_mu_residual(spectra_slice, gaussian_model)]) # Use 9999.0 as maximum to avoid overflows
        except Exception as e:
            rms_gaussian = 9999.0
            mu_residual_gaussian = 9999.0
            discard_gaussian = True
    
    if not discard_voigt:
        try:
            voigt_model = fit_voigt(spectra_slice, continuum_model, mu, sig=sig, A=A, gamma=gamma, prioritize_deeper_fluxes=prioritize_deeper_fluxes)
            rms_voigt = np.sqrt(np.sum(np.power(voigt_model.residuals(), 2)) / len(voigt_model.residuals()))
            mu_residual_voigt = np.min([9999.0, calculate_mu_residual(spectra_slice, voigt_model)]) # Use 9999.0 as maximum to avoid overflows
        except Exception as e:
            rms_voigt = 9999.0
            mu_residual_voigt = 9999.0
            discard_voigt = True
        
    # Mu (peak) inside the spectra region & Coherent parameters for an absorption line
    if not discard_gaussian and (gaussian_model.mu < spectra_slice['waveobs'][0] or gaussian_model.mu > spectra_slice['waveobs'][-1] or gaussian_model.A >= 0 or gaussian_model.sig < 0):
        discard_gaussian = True
    
    # Mu (peak) inside the spectra region & Coherent parameters for an absorption line
    if not discard_voigt and (voigt_model.mu < spectra_slice['waveobs'][0] or voigt_model.mu > spectra_slice['waveobs'][-1] or voigt_model.A >= 0 or voigt_model.sig < 0 or voigt_model.gamma < 0):
        discard_voigt = True
    
    if (not discard_gaussian and not discard_voigt and rms_gaussian <= rms_voigt) or (not discard_gaussian and discard_voigt):
        return gaussian_model, rms_gaussian, mu_residual_gaussian
    elif (not discard_gaussian and not discard_voigt and rms_gaussian > rms_voigt) or (discard_gaussian and not discard_voigt):
        return voigt_model, rms_voigt, mu_residual_voigt
    else:
        raise Exception("Gaussian or Voigt fit for absorption line not possible.")
    
    #if (not discard_gaussian and not discard_voigt and rms_gaussian + mu_residual_gaussian <= rms_voigt + mu_residual_voigt) or (not discard_gaussian and discard_voigt):
        #return gaussian_model, rms_gaussian, mu_residual_gaussian
    #elif (not discard_gaussian and not discard_voigt and rms_gaussian + mu_residual_gaussian > rms_voigt + mu_residual_voigt) or (discard_gaussian and not discard_voigt):
        #return voigt_model, rms_voigt, mu_residual_voigt
    #else:
        #raise Exception("Gaussian or Voigt fit for absorption line not possible.")
    

def calculate_mu_residual(spectra, line_model):
    ## Mu residual (sometimes RMS is very good but the fit is bad and it can be detected by
    ##              evaluating the flux at the mu position)
    ##            * Maximum value: 9999.0
    # - Modeled flux at mu
    modeled_mu_flux = line_model(line_model.mu)
    # - Observed flux at mu (interpolate if needed)
    observed_mu_flux = np.interp(line_model.mu, spectra['waveobs'], spectra['flux'] - continuum_model(spectra['waveobs']))
    return np.abs(modeled_mu_flux - observed_mu_flux)

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
def assert_structure(spectra, peaks, base_points):
    # Limit the base_points array to the ones that are useful, considering that
    # the first and last peak are always removed
    first_wave_peak = spectra['waveobs'][peaks][0]
    first_wave_base = spectra['waveobs'][base_points][0]
    if first_wave_peak > first_wave_base:
        if len(base_points) - len(peaks) == 1:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: base point
            base_points = base_points
            peaks = peaks
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: peak (this last one will be removed)
            base_points = base_points
            peaks = peaks[:-1]
        else:
            raise Exception("This should not happen")
    else:
        if len(base_points) - len(peaks) == -1:
            ## First feature found in spectra: peak (this first one will be removed)
            ## Last feature found in spectra: peak (this last one will be removed)
            base_points = base_points
            peaks = peaks[1:-1]
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: peak (this first one will be removed)
            ## Last feature found in spectra: base point
            base_points = base_points
            peaks = peaks[1:]
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
# If prioritize_deeper_fluxes is True, values with lower flux values will
# have more weight in the fitting process
# To save computation time, min and max depth can be indicated and all the lines out
# of this range will not be considered for fit process (save CPU time) although
# the rest of the information of the line will be conserved in the output
# Returns a complete structure with all the necessary information to
# determine if it is a line of interest
def generate_linemasks(spectra, peaks, base_points, continuum_model, minimum_depth=None, maximum_depth=None, smoothed_spectra=None, vald_linelist_file="input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", discard_gaussian = False, discard_voigt = False, prioritize_deeper_fluxes=False):
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

    num_peaks = len(peaks)
    linemasks = np.recarray((num_peaks, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100'), ('peak', int), ('base', int), ('top', int), ('depth', float), ('relative_depth', float), ('wave_base_fit', float), ('wave_top_fit', float), ('base_fit', int), ('top_fit', int), ('mu', float), ('sig', float), ('A', float), ('gamma', float), ('depth_fit', float), ('relative_depth_fit', float), ('integrated_flux', float), ('ew', float), ('rms', float), ('mu_residual', float), ('wave_peak_diff', float), ('element', '|S4'), ('lower_state(eV)', float), ('log(gf)', float), ('solar_depth', float), ('discarded', bool)])
    
    for i in np.arange(num_peaks):
        linemasks['discarded'][i] = False
        # Line mask
        linemasks['wave_peak'][i] = spectra['waveobs'][peaks[i]]
        linemasks['wave_base'][i] = spectra['waveobs'][base_points[i]]
        linemasks['wave_top'][i] = spectra['waveobs'][base_points[i+1]]
        linemasks['note'][i] = ""
        # Position in the spectra, usefull for noise detection
        linemasks['peak'][i] = peaks[i]
        linemasks['base'][i] = base_points[i]
        linemasks['top'][i] = base_points[i+1]
        linemasks['depth'][i] = depth[i]
        # Relative depth is "peak - mean_base_point" with respect to the total continuum
        # - In case that the mean base point is higher than the continuum, relative_depth < 0
        # - relative_depth < depth is always true
        flux_from_top_base_point_to_continuum = np.abs(continuum_at_peak[i] - np.mean([spectra['flux'][base_points[i]], spectra['flux'][base_points[i+1]]]))
        linemasks['relative_depth'][i] = ((continuum_at_peak[i] - (flux_at_peak[i] + flux_from_top_base_point_to_continuum)) / continuum_at_peak[i])
        # Model: fit gaussian
        # Adjust edges
        new_base, new_top = improve_linemask_edges(smoothed_spectra, linemasks['base'][i], linemasks['top'][i], linemasks['peak'][i])
        linemasks['base_fit'][i] = new_base
        linemasks['top_fit'][i] = new_top
        linemasks['wave_base_fit'][i] = spectra['waveobs'][new_base]
        linemasks['wave_top_fit'][i] = spectra['waveobs'][new_top]
       
        fitting_not_possible = False
        if accepted_for_fitting[i]:
            try:
                line_model, rms, mu_residual = fit_line(spectra[new_base:new_top+1], continuum_model, linemasks['wave_peak'][i], sig=0.02, A=-0.025, gamma=0.025, discard_gaussian = discard_gaussian, discard_voigt = discard_voigt, prioritize_deeper_fluxes = prioritize_deeper_fluxes)
                linemasks['mu'][i] = line_model.mu
                linemasks['sig'][i] = line_model.sig
                linemasks['A'][i] = line_model.A
                if type(line_model) == VoigtModel:
                    linemasks['gamma'][i] = line_model.gamma
                else:
                    # The model is Gaussian, do not use 'gamma'
                    linemasks['gamma'][i] = 9999.0
                
                # Depth of the peak with respect to the total continuum in % over the total continuum
                # - In case that the peak is higher than the continuum, depth < 0
                continuum = continuum_model(line_model.mu)
                flux = line_model(line_model.mu) + continuum
                linemasks['depth_fit'][i] = ((continuum - flux) / continuum)
                # Relative depth is "peak - mean_base_point" with respect to the total continuum
                # - In case that the mean base point is higher than the continuum, relative_depth < 0
                # - relative_depth < depth is always true
                flux_from_top_base_point_to_continuum = np.abs(continuum - np.mean([spectra['flux'][base_points[i]], spectra['flux'][base_points[i+1]]]))
                linemasks['relative_depth_fit'][i] = ((continuum - (flux + flux_from_top_base_point_to_continuum)) / continuum)
                
                # Equivalent Width
                linemasks['integrated_flux'][i] = -1 * line_model.integrate(spectra['waveobs'][new_base], spectra['waveobs'][new_top])
                linemasks['ew'][i] = linemasks['integrated_flux'][i] / np.mean(continuum_model(spectra['waveobs'][new_base:new_top+1]))
                # RMS
                linemasks['rms'][i] = rms
                # Mu residual
                # - Residual of the flux at the mu position, for observed spectra interpolation is used if needed
                # - Maximum value (to avoid overflows): 9999.0
                linemasks['mu_residual'][i] = mu_residual
            except Exception as e:
                #print "WARNING: Bad line fit (", i, ") - ", e.message
                fitting_not_possible = True
        
        if not accepted_for_fitting[i] or fitting_not_possible:
            # Values that indicate that the line has not been fitted
            linemasks['mu'][i] = 0.0
            linemasks['sig'][i] = 0.0
            linemasks['A'][i] = 0.0
            linemasks['gamma'][i] = 0.0
            linemasks['depth_fit'][i] = 0.0
            linemasks['relative_depth_fit'][i] = 0.0
            linemasks['integrated_flux'][i] = 0.0
            linemasks['ew'][i] = 0.0

        if not accepted_for_fitting[i] or fitting_not_possible:
            linemasks['rms'][i] = 0.0
            linemasks['mu_residual'][i] = 0.0
        elif fitting_not_possible:
            linemasks['rms'][i] = 9999.0
            linemasks['mu_residual'][i] = 9999.0
        
        if (i % 100) == 0:
            print "%.2f%%" % (((i*1.0)/num_peaks) * 100)
        
    linemasks = fill_with_VALD_info(linemasks, vald_linelist_file=vald_linelist_file)
    
    return linemasks

# Cross-match linemasks with a VALD linelist in order to find
# the nearest lines and copy the information into the linemasks structure
def fill_with_VALD_info(linemasks, vald_linelist_file="input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst"):
        ## Load original VALD linelist
        vald_linelist = load_and_filter_VALD(vald_linelist_file, depth_limit=0.0)
        
        ## Convert wavelengths from armstrong to nm
        vald_linelist['wave (A)'] = vald_linelist['wave (A)'] / 10.0
        vald_linelist = rfn.rename_fields(vald_linelist, {'wave (A)':'wave_peak',})
        
        ## Detect duplicates
        # Sort by wave_peak and descending depth (for that we create a temporary field)
        vald_linelist = rfn.append_fields(vald_linelist, "reverse_depth", dtypes=float, data=1-vald_linelist['depth'])
        vald_linelist.sort(order=['wave_peak', 'reverse_depth'])
        vald_linelist = rfn.drop_fields(vald_linelist, ['reverse_depth'])
        # Find duplicates
        dups, dups_index = find_duplicates(vald_linelist, 'wave_peak')
        # Filter all duplicates except the first one (which corresponds to the biggest depth)
        valid = ~np.isnan(vald_linelist['wave_peak'])
        last_wave = None
        for i in np.arange(len(dups)):
            current_wave = dups[i]['wave_peak']
            if last_wave == None:
                last_wave = current_wave
                continue
            if last_wave == current_wave:
                pos = dups_index[i]
                valid[pos] = False
            else:
                # Do not filter the first duplicated value
                last_wave = dups[i]['wave_peak']
        # Remove duplicates, leaving only those with the biggest depth
        vald_linelist = vald_linelist[valid]
        
        if vald_linelist['wave_peak'][0] > linemasks['wave_peak'][0] or vald_linelist['wave_peak'][-1] < linemasks['wave_peak'][-1]:
            print "WARNING: VALD linelist does not cover the linemask wavelength range"
            print "- VALD range from", vald_linelist['wave_peak'][0], "to", vald_linelist['wave_peak'][-1], "nm"
            print "- Linemask range from", linemasks['wave_peak'][0], "to", linemasks['wave_peak'][-1], "nm"
        
        # Cross-match the linemask with the VALD linelist by finding the closest line
        # - There is no need to sort, linemask is already sorted by wavelength because of the way it is created
        #linemasks.sort(order=['wave_peak'])
        last_diff = None
        i = 0
        j = 0
        while i < len(linemasks) and j < len(vald_linelist):
            diff = linemasks['wave_peak'][i] - vald_linelist['wave_peak'][j]
            if last_diff == None:
                last_diff = diff
                j += 1
                continue
            if np.abs(last_diff) < np.abs(diff):
                # Save the information
                linemasks["wave_peak_diff"][i] = last_diff
                linemasks["element"][i] = vald_linelist["element"][j-1]
                linemasks["lower_state(eV)"][i] = vald_linelist["lower state (eV)"][j-1]
                linemasks["log(gf)"][i] = vald_linelist["log(gf)"][j-1]
                linemasks["solar_depth"][i] = vald_linelist["depth"][j-1]
                last_diff = None
                i += 1 # Next line in the linemask
                j = np.max([0, j-10]) # Jump back in the VALD linelist, just in case
            else:
                last_diff = diff
                j += 1
        
        return linemasks

# Works better with a smoothed spectra (i.e. convolved using 2*resolution)
def find_peaks_and_base_points(spectra, use_EW_SNR=False, resolution=None):
    # Determine peaks and base points (also known as continuum points)
    if use_EW_SNR:
        if resolution == None:
            raise Exception("Resolution is needed")
        # - SNR: Schneider et al (1993). The Hubble Space Telescope quasar absorption line key project. II - Data calibration and absorption-line selection
        corrected_convolved_ew_snr = get_corrected_convolved_ew_snr(spectra, resolution)
        peaks = find_local_max_values(corrected_convolved_ew_snr)
        base_points = find_local_min_values(corrected_convolved_ew_snr)
    else:
        ## - Convolved spectrum
        peaks = find_local_min_values(spectra['flux'])
        base_points = find_local_max_values(spectra['flux'])
    
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
    peaks, base_points = assert_structure(spectra, peaks, base_points)
    
    if len(peaks[peaks - base_points[:-1] <= 0]):
        raise Exception("This should not happen")
    
    if use_EW_SNR:
        return peaks, base_points, corrected_convolved_ew_snr
    else:
        return peaks, base_points


def build_fitted_spectrum(waveobs, continuum_model, lines):
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
                line_flux = A*gamma/pi/(waveobs*waveobs - 2*waveobs*mu+mu*mu+gamma*gamma)
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
def improve_linemask_edges(spectra, base, top, peak):
    # Try to use two additional position, since we are going to lose them
    # by doing the first and second derivative
    original_top = top
    top = np.min([top+2, len(spectra)]) 
    flux = spectra[base:top+1]['flux']
    waveobs = spectra[base:top+1]['waveobs']
    # First derivative (positive => flux increases, negative => flux decreases)
    dflux_dwave = (flux[:-1] - flux[1:]) / (waveobs[:-1] - waveobs[1:])
    # Second derivative (positive => convex, negative => concave)
    d2flux_dwave2 = (dflux_dwave[:-1] - dflux_dwave[1:])  / (waveobs[:-2] - waveobs[2:])
    # Peak position inside the linemask region
    peak_relative_pos = peak - base
    # The peak should be in a convex region => the second derivative should be positive
    # - It may happen that the peak falls in the beginning/end of a concave region, accept also this cases
    if peak_relative_pos < len(d2flux_dwave2) and (d2flux_dwave2[peak_relative_pos-1] > 0 or d2flux_dwave2[peak_relative_pos] > 0 or d2flux_dwave2[peak_relative_pos+1] > 0):
        # Find the concave positions at both sides of the peak
        concave_pos = np.where(d2flux_dwave2<0)[0]
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
        raise Exception("This should not happen")
        #plt.plot(waveobs, flux)
        #l = plt.axvline(x = waveobs[peak_relative_pos], linewidth=1, color='red')
        #print d2flux_dwave2, d2flux_dwave2[peak_relative_pos]
        #plt.show()
    
    return new_base, new_top

def fwhm_and_resolution(linemasks):
    # Profile types
    gaussian = linemasks['gamma'] == 9999.0
    voigt = linemasks['gamma'] != 9999.0
    ## Resolution
    # Light speed in vacuum
    c = 299792458.0 # m/s
    sigma = linemasks['sig']
    # FWHM considering all gaussians
    fwhm = sigma * (2*np.sqrt(2*np.log(2))) # nm
    # FWHM for lorentzian (temporary value needed for voigt FWHM calculation)
    fwhm_lorentzian = 2*linemasks['gamma'][voigt]
    # FWHM for voigt
    # Formula from Olivero et al (1977) "Empirical fits to the Voigt line width: A brief review"
    # http://www.sciencedirect.com/science/article/pii/0022407377901613
    # http://en.wikipedia.org/wiki/Voigt_profile#The_width_of_the_Voigt_profile
    fwhm[voigt] = 0.5346*(2*fwhm_lorentzian) + np.sqrt(0.2166*np.power(fwhm_lorentzian, 2) + np.power(fwhm[voigt], 2))
    
    R = np.zeros(len(fwhm)) # Avoid zero division
    R[fwhm != 0] = linemasks['wave_peak'][fwhm != 0] / fwhm[fwhm != 0]
    # In m/s:
    #fwhm = c / R # m/s
    # ... or ...
    #fwhm = c * (fwhm / linemasks['wave_peak'][~discarded]) # m/s
    fwhm_kms = np.zeros(len(fwhm)) # Avoid zero division
    fwhm_kms[R != 0] = c / R[R != 0] # m/s
    fwhm_kms = fwhm_kms / 1000 # km/s

    return fwhm, fwhm_kms, R
 


def print_linemasks_stats(linemasks, discarded):
    # Profile types
    gaussian = linemasks['gamma'] == 9999.0
    voigt = linemasks['gamma'] != 9999.0
    fwhm, fwhm_kms, R = fwhm_and_resolution(linemasks)
    
    print "---------------GENERAL----------------"
    print "\t\tMean\tMedian\tStdev"
    print "RMS:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['rms'][~discarded]), np.median(linemasks['rms'][~discarded]), np.std(linemasks['rms'][~discarded]))
    print "Mu residual:\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['mu_residual'][~discarded]), np.median(linemasks['mu_residual'][~discarded]), np.std(linemasks['mu_residual'][~discarded]))
    print "A:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['A'][~discarded]), np.median(linemasks['A'][~discarded]), np.std(linemasks['A'][~discarded]))
    print "sig:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['sig'][~discarded]), np.median(linemasks['sig'][~discarded]), np.std(linemasks['sig'][~discarded]))
    print "gamma:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['gamma'][voigt & ~discarded]), np.median(linemasks['gamma'][voigt & ~discarded]), np.std(linemasks['gamma'][voigt & ~discarded]))
    print "R:\t\t%i\t%i\t%i" % (np.median(R[~discarded]), np.mean(R[~discarded]), np.std(R[~discarded]))
    print "FWHM (nm):\t%.3f\t%.3f\t%.3f" % (np.median(fwhm[~discarded]), np.mean(fwhm[~discarded]), np.std(fwhm[~discarded]))
    print "FWHM (km/s):\t%.3f\t%.3f\t%.3f" % (np.median(fwhm_kms[~discarded]), np.mean(fwhm_kms[~discarded]), np.std(fwhm_kms[~discarded]))
    if len(linemasks[gaussian & ~discarded]) > 0:
        print "---------------GAUSSIAN----------------"
        print "\t\tMean\tMedian\tStdev"
        print "RMS:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['rms'][gaussian & ~discarded]), np.median(linemasks['rms'][gaussian & ~discarded]), np.std(linemasks['rms'][gaussian & ~discarded]))
        print "Mu residual:\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['mu_residual'][gaussian & ~discarded]), np.median(linemasks['mu_residual'][gaussian & ~discarded]), np.std(linemasks['mu_residual'][gaussian & ~discarded]))
        print "A:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['A'][gaussian & ~discarded]), np.median(linemasks['A'][gaussian & ~discarded]), np.std(linemasks['A'][gaussian & ~discarded]))
        print "sig:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['sig'][gaussian & ~discarded]), np.median(linemasks['sig'][gaussian & ~discarded]), np.std(linemasks['sig'][gaussian & ~discarded]))
        print "R:\t\t%i\t%i\t%i" % (np.median(R[gaussian & ~discarded]), np.mean(R[gaussian & ~discarded]), np.std(R[gaussian & ~discarded]))
        print "FWHM (nm):\t%.3f\t%.3f\t%.3f" % (np.median(fwhm[gaussian & ~discarded]), np.mean(fwhm[gaussian & ~discarded]), np.std(fwhm[gaussian & ~discarded]))
        print "FWHM (km/s):\t%.3f\t%.3f\t%.3f" % (np.median(fwhm_kms[gaussian & ~discarded]), np.mean(fwhm_kms[gaussian & ~discarded]), np.std(fwhm_kms[gaussian & ~discarded]))
    if len(linemasks[voigt & ~discarded]) > 0:
        print "---------------VOIGT-------------------"
        print "\t\tMean\tMedian\tStdev"
        print "RMS:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['rms'][voigt & ~discarded]), np.median(linemasks['rms'][voigt & ~discarded]), np.std(linemasks['rms'][voigt & ~discarded]))
        print "Mu residual:\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['mu_residual'][voigt & ~discarded]), np.median(linemasks['mu_residual'][voigt & ~discarded]), np.std(linemasks['mu_residual'][voigt & ~discarded]))
        print "A:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['A'][voigt & ~discarded]), np.median(linemasks['A'][voigt & ~discarded]), np.std(linemasks['A'][voigt & ~discarded]))
        print "sig:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['sig'][voigt & ~discarded]), np.median(linemasks['sig'][voigt & ~discarded]), np.std(linemasks['sig'][voigt & ~discarded]))
        print "gamma:\t\t%.3f\t%.3f\t%.3f" % (np.mean(linemasks['gamma'][voigt & ~discarded]), np.median(linemasks['gamma'][voigt & ~discarded]), np.std(linemasks['gamma'][voigt & ~discarded]))
        print "R:\t\t%i\t%i\t%i" % (np.median(R[voigt & ~discarded]), np.mean(R[voigt & ~discarded]), np.std(R[voigt & ~discarded]))
        print "FWHM (nm):\t%.3f\t%.3f\t%.3f" % (np.median(fwhm[voigt & ~discarded]), np.mean(fwhm[voigt & ~discarded]), np.std(fwhm[voigt & ~discarded]))
        print "FWHM (km/s):\t%.3f\t%.3f\t%.3f" % (np.median(fwhm_kms[voigt & ~discarded]), np.mean(fwhm_kms[voigt & ~discarded]), np.std(fwhm_kms[voigt & ~discarded]))
    

if __name__ == '__main__':
    ### VALD line list
    #VALD_to_SPECTRUM_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/linelists/VALD.300_1100nm.lst", depth_limit=0.0)

    #VALD_top_3_to_RV_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/rv/VALD.300_1100nm.rv.lst", top=1, wave_step=10)
    
    ### Line detection and fitting example:
    print "Reading spectrum..."
    #star, resolution = "/home/marble/Downloads/HD146233/Narval/narval_hd146233_100312.s", 75000
    #star, resolution = "input/L082N03_spec_norm/05oct08/sp2_Normal/004_vesta_001.s", 75000
    #star, resolution = "/home/marble/Downloads/HD146233/UVES/hd146233_uves.txt", 47000
    #star, resolution = "input/test/observed_arcturus.s.gz", 47000
    #star, resolution = "input/test/observed_mu_cas_a.s.gz", 47000
    #star, resolution = "input/test/observed_mu_leo.s.gz", 47000
    #star, resolution = "input/test/observed_sun.s.gz", 47000
    ## Test
    #star, resolution = "/home/marble/Downloads/HD146233/Narval/narval_hd146233_100312_segment.s", 65000
    #star, resolution = "input/test/observed_sun.s.gz", 47000
    #star, resolution = "/home/marble/Downloads/HD146233/UVES/hd146233_uves_segment.txt", 47000
    #star, resolution = "/home/marble/Downloads/HD146233/UVES/hd146233_uves_segment_not_normalized.txt", 47000
    
    ## Smooth spectra using the instrumental resolution
    #normalize = True
    normalize = False
    smooth_spectra = True
    nknots_factor = 1
    
    find_continuum_by_interpolating_base_points = True
    smooth_continuum_interpolation = True
    find_continuum_by_spline_fitting = np.logical_not(find_continuum_by_interpolating_base_points)
    
    #############################
    # Telluric lines: constant continuum at 1, no filtering
    star, resolution = "input/telluric/standard_atm_air.s.gz", 100000
    ## Smooth spectra using the instrumental resolution
    normalize = False
    smooth_spectra = False
    nknots_factor = 3
    
    find_continuum_by_interpolating_base_points = True
    smooth_continuum_interpolation = False
    find_continuum_by_spline_fitting = np.logical_not(find_continuum_by_interpolating_base_points)
    #############################
    
    ## Reading spectra
    spectra = read_spectra(star)
    
    if normalize:
        print "Continuum normalization..."
        continuum_base_points = determine_continuum_base_points(spectra)
        if find_continuum_by_interpolating_base_points:
            # OPTION A: DIRECT LINEAR INTERPOLATION
            continuum_model_for_normalization = interpolate_continuum(spectra, continuum_base_points, smooth_continuum_interpolation)
        else:
            # OPTION B: SPLINE FITTING
            # * 1 knot every 10 nm in average
            nknots = np.max([1, int((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 10)])
            nknots = nknots_factor * nknots # N knots every 10nm in average
            continuum_model_for_normalization = fit_continuum(spectra, nknots=nknots)
        spectra['flux'] /= continuum_model_for_normalization(spectra['waveobs'])
    
    # Use method from Schneider et al (1993).
    #   The Hubble Space Telescope quasar absorption line key project. 
    #   II - Data calibration and absorption-line selection
    use_EW_SNR = False

    ### Fitting parameters
    discard_gaussian = False
    discard_voigt = False
    prioritize_deeper_fluxes_for_fit = False # False needed for telluric
    generate_fitted_spectra = True
    
    ### Filtering parameters
    # Discard lines that do not have at least a given depth
    minimum_depth = 0.00 # (% of the continuum)
    maximum_depth = 1.00 # (% of the continuum)
    # Discard outliers of the VALD cross-matched lines
    discard_bad_VALD_crossmatch = False
    # Discard potential gaps in the spectra that are identified as lines
    discard_too_big_wavelength_range = False
    # Discard outliers (biggest RMS or mu_residual)
    discard_outlier_fit = False
    # Discard outliers in resolving power
    discard_outlier_R = False 
    
    print "Convolving (smoothing)..."
    original_spectra = spectra
    if smooth_spectra:
        # Use 2 times the resolution to smooth the spectra without loosing too much details
        spectra = convolve_spectra(original_spectra, 2*resolution)
    
    print "Finding peaks and base points..."
    if use_EW_SNR:
        peaks, base_points, corrected_convolved_ew_snr = find_peaks_and_base_points(spectra, use_EW_SNR=True)
    else:
        peaks, base_points = find_peaks_and_base_points(spectra, use_EW_SNR=False)
    
    # Fit continuum
    print "Fitting continuum..."
    continuum_base_points = determine_continuum_base_points(spectra)
    if find_continuum_by_interpolating_base_points:
        # OPTION A: DIRECT LINEAR INTERPOLATION
        continuum_model = interpolate_continuum(spectra, continuum_base_points, smooth_continuum_interpolation)
    else:
        # OPTION B: SPLINE FITTING
        # * 1 knot every 10 nm in average
        nknots = np.max([1, int((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 10)])
        nknots = nknots_factor * nknots # N knots every 10nm in average
        continuum_model = fit_continuum(spectra, nknots=nknots)
    
    
    print "Generating linemasks and fitting gaussians..."
    #linemasks = generate_linemasks(original_spectra, peaks, base_points, continuum_model, smoothed_spectra=spectra ,vald_linelist_file="input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", discard_gaussian = discard_gaussian, discard_voigt = discard_voigt, prioritize_deeper_fluxes = prioritize_deeper_fluxes_for_fit)
    linemasks = generate_linemasks(original_spectra, peaks, base_points, continuum_model, minimum_depth=minimum_depth, maximum_depth=maximum_depth, smoothed_spectra=spectra ,vald_linelist_file="input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", discard_gaussian = discard_gaussian, discard_voigt = discard_voigt, prioritize_deeper_fluxes = prioritize_deeper_fluxes_for_fit)
    
    ##### Filters
    rejected_by_noise = detect_false_positives_and_noise(spectra, linemasks)

    if use_EW_SNR:
        # Filter low SNR peaks
        rejected_by_low_snr = (corrected_convolved_ew_snr[peaks] < np.median(corrected_convolved_ew_snr))
        rejected_by_noise = np.logical_or(rejected_by_noise, rejected_by_low_snr)
    
    # Identify outliers in the cross-matching with the VALD line list
    wave_diff_selected, wave_diff_selected_filter = sigma_clipping(linemasks["wave_peak_diff"], sig=3, meanfunc=np.median) # Discard outliers
    rejected_by_bad_VALD_crossmatch = np.logical_not(wave_diff_selected_filter)

    # Identify linemasks with too big wavelength range (outliers)
    # WARNING: If the spectrum has a small wavelength range with some few strong features
    #          they may be erroneously discarded by this filter
    wave_diff = (linemasks['wave_top'] - linemasks['wave_base']) / (linemasks['top'] - linemasks['base'])
    wave_diff_selected, wave_diff_selected_filter = sigma_clipping(wave_diff, sig=3, meanfunc=np.median) # Discard outliers
    rejected_by_wave_gaps = np.logical_not(wave_diff_selected_filter)
    
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
    rejected_by_bad_fit = np.logical_or(linemasks['rms'] >= 9999.0, linemasks['mu_residual'] >= 9999.0)
    # Also it is a bad fit if the mu residual is bigger than the flux difference between
    # the highest and lowest point
    lower_flux = np.min(original_spectra['flux'][linemasks['peak']])
    higher_flux = np.max(original_spectra['flux'][linemasks['base']])
    maximum_flux_difference = np.abs(higher_flux - lower_flux)
    rejected_by_bad_fit = np.logical_or(rejected_by_bad_fit, np.abs(linemasks['mu_residual']) > maximum_flux_difference)
    
    
    discarded = rejected_by_noise
    if discard_bad_VALD_crossmatch:
        discarded = np.logical_or(discarded, rejected_by_bad_VALD_crossmatch)
    if discard_too_big_wavelength_range:
        discarded = np.logical_or(discarded, rejected_by_wave_gaps)
    discarded = np.logical_or(discarded, rejected_by_depth_higher_than_continuum)
    discarded = np.logical_or(discarded, rejected_by_depth_limits)
    discarded = np.logical_or(discarded, rejected_by_bad_fit)

    ### Outliers identification should be done avoiding the already discarded points
    ### otherwise, we will have noise contamination
    # Identify outliers considering RMS and mu_residual
    rms = (linemasks['rms'] + linemasks['mu_residual']) / 2
    rms_selected, rms_selected_filter = sigma_clipping(rms[~discarded], sig=3, meanfunc=np.median) # Discard outliers
    accepted_index = np.arange(len(rms))[~discarded][rms_selected_filter]
    rejected_by_outlier_fit = rms < 0 # Create an array of booleans
    rejected_by_outlier_fit[:] = True # Initialize
    rejected_by_outlier_fit[accepted_index] = False
    
    if discard_outlier_fit:
        discarded = np.logical_or(discarded, rejected_by_outlier_fit)
    
    # Identify outliers considering the resolution
    fwhm, fwhm_kms, R = fwhm_and_resolution(linemasks)
    r_selected, r_selected_filter = sigma_clipping(R[~discarded], sig=3, meanfunc=np.median) # Discard outliers
    accepted_index = np.arange(len(R))[~discarded][r_selected_filter]
    rejected_by_outlier_R = R < 0 # Create an array of booleans
    rejected_by_outlier_R[:] = True # Initialize
    rejected_by_outlier_R[accepted_index] = False

    if discard_outlier_R:
        discarded = np.logical_or(discarded, rejected_by_outlier_R)
    
    linemasks['discarded'][discarded] = True
    
    # Profile types
    gaussian = linemasks['gamma'] == 9999.0
    voigt = linemasks['gamma'] != 9999.0


    ## Peaks
    print "--------------------------------------"
    print "Number of peak candidates:\t", len(linemasks)
    print "- Noise and false positive:\t", len(linemasks[rejected_by_noise])
    print "- Peaks higher than continuum:\t", len(linemasks[rejected_by_depth_higher_than_continuum])
    if prioritize_deeper_fluxes_for_fit:
        print "- Bad fits:\t\t\t", len(linemasks[rejected_by_bad_fit]), "\t[Deeper fluxes prioritized]"
    else:
        print "- Bad fits:\t\t\t", len(linemasks[rejected_by_bad_fit])
    print "- Out of depth limits:\t\t", len(linemasks[rejected_by_depth_limits]), "\t[%.2f - %.2f]" % (minimum_depth, maximum_depth)
    if discard_bad_VALD_crossmatch:
        print "- Bad VALD cross-match:\t\t", len(linemasks[rejected_by_bad_VALD_crossmatch])
    else:
        print "- Bad VALD cross-match:\t\t", len(linemasks[rejected_by_bad_VALD_crossmatch]), "\t[Disabled]"
    if discard_too_big_wavelength_range:
        print "- Too big wavelength range:\t", len(linemasks[rejected_by_wave_gaps])
    else:
        print "- Too big wavelength range:\t", len(linemasks[rejected_by_wave_gaps]), "\t[Disabled]"
    if discard_outlier_fit:
        print "- Outlier fits:\t\t\t", len(linemasks[rejected_by_outlier_fit])
    else:
        print "- Outlier fits:\t\t\t", len(linemasks[rejected_by_outlier_fit]), "\t[Disabled]"
    if discard_outlier_R:
        print "- Resolving power outliers:\t", len(linemasks[rejected_by_outlier_R])
    else:
        print "- Resolving power outliers:\t", len(linemasks[rejected_by_outlier_R]), "\t[Disabled]"
    print "Final number of peaks:\t\t", len(linemasks[~discarded])
    print "- Gaussian profile:\t\t", len(linemasks[gaussian & ~discarded])
    print "- Voigt profile:\t\t", len(linemasks[voigt & ~discarded])
    print_linemasks_stats(linemasks, discarded)
    
    if generate_fitted_spectra:
        print "Building a fitted spectrum..."
        #waveobs = generate_wavelength_grid(np.min(spectra['waveobs']), np.max(spectra['waveobs']), resolution, points_per_fwhm = 3)
        waveobs = spectra['waveobs']
        fitted_spectra = build_fitted_spectrum(waveobs, continuum_model, linemasks[~discarded])
        
        # Plot
        if use_EW_SNR:
            fig = plt.figure()
            axes = fig.add_subplot(211)
            plt.plot(spectra['waveobs'], original_spectra['flux'])
            plt.plot(spectra['waveobs'], continuum_model(spectra['waveobs']))
            plt.plot(fitted_spectra['waveobs'], fitted_spectra['flux'])
            plt.scatter(spectra['waveobs'][peaks], spectra['flux'][peaks], c='red')
            plt.scatter(spectra['waveobs'][base_points], spectra['flux'][base_points], c='green')
            axes = fig.add_subplot(212, sharex=axes)
            plt.plot(spectra['waveobs'], corrected_convolved_ew_snr) # green
            plt.show()
        else:
            fig = plt.figure()
            plt.plot(spectra['waveobs'], original_spectra['flux'])
            plt.plot(spectra['waveobs'], continuum_model(spectra['waveobs']))
            plt.plot(fitted_spectra['waveobs'], fitted_spectra['flux'])
            plt.scatter(spectra['waveobs'][peaks], original_spectra['flux'][peaks], c='red')
            plt.scatter(spectra['waveobs'][base_points], original_spectra['flux'][base_points], c='green')
            plt.show()
    else:
        fitted_spectra = None

    ##### Saving...
    dump_filename = "output/" + os.path.basename(star) + ".dump"
    # Complete dump
    version = 20120402
    data = (os.path.abspath(star), resolution, continuum_model.data, continuum_model.pardict, linemasks, fitted_spectra)
    pickle.dump((version, data), gzip.open(dump_filename, "wb", compresslevel=3), protocol=2)

    ## Complete txt line list
    #asciitable.write(linemasks, output=output_filename+"_complete.txt", delimiter="\t")
    ## Masks regions for visualization
    #linemasks['note'][~discarded] = linemasks['element'][~discarded]
    #linemasks['wave_peak'][~discarded] = linemasks['mu'][~discarded]
    #linemasks['wave_base'][~discarded] = linemasks['wave_base_fit'][~discarded]
    #linemasks['wave_top'][~discarded] = linemasks['wave_top_fit'][~discarded]
    #asciitable.write(linemasks[~discarded], output=output_filename+".txt", delimiter="\t", include_names=['wave_peak', 'wave_base', 'wave_top', 'note'])
    ##### Restoring...
    dump_filename = "output/" + os.path.basename(star) + ".dump"
    version, data = pickle.load(gzip.open(dump_filename, "rb"))
    continuum_model = UniformCDFKnotSplineModel()
    star, resolution, continuum_model.data, continuum_model.pardict, linemasks, fitted_spectra = data
    ##################

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['wave_peak'], fwhm[~discarded])
    #plt.scatter(linemasks[~discarded]['wave_peak'], fwhm[~discarded], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('FWHM (nm)')
    #plt.show()
    
    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['wave_peak'], fwhm_kms[~discarded])
    #plt.scatter(linemasks[~discarded]['wave_peak'], fwhm_kms[~discarded], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('FWHM (km/s)')
    #plt.show()
    
    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[~discarded]['wave_peak'], R[~discarded])
    #plt.scatter(linemasks[~discarded]['wave_peak'], R[~discarded], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()
    
    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[gaussian & ~discarded]['wave_peak'], R[gaussian & ~discarded])
    #plt.scatter(linemasks[gaussian & ~discarded]['wave_peak'], R[gaussian & ~discarded], s=4)
    #plt.plot(linemasks[gaussian & ~discarded]['wave_peak'], trend(linemasks[gaussian & ~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()

    #show_histogram(R[~discarded][R[~discarded] < 100000])

    #from pymodelfit import LinearModel
    #trend = LinearModel()
    #trend.fitData(linemasks[gaussian & ~discarded]['depth'], R[gaussian & ~discarded])
    #plt.scatter(linemasks[gaussian & ~discarded]['depth'], R[gaussian & ~discarded], s=4)
    #plt.plot(linemasks[gaussian & ~discarded]['depth'], trend(linemasks[gaussian & ~discarded]['depth']), color="red", linewidth=3)
    #plt.xlabel('Depth')
    #plt.ylabel('Resolving power')
    #plt.show()

    from pymodelfit import LinearModel
    plt.scatter(R[gaussian & ~discarded], linemasks[gaussian & ~discarded]['depth'], s=4)
    plt.ylabel('Depth')
    plt.xlabel('Resolving power')
    plt.show()

    ############ STATS
    #show_histogram(np.abs(linemasks["wave_peak_diff"]), nbins=100)
    #show_histogram(np.abs(linemasks["depth"] - linemasks["relative_depth"]), nbins=100)
    #show_histogram(linemasks["rms"], nbins=100)
    #show_histogram(linemasks["A"], nbins=100)
    #show_histogram(linemasks["sig"], nbins=100)
    #
    #show_histogram(linemasks["depth"], nbins=100)
    #show_histogram(linemasks["relative_depth"], nbins=100)
    #show_histogram(np.abs(linemasks["depth"] - linemasks["solar_depth"]), nbins=100)
    #show_histogram(linemasks[~discarded]['rms'], nbins=100)
    #show_histogram(linemasks[~discarded]['mu_residual'], xlabel='mu residual', nbins=100)
    
    #show_histogram(linemasks[~discarded]['sig'], xlabel='sig', nbins=100)
    #show_histogram(linemasks[~discarded]['A'], xlabel='A', nbins=100)
    #show_histogram(linemasks[~discarded]['gamma'], xlabel='gamma', nbins=100)
    
    #show_histogram(R[~discarded][R[~discarded] < 200000], xlabel='R', nbins=100)
    
    
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
    
    #plt.plot(linemasks[~discarded]['mu'], np.abs(linemasks[~discarded]['wave_peak_diff']))
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Wavelength diff. with VALD')
    #plt.show()
    
    #plt.scatter(linemasks[~discarded]['depth'], linemasks[~discarded]['log(gf)'], s=4)
    #plt.xlabel('Depth')
    #plt.ylabel('log(gf)')
    #plt.show()

    #plt.scatter(linemasks[~discarded]['depth'], linemasks[~discarded]['lower_state(eV)'], s=4)
    #plt.xlabel('Depth')
    #plt.ylabel('lower_state(eV)')
    #plt.show()
    
    #plt.scatter(linemasks[~discarded]['log(gf)'], linemasks[~discarded]['lower_state(eV)'], s=4)
    #plt.xlabel('log(gf)')
    #plt.ylabel('lower_state(eV)')
    #plt.show()
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot_wireframe(, rstride=10, cstride=10)
    #plt.show()
    
    ##import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization! 
    #fig = plt.figure()
    #ax = Axes3D(fig) #<-- Note the difference from your original code...
    #X, Y, Z = linemasks[~discarded]['log(gf)'], linemasks[~discarded]['lower_state(eV)'], linemasks[~discarded]['depth']
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
