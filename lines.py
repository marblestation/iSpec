import asciitable
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions
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
def fit_gaussian(spectra_slice, continuum_model, mu, sig=0.02, A=-0.025):
    model = GaussianModel()
    model.mu = mu
    model.sig = sig
    model.A = A
    cont = continuum_model(spectra_slice['waveobs'])
    conterr = 0
    
    # More weight to the deeper fluxes
    min_flux = np.min(spectra_slice['flux'])
    if min_flux < 0:
        weights = spectra_slice['flux'] + -1*(min_flux) + 0.01 # Above zero
        weights = np.min(weights) / weights
    else:
        weights = min_flux / spectra_slice['flux']
    
    # Without weigths:
    #model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, fixedpars=[])
    # Priorizing deeper fluxes:
    model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=weights, fixedpars=[])
    # Priorizing fluxes with lower errors:
    #model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=1/spectra_slice['err'], fixedpars=['mu','sig','A'])
    
    # A gaussian with a positive A and a negative Sigma is exactly the same as the inverse
    # - We inverse the signs since it is more intuitive for absorption lines (negative amplitudes)
    #   and it will be easier to apply filtering rules after (such as, discard all gaussian/voigt with positive A)
    if model.A > 0 and model.sig < 0:
        model.A = -1 * model.A
        model.sig = -1 * model.sig
    
    return model

# Fits a voigt at a given wavelength location using a fitted continuum model
# - For absorption lines, it will alway be true:
#      model.A < 0 and model.sig > 0
#   The oposite indicates a potential bad fit
# - For absorption lines, model.gamma < 0 indicates strange wings and probably a bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_voigt(spectra_slice, continuum_model, mu, sig=0.02, A=-0.025, gamma=0.025):
    model = VoigtModel()
    model.gamma = gamma
    model.mu = mu
    model.sig = sig
    model.A = A
    cont = continuum_model(spectra_slice['waveobs'])
    conterr = 0
    
    # More weight to the deeper fluxes
    min_flux = np.min(spectra_slice['flux'])
    if min_flux < 0:
        weights = spectra_slice['flux'] + -1*(min_flux) + 0.01 # Above zero
        weights = np.min(weights) / weights
    else:
        weights = min_flux / spectra_slice['flux']
    
    # Without weigths:
    #model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, fixedpars=[])
    # Priorizing deeper fluxes:
    model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=weights, fixedpars=[])
    # Priorizing fluxes with lower errors:
    #model.fitData(spectra_slice['waveobs'], spectra_slice['flux'] - cont, weights=1/spectra_slice['err'], fixedpars=['mu','sig','A'])
    
    # A voigt with a positive A, a negative Sigma and a negative gamma is exactly the same as the inverse
    # - We inverse the signs since it is more intuitive for absorption lines (negative amplitudes)
    #   and it will be easier to apply filtering rules after (such as, discard all gaussian/voigt with positive A)
    if model.A > 0 and model.sig < 0 and model.gamma < 0:
        model.A = -1 * model.A
        model.sig = -1 * model.sig
        model.gamma = -1 * model.gamma
    
    return model
    

# Fits a gaussian or a voigt at a given wavelength location using a fitted continuum model
# - For absorption lines, it will alway be true:
#      model.A < 0 and model.sig > 0
#   The oposite indicates a potential bad fit
# - For absorption lines, model.gamma < 0 indicates strange wings and probably a bad fit
# - model.mu outside the region used for the fitting it is also a symptom of bad fit
def fit_line(spectra_slice, continuum_model, mu, sig=0.02, A=-0.025, gamma=0.025):
    discard_gaussian = False
    discard_voigt = False
    
    try:
        gaussian_model = fit_gaussian(spectra_slice, continuum_model, mu, sig=sig, A=A)
        rms_gaussian = np.sqrt(np.sum(np.power(gaussian_model.residuals(), 2)) / len(gaussian_model.residuals()))
        mu_residual_gaussian = calculate_mu_residual(spectra_slice, gaussian_model)
    except Exception as e:
        rms_gaussian = 9999.0
        mu_residual_gaussian = 9999.0
        discard_gaussian = True
    
    try:
        voigt_model = fit_voigt(spectra_slice, continuum_model, mu, sig=sig, A=A, gamma=gamma)
        rms_voigt = np.sqrt(np.sum(np.power(voigt_model.residuals(), 2)) / len(voigt_model.residuals()))
        mu_residual_voigt = calculate_mu_residual(spectra_slice, voigt_model)
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
    
    if (not discard_gaussian and not discard_voigt and rms_gaussian + mu_residual_gaussian <= rms_voigt + mu_residual_voigt) or (not discard_gaussian and discard_voigt):
        return gaussian_model
    elif (not discard_gaussian and not discard_voigt and rms_gaussian + mu_residual_gaussian > rms_voigt + mu_residual_voigt) or (discard_gaussian and not discard_voigt):
        return voigt_model
    else:
        raise Exception("Gaussian or Voigt fit for absorption line not possible.")
    

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
# Returns a complete structure with all the necessary information to
# determine if it is a line of interest
def generate_linemasks(spectra, peaks, base_points, continuum_model, smoothed_spectra=None, vald_linelist_file="input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst"):
    if smoothed_spectra == None:
        smoothed_spectra = spectra
    num_peaks = len(peaks)
    linemasks = np.recarray((num_peaks, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100'), ('peak', int), ('base', int), ('top', int), ('wave_base_fit', float), ('wave_top_fit', float), ('base_fit', int), ('top_fit', int), ('mu', float), ('sig', float), ('A', float), ('gamma', float), ('depth', float), ('relative_depth', float), ('flux', float), ('ew', float), ('rms', float), ('mu_residual', float), ('wave_peak_diff', float), ('element', '|S4'), ('lower_state(eV)', float), ('log(gf)', float), ('solar_depth', float), ('discarded', bool)])
    
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
        # Depth of the peak with respect to the total continuum in % over the total continuum
        # - In case that the peak is higher than the continuum, depth < 0
        flux = spectra['flux'][peaks][i]
        continuum = continuum_model(linemasks['wave_peak'][i])
        linemasks['depth'][i] = ((continuum - flux) / continuum)
        # Relative depth is "peak - mean_base_point" with respect to the total continuum
        # - In case that the mean base point is higher than the continuum, relative_depth < 0
        # - relative_depth < depth is always true
        flux_from_top_base_point_to_continuum = np.abs(continuum - np.mean([spectra['flux'][base_points[i]], spectra['flux'][base_points[i+1]]]))
        linemasks['relative_depth'][i] = ((continuum - (flux + flux_from_top_base_point_to_continuum)) / continuum)
        # Model: fit gaussian
        # Adjust edges
        new_base, new_top = improve_linemask_edges(smoothed_spectra, linemasks['base'][i], linemasks['top'][i], linemasks['peak'][i])
        linemasks['base_fit'][i] = new_base
        linemasks['top_fit'][i] = new_top
        linemasks['wave_base_fit'][i] = spectra['waveobs'][new_base]
        linemasks['wave_top_fit'][i] = spectra['waveobs'][new_top]
        try:
            #line_model = fit_line(spectra[base_points[i]:base_points[i+1]+1], continuum_model, linemasks['wave_peak'][i])
            line_model = fit_line(spectra[new_base:new_top+1], continuum_model, linemasks['wave_peak'][i])
            linemasks['mu'][i] = line_model.mu
            linemasks['sig'][i] = line_model.sig
            linemasks['A'][i] = line_model.A
            if type(line_model) == VoigtModel:
                linemasks['gamma'][i] = line_model.gamma
            else:
                # The model is Gaussian, do not use 'gamma'
                linemasks['gamma'][i] = 9999.0
            # Equivalent Width
            linemasks['flux'][i] = -1 * line_model.integrate(spectra['waveobs'][new_base], spectra['waveobs'][new_top])
            linemasks['ew'][i] = linemasks['flux'][i] / np.mean(continuum_model(spectra['waveobs'][new_base:new_top+1]))
            # RMS
            linemasks['rms'][i] = np.sqrt(np.sum(np.power(line_model.residuals(), 2)) / len(line_model.residuals()))
            ## Mu residual (sometimes RMS is very good but the fit is bad and it can be detected by
            ##              evaluating the flux at the mu position)
            ##            * Maximum value: 9999.0
            # - Modeled flux at mu
            modeled_mu_flux = line_model(line_model.mu)
            # - Observed flux at mu (interpolate if needed)
            observed_mu_flux = np.interp(line_model.mu, spectra[new_base:new_top+1] ['waveobs'], spectra[new_base:new_top+1] ['flux'] - continuum_model(spectra[new_base:new_top+1] ['waveobs']))
            linemasks['mu_residual'][i] = np.min([9999.0, np.abs(modeled_mu_flux - observed_mu_flux)])
        except Exception as e:
            # This should not happen, but just in case...
            linemasks['mu'][i] = 0.0
            linemasks['sig'][i] = 0.0
            linemasks['A'][i] = 0.0
            linemasks['gamma'][i] = 0.0
            linemasks['flux'][i] = 0.0
            linemasks['ew'][i] = 0.0
            linemasks['rms'][i] = 9999.0
            linemasks['mu_residual'][i] = 9999.0
            print "WARNING: Bad line fit (", i, ") - ", e.message
        
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
        # Fin duplicates
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

if __name__ == '__main__':
    ### VALD line list
    #VALD_to_SPECTRUM_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/linelists/VALD.300_1100nm.lst", depth_limit=0.0)

    #VALD_top_3_to_RV_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/rv/VALD.300_1100nm.rv.lst", top=1, wave_step=10)
    
    ### Line detection and fitting example:
    print "Reading spectrum..."
    spectra = read_spectra("input/L082N03_spec_norm/05oct08/sp2_Normal/004_vesta_001.s")
    resolution = 65000
    #spectra = read_spectra("input/test/observed_arcturus.s.gz")
    #spectra = read_spectra("input/test/observed_mu_cas_a.s.gz")
    #spectra = read_spectra("input/test/observed_mu_leo.s.gz")
    #spectra = read_spectra("input/test/observed_sun.s.gz")
    #resolution = 47000
    
    ## Test
    #spectra = read_spectra("/home/marble/Downloads/HD146233/Narval/narval_hd146233_100312_segment.s")
    #resolution = 65000
    #spectra = read_spectra("/home/marble/Downloads/HD146233/UVES/hd146233_uves_segment.txt")
    #spectra = read_spectra("/home/marble/Downloads/HD146233/UVES/hd146233_uves_segment_not_normalized.txt")
    #resolution = 47000
    ## Smooth spectra using the instrumental resolution
    smooth_spectra = True
    nknots_factor = 1
    
    # Telluric lines: constant continuum at 1, no filtering
    #spectra = read_spectra("input/telluric/standard_atm_air.s.gz")
    #resolution = 40000
    ## Smooth spectra using the instrumental resolution
    #smooth_spectra = True
    #nknots_factor = 3
    
    # Use method from Schneider et al (1993).
    #   The Hubble Space Telescope quasar absorption line key project. 
    #   II - Data calibration and absorption-line selection
    use_EW_SNR = False

    # Discard lines that do not have at least a given depth
    minimum_depth = 0.00 # (% of the continuum)
    maximum_depth = 1.00 # (% of the continuum)
    # Discard fitted lines with RMS higher than the median RMS
    discard_high_rms = False
    discard_bad_VALD_crossmatch = False
    discard_too_big_wavelength_range = True
    
    print "Convolving (smoothing)..."
    original_spectra = spectra
    if smooth_spectra:
        spectra = convolve_spectra(original_spectra, resolution)
    
    print "Finding peaks and base points..."
    if use_EW_SNR:
        peaks, base_points, corrected_convolved_ew_snr = find_peaks_and_base_points(spectra, use_EW_SNR=True)
    else:
        peaks, base_points = find_peaks_and_base_points(spectra, use_EW_SNR=False)
    
    # Fit continuum
    print "Fitting continuum..."
    # * 1 knot every 10 nm in average
    nknots = np.max([1, int((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 10)])
    nknots = nknots_factor * nknots # N knots every 10nm in average
    continuum_model = fit_continuum(spectra, nknots=nknots)


    print "Generating linemasks and fitting gaussians..."
    linemasks = generate_linemasks(original_spectra, peaks, base_points, continuum_model, smoothed_spectra=spectra ,vald_linelist_file="input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
    
    
    ##### Filters
    rejected_by_noise = detect_false_positives_and_noise(spectra, linemasks)

    if use_EW_SNR:
        # Filter low SNR peaks
        rejected_by_low_snr = (corrected_convolved_ew_snr[peaks] < np.median(corrected_convolved_ew_snr))
        rejected_by_noise = np.logical_or(rejected_by_noise, rejected_by_low_snr)
    
    if discard_bad_VALD_crossmatch:
        wave_diff_selected, wave_diff_selected_filter = sigma_clip(linemasks["wave_peak_diff"], sig=3, iters=1, varfunc=np.var,meanfunc=np.median) # Discard outliers
        rejected_by_bad_VALD_crossmatch = np.logical_not(wave_diff_selected_filter)
        #rejected_by_bad_VALD_crossmatch = np.abs(linemasks["wave_peak_diff"]) < np.mean(np.abs(linemasks["wave_peak_diff"]))
    else:
        rejected_by_bad_VALD_crossmatch = np.asarray([False]*len(linemasks))
    
    if discard_too_big_wavelength_range:
        # Reject linemasks with too big wavelength range (outliers)
        # WARNING: If the spectrum has a small wavelength range with some few strong features
        #          they may be erroneously discarded by this filter
        wave_diff = (linemasks['wave_top'] - linemasks['wave_base']) / (linemasks['top'] - linemasks['base'])
        wave_diff_selected, wave_diff_selected_filter = sigma_clip(wave_diff, sig=3, iters=1, varfunc=np.var,meanfunc=np.median) # Discard outliers
        rejected_by_wave_gaps = np.logical_not(wave_diff_selected_filter)
    else:
        rejected_by_wave_gaps = np.asarray([False]*len(linemasks))
    
   
    # Reject base points or peaks higher than continuum
    # - Depth is negative if the peak is higher than the continuum
    # - Relative depth is negative if the mean base point is higher than the continuum
    # DISCARDED: it is possible that a good line have base points a little bit higher than the continuum,
    # it is better to only consider the peak
    #rejected_by_depth_higher_than_continuum = np.logical_or((linemasks['relative_depth'] < 0), (linemasks['depth'] < 0))
    rejected_by_depth_higher_than_continuum = (linemasks['depth'] < 0)
    
    # Discard peaks with a depth inferior/superior to a given limit (% of the continuum)
    rejected_by_depth_limits = np.logical_or((linemasks['depth'] <= minimum_depth), (linemasks['depth'] >= maximum_depth))
    
    # Filter lines with a bad fit and high RMS (outliers)
    rejected_by_bad_mu_fit = np.logical_or((linemasks['mu'] < linemasks['wave_base']), (linemasks['mu'] > linemasks['wave_top']))
    
    rejected_by_bad_A_fit = linemasks["A"] > 0
    A_selected, A_selected_filter = sigma_clip(linemasks['A'], sig=3, iters=1, varfunc=np.var,meanfunc=np.median) # Discard outliers
    rejected_by_high_A = np.logical_or(rejected_by_bad_A_fit, np.logical_not(A_selected_filter))
    
    rejected_by_bad_sig_fit = linemasks["sig"] < 0
    sig_selected, sig_selected_filter = sigma_clip(linemasks['sig'], sig=3, iters=1, varfunc=np.var,meanfunc=np.median) # Discard outliers
    rejected_by_bad_sig_fit = np.logical_or(rejected_by_bad_sig_fit, np.logical_not(sig_selected_filter))
    
    # Do not reject gamma values of 9999.0, they indicate that the line is a Gaussian and not a Voigt
    rejected_by_bad_gamma_fit = linemasks["gamma"] != 9999.0
    gamma_selected, gamma_selected_filter = sigma_clip(linemasks['gamma'], sig=3, iters=1, varfunc=np.var,meanfunc=np.median) # Discard outliers
    rejected_by_bad_gamma_fit = np.logical_and(rejected_by_bad_gamma_fit, np.logical_not(gamma_selected_filter))
    # Discard negative gamma (strange wings that are not coherent with absoption line profiles)
    rejected_by_bad_gamma_fit = np.logical_or(rejected_by_bad_gamma_fit, linemasks["gamma"] < 0)
    
    if discard_high_rms:
        rms_selected, rms_selected_filter = sigma_clip(linemasks['rms'], sig=3, iters=1, varfunc=np.var,meanfunc=np.median) # Discard outliers
        rejected_by_high_rms = np.logical_not(rms_selected_filter)
        #rejected_by_high_rms = (linemasks['rms'] > np.mean(linemasks['rms']))
    else:
        rejected_by_high_rms = (linemasks['rms'] >= 9999.0) # 9999.0: Cases where has not been possible to fit a gaussian
    
    # Mu residual
    mu_residual_selected, mu_residual_selected_filter = sigma_clip(linemasks["mu_residual"], sig=3, iters=1, varfunc=np.var, meanfunc=np.median) # Discard outliers
    rejected_by_mu_residual = np.logical_not(mu_residual_selected_filter)
    rejected_by_mu_residual = np.logical_or(rejected_by_mu_residual, linemasks["mu_residual"] > 1)
    
    discarded = np.logical_or(rejected_by_noise, rejected_by_bad_VALD_crossmatch)
    discarded = np.logical_or(discarded, rejected_by_wave_gaps)
    discarded = np.logical_or(discarded, rejected_by_depth_higher_than_continuum)
    discarded = np.logical_or(discarded, rejected_by_depth_limits)
    discarded = np.logical_or(discarded, rejected_by_bad_mu_fit)
    discarded = np.logical_or(discarded, rejected_by_bad_A_fit)
    discarded = np.logical_or(discarded, rejected_by_bad_sig_fit)
    discarded = np.logical_or(discarded, rejected_by_bad_gamma_fit)
    discarded = np.logical_or(discarded, rejected_by_high_rms)
    discarded = np.logical_or(discarded, rejected_by_mu_residual)
    linemasks['discarded'][discarded] = True
    
    ## Peaks
    print "Number of peak candidates:", len(linemasks)
    print "- Noise and false positive:", len(linemasks[rejected_by_noise])
    print "- Bad VALD cross-match:", len(linemasks[rejected_by_bad_VALD_crossmatch])
    print "- Too big wavelength range (outliers):", len(linemasks[rejected_by_wave_gaps])
    print "- Peaks higher than continuum:", len(linemasks[rejected_by_depth_higher_than_continuum])
    print "- Out of depth limits:", len(linemasks[rejected_by_depth_limits])
    print "- Bad Gaussian fits - mu out of range:", len(linemasks[rejected_by_bad_mu_fit])
    print "- Bad Gaussian fits - A > 0 or outlier:", len(linemasks[rejected_by_bad_A_fit])
    print "- Bad Gaussian fits - sig < 0 or outlier:", len(linemasks[rejected_by_bad_sig_fit])
    print "- Bad Gaussian fits - gamma outlier:", len(linemasks[rejected_by_bad_gamma_fit])
    print "- Bad Gaussian fits - RMS outliers:", len(linemasks[rejected_by_high_rms])
    print "- Bad Gaussian fits - mu residual outliers:", len(linemasks[rejected_by_mu_residual])
    print "Final number of peaks:", len(linemasks[~discarded])
    print "- Voigt profile:", len(linemasks[~discarded][linemasks[~discarded]['gamma'] != 9999.0])
    print "- Gaussian profile:", len(linemasks[~discarded][linemasks[~discarded]['gamma'] == 9999.0])

    
    ##### Saving...
    import cPickle as pickle
    output_filename = "output/linemasks"
    # Complete
    asciitable.write(linemasks, output=output_filename+"_complete.txt", delimiter="\t")
    # Masks regions for visualization
    linemasks['note'][~discarded] = linemasks['element'][~discarded]
    linemasks['wave_peak'][~discarded] = linemasks['mu'][~discarded]
    linemasks['wave_base'][~discarded] = linemasks['wave_base_fit'][~discarded]
    linemasks['wave_top'][~discarded] = linemasks['wave_top_fit'][~discarded]
    asciitable.write(linemasks[~discarded], output=output_filename+".txt", delimiter="\t", include_names=['wave_peak', 'wave_base', 'wave_top', 'note'])
    # Continuum model
    pickle.dump((continuum_model.data, continuum_model.pardict), open(output_filename+"_continuum.dump", 'w'))
    ##### Restoring...
    import cPickle as pickle
    output_filename = "output/linemasks"
    linemasks = asciitable.read(output_filename+"_complete.txt", delimiter="\t")
    continuum_model = UniformCDFKnotSplineModel()
    continuum_model.data, continuum_model.pardict = pickle.load(open(output_filename+"_continuum.dump"))
    
    #print "Building a fitted spectrum..."
    waveobs = generate_wavelength_grid(np.min(spectra['waveobs']), np.max(spectra['waveobs']), resolution, points_per_fwhm = 3)
    fitted_spectra = build_fitted_spectrum(waveobs, continuum_model, linemasks[~discarded])
    
    # Plot
    used_peaks = peaks

    if use_EW_SNR:
        fig = plt.figure()
        axes = fig.add_subplot(211)
        plt.plot(spectra['waveobs'], spectra['flux'])
        plt.plot(spectra['waveobs'], continuum_model(spectra['waveobs']))
        plt.plot(fitted_spectra['waveobs'], fitted_spectra['flux'])
        plt.scatter(spectra['waveobs'][used_peaks], spectra['flux'][used_peaks], c='red')
        plt.scatter(spectra['waveobs'][base_points], spectra['flux'][base_points], c='green')
        axes = fig.add_subplot(212, sharex=axes)
        plt.plot(spectra['waveobs'], corrected_convolved_ew_snr) # green
        plt.show()
    else:
        fig = plt.figure()
        plt.plot(spectra['waveobs'], original_spectra['flux'])
        plt.plot(spectra['waveobs'], continuum_model(spectra['waveobs']))
        plt.plot(fitted_spectra['waveobs'], fitted_spectra['flux'])
        plt.scatter(spectra['waveobs'][used_peaks], original_spectra['flux'][used_peaks], c='red')
        plt.scatter(spectra['waveobs'][base_points], original_spectra['flux'][base_points], c='green')
        #for line in a:
            #a = plt.axvspan(line['wave_base'], line['wave_top'], facecolor='yellow', alpha=0.30)
            #l = plt.axvline(x = line['wave_peak'], linewidth=1, color='orange')
        plt.show()

    ## Resolution
    # Light speed in vacuum
    c = 299792458.0 # m/s
    sigma = linemasks['sig']
    fwhm = sigma * (2*np.sqrt(2*np.log(2))) # nm
    R = linemasks['wave_peak'] / fwhm
    
    # In m/s:
    #fwhm = c / R # m/s
    # ... or ...
    #fwhm = c * (fwhm / linemasks['wave_peak'][~discarded]) # m/s
    fwhm_kms = c / R # m/s
    fwhm_kms = fwhm_kms / 1000 # km/s

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
    #trend.fitData(linemasks['wave_peak'][~discarded][R[~discarded]<100000], R[~discarded][R[~discarded]<100000])
    #plt.scatter(linemasks['wave_peak'][~discarded][R[~discarded]<100000], R[~discarded][R[~discarded]<100000], s=4)
    #plt.plot(linemasks[~discarded]['wave_peak'], trend(linemasks[~discarded]['wave_peak']), color="red", linewidth=3)
    #plt.xlabel('Wavelength peak')
    #plt.ylabel('Resolving power')
    #plt.show()

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
