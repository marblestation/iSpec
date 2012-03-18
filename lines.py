import asciitable
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions
import ipdb
from common import *
from continuum import *
from fitting import *
from radial_velocity import *
from convolve import *
import matplotlib.pyplot as plt
from pymodelfit import UniformKnotSplineModel

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
        except ValueError:
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

# For an array of values, determine if a value is the maximum of a window
# of "span" elements
def find_max_win(x, span=3):
    ret = []
    n = len(x)
    dist = (span + 1) / 2;
    m = 0;
    for i in np.arange(n):
        l_min = np.max([i-dist+1, 0])
        l_max = i-1
        r_min = i+1
        r_max = np.min([i+dist-1, n-1])
        is_max = 1;
        # left side
        j = l_min
        while j <= l_max:
            if (x[j] > x[i]):
                is_max = 0;
                break
            j += 1
        
        # right side
        if (is_max == 1):
            j = r_min
            while j <= r_max:
                if (x[j] > x[i]):
                    is_max = 0;
                    break
                j += 1
        if (is_max == 1):
            ret.append(i)
    return np.asarray(ret)

# For an array of values, determine if a value is the minimum of a window
# of "span" elements
def find_min_win(x, span=3):
    ret = []
    n = len(x)
    dist = (span + 1) / 2;
    m = 0;
    for i in np.arange(n):
        l_min = np.max([i-dist+1, 0])
        l_max = i-1
        r_min = i+1
        r_max = np.min([i+dist-1, n-1])
        is_min = 1;
        # left side
        j = l_min
        while j <= l_max:
            if (x[j] < x[i]):
                is_min = 0;
                break
            j += 1
        
        # right side
        if (is_min == 1):
            j = r_min
            while j <= r_max:
                if (x[j] < x[i]):
                    is_min = 0;
                    break
                j += 1
        if (is_min == 1):
            ret.append(i)
    return np.asarray(ret)



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

# Remove consecutive features (i.e. peaks or base points)
def remove_consecutives_features(features):
    duplicated_features = (np.abs(features[1:-1] - features[2:]) == 1)
    duplicated_features = np.array([False] + duplicated_features.tolist() + [False])
    cleaned_features = features[~duplicated_features]
    return cleaned_features


def detect_false_positives_and_noise(spectra, peaks, base_points):
    # WARNING: First and last peak will be always ignored in order to simplify the algorithm
    #          but in any case, it is not a big deal
    # Find
    # - Peaks that are less deep than it nearby base points (false positives)
    # - Peaks too close (2 or less positions) to the next and previous base point (noise)
    first_wave_peak = spectra['waveobs'][peaks][0]
    first_wave_base = spectra['waveobs'][base_points][0]
    if first_wave_peak > first_wave_base:
        if len(base_points) - len(peaks) == 1:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: base point
            # Left
            peak_base_diff_left = spectra['flux'][base_points[:-1]] - spectra['flux'][peaks]
            peak_base_index_diff_left = base_points[:-1] - peaks
            # Right
            peak_base_diff_right = spectra['flux'][base_points[1:]] - spectra['flux'][peaks]
            peak_base_index_diff_right = base_points[1:] - peaks
            # Ignore first and last peak
            peak_base_diff_left = peak_base_diff_left[1:-1]
            peak_base_index_diff_left = peak_base_index_diff_left[1:-1] 
            peak_base_diff_right = peak_base_diff_right[1:-1]
            peak_base_index_diff_right = peak_base_index_diff_right[1:-1]
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: peak (this last one will be ignored)
            # Left
            peak_base_diff_left = spectra['flux'][base_points[:-1]] - spectra['flux'][peaks[:-1]]
            peak_base_index_diff_left = base_points[:-1] - peaks[:-1]
            # Right
            peak_base_diff_right = spectra['flux'][base_points[1:]] - spectra['flux'][peaks[:-1]]
            peak_base_index_diff_right = base_points[1:] - peaks[:-1]
            # Ignore also first peak
            peak_base_diff_left = peak_base_diff_left[1:]
            peak_base_index_diff_left = peak_base_index_diff_left[1:] 
            peak_base_diff_right = peak_base_diff_right[1:]
            peak_base_index_diff_right = peak_base_index_diff_right[1:]
        else:
            raise Exception("This should not happen")
    else:
        if len(base_points) - len(peaks) == -1:
            ## First feature found in spectra: peak (this first one will be ignored)
            ## Last feature found in spectra: peak (this last one will be ignored)
            # Left
            peak_base_diff_left = spectra['flux'][base_points[:-1]] - spectra['flux'][peaks[1:-1]]
            peak_base_index_diff_left = base_points[:-1] - peaks[1:-1]
            peak_base_diff_right = spectra['flux'][base_points[1:]] - spectra['flux'][peaks[1:-1]]
            peak_base_index_diff_right = base_points[1:] - peaks[1:-1]
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: peak (this first one will be ignored)
            ## Last feature found in spectra: base point
            # Left
            peak_base_diff_left = spectra['flux'][base_points[:-1]] - spectra['flux'][peaks[1:]]
            peak_base_index_diff_left = base_points[:-1] - peaks[1:]
            # Right
            peak_base_diff_right = spectra['flux'][base_points[1:]] - spectra['flux'][peaks[1:]]
            peak_base_index_diff_right = base_points[1:] - peaks[1:]
            # Ignore also last peak
            peak_base_diff_left = peak_base_diff_left[:-1]
            peak_base_index_diff_left = peak_base_index_diff_left[:-1] 
            peak_base_diff_right = peak_base_diff_right[:-1]
            peak_base_index_diff_right = peak_base_index_diff_right[:-1]
        else:
            raise Exception("This should not happen")
    
    # Find peaks too close (2 positions) to the next or previous peak (noise)
    # First and last peaks are ignored
    peak_peak_index_diff_left = peaks[1:-1] - peaks[2:]
    peak_peak_index_diff_right = peaks[1:-1] - peaks[:-2]
    
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

# Determine lines that have a depth less than a tolerance level
# Depth is calculated by estimating a fixed continuum value
def reject_by_depth(spectra, peaks, base_points, minimum_depth, continuum_model=None):
    # WARNING: First and last peak will be always ignored in order to simplify the algorithm
    #          but in any case, it is not a big deal
    if continuum_model==None:
        continuum = np.median(spectra['flux'][base_points])
    else:
        continuum = continuum_model(spectra['waveobs'][peaks[1:-1]])
    flux = spectra['flux'][peaks[1:-1]]
    depth = (continuum - flux) / continuum
    rejected_by_depth = (depth < minimum_depth)
    return rejected_by_depth

def generate_linemasks_and_models(spectra, peaks, base_points, discarded, continuum_model):
    # Limit the base_points array to the ones that are useful, considering that
    # the first and last peak are always removed
    first_wave_peak = spectra['waveobs'][peaks][0]
    first_wave_base = spectra['waveobs'][base_points][0]
    if first_wave_peak > first_wave_base:
        if len(base_points) - len(peaks) == 1:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: base point
            base_points = base_points[1:-1]
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: base point
            ## Last feature found in spectra: peak (this last one will be ignored)
            base_points = base_points[1:]
        else:
            raise Exception("This should not happen")
    else:
        if len(base_points) - len(peaks) == -1:
            ## First feature found in spectra: peak (this first one will be ignored)
            ## Last feature found in spectra: peak (this last one will be ignored)
            base_points = base_points
        elif len(base_points) - len(peaks) == 0:
            ## First feature found in spectra: peak (this first one will be ignored)
            ## Last feature found in spectra: base point
            base_points = base_points[:-1]
        else:
            raise Exception("This should not happen")
    total_regions = len(peaks[1:-1][~discarded])
    lines = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100'), ('mu', float), ('sig', float), ('A', float), ('depth', float), ('rms', float)])
    
    if continuum_model == None:
        continuum = np.median(spectra['flux'][base_points])
    
    num_models = len(peaks[1:-1][~discarded])
    num_peaks = len(peaks[1:-1])
    models = np.empty(num_peaks, dtype=object) # As big as "discarded"
    model_rms = np.ones(num_peaks) # As big as "discarded"
    model_rms *= -9999 # If no model is fitted, the RMS will be "infinite" (represented by -9999)
    j = 0
    for i in np.arange(num_peaks):
        if not discarded[i]:
            # Line mask
            lines['wave_peak'][j] = spectra['waveobs'][peaks[1:-1][i]]
            lines['wave_base'][j] = spectra['waveobs'][base_points[i]]
            lines['wave_top'][j] = spectra['waveobs'][base_points[i+1]]
            flux = spectra['flux'][peaks[1:-1]][i]
            continuum = continuum_model(lines['wave_peak'][j])
            lines['note'][j] = ""
            lines['depth'][j] = ((continuum - flux) / continuum)
            ## Model: fit gaussian
            # Initial parameters
            # pymodelfit uses (self.mu-abs(sig)*4,self.mu+abs(sig)*4) as initial limits to search for the peak
            # so we should specify a coherent sigma for the wavelength range of the line mask
            max_separation = np.max([lines['wave_peak'][j] - lines['wave_base'][j], lines['wave_top'][j] - lines['wave_peak'][j]])
            sig = max_separation / 4
            line_model, continuum_value, lineflux, ew = fit_line(spectra[base_points[i]:base_points[i+1]], lines['wave_peak'][j], sig=sig, continuum_model=continuum_model)
            lines['mu'][j] = line_model.mu
            lines['sig'][j] = line_model.sig
            lines['A'][j] = line_model.A
            lines['rms'][j] = np.sqrt(np.sum(np.power(line_model.residuals(), 2)) / len(line_model.residuals()))
            # Only accept the fit if the peak (mu) is in the delimited region
            if line_model.mu >= lines['wave_base'][j] and line_model.mu <= lines['wave_top'][j]:
                models[i] = (line_model)
                model_rms[i] = lines['rms'][j]
            if (j % 100) == 0:
                print "%.2f%%" % (((i*1.0)/num_peaks) * 100)
            j += 1
    return lines, models, model_rms

def find_peaks_and_base_points(spectra, use_EW_SNR=False, resolution=None):
    # Determine peaks and base points (also known as continuum points)
    if use_EW_SNR:
        if resolution == None:
            raise Exception("Resolution is needed")
        # - SNR: Schneider et al (1993). The Hubble Space Telescope quasar absorption line key project. II - Data calibration and absorption-line selection
        # WARNING: Span should be 3 in order to alway have the same number of peaks and base_points +/-1
        #          which is a necessary condition for the rest of the algorithm
        corrected_convolved_ew_snr = get_corrected_convolved_ew_snr(spectra, resolution)
        peaks = find_max_win(corrected_convolved_ew_snr, span=3)
        base_points = find_min_win(corrected_convolved_ew_snr, span=3)
    else:
        ## - Convolved spectrum
        # WARNING: Span should be 3 in order to alway have the same number of peaks and base_points +/-1
        #          which is a necessary condition for the rest of the algorithm
        peaks = find_min_win(spectra['flux'], span=3)
        base_points = find_max_win(spectra['flux'], span=3)

    # WARNING: Due to three or more consecutive values with exactly the same flux
    # find_max_win or find_min_win will identify all of them as peaks or bases,
    # where only one of the should be marked as peak or base.
    # These cases break the necessary condition of having the same number of 
    # peaks and base_points +/-1
    # It is necessary to find those "duplicates" and remove them:
    peaks = remove_consecutives_features(peaks)
    base_points = remove_consecutives_features(base_points)

    if not (len(peaks) - len(base_points)) in [-1, 0, 1]:
        raise Exception("This should not happen")
    
    return peaks, base_points


# Determine max points by using a moving window of 3 elements (also used for line determination)
# Remove outliers and fits a uniform knot spline model with n knots
# Returns the model
def fit_continuum(spectra, nknots=None):
    points_per_bin = 20
    # Find max points in windows of 3 measures
    base_points = find_max_win(spectra['flux'], span=3)
    # Find outliers considering 3 sigma around the median
    flux_selected, filter_outliers = sigma_clip(spectra['flux'][base_points], sig=3, iters=10)
    base_points = base_points[filter_outliers]
    # If there are enough points
    num_base_points = len(base_points)
    if num_base_points > points_per_bin:
        # Group points in bins and use only the one with the higher flux
        num_final_base_points = int(num_base_points/points_per_bin)
        final_base_points = np.empty(num_final_base_points, dtype=int)
        for i in np.arange(num_final_base_points):
            ini = i*points_per_bin
            end = np.min([(i+1)*points_per_bin, num_base_points])
            # Find max for this bin
            max_flux = np.max(spectra['flux'][base_points[ini:end]])
            index = np.where(spectra['flux'][base_points[ini:end]] == max_flux)
            final_base_points[i] = base_points[ini:end][index]
    else:
        num_final_base_points = num_base_points
        final_base_points = base_points
    if nknots == None:
        # * 1 knot every 10 nm
        nknots = np.max([1, int((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 10)])
    # Fit a uniformly divided splines (n knots)
    continuum_model = UniformKnotSplineModel(nknots=nknots)
    #continuum_model.fitData(spectra['waveobs'][base_points][filter_outliers], spectra['flux'][base_points][filter_outliers])
    continuum_model.fitData(spectra['waveobs'][final_base_points], spectra['flux'][final_base_points])
    return continuum_model

def build_fitted_spectrum(waveobs, line_models):
    # Build a fitted spectrum
    fluxes = np.zeros(len(waveobs))
    num_line_models = len(line_models)
    i = 0
    for line_model in line_models:
        if (i % 100) == 0:
            print "%.2f%%" % (((i*1.0)/num_line_models) * 100)
        line_flux = line_model(waveobs)
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


#if __name__ == '__main__':
    #pass
    ### VALD lines
    ##VALD_to_SPECTRUM_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/linelists/VALD.300_1100nm.lst", depth_limit=0.0)

    ##VALD_top_3_to_RV_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/rv/VALD.300_1100nm.rv.lst", top=1, wave_step=10)
    
##############################################

print "Reading spectrum..."
#spectra = read_spectra("input/L082N03_spec_norm/05oct08/sp2_Normal/004_vesta_001.s")
spectra = read_spectra("input/test/observed_arcturus.s.gz")
#spectra = read_spectra("input/test/observed_mu_cas_a.s.gz")
#spectra = read_spectra("input/test/observed_mu_leo.s.gz")
#spectra = read_spectra("input/test/observed_sun.s.gz")
resolution = 65000
minimum_depth = None # (% of the continuum)
minimum_depth = 0.01

print "Convolving..."
# Smooth spectra using the instrumental resolution
convolved_spectra = convolve_spectra(spectra, resolution)

print "Finding peaks and base points..."
use_EW_SNR = False
peaks, base_points = find_peaks_and_base_points(convolved_spectra, use_EW_SNR=False)

# Fit continuum
print "Fitting continuum..."
continuum_model = fit_continuum(convolved_spectra)

rejected_by_noise = detect_false_positives_and_noise(convolved_spectra, peaks, base_points)

if use_EW_SNR:
    # Filter low SNR peaks
    rejected_by_low_snr = (corrected_convolved_ew_snr[peaks[1:-1]] < np.median(corrected_convolved_ew_snr))
    rejected = np.logical_or(rejected_by_noise, rejected_by_low_snr)

# Discard peaks with a depth inferior to a given limit (% of the continuum)
# - The higher the tolerance, the more rejected points
if minimum_depth != None:
    rejected_by_depth = reject_by_depth(convolved_spectra, peaks, base_points, minimum_depth, continuum_model=continuum_model)
    rejected = np.logical_or(rejected_by_noise, rejected_by_depth)
else:
    rejected = rejected_by_noise

print "Generating linemasks and fitting gaussians..."
linemasks, line_models, model_rms = generate_linemasks_and_models(spectra, peaks, base_points, rejected, continuum_model)
#asciitable.write(lines, output="x.txt", delimiter="\t")

# Filter lines with a bad fit/RMS
none_rms = (model_rms == -9999) # lines without model
median_rms = np.median(np.abs(model_rms[~none_rms]))
rejected_by_badfit = np.logical_or((model_rms > median_rms), (model_rms < -1*median_rms))
rejected_by_badfit = np.logical_or(rejected_by_badfit, none_rms)
discarded = np.logical_or(rejected, rejected_by_badfit)

print "Number of peak candidates:", len(peaks[1:-1])
print "After filtering noise and false positive:", len(peaks[1:-1][~rejected_by_noise])
print "After tolerance:", len(peaks[1:-1][~rejected])
print "After Gaussian fits:", len(peaks[1:-1][~discarded])

print "Building a fitted spectrum..."
waveobs = generate_wavelength_grid(np.min(spectra['waveobs']), np.max(spectra['waveobs']), resolution, points_per_fwhm = 3)
fitted_spectra = build_fitted_spectrum(waveobs, line_models[~discarded])

# Plot
final_peaks = peaks[1:-1][~discarded]

if use_EW_SNR:
    fig = plt.figure()
    axes = fig.add_subplot(211)
    plt.plot(spectra['waveobs'], spectra['flux'])
    plt.scatter(spectra['waveobs'][final_peaks], spectra['flux'][final_peaks], c='red')
    plt.scatter(spectra['waveobs'][base_points], spectra['flux'][base_points], c='green')
    axes = fig.add_subplot(212, sharex=axes)
    plt.plot(spectra['waveobs'], corrected_convolved_ew_snr) # green
    plt.show()
else:
    fig = plt.figure()
    plt.plot(spectra['waveobs'], spectra['flux'])
    plt.plot(spectra['waveobs'], continuum_model(spectra['waveobs']))
    plt.plot(fitted_spectra['waveobs'], fitted_spectra['flux'])
    plt.scatter(spectra['waveobs'][final_peaks], spectra['flux'][final_peaks], c='red')
    plt.scatter(spectra['waveobs'][base_points], spectra['flux'][base_points], c='green')
    plt.show()
