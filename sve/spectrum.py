#
#    This file is part of Spectra Visual Editor (SVE).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    SVE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SVE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with SVE. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
import pyfits
import scipy.ndimage as ndi
from spectrum import *
from plotting import *
from common import *
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt
import log
import logging

def __read_spectrum(spectrum_filename):
    try:
        spectrum = np.array([tuple(line.rstrip('\r\n').split("\t")) for line in open(spectrum_filename,)][1:], dtype=[('waveobs', float),('flux', float),('err', float)])
    except Exception as err:
        # Try without error column
        spectrum_tmp = np.array([tuple(line.rstrip('\r\n').split("\t")) for line in open(spectrum_filename,)][1:], dtype=[('waveobs', float),('flux', float)])
        spectrum = np.recarray((len(spectrum_tmp), ), dtype=[('waveobs', float),('flux', float),('err', float)])
        spectrum['waveobs'] = spectrum_tmp['waveobs']
        spectrum['flux'] = spectrum_tmp['flux']
        spectrum['err'] = 0.0
    return spectrum

def read_spectrum(spectrum_filename):
    """
    Return spectrum recarray structure from a filename.
    The file format shouldd be plain text files with **tab** character as column delimiter.
    Three columns should exists: wavelength, flux and error (although this last one is not a relevant value
    for the editor and it can be set all to zero).
    The first line should contain the header names 'waveobs', 'flux' and 'err' such as in the following example:
    ::

        waveobs       flux          err
        370.000000000 1.26095742505 1.53596736433
        370.001897436 1.22468868618 1.55692475754
        370.003794872 1.18323884263 1.47304952231
        370.005692308 1.16766911881 1.49393329036

    To save space, the file can be compressed in gzip format.

    If the specified file does not exists, it checks if there is a compressed version
    with the extension '.gz' (gzip) and if it exists, it will be automatically uncompressed.
    """
    # If it is not compressed
    if os.path.exists(spectrum_filename) and spectrum_filename[-3:] != ".gz":
        spectrum = __read_spectrum(spectrum_filename)
    elif (os.path.exists(spectrum_filename) and spectrum_filename[-3:] == ".gz") or (os.path.exists(spectrum_filename + ".gz")):
        if spectrum_filename[-3:] != ".gz":
            spectrum_filename = spectrum_filename + ".gz"

        tmp_spec = tempfile.mktemp() + str(int(random.random() * 100000000))
        # Uncompress to a temporary file
        f_out = open(tmp_spec, 'wb')
        f_in = gzip.open(spectrum_filename, 'rb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()

        spectrum = __read_spectrum(tmp_spec)
        os.remove(tmp_spec)

    # Filtering...
    valid = ~np.isnan(spectrum['flux'])

    # Find duplicate wavelengths
    dups, dups_index = find_duplicates(spectrum, 'waveobs')

    # Filter all duplicates except the first one
    last_wave = None
    for i in np.arange(len(dups)):
        if last_wave == None:
            last_wave = dups[i]['waveobs']
            continue
        if last_wave == dups[i]['waveobs']:
            pos = dups_index[i]
            valid[pos] = False
        else:
            # Do not filter the first duplicated value
            last_wave = dups[i]['waveobs']

    # Filter invalid and duplicated values
    spectrum = spectrum[valid]

    spectrum.sort(order='waveobs') # Make sure it is ordered by wavelength

    return spectrum

def write_spectrum(spectrum, spectrum_filename, compress=True):
    """
    Write spectrum to a file with the following file format:
    ::

        waveobs       flux          err
        370.000000000 1.26095742505 1.53596736433
        370.001897436 1.22468868618 1.55692475754
        370.003794872 1.18323884263 1.47304952231
        370.005692308 1.16766911881 1.49393329036

    """
    if compress:
        if spectrum_filename[-3:] != ".gz":
            spectrum_filename = spectrum_filename + ".gz"

        tmp_spec = tempfile.mktemp() + str(int(random.random() * 100000000))
        out = open(tmp_spec, "w")
        out.write("waveobs\tflux\terr\n")
        out.write("\n".join(["\t".join(map(str, (line['waveobs'], line['flux'], line['err']))) for line in spectrum]))
        out.close()

        # Compress the temporary file
        f_in = open(tmp_spec, 'rb')
        f_out = gzip.open(spectrum_filename, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        os.remove(tmp_spec)
    else:
        out = open(spectrum_filename, "w")
        out.write("waveobs\tflux\terr\n")
        out.write("\n".join(["\t".join(map(str, (line['waveobs'], line['flux'], line['err']))) for line in spectrum]))
        out.close()


def estimate_snr(flux, num_points=10, frame=None):
    """
    Estimate the Signal-to-Noise ratio for a given spectrum calculating the
    signal over standard deviation in blocks of N points and returning the average.
    """
    # Avoid negative values and outliers
    flux = flux[flux > 0.0]
    #flux, f = sigma_clipping(flux, sig=3, meanfunc=np.median)
    #flux, f = interquartile_range_filtering(flux, k=1.5)
    last_reported_progress = -1
    if num_points == 1:
        snr = np.mean(flux) / np.std(flux)
    else:
        snr = []
        total_num_blocks = len(flux)-num_points
        for i in np.arange(total_num_blocks):
            values = flux[i:i+num_points]
            stdev = np.std(values)
            if stdev != 0:
                snr.append(np.mean(values) / stdev)

            current_work_progress = ((i*1.0 / total_num_blocks) * 100.0)
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                logging.info("%.2f%%" % current_work_progress)
                if frame != None:
                    frame.update_progress(current_work_progress)
        snr = np.asarray(snr)
    #snr, s = sigma_clipping(snr, sig=3, meanfunc=np.median)
    snr, s = interquartile_range_filtering(snr, k=1.5)
    estimated_snr = np.mean(snr)
    logging.info("SNR = %.2f" % estimated_snr)
    return estimated_snr


def __get_fwhm(lambda_peak, from_resolution, to_resolution):
    """
    Calculate the FWHM of the gaussian needed to convert
    a spectrum from one resolution to another at a given wavelength point.
    """
    if from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are equal or bigger than original")
    from_delta_lambda = (1.0*lambda_peak) / from_resolution
    to_delta_lambda = (1.0*lambda_peak) / to_resolution
    fwhm = np.sqrt(to_delta_lambda**2 - from_delta_lambda**2)
    return fwhm

def __fwhm_to_sigma(fwhm):
    """
    Calculate the sigma value from the FWHM.
    """
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    return sigma


def convolve_spectrum(spectrum, from_resolution, to_resolution=None, frame=None):
    """
    Spectra resolution smoothness/degradation. Procedure:

    1) Define a bin per measure which marks the wavelength range that it covers.
    2) For each point, identify the window segment to convolve by using the bin widths and the FWHM.
    3) Build a gaussian using the sigma value and the wavelength values of the spectrum window.
    4) Convolve the spectrum window with the gaussian and save the convolved value.

    If "to_resolution" is not specified or its equal to "from_resolution", then the spectrum
    is convolved with the instrumental gaussian defined by "from_resolution".

    If "to_resolution" is specified, the convolution is made with the difference of
    both resolutions in order to degrade the spectrum.
    """
    if to_resolution != None and from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")

    total_points = len(spectrum)
    convolved_spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    convolved_spectrum['waveobs'] = spectrum['waveobs']
    convolved_spectrum['err'] = spectrum['err']

    last_reported_progress = -1
    if frame != None:
        frame.update_progress(0)

    flux = spectrum['flux']
    err = spectrum['err']
    # Consider the wavelength of the measurements as the center of the bins
    waveobs = spectrum['waveobs']
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

    # FWHM of the gaussian for the given resolution
    if to_resolution == None:
        # Convolve using instrumental resolution (smooth but not degrade)
        fwhm = waveobs / from_resolution
    else:
        # Degrade resolution
        fwhm = __get_fwhm(waveobs, from_resolution, to_resolution)
    sigma = __fwhm_to_sigma(fwhm)
    # Convert from wavelength units to bins
    fwhm_bin = fwhm / bin_width

    # Round number of bins per FWHM
    nbins = np.ceil(fwhm_bin) #npixels

    # Number of measures
    nwaveobs = len(waveobs)

    # In theory, len(nbins) == len(spectrum)
    for i in np.arange(len(nbins)):
        current_nbins = 2 * nbins[i] # Each side
        current_center = waveobs[i] # Center
        current_sigma = sigma[i]

        # Find lower and uper index for the gaussian, taking care of the current spectrum limits
        lower_pos = int(max(0, i - current_nbins))
        upper_pos = int(min(nwaveobs, i + current_nbins + 1))

        # Select only the flux values for the segment that we are going to convolve
        flux_segment = flux[lower_pos:upper_pos+1]
        err_segment = err[lower_pos:upper_pos+1]
        waveobs_segment = waveobs[lower_pos:upper_pos+1]

        nsegments = len(flux_segment)

        # Build the gaussian corresponding to the instrumental spread function
        gaussian = np.exp(- ((waveobs_segment - current_center)**2) / (2*current_sigma**2)) / np.sqrt(2*np.pi*current_sigma**2)
        gaussian = gaussian / np.sum(gaussian)

        # Convolve the current position by using the segment and the gaussian
        if flux[i] > 0:
            weighted_flux = flux_segment * gaussian
            current_convolved_flux = weighted_flux.sum()
            convolved_spectrum['flux'][i] = current_convolved_flux
        else:
            convolved_spectrum['flux'][i] = 0.0

        if err[i] > 0:
            # * Propagate error Only if the current value has a valid error value assigned
            #
            # Error propagation considering that measures are dependent (more conservative approach)
            # because it is common to find spectra with errors calculated from a SNR which
            # at the same time has been estimated from all the measurements in the same spectra
            #
            weighted_err = err_segment * gaussian
            current_convolved_err = weighted_err.sum()
            #current_convolved_err = np.sqrt(np.power(weighted_err, 2).sum()) # Case for independent errors
            convolved_spectrum['err'][i] = current_convolved_err
        else:
            convolved_spectrum['err'][i] = 0.0

        current_work_progress = (i*1.0 / total_points) * 100
        if report_progress(current_work_progress, last_reported_progress):
            last_reported_progress = current_work_progress
            logging.info("%.2f%%" % current_work_progress)
            if frame != None:
                frame.update_progress(current_work_progress)
    logging.info("Spectra convolved!")

    return convolved_spectrum


def __interpolate_flux(spectrum, wavelength):
    """
    Interpolate flux for a given wavelength by using Bessel's Central-Difference Interpolation.
    It considers:

    - 4 points in general
    - 2 when there are not more (i.e. at the beginning of the array or outside)
    """
    # Target wavelength
    objective_wavelength = wavelength
    fluxes = spectrum['flux']
    waveobs = spectrum['waveobs']

    # Find the index position of the first wave length equal or higher than the objective
#    index = np.where(waveobs >= objective_wavelength)[0][0]
    index = waveobs.searchsorted(objective_wavelength)

    total_points = len(spectrum)
    if index == total_points:
        # DISCARD: Linear extrapolation using index-1 and index-2
        # flux = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index-1]-fluxes[index-2])/(waveobs[index-1]-waveobs[index-2]))
        # JUST DUPLICATE:
        flux = fluxes[index-1]
    elif index == 1 or index == total_points-1:
        # Linear interpolation between index and index-1
        # http://en.wikipedia.org/wiki/Linear_interpolation#Linear_interpolation_between_two_known_points
        flux = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index]-fluxes[index-1])/(waveobs[index]-waveobs[index-1]))
    elif index == 0 and waveobs[index] != objective_wavelength:
        # DISCARD: Linear extrapolation using index+1 and index
        # flux = fluxes[index] + (objective_wavelength - waveobs[index]) * ((fluxes[index+1]-fluxes[index])/(waveobs[index+1]-waveobs[index]))
        # JUST DUPLICATE:
        flux = fluxes[index]
    elif waveobs[index] == objective_wavelength:
        flux = fluxes[index]
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

        p = (objective_wavelength - wave_x0) / (wave_x1 - wave_x0)
        flux = flux_x0 + p * (flux_x1 - flux_x0) + (p * (p - 1) / 4) * (flux_x2 - flux_x1 - flux_x0 + flux_x_1)


#    print flux, fluxes[index], wavelength
    return flux, index


def resample_spectrum(spectrum, xaxis, linear=True, frame=None):
    """
    Returns a new spectrum with measures at the given xaxis wavelength
    Interpolation is completely linear by default (fastest option) but a Bessel's
    Central-Difference Interpolation with 4 points can be activated by
    specifying "linear=False", in that case interpolation is linear only
    when there are not enough points (i.e. beginning/end of spectrum).
    """
    total_points = len(xaxis)
    last_reported_progress = -1

    if linear:
        current_work_progress = 10.0
        logging.info("%.2f%%" % current_work_progress)
        if frame != None:
            frame.update_progress(current_work_progress)
        resampled_spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        resampled_spectrum['waveobs'] = xaxis
        resampled_spectrum['flux'] = np.interp(xaxis, spectrum['waveobs'], spectrum['flux'], left=0.0, right=0.0) # No extrapolation, just returns zeros

        current_work_progress = 90.0
        logging.info("%.2f%%" % current_work_progress)
        if frame != None:
            frame.update_progress(current_work_progress)
    else:
        resampled_spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        resampled_spectrum['waveobs'] = xaxis

        from_index = 0 # Optimization: discard regions already processed
        for i in np.arange(total_points):
            resampled_spectrum['flux'][i], index = __interpolate_flux(spectrum[from_index:], resampled_spectrum['waveobs'][i])
            if index > 4:
                from_index = index - 4
            current_work_progress = np.min([(i*1.0 / total_points) * 100, 90.0])
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                logging.info("%.2f%%" % current_work_progress)
                if frame != None:
                    frame.update_progress(current_work_progress)

    resampled_spectrum['err'] = np.interp(xaxis, spectrum['waveobs'], spectrum['err'])
    return resampled_spectrum

def correct_velocity(spectrum, velocity):
    """
    Correct velocity in km/s.
    """
    # Speed of light in m/s
    c = 299792458.0
    # Radial/barycentric velocity from km/s to m/s
    velocity = velocity * 1000

    # Correct wavelength scale for radial velocity
    spectrum['waveobs'] = spectrum['waveobs'] / ((velocity / c) + 1)
    return spectrum

def correct_velocity_regions(regions, velocity, with_peak=False):
    """
    Correct regions' velocity in km/s.
    """
    # Speed of light in m/s
    c = 299792458.0
    # Radial/barycentric velocity from km/s to m/s
    velocity = velocity * 1000

    regions['wave_base'] = regions['wave_base'] / ((velocity / c) + 1)
    regions['wave_top'] = regions['wave_top'] / ((velocity / c) + 1)
    if with_peak:
        regions['wave_peak'] = regions['wave_peak'] / ((velocity / c) + 1)
    return regions


