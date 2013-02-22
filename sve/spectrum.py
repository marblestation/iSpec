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
import numpy.lib.recfunctions as rfn # Extra functions
import pyfits
import scipy.ndimage as ndi
from spectrum import *
from plotting import *
from common import *
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import log
import logging

def __read_fits_spectrum(spectrum_filename, fluxhdu="PRIMARY", errorhdu=None):
    """
    Reads the 'PRIMARY' HDU of the FITS file, considering that it contains the fluxes.

    The wavelength are derived from the headers, if not possible it checks if the
    data from the HDU contains 2 axes and takes the first as the wavelength
    and the second as the flux.

    It tries to find the errors in other HDU by checking the names and the length,
    if none are found then they are set to zero.

    Inspired by pyspeckit:
        https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/src/ae1e0714410b58905466740b04b54318d5f318f8/pyspeckit/spectrum/readers/fits_reader.py?at=default
    """
    hdulist = pyfits.open(spectrum_filename)

    data = hdulist[fluxhdu].data
    hdr = hdulist[fluxhdu].header

    axis = 1
    specaxis = str(axis)
    # Try to determine if wavelength is in Angstrom (by default) or nm
    ctype = hdr.get('CTYPE%i' % axis)
    cunit = hdr.get('CUNIT%i' % axis)
    if str(cunit).upper() in ['NM']:
        unit = "nm"
    else:
        #if str(ctype).upper() in ['AWAV', 'ANGSTROM', '0.1 NM'] or str(cunit).upper() in ['AWAV', 'ANGSTROM', '0.1 NM']:
        unit = "Angstrom"

    flux = data.flatten()
    waveobs = None
    # NOTE: One can use hdr.get('ORIGIN') for special treatments
    if hdr.get(str('CD%s_%s' % (specaxis,specaxis))) != None:
        wave_step = hdr['CD%s_%s' % (specaxis,specaxis)]
        wave_base = hdr['CRVAL%s' % (specaxis)]
        reference_pixel = hdr['CRPIX%s' % (specaxis)]
        logging.info("Using the FITS CD matrix.  PIX=%f VAL=%f DELT=%f UNIT=%s" % (reference_pixel,wave_base,wave_step,unit))
    elif hdr.get(str('CDELT%s' % (specaxis))) != None:
        wave_step = hdr['CDELT%s' % (specaxis)]
        wave_base = hdr['CRVAL%s' % (specaxis)]
        reference_pixel = hdr['CRPIX%s' % (specaxis)]
        logging.info("Using the FITS CDELT value.  PIX=%f VAL=%f DELT=%f UNIT=%s" % (reference_pixel,wave_base,wave_step,unit))
    elif len(data.shape) > 1:
        logging.info("No CDELT or CD in header.  Assuming 2D input with 1st line representing the spectral axis.")
        # try assuming first axis is X axis
        if hdr.get('CUNIT%s' % (specaxis)) != None:
            waveobs = data[0,:]
            flux = data[1,:]
            if data.shape[0] > 2:
                errspec = data[2,:]
        else:
            raise Exception("Unknown FITS file format")
    else:
        raise Exception("Unknown FITS file format")

    # Angstrom to nm
    if unit != "nm":
        wave_base /= 10
        wave_step /= 10

    # Deal with logarithmic wavelength binning if necessary
    if waveobs is None:
        if hdr.get('WFITTYPE') == 'LOG-LINEAR':
            xconv = lambda v: 10**((v-reference_pixel+1)*wave_step+wave_base)
            waveobs = xconv(np.arange(len(flux)))
        else:
            xconv = lambda v: ((v-reference_pixel+1)*wave_step+wave_base)
            waveobs = xconv(np.arange(len(flux)))


    num_measures = len(flux)
    spectrum = np.recarray((num_measures, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    spectrum['waveobs'] = waveobs
    spectrum['flux'] = flux

    if errorhdu == None:
        spectrum['err'] = np.zeros(len(flux))
        # Try to find the errors in the extensions (HDU different than the PRIMARY):
        for i in xrange(len(hdulist)):
            name = hdulist[i].name.upper()
            if name == str(fluxhdu) or len(hdulist[i].data.flatten()) != len(flux):
                continue
            if 'IVAR' in name or 'VARIANCE' in name:
                spectrum['err'] = 1. / hdulist[i].data.flatten()
                break
            if 'NOISE' in name or 'ERR' in name or 'SIGMA' in name:
                spectrum['err'] = hdulist[i].data.flatten()
                break
    else:
        spectrum['err'] = hdulist[errorhdu].data.flatten()

    hdulist.close()

    return spectrum

def __read_spectrum(spectrum_filename):
    try:
        spectrum = np.array([tuple(line.rstrip('\r\n').split("\t")) for line in open(spectrum_filename,)][1:], dtype=[('waveobs', float),('flux', float),('err', float)])
        if len(spectrum) == 0:
            raise Exception("empty spectrum or incompatible format")
    except Exception as err:
        # try narval plain text format:
        # - ignores 2 first lines (header)
        # - ignores last line (empty)
        # - lines separated by \r
        # - columns separated by space
        narval = open(spectrum_filename,).readlines()[0].split('\r')
        spectrum = np.array([tuple(line.rstrip('\r').split()) for line in narval[2:-1]], dtype=[('waveobs', float),('flux', float),('err', float)])
    return spectrum

def read_spectrum(spectrum_filename, fits_options={"fluxhdu": "PRIMARY", "errorhdu": None}):
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

    ** It can recognise FITS files by the filename (extensions .FITS or .FIT), if this is
    the case, then it tries to load the PRIMARY spectra by default and tries to search the errors
    in the extensions of the FITS file (this behaviour can be modified by specifying
    the HDU for the flux and the errors via "fits_options").
    """
    # If it is not compressed
    if os.path.exists(spectrum_filename) and (spectrum_filename[-4:].lower() == ".fit" or spectrum_filename[-5:].lower() == ".fits") :
        # Make sure that "fluxhdu" and "errorhdu" exist:
        default = {"fluxhdu": "PRIMARY", "errorhdu": None}
        default.update(fits_options) # Overwrite with user input values
        fits_options = default
        spectrum = __read_fits_spectrum(spectrum_filename, fluxhdu = fits_options["fluxhdu"], errorhdu = fits_options["errorhdu"])
    elif os.path.exists(spectrum_filename) and spectrum_filename[-3:].lower() != ".gz":
        spectrum = __read_spectrum(spectrum_filename)
    elif (os.path.exists(spectrum_filename) and spectrum_filename[-3:].lower() == ".gz") or (os.path.exists(spectrum_filename.lower() + ".gz")):
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
    else:
        raise Exception("Spectrum file does not exists!")

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

def write_spectrum(spectrum, spectrum_filename):
    """
    Write spectrum to a file with the following file format:
    ::

        waveobs       flux          err
        370.000000000 1.26095742505 1.53596736433
        370.001897436 1.22468868618 1.55692475754
        370.003794872 1.18323884263 1.47304952231
        370.005692308 1.16766911881 1.49393329036

    ** If the filename has the extension ".FIT" or ".FITS", then the file is saved
    in FITS format. If the spectrum is not regularly sampled, then it will save
    the flux and the wavelengths as a matrix in the primary HDU. If the errors
    are different from zero, they will be saved as an extension.
    """
    if spectrum_filename[-4:].lower() == ".fit" or spectrum_filename[-5:].lower() == ".fits":
        wave_diff = spectrum['waveobs'][1:] - spectrum['waveobs'][:-1]
        if np.all(np.abs(wave_diff - np.median(wave_diff)) < 0.0000001):
            # Regularly sampled spectrum
            data = spectrum['flux']
            header = pyfits.Header()
            header.update('TTYPE1', "FLUX")
            header.update('TUNIT1', "COUNTS")
            header.update('CD1_1', spectrum['waveobs'][1] - spectrum['waveobs'][0])
            header.update('CRVAL1', spectrum['waveobs'][0])
            header.update('CRPIX1', 1)

            header.update('NAXIS', 1)
            header.update('NAXIS1', len(spectrum['flux']))
        else:
            # Since it is not regularly sampled, add the wavelength and flux together
            # in a matrix
            data = np.vstack([spectrum['waveobs'], spectrum['flux']])
            header = pyfits.Header()
            header.update('TTYPE1', "WAVELENGTH")
            header.update('TUNIT1', "NM")
            header.update('TTYPE2', "FLUX")
            header.update('TUNIT2', "COUNTS")
            header.update('NAXIS', 2)
            header.update('NAXIS1', len(spectrum['waveobs']))
            header.update('NAXIS2', 2) # waveobs and flux
        header.update('CUNIT1', "NM")
        header.update('CTYPE1', "WAVELENGTH")
        header.update('ORIGIN', "SVE")
        #header.update('VERSION', "SVE")
        header.update('UTCSAVED', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        primary_hdu = pyfits.PrimaryHDU(data=data, header=header)

        # Add an HDU extension with errors if they exist
        if np.any(spectrum['err'] != 0):
            # Error extension
            data = spectrum['err']
            header = pyfits.Header()
            header.update('TTYPE1', "FLUX")
            header.update('TUNIT1', "COUNTS")
            header.update('NAXIS', 1)
            header.update('NAXIS1', len(spectrum['err']))

            sigma_hdu = pyfits.ImageHDU(data=data, header=header, name="SIGMA")
            fits = pyfits.HDUList([primary_hdu, sigma_hdu])
        else:
            fits = pyfits.HDUList(primary_hdu)

        fits.writeto(spectrum_filename, clobber=True)
    elif spectrum_filename[-3:].lower() == ".gz":
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


def resample_spectrum(spectrum, xaxis, method="linear", frame=None):
    """
    Returns a new spectrum with measures at the given xaxis wavelength
    Interpolation is completely linear by default (fastest option) but
    it can also be:

    * method = "bessel": A Bessel's Central-Difference Interpolation with
    4 points. In this case interpolation is linear only when there are
    not enough points (i.e. beginning/end of spectrum).
    * method = "spline": A spline interpolation that may may be usefull
    for obtaining a smooth oversampled spectra by following this simple
    procedure:
        - If the original spectrum is not uniformly sampled, apply first
        the linear interpolation with a wave step equal to the median step
        of the spectrum.
        - Once the spectrum is uniformly sampled, apply the spline interpolation
        with a wave step smaller than the current wave step.

    """
    total_points = len(xaxis)
    last_reported_progress = -1

    current_work_progress = 10.0
    logging.info("%.2f%%" % current_work_progress)
    if frame != None:
        frame.update_progress(current_work_progress)
    resampled_spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    resampled_spectrum['waveobs'] = xaxis
    if method.lower() == "linear":
        ## Scipy linear interpolation (same result as numpy):
        #f = interpolate.interp1d(spectrum['waveobs'], spectrum['flux'], kind='linear', bounds_error=False, fill_value=0.0)
        #resampled_spectrum['flux'] = f(xaxis)
        ## Numpy linear interpolation:
        resampled_spectrum['flux'] = np.interp(xaxis, spectrum['waveobs'], spectrum['flux'], left=0.0, right=0.0) # No extrapolation, just returns zeros
        current_work_progress = 90.0
        logging.info("%.2f%%" % current_work_progress)
        if frame != None:
            frame.update_progress(current_work_progress)
    elif method.lower() == "spline":
        f = interpolate.InterpolatedUnivariateSpline(spectrum['waveobs'], spectrum['flux'], k=3)
        resampled_spectrum['flux'] = f(xaxis)
        current_work_progress = 90.0
        logging.info("%.2f%%" % current_work_progress)
        if frame != None:
            frame.update_progress(current_work_progress)
    elif method.lower() == "bessel":
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
    else:
        raise Exception("Unknown method")

    resampled_spectrum['err'] = np.interp(xaxis, spectrum['waveobs'], spectrum['err'])
    return resampled_spectrum

def correct_velocity(spectrum, velocity):
    """
    Correct velocity in km/s.
    """
    # Speed of light in m/s
    c = 299792458.0

    # Correct wavelength scale for radial velocity
    # - Newtonian version:
    ##Radial/barycentric velocity from km/s to m/s
    ##velocity = velocity * 1000
    #spectrum['waveobs'] = spectrum['waveobs'] / ((velocity / c) + 1)
    # - Relativistic version:
    spectrum['waveobs'] = spectrum['waveobs'] * np.sqrt((1.-(velocity*1000.)/c)/(1.+(velocity*1000.)/c))
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


