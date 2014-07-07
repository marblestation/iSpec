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
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions
from astropy.io import fits as pyfits
import scipy.ndimage as ndi
from spectrum import *
from plotting import *
from common import *
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import log
import logging

def __read_fits_spectrum(spectrum_filename):
    """
    Reads the 'PRIMARY' HDU of the FITS file, considering that it contains the fluxes.

    The wavelength are derived from the headers, if not possible it checks if the
    data from the HDU contains 2 axes and takes the first as the wavelength
    and the second as the flux.

    It tries to find the errors in other HDU by checking the names and the length,
    if none are found then they are set to zero.

    Inspired by pyspeckit:
        https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/src/ae1e0714410b58905466740b04b54318d5f318f8/pyspeckit/spectrum/readers/fits_reader.py?at=default

    Finally, if nothing has worked, it searches for a binary table with 3 columns
    'AWAV', 'FLUXES' and 'SIGMA'.
    """
    hdulist = pyfits.open(spectrum_filename)

    data = hdulist['PRIMARY'].data
    hdr = hdulist['PRIMARY'].header

    if data is not None and (type(hdulist['PRIMARY']) is pyfits.hdu.image.PrimaryHDU or \
                                type(hdulist['PRIMARY']) is pyfits.hdu.image.ImageHDU):
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

        # Try to read World Coordinate System (WCS) that defines the wavelength grid
        if hdr.get(str('CD%s_%s' % (specaxis,specaxis))) is not None:
            wave_step = hdr['CD%s_%s' % (specaxis,specaxis)]
            wave_base = hdr['CRVAL%s' % (specaxis)]
            reference_pixel = hdr['CRPIX%s' % (specaxis)]
            logging.info("Using the FITS CD matrix.  PIX=%f VAL=%f DELT=%f UNIT=%s" % (reference_pixel,wave_base,wave_step,unit))
        elif hdr.get(str('CDELT%s' % (specaxis))) is not None:
            wave_step = hdr['CDELT%s' % (specaxis)]
            wave_base = hdr['CRVAL%s' % (specaxis)]
            reference_pixel = hdr['CRPIX%s' % (specaxis)]
            logging.info("Using the FITS CDELT value.  PIX=%f VAL=%f DELT=%f UNIT=%s" % (reference_pixel,wave_base,wave_step,unit))
        elif len(data.shape) > 1:
            logging.info("No CDELT or CD in header.  Assuming 2D input with 1st line representing the spectral axis.")
            # No valid WCS, try assuming first axis is the wavelength axis
            if hdr.get('CUNIT%s' % (specaxis)) is not None:
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
                print "Log scale"
            else:
                xconv = lambda v: ((v-reference_pixel+1)*wave_step+wave_base)
                waveobs = xconv(np.arange(len(flux)))

        num_measures = len(flux)
        spectrum = create_spectrum_structure(waveobs, flux)

        # Try to find the errors in the extensions (HDU different than the PRIMARY):
        spectrum['err'] = np.zeros(len(flux))
        for i in xrange(len(hdulist)):
            name = hdulist[i].name.upper()
            if name == str('PRIMARY') or len(hdulist[i].data.flatten()) != len(flux) or type(hdulist[i]) is pyfits.hdu.table.BinTableHDU:
                continue
            if 'IVAR' in name or 'IVARIANCE' in name:
                spectrum['err'] = np.sqrt(1. / hdulist[i].data.flatten()) # Not sure
                #spectrum['err'] = 1. / hdulist[i].data.flatten()
                break
            elif 'VAR' in name or 'VARIANCE' in name:
                spectrum['err'] = np.sqrt(hdulist[i].data.flatten()) # Not sure
                #spectrum['err'] = hdulist[i].data.flatten()
                break
            elif 'NOISE' in name or 'ERR' in name or 'SIGMA' in name:
                spectrum['err'] = hdulist[i].data.flatten()
                break

    elif data is None:
        # Try to find a binary table with an irregular spectra and 3 columns
        spectrum = None
        for i in xrange(len(hdulist)):
            if type(hdulist[i]) is pyfits.hdu.table.BinTableHDU:
                data = hdulist[i].data
                # iSpec binary table for irregular spectra
                try:
                    spectrum = create_spectrum_structure(data['AWAV'], data['FLUX'], data['SIGMA'])
                except:
                    continue
        if spectrum is None:
            raise Exception("Unknown FITS file format")
    else:
        raise Exception("Unknown FITS file format")

    hdulist.close()

    return spectrum

def __read_spectrum(spectrum_filename):
    try:
        spectrum = np.array([tuple(line.rstrip('\r\n').split("\t")) for line in open(spectrum_filename,)][1:], dtype=[('waveobs', float),('flux', float),('err', float)])
        if len(spectrum) == 0:
            raise Exception("Empty spectrum or incompatible format")
    except Exception as err:
        # try narval plain text format:
        # - ignores 2 first lines (header)
        # - ignores last line (empty)
        # - lines separated by \r
        # - columns separated by space
        try:
            narval = open(spectrum_filename,).readlines()[0].split('\r')
            spectrum = np.array([tuple(line.rstrip('\r').split()) for line in narval[2:-1]], dtype=[('waveobs', float),('flux', float),('err', float)])
            if len(spectrum) == 0:
                raise Exception("Empty spectrum or incompatible format")
        except:
            # try espadons plain text format:
            # - ignores 2 first lines (header)
            # - columns separated by space
            espadons = open(spectrum_filename,).readlines()
            spectrum = np.array([tuple(line.rstrip('\r').split()) for line in espadons[2:]], dtype=[('waveobs', float),('flux', float),('err', float)])

    if len(spectrum) == 0:
        raise Exception("Empty spectrum or incompatible format")
    return spectrum

def read_spectrum(spectrum_filename, apply_filters=True, sort=True):
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
    in the extensions of the FITS file. In case the PRIMARY is empty, it searches for binary tables
    with wavelengths, fluxes and errors.
    """
    # If it is not compressed
    if os.path.exists(spectrum_filename) and (spectrum_filename[-4:].lower() == ".fit" or spectrum_filename[-5:].lower() == ".fits") :
        spectrum = __read_fits_spectrum(spectrum_filename)
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

    if apply_filters:
        # Filtering...
        valid = ~np.isnan(spectrum['flux'])

        if len(spectrum[valid]) > 2:
            # Find duplicate wavelengths
            dups, dups_index = find_duplicates(spectrum, 'waveobs')

            # Filter all duplicates except the first one
            last_wave = None
            for i in np.arange(len(dups)):
                if last_wave is None:
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

    if sort:
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
        header = pyfits.Header()
        header.set('ORIGIN', "iSpec")
        #header.set('VERSION', "iSpec")
        header.set('UTCSAVED', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

        wave_diff = spectrum['waveobs'][1:] - spectrum['waveobs'][:-1]
        median_wave_step = np.median(wave_diff)
        if np.all(np.abs(wave_diff - median_wave_step) < 0.0000001):
            ### Regularly sampled spectrum
            primary_data = np.asarray(spectrum['flux'], dtype='float32')

            # Coordinates
            header.set('CUNIT1', "NM")
            header.set('CTYPE1', "AWAV") # Air wavelength
            #header.set('CD1_1', spectrum['waveobs'][1] - spectrum['waveobs'][0])
            header.set('CDELT1', spectrum['waveobs'][1] - spectrum['waveobs'][0])
            header.set('CRVAL1', spectrum['waveobs'][0])
            header.set('CRPIX1', 1)
            #
            header.set('NAXIS', 1)
            header.set('NAXIS1', len(spectrum['flux']))
            header.set('TTYPE1', "FLUX")
            header.set('TUNIT1', "COUNTS")

            primary_hdu = pyfits.PrimaryHDU(data=primary_data, header=header)

            # Add an HDU extension (image) with errors if they exist
            if np.any(spectrum['err'] != 0):
                # Error extension
                ext_data = np.asarray(spectrum['err'], dtype='float32')
                extheader = header.copy()
                extheader.set('TTYPE1', "FLUX")
                extheader.set('TUNIT1', "COUNTS")
                extheader.set('NAXIS', 1)
                extheader.set('NAXIS1', len(spectrum['err']))
                extension_hdu = pyfits.ImageHDU(data=ext_data, header=extheader, name="SIGMA")
                fits_format = pyfits.HDUList([primary_hdu, extension_hdu])
            else:
                fits_format = pyfits.HDUList(primary_hdu)
        else:
            # Since it is not regularly sampled, we use a BinTable extension with an empty primary
            primary_hdu = pyfits.PrimaryHDU(data=None, header=header)
            bintable_hdu = pyfits.BinTableHDU(spectrum)

            #bintable_hdu.columns.change_name('waveobs', 'AWAV')
            bintable_hdu.data.columns[0].name = "AWAV"
            bintable_hdu.header.set('TTYPE1', "AWAV")
            bintable_hdu.header.set('TUNIT1', "NM")
            #bintable_hdu.columns.change_name('flux', 'FLUX')
            bintable_hdu.data.columns[1].name = "FLUX"
            bintable_hdu.header.set('TTYPE2', "FLUX")
            bintable_hdu.header.set('TUNIT2', "COUNTS")
            #bintable_hdu.columns.change_name('err', 'SIGMA')
            bintable_hdu.data.columns[2].name = "SIGMA"
            bintable_hdu.header.set('TTYPE3', "SIGMA")
            bintable_hdu.header.set('TUNIT3', "COUNTS")

            fits_format = pyfits.HDUList([primary_hdu, bintable_hdu])


        fits_format.writeto(spectrum_filename, clobber=True)
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

def normalize_spectrum(spectrum, continuum_model, consider_continuum_errors=True):
    """
    Normalizes a spectrum given a continuum fit
    """
    continuum_flux = continuum_model(spectrum['waveobs'])
    continuum_errors = continuum_model.placement_errors(spectrum['waveobs'])
    if np.all(continuum_flux == 1) and (not consider_continuum_errors or np.all(continuum_errors == 0)):
        return spectrum.copy()

    normalized_spectrum = create_spectrum_structure(spectrum['waveobs'])
    zeros = continuum_flux == 0
    normalized_spectrum['flux'][zeros] = 1.
    normalized_spectrum['err'][zeros] = 0.
    normalized_spectrum['flux'][~zeros] = spectrum['flux'][~zeros] / continuum_flux[~zeros]

    # Error propagation considering errors in continuum (most conservative operation)
    if (consider_continuum_errors and np.all(continuum_errors == 0)) or np.all(spectrum['err'] == 0):
        normalized_spectrum['err'] = 0.
    else:
        if consider_continuum_errors:
            zeros = np.logical_or(spectrum['flux'] == 0, continuum_flux == 0)
            normalized_spectrum['err'][~zeros] = normalized_spectrum['flux'][~zeros] * ((spectrum['err'][~zeros] / spectrum['flux'][~zeros]) + (continuum_errors[~zeros] / continuum_flux[~zeros]))
            # Non-correlated errors:
            #normalized_spectrum['err'][~zeros] = normalized_spectrum['flux'][~zeros] * np.sqrt(np.power(spectrum['err'][~zeros] / spectrum['flux'][~zeros], 2) + np.power(continuum_errors[~zeros] / continuum_flux[~zeros], 2))
            normalized_spectrum['err'][zeros] = 0.
        else:
            normalized_spectrum['err'][~zeros] = spectrum['err'][~zeros] / continuum_flux[~zeros]


    return normalized_spectrum

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
                if frame is not None:
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

try:
    import pyximport
    import numpy as np
    pyximport.install(setup_args={'include_dirs':[np.get_include()]})
    from spectrum_c import convolve_spectrum as __convolve_spectrum
    from spectrum_c import interpolation as __interpolation
except:
    print "*********************************************************************"
    print "Not optimized version loaded!"
    print "*********************************************************************"

    def __interpolation(waveobs, fluxes, err, resampled_waveobs, bessel=False, zero_edges=True, frame=None):
        """
        Interpolate flux for a given wavelength by using Bessel's Central-Difference Interpolation.
        It considers:

        - 4 points in general
        - 2 when there are not more (i.e. at the beginning of the array or outside)

        * It does not interpolate if any of the fluxes used for interpolation is zero or negative
        this way it can respect gaps in the spectrum
        """
        last_reported_progress = -1
        current_work_progress = 10.0

        total_points = len(waveobs)
        new_total_points = len(resampled_waveobs)
        resampled_flux = np.zeros(new_total_points)
        resampled_err = np.zeros(new_total_points)
        from_index = 0 # Optimization: discard regions already processed
        for i in np.arange(new_total_points):
            # Target wavelength
            objective_wavelength = resampled_waveobs[i]

            # Find the index position of the first wave length equal or higher than the objective
            index = waveobs[from_index:].searchsorted(objective_wavelength)
            index += from_index

            if index == total_points:
                # DISCARD: Linear extrapolation using index-1 and index-2
                # flux = fluxes[index-1] + (objective_wavelength - waveobs[index-1]) * ((fluxes[index-1]-fluxes[index-2])/(waveobs[index-1]-waveobs[index-2]))
                if zero_edges:
                    # JUST ZERO:
                    resampled_flux[i] = 0.0
                    resampled_err[i] = 0.0
                else:
                    # JUST DUPLICATE:
                    resampled_flux[i] = fluxes[index-1]
                    resampled_err[i] = err[index-1]
            #elif index == 0 and waveobs[index] != objective_wavelength:
            elif index == 0:
                # DISCARD: Linear extrapolation using index+1 and index
                # flux = fluxes[index] + (objective_wavelength - waveobs[index]) * ((fluxes[index+1]-fluxes[index])/(waveobs[index+1]-waveobs[index]))
                if zero_edges:
                    # JUST ZERO:
                    resampled_flux[i] = 0.0
                    resampled_err[i] = 0.0
                else:
                    # JUST DUPLICATE:
                    resampled_flux[i] = fluxes[index]
                    resampled_err[i] = err[index]
            # Do not do this optimization because it can produce a value surounded
            # by zeros because of the condition "Do not interpolate if any of the
            # fluxes is zero or negative" implemented in the rest of the cases
            #elif waveobs[index] == objective_wavelength:
                #resampled_flux[i] = fluxes[index]
                #resampled_err[i] = err[index]
            else:
                if not bessel or index == 1 or index == total_points-1:
                    # Do not interpolate if any of the fluxes is zero or negative
                    if fluxes[index-1] <= 1e-10 or fluxes[index] <= 1e-10:
                        resampled_flux[i] = 0.0
                        resampled_err[i] = 0.0
                    else:
                        # Linear interpolation between index and index-1
                        # http://en.wikipedia.org/wiki/Linear_interpolation#Linear_interpolation_between_two_known_points
                        d1 = (objective_wavelength - waveobs[index-1])
                        d2 = (waveobs[index]-waveobs[index-1])
                        resampled_flux[i] = fluxes[index-1] + d1 * ((fluxes[index]-fluxes[index-1])/d2)
                        # Same formula as for interpolation but I have re-arranged the terms to make
                        # clear that it is valid for error propagation (sum of errors multiplied by constant values)
                        resampled_err[i] = (err[index-1] * (d2 - d1)  + (err[index] * d1)) / d2
                        # Do not allow negative fluxes or errors
                        if resampled_err[i] < 0:
                            resampled_err[i] = 1e-10
                        if resampled_flux[i] < 0:
                            resampled_flux[i] = 0
                            resampled_err[i] = 0
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

                    err_x_1 = err[index - 2]
                    err_x0 = err[index - 1]
                    err_x1 = err[index]
                    err_x2 = err[index + 1]

                    # Do not interpolate if any of the fluxes is zero or negative
                    if flux_x_1 <= 1e-10 or flux_x0 <= 1e-10 or flux_x1 <= 1e-10 or flux_x2 <= 1e-10:
                        resampled_flux[i] = 0.0
                        resampled_err[i] = 0.0
                    else:
                        p = (objective_wavelength - wave_x0) / (wave_x1 - wave_x0)
                        factor = (p * (p - 1) / 4)
                        resampled_flux[i] = flux_x0 + p * (flux_x1 - flux_x0) + factor * (flux_x2 - flux_x1 - flux_x0 + flux_x_1)
                        # Same formula as for interpolation but I have re-arranged the terms to make
                        # clear that it is valid for error propagation (sum of errors multiplied by constant values)
                        resampled_err[i] = err_x_1 * factor + err_x0 * (1 - p - factor) + err_x1 * (p - factor) + err_x2 * factor
                        # Do not allow negative fluxes or errors
                        if resampled_err[i] < 0:
                            resampled_err[i] = 1e-10
                        if resampled_flux[i] < 0:
                            resampled_flux[i] = 0
                            resampled_err[i] = 0

            if index > 4:
                from_index = index - 4

            current_work_progress = np.min([(i*1.0 / new_total_points) * 100, 90.0])
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                logging.info("%.2f%%" % current_work_progress)
                if frame is not None:
                    frame.update_progress(current_work_progress)

        return resampled_waveobs, resampled_flux, resampled_err

    def __convolve_spectrum_slow(waveobs, flux, err, to_resolution, from_resolution=None, frame=None):
        """
        Slower implementation of the convolution but easier to understand
        """
        total_points = len(waveobs)
        convolved_flux = np.zeros(total_points)
        convolved_err = np.zeros(total_points)

        last_reported_progress = -1
        if frame is not None:
            frame.update_progress(0)

        # FWHM of the gaussian for the given resolution
        if from_resolution is None:
            # Convolve using instrumental resolution (smooth but not degrade)
            fwhm = waveobs / to_resolution
        else:
            # Degrade resolution
            fwhm = __get_fwhm(waveobs, from_resolution, to_resolution)
        sigma = __fwhm_to_sigma(fwhm)

        for i in np.arange(total_points):
            if flux[i] <= 1e-10:
                continue

            lambda_peak = waveobs[i] # Current lambda (wavelength) to be modified

            # Only work with a limited window considering 2 times the fwhm in each side of the current
            # position to be modified and saved in the convolved spectra
            wave_filter_base = waveobs.searchsorted(lambda_peak - 2*fwhm[i])
            wave_filter_top = waveobs.searchsorted(lambda_peak + 2*fwhm[i])
            waveobs_window = waveobs[wave_filter_base:wave_filter_top]
            flux_window = flux[wave_filter_base:wave_filter_top]

            # Construct the gaussian
            gaussian = np.exp(- ((waveobs_window - lambda_peak)**2) / (2*sigma[i]**2)) / np.sqrt(2*np.pi*sigma[i]**2)
            gaussian = gaussian / np.sum(gaussian)
            # Convolve
            convolved_flux[i] = np.sum(flux_window * gaussian)
            if err[i] > 0:
                # * Propagate error Only if the current value has a valid error value assigned
                #
                # Error propagation considering that measures are dependent (more conservative approach)
                # because it is common to find spectra with errors calculated from a SNR which
                # at the same time has been estimated from all the measurements in the same spectra
                #
                err_window = err[wave_filter_base:wave_filter_top]
                convolved_err[i] = np.sum(err_window * gaussian)

            current_work_progress = (i*1.0 / total_points) * 100
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                logging.info("%.2f%%" % current_work_progress)
                if frame is not None:
                    frame.update_progress(current_work_progress)

        return waveobs, convolved_flux, convolved_err

    def __convolve_spectrum(waveobs, flux, err, to_resolution, from_resolution=None, frame=None):
        """
        Spectra resolution smoothness/degradation. Procedure:

        1) Define a bin per measure which marks the wavelength range that it covers.
        2) For each point, identify the window segment to convolve by using the bin widths and the FWHM.
        3) Build a gaussian using the sigma value and the wavelength values of the spectrum window.
        4) Convolve the spectrum window with the gaussian and save the convolved value.

        If "from_resolution" is not specified or its equal to "to_resolution", then the spectrum
        is convolved with the instrumental gaussian defined by "to_resolution".

        If "to_resolution" is specified, the convolution is made with the difference of
        both resolutions in order to degrade the spectrum.
        """
        if from_resolution is not None and from_resolution <= to_resolution:
            raise Exception("This method cannot deal with final resolutions that are bigger than original")

        total_points = len(waveobs)
        convolved_flux = np.zeros(total_points)
        convolved_err = np.zeros(total_points)

        last_reported_progress = -1
        if frame is not None:
            frame.update_progress(0)

        # Consider the wavelength of the measurements as the center of the bins
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
        if from_resolution is None:
            # Convolve using instrumental resolution (smooth but not degrade)
            fwhm = waveobs / to_resolution
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

        # In theory, len(nbins) == len(waveobs)
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

            # Build the gaussian corresponding to the instrumental spread function
            gaussian = np.exp(- ((waveobs_segment - current_center)**2) / (2*current_sigma**2)) / np.sqrt(2*np.pi*current_sigma**2)
            gaussian = gaussian / np.sum(gaussian)

            # Convolve the current position by using the segment and the gaussian
            if flux[i] > 0:
                # Zero or negative values are considered as gaps in the spectrum
                only_positive_fluxes = flux_segment > 0
                weighted_flux = flux_segment[only_positive_fluxes] * gaussian[only_positive_fluxes]
                current_convolved_flux = weighted_flux.sum()
                convolved_flux[i] = current_convolved_flux
            else:
                convolved_err[i] = 0.0

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
                convolved_err[i] = current_convolved_err
            else:
                convolved_err[i] = 0.0

            current_work_progress = (i*1.0 / total_points) * 100
            if report_progress(current_work_progress, last_reported_progress):
                last_reported_progress = current_work_progress
                #logging.info("%.2f%%" % current_work_progress)
                if frame is not None:
                    frame.update_progress(current_work_progress)
        logging.info("Spectra convolved!")

        return waveobs, convolved_flux, convolved_err

def convolve_spectrum(spectrum, to_resolution, from_resolution=None, frame=None):
    """
    Spectra resolution smoothness/degradation.

    If "from_resolution" is not specified or its equal to "to_resolution", then the spectrum
    is convolved with the instrumental gaussian defined by "to_resolution".

    If "from_resolution" is specified, the convolution is made with the difference of
    both resolutions in order to degrade the spectrum.
    """
    if from_resolution is not None and from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")

    waveobs, flux, err = __convolve_spectrum(spectrum['waveobs'], spectrum['flux'], spectrum['err'], to_resolution, from_resolution=from_resolution, frame=frame)
    convolved_spectrum = create_spectrum_structure(waveobs, flux, err)
    return convolved_spectrum

def create_spectrum_structure(waveobs, flux=None, err=None):
    """
    Create spectrum structure
    """
    spectrum = np.recarray((len(waveobs), ), dtype=[('waveobs', float),('flux', float),('err', float)])
    spectrum['waveobs'] = waveobs

    if flux is not None:
        spectrum['flux'] = flux
    else:
        spectrum['flux'] = 0.0

    if err is not None:
        spectrum['err'] = err
    else:
        spectrum['err'] = 0.0

    return spectrum


def resample_spectrum(spectrum, xaxis, method="bessel", zero_edges=True, frame=None):
    """
    Returns a new spectrum with measures at the given xaxis wavelength
    Interpolation method can be:

    - method = "bessel": A Bessel's Central-Difference Interpolation with
      4 points. In this case interpolation is linear only when there are
      not enough points (i.e. beginning/end of spectrum).
    - method = "linear": Linear interpolation using 2 points
    - method = "spline": A spline interpolation that may may be usefull
      for obtaining a smooth oversampled spectra by following this simple
      procedure:

        - If the original spectrum is not uniformly sampled, apply first
          the linear interpolation with a wave step equal to the median step
          of the spectrum.
        - Once the spectrum is uniformly sampled, apply the spline interpolation
          with a wave step smaller than the current wave step.

        WARNING: spline method does not propagate correctly the errors
    - zero_edges = False, the first and last value in the edge will be propagated
      if the xaxis is bigger than spectrum['waveobs']. If it is set to True, then
      zero values will be assigned instead.


    """
    total_points = len(xaxis)
    last_reported_progress = -1

    current_work_progress = 10.0
    logging.info("%.2f%%" % current_work_progress)
    if frame is not None:
        frame.update_progress(current_work_progress)

    if method.lower() == "linear":
        ## Scipy linear interpolation (same result as numpy):
        #f = interpolate.interp1d(spectrum['waveobs'], spectrum['flux'], kind='linear', bounds_error=False, fill_value=0.0)
        #flux = f(xaxis)
        ## Numpy linear interpolation:
        #flux = np.interp(xaxis, spectrum['waveobs'], spectrum['flux'], left=0.0, right=0.0) # No extrapolation, just returns zeros
        #err = np.interp(xaxis, spectrum['waveobs'], spectrum['err'], left=0.0, right=0.0) # No extrapolation, just returns zeros
        #current_work_progress = 90.0
        #logging.info("%.2f%%" % current_work_progress)
        #if frame is not None:
            #frame.update_progress(current_work_progress)

        # iSpec linear interpolation:
        waveobs, flux, err = __interpolation(spectrum['waveobs'], spectrum['flux'], spectrum['err'], xaxis, bessel=False, zero_edges=zero_edges, frame=frame)
    elif method.lower() == "spline":
        logging.warn("Spline interpolation not recommended.")
        f = interpolate.InterpolatedUnivariateSpline(spectrum['waveobs'], spectrum['flux'], k=3)
        e = interpolate.InterpolatedUnivariateSpline(spectrum['waveobs'], spectrum['err'], k=3)
        flux = f(xaxis)
        err = e(xaxis)
        current_work_progress = 90.0
        logging.info("%.2f%%" % current_work_progress)
        if frame is not None:
            frame.update_progress(current_work_progress)
    elif method.lower() == "bessel":
        waveobs, flux, err = __interpolation(spectrum['waveobs'], spectrum['flux'], spectrum['err'], xaxis, bessel=True, zero_edges=zero_edges, frame=frame)
    else:
        raise Exception("Unknown method")

    resampled_spectrum = create_spectrum_structure(xaxis, flux, err)
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

def add_noise(spectrum, snr, distribution="poisson"):
    """
    Add noise to the spectrum fluxes to simulate a given SNR (Signal to noise Ratio).
    The distribution can be "poisson" or "gaussian"
    """
    noisy_spectrum = create_spectrum_structure(spectrum['waveobs'], spectrum['flux'], spectrum['err'])
    if distribution.lower() == "gaussian":
        sigma = spectrum['flux']/snr
        sigma[sigma <= 0.0] = 1.0e-10
        noisy_spectrum['flux'] += np.random.normal(0, sigma, len(spectrum))
        noisy_spectrum['err'] += sigma
    else:
        # poison
        sigma = 1./snr
        lamb = spectrum['flux'] / (sigma * sigma)
        lamb[lamb < 0.0] = 0.0
        noisy_spectrum['flux'] = np.random.poisson(lamb)
        noisy_spectrum['flux'] *= sigma*sigma
        noisy_spectrum['err'] += noisy_spectrum['flux'] / np.sqrt(lamb)
    return noisy_spectrum

def random_realizations(spectrum, number, distribution="poisson"):
    """
    Derive a group of spectra from a single spectrum by considering fluxes as
    mean (mu), errors as standard deviations (sigma) and a distribution.
    The distribution can be "poisson" or "gaussian"
    """
    realizations = []
    for i in xrange(number):
        new_derived_spectrum = create_spectrum_structure(spectrum['waveobs'], spectrum['flux'], spectrum['err'])
        if distribution.lower() == "gaussian":
            sigma = spectrum['err']
            sigma[sigma <= 0.0] = 1.0e-10
            new_derived_spectrum['flux'] += np.random.normal(0, sigma, len(spectrum))
            #new_derived_spectrum['err'] += sigma
        else:
            # poison
            sigma = spectrum['err']
            lamb = spectrum['flux'] / np.power(sigma, 2)
            lamb[lamb < 0.0] = 0.0
            new_derived_spectrum['flux'] = np.random.poisson(lamb)
            new_derived_spectrum['flux'] *= sigma*sigma
            #new_derived_spectrum['err'] += new_derived_spectrum['flux'] / np.sqrt(lamb)
        realizations.append(new_derived_spectrum)
    return realizations

def create_wavelength_filter(spectrum, wave_base=None, wave_top=None, regions=None):
    """
    Filter out the measurement of a spectrum between two given wavelengths limits
    or outside the regions if it is not "None".
    """
    wfilter = None
    if regions is None:
        if wave_base is None:
            wave_base = np.min(spectrum['waveobs'])
        if wave_top is None:
            wave_top = np.max(spectrum['waveobs'])
        wfilter = (spectrum['waveobs'] >= wave_base) & (spectrum['waveobs'] <= wave_top)
    elif len(regions) == 0:
        wfilter = spectrum['waveobs'] == np.min(spectrum['waveobs']) - 1.0 # Make all false
    else:
        # Build wavelength points from regions
        for region in regions:
            wave_base = region['wave_base']
            wave_top = region['wave_top']

            if wfilter is None:
                wfilter = np.logical_and(spectrum['waveobs'] >= wave_base, spectrum['waveobs'] <= wave_top)
            else:
                wfilter = np.logical_or(wfilter, np.logical_and(spectrum['waveobs'] >= wave_base, spectrum['waveobs'] <= wave_top))
    return wfilter


def air_to_vacuum(spectrum):
    """
    It converts spectrum's wavelengths (nm) from air to vacuum
    """
    converted_spectrum = create_spectrum_structure(spectrum['waveobs'], spectrum['flux'], spectrum['err'])
    sigma2 = np.power(1.e3/spectrum['waveobs'], 2) # nm
    # Compute conversion factor
    fact = 1. + 6.4328e-5 + 2.94981e-2/(146.-sigma2) + 2.5540e-4/(41.-sigma2)
    fact = fact*(spectrum['waveobs'] >= 200.) + 1.*(spectrum['waveobs'] < 200.)
    converted_spectrum['waveobs'] = spectrum['waveobs']*fact
    return converted_spectrum


def vacuum_to_air(spectrum):
    """
    It converts spectrum's wavelengths from vacuum to air
    """
    converted_spectrum = create_spectrum_structure(spectrum['waveobs'], spectrum['flux'], spectrum['err'])
    wave2 = np.power(spectrum['waveobs'], 2)
    # Compute conversion factor
    fact = 1. + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2**2.)
    fact = fact * ( spectrum['waveobs'] >= 200. ) + 1.*( spectrum['waveobs'] < 200. )
    converted_spectrum['waveobs'] = spectrum['waveobs'] / fact
    return converted_spectrum


def create_filter_cosmic_rays(spectrum, continuum_model, resampling_wave_step=0.001, window_size=15, variation_limit=0.01):
    """
    It uses a median filter to smooth out single-measurement deviations. Then it uses
    sigma-clipping to remove large variations between the actual and smoothed image.

    For doing the comparison, the original spectrum should be resampled to have
    homonenous wave step.

    Only those detected cosmics above the continuum will be discarded.
    """
    import scipy.signal
    wavelengths = np.arange(np.min(spectrum['waveobs']), np.max(spectrum['waveobs']), resampling_wave_step)
    resampled_spectrum = resample_spectrum(spectrum, wavelengths)

    resampled_smooth = create_spectrum_structure(resampled_spectrum['waveobs'])
    resampled_smooth['flux'] = scipy.signal.medfilt(resampled_spectrum['flux'], 15)
    smooth = resample_spectrum(resampled_smooth, spectrum['waveobs'])

    cosmics = (spectrum['flux'] - smooth['flux'])/ continuum_model(spectrum['waveobs']) > variation_limit
    cosmics = np.logical_and(cosmics, (spectrum['flux'] / continuum_model(spectrum['waveobs'])) > 1.0)
    return cosmics




