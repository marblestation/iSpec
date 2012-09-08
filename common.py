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
from scipy.interpolate import UnivariateSpline
import numpy.lib.recfunctions as rfn # Extra functions
import numpy as np
import matplotlib.pyplot as plt
import calendar
import re
import shutil
import gzip
import tempfile
import os, errno
#import ipdb
import random
import sys


# Estimate the Signal-to-Noise ratio for a given spectrum
# - If the spectra is normalized, we can estimate the SNR by considering
#   fluxes near the continuum
def estimate_snr(flux, num_points=10, frame=None):
    # Avoid negative values and outliers
    flux = flux[flux > 0.0]
    #flux, f = sigma_clipping(flux, sig=3, meanfunc=np.median)
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
            if i % 2000 == 0:
                progress = ((i*1.0 / total_num_blocks) * 100.0)
                print "%.2f%%" % progress
                if frame != None:
                    frame.update_progress(progress)
        snr = np.asarray(snr)
    snr, s = sigma_clipping(snr, sig=3, meanfunc=np.median)
    estimated_snr = np.median(snr)
    return estimated_snr


## Select spectra with identified file path doing a inner join between narval.vr and liste_spectre.
## - It uses the name of the star, the date (true and previous day) and the SNR similarity
## - Writes the result to a txt file
## - Returns a numpy array with the selection
def select_and_create_spectra_list(rv_name = "input/narval.vr", spectra_list_name = "input/liste_spectres", output_file = "output/spectra_list.txt"):
    # Dictionaries for month transformation
    month_name2num = dict((v.lower(),k) for k,v in enumerate(calendar.month_abbr))
    month_num2name = dict((k,v.lower()) for k,v in enumerate(calendar.month_abbr))

    # Read radial velocity file and spectra paths file
    rv = asciitable.read(rv_name, names=['name', 'julian_date', 'radial_velocity', 'err', 'snr'])
    spectra_list = asciitable.read(spectra_list_name, names=['name', 'file', 'snr', 'ignore1', 'ignore2'], exclude_names=['ignore1', 'ignore2'])

    # Join the information of radial velocity with the path where the spectra is located
    # by using the star's name and the date (which correspond to a subdirectory of the file's path)
    selected_spectra = []
    discarded_spectra = []
    already_matched_paths = []
    for i in np.arange(len(rv)):
        # Transform from julian to normal python dates
        d = jd_to_calendar(rv['julian_date'][i])

        # Transform date to format: 05oct07
        date_string1 = "%.2i" % d.day + month_num2name[d.month] + str(d.year)[2:]
        date_string2 = "%.2i" % (int(d.day)-1) + month_num2name[d.month] + str(d.year)[2:]
        path = ""
        min_diff_snr = 9999.9
        # For every spectra with the same name
        for spec in spectra_list[spectra_list['name'] == rv['name'][i]]:
            # Ignore already matched paths
            if spec['file'] in already_matched_paths:
                continue
            # Match the date looking inside the file's path
            match1 = re.search(".*/%s/.*/.*\.s" % date_string1, spec['file']) != None
            # Some of them are in a directory named with the previous day
            match2 = re.search(".*/%s/.*/.*\.s" % date_string2, spec['file']) != None
            diff_snr = abs(float(spec['snr'] - rv['snr'][i]) / spec['snr'])
            if (match1 or match2):
                # If there is more than one match, get the one with a smaller snr difference
                # - I do not know why, but SNR are not identical in narval.vr and liste_spectres
                if diff_snr < min_diff_snr:
                    min_diff_snr = diff_snr
                    path = spec['file']

        # If min_diff_snr has not been set, there is no match for this entry
        if min_diff_snr == 9999.9:
            #print "File path not found:", "\t", rv['name'][i], "\t", date_string1
            discarded_spectra.append([rv['name'][i], date_string1, str(d), rv['radial_velocity'][i], rv['err'][i], rv['snr'][i]])
        else:
            if min_diff_snr > 0.05:
                print "Warning: File matched but with 5% of difference in SNR\n", "\t", rv['name'][i], "\t", date_string1, "%.2f" % min_diff_snr
            already_matched_paths.append(path)
            # Separate filename from dirname
            filename = path.split('/')[-1]
            filename_length = len(filename)
            dirname = path[:-filename_length]
            selected_spectra.append((rv['name'][i], date_string1, str(d), rv['radial_velocity'][i], rv['err'][i], rv['snr'][i], dirname, filename))


    if output_file != None:
        asciitable.write(selected_spectra, output=output_file, delimiter='\t', names=['name', 'date', 'datetime', 'radial_velocity', 'err', 'snr', 'dirname', 'filename'])
        asciitable.write(discarded_spectra, output=output_file+".discard", delimiter='\t', names=['name', 'date', 'datetime', 'radial_velocity', 'err', 'snr'])

    return np.array(selected_spectra, dtype=[('name', '|S20'), ('date', '|S8'), ('datetime', '|S26'), ('radial_velocity', float), ('err', float), ('snr', int), ('dirname', '|S100'), ('filename', '|S100')])


## Read the selected spectra list
def read_spectra_list(input_file = "output/spectra_list.txt"):
    return asciitable.read(input_file, delimiter='\t', names=['name', 'date', 'datetime', 'radial_velocity', 'err', 'snr', 'dirname', 'filename'])


## Make dir (mkdir -p)
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

## Copy and optionally compress the selected spectra to the current working directory,
## respecting the subdirectory structure.
## - It should be executed in vanoise
def copy_input_spectra(selected_spectra, compress = True, base_path = "/m2a/soubiran/SPECTRES/NARVAL/"):
    for spec in selected_spectra:
        # Ignore if the file has been already copied and compressed
        if os.path.exists("input/" + spec['dirname'] + spec['filename'] + ".gz"):
            continue

        if not os.path.exists("input/" + spec['dirname']):
            mkdir_p("input/" + spec['dirname'])

        print base_path + spec['dirname'] + spec['filename'], "=>", "input/" + spec['dirname']

        if not compress:
            shutil.copy(base_path + spec['dirname'] + spec['filename'], "input/" + spec['dirname'])
        else:
            f_in = open(base_path + spec['dirname'] + spec['filename'], 'rb')
            f_out = gzip.open("input/" + spec['dirname'] + spec['filename'] + ".gz", 'wb')
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()

## Create output dir structure
def prepare_output_dirs(selected_spectra):
    for spec in selected_spectra:
        if not os.path.exists("output/" + spec['dirname']):
            mkdir_p("output/" + spec['dirname'])


# Find duplicates in a column of a recarray
# This is a simplified version of:
#   import numpy.lib.recfunctions as rfn
#   rfn.find_duplicates(...)
def find_duplicates(a, key):
    a = np.asanyarray(a).ravel()
    # Get the sorting data (by selecting the corresponding field)
    base = a[key]
    # Get the sorting indices and the sorted data
    sortidx = base.argsort()
    sorteddata = base[sortidx]
    # Compare the sorting data
    flag = (sorteddata[:-1] == sorteddata[1:])
    flag = np.concatenate(([False], flag))
    # We need to take the point on the left as well (else we're missing it)
    flag[:-1] = flag[:-1] + flag[1:]
    duplicates = a[sortidx][flag]
    duplicates_index = sortidx[flag]
    return (duplicates, duplicates_index)


# Reduces outliers considering the mean value and 3 sigma (3*stdev),
# iterating until convergence
def sigma_clipping(data, sig=3, meanfunc=np.mean):
    last_total = len(data)

    # First iteration
    stdev = np.std(data)
    diff = data - meanfunc(data)
    sfilter = np.abs(diff) < sig*stdev
    current_total = len(data[sfilter])
    # Continue iterating until convergence (no more points are removed)
    while last_total > current_total:
        #print current_total, stdev
        last_total = current_total

        stdev = np.std(data[sfilter])
        diff = data - meanfunc(data[sfilter])
        sfilter = np.abs(diff) < sig*stdev

        current_total = len(data[sfilter])

    return data[sfilter], sfilter


# For an array of values, find local maximum values considering a window
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

# For an array of values, find local minimum values considering a window
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


# For an array of values, find the position of local maximum values considering only
# the next and previous elements, except they have the same value.
# In that case, the next/previous different value is checked. Therefore,
# find_local_max([1,2,3,3,2,1,4,3]) would return [2, 3, 6]
def find_local_max_values(x):
    ret = []
    n = len(x)
    m = 0;
    for i in np.arange(n):
        l_min = np.max([i-1, 0])
        #l_max = i-1
        #r_min = i+1
        #r_max = np.min([i+1, n-1])
        r_min = np.min([i+1, n-1])
        is_max = True

        # left side
        j = l_min
        # If value is equal, search for the last different value
        while j >= 0 and x[j] == x[i]:
            j -= 1

        if j < 0 or x[j] > x[i]:
            is_max = False

        # right side
        if is_max:
            j = r_min
            # If value is equal, search for the next different value
            while j < n and x[j] == x[i]:
                j += 1
            if j >= n or x[j] > x[i]:
                is_max = False

        if is_max:
            ret.append(i)
    return np.asarray(ret)

# For an array of values, find the position of local maximum values considering only
# the next and previous elements, except they have the same value.
# In that case, the next/previous different value is checked. Therefore,
# find_local_max([10,9,3,3,9,10,4,30]) would return [2, 3, 6]
def find_local_min_values(x):
    ret = []
    n = len(x)
    m = 0;
    for i in np.arange(n):
        l_min = np.max([i-1, 0])
        #l_max = i-1
        #r_min = i+1
        #r_max = np.min([i+1, n-1])
        r_min = np.min([i+1, n-1])
        is_min = True
        # left side
        j = l_min
        # If value is equal, search for the last different value
        while j >= 0 and x[j] == x[i]:
            j -= 1

        if j < 0 or x[j] < x[i]:
            is_min = False

        # right side
        if is_min:
            j = r_min
            # If value is equal, search for the next different value
            while j < n and x[j] == x[i]:
                j += 1

            if j >= n or x[j] < x[i]:
                is_min = False

        if is_min:
            ret.append(i)
    return np.asarray(ret)

## Returns spectra from a filename:
## - if the file does not exists, checks if it exists a compressed version (gzip)
## - it the file exists, it can be automatically uncompressed (gzip)
def read_spectra(spectra_filename, estimate_errors_if_not_present=False):
    #spectra_filename = "input/" + spectra_filename
    # If it is not compressed
    if os.path.exists(spectra_filename) and spectra_filename[-3:] != ".gz":
        try:
            spectra = asciitable.read(table=spectra_filename, names=['waveobs', 'flux', 'err'])
        except asciitable.core.InconsistentTableError as err:
            try:
                # If it fails, try indicating that data starts at line 2 (original NARVAL spectra need this)
                spectra = asciitable.read(table=spectra_filename, data_start=2, names=['waveobs', 'flux', 'err'])
            except asciitable.core.InconsistentTableError as err:
                # Try without error column
                spectra_tmp = asciitable.read(table=spectra_filename, names=['waveobs', 'flux'])
                spectra = np.recarray((len(spectra_tmp), ), dtype=[('waveobs', float),('flux', float),('err', float)])
                spectra['waveobs'] = spectra_tmp['waveobs']
                spectra['flux'] = spectra_tmp['flux']
                if estimate_errors_if_not_present:
                    print "Estimating errors based on estimated SNR..."
                    snr = estimate_snr(spectra['flux'])
                    spectra['err'] = spectra['flux'] / snr
                else:
                    spectra['err'] = np.zeros(len(spectra)) # Add a zeroed error column

    elif (os.path.exists(spectra_filename) and spectra_filename[-3:] == ".gz") or (os.path.exists(spectra_filename + ".gz")):
        if spectra_filename[-3:] != ".gz":
            spectra_filename = spectra_filename + ".gz"

        tmp_spec = tempfile.mktemp() + str(int(random.random() * 100000000))
        # Uncompress to a temporary file
        f_out = open(tmp_spec, 'wb')
        f_in = gzip.open(spectra_filename, 'rb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()

        try:
            spectra = asciitable.read(table=tmp_spec, names=['waveobs', 'flux', 'err'])
        except asciitable.core.InconsistentTableError as err:
            try:
                # If it fails, try indicating that data starts at line 2 (original NARVAL spectra need this)
                spectra = asciitable.read(table=tmp_spec, data_start=2, names=['waveobs', 'flux', 'err'])
            except asciitable.core.InconsistentTableError as err:
                # Try without error column
                spectra_tmp = asciitable.read(table=tmp_spec, names=['waveobs', 'flux'])
                spectra = np.recarray((len(spectra_tmp), ), dtype=[('waveobs', float),('flux', float),('err', float)])
                spectra['waveobs'] = spectra_tmp['waveobs']
                spectra['flux'] = spectra_tmp['flux']
                print "Estimating errors based on estimated SNR..."
                snr = estimate_snr(spectra['flux'])
                spectra['err'] = spectra['flux'] / snr
                #spectra['err'] = np.zeros(len(spectra)) # Add a zeroed error column
        os.remove(tmp_spec)

    # Filter invalid errors and fluxes
    # TODO: Decide if we should request a valid 'error' column
    #valid = (spectra['err'] > 0) & ~np.isnan(spectra['err']) & (spectra['flux'] > 0) & ~np.isnan(spectra['flux'])
    #valid = (spectra['flux'] > 0) & ~np.isnan(spectra['flux'])
    valid = ~np.isnan(spectra['flux'])

    # Find duplicate wavelengths
    dups, dups_index = find_duplicates(spectra, 'waveobs')

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
    spectra = spectra[valid]

    spectra.sort(order='waveobs') # Make sure it is ordered by wavelength

    return spectra

## Write spectra to file
def write_spectra(spectra, spectra_filename, compress=True):
    #spectra_filename = "output/" + spectra_filename
    if compress:
        if spectra_filename[-3:] != ".gz":
            spectra_filename = spectra_filename + ".gz"

        tmp_spec = tempfile.mktemp() + str(int(random.random() * 100000000))
        asciitable.write(spectra, output=tmp_spec, delimiter='\t')

        # Compress the temporary file
        f_in = open(tmp_spec, 'rb')
        f_out = gzip.open(spectra_filename, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        os.remove(tmp_spec)
    else:
        asciitable.write(spectra, output=spectra_filename, delimiter='\t')


#### Read & write regions
# Continuum
def read_continuum_regions(continuum_regions_filename):
    continuum_regions = asciitable.read(table=continuum_regions_filename, comment='#', names=['wave_base', 'wave_top'])
    return continuum_regions

def write_continuum_regions(continuum_regions, continuum_regions_filename):
    asciitable.write(continuum_regions, output=continuum_regions_filename, delimiter='\t')

# Lines
def read_line_regions(line_regions_filename):
    line_regions = asciitable.read(table=line_regions_filename, comment='#', names=['wave_peak', 'wave_base', 'wave_top', 'note'], quotechar="\"")
    return line_regions

def write_line_regions(line_regions, line_regions_filename):
    asciitable.write(line_regions, output=line_regions_filename, delimiter='\t', quotechar="\"")

# Segments
def read_segment_regions(segment_regions_filename):
    segment_regions = asciitable.read(table=segment_regions_filename, comment='#', names=['wave_base', 'wave_top'])
    return segment_regions

def write_segment_regions(segment_regions, segment_regions_filename):
    asciitable.write(segment_regions, output=segment_regions_filename, delimiter='\t')


def show_histogram(x, xlabel='Units', nbins=50):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # the histogram of the data
    n, bins, patches = ax.hist(x, nbins, normed=0, facecolor='green', alpha=0.75)
    # Mean
    l = plt.axvline(x = np.mean(x), linewidth=1, color='red')
    ax.annotate('Mean', xy=(np.mean(x), np.max(n)),  xycoords='data',
        xytext=(10, 20), textcoords='offset points',
        size=8,
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10",
                        edgecolor='black'),
        horizontalalignment='right', verticalalignment='top',
        )
    # Median
    l = plt.axvline(x = np.median(x), linewidth=1, color='orange')
    ax.annotate('Median', xy=(np.median(x), np.max(n)),  xycoords='data',
        xytext=(10, 35), textcoords='offset points',
        size=8,
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10",
                        edgecolor='black'),
        horizontalalignment='right', verticalalignment='top',
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    ax.grid(True)
    plt.show()




################################################################################
#### [start] Copied from astropysics.obsutils (if not, pyinstaller fails)
# http://packages.python.org/Astropysics/
# https://github.com/eteq/astropysics/blob/master/astropysics/obstools.py
"""
Offset between Julian Date and Modified Julian Date - e.g. mjd = jd - mjdoffset
"""

def jd_to_calendar(jd,rounding=1000000,output='datetime',gregorian=None,mjd=False):
    """
    Converts a julian date to a calendar date and time.

    :param jd:
        The Julian Date at which to compute the calendar date/time, a sequence
        of JDs, or None for the current date/time at the moment the function is
        called.
    :type jd: scalar, array-like, or None
    :param rounding:
        If non-0, Performs a fix for floating-point errors. It specifies the
        number of milliseconds by which to round the result to the nearest
        second. If 1000000 (one second), no milliseconds are recorded. If
        larger, a ValueError is raised.
    :type rounding: scalar
    :param output:
        Determines the format of the returned object and can be:

            * 'datetime'
                A list of :class:`datetime.datetime` objects in UTC will be
                returned. If the input is a scalar, a single object will be
                returned.
            * 'array'
                A Nx7 array will be returned of the form
                [(year,month,day,hr,min,sec,msec),...] unless the input was a
                scalar, in which case it will be a length-7 array.
            * 'fracarray'
                An Nx3 array (year,month,day) where day includes the decimal
                portion.

    :param gregorian:
        If True, the output will be in the Gregorian calendar. Otherwise, it
        will be Julian. If None, it will be assumed to switch over on October
        4/15 1582.
    :type gregorian: bool or None
    :param bool mjd:
        If True, the input is interpreted as a modified julian date instead of a
        standard julian date.

    :returns:
        The calendar date and time in a format determined by the `output`
        parameter (see above).

    :except ValueError:
        If `rounding` is larger than one second, or `output` is invalid.


    **Examples**

    >>> jd_to_calendar(2451545)
    datetime.datetime(2000, 1, 1, 12, 0, tzinfo=tzutc())
    >>> jd_to_calendar(2305812.5)
    datetime.datetime(1600, 12, 31, 0, 0, tzinfo=tzutc())
    >>> jd_to_calendar([2415020.5,2305447.5],output='array')
    array([[1900,    1,    1,    0,    0,    0,    0],
           [1600,    1,    1,    0,    0,    0,    0]])
    >>> jd_to_calendar(0.0,output='fracarray')
    array([[ -4.71200000e+03,   1.00000000e+00,   1.50000000e+00]])

    """
    import datetime
    from dateutil import tz

    if jd is None:
        jd = calendar_to_jd(datetime.datetime.now(tz.tzlocal()))

    jd = np.array(jd,copy=True,dtype=float)
    scalar = jd.shape == ()
    jd = jd.ravel()

    if mjd:
        jd += mjdoffset

    if rounding > 1000000:
        raise ValueError('rounding cannot exceed a second')
    elif rounding <= 0:
        jd += .5
    else:
        rounding = int(rounding)
        roundingfrac = rounding/86400000000
        jd += .5 + roundingfrac

    z = np.floor(jd).astype(int)
    dec = jd - z #fractional piece

    #fix slight floating-point errors if they hapepn TOOD:check
    dgtr1 = dec>=1.0
    dec[dgtr1] -= 1.0
    z[dgtr1] += 1


    if gregorian is None:
        gregorian = 2299161

    if gregorian is True:
        alpha = ((z-1867216.25)/36524.25).astype(int)
        z += 1 + alpha - alpha//4
    elif gregorian is False:
        pass
    else:
        gmask = z >= gregorian
        alpha = ((z[gmask]-1867216.25)/36524.25).astype(int)
        z[gmask] += 1 + alpha - alpha//4

    b = z + 1524
    c = ((b-122.1)/365.25).astype(int)
    d = (365.25*c).astype(int)
    e = ((b-d)/30.6001).astype(int)

    day = b - d - (30.6001*e).astype(int)

    mmask = e<14
    month = e
    month[mmask] -= 1
    month[~mmask] -= 13
    year = c
    year[month>2] -= 4716
    year[month<=2] -= 4715

    if output == 'fracarray':
        dec = dec-roundingfrac
        dec[dec<0]=0
        return np.array((year,month,day+dec)).T

    if rounding == 1000000:
        secdec = dec*86400
        sec = secdec.astype(int)
        min = sec//60
        sec -= 60*min
        hr = min//60
        min -= 60*hr
        #sec[sec==secdec] -= 1
        msec = None
    else:
        msec = (dec*86400000000.).astype('int64')
        if rounding > 0:
            div = (msec//1000000)*1000000
            toround = (msec - div)<(2*rounding)
            msec[toround] = div + rounding
            msec  -= rounding

        sec = msec//1000000
        msec -= 1000000*sec

        min = sec//60
        sec -= 60*min
        hr = min//60
        min -= 60*hr

    if output == 'datetime':
        tzi = tz.tzutc()
        if msec is None:
            ts = (year,month,day,hr%24,min%60,sec%60)
        else:
            ts = (year,month,day,hr%24,min%60,sec%60,msec%1000000)
        res = [datetime.datetime(*t,**dict(tzinfo=tzi)) for t in zip(*ts)]
    elif output == 'array':
        msec = np.zeros_like(sec) if msec is None else msec
        res = np.array([year,month,day,hr%24,min%60,sec%60,msec]).T
    else:
        raise ValueError('invlid output form '+str(output))
    if scalar:
        return res[0]
    else:
        return res




def calendar_to_jd(caltime,tz=None,gregorian=True,mjd=False):

    """
    Convert a calendar date and time to julian date.

    :param caltime:
        The date and time to compute the JD.  Can be in one of these forms:

            * A sequence of floats in the order (yr,month,day,[hr,min,sec]).
            * A sequence in the order (yr,month,day,[hr,min,sec]) where at least
               one of the elements is a sequence (a sequence will be returned).
            * A :class:`datetime.datetime` or :class:`datetime.date` object
            * A sequence of :class:`datetime.datetime` or :class:`datetime.date`
              objects (a sequence will be returned).
            * None : returns the JD at the moment the function is called.

        If the time is unspecified, it is taken to be noon (i.e. Julian Date =
        Julian Day Number)

    :param tz:
        Sets the time zone to assume for the inputs for conversion to UTC. Can
        be any of the following:

            * None
                No time zone conversion will occur unless `caltime` is given as
                :class:`datetime.datetime` or :class:`datetime.date` objects
                with `tzinfo`, in which case they will be converted to UTC using
                their own `tzinfo`.
            * a string
                Specifies a timezone name (resolved into a timezone using the
                :func:`dateutil.tz.gettz` function).
            * a scalar
                The hour offset of the timezone.
            * a :class:`datetime.tzinfo` object,
                This object will be used for timezone information.

    :param gregorian:
        If True, the input will be interpreted as in the Gregorian calendar.
        Otherwise, it will be Julian. If None, it will be assumed to switch over
        on October 4/15, 1582.
    :type gregorian: bool or None
    :param bool mjd:
        If True, a modified julian date is returned instead of the standard
        julian date.

    :returns: JD as a float, or a sequence of JDs if sequences were input.


    **Examples**

    >>> import datetime,dateutil
    >>> calendar_to_jd((2010,1,1))
    2455198.0
    >>> calendar_to_jd(datetime.datetime(2000,12,21,3,0,0))
    2451899.625
    >>> calendar_to_jd([2004,3,(5,6)])
    array([ 2453070.,  2453071.])
    >>> dates = [datetime.datetime(2004,3,5),datetime.datetime(2004,3,9)]
    >>> calendar_to_jd(dates)
    array([ 2453069.5,  2453073.5])
    >>> tz = dateutil.tz.tzoffset('2',3*3600)
    >>> calendar_to_jd((2010,1,1),tz)
    2455197.875


    """
    #Adapted from xidl  jdcnv.pro
    from datetime import datetime,date,tzinfo

    if caltime is None:
        from dateutil.tz import tzlocal
        datetimes = [datetime.now(tzlocal())]
        scalarout = True
    elif isinstance(caltime,datetime) or isinstance(caltime,date):
        datetimes = [caltime]
        scalarout = True
    elif all([isinstance(ct,datetime) or isinstance(ct,date) for ct in caltime]):
        datetimes = caltime
        scalarout = False
    else:
        datetimes = None
        caltime = list(caltime)
        if not (3 <= len(caltime) < 8):
            raise ValueError('caltime input sequence is invalid size')
        while len(caltime) < 7:
            if len(caltime) == 3:
                #make hours 12
                caltime.append(12*np.ones_like(caltime[-1]))
            else:
                caltime.append(np.zeros_like(caltime[-1]))
        yr,month,day,hr,min,sec,msec = caltime
        scalarout = all([np.shape(v) is tuple() for v in caltime])

    #if input objects are datetime objects, generate arrays
    if datetimes is not None:
        yr,month,day,hr,min,sec,msec = [],[],[],[],[],[],[]
        for dt in datetimes:
            if not hasattr(dt,'hour'):
                dt = datetime(dt.year,dt.month,dt.day,12)

            if tz is None:
                off = dt.utcoffset()
                if off is not None:
                    dt = dt - off

            yr.append(dt.year)
            month.append(dt.month)
            day.append(dt.day)
            hr.append(dt.hour)
            min.append(dt.minute)
            sec.append(dt.second)
            msec.append(dt.microsecond)



    yr = np.array(yr,dtype='int64',copy=False).ravel()
    month = np.array(month,dtype='int64',copy=False).ravel()
    day = np.array(day,dtype='int64',copy=False).ravel()
    hr = np.array(hr,dtype=float,copy=False).ravel()
    min = np.array(min,dtype=float,copy=False).ravel()
    sec = np.array(sec,dtype=float,copy=False).ravel()
    msec = np.array(msec,dtype=float,copy=False).ravel()

    #do tz conversion if tz is provided
    if isinstance(tz,basestring) or isinstance(tz,tzinfo):
        if isinstance(tz,basestring):
            from dateutil import tz
            tzi = tz.gettz(tz)
        else:
            tzi = tz

        utcoffset = []
        for t in zip(yr,month,day,hr,min,sec,msec):
            #microsecond from float component of seconds

            dt = datetime(*[int(ti) for ti in t],**dict(tzinfo=tzi))
            utcdt = dt.utcoffset()
            if utcdt is None:
                utcoffset.append(0)
            else:
                utcoffset.append(utcdt.days*24 + (utcdt.seconds + utcdt.microseconds*1e-6)/3600)
    else:
        utcoffset = tz

#    ly = ((month-14)/12).astype(int) #In leap years, -1 for Jan, Feb, else 0
#    jdn = day - 32075l + 1461l*(yr+4800l+ly)//4

#    jdn += 367l*(month - 2-ly*12)//12 - 3*((yr+4900l+ly)//100)//4

#    res = jdn + (hr/24.0) + min/1440.0 + sec/86400.0 - 0.5

    #this algorithm from meeus 2ed
    m3 = month < 3
    yr[m3] -= 1
    month[m3] += 12

    cen = yr//100

    if gregorian is None:
        gregorian = (1582,10,4)
    if gregorian is True:
        gregoffset = 2 - cen + cen//4
    elif gregorian is False:
        gregoffset = 0
    else:
        gregoffset = 2 - cen + cen//4
        gmask = (yr>gregorian[0])&(month>gregorian[1])&(day>gregorian[2])
        gregoffset[~gmask] = 0


    jdn = (365.25*(yr+4716)).astype(int) + \
          (30.6001*(month + 1)).astype(int) + \
               day + gregoffset - 1524.5
    res = jdn + hr/24.0 + min/1440.0 + sec/86400.0

    if mjd:
        res -= mjdoffset

    if np.any(utcoffset):
        res -= np.array(utcoffset)/24.0

    if scalarout:
        return res[0]
    else:
        return res
#### [end]   Astropysics


