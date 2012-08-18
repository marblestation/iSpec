"""
    This file is part of Spectra.
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com

    Spectra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spectra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Spectra.  If not, see <http://www.gnu.org/licenses/>.
"""
import asciitable
from scipy.interpolate import UnivariateSpline
import numpy.lib.recfunctions as rfn # Extra functions
import numpy as np
import matplotlib.pyplot as plt
from astropysics import obstools
import calendar
import re
import shutil
import gzip
import tempfile
import os, errno
#import ipdb
import random

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
    estimated_snr = np.mean(snr)
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
        d = obstools.jd_to_calendar(rv['julian_date'][i])

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



