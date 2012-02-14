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

## Returns spectra from a filename:
## - if the file does not exists, checks if it exists a compressed version (gzip)
## - it the file exists, it can be automatically uncompressed (gzip)
def read_spectra(spectra_filename):    
    #spectra_filename = "input/" + spectra_filename
    # If it is not compressed
    if os.path.exists(spectra_filename) and spectra_filename[-3:] != ".gz":
        try:
            spectra = asciitable.read(table=spectra_filename, delimiter=' ', names=['waveobs', 'flux', 'err'])
        except asciitable.core.InconsistentTableError, err:
            # If it fails, try indicating that data starts at line 2 (original NARVAL spectra need this)
            spectra = asciitable.read(table=spectra_filename, delimiter=' ', data_start=2, names=['waveobs', 'flux', 'err'])
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
            spectra = asciitable.read(table=tmp_spec, delimiter=' ', names=['waveobs', 'flux', 'err'])
        except asciitable.core.InconsistentTableError, err:
            # If it fails, try indicating that data starts at line 2 (original NARVAL spectra need this)
            spectra = asciitable.read(table=tmp_spec, delimiter=' ', data_start=2, names=['waveobs', 'flux', 'err'])
        os.remove(tmp_spec)
    
    spectra.sort(order='waveobs') # Make sure it is ordered by wavelength
    
    return spectra


## Write spectra to file
def write_spectra(spectra, spectra_filename, compress=True):
    #spectra_filename = "output/" + spectra_filename
    if compress:
        if spectra_filename[-3:] != ".gz":
            spectra_filename = spectra_filename + ".gz"
        
        tmp_spec = tempfile.mktemp() + str(int(random.random() * 100000000))
        asciitable.write(spectra, output=tmp_spec, delimiter=' ')
        
        # Compress the temporary file
        f_in = open(tmp_spec, 'rb')
        f_out = gzip.open(spectra_filename, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        os.remove(tmp_spec)
    else:
        asciitable.write(spectra, output=spectra_filename, delimiter=' ')


#### Read & write regions
# Continuum
def read_continuum_regions(continuum_regions_filename):
    continuum_regions = asciitable.read(table=continuum_regions_filename, delimiter='\t', comment='#', names=['wave_base', 'wave_top'])
    return continuum_regions

def write_continuum_regions(continuum_regions, continuum_regions_filename):
    asciitable.write(continuum_regions, output=continuum_regions_filename, delimiter='\t')

# Lines
def read_line_regions(line_regions_filename):
    line_regions = asciitable.read(table=line_regions_filename, delimiter='\t', comment='#', names=['wave_peak', 'wave_base', 'wave_top', 'note'])
    return line_regions

def write_line_regions(line_regions, line_regions_filename):
    asciitable.write(line_regions, output=line_regions_filename, delimiter='\t')

# Segments
def read_segment_regions(segment_regions_filename):
    segment_regions = asciitable.read(table=segment_regions_filename, delimiter='\t', comment='#', names=['wave_base', 'wave_top'])
    return segment_regions

def write_segment_regions(segment_regions, segment_regions_filename):
    asciitable.write(segment_regions, output=segment_regions_filename, delimiter='\t')





