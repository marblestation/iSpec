#!/usr/bin/env python
#################
# Run with ipython -pdb -c "%run process.py"
#################
import ipdb
import asciitable
import numpy as np
from multiprocessing import Pool
from plotting import *
from continuum import *
from common import *
from radial_velocity import *
from interpolate import *
from fitting import *

## Homogenize spectra
# - Corrects radial velocity
# - Converts spectra to a common wavelength grid
def homogenize_spectra(spectra_list, base_wave = 370.0, top_wave = 1050.0, resolution = 65000):
    # Prepare a common wavelength grid for a given resolution
    xaxis = generate_wavelength_grid(base_wave, top_wave, resolution)
    
    for spec in spectra_list:
        rv_filename = "output/" + spec['dirname'] + spec['filename'] + ".rv" + ".gz"
        uniform_filename = "output/" + spec['dirname'] + spec['filename'] + ".common_grid" + ".gz"
        if not os.path.exists(uniform_filename):
            print "Reading %s" % "input/" + spec['dirname'] + spec['filename']
            spectra_data = read_spectra("input/" + spec['dirname'] + spec['filename'])

            print "Radial velocity correction..."
            radial_vel = spec['radial_velocity']
            spectra_data = correct_radial_velocity(spectra_data, radial_vel)
            #~ write_spectra(spectra_data, rv_filename, compress=True)
            #~ print "... saved %s" % rv_filename

            print "Resampling for a comon wavelength grid..."
            resampled_spectra_data = resample_spectra(spectra_data, xaxis)
            write_spectra(resampled_spectra_data, uniform_filename, compress=True)
            print "... saved %s" % uniform_filename
        else:
            print "Skip already existing file: %s" % uniform_filename
        #~ break

## Execute homogenize_spectra in several parallel processes
def concurrent_homogenize_spectra(spectra_list, num_processes = 8):
    p = Pool(num_processes)

    spec_per_process = int(np.round(len(spectra_list) / num_processes))
    arguments = []
    from_s = 0
    to_s = spec_per_process - 1
    for i in np.arange(num_processes):
        if i != num_processes - 1:
            arguments.append(spectra_list[from_s:to_s])
            from_s = to_s
            to_s = from_s + spec_per_process
        else:
            # For the last, include everything until the end
            arguments.append(spectra_list[from_s:])

    p.map(homogenize_spectra, arguments)


# Combine spectra which has been previously homogenize
# - It only consider those with SNR >= 100
def combine_spectra(spectra_list):
    ## Consider only spectra with SNR >= 100
    spectra_list = spectra_list[spectra_list["snr"] >= 100]

    ## Load uniform spectra
    # Read the first spectra as a reference
    spectra_ref = read_spectra("output/" + spectra_list[0]['dirname'] + spectra_list[0]['filename'] + ".common_grid.gz")

    total_spectra = len(spectra_list)
    total_wavelengths = len(spectra_ref)
    matrix = np.empty((total_spectra, total_wavelengths))

    for i in np.arange(total_spectra):
        print "Reading", spectra_list[i]['filename'], "(", i+1, "of", total_spectra, ")..."
        spectra_data = read_spectra("output/" + spectra_list[i]['dirname'] + spectra_list[i]['filename'] + ".common_grid.gz")
        matrix[i] = spectra_data['flux']

    ## Create a spectra as a result of a global combination
    print "Creating a cumulative spectra..."
    cumulative_mean_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    cumulative_median_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
    cumulative_mean_spectra['waveobs'] = spectra_ref['waveobs'] # Wavelengths from the reference spectra
    cumulative_median_spectra['waveobs'] = spectra_ref['waveobs']

    cumulative_mean_spectra['flux'] = np.mean(matrix, axis=0)
    cumulative_median_spectra['flux'] = np.median(matrix, axis=0)
    std = np.std(matrix, axis=0)
    cumulative_mean_spectra['err'] = std
    cumulative_median_spectra['err'] = std

    return cumulative_mean_spectra, cumulative_median_spectra


#~ if __name__ == '__main__':

#### Select which spectra to use
#spectra_list = select_and_create_spectra_list(output_file = "output/spectra_list.txt") # Run just once

#### Read spectra list previously selected
spectra_list = read_spectra_list(input_file = "output/spectra_list.txt")

#### Input and output preparations
#copy_input_spectra(spectra_list) # Needs to be executed in vanoise (but only once)
prepare_output_dirs(spectra_list)

#### Homogenize correcting radial velocity and imposing a common wavelength grid
#~ concurrent_homogenize_spectra(spectra_list, num_processes = 8)
#homogenize_spectra(spectra_list)


#### Combine homogenized spectra with SNR >= 100
#~ cumulative_mean_spectra, cumulative_median_spectra = combine_spectra(spectra_list)
#~ write_spectra(cumulative_mean_spectra, "output/cumulative_mean_spectra.s.common_grid.gz", compress=True)
cumulative_mean_spectra = read_spectra("output/cumulative_mean_spectra.s.common_grid.gz")
#~ write_spectra(cumulative_median_spectra, "output/cumulative_median_spectra.s.common_grid.gz", compress=True)
cumulative_median_spectra = read_spectra("output/cumulative_median_spectra.s.common_grid.gz")

#~ plot_spectra([cumulative_mean_spectra, cumulative_median_spectra])


#### Fit continuum model
continuum_model = fit_continuum(cumulative_mean_spectra)
#~ spectra_continuum_mean = get_spectra_from_model(continuum_model, cumulative_mean_spectra['waveobs'])
#~ write_spectra(spectra_continuum_mean, "output/cumulative_spectra_continuum_mean.gz", compress=True)
spectra_continuum_mean = read_spectra("output/cumulative_spectra_continuum_mean.gz")

continuum_model = fit_continuum(cumulative_median_spectra)
#~ spectra_continuum_median = get_spectra_from_model(continuum_model, cumulative_median_spectra['waveobs'])
#~ write_spectra(spectra_continuum_median, "output/cumulative_spectra_continuum_median.gz", compress=True)
spectra_continuum_median = read_spectra("output/cumulative_spectra_continuum_median.gz")


#### Find continuum blocks
#~ resolution = 65000
#~ continuum_mean_blocks = find_continuum(cumulative_mean_spectra, resolution, continuum_model=continuum_model)
#~ write_continuum_blocks(continuum_mean_blocks, "output/cumulative_mean.continuum.txt")
continuum_mean_blocks = read_continuum_blocks("output/cumulative_mean.continuum.txt")
#~ #continuum_mean_blocks = merge_blocks(cumulative_mean_spectra, continuum_mean_blocks)
#~ plot_spectra([cumulative_mean_spectra, spectra_continuum_mean], continuum=continuum_mean_blocks)

#~ resolution = 65000
#~ continuum_median_blocks = find_continuum(cumulative_median_spectra, resolution, continuum_model=continuum_model)
#~ write_continuum_blocks(continuum_median_blocks, "output/cumulative_median.continuum.txt")
continuum_median_blocks = read_continuum_blocks("output/cumulative_median.continuum.txt")
#~ #continuum_median_blocks = merge_blocks(cumulative_median_spectra, continuum_median_blocks)
#~ plot_spectra([cumulative_median_spectra, spectra_continuum_median], continuum=continuum_median_blocks)


#### Find continuum blocks from "stable" regions considering low variation among the combined spectra
#~ stable_blocks = find_stable_blocks(cumulative_mean_spectra)
#~ merged_blocks = merge_blocks(cumulative_mean_spectra, stable_blocks)
#~ # At least 12 measures per block
#~ min_blocks = merged_blocks[merged_blocks['num_measures'] >= 12]
#~ plot_spectra([cumulative_mean_spectra, spectra_continuum_mean], continuum=min_blocks)

#~ resolution = 65000
#~ continuum_mean_blocks = find_continuum_block_limited(cumulative_mean_spectra, resolution, min_blocks, continuum_model=continuum_model)
#~ write_continuum_blocks(continuum_mean_blocks, "output/cumulative_mean.stable.continuum.txt")
continuum_mean_blocks = read_continuum_blocks("output/cumulative_mean.stable.continuum.txt")
plot_spectra([cumulative_mean_spectra, spectra_continuum_mean], continuum=continuum_mean_blocks)


#~ stable_blocks = find_stable_blocks(cumulative_median_spectra)
#~ merged_blocks = merge_blocks(cumulative_median_spectra, stable_blocks)
#~ # At least 12 measures per block
#~ min_blocks = merged_blocks[merged_blocks['num_measures'] >= 12]
#~ plot_spectra([cumulative_median_spectra, spectra_continuum_median], continuum=min_blocks)

#~ resolution = 65000
#~ continuum_blocks = find_continuum_block_limited(cumulative_median_spectra, resolution, min_blocks, continuum_model=continuum_model)
#~ write_continuum_blocks(continuum_median_blocks, "output/cumulative_median.stable.continuum.txt")
continuum_median_blocks = read_continuum_blocks("output/cumulative_median.stable.continuum.txt")
plot_spectra([cumulative_median_spectra, spectra_continuum_median], continuum=continuum_median_blocks)
