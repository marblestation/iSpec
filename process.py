#!/usr/bin/env python
#################
# Run with ipython -pdb -c "%run process.py"
#################
#import ipdb
import asciitable
import numpy as np
from multiprocessing import Pool
from plotting import *
from continuum import *
from common import *
from radial_velocity import *
from interpolate import *

def process_spectra(spectra_list):
    for spec in spectra_list:
        rv_filename = "output/" + spec['dirname'] + spec['filename'] + ".rv" + ".gz"
        uniform_filename = "output/" + spec['dirname'] + spec['filename'] + ".rv.uniform" + ".gz"
        if not os.path.exists(uniform_filename):
            print "Reading %s" % "input/" + spec['dirname'] + spec['filename']
            spectra_data = read_spectra("input/" + spec['dirname'] + spec['filename'])

            print "Radial velocity correction..."
            radial_vel = spec['radial_velocity']
            spectra_data = correct_radial_velocity(spectra_data, radial_vel)
            write_spectra(spectra_data, rv_filename, compress=True)
            print "... saved %s" % rv_filename

            print "Resampling for uniform spacing..."
            resampled_spectra_data = resample_spectra(spectra_data)
            write_spectra(resampled_spectra_data, uniform_filename, compress=True)
            print "... saved %s" % uniform_filename
        else:
            print "Skip already existing file: %s" % uniform_filename
        #~ break

def concurrent_process_spectra(spectra_list):
    num_processes = 8
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

    p.map(process_spectra, arguments)


#~ if __name__ == '__main__':

## Select which spectra to use
#pectra_list = select_and_create_spectra_list(output_file = "output/spectra_list.txt") # Run just once
## Read spectra list previously selected
spectra_list = read_spectra_list(input_file = "output/spectra_list.txt")

## Input and output preparations
#copy_input_spectra(spectra_list) # Needs to be executed in vanoise (but only once)
prepare_output_dirs(spectra_list)

##
#concurrent_process_spectra(spectra_list)

## Read spectra list previously selected
spectra_list = read_spectra_list(input_file = "output/spectra_list.txt")
## Consider only spectra with SNR > 100
spectra_list = spectra_list[spectra_list["snr"] >= 100]

## Load uniform spectra
# Read the first spectra as a reference
spectra_ref = read_spectra("output/" + spectra_list[0]['dirname'] + spectra_list[0]['filename'] + ".rv.uniform.gz")

total_spectra = len(spectra_list)
total_wavelengths = len(spectra_ref)
matrix = np.empty((total_spectra, total_wavelengths))

for i in np.arange(total_spectra):
    print "Reading", spectra_list[i]['filename'], "(", i+1, "of", total_spectra, ")..."
    spectra_data = read_spectra("output/" + spectra_list[i]['dirname'] + spectra_list[i]['filename'] + ".rv.uniform.gz")
    matrix[i] = spectra_data['flux']

## Create a spectra as a result of a global combination
print "Creating a cumulative spectra..."
cumulative_mean_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
cumulative_median_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
cumulative_mean_spectra['waveobs'] = spectra_ref['waveobs'] # Wavelengths from the reference spectra
cumulative_median_spectra['waveobs'] = spectra_ref['waveobs']

#matrix = np.vstack( (spectra_data[0]['flux'], spectra_data[1]['flux']) )
#matrix = np.vstack( (spectra_data[i]['flux'] for i in np.arange(total_spectra)) )

cumulative_mean_spectra['flux'] = np.mean(matrix, axis=0)
cumulative_median_spectra['flux'] = np.median(matrix, axis=0)
#std = np.std(matrix, axis=0) # Not used because it returns a MemoryError, replaced by a loop
for i in np.arange(total_spectra):
    std = np.std(matrix[i])
    cumulative_mean_spectra['err'] = std
    cumulative_median_spectra['err'] = std

print "Saving cumulative spectra..."
write_spectra(cumulative_mean_spectra, "output/cumulative_mean_spectra.s.rv.uniform.gz", compress=True)
write_spectra(cumulative_median_spectra, "output/cumulative_median_spectra.s.rv.uniform.gz", compress=True)

#~ plot_spectra([cumulative_mean_spectra, cumulative_median_spectra])

## Continuum
continuum_blocks1 = find_continuum(cumulative_mean_spectra, filename="output/cumulative_mean_spectra.s.rv.uniform")
continuum_blocks2 = find_continuum(cumulative_mean_spectra, filename="output/cumulative_median_spectra.s.rv.uniform")
#~ #plot_spectra(cumulative_mean_spectra, filename="output/cumulative_mean_spectra.s.rv.uniform", continuum=continuum_blocks1)
#~ plot_spectra([cumulative_mean_spectra], continuum=continuum_blocks1)


################ Fit continuum and lines
from pymodelfit import GaussianModel
from pymodelfit import QuadraticModel
spectra = read_spectra("input/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s")
#~ wave_base = 486
#~ wave_top = 486.2
wave_base = 485
wave_top = 487
wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)    
spectra_visible = spectra[wave_filter]

## Find continuum using a quadratic model and excluding outliers for the fit
continuum = QuadraticModel()
weights = 1.0 / spectra_visible['err']
mean = np.mean(spectra_visible['flux'])
# TODO: Consider continuum regions
filter_outliers = spectra_visible['flux'] < mean
weights[filter_outliers] = 0.01 # Downgrade importance to points that probably belong to absorption lines
continuum.fitData(spectra_visible['waveobs'], spectra_visible['flux'], weights=weights)
print continuum.pardict
#~ q.plot()

## Fit line
g = GaussianModel()
g.mu = 486.11
g.sig = 0.02
g.A = -0.025
#~ g.mu = 486.34
#~ g.sig = 0.02
#~ g.A = -0.025

weights = 1.0 / spectra_visible['err']
#~ g.fitData(spectra_visible['waveobs'], spectra_visible['flux'] - continuum(spectra_visible['waveobs']), weights=weights, fixedpars=['mu','sig','A'])
g.fitData(spectra_visible['waveobs'], spectra_visible['flux'] - continuum(spectra_visible['waveobs']), weights=weights)

print g.pardict
g.plot()



#~ from common import *
#~ from plotting import *
#~ spectra = []
#~ spectra.append(read_spectra("input/L111N03_spec_norm/09may11/sp2_Normal/4_vesta_001.s.gz"))
#~ spectra.append(read_spectra("output/L111N03_spec_norm/09may11/sp2_Normal/4_vesta_001.s.rv.gz"))
#~ spectra.append(read_spectra("output/L111N03_spec_norm/09may11/sp2_Normal/4_vesta_001.s.rv.uniform.gz"))
#~ spectra = []
#~ spectra.append(read_spectra("input/L101N05_spec_norm/06mar10/sp2_Normal/hd025329_001.s.gz"))
#~ spectra.append(read_spectra("output/L101N05_spec_norm/06mar10/sp2_Normal/hd025329_001.s.rv.gz"))
#~ spectra.append(read_spectra("output/L101N05_spec_norm/06mar10/sp2_Normal/hd025329_001.s.rv.uniform.gz"))
# -25.478
#~ spectra.append(read_spectra("output/L111N03_spec_norm/15may11/sp2_Normal/hip077348_001.s.rv.uniform.gz"))
# 1.868 HIP077348 15may11
#~ spectra.append(read_spectra(""))
#~ spectra.append(read_spectra(""))
#~ spectra = []
#~ spectra.append(read_spectra("output/L111N03_spec_norm/09may11/sp2_Normal/4_vesta_001.s.rv.uniform.gz"))
#~ spectra.append(read_spectra("output/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s.rv.uniform.gz"))
#~ ./L101N05_spec_norm/12mar10/sp2_Normal/hd122563_001.s.rv.uniform.gz
#~ spectra = []
#~ spectra.append(read_spectra("input/L111N03_spec_norm/09may11/sp2_Normal/4_vesta_001.s.gz"))
#~ spectra.append(read_spectra("input/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s.gz"))
#~ 
#~ spectra = []
#~ spectra.append(read_spectra("output/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s.rv.uniform.gz")) # blue -12.426 0.023 139
#~ spectra.append(read_spectra("output/L091N03_spec_norm/11mar09/sp2_Normal/hd140538_001.s.rv.uniform.gz")) # green 19.051 0.012 273
#~ spectra = []
#~ spectra.append(read_spectra("input/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s.gz"))
#~ spectra.append(read_spectra("input/./L091N03_spec_norm/11mar09/sp2_Normal/hd140538_001.s.gz"))
#~ 
#~ spectra = []
#~ spectra.append(read_spectra("input/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s.gz")) # blue -12.426 0.023 139
#~ spectra.append(read_spectra("output/L091N03_spec_norm/11mar09/sp2_Normal/hd140538_001.s.rv.gz")) # green 19.051 0.012 273
#~ plot_spectra(spectra)

#~ spectra_data1 = read_spectra("input/L091N03_spec_norm/23apr09/sp2_Normal/hd125184_001.s.gz")
#~ radial_vel = -12.426
#~ spectra_data1 = correct_radial_velocity(spectra_data1, radial_vel)
#~ write_spectra(spectra_data1, "output/hd125184_001.s.rv.gz", compress=True)
#~ resampled_spectra_data1 = resample_spectra(spectra_data1)
#~ write_spectra(resampled_spectra_data1, "output/hd125184_001.s.rv.uniform.gz", compress=True)
#~ 
#~ spectra_data2 = read_spectra("input/L091N03_spec_norm/11mar09/sp2_Normal/hd140538_001.s.gz")
#~ radial_vel = 19.051
#~ spectra_data2 = correct_radial_velocity(spectra_data2, radial_vel)
#~ write_spectra(spectra_data2, "output/hd140538_001.s.rv.gz", compress=True)
#~ resampled_spectra_data2 = resample_spectra(spectra_data2)
#~ write_spectra(resampled_spectra_data2, "output/hd140538_001.s.rv.uniform.gz", compress=True)
#~ 
#~ plot_spectra([resampled_spectra_data1, resampled_spectra_data2])

#~ base_dir = "/media/asteroid/Data/Universidad/Bordeaux1/Spectra/process/input/"
#~ #filename = "arcturus_001.s"
#~ filename = "hd4614_001.s"
#~ 
#~ spectra = asciitable.read(table=base_dir+filename, delimiter=' ', data_start=2, names=['waveobs', 'flux', 'err'])
#~ spectra.sort(order='waveobs') # Make sure it is ordered by wavelength
#~ 
#~ wave_initial_base = spectra['waveobs'].min()
#~ wave_final_top = spectra['waveobs'].max()
#~ wave_increase = (wave_final_top - wave_initial_base) / 4
#~ 
#~ wave_base = wave_initial_base
#~ wave_top = wave_base + wave_increase
#~ num = 0
#~ while wave_base < wave_final_top:
    #~ wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] < wave_top)
    #~ find_continuum(spectra[wave_filter], filename=filename+"."+str(num), plot=True)
    #~ wave_base = wave_top
    #~ wave_top = wave_base + wave_increase
    #~ num += 1

