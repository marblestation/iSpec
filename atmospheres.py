import os
import sys
import numpy as np
from subprocess import Popen, PIPE
import asciitable
import math
from datetime import datetime
from scipy import interpolate
import matplotlib.pyplot as plt
import tempfile
import cPickle as pickle

# SPECTRUM is compatible only with the plane-parallel atmospheres.
# The first layer represents the surface. 
# Elemental abundances in the stdatom.dat file are used (and scaled with the [M/H] value)

# Constant class used for microturbulent velocities because they are constant for all layers and atmospheres
class ConstantValue:
    def __init__(self, value):
        self.value = value
    
    def __call__(self, x, y):
        return self.value

# Read castelli et kurucz atmospheres and builds the following structure:
#  atmospheres_params (as many as different metallicities):
#   - atmospheres_params_with_same_metallicity
#       - atmospheres
#           -> params
#  atmospheres (as many as different metallicities):
#   - atmospheres_with_same_metallicity
#       - atmospheres
#           -> values
#               * Values' order:
#                   mass_depth = 0
#                   temperature_K = 1
#                   gas_preassure = 2
#                   electron_density = 3
#                   mean_absorption_coeff = 4
#                   radiation_preassure = 5
#                   microturbulent_vel = 6 # In m/s
def read_castelli_kurucz_atmospheres(atmosphere_models, atmosphere_models_metalicities):
    read_atmosphere_data = False
    atmospheres_params = []
    atmospheres = []
    for i in np.arange(len(atmosphere_models)):
        atmosphere_model = atmosphere_models[i]
        f = open(atmosphere_model)
        
        atmospheres_with_same_metallicity = []
        atmospheres_params_with_same_metallicity = []
        for line in f.readlines():
            vline = line.split()
            if vline[0] == "TEFF":
                teff = np.round(float(vline[1]), 2) # K
                logg = np.round(float(vline[3]), 2)
            elif vline[0] == "TITLE":
                metallicity = float(vline[1].strip("[").strip("]").strip("A")) # [M/H]
            elif vline[0] == "READ":
                # New atmosphere
                current_atmosphere = []
                read_atmosphere_data = True
            elif read_atmosphere_data:
                if vline[0] == "PRADK":
                    read_atmosphere_data = False
                    if len(current_atmosphere) == 71:
                        # There is one model with 71 layers instead of 72
                        #  - ap05k2odfnew.dat teff 5250.0 logg 0.0 metallicity 0.5
                        # So we duplicate the last layer in order to have 72
                        current_atmosphere.append(layer[0:7])
                    atmospheres_params_with_same_metallicity.append([teff, logg, metallicity])
                    atmospheres_with_same_metallicity.append(current_atmosphere)
                    
                else:
                    layer = line.split(" ")
                    # Clean empty values due to consecutive spaces
                    while True:
                        try:
                            layer.remove('')
                        except ValueError:
                            break
                    # Only use the 7 first values
                    current_atmosphere.append(layer[0:7])
        atmospheres.append(atmospheres_with_same_metallicity)
        atmospheres_params.append(atmospheres_params_with_same_metallicity)

    atmospheres_params = np.array(atmospheres_params)
    return atmospheres, atmospheres_params

# Calculate ranges, total number and min/max values for atmospheric models
def atmospheres_statistics(atmospheres_params):
    # Stats for the first atmosphere (but all are going to be the same)
    teff_range = np.unique(atmospheres_params[0][:,0])
    logg_range = np.unique(atmospheres_params[0][:,1])

    MH_range = []
    for atmospheres_params_with_same_metallicity in atmospheres_params:
        MH_range.append(float(atmospheres_params_with_same_metallicity[0][2]))

    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)
    nlayers = len(atmospheres[0][0])
    nvalues = len(atmospheres[0][0][0])

    teff_min = np.min(teff_range)
    teff_max = np.max(teff_range)
    logg_min = np.min(logg_range)
    logg_max = np.max(logg_range)
    MH_min = np.min(MH_range)
    MH_max = np.max(MH_range)
    
    return teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max


# Looks for what combinations of teff-logg we have a model and builds a structure:
#  - model_combinations (as many as nMH)
#     - model combination matrix
#       -> True  : It exists
#       -> False
def check_teff_logg_combinations(atmospheres_params, nteff, nlogg, nMH):
    model_combination_matrix = np.array([-1]*nteff*nlogg)
    model_combination_matrix = model_combination_matrix.reshape(nteff, nlogg)
    model_combinations = [model_combination_matrix] * nMH

    for metal_num in np.arange(nMH):
        atmospheres_params_with_same_metallicity = atmospheres_params[metal_num]
        for i in np.arange(len(atmospheres_params_with_same_metallicity)):
            atm_params = atmospheres_params_with_same_metallicity[i]
            # Mark model as read
            teff_index = np.where(teff_range==atm_params[0])[0][0]
            logg_index = np.where(logg_range==atm_params[1])[0][0]
            model_combinations[metal_num][teff_index][logg_index] = i
    return model_combinations


# Builds an structure where each value of each layer has a RectBivariateSpline (based on the values
# read from atmospheric models) that can be used for interpolation. The structure is:
# - Modeled layers (as many as different metallicities)
#    - Layers
#       - Modeled values
#           -> [Matrix value (for each teff-logg), Model for interpolation]
def modeled_interpolated_layer_values(atmospheres, model_combinations, nteff, nlogg, nMH, nlayers, nvalues, teff_range, logg_range):
    ### Construct base grid for interpolation
    grid_teff, grid_logg = np.indices((nteff, nlogg), dtype=float)
    # Assign values to each grid position
    for i in np.arange(nteff):
        grid_teff[i] = teff_range[i]
    for i in np.arange(nlogg):
        grid_logg[:,i] = logg_range[i]
    
    modeled_layers = []
    # For each atmosphere with the same metallicity
    for metal_num in np.arange(nMH):
        atmospheres_with_same_metallicity = atmospheres[metal_num]
        modeled_layers_with_same_metallicity = []
        print "\nAtmosphere models", metal_num
        print "\tLayer:",
        # For each layer, group and modelize the different atmospheres (teff, logg)
        for layer_num in np.arange(nlayers):
            modeled_values = []
            print layer_num,
            sys.stdout.flush()
            for value_num in np.arange(nvalues):
                # Prepare structure
                val = np.zeros(nteff*nlogg)
                val = val.reshape(nteff, nlogg)
                # Fill values for every teff-logg combination
                for j in np.arange(nlogg):
                    logg = logg_range[j]
                    last_valid_pos = -1
                    for i in np.arange(nteff):
                        teff = teff_range[i]
                        teff_index = np.where(teff_range==teff)[0][0]
                        logg_index = np.where(logg_range==logg)[0][0]
                        pos = model_combinations[metal_num][teff_index][logg_index]
                        if pos != -1:
                            # For this teff-logg combination, there is an atmosphere model
                            val[i][j] = atmospheres_with_same_metallicity[pos][layer_num][value_num]
                            last_valid_pos = pos
                        else:
                            if last_valid_pos == -1:
                                # There should be at least one model for each logg of the first temperature
                                raise Exception("Bad model format")
                            else:
                                # Since there is no atmosphere model for this teff-logg combination
                                # we replicate the last one used with the same logg but lower temperature
                                # we do this to minimize strange effects in the posterior cubic spline interpolation
                                val[i][j] = atmospheres_with_same_metallicity[last_valid_pos][layer_num][value_num]
                if value_num == 6: # microturbulent_vel value
                    # Microturbulence velocity is constant for all layers of the castelli models
                    model_val = ConstantValue(val[0][0])
                else:
                    ### Modelize for a finer grid
                    model_val = interpolate.RectBivariateSpline(grid_teff[:,0], grid_logg[1], val)
                # Add to current layer values
                modeled_values.append([val, model_val])
            modeled_layers_with_same_metallicity.append(modeled_values)
        modeled_layers.append(modeled_layers_with_same_metallicity)

    return modeled_layers


# Plots all the values (except microturbulence vel which is constant) for 3 layers (0, 35 and 70) for
# all the different original metallicities:
#  - Values for the original teff-logg grid
#  - Interpolated values in a finer teff-logg grid
def plot_original_and_finer_grid(modeled_layers, teff_range, logg_range, teff_min, teff_max, logg_min, logg_max, nMH, nlayers, nvalues, nteff, nlogg, teff_step=100, logg_step=0.10):
    ### Plotting
    plot_output_dir = "output/plots/"
    
    ### Construct original grid
    grid_teff, grid_logg = np.indices((nteff, nlogg), dtype=float)
    # Assign values to each grid position
    for i in np.arange(nteff):
        grid_teff[i] = teff_range[i]
    for i in np.arange(nlogg):
        grid_logg[:,i] = logg_range[i]
    

    # Construct a finer grid for plotting
    new_grid_teff, new_grid_logg = np.mgrid[teff_min:teff_max:teff_step,logg_min:logg_max:logg_step]
    new_nteff = len(new_grid_teff[:,0])
    new_nlogg = len(new_grid_logg[:,0])
    
    
    for metal_num in np.arange(nMH):
        print "\nAtmosphere models", metal_num
        for layer_num in np.arange(0, nlayers, 35): # Only 3 layers (0, 35 and 70)
            print "\t Layer:", layer_num
            print "\t\t Value:",
            for value_num in np.arange(nvalues-1): # -1: do not plot last value because it is always constant (microturbulence velocity)
                print value_num,
                sys.stdout.flush()
                model_val = modeled_layers[metal_num][layer_num][value_num][1]
                val = modeled_layers[metal_num][layer_num][value_num][0]
                refined_val = model_val(new_grid_teff[:,0], new_grid_logg[0,:])
                ## Draw
                ##plt.figure(figsize=(9,6))
                plt.subplot(211)
                plt.pcolor(grid_teff, grid_logg, val)
                plt.xlabel("Temperature (K)")
                plt.ylabel("Log g (dex)")
                plt.title("Original")
                
                plt.subplot(212)
                plt.pcolor(new_grid_teff, new_grid_logg, refined_val)
                plt.xlabel("Temperature (K)")
                plt.ylabel("Log g (dex)")
                plt.title("Interpolated")
                
                plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
                cax = plt.axes([0.85, 0.1, 0.075, 0.8])
                plt.colorbar(cax=cax)
                plt.savefig(plot_output_dir + "metal_" + str(metal_num) + "_value_" + str(value_num) + "_layer_" + str(layer_num) + ".png")

# Checks if the objectif teff, logg and metallicity can be obtained by using the loaded model
#  - Validates limits http://wwwuser.oat.ts.astro.it/castelli/grids/gridp00k2odfnew/ap00k2tab.html
def valid_objective(modeled_layers_pack, teff_obj, logg_obj, MH_obj):
    modeled_layers, model_combinations, atmosphere_models, atmosphere_models_metalicities, teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max = modeled_layers_pack
    
    teff_index = np.searchsorted(teff_range, teff_obj)
    if teff_index == 0 and teff_obj != teff_range[0]:
        #raise Exception("Out of range: low teff value")
        return False
    if teff_index >= nteff:
        #raise Exception("Out of range: high teff value")
        return False
    
    logg_index = np.searchsorted(logg_range, logg_obj)
    if logg_index == 0 and logg_obj != logg_range[0]:
        #raise Exception("Out of range: low logg value")
        return False
    if logg_index >= nlogg:
        #raise Exception("Out of range: high logg value")
        return False
    
    valid = True
    for metal_num in np.arange(nMH):
        valid = valid & (model_combinations[metal_num][teff_index][logg_index] != -1)
        if teff_obj != teff_range[logg_index]:
            # teff_obj is between teff_index-1 and teff_index
            valid = valid & (model_combinations[metal_num][teff_index-1][logg_index] != -1)
        if logg_obj != logg_range[logg_index]:
            # logg_obj is between logg_index-1 and logg_index
            valid = valid & (model_combinations[metal_num][teff_index][logg_index-1] != -1)
        if logg_obj != logg_range[logg_index] and teff_obj != teff_range[logg_index]:
            # One additional check
            valid = valid & (model_combinations[metal_num][teff_index-1][logg_index-1] != -1)
    
    return valid


# Generates an interpolated atmosphere for a given teff, logg and metallicity
def interpolate_atmosphere_layers(modeled_layers_pack,  teff_obj, logg_obj, MH_obj):
    modeled_layers, model_combinations, atmosphere_models, atmosphere_models_metalicities, teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max = modeled_layers_pack
    
    ## For each layer
    layers = []
    for layer_num in np.arange(nlayers):
        values = []
        for value_num in np.arange(nvalues):
            base_values = []
            # For each metallicity
            for metal_num in np.arange(nMH):
                # Value for a given teff and logg
                model_val = modeled_layers[metal_num][layer_num][value_num][1]
                val = model_val(teff_obj, logg_obj)
                base_values.append(val)
            # Modelize for obtaining the interpolated value for a given metallicity
            model_MH_val = interpolate.UnivariateSpline(atmosphere_models_metalicities, base_values)
            val = model_MH_val(MH_obj)
            values.append(val)
        layers.append(values)
    return layers

def write_atmosphere(teff, logg, MH, layers):
    atm_file = tempfile.NamedTemporaryFile(delete=False)
    atm_file.write("%.1f  %.5f  %.2f  %i\n" % (teff, logg, MH, len(layers)) )
    for layer in layers:
        atm_file.write("%.8e   %.1f %.3e %.3e %.3e %.3e %.3e\n" % (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], layer[6]) )
    atm_file.close()
    return atm_file.name


# Serialize modeled layers and stats
def save_modeled_layers_pack(modeled_layers, model_combinations, atmosphere_models, atmosphere_models_metalicities, teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max, filename='input/atmospheres/modeled_layers_pack.dump'):
    modeled_layers_pack = (modeled_layers, model_combinations, atmosphere_models, atmosphere_models_metalicities, teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max)
    pickle.dump(modeled_layers_pack, open(filename, 'w'))

# Restore modeled layers and stats
def load_modeled_layers_pack(filename='input/atmospheres/modeled_layers_pack.dump'):
    sys.modules['__main__'].ConstantValue = ConstantValue
    modeled_layers_pack = pickle.load(open(filename))
    return modeled_layers_pack

if __name__ == '__main__':
    # Model file and metallicity
    #atmosphere_models = ["input/atmospheres/am25k2odfnew.dat", "input/atmospheres/am20k2odfnew.dat", "input/atmospheres/am15k2odfnew.dat", "input/atmospheres/am10k2odfnew.dat", "input/atmospheres/am05k2odfnew.dat", "input/atmospheres/ap00k2odfnew.dat", "input/atmospheres/ap02k2odfnew.dat", "input/atmospheres/ap05k2odfnew.dat"]
    #atmosphere_models_metalicities = [0.5, 0.2, 0.0, -0.5, -1.0, -1.5, -2.0, -2.5]
    #atmospheres, atmospheres_params = read_castelli_kurucz_atmospheres(atmosphere_models, atmosphere_models_metalicities)

    ### Stats
    #teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max = atmospheres_statistics(atmospheres_params)

    ## Check teff-logg combinations and register its position in atmospheres array (-1 if it does not exist)
    #model_combinations = check_teff_logg_combinations(atmospheres_params, nteff, nlogg, nMH)

    #modeled_layers = modeled_interpolated_layer_values(atmospheres, model_combinations, nteff, nlogg, nMH, nlayers, nvalues, teff_range, logg_range)

    #plot_original_and_finer_grid(modeled_layers, teff_range, logg_range, teff_min, teff_max, logg_min, logg_max, nMH, nlayers, nvalues, nteff, nlogg, teff_step=100, logg_step=0.10)

    ### Serialize
    #save_modeled_layers_pack(modeled_layers, model_combinations, atmosphere_models, atmosphere_models_metalicities, teff_range, logg_range, MH_range, nteff, nlogg, nMH, nlayers, nvalues, teff_min, teff_max, logg_min, logg_max, MH_min, MH_max)
    modeled_layers_pack = load_modeled_layers_pack(filename='input/atmospheres/modeled_layers_pack.dump')

    teff_obj = 5500
    logg_obj = 0.0
    MH_obj = 0.0

    valid_objective(modeled_layers_pack, teff_obj, logg_obj, MH_obj)
    layers = interpolate_atmosphere_layers(modeled_layers_pack, teff_obj, logg_obj, MH_obj)
    atm_filename = write_atmosphere(teff_obj, logg_obj, MH_obj, layers)

    #os.remove(atm_filename)
