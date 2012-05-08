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
import ipdb

# SPECTRUM is compatible only with the plane-parallel atmospheres.
# The first layer represents the surface. 
# Elemental abundances in the stdatom.dat file are used (and scaled with the [M/H] value)

# Constant class used for microturbulent velocities because they are constant for all layers and atmospheres
class ConstantValue:
    def __init__(self, value):
        self.value = value
    
    def __call__(self, x, y):
        return [[self.value]]

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
def read_stats_kurucz_atmospheres(atmosphere_models, required_layers=72):
    min_values_per_layer = 7
    read_atmosphere_data = False
    atmospheres_params = []
    atmospheres = []
    
    metallicities = []
    temperatures = []
    gravities = []
    for i in np.arange(len(atmosphere_models)):
        atmosphere_model = atmosphere_models[i]
        f = open(atmosphere_model)
        
        metallicity = float(atmosphere_model.split("/")[-1][2] + "." + atmosphere_model.split("/")[-1][3])
        if atmosphere_model.split("/")[-1][1] == "m":
            metallicity *= -1
        metallicities.append(metallicity)
        
        atmospheres_with_same_metallicity = []
        atmospheres_params_with_same_metallicity = []
        for line in f.readlines():
            vline = line.split()
            if vline[0] == "TEFF":
                teff = np.round(float(vline[1]), 2) # K
                logg = np.round(float(vline[3]), 2)
            elif vline[0] == "READ":
                # New atmosphere
                current_atmosphere = []
                num_layers = 0
                read_atmosphere_data = True
            elif read_atmosphere_data:
                if vline[0] == "PRADK":
                    # Only consider atmospheres with the required number of layers
                    if num_layers == required_layers:
                        temperatures.append(teff)
                        gravities.append(logg)
                        atmospheres_params_with_same_metallicity.append([teff, logg, metallicity])
                        atmospheres_with_same_metallicity.append(current_atmosphere)
                    read_atmosphere_data = False
                else:
                    num_layers += 1
                    layer = line.split(" ")
                    # Clean empty values due to consecutive spaces
                    while True:
                        try:
                            layer.remove('')
                        except ValueError:
                            break
                    
                    if len(layer) < min_values_per_layer:
                        raise Exception("FORMAT ERROR: Not enough values")
                    
                    # Only use the 7 first values
                    current_atmosphere.append(layer[0:7])
        atmospheres.append(atmospheres_with_same_metallicity)
        atmospheres_params.append(atmospheres_params_with_same_metallicity)
        f.close()
    
    MH_range = np.unique(metallicities)
    teff_range = np.unique(temperatures)
    logg_range = np.unique(gravities)
    
    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)
    
    
    # Build an array with as many matrices as metallicities
    model_combinations = []
    for metal_num in np.arange(nMH):
        # Build a matrix Teff x Log(g) filled with 'None'
        model_combination_matrix = np.array([None]*nteff*nlogg)
        model_combination_matrix = model_combination_matrix.reshape(nteff, nlogg)
        model_combinations.append(model_combination_matrix)

    ## Fill the array of matrices with the corresponding atmospheric layers
    # For each group of athmospheres with the same metallicity
    for metal_num in np.arange(nMH):
        atmospheres_with_same_metallicity = atmospheres[metal_num]
        atmospheres_params_with_same_metallicity = atmospheres_params[metal_num]
        # For each atmosphere
        for i in np.arange(len(atmospheres_params_with_same_metallicity)):
            # Current atmosphere and its parameters
            atm = atmospheres_with_same_metallicity[i]
            atm_teff = atmospheres_params_with_same_metallicity[i][0]
            atm_logg = atmospheres_params_with_same_metallicity[i][1]
            # Assign the atmosphere to the correct position inside the array of matrices
            teff_index = np.where(teff_range==atm_teff)[0][0]
            logg_index = np.where(logg_range==atm_logg)[0][0]
            model_combinations[metal_num][teff_index][logg_index] = atm
    
    # Output:
    #  - model_combinations (as many as nMH)
    #     - model combination matrix
    #       -> Atmosphere  : It exists
    #       -> None
    return model_combinations, teff_range, logg_range, MH_range
 

# Builds an structure where each value of each layer has a RectBivariateSpline (based on the values
# read from atmospheric models) that can be used for interpolation. The structure is:
# - Modeled layers (as many as different metallicities)
#    - Layers
#       - Modeled values
#           -> Model for interpolation
# - Values used for modeled layers (as many as different metallicities) => Only useful for posterior plotting
#    - Layers
#       - Used values
#           -> Matrix value (for each teff-logg)
def modeled_interpolated_layer_values(model_combinations, teff_range, logg_range, MH_range, required_layers=72):
    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)
    nlayers = required_layers
    nvalues = 7 # Only use the 7 first values

    ### Construct base grid for interpolation
    grid_teff, grid_logg = np.indices((nteff, nlogg), dtype=float)
    # Assign values to each grid position
    for i in np.arange(nteff):
        grid_teff[i] = teff_range[i]
    for i in np.arange(nlogg):
        grid_logg[:,i] = logg_range[i]
    
    modeled_layers = []
    used_values_for_layers = [] # Useful only for plotting
    # For each atmosphere with the same metallicity
    for metal_num in np.arange(nMH):
        modeled_layers_with_same_metallicity = []
        used_values_for_layers_with_same_metallicity = []  # Useful only for plotting
        print "\nAtmosphere models", metal_num
        print "\tLayer:",
        # For each layer, group and modelize the different atmospheres (teff, logg)
        for layer_num in np.arange(nlayers):
            modeled_values = []
            used_values = []  # Useful only for plotting
            print layer_num,
            sys.stdout.flush()
            for value_num in np.arange(nvalues):
                # Prepare structure
                val = np.zeros(nteff*nlogg)
                val = val.reshape(nteff, nlogg)
                # Fill values for every teff-logg combination
                for j in np.arange(nlogg):
                    logg = logg_range[j]
                    last_valid_atm = None
                    for i in np.arange(nteff):
                        teff = teff_range[i]
                        teff_index = np.where(teff_range==teff)[0][0]
                        logg_index = np.where(logg_range==logg)[0][0]
                        atm = model_combinations[metal_num][teff_index][logg_index]
                        if atm != None:
                            # For this teff-logg combination, there is an atmosphere model
                            val[i][j] = atm[layer_num][value_num]
                            last_valid_atm = atm
                        else:
                            if last_valid_atm != None:
                                # Since there is no atmosphere model for this teff-logg combination
                                # we replicate the last one used with the same logg but lower temperature
                                # we do this to minimize strange effects in the posterior cubic spline interpolation
                                val[i][j] = last_valid_atm[layer_num][value_num]
                            else:
                                # Since there is no atmosphere model for this teff-logg combination
                                # and there is not valid alternative with the same logg but lower temperature
                                # we search for the next valid atmosphere with higher temperature
                                # we do this to minimize strange effects in the posterior cubic spline interpolation
                                t = teff_index + 1
                                while t < nteff:
                                    future_valid_atm = model_combinations[metal_num][t][logg_index]
                                    if future_valid_atm != None:
                                        break 
                                    t += 1
                                
                                # If atmosphere found
                                if future_valid_atm != None:
                                    val[i][j] = future_valid_atm[layer_num][value_num]
                                else:
                                    # There should be at least one valid model for each logg
                                    raise Exception("Bad model format")
                if value_num == 6: # microturbulent_vel value
                    # Microturbulence velocity is constant for all layers of the castelli models
                    model_val = ConstantValue(val[0][0])
                else:
                    ### Modelize for a finer grid
                    model_val = interpolate.RectBivariateSpline(grid_teff[:,0], grid_logg[0], val)
                # Add to current layer values
                modeled_values.append(model_val)
                used_values.append(val)
            modeled_layers_with_same_metallicity.append(modeled_values)
            used_values_for_layers_with_same_metallicity.append(used_values)
        modeled_layers.append(modeled_layers_with_same_metallicity)
        used_values_for_layers.append(modeled_layers_with_same_metallicity)
    return modeled_layers, used_values_for_layers


# Plots all the values (except microturbulence vel which is constant) for 3 layers (0, 35 and 70) for
# all the different original metallicities:
#  - Values for the original teff-logg grid
#  - Interpolated values in a finer teff-logg grid
def plot_original_and_finer_grid(modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range, required_layers=72, teff_step=100, logg_step=0.10):
    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)
    nlayers = required_layers
    nvalues = min_values_per_layer # Only use the 7 first values

    teff_min = np.min(teff_range)
    teff_max = np.max(teff_range)
    logg_min = np.min(logg_range)
    logg_max = np.max(logg_range)
    MH_min = np.min(MH_range)
    MH_max = np.max(MH_range)
    
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
                model_val = modeled_layers[metal_num][layer_num][value_num]
                val = used_values_for_layers[metal_num][layer_num][value_num]
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
#  - Validates limits (i.e http://wwwuser.oat.ts.astro.it/castelli/grids/gridp00k2odfnew/ap00k2tab.html)
def valid_objective(modeled_layers_pack, teff_obj, logg_obj, MH_obj):
    modeled_layers, model_combinations, teff_range, logg_range, MH_range, nlayers = modeled_layers_pack
    
    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)
    
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
    
    MH_index = np.searchsorted(MH_range, MH_obj)
    if MH_index == 0 and MH_obj != MH_range[0]:
        #raise Exception("Out of range: low MH value")
        return False
    if MH_index >= nMH:
        #raise Exception("Out of range: high MH value")
        return False
    
    valid = True
    for metal_num in [MH_index-1, MH_index]:
        valid = valid & (model_combinations[metal_num][teff_index][logg_index] != None)
        if teff_obj != teff_range[logg_index]:
            # teff_obj is between teff_index-1 and teff_index
            valid = valid & (model_combinations[metal_num][teff_index-1][logg_index] != None)
        if logg_obj != logg_range[logg_index]:
            # logg_obj is between logg_index-1 and logg_index
            valid = valid & (model_combinations[metal_num][teff_index][logg_index-1] != None)
        if logg_obj != logg_range[logg_index] and teff_obj != teff_range[logg_index]:
            # One additional check
            valid = valid & (model_combinations[metal_num][teff_index-1][logg_index-1] != None)
    
    return valid

# Generates an interpolated atmosphere for a given teff, logg and metallicity
def interpolate_atmosphere_layers(modeled_layers_pack,  teff_obj, logg_obj, MH_obj):
    modeled_layers, model_combinations, teff_range, logg_range, MH_range, nlayers = modeled_layers_pack
    
    nMH  = len(MH_range)
    nvalues = 7
    MH_index = MH_range.searchsorted(MH_obj)
    if MH_index == 0 and MH_obj != MH_range[0]:
        raise Exception("Out of range: low MH value")
    if MH_index >= nMH:
        raise Exception("Out of range: high MH value")
    
    ## For each layer
    layers = []
    for layer_num in np.arange(nlayers):
        values = []
        
        for value_num in np.arange(nvalues):
            if MH_obj == MH_range[MH_index]:
                model_val = modeled_layers[MH_index][layer_num][value_num]
                val = model_val(teff_obj, logg_obj)[0][0]
            else:
                # In between two known metallicities
                MH_xcoord = []
                MH_xcoord.append(MH_range[MH_index-1])
                MH_xcoord.append(MH_range[MH_index])
                
                MH_ycoord = []
                # Value for a given teff and logg
                model_val = modeled_layers[MH_index-1][layer_num][value_num]
                val = model_val(teff_obj, logg_obj)
                MH_ycoord.append(val[0][0])
                model_val = modeled_layers[MH_index][layer_num][value_num]
                val = model_val(teff_obj, logg_obj)
                MH_ycoord.append(val[0][0])
                
                # Modelize for obtaining the interpolated value for a given metallicity
                model_MH_val = interpolate.interp1d(MH_xcoord, MH_ycoord)
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
def save_modeled_layers_pack(modeled_layers, model_combinations, teff_range, logg_range, MH_range, required_layers=72, filename='input/atmospheres/modeled_layers_pack.dump'):
    nlayers = required_layers
    
    modeled_layers_pack = (modeled_layers, model_combinations, teff_range, logg_range, MH_range, nlayers)
    pickle.dump(modeled_layers_pack, open(filename, 'w'))

# Restore modeled layers and stats
def load_modeled_layers_pack(filename='input/atmospheres/modeled_layers_pack.dump'):
    sys.modules['__main__'].ConstantValue = ConstantValue
    modeled_layers_pack = pickle.load(open(filename))
    return modeled_layers_pack

if __name__ == '__main__':
    # Model file and metallicity
    kurutz = True
    if kurutz:
        ###### Kurucz
        # Files should be listed in increasing order of metallicity
        atmosphere_models = ["input/atmospheres/kurucz/am50k2.dat", 
        "input/atmospheres/kurucz/am45k2.dat", 
        "input/atmospheres/kurucz/am40k2.dat", 
        "input/atmospheres/kurucz/am35k2.dat", 
        "input/atmospheres/kurucz/am30k2.dat", 
        "input/atmospheres/kurucz/am25k2.dat", 
        "input/atmospheres/kurucz/am20k2.dat", 
        "input/atmospheres/kurucz/am15k2.dat", 
        "input/atmospheres/kurucz/am10k2.dat", 
        "input/atmospheres/kurucz/am05k2.dat", 
        "input/atmospheres/kurucz/am03k2.dat", 
        "input/atmospheres/kurucz/am02k2.dat", 
        "input/atmospheres/kurucz/am01k2.dat", 
        "input/atmospheres/kurucz/ap00k2.dat", 
        "input/atmospheres/kurucz/ap01k2.dat", 
        "input/atmospheres/kurucz/ap02k2.dat", 
        "input/atmospheres/kurucz/ap03k2.dat", 
        "input/atmospheres/kurucz/ap05k2.dat", 
        "input/atmospheres/kurucz/ap10k2.dat"]
        dump_filename = 'input/atmospheres/kurucz/modeled_layers_pack.dump'
    else:
        ###### Castelli atmospheres
        # Files should be listed in increasing order of metallicity
        atmosphere_models = ["input/atmospheres/castelli/am25k2odfnew.dat", "input/atmospheres/castelli/am20k2odfnew.dat", "input/atmospheres/castelli/am15k2odfnew.dat", "input/atmospheres/castelli/am10k2odfnew.dat", "input/atmospheres/castelli/am05k2odfnew.dat", "input/atmospheres/castelli/ap00k2odfnew.dat", "input/atmospheres/castelli/ap02k2odfnew.dat", "input/atmospheres/castelli/ap05k2odfnew.dat"]
        dump_filename = 'input/atmospheres/castelli/modeled_layers_pack.dump'
    

    model_combinations, teff_range, logg_range, MH_range = read_stats_kurucz_atmospheres(atmosphere_models, required_layers=72)

    modeled_layers, used_values_for_layers = modeled_interpolated_layer_values(model_combinations, teff_range, logg_range, MH_range)
    
    #plot_original_and_finer_grid(modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range, teff_step=100, logg_step=0.10)
    
    #### Serialize
    save_modeled_layers_pack(modeled_layers, model_combinations, teff_range, logg_range, MH_range, required_layers=72, filename=dump_filename)
    
    modeled_layers_pack = load_modeled_layers_pack(filename=dump_filename)
    
    #modeled_layers_pack = (modeled_layers, model_combinations, teff_range, logg_range, MH_range, 72)
    
    teff_obj = 5500
    logg_obj = 0.0
    MH_obj = 0.0

    valid_objective(modeled_layers_pack, teff_obj, logg_obj, MH_obj)
    layers = interpolate_atmosphere_layers(modeled_layers_pack, teff_obj, logg_obj, MH_obj)
    atm_filename = write_atmosphere(teff_obj, logg_obj, MH_obj, layers)

    os.remove(atm_filename)
    
