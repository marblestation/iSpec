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
import os
import sys
import numpy as np
from subprocess import Popen, PIPE
import asciitable
import math
from datetime import datetime
from scipy import interpolate
import tempfile
import cPickle as pickle

# SPECTRUM is compatible only with the plane-parallel atmospheres.
# The first layer represents the surface.
# Elemental abundances in the stdatom.dat file are used (and scaled with the [M/H] value)

class ConstantValue:
    """ Constant class used for microturbulent velocities because they are
        constant for all layers and atmospheres """
    def __init__(self, value):
        self.value = value

    def __call__(self, x, y):
        return self.value

def read_kurucz_atmospheres(atmosphere_models, required_layers=56):
    """
    Read castelli and kurucz atmospheres.

    :param atmosphere_models:
        List or array of files with Kurucz atmospheres ordered from lower to
        higher metallicity. Example: atmosphere_models = ["input/atmospheres/kurucz/am50k2.dat",]
    :type atmosphere_models: array

    :returns:
        Atmosphere's parameters and values

        teff_range, logg_range, MH_range are arrays with the list of effective
        temperature, gravity and metallicity

    """
    min_values_per_layer = 7
    read_atmosphere_data = False
    atmospheres_params = []
    atmospheres = []

    # Build the following structures:
     #atmospheres_params (as many as different metallicities):
      #- atmospheres_params_with_same_metallicity
          #- atmospheres
              #-> params
     #atmospheres (as many as different metallicities):
      #- atmospheres_with_same_metallicity
          #- atmospheres
              #-> values
                  #* Values' order:
                      #mass_depth = 0
                      #temperature_K = 1
                      #gas_preassure = 2
                      #electron_density = 3
                      #mean_absorption_coeff = 4
                      #radiation_preassure = 5
                      #microturbulent_vel = 6 # In m/s

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
                    current_atmosphere.append(map(float, layer[0:7]))
        atmospheres.append(atmospheres_with_same_metallicity)
        atmospheres_params.append(atmospheres_params_with_same_metallicity)
        f.close()

    MH_range = np.unique(metallicities)
    teff_range = np.unique(temperatures)
    logg_range = np.unique(gravities)

    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)


    return atmospheres_params, atmospheres, teff_range, logg_range, MH_range

def __extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y = np.where(x<xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x>xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
    return y

def build_modeled_interpolated_layer_values(atmospheres_params, atmospheres, teff_range, logg_range, MH_range, required_layers=56):
    """
    Builds an structure where each value of each layer has a RectBivariateSpline (based on the values
    read from atmospheric models) that can be used for interpolation.

    :param atmospheres:
        Output from read_kurucz_atmospheres method
    :type atmospheres: array

    :param teff_range:
        Output from read_kurucz_atmospheres method
    :type teff_range: array

    :param logg_range:
        Output from read_kurucz_atmospheres method
    :type teff_range: array

    :param MH_range:
        Output from read_kurucz_atmospheres method
    :type teff_range: array

    :returns:
        modeled_layers is an array with as many elements as different metallicities

        Layers - Modeled values - Model for interpolation

        used_values_for_layers is an array with as many elements as different metallicities
        which is basicly useful only for plotting the values used for building the models:

        Layers - Used values - Matrix value (for each teff-logg)
    """
    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)
    nlayers = required_layers
    nvalues = 7 # Only use the 7 first values

    models_atm_same_metallicity = []
    values_atm_same_metallicity = [] # Useful only for plotting
    structured_values_atm_same_metallicity = [] # Useful only for plotting
    params_atm_same_metallicity = []

    print "\n1) Searching for minimum/maximum values in each atmosphere/layer"
    print   "----------------------------------------------------------------"
    min_value = np.array([np.inf]*(nvalues))
    max_value = np.array([-np.inf]*(nvalues))
    for metal_num in np.arange(nMH):
        print "\nAtmosphere models", metal_num
        print "\tLayer:",
        for layer_num in np.arange(nlayers):
            print layer_num,
            sys.stdout.flush()
            for value_num in np.arange(nvalues):
                for atm_num in xrange(len(atmospheres_params[metal_num])):
                    single_value = float(atmospheres[metal_num][atm_num][layer_num][value_num])
                    max_value[value_num] = np.max([single_value, max_value[value_num]])
                    min_value[value_num] = np.min([single_value, min_value[value_num]])


    print "\n\n2) Building models for interpolation/extrapolation"
    print   "----------------------------------------------------------------"
    # For each atmosphere with the same metallicity
    for metal_num in np.arange(nMH):
        print "\nAtmosphere models", metal_num
        print "\tLayer:",

        models_atm = []
        values_atm = []  # Useful only for plotting
        structured_values_atm = []
        params_atm = []
        atm_teff = []
        atm_logg = []
        # For this metalicity what teff-logg combinations do we have:
        for atm_num in xrange(len(atmospheres_params[metal_num])):
            atm_teff.append(atmospheres_params[metal_num][atm_num][0])
            atm_logg.append(atmospheres_params[metal_num][atm_num][1])
            #metallicity = atmospheres_params[metal_num][i][3]
        params_atm.append((atm_teff, atm_logg))

        # For each layer, group and modelize the different atmospheres (teff, logg)
        for layer_num in np.arange(nlayers):
            print layer_num,
            sys.stdout.flush()

            models_single_layer = []
            values_single_layer = [] # Useful only for plotting
            structured_values_single_layer = [] # Useful only for plotting
            for value_num in np.arange(nvalues):
                # Prepare structure
                total = nteff*nlogg
                structured_single_values = np.ma.array(np.array([np.nan]*total), mask=[i<0 for i in xrange(total)])
                structured_single_values = structured_single_values.reshape(nteff, nlogg)

                # Recollect values and save them in a list and a structured container
                single_values = []
                for atm_num in xrange(len(atmospheres_params[metal_num])):
                    single_value = float(atmospheres[metal_num][atm_num][layer_num][value_num])
                    single_values.append(single_value)
                    # Save
                    teff = atmospheres_params[metal_num][atm_num][0]
                    logg = atmospheres_params[metal_num][atm_num][1]
                    teff_index = np.where(teff_range==teff)[0][0]
                    logg_index = np.where(logg_range==logg)[0][0]
                    structured_single_values[teff_index][logg_index] = single_value


                if value_num == 6: # microturbulent_vel value
                    # Microturbulence velocity is constant for all layers
                    single_model = ConstantValue(single_values[0])
                else:
                    # For those teff-logg combination that we do not have an atmosphere:
                    # - Linearly interpolate each value from the existing atmospheres:
                    #   a. Interpolate using values for all the available temperatures with the fixed current logg
                    #   b. Interpolate using values for all the available logg with the fixed current temperature
                    #   c. If interpolation has been possible in (a) and (b):
                    #      - Average the two results
                    #      If interpolation has been possible in (a) or (b):
                    #      - Use the interpolated result from (a) or (b)
                    #      If no interpolation has been possible, **linear extrapolation** is needed:
                    #         1. Extrapolate using the two nearest real values in the logg axis
                    #         2. Extrapolate using the two nearest real value in the temperature axis
                    #         3. If a extrapolated value has been derived in (1) and (2):
                    #            - Do a weighted average of the two results by giving more
                    #              weight to the value that has used the nearest real values
                    #            If a value has been found in (1) or (2), go to second extrapolation process
                    # - Second extrapolation process:
                    #   a. Extrapolate using the two nearest real values in the logg axis
                    #   b. Extrapolate using the two nearest real value in the temperature axis
                    #   c. Use the extrapolated result from (a) or (b)
                    #      - It is expected to have at least two atmospheres with
                    #        the same logg or teff in the grid of combinations.
                    #        Therefore, extrapolation (a) or (b) should be always possible.
                    #        If not, an exception (error message) will occur.
                    #
                    # NOTE: In all cases, to avoid unphysical results, extrapolated
                    #       values are limited by the maximum and minimum values
                    #       found in all the real atmospheres.
                    structured_single_values_interpolated = structured_single_values.copy() # Original values + interpolated/duplicated ones
                    for logg_index in np.arange(nlogg):
                        logg = logg_range[logg_index]
                        for teff_index in np.arange(nteff):
                            teff = teff_range[teff_index]
                            if np.isnan(structured_single_values[teff_index][logg_index]):
                                # Mark this value to keep track of those that are not original
                                structured_single_values.mask[teff_index][logg_index] = True
                                structured_single_values_interpolated.mask[teff_index][logg_index] = True
                                # For interpolation/extrapolation, ignore the np.nan values:
                                valid_logg = np.where(np.logical_not(np.isnan(structured_single_values.data[teff_index])))
                                valid_teff = np.where(np.logical_not(np.isnan(structured_single_values.data.T[logg_index])))

                                # Interpolate/extrapolate by using logg values for a fixed teff
                                interp_logg_value = np.nan
                                extrap_logg_value = np.nan
                                if len(logg_range[valid_logg]) >= 2:
                                    if logg_index < valid_logg[0][0]:
                                        # Extrapolate and do not allow to go beyong the max/min values found in the real atmospheres
                                        extrap_logg_value = __extrap(logg_range[logg_index], logg_range[valid_logg][:2], structured_single_values.data[teff_index][valid_logg][:2])
                                        extrap_logg_value = np.min((max_value[value_num], extrap_logg_value))
                                        extrap_logg_value = np.max((min_value[value_num], extrap_logg_value))
                                        logg_limit_index = valid_logg[0][0]
                                    elif logg_index > valid_logg[0][-1]:
                                        # Extrapolate and do not allow to go beyong the max/min values found in the real atmospheres
                                        extrap_logg_value = __extrap(logg_range[logg_index], logg_range[valid_logg][len(logg_range[valid_logg])-2:], structured_single_values.data[teff_index][valid_logg][len(logg_range[valid_logg])-2:])
                                        extrap_logg_value = np.min((max_value[value_num], extrap_logg_value))
                                        extrap_logg_value = np.max((min_value[value_num], extrap_logg_value))
                                        logg_limit_index = valid_logg[0][-1]
                                    else:
                                        # Interpolate
                                        interp_logg_value = np.interp(logg_range[logg_index], logg_range[valid_logg], structured_single_values.data[teff_index][valid_logg], left=np.nan, right=np.nan)

                                # Interpolate/extrapolate by using teff values for a fixed logg
                                interp_teff_value = np.nan
                                extrap_teff_value = np.nan
                                if len(teff_range[valid_teff]) >= 2:
                                    if teff_index < valid_teff[0][0]:
                                        # Extrapolate and do not allow to go beyong the max/min values found in the real atmospheres
                                        extrap_teff_value = __extrap(teff_range[teff_index], teff_range[valid_teff][:2], structured_single_values.data.T[logg_index][valid_teff][:2])
                                        extrap_teff_value = np.min((max_value[value_num], extrap_teff_value))
                                        extrap_teff_value = np.max((min_value[value_num], extrap_teff_value))
                                        teff_limit_index = valid_teff[0][0]
                                        interp_teff_value = np.nan
                                    elif teff_index > valid_teff[0][-1]:
                                        # Extrapolate and do not allow to go beyong the max/min values found in the real atmospheres
                                        extrap_teff_value = __extrap(teff_range[teff_index], teff_range[valid_teff][len(teff_range[valid_teff])-2:], structured_single_values.data.T[logg_index][valid_teff][len(teff_range[valid_teff])-2:])
                                        extrap_teff_value = np.min((max_value[value_num], extrap_teff_value))
                                        extrap_teff_value = np.max((min_value[value_num], extrap_teff_value))
                                        teff_limit_index = valid_teff[0][-1]
                                        interp_teff_value = np.nan
                                    else:
                                        # Interpolate
                                        interp_teff_value = np.interp(teff_range[teff_index], teff_range[valid_teff], structured_single_values.data.T[logg_index][valid_teff], left=np.nan, right=np.nan)

                                # Interpolated values have priority over extrapolated ones
                                if not np.isnan(interp_logg_value) and not np.isnan(interp_teff_value):
                                    derived_value = (interp_logg_value + interp_teff_value) / 2.0
                                elif not np.isnan(interp_logg_value):
                                    derived_value = interp_logg_value
                                elif not np.isnan(interp_teff_value):
                                    derived_value = interp_teff_value
                                else:
                                    if not np.isnan(extrap_logg_value) and not np.isnan(extrap_teff_value):
                                        # The final extrapolated value is a weighted average:
                                        # - The value extrapolated with closer real values has more weight
                                        # WARNING: It is considered that each step in logg/teff are approximatelly of the same order of magnitude
                                        logg_jumps = float(np.abs(logg_index - logg_limit_index))
                                        teff_jumps = float(np.abs(teff_index - teff_limit_index))
                                        total_jumps = logg_jumps + teff_jumps
                                        derived_value = extrap_logg_value*(1.-(logg_jumps/total_jumps)) + extrap_teff_value*(1.-(teff_jumps/total_jumps))
                                    else:
                                        # Cases delegated to the second extrapolation process:
                                        # 1. Cases where there are only one extrapolated values
                                        # 2. Cases where it has not been possible to extrapolate values due to
                                        #    a lack of atmospheres (there should be at least two atmospheres
                                        #    with the same logg or teff in the grid of combinations)
                                        derived_value = np.nan
                                structured_single_values_interpolated[teff_index][logg_index] = derived_value

                    ############################################################
                    # Second extrapolation process
                    # - In this process, real values + previous averaged extrapoled values will
                    #   be considered to derive values for the cases where only one extrapolation
                    #   can be done (case 1)
                    ############################################################
                    structured_single_values_interpolated2 = structured_single_values_interpolated.copy()
                    for logg_index in np.arange(nlogg):
                        logg = logg_range[logg_index]
                        for teff_index in np.arange(nteff):
                            teff = teff_range[teff_index]
                            if np.isnan(structured_single_values_interpolated.data[teff_index][logg_index]):
                                # Mark this value to keep track of those that are not original
                                structured_single_values_interpolated2.mask[teff_index][logg_index] = True
                                # For extrapolation, ignore the np.nan values:
                                valid_logg = np.where(np.logical_not(np.isnan(structured_single_values_interpolated.data[teff_index])))
                                valid_teff = np.where(np.logical_not(np.isnan(structured_single_values_interpolated.data.T[logg_index])))

                                if len(logg_range[valid_logg]) >= 2:
                                    # Extrapolate by using logg values for a fixed teff
                                    if logg_index < valid_logg[0][0]:
                                        extrap_logg_value = __extrap(logg_range[logg_index], logg_range[valid_logg][:2], structured_single_values_interpolated.data[teff_index][valid_logg][:2])
                                        extrap_logg_value = np.min((max_value[value_num], extrap_logg_value))
                                        extrap_logg_value = np.max((min_value[value_num], extrap_logg_value))
                                    elif logg_index > valid_logg[0][-1]:
                                        extrap_logg_value = __extrap(logg_range[logg_index], logg_range[valid_logg][len(logg_range[valid_logg])-2:], structured_single_values_interpolated.data[teff_index][valid_logg][len(logg_range[valid_logg])-2:])
                                        extrap_logg_value = np.min((max_value[value_num], extrap_logg_value))
                                        extrap_logg_value = np.max((min_value[value_num], extrap_logg_value))
                                    else:
                                        raise Exception("Only values that should be extrapolated are expected at this point")
                                    structured_single_values_interpolated2[teff_index][logg_index] = extrap_logg_value
                                    continue
                                elif len(teff_range[valid_teff]) >= 2:
                                    # Extrapolate by using teff values for a fixed logg
                                    if teff_index < valid_teff[0][0]:
                                        extrap_teff_value = __extrap(teff_range[teff_index], teff_range[valid_teff][:2], structured_single_values_interpolated.data.T[logg_index][valid_teff][:2])
                                        extrap_teff_value = np.min((max_value[value_num], extrap_teff_value))
                                        extrap_teff_value = np.max((min_value[value_num], extrap_teff_value))
                                    elif teff_index > valid_teff[0][-1]:
                                        extrap_teff_value = __extrap(teff_range[teff_index], teff_range[valid_teff][len(teff_range[valid_teff])-2:], structured_single_values_interpolated.data.T[logg_index][valid_teff][len(teff_range[valid_teff])-2:])
                                        extrap_teff_value = np.min((max_value[value_num], extrap_teff_value))
                                        extrap_teff_value = np.max((min_value[value_num], extrap_teff_value))
                                    else:
                                        raise Exception("Only values that should be extrapolated are expected at this point")
                                    structured_single_values_interpolated2[teff_index][logg_index] = extrap_teff_value
                                    continue
                                else:
                                    raise Exception("There should be at least two atmospheres with the same logg or teff in the grid of combinations")

                    ############################################################
                    # Build a spline model to be able to estimate a value for a given teff/logg combination
                    single_model = interpolate.RectBivariateSpline(teff_range, logg_range, structured_single_values_interpolated2.data)
                # Add to current model and used values
                models_single_layer.append(single_model)
                values_single_layer.append(single_values)
                structured_values_single_layer.append(structured_single_values_interpolated2)
            # Save models and used values for the whole layer
            models_atm.append(models_single_layer)
            values_atm.append(values_single_layer)
            structured_values_atm.append(structured_values_single_layer)

        models_atm_same_metallicity.append(models_atm)
        values_atm_same_metallicity.append(values_atm)
        structured_values_atm_same_metallicity.append(structured_values_atm)
        params_atm_same_metallicity.append(params_atm)

    return models_atm_same_metallicity, structured_values_atm_same_metallicity


def valid_atmosphere_target(modeled_layers_pack, teff_target, logg_target, MH_target):
    """
    Checks if the objectif teff, logg and metallicity can be obtained by using the loaded model

    :param modeled_layers_pack:
        Output from load_modeled_layers_pack
    :type modeled_layers_pack: array

    :returns:
        True if the target teff, logg and metallicity can be obtained with the
        models
    """
    modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range, nlayers = modeled_layers_pack

    nteff = len(teff_range)
    nlogg = len(logg_range)
    nMH  = len(MH_range)

    teff_index = np.searchsorted(teff_range, teff_target)
    if teff_index == 0 and teff_target != teff_range[0]:
        #raise Exception("Out of range: low teff value")
        return False
    if teff_index >= nteff:
        #raise Exception("Out of range: high teff value")
        return False

    logg_index = np.searchsorted(logg_range, logg_target)
    if logg_index == 0 and logg_target != logg_range[0]:
        #raise Exception("Out of range: low logg value")
        return False
    if logg_index >= nlogg:
        #raise Exception("Out of range: high logg value")
        return False

    MH_index = np.searchsorted(MH_range, MH_target)
    if MH_index == 0 and MH_target != MH_range[0]:
        #raise Exception("Out of range: low MH value")
        return False
    if MH_index >= nMH:
        #raise Exception("Out of range: high MH value")
        return False

    return True


def interpolate_atmosphere_layers(modeled_layers_pack,  teff_target, logg_target, MH_target):
    """
    Generates an interpolated atmosphere for a given teff, logg and metallicity

    :param modeled_layers_pack:
        Output from load_modeled_layers_pack
    :type modeled_layers_pack: array

    :returns:
        Interpolated model atmosphere
    """
    modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range, nlayers = modeled_layers_pack

    nMH  = len(MH_range)
    nvalues = 7
    MH_index = MH_range.searchsorted(MH_target)
    if MH_index == 0 and MH_target != MH_range[0]:
        raise Exception("Out of range: low MH value")
    if MH_index >= nMH:
        raise Exception("Out of range: high MH value")

    # Prepare structure
    total = nMH*nlayers
    structured_single_values = np.ma.array(np.array([np.nan]*total), mask=[i<0 for i in xrange(total)])
    structured_single_values = structured_single_values.reshape(nMH, nlayers)

    values = []
    for value_num in np.arange(nvalues):
        # Calculate this value for all metallicities and layers
        for metal_num in np.arange(nMH):
            for layer_num in np.arange(nlayers):
                model_val = modeled_layers[metal_num][layer_num][value_num]
                val = model_val(teff_target, logg_target)
                structured_single_values[metal_num][layer_num] = val
        # Fit Bivariate Spline to consider not only variation in metallicity but also in layers
        model_val = interpolate.RectBivariateSpline(MH_range, np.arange(nlayers), structured_single_values.data)
        # Interpolate the current value for all layers
        values.append(model_val(MH_target, np.arange(nlayers))[0])
    # Transpose to have a layer in each position instead than value (easier to write to disk later)
    layers = np.array(values).T

    return layers


def write_atmosphere(teff, logg, MH, layers):
    """
    Write a model atmosphere to a temporary file

    :param layers:
        Output from interpolate_atmosphere_layers
    :type modeled_layers_pack: array

    :returns:
        Name of the temporary file
    """
    atm_file = tempfile.NamedTemporaryFile(delete=False)
    atm_file.write("%.1f  %.5f  %.2f  %i\n" % (teff, logg, MH, len(layers)) )
    for layer in layers:
        atm_file.write("%.8e   %.1f %.3e %.3e %.3e %.3e %.3e\n" % (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], layer[6]) )
    atm_file.close()
    return atm_file.name


# Serialize modeled layers and stats
def dump_modeled_layers_pack(modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range, filename, required_layers=56):
    """
    Build a list of modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range and nlayers
    in order to serialize it to disk for easier later recovery.

    :param filename:
        Name of the output file (i.e. models.dump)
    :type filename: string

    """
    nlayers = required_layers

    modeled_layers_pack = (modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range, nlayers)
    pickle.dump(modeled_layers_pack, open(filename, 'w'))

def load_modeled_layers_pack(filename):
    """
    Restore modeled layers and stats saved previously with save_modeled_layers_pack

    :param filename:
        Name of the input file (i.e. models.dump)
    :type filename: string

    :returns:
        List of modeled_layers, used_values_for_layers, teff_range, logg_range, MH_range and nlayers
    """
    sys.modules['__main__'].ConstantValue = ConstantValue
    modeled_layers_pack = pickle.load(open(filename))
    return modeled_layers_pack

