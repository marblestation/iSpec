#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
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
import os
import sys
import numpy as np
from subprocess import Popen, PIPE
import math
from datetime import datetime
import tempfile
import pickle as pickle
from . import log
import logging
import subprocess
import shutil
import pandas as pd
from astropy.io import fits
from scipy import spatial
from scipy.interpolate import LinearNDInterpolator
import glob
import h5py
import hdf5plugin # required to be able to read hdf5 files compressed with lz4 (NLTE departure coefficient grids), even if it is not explicitly used in the reading code here
from .common import is_turbospectrum_support_enabled, is_spectrum_support_enabled, enhance_solar_abundances
from .synth.effects import _filter_linelist


# SPECTRUM is compatible only with the plane-parallel atmospheres.
# The first layer represents the surface.
# Elemental abundances in the stdatom.dat file are used (and scaled with the [M/H] value)

class ConstantValue(object):
    """ Constant class used for microturbulent velocities because they are
        constant for all layers and atmospheres """
    def __init__(self, value):
        self.value = value

    def __call__(self, x, y):
        return self.value

def valid_atmosphere_target(modeled_layers_pack, target):
    """
    Checks if the objectif teff, logg and metallicity can be obtained by using the loaded model

    :param modeled_layers_pack:
        Output from load_modeled_layers_pack
    :type modeled_layers_pack: array

    :returns:
        True if the target teff, logg and metallicity can be obtained with the
        models
    """
    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, nlte_dep_grid, base_dirname = modeled_layers_pack

    for param in free_parameters:
        param_target = target[param]
        param_range = ranges[param]
        ntotal = len(param_range)

        param_index = np.searchsorted(param_range, param_target)
        if param_index == 0 and param_target != param_range[0]:
            return False
        if param_index >= ntotal:
            return False
    return True

def __closest(kdtree, existing_points, filenames, read_point_value, target_point):
    """
    If there is no model in the extreme points, copy the closest one.
    input:
    - existing_points: points in the parameters space such as [( 3000.,  3.5,  0.  ,  0. ), ( 3000.,  3.5,  0.25,  0. ), ...]
    - target_point: parameters for the target model
    output:
    - created extreme points such as [(3000, 0.0, -5.0, 0.0), (3000, 0.0, -5.0, 0.40000000000000002),]
    """
    distance, index = kdtree.query(target_point, k=1)
    closest_existing_point = existing_points[index]
    closest_existing_point_filename = filenames[index]
    value = read_point_value(closest_existing_point_filename)
    logging.info("Closest to target point '{}' is '{}'".format(" ".join(map(str, target_point)), " ".join(map(str, closest_existing_point))))
    return value

def _interpolate(delaunay_triangulations, kdtree, existing_points, filenames, read_point_value, value_fields, target_point):
    """
    input:
    - existing_points: points in the parameters space such as [( 3000.,  3.5,  0.  ,  0. ), ( 3000.,  3.5,  0.25,  0. ), ...]
    - read_point_value: function to read the value from the model filename
    - target_point: parameters for the target model
    output:
    - interpolated value
    """
    target_point_cannot_be_interpolated = True
    for subset, delaunay_triangulation in zip(delaunay_triangulations['subsets'], delaunay_triangulations['precomputed']):
        if delaunay_triangulation is None:
            continue
        simplex = delaunay_triangulation.find_simplex(target_point)
        if not np.any(simplex == -1):
            target_point_cannot_be_interpolated = False
            break
    if target_point_cannot_be_interpolated:
        logging.warning("Target point '{}' is out of bound, using the closest".format(" ".join(map(str, target_point))))
        return __closest(kdtree, existing_points, filenames, read_point_value, target_point)
    index = delaunay_triangulation.simplices[simplex]
    points = []
    values = {}
    for point, filename in zip(existing_points[subset][index], filenames[subset][index]):
        value = read_point_value(filename)
        points.append(point.tolist())
        for field in value_fields:
            if field not in values:
                values[field] = value[field]
            else:
                values[field] = np.vstack((values[field], value[field]))
    #
    interpolated_values = None
    for field in value_fields:
        interpolator = LinearNDInterpolator(points, values[field])
        interpolated_value = interpolator((target_point))[0]
        if interpolated_values is None:
            interpolated_values = pd.DataFrame(interpolated_value, columns=[field])
        else:
            interpolated_values[field] = interpolated_value
    return interpolated_values.to_records(index=False)



def interpolate_atmosphere_layers(modeled_layers_pack, target, code="spectrum"):
    """
    Generates an interpolated atmosphere for a given teff, logg and metallicity

    :param modeled_layers_pack:
        Output from load_modeled_layers_pack
    :type modeled_layers_pack: array

    :returns:
        Interpolated model atmosphere
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'moog-scat', 'width', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    if code == 'moog-scat':
        # MOOG-SCAT is backward compatible with MOOG
        code = 'moog'

    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, nlte_dep_grid, base_dirname = modeled_layers_pack

    if not valid_atmosphere_target(modeled_layers_pack, target):
        raise Exception("Target parameters '{}' are out of range.".format(target))

    target_point = []
    for param in free_parameters:
        target_point.append(target[param])
    interpolated_atm = _interpolate(delaunay_triangulations, kdtree, existing_points, filenames, read_point_value, value_fields, target_point)

    compatible_fields = ["rhox", "temperature", "pgas", "xne", "abross", "accrad", "vturb", "logtau5", "depth", "pelectron"]
    if "MARCS" in base_dirname:
        compatible_fields.append("radius")
    #interpolated_atm_compatible_format = interpolated_atm[compatible_fields]

    try:
        interpolated_atm_compatible_format = pd.DataFrame(interpolated_atm)[compatible_fields].to_records(index=False)
    except ValueError:
        # Pandas raises exception when interpolated_atm is a FITS table read with astropy (case where there wasn't an interpolation but just the closest point was read)
        #  - ValueError: Big-endian buffer not supported on little-endian compiler
        #
        # The issue: interpolated_atm has big-endian byte order ('>f8' in dtype)
        # but pandas/numpy on this system expects little-endian data
        #
        # newbyteorder('=') converts the data to the native byte order of the current system:
        # - '>' means big-endian (most significant byte first)
        # - '<' means little-endian (least significant byte first)
        # - '=' means native byte order (whatever the current system uses)
        #
        # This conversion ensures the data is in the format expected by pandas
        interpolated_atm_native = interpolated_atm.astype(interpolated_atm.dtype.newbyteorder('='))

        # Now we can safely create a DataFrame and convert to records without endianness errors
        interpolated_atm_compatible_format = pd.DataFrame(interpolated_atm_native)[compatible_fields].to_records(index=False)
    interpolated_atm_compatible_format = interpolated_atm_compatible_format.view(float).reshape(interpolated_atm_compatible_format.shape + (-1,))
    # SME fails if it is not a np.ndarray (it does not accept views either)
    # built like this:
    interpolated_atm_compatible_format_ndarray = np.array(interpolated_atm_compatible_format.T.tolist()).T
    return interpolated_atm_compatible_format_ndarray


def write_atmosphere(atmosphere_layers, teff, logg, MH, atmosphere_filename=None, code='spectrum', tmp_dir=None):
    """
    Write a model atmosphere to file
    If filename is not specified, a temporary file is created and the name is returned.

    :param layers:
        Output from interpolate_atmosphere_layers
    :type modeled_layers_pack: array

    :returns:
        Name of the temporary file
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'moog-scat', 'width', 'synthe']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    if code == 'moog-scat':
        # MOOG-SCAT is backward compatible with MOOG
        code = 'moog'

    if atmosphere_filename is not None:
        atm_file = open(atmosphere_filename, "w")
    else:
        # Temporary file
        atm_file = tempfile.NamedTemporaryFile(mode="wt", delete=False, dir=tmp_dir, encoding='utf-8')

    if code == "moog":
        atm_file.write("KURUCZ\n")
        atm_file.write("TEFF=%i,LOGG=%.2f,[FE/H]=%.2f,marcs\n" % (teff, logg, MH))
        atm_file.write("ND=            %i\n" % (len(atmosphere_layers)))
                                #0.62483359E-03   3636.6 1.937E+01 2.604E+09
        for i, layer in enumerate(atmosphere_layers):
            atm_file.write("%.8E %.1f %.4E %.4E %.4E\n" % (layer[0], layer[1], layer[2], layer[3], layer[4]))
        #vmic = 1.0
        #atm_file.write("  %.2f\n" % (vmic))
        #atm_file.write("NATOMS = 0  %.2f\n" % (MH))
        #atm_file.write("NMOL 0")
    elif code == "turbospectrum":
        is_marcs_model = len(atmosphere_layers[0]) == 11
        if not is_marcs_model:
            atm_file.write("TEFF   %i.  GRAVITY %.5f LTE\n" % (teff, logg))
            atm_file.write("READ DECK6 %i RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB, FLXCNV,VCONV,VELSND\n" % (len(atmosphere_layers)))
            for i, layer in enumerate(atmosphere_layers):
                atm_file.write("%.8E %.1f %.4E %.4E %.4E\n" % (layer[0], layer[1], layer[2], layer[3], layer[4]))
        else:
            for i, layer in enumerate(atmosphere_layers):
                atm_file.write("%.8E %.1f %.4E %.4E %.4E\n" % (layer[0], layer[1], layer[2], layer[3], layer[4]))
            radius = atmosphere_layers[0][-1]
            if radius > 1.0:
                atm_file.write("spherical model\n")
                atm_file.write("  1.0        Mass [Msun]\n")
                atm_file.write("  %.4E Radius [cm] at Tau(Rosseland)=1.0\n" % (radius))
            else:
                atm_file.write("plane-parallel model\n")
                atm_file.write("  0.0        No mass for plane-parallel models\n")
                atm_file.write("  1.0000E+00 1 cm radius for plane-parallel models\n")
            atm_file.write("  %i Number of depth points\n" % (len(atmosphere_layers)))
            atm_file.write("Model structure\n")
            atm_file.write(" k lgTauR  lgTau5    Depth     T        Pe          Pg         Prad       Pturb\n")
            #  temperature structure measured at the continuum optical depth at 5000 A
            #1 -5.00 -4.9174 -6.931E+07  4066.8  2.1166E-02  2.6699E+02  1.4884E+00  0.0000E+00
            ## Gustafsson et al. 2008:
            #  http://marcs.astro.uu.se/GEEJNP08.pdf
            # We have tested the use of Eq. (7) to simulate
            # the effects of turbulent pressure for a number of models at
            # various points in the grid and find that it leads to very small errors
            # in the temperature structure (less than 5 K in the temperature
            # throughout the model for a depth independent vt in the interval 0
            # to 10 km s-1). We, therefore, have chosen to set vt = 0 for all
            # grid models, and advise those who would have liked a different
            # choice to use models with a different mass or g, according to
            # the recipe given in Eq. (7).
            lgTauR = -5.00 # lgTauR is only needed by Turbospectrum to read the model atmosphere but it has no effect
            for i, layer in enumerate(atmosphere_layers):
                atm_file.write("%i %.2f %.4f %.3E %.1f %.4E %.4E %.4E %.4E\n" % (i+1, lgTauR, layer[7], layer[8], layer[1], layer[9], layer[2], layer[5], 0.))
                if lgTauR <= -3.00 or lgTauR >= 1.00:
                    lgTauR += 0.20
                else:
                    lgTauR += 0.10
    elif code == "spectrum":
        # Spectrum
        # mass depth, temperature in kelvin, gas pressure, electron density, Rosseland mean absorption coefficient, radiation pressure, microturbulent velocity in meters/second.
        atm_file.write("%.1f  %.5f  %.2f  %i\n" % (teff, logg, MH, len(atmosphere_layers)))
        #atm_file.write("\n".join(["  ".join(map(str, (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], layer[6]))) for layer in atmosphere_layers]))
        atm_file.write("\n".join(["  ".join(map(str, (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], 1.0))) for layer in atmosphere_layers]))
        #for layer in layers:
            #atm_file.write("%.8e   %.1f %.3e %.3e %.3e %.3e %.3e\n" % (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], layer[6]) )
    elif code == "width" or code == "synthe":
        command_input = ""
        command_input += "READ DECK6 %i RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB\n" % (len(atmosphere_layers))
        #command_input += " 6.12960183E-04   3686.1 1.679E+01 2.580E+09 2.175E-04 4.386E-02 1.000E+05\n"
        #atm_kurucz.write("%.8e   %.1f %.3e %.3e %.3e %.3e %.3e" % (rhox[i], temperature[i], pgas[i], xne[i], abross[i], accrad[i], vturb[i]) )
        #command_input += "\n".join([" %.8E %8.1f %.3E %.3E %.3E %.3E %.3E" % (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], layer[6]) for layer in atmosphere_layers])
        #command_input += "\n".join([" %.8E %8.1f %.3E %.3E %.3E %.3E %.3E" % (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], 1.0e5) for layer in atmosphere_layers])
        # Force microturbulence in model to zero because later it will be used to add to the real microturbulence that we want:
        command_input += "\n".join([" %.8E %8.1f %.3E %.3E %.3E %.3E %.3E" % (layer[0], layer[1], layer[2], layer[3], layer[4], layer[5], 0.0e5) for layer in atmosphere_layers])
        command_input += "\nPRADK 1.4878E+00\n"
        command_input += "READ MOLECULES\n"
        command_input += "MOLECULES ON\n"
        command_input += "BEGIN                    ITERATION  15 COMPLETED\n"
        command_input += "END\n"
        command_input += "STOP\n"
        atm_file.write(command_input)
    atm_file.close()
    return atm_file.name

def model_atmosphere_is_closest_copy(modeled_layers_pack, target):
    """
        Returns True if model could not be interpolated
    """
    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, nlte_dep_grid, base_dirname = modeled_layers_pack
    target_point = []
    for param in free_parameters:
        target_point.append(target[param])
    target_point_cannot_be_interpolated = None
    for delaunay_triangulation in delaunay_triangulations['precomputed']:
        if delaunay_triangulation is None:
            continue
        simplex = delaunay_triangulation.find_simplex(target_point)
        if target_point_cannot_be_interpolated is None:
            target_point_cannot_be_interpolated = np.any(simplex == -1)
        else:
            target_point_cannot_be_interpolated = target_point_cannot_be_interpolated and np.any(simplex == -1)
    return target_point_cannot_be_interpolated if target_point_cannot_be_interpolated is not None else True



def load_modeled_layers_pack(input_path):
    """
    Restore modeled atmospheric layers and statistics, previously processed by
    iSpec (i.e. save_modeled_layers_pack). By default, iSpec is distributed with
    several pre-processed atmospheric models that can be loaded with this function.

    :param input_path:
        Name of the input file (i.e. models.dump) or directory (new interpolator)
    :type input_path: string

    :returns:
        List of modeled_layers, used_values_for_layers, proximity, teff_range, logg_range, MH_range, alpha_range and nlayers
    """
    if not os.path.isdir(input_path):
        raise Exception("Input path '{}' is not a directory".format(input_path))

    base_dirname = input_path
    atm_dirname = os.path.join(base_dirname, "grid")
    params_filename = os.path.join(base_dirname, "parameters.tsv")
    cache_filename = os.path.join(base_dirname, "cache.dump")

    if not os.path.exists(atm_dirname):
        raise Exception("Grid path '{}' does not exist".format(atm_dirname))
    if not os.path.exists(params_filename):
        raise Exception("Parameters file '{}' does not exist".format(params_filename))

    parameters = pd.read_csv(params_filename, sep="\t")
    # Subsets are used in spectral grids to spread the number of data points among
    # several Delaunay Triangulations, but this has not been necessary for model atmospheres
    parameters_subsets = [np.array([True]*len(parameters))]
    filenames = base_dirname + "/" + parameters['filename']
    filenames = np.asarray(filenames)
    ranges = {}
    free_parameters = parameters.columns.to_numpy().tolist()
    free_parameters = [x for x in free_parameters if not x.startswith("fixed_") and x != "filename"]
    for free_param in free_parameters:
        free_param_range = np.unique(parameters[free_param])
        ranges[free_param] = free_param_range
    existing_points = parameters[free_parameters]
    existing_points = np.array(existing_points)

    # The delaunay triangulation and kdtree can be computationally expensive,
    # do it once and save a pickle dump
    # But first verify that the computed one (if exists) matches what is expected
    use_dump = False
    if os.path.exists(cache_filename):
        use_dump = True
        delaunay_triangulations, kdtree = pickle.load(open(cache_filename, 'rb'))
        for delaunay_triangulation, parameters_subset in zip(delaunay_triangulations['precomputed'], parameters_subsets):
            if delaunay_triangulation is None:
                continue
            if not (delaunay_triangulation.ndim == len(free_parameters) and delaunay_triangulation.npoints == len(existing_points[parameters_subset])):
                use_dump = False
                break

    if not os.path.exists(cache_filename) or not use_dump:
        delaunay_triangulations = {'subsets': parameters_subsets, 'precomputed': []}
        for i, parameters_subset in enumerate(parameters_subsets):
            logging.info("Pre-computing model atmosphere triangulation [{}/{}]...".format(i+1, len(parameters_subsets)))
            if len(existing_points[parameters_subset]) > 0:
                delaunay_triangulations['precomputed'].append(spatial.Delaunay(existing_points[parameters_subset]))
            else:
                delaunay_triangulations['precomputed'].append(None)
        kdtree = spatial.cKDTree(existing_points)
        pickle.dump((delaunay_triangulations, kdtree), open(cache_filename, 'wb'), protocol=2)

    # Functions will receive the parameters in the same order
    read_point_value = lambda f: fits.open(f)[1].data

    value_fields = ["rhox", "temperature", "pgas", "xne", "abross", "accrad", "vturb", "logtau5", "depth", "pelectron"]
    value_fields += ["alpha_enhancement", "c_enhancement", "n_enhancement", "o_enhancement", "rapid_neutron_capture_enhancement", "slow_neutron_capture_enhancement"]
    value_fields += ["radius", "mass", "vmic"]

    # NLTE departure coefficient grids
    departure_grid_dirname = os.path.join(input_path, "dep-grid")
    nlte_dep_grid = {}
    if os.path.exists(departure_grid_dirname):
        available_elements_for_nlte_departures = get_available_elements_for_nlte_departures(departure_grid_dirname)
        #available_elements_for_nlte_departures = ('Fe',) # TODO: remove
        for i, available_element in enumerate(available_elements_for_nlte_departures):
            dep_grid_h5_filename = os.path.join(departure_grid_dirname, f"{available_element}_nlte_grid_data.h5")

            dep_cache_filename = os.path.join(departure_grid_dirname, os.path.splitext(os.path.basename(dep_grid_h5_filename))[0] + "_cache.dump")
            if not os.path.exists(dep_cache_filename) or not use_dump:
                with h5py.File(dep_grid_h5_filename, 'r') as f:
                    dep_existing_points = f['points'][:]
                dep_teff_logg_MH_existing_points = np.unique(dep_existing_points[:, :3], axis=0) # Ignore absolute abundances (A(X)) in the 4th position
                assert existing_points.shape[0] == dep_teff_logg_MH_existing_points.shape[0], "Departure coefficient should have the exact parameter coverage of the model atmosphere"

                logging.info(f"Pre-computing departure coefficient triangulation for '{available_element}' [{i+1}/{len(available_elements_for_nlte_departures)}]...")
                dep_delaunay_triangulations = spatial.Delaunay(dep_teff_logg_MH_existing_points)
                dep_kdtree = spatial.cKDTree(dep_teff_logg_MH_existing_points)
                pickle.dump((dep_delaunay_triangulations, dep_kdtree), open(dep_cache_filename, 'wb'), protocol=2)
            else:
                dep_delaunay_triangulations, dep_kdtree = pickle.load(open(dep_cache_filename, 'rb'))

            nlte_dep_grid[available_element] = (dep_delaunay_triangulations, dep_kdtree, os.path.abspath(dep_grid_h5_filename))

    return existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, nlte_dep_grid, base_dirname



def calculate_opacities(atmosphere_layers_file, abundances, MH, microturbulence_vel, wave_base, wave_top, wave_step, verbose=0, opacities_filename=None, tmp_dir=None, is_marcs_model=True):
    """
    abundances should have been already modified acording to MH
    """
    if not is_turbospectrum_support_enabled():
        raise Exception("Turbospectrum support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    turbospectrum_dir = ispec_dir + "/synthesizer/turbospectrum/"
    turbospectrum_data = turbospectrum_dir + "/DATA/"
    turbospectrum_babsma_lu = turbospectrum_dir + "bin/babsma_lu"

    if opacities_filename is None:
        # Temporary file
        out = tempfile.NamedTemporaryFile(mode="wt", delete=False, dir=tmp_dir, encoding='utf-8')
        out.close()
        opacities_filename = out.name

    command = turbospectrum_babsma_lu
    command_input = "'LAMBDA_MIN:'  '"+str(wave_base*10.)+"'\n"
    command_input += "'LAMBDA_MAX:'  '"+str(wave_top*10.)+"'\n"
    command_input += "'LAMBDA_STEP:' '"+str(wave_step*10.)+"'\n"
    command_input += "'MODELINPUT:' '"+atmosphere_layers_file+"'\n"
    if is_marcs_model:
        command_input += "'MARCS-FILE    :' '.true.'\n"
    else:
        command_input += "'MARCS-FILE    :' '.false.'\n"
    command_input += "'MODELOPAC:' '"+opacities_filename+"'\n"
    #command_input += "'METALLICITY:'    '"+str(MH)+"'\n"
    command_input += "'METALLICITY:'    '0.00'\n" # We have done the abundance changes already
    command_input += "'ALPHA/Fe   :'    '0.00'\n"
    command_input += "'HELIUM     :'    '0.00'\n"
    command_input += "'R-PROCESS  :'    '0.00'\n"
    command_input += "'S-PROCESS  :'    '0.00'\n"
    #command_input += "'INDIVIDUAL ABUNDANCES:'   '0'\n"
    #command_input += "'INDIVIDUAL ABUNDANCES:'   '1'\n"
    #command_input += "3  1.05\n"
    atom_abundances = abundances[abundances['code'] <= 92]
    if len(atom_abundances) != 92:
        raise Exception("No abundances for all 92 elements!")
    command_input += "'INDIVIDUAL ABUNDANCES:'   '"+str(len(atom_abundances))+"'\n"
    for atom_abundance in atom_abundances:
        abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
        command_input +=  "%i  %.2f\n" % (atom_abundance['code'], abund)
    command_input += "'XIFIX:' 'T'\n"
    command_input += str(microturbulence_vel)+"\n"

    tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
    os.symlink(turbospectrum_data, tmp_execution_dir+"/DATA")
    previous_cwd = os.getcwd()
    os.chdir(tmp_execution_dir)

    if verbose == 1:
        proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
    else:
        proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    # wait for the process to terminate
    out, err = proc.communicate(input=command_input.encode('utf-8'))
    errcode = proc.returncode

    os.chdir(previous_cwd)

    shutil.rmtree(tmp_execution_dir)

    return opacities_filename

#--------------------------------------------------------------------------------
# NLTE related:
#--------------------------------------------------------------------------------

def interpolate_nlte_departure_coefficients(modeled_layers_pack, abundances, target, fixed_abundances=None, linelist=None, regions=None, code="turbospectrum"):
    """
    Interpolates a quantity from a sparse HDF5 grid using Delaunay triangulation.

    If the target point is within the convex hull of the grid, the function performs
    linear barycentric interpolation within the enclosing tetrahedron.

    If the target point is outside the convex hull, it falls back to a nearest-neighbor
    approach using a KDTree to find the closest grid point in the (Teff, logg, [M/H])
    space and issues a warning.

    This function correctly handles irregular grids and maintains memory efficiency
    by only loading the necessary model data from the HDF5 file.

    Args:
        points_arr (np.ndarray): A 2D array of shape (N, 4) where each row is
                                 (Teff, logg, [M/H], A(X)).
        values_dset (h5py.Dataset): The HDF5 dataset containing the values to be
                                    interpolated. Must have shape (N, ...).
        target_teff (float): The target effective temperature.
        target_logg (float): The target surface gravity.
        target_MH (float): The target metallicity [M/H].
        target_absolute_abundance (float): The target absolute abundance A(X).

    Returns:
        np.ndarray: The interpolated (or nearest-neighbor) data array.
    """
    assert code == "turbospectrum", "NLTE calculations are only available with Turbospectrum"

    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, nlte_dep_grid, base_dirname = modeled_layers_pack
    nlte_departure_coefficients = {}

    if linelist is not None:
        # Ignore lines outside of the regions to be synthesize
        if regions is not None:
            linelist = _filter_linelist(linelist, regions)

        # Filter out lines not supported by a given synthesizer
        if code == 'moog-scat':
            # MOOG-SCAT is backward compatible with MOOG
            linelist_code = 'moog'
        else:
            linelist_code = code
        lcode = linelist[linelist_code+'_support'] == "T"
        linelist = linelist[lcode]

    target_teff = target['teff']
    target_logg = target['logg']
    target_MH = target['MH']
    target_alpha = target['alpha']

    #--------------------------------------------------------------------------------
    # Scale abundances to get the abundance used for departure coefficient interpolation
    # matching the same value that will be used for synthesis later too
    atom_abundances = enhance_solar_abundances(abundances, target_alpha)

    # Turbospectrum is not going to scale the abundances because we are indicating
    # our abundances in the input and that overrides any other prescription, thus
    # we have to manually scale (but do not change Hydrogen and Helium!)
    atom_abundances = atom_abundances[atom_abundances['code'] <= 92]
    if len(atom_abundances) != 92:
        raise Exception("No abundances for all 92 elements!")
    efilter = np.logical_and(atom_abundances['code'] != 1, atom_abundances['code'] != 2)
    atom_abundances['Abund'][efilter] += target_MH

    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        atom_abundances = atom_abundances.copy()
        for fixed_abundance in fixed_abundances:
            index = np.where(atom_abundances['code'] == fixed_abundance['code'])[0]
            atom_abundances['Abund'][index] = fixed_abundance['Abund']
    #--------------------------------------------------------------------------------


    for element_name, (dep_delaunay_triangulations, dep_kdtree, dep_grid_h5_filename) in nlte_dep_grid.items():
        with h5py.File(dep_grid_h5_filename, 'r') as f:
            dep_existing_points = f['points'][:]
            # Loading sparse grid parameters from HDF5 file...
            depart_values_dset = f['depart_values']
            tau_values_dset = f['tau_values']

            assert element_name == f.attrs['element_name']
            atom_model_filename = f.attrs['atom_model_filename']
            grid_filename = f.attrs['grid_filename']
            index_filename = f.attrs['index_filename']
            ndep = f.attrs['ndep']
            nk = f.attrs['nk']
            absolute_abundance = f.attrs['absolute_abundance'] # solar abundance used in the model atmosphere
            atomic_number = f.attrs['atomic_number']

            if atomic_number > 1: # hydrogen (atomic number == 1) is usually encoded into radiative transfer codes and not necessarily present in the atomic linelist
                # Check if there is any atomic line related to this element that has NLTE data
                element_lines_filter = np.floor(np.asarray(linelist['turbospectrum_species'], dtype=float)) == atomic_number
                element_linelist = linelist[element_lines_filter]
                linelist_has_nlte_data_for_element = 'T' in element_linelist['nlte'] and (np.any(element_linelist[element_linelist['nlte'] == 'T']['nlte_label_low'] != 'none') or np.any(element_linelist[element_linelist['nlte'] == 'T']['nlte_label_up'] != 'none'))
                if not linelist_has_nlte_data_for_element:
                    # Skip interpolattion of departure coefficient for elements that do not have any line with NLTE data in the linelist
                    continue

            target_absolute_abundance = atom_abundances[atom_abundances['code'] == atomic_number]['Abund'][0] + 12.036 # + 0.004 # turbospectrum check (bsyn.f) seems to rely on + 12.036 + 0.004

            interpolated_departures = _interpolate_from_h5_grid_delaunay(
                dep_existing_points,
                dep_delaunay_triangulations,
                dep_kdtree,
                depart_values_dset,
                target_teff,
                target_logg,
                target_MH,
                target_absolute_abundance
            )

            interpolated_tau = _interpolate_from_h5_grid_delaunay(
                dep_existing_points,
                dep_delaunay_triangulations,
                dep_kdtree,
                tau_values_dset,
                target_teff,
                target_logg,
                target_MH,
                target_absolute_abundance
            )

            if np.any(np.isnan(interpolated_departures)) or np.any(np.isinf(interpolated_departures)):
                # Some original departure coefficient are nan and/or inf (e.g., teff = 3004 logg = 2.44 MH = -0.1)
                # thus the interpolated coefficients also contain nan and/or inf
                # but then turbospectrum may fail (segmentation fault) if it encounter lines that require those values
                # hence we fill them with either interpolating or extrapolating on the tau dimension (depth)
                interpolated_departures = _fill_missing_departures(interpolated_tau, interpolated_departures)
                max_value = np.nanmax(interpolated_departures[~np.isinf(interpolated_departures)])
                interpolated_departures[np.isnan(interpolated_departures)] = max_value
                interpolated_departures[np.isinf(interpolated_departures)] = max_value
            #

        absolute_atom_model_filename = os.path.join(os.path.dirname(dep_grid_h5_filename), "atoms", os.path.basename(atom_model_filename))
        #nlte_departure_coefficients[element_name] = (interpolated_departures, interpolated_tau, ndep, nk, atomic_number, absolute_abundance, absolute_atom_model_filename, )
        nlte_departure_coefficients[element_name] = (interpolated_departures, interpolated_tau, ndep, nk, atomic_number, target_absolute_abundance, absolute_atom_model_filename, )
        logging.debug(f"Interpolated NLTE departures coefficients for '{element_name}'")
    if len(nlte_departure_coefficients) == 0:
        logging.warning("No NLTE departures coefficients in this model atmosphere or no atomic lines with NLTE data in the considered regions")
    return nlte_departure_coefficients


def _fill_missing_departures(tau, departures):
    filled_departures = departures.copy()
    nk, ndep = filled_departures.shape
    filled_departures[np.isinf(departures)] = np.nan

    # --- 2. Iterate and interpolate with np.interp ---
    for i in range(nk):
        departure_row = filled_departures[i, :]
        is_finite = np.isfinite(departure_row)

        if not np.any(is_finite):
            print(f"Warning: Row {i} contains no valid data. Filling with 1.0.")
            filled_departures[i, :] = 1.0
            continue

        # np.interp needs the x-coordinates of the points to be interpolated,
        # and the (x, y) coordinates of the known data points.
        # It automatically handles extrapolation by repeating the first/last y-value.
        filled_departures[i, :] = np.interp(
            tau,                        # x-values to evaluate at (the full grid)
            tau[is_finite],             # x-coordinates of known points
            departure_row[is_finite]    # y-coordinates of known points
        )
    return filled_departures


def _interpolate_from_h5_grid_delaunay(points_arr, delaunay_grid, kdtree, values_dset, target_teff, target_logg, target_MH, target_absolute_abundance):

    # Isolate the 3D grid points for triangulation/search
    grid_points_3d = np.unique(points_arr[:, :3], axis=0)
    target_point_3d = np.array([target_teff, target_logg, target_MH])

    # Find the enclosing tetrahedron
    simplex_index = delaunay_grid.find_simplex(target_point_3d)

    # Check if point is inside or outside the convex hull
    if simplex_index == -1:
        # The point is outside the convex hull. Use KDTree to find the nearest neighbor.
        # - Query the tree for the single closest point (k=1)
        distance, closest_vertex_idx_in_3d_grid = kdtree.query(target_point_3d, k=1)
        closest_vertex_coord = grid_points_3d[closest_vertex_idx_in_3d_grid]

        # - Find the best-matching model from the full `points_arr` based on the target abundance `A(X)`.
        mask = (points_arr[:, 0] == closest_vertex_coord[0]) & \
               (points_arr[:, 1] == closest_vertex_coord[1]) & \
               (points_arr[:, 2] == closest_vertex_coord[2])
        candidate_indices = np.where(mask)[0]

        # Of these candidates, find the one with the A(X) closest to the target A(X)
        candidate_abundances = points_arr[candidate_indices, 3]
        closest_ax_idx_in_candidates = np.argmin(np.abs(candidate_abundances - target_absolute_abundance))
        final_grid_index_to_load = candidate_indices[closest_ax_idx_in_candidates]

        logging.warning(
            f"Target point {target_point_3d} is outside the convex hull of the grid. "
            f"Falling back to the nearest grid point: {final_grid_index_to_load} with A(X) = {points_arr[final_grid_index_to_load, 3]:.2f}"
        )

        result = values_dset[final_grid_index_to_load]
    else:
        # Point is inside the hull

        # - Identify the 4 vertices of the tetrahedron and resolve the A(X) dimension
        # These are the indices of the vertices in `grid_points_3d`
        vertex_indices_in_3d_grid = delaunay_grid.simplices[simplex_index]
        # These are the coordinates of the 4 vertices
        vertices_coords = grid_points_3d[vertex_indices_in_3d_grid]

        # - For each vertex, find the best-matching model from the full `points_arr` based on the target abundance `A(X)`.
        indices_to_load = np.zeros(4, dtype=int)
        for i, vertex_coord in enumerate(vertices_coords):
            # Find all models in the original grid that match this vertex's (T,g,M)
            mask = (points_arr[:, 0] == vertex_coord[0]) & \
                   (points_arr[:, 1] == vertex_coord[1]) & \
                   (points_arr[:, 2] == vertex_coord[2])
            candidate_indices = np.where(mask)[0]
            # Of these candidates, find the one with the A(X) closest to the target A(X)
            candidate_abundances = points_arr[candidate_indices, 3]
            closest_ax_idx_in_candidates = np.argmin(np.abs(candidate_abundances - target_absolute_abundance))
            final_grid_index = candidate_indices[closest_ax_idx_in_candidates]

            indices_to_load[i] = final_grid_index
            #print(f"  Vertex {i} ({vertex_coord[0]:.0f}, {vertex_coord[1]:.2f}, {vertex_coord[2]:.2f}): "
            #      f"Using model index {final_grid_index} with A(X) = {points_arr[final_grid_index, 3]:.2f}")

        # Read ONLY the required data from HDF5
        sort_order = np.argsort(indices_to_load)
        sorted_indices_to_load = indices_to_load[sort_order]
        vertex_values_sorted_flat = values_dset[sorted_indices_to_load] # This is the slowest part, limited by file disk read speed

        value_shape = values_dset.shape[1:]
        vertex_values_flat = np.empty((4,) + value_shape, dtype=vertex_values_sorted_flat.dtype)
        vertex_values_flat[sort_order] = vertex_values_sorted_flat

        # Perform interpolation using barycentric coordinates
        # - This is the standard way to do linear interpolation inside a simplex.
        # - The transform matrix helps convert Cartesian coordinates to barycentric coordinates.
        # - It is structured as (ndim+1, ndim), where the first ndim rows are the inverse transform matrix and the last row is the offset vector.
        transform_matrix = delaunay_grid.transform[simplex_index]
        # Separate the inverse transform (3x3) and the offset vector (3,)
        T_inv = transform_matrix[:3, :3]
        offset_r = transform_matrix[3, :]
        # Calculate the first 3 barycentric coordinates
        b = T_inv @ (target_point_3d - offset_r)
        # The 4 barycentric coordinates must sum to 1. The fourth is 1 - sum(others).
        bary_coords = np.append(b, 1 - b.sum())

        # The interpolated result is the weighted average of the vertex values, where the weights are the barycentric coordinates.
        # We need to reshape the barycentric coords to allow broadcasting with multi-dimensional values.
        bary_coords_reshaped = bary_coords.reshape((4,) + (1,) * len(value_shape))

        # Performing interpolation using barycentric coordinates...
        with np.errstate(invalid='ignore'):
            # NOTE: Omit warning "RuntimeWarning: invalid value encountered in multiply"
            # which happens because there can be NaN and/or inifites involved
            # but they will be propagated and fixed later on before using themm with turbospectrum
            result = np.sum(vertex_values_flat * bary_coords_reshaped, axis=0)

    return result


def get_available_elements_for_nlte_departures(grid_directory):
    """
    Scans a directory to find available elements based on file naming conventions.

    An element is considered "available" if both a data file (`{El}_..._data.h5`)
    and an index file (`{El}_..._index.tsv`) exist for it.

    Args:
        grid_directory (str): The path to the directory containing grid files.

    Returns:
        list[str]: A sorted list of available element symbols (e.g., ['Al', 'Fe', 'Na']).
                   Returns an empty list if the directory is not found or is empty.
    """
    if not os.path.isdir(grid_directory):
        print(f"Warning: Grid directory '{grid_directory}' not found.")
        return []

    try:
        all_files = os.listdir(grid_directory)
    except OSError as e:
        print(f"Warning: Could not read directory '{grid_directory}'. Error: {e}")
        return []

    # Use a set for efficient storage and automatic deduplication
    found_elements = set()

    for filename in all_files:
        # Check for the primary data file
        if filename.endswith("_nlte_grid_data.h5"):
            # The element symbol is the part before the first underscore
            element_symbol = filename.split('_')[0]

            # Construct the name of the required companion index file
            index_filename = f"{element_symbol}_nlte_grid_index.tsv"

            # Check if the companion file also exists in the directory
            if index_filename in all_files:
                found_elements.add(element_symbol)

    # Return a sorted list for consistent order
    return sorted(list(found_elements))

