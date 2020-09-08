from __future__ import absolute_import
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
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import zip
from builtins import map
from builtins import object
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
from .common import is_turbospectrum_support_enabled, is_spectrum_support_enabled

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
    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname = modeled_layers_pack

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
    if code not in ['spectrum', 'turbospectrum', 'moog', 'width', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname = modeled_layers_pack

    if code == "turbospectrum" and "MARCS" not in base_dirname:
        # Spherical models in turbospectrum require a parameters that is only provided in MARCS model atmosphere
        raise Exception("Turbospectrum can only be used with MARCS model atmospheres.")

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
        # Bytes should be swapped and FITS converted to numpy array:
        interpolated_atm_compatible_format = pd.DataFrame(np.array(interpolated_atm.byteswap().newbyteorder()))[compatible_fields].to_records(index=False)
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
    if code not in ['spectrum', 'turbospectrum', 'moog', 'width', 'synthe']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

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
        nvalues = len(atmosphere_layers[0])
        if nvalues != 11:
            raise Exception("Turbospectrum can only be used with MARCS model atmospheres.")
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
    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname = modeled_layers_pack
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
            logging.info("Pre-computing [{}/{}]...".format(i+1, len(parameters_subsets)))
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
    return existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname


def calculate_opacities(atmosphere_layers_file, abundances, MH, microturbulence_vel, wave_base, wave_top, wave_step, verbose=0, opacities_filename=None, tmp_dir=None):
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
    command_input += "'MARCS-FILE:' '.true.'\n"
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
