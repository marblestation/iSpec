#
#    This file is part of the Integrated Spectroscopic Framework (iSpec).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
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
import cPickle as pickle
import log
import logging
import subprocess
import shutil
import pandas as pd
from astropy.io import fits
from scipy import spatial
from scipy.interpolate import LinearNDInterpolator
import glob
from common import is_turbospectrum_support_enabled, is_spectrum_support_enabled

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
    existing_points, existing_point_filename_pattern_builder, read_point_value, value_fields, delaunay_triangulation, kdtree, teff_range, logg_range, MH_range, base_dirname = modeled_layers_pack

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

def __find_matching_filename(filename_pattern):
    filenames = glob.glob(filename_pattern)
    if len(filenames) > 1:
        raise Exception("Ambigous filename pattern '{}'".format(filename_pattern))
    elif len(filenames) == 0:
        raise Exception("Non existent filenames with pattern '{}'".format(filename_pattern))
    return filenames[0]

def __closest(kdtree, existing_points, existing_point_filename_pattern_builder, read_point_value, target_point):
    """
    If there is no model in the extreme points, copy the closest one.
    input:
    - existing_points: points in the parameters space such as [( 3000.,  3.5,  0.  ,  0. ), ( 3000.,  3.5,  0.25,  0. ), ...]
    - existing_point_filename_pattern_builder: function that returns a filename pattern to find the matching model
    - target_point: parameters for the target model
    output:
    - created extreme points such as [(3000, 0.0, -5.0, 0.0), (3000, 0.0, -5.0, 0.40000000000000002),]
    """
    distance, index = kdtree.query(target_point, k=1)
    closest_existing_point = existing_points[index]
    closest_existing_point_filename = __find_matching_filename(existing_point_filename_pattern_builder(closest_existing_point))
    value = read_point_value(closest_existing_point_filename)
    logging.info("Closest to target point '{}' is '{}'".format(" ".join(map(str, target_point)), " ".join(map(str, closest_existing_point))))
    return value

def _interpolate(delaunay_triangulation, kdtree, existing_points, existing_point_filename_pattern_builder, read_point_value, value_fields, target_point):
    """
    input:
    - existing_points: points in the parameters space such as [( 3000.,  3.5,  0.  ,  0. ), ( 3000.,  3.5,  0.25,  0. ), ...]
    - existing_point_filename_pattern_builder: function that returns a filename pattern to find the matching model
    - read_point_value: function to read the value from the model filename
    - target_point: parameters for the target model
    output:
    - interpolated value
    """
    simplex = delaunay_triangulation.find_simplex(target_point)
    if np.any(simplex == -1):
        logging.warn("Target point '{}' is out of bound, using the closest".format(" ".join(map(str, target_point))))
        return __closest(kdtree, existing_points, existing_point_filename_pattern_builder, read_point_value, target_point)
    index = delaunay_triangulation.simplices[simplex]
    points = []
    values = {}
    existing_points = np.array(existing_points)
    for point in existing_points[index]:
        filename = __find_matching_filename(existing_point_filename_pattern_builder(point))
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
        interpolated_value = interpolator((target_point))
        if interpolated_values is None:
            interpolated_values = pd.DataFrame(interpolated_value, columns=[field])
        else:
            interpolated_values[field] = interpolated_value
    return interpolated_values.to_records(index=False)



def interpolate_atmosphere_layers(modeled_layers_pack,  teff_target, logg_target, MH_target, code="spectrum"):
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

    if code == "turbospectrum" and "MARCS" not in base_dirname:
        # Spherical models in turbospectrum require a parameters that is only provided in MARCS model atmosphere
        raise Exception("Turbospectrum can only be used with MARCS model atmospheres.")

    if not valid_atmosphere_target(modeled_layers_pack, teff_target, logg_target, MH_target):
        raise Exception("Target parameters '{} {} {}' are out of range.".format(teff_target, logg_target, MH_target))

    existing_points, existing_point_filename_pattern_builder, read_point_value, value_fields, delaunay_triangulation, kdtree, teff_range, logg_range, MH_range, base_dirname = modeled_layers_pack
    target_point = (teff_target, logg_target, MH_target)
    interpolated_atm = _interpolate(delaunay_triangulation, kdtree, existing_points, existing_point_filename_pattern_builder, read_point_value, value_fields, target_point)

    interpolated_atm_compatible_format = interpolated_atm[["rhox", "temperature", "pgas", "xne", "abross", "accrad", "vturb", "logtau5", "depth", "pelectron"]]
    interpolated_atm_compatible_format = interpolated_atm_compatible_format.view(float).reshape(interpolated_atm_compatible_format.shape + (-1,))
    return interpolated_atm_compatible_format


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
        atm_file = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)

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
        lgTauR = -5.00
        for i, layer in enumerate(atmosphere_layers):
            atm_file.write("%i %.2f %.4f %.3E %.1f %.4E %.4E %.4E %.4E\n" % (i+1, lgTauR, layer[7], layer[8], layer[1], layer[9], layer[2], layer[5], 0.))
            lgTauR += 0.20
    elif code == "spectrum":
        # Spectrum
        # mass depth, temperature in kelvin, gas pressure, electron density, Rosseland mean absorption coefficient, radiation pressure, microturbulent velocity in meters/second.
        atm_file.write("%.1f  %.5f  %.2f  %i\n" % (teff, logg, MH, len(atmosphere_layers)) )
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

def model_atmosphere_is_closest_copy(modeled_layers_pack, teff_target, logg_target, MH_target):
    """
        Returns True if model could not be interpolated
    """
    existing_points, existing_point_filename_pattern_builder, read_point_value, value_fields, delaunay_triangulation, kdtree, teff_range, logg_range, MH_range, base_dirname = modeled_layers_pack
    target_point = (teff_target, logg_target, MH_target)
    simplex = delaunay_triangulation.find_simplex(target_point)
    target_point_cannot_be_interpolated = np.any(simplex == -1)
    return target_point_cannot_be_interpolated



def load_modeled_layers_pack(input_path):
    """
    Restore modeled atmospheric layers and statistics, previously processed by
    iSpec (i.e. save_modeled_layers_pack). By default, iSpec is distributed with
    several pre-processed atmospheric models that can be loaded with this function.

    :param input_path:
        Name of the input file (i.e. models.dump) or directory (new interpolator)
    :type input_path: string

    :returns:
        List of modeled_layers, used_values_for_layers, proximity, teff_range, logg_range, MH_range and nlayers
    """
    if not os.path.isdir(input_path):
        raise Exception("Input path '{}' is not a directory".format(input_path))

    base_dirname = input_path
    atm_dirname = os.path.join(base_dirname, "grid")
    params_filename = os.path.join(base_dirname, "parameters.tsv")

    if not os.path.exists(atm_dirname):
        raise Exception("Grid path '{}' does not exist".format(atm_dirname))
    if not os.path.exists(params_filename):
        raise Exception("Parameters file '{}' does not exist".format(params_filename))

    params = pd.read_csv(params_filename, sep="\t", names=["teff", "logg", "mh"])
    teff_range = np.unique(params['teff'])
    logg_range = np.unique(params['logg'])
    MH_range = np.unique(params['mh'])
    existing_points = zip(params["teff"], params["logg"], params["mh"])

    delaunay_triangulation = spatial.Delaunay(existing_points)
    kdtree = spatial.KDTree(existing_points)

    # Functions will receive the parameters in the same order
    existing_point_filename_pattern_builder = lambda p: os.path.join(atm_dirname, "[ps]{0:0.0f}_g{1:+0.1f}*_m?.?_t??_st_z{2:+0.2f}_a*.fits.gz".format(*p))
    read_point_value = lambda f: fits.open(f)[1].data

    value_fields = ["rhox", "temperature", "pgas", "xne", "abross", "accrad", "vturb", "logtau5", "depth", "pelectron"]
    value_fields += ["alpha_enhancement", "c_enhancement", "n_enhancement", "o_enhancement", "rapid_neutron_capture_enhancement", "slow_neutron_capture_enhancement"]
    value_fields += ["radius", "mass", "vmic"]
    return existing_points, existing_point_filename_pattern_builder, read_point_value, value_fields, delaunay_triangulation, kdtree, teff_range, logg_range, MH_range, base_dirname


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
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
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
    out, err = proc.communicate(input=command_input)
    errcode = proc.returncode

    os.chdir(previous_cwd)

    shutil.rmtree(tmp_execution_dir)

    return opacities_filename
