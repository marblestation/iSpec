"""
    This file is part of the Integrated Spectroscopic Framework (iSpec).
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com

    iSpec is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    iSpec is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division

import os
import sys

if sys.hexversion < 0x02050000:
    raise RuntimeError("iSpec requires at least Python 2.5")

## PyInstaller resource access
def resource_path(relative):
    if getattr(sys, 'frozen', None):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(__file__)
    return os.path.join(basedir, relative)


from abundances import determine_abundances
from abundances import read_solar_abundances
from abundances import enhance_solar_abundances
from abundances import determine_abundance_enchancements
from abundances import create_free_abundances_structure
from atmospheres import interpolate_atmosphere_layers
from atmospheres import load_modeled_layers_pack
from atmospheres import valid_atmosphere_target
from atmospheres import write_atmosphere
from common import is_turbospectrum_support_enabled, is_spectrum_support_enabled, is_moog_support_enabled, is_width_support_enabled, is_ares_support_enabled, is_synthe_support_enabled, is_sme_support_enabled
from synth import generate_spectrum
from synth import load_spectral_grid, valid_interpolated_spectrum_target
from synth import model_spectrum, precompute_synthetic_grid, estimate_initial_ap
from common import estimate_vmic, estimate_vmac, estimate_mass_radius
from synth import model_spectrum_from_ew
from synth import calculate_theoretical_ew_and_depth

from common import calculate_barycentric_velocity_correction, sigma_clipping, interquartile_range_filtering, save_results, restore_results, mkdir_p
from spectrum import *
from continuum import *
from lines import *
from segments import *
from mpfitmodels import GaussianModel, VoigtModel
import log
import logging
