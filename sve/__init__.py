"""
    This file is part of Spectra Visual Editor (SVE).
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com

    SVE is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SVE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with SVE. If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division

import os
import sys

if sys.hexversion < 0x02050000:
    raise RuntimeError("SVE requires at least Python 2.5")

## PyInstaller resource access
def resource_path(relative):
    if getattr(sys, 'frozen', None):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(__file__)
    return os.path.join(basedir, relative)

#if os.path.exists(resource_path("synthesizer.so")):
try:
    from abundances import determine_abundances
    from abundances import read_SPECTRUM_abundances
    from atmospheres import interpolate_atmosphere_layers
    from atmospheres import load_modeled_layers_pack
    from atmospheres import valid_atmosphere_target
    from synth import generate_spectrum
    from synth import modelize_spectrum
    from synth import read_SPECTRUM_linelist
except ImportError as e:
    pass
from common import calculate_barycentric_velocity_correction, sigma_clipping, interquartile_range_filtering
from spectrum import *
from continuum import *
from lines import *
from segments import *
from mpfitmodels import GaussianModel, VoigtModel
import log
import logging
