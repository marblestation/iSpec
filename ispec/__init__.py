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
from common import is_turbospectrum_support_enabled
from common import is_spectrum_support_enabled
from common import is_moog_support_enabled
from common import is_width_support_enabled
from common import is_ares_support_enabled
from common import is_synthe_support_enabled
from common import is_sme_support_enabled
from common import estimate_vmic
from common import estimate_vmac
from common import estimate_mass_radius
from common import calculate_barycentric_velocity_correction
from common import sigma_clipping
from common import interquartile_range_filtering
from common import save_results
from common import restore_results
from common import mkdir_p
from spectrum import read_spectrum
from spectrum import write_spectrum
from spectrum import normalize_spectrum
from spectrum import estimate_snr
from spectrum import convolve_spectrum
from spectrum import create_spectrum_structure
from spectrum import resample_spectrum
from spectrum import correct_velocity
from spectrum import correct_velocity_regions
from spectrum import add_noise
from spectrum import random_realizations
from spectrum import create_wavelength_filter
from spectrum import air_to_vacuum
from spectrum import vacuum_to_air
from spectrum import create_filter_cosmic_rays
from continuum import read_continuum_regions
from continuum import write_continuum_regions
from continuum import fit_continuum
from continuum import find_continuum
from continuum import generate_random_continuum
from lines import read_telluric_linelist
from lines import read_isotope_data
from lines import write_isotope_data
from lines import read_cross_correlation_mask
from lines import read_chemical_elements
from lines import read_molecular_symbols
from lines import read_line_regions
from lines import write_line_regions
from lines import create_linemasks_structure
from lines import read_atomic_linelist
from lines import write_atomic_linelist
from lines import find_linemasks
from lines import reset_fitted_data_fields
from lines import fit_lines
from lines import adjust_linemasks
from lines import create_filter_for_regions_affected_by_tellurics
from lines import cross_correlate_with_mask
from lines import cross_correlate_with_template
from lines import select_good_velocity_profile_models
from lines import update_ew_with_ares
from lines import van_der_Waals_ABO_to_single_gamma_format
from segments import read_segment_regions
from segments import write_segment_regions
from segments import create_segments_around_lines
from modeling.mpfitmodels import GaussianModel
from modeling.mpfitmodels import VoigtModel
from modeling.ew import model_spectrum_from_ew
from modeling.ssf import model_spectrum
from synth.common import generate_fundamental_spectrum
from synth.common import generate_spectrum
from synth.effects import apply_post_fundamental_effects
from synth.spectrum import calculate_theoretical_ew_and_depth
from synth.grid import load_spectral_grid
from synth.grid import valid_interpolated_spectrum_target
from synth.grid import precompute_synthetic_grid
from synth.grid import estimate_initial_ap
import log
import logging
