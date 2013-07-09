#!/usr/bin/env python
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

class Spectrum():
    def __init__(self, data, name, color='b', not_saved=False, plot_id=None, continuum_model=None, continuum_data=None, continuum_plot_id=None, path=None):
        self.data = data
        self.name = name
        self.color = color
        self.not_saved = not_saved
        self.plot_id = plot_id
        self.continuum_model = continuum_model
        self.continuum_data = continuum_data
        self.continuum_plot_id = continuum_plot_id
        self.errors_plot_id1 = None
        self.errors_plot_id2 = None
        self.velocity_atomic = 0.0 # Radial velocity (km/s)
        self.velocity_telluric = 0.0 # Barycentric velocity (km/s)
        self.velocity_template = 0.0 # Barycentric velocity (km/s)
        self.path = path
        # Resolution power from the velocity profile relative to atomic and telluric lines
        self.resolution_atomic = 0.0
        self.resolution_telluric = 0.0
        self.velocity_profile_atomic_xcoord = None
        self.velocity_profile_atomic_fluxes = None
        self.velocity_profile_atomic_errors = None
        self.velocity_profile_atomic_models = None
        self.velocity_profile_atomic_num_used_lines = None
        self.velocity_profile_atomic_rv_step = None
        self.velocity_profile_telluric_xcoord = None
        self.velocity_profile_telluric_fluxes = None
        self.velocity_profile_telluric_errors = None
        self.velocity_profile_telluric_models = None
        self.velocity_profile_telluric_num_used_lines = None
        self.velocity_profile_telluric_rv_step = None
        self.velocity_profile_telluric_fwhm_correction = 0.0
        self.velocity_profile_template_xcoord = None
        self.velocity_profile_template_fluxes = None
        self.velocity_profile_template_errors = None
        self.velocity_profile_template_models = None
        self.velocity_profile_template_num_used_lines = None
        self.velocity_profile_template_rv_step = None
        self.velocity_profile_template = None
        self.velocity_profile_internal_template = None
        self.snr = None
        self.linemasks = None # Linemasks that has been fitted and cross-matched with atomic & telluric data
        self.abundances = None
        self.abundances_teff = 5777.0
        self.abundances_logg = 4.44
        self.abundances_MH = 0.02
        self.abundances_microturbulence_vel = 2.0
        self.dialog = {}


