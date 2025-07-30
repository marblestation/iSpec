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
import time
from datetime import datetime, timedelta
import numpy as np
import logging

from .mpfitmodels import MPFitModel
from ispec.abundances import write_solar_abundances
from ispec.abundances import determine_abundance_enchancements
from ispec.atmospheres import write_atmosphere, interpolate_atmosphere_layers, model_atmosphere_is_closest_copy, interpolate_nlte_departure_coefficients
from ispec.lines import write_atomic_linelist, write_isotope_data, _get_atomic_linelist_definition
from ispec.common import estimate_vmic, estimate_vmac
from ispec.spectrum import create_spectrum_structure, convolve_spectrum, resample_spectrum, read_spectrum, create_wavelength_filter, read_spectrum, normalize_spectrum
from ispec.synth.effects import _filter_linelist, apply_post_fundamental_effects
from .common import Constants, _filter_linemasks_not_in_segments, _create_comparing_mask, _get_stats_per_linemask
from ispec.synth.spectrum import _create_waveobs_mask
from ispec.synth.common import generate_fundamental_spectrum
from ispec.synth.grid import load_spectral_grid, valid_interpolated_spectrum_target
from ispec.segments import merge_overlapping_regions

PARINFO_BASE = 10

class SynthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, isotopes, linelist_free_loggf, abundances, enhance_abundances=True, scale=None, teff=5000, logg=3.0, MH=0.0, alpha=0.0, vmic=2.0, vmac=0.0, vsini=2.0, limb_darkening_coeff=0.0, R=0, precomputed_grid_dir=None, grid=None, normalize_func=None):
        self.precomputed_grid_dir = precomputed_grid_dir
        self.grid = grid
        self.normalize_func = normalize_func
        self.elements = {}
        #self.elements["1"] = "H"
        #self.elements["2"] = "He"
        self.elements["3"] = "Li"
        self.elements["4"] = "Be"
        self.elements["5"] = "B"
        self.elements["6"] = "C"
        self.elements["7"] = "N"
        self.elements["8"] = "O"
        self.elements["9"] = "F"
        self.elements["10"] = "Ne"
        self.elements["11"] = "Na"
        self.elements["12"] = "Mg"
        self.elements["13"] = "Al"
        self.elements["14"] = "Si"
        self.elements["15"] = "P"
        self.elements["16"] = "S"
        self.elements["17"] = "Cl"
        self.elements["18"] = "Ar"
        self.elements["19"] = "K"
        self.elements["20"] = "Ca"
        self.elements["21"] = "Sc"
        self.elements["22"] = "Ti"
        self.elements["23"] = "V"
        self.elements["24"] = "Cr"
        self.elements["25"] = "Mn"
        self.elements["26"] = "Fe"
        self.elements["27"] = "Co"
        self.elements["28"] = "Ni"
        self.elements["29"] = "Cu"
        self.elements["30"] = "Zn"
        self.elements["31"] = "Ga"
        self.elements["32"] = "Ge"
        self.elements["33"] = "As"
        self.elements["34"] = "Se"
        self.elements["35"] = "Br"
        self.elements["36"] = "Kr"
        self.elements["37"] = "Rb"
        self.elements["38"] = "Sr"
        self.elements["39"] = "Y"
        self.elements["40"] = "Zr"
        self.elements["41"] = "Nb"
        self.elements["42"] = "Mo"
        self.elements["43"] = "Tc"
        self.elements["44"] = "Ru"
        self.elements["45"] = "Rh"
        self.elements["46"] = "Pd"
        self.elements["47"] = "Ag"
        self.elements["48"] = "Cd"
        self.elements["49"] = "In"
        self.elements["50"] = "Sn"
        self.elements["51"] = "Sb"
        self.elements["52"] = "Te"
        self.elements["53"] = "I"
        self.elements["54"] = "Xe"
        self.elements["55"] = "Cs"
        self.elements["56"] = "Ba"
        self.elements["57"] = "La"
        self.elements["58"] = "Ce"
        self.elements["59"] = "Pr"
        self.elements["60"] = "Nd"
        self.elements["61"] = "Pm"
        self.elements["62"] = "Sm"
        self.elements["63"] = "Eu"
        self.elements["64"] = "Gd"
        self.elements["65"] = "Tb"
        self.elements["66"] = "Dy"
        self.elements["67"] = "Ho"
        self.elements["68"] = "Er"
        self.elements["69"] = "Tm"
        self.elements["70"] = "Yb"
        self.elements["71"] = "Lu"
        self.elements["72"] = "Hf"
        self.elements["73"] = "Ta"
        self.elements["74"] = "W"
        self.elements["75"] = "Re"
        self.elements["76"] = "Os"
        self.elements["77"] = "Ir"
        self.elements["78"] = "Pt"
        self.elements["79"] = "Au"
        self.elements["80"] = "Hg"
        self.elements["81"] = "Tl"
        self.elements["82"] = "Pb"
        self.elements["83"] = "Bi"
        self.elements["84"] = "Po"
        self.elements["85"] = "At"
        self.elements["86"] = "Rn"
        self.elements["87"] = "Fr"
        self.elements["88"] = "Ra"
        self.elements["89"] = "Ac"
        self.elements["90"] = "Th"
        self.elements["91"] = "Pa"
        self.elements["92"] = "U"
        self.elements["101"] = "Md"
        self.elements["106"] = "Sg"
        self.elements["107"] = "Bh"
        self.elements["108"] = "Hs"
        self.elements["112"] = "Cn"
        self.elements["113"] = "Uut"
        self.elements["114"] = "Uuq"

        self.modeled_layers_pack = modeled_layers_pack
        self.linelist = linelist
        self.isotopes = isotopes
        self.linelist_free_loggf = linelist_free_loggf
        self.abundances = abundances
        self.enhance_abundances = enhance_abundances
        self.scale = scale
        #
        self.calculation_time = 0
        self.waveobs = None
        self.segments = None
        self.waveobs_mask = None
        self.cache = {}
        p = [teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, R ]
        #
        self.abundances_files = []
        self.linelist_file = None
        self.isotope_file = None
        self.molecules_files = None
        self.atmosphere_layers_files = []
        super(SynthModel, self).__init__(p)

    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in range(len(p)):
                self._parinfo[i]['value'] = p[i]

        self.last_final_fluxes = [None] * self.n_stellar_components
        for component in range(self.n_stellar_components):
            key = "%.0f %.2f %.2f %.2f %.2f " % (self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), self.vmic(component=component))
            complete_key = "%.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %i %.2f" % (self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), self.vmic(component=component), self.vmac(component=component), self.vsini(component=component), self.limb_darkening_coeff(component=component), int(self.R(component=component)), self.lf(component=component))

            # Consider new loggf
            linelist_free_loggf = self.generate_linelist_free_loggf(component=component)

            loggf_key = " ".join(["%.3f" % (x) for x in linelist_free_loggf['loggf']])
            complete_key += " loggf [" + loggf_key + "]"
            key += loggf_key

            if len(linelist_free_loggf) > 0:
                linelist = np.hstack((self.linelist, linelist_free_loggf))
                linelist.sort(order=['wave_nm'])
            else:
                linelist = self.linelist

            # Consider new abundances as fixed
            fixed_abundances = self.free_abundances(component=component)

            abundances_key = " ".join(["%.2f" % (x) for x in fixed_abundances['Abund']])
            complete_key += " abund [" + abundances_key + "]"
            key += abundances_key

            # vrad
            vrad = self.vrad(component=component)
            if np.all(vrad[0] == vrad):
                # When fitting a spectroscopic binary, the second component has a fixed vrad for all segments, simplify by showing just one value
                vrad = np.unique(vrad)
            vrad_key = " ".join(["%.2f" % (x) for x in vrad])
            complete_key += " vrad [" + vrad_key + "]"

            ##### [start] Check precomputed (solar abundance)
            filename = "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}_{7:.2f}.fits.gz".format(int(self.teff(component=component)), self.logg(component=component), self.MH(component=component), self.alpha(component=component), self.vmic(component=component), self.vmac(component=component), self.vsini(component=component), self.limb_darkening_coeff(component=component))
            fundamental_filename = "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}_{7:.2f}.fits.gz".format(int(self.teff(component=component)), self.logg(component=component), self.MH(component=component), self.alpha(component=component), self.vmic(component=component), 0., 0., 0.)
            precomputed_file = os.path.join(str(self.precomputed_grid_dir), "grid", filename)
            precomputed_step_file = os.path.join(str(self.precomputed_grid_dir), "steps", filename)
            fundamental_precomputed_file = os.path.join(str(self.precomputed_grid_dir), "grid", fundamental_filename)
            fundamental_precomputed_step_file = os.path.join(str(self.precomputed_grid_dir), "steps", fundamental_filename)
            if self.precomputed_grid_dir is not None and abundances_key == "" and (os.path.exists(precomputed_file) or os.path.exists(precomputed_step_file)):
                if not self.quiet:
                    print("Pre-computed:", complete_key)
                if os.path.exists(precomputed_file):
                    precomputed = read_spectrum(precomputed_file)
                else:
                    precomputed = read_spectrum(precomputed_step_file)
                convolved_precomputed = convolve_spectrum(precomputed, self.R(component=component))

                convolved_precomputed = resample_spectrum(convolved_precomputed, self.waveobs, method="linear", zero_edges=True)
                convolved_precomputed['flux'][self.waveobs_mask == 0] = 1.
                self.last_fluxes = convolved_precomputed['flux'].copy()
                self.last_final_fluxes[component] = convolved_precomputed['flux'].copy()

            elif self.precomputed_grid_dir is not None and abundances_key == "" and (os.path.exists(fundamental_precomputed_file) or os.path.exists(fundamental_precomputed_step_file)):
                if not self.quiet:
                    print("Pre-computed (fundamental):", complete_key)
                if os.path.exists(fundamental_precomputed_file):
                    fundamental_precomputed = read_spectrum(fundamental_precomputed_file)
                else:
                    fundamental_precomputed = read_spectrum(fundamental_precomputed_step_file)
                fundamental_precomputed = resample_spectrum(fundamental_precomputed, self.waveobs, method="linear", zero_edges=True)
                fundamental_precomputed['flux'][self.waveobs_mask == 0] = 1.
                self.last_fluxes = fundamental_precomputed['flux']
                # Optimization to avoid too small changes in parameters or repetition
                self.cache[key] = self.last_fluxes.copy()
                self.last_final_fluxes[component] = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, self.segments, macroturbulence=self.vmac(component=component), vsini=self.vsini(component=component), limb_darkening_coeff=self.limb_darkening_coeff(component=component), R=self.R(component=component), vrad=self.vrad(component=component), verbose=0)
                self.last_final_fluxes[component][self.waveobs_mask == 0] = 1.
            else:
                if key in self.cache:
                    if not self.quiet:
                        print("Cache:", complete_key)
                    self.last_fluxes = self.cache[key].copy()
                else:
                    if not self.quiet:
                        print("Generating:", complete_key)

                    if self.code != "grid":

                        # Atmosphere
                        atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, {'teff':self.teff(component=component), 'logg':self.logg(component=component), 'MH':self.MH(component=component), 'alpha':self.alpha(component=component)}, code=self.code)
                        # Fundamental synthetic fluxes
                        if self.code == "turbospectrum":
                            if len(self.nlte_originally_available) > 0:
                                # NLTE departure coefficient interpolation
                                nlte_departure_coefficients = interpolate_nlte_departure_coefficients(self.modeled_layers_pack, self.abundances, {'teff':self.teff(component=component), 'logg': self.logg(component=component), 'MH': self.MH(component=component), 'alpha': self.alpha(component=component)}, fixed_abundances=fixed_abundances, linelist=linelist, regions=self.segments, code=self.code)
                                if len(nlte_departure_coefficients) != len(self.nlte_originally_available):
                                    # Some elements may be ignored if there are no absorption lines in the considered regions
                                    self.nlte_available = list(nlte_departure_coefficients.keys())
                                    self.nlte_ignored = list(set(self.nlte_originally_available) - set(nlte_departure_coefficients.keys()))
                            else:
                                nlte_departure_coefficients = None
                            #
                            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), linelist, self.isotopes, self.abundances, fixed_abundances, self.vmic(component=component), atmosphere_layers_file=self.atmosphere_layers_files[component], abundances_file=self.abundances_files[component], linelist_file=self.linelist_file, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, use_molecules=self.use_molecules, nlte_departure_coefficients=nlte_departure_coefficients, tmp_dir=self.tmp_dir, timeout=self.timeout)
                        elif self.code in ("moog", "moog-scat"):
                            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), linelist, self.isotopes, self.abundances, fixed_abundances, self.vmic(component=component), atmosphere_layers_file=self.atmosphere_layers_files[component], abundances_file=self.abundances_files[component], linelist_file=self.linelist_file, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, tmp_dir=self.tmp_dir, timeout=self.timeout)
                        elif self.code == "synthe":
                            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), linelist, self.isotopes, self.abundances, fixed_abundances, self.vmic(component=component), atmosphere_layers_file=self.atmosphere_layers_files[component], abundances_file=self.abundances_files[component], linelist_file=self.linelist_file, molecules_files=self.molecules_files, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, tmp_dir=self.tmp_dir, timeout=self.timeout)
                        elif self.code == "sme":
                            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), linelist, self.isotopes, self.abundances, fixed_abundances, self.vmic(component=component), regions=self.segments, verbose=0, code=self.code, timeout=self.timeout)
                            ## Do not abort failed synthesis, the minimization algorithm will just consider this point as a bad one
                            #if np.all(self.last_fluxes == 0):
                                #raise Exception("SME has failed.")
                        elif self.code == "spectrum":
                            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), linelist, self.isotopes, self.abundances, fixed_abundances, self.vmic(component=component),  atmosphere_layers_file=self.atmosphere_layers_files[component], abundances_file=self.abundances_files[component], linelist_file=self.linelist_file, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, tmp_dir=self.tmp_dir, timeout=self.timeout)

                            ## Do not abort failed synthesis, the minimization algorithm will just consider this point as a bad one
                            #if np.all(self.last_fluxes == 0):
                                #raise Exception("SPECTRUM has failed.")
                        else:
                            raise Exception("Unknown code: %s" % (self.code))
                    else:
                        atmosphere_layers = None
                        linelist = None
                        abundances = None
                        fixed_abundances = None
                        self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), linelist, self.isotopes, abundances, fixed_abundances, self.vmic(component=component), code=self.code, grid=self.grid, regions=self.segments)

                    # Optimization to avoid too small changes in parameters or repetition
                    self.cache[key] = self.last_fluxes.copy()

                if not np.all(self.last_fluxes == 0):
                    # If synthesis did not fail
                    self.last_final_fluxes[component] = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, self.segments, macroturbulence=self.vmac(component=component), vsini=self.vsini(component=component), limb_darkening_coeff=self.limb_darkening_coeff(component=component), R=self.R(component=component), vrad=self.vrad(component=component), verbose=0)

                    if self.normalize_func is not None:
                        self.last_final_fluxes[component] = self.normalize_func(create_spectrum_structure(self.waveobs, self.last_final_fluxes[component]))['flux']
                else:
                    self.last_final_fluxes[component] = self.last_fluxes

                if (self.code != "grid" and model_atmosphere_is_closest_copy(self.modeled_layers_pack, {'teff':self.teff(component=component), 'logg':self.logg(component=component), 'MH':self.MH(component=component), 'alpha':self.alpha(component=component), 'vmic': self.vmic(component=component)})) \
                   or (self.code == "grid" and not valid_interpolated_spectrum_target(self.grid, {'teff':self.teff(component=component), 'logg':self.logg(component=component), 'MH':self.MH(component=component), 'alpha':self.alpha(component=component), 'vmic': self.vmic(component=component)})):
                       self.last_final_fluxes[component] *= 0.

        self.combined_final_fluxes = np.zeros_like(self.last_final_fluxes[0])
        for component in range(self.n_stellar_components):
            self.combined_final_fluxes += self.last_final_fluxes[component] * self.lf(component=component)

        return self.combined_final_fluxes[self.comparing_mask]

    def fitData(self, spectrum, segments, comparing_mask, weights, parinfo=None, parinfo_bases=None, use_errors=False, max_iterations=20, quiet=True, code="spectrum", use_molecules=False, vmic_from_empirical_relation=True, vmac_from_empirical_relation=True, n_stellar_components=1, tmp_dir=None, timeout=1800):
        code = code.lower()
        if code not in ['spectrum', 'turbospectrum', 'moog', 'moog-scat', 'synthe', 'sme', 'grid']:
            raise Exception("Unknown radiative transfer code: %s" % (code))

        self.timeout = timeout
        self.use_errors = use_errors

        base = PARINFO_BASE
        if len(parinfo) < base*n_stellar_components:
            raise Exception("Wrong number of parameters!")
        self.parinfo_bases = parinfo_bases
        self.n_stellar_components = n_stellar_components

        if sys.platform == "win32":
            # On Windows, the best timer is time.clock()
            default_timer = time.clock
        else:
            # On most other platforms the best timer is time.time()
            default_timer = time.time
        self.spectrum = spectrum
        if weights is None:
            weights = np.ones(len(spectrum['waveobs']))
        self.weights = weights
        self.code = code
        if self.code == "grid" and self.grid is None:
            self.grid = load_spectral_grid(self.precomputed_grid_dir)
        self.nlte_originally_available = []
        if self.code == "turbospectrum":
            nlte_dep_grid = self.modeled_layers_pack[8]
            if nlte_dep_grid is not None and len(nlte_dep_grid) > 0:
                self.nlte_originally_available = list(nlte_dep_grid.keys())
        self.nlte_available = self.nlte_originally_available
        self.nlte_ignored = []
        self.use_molecules = use_molecules
        self.vmic_from_empirical_relation = vmic_from_empirical_relation
        self.vmac_from_empirical_relation = vmac_from_empirical_relation
        self.tmp_dir = tmp_dir


        self.segments = segments
        self.waveobs = spectrum['waveobs']
        self.waveobs_mask = _create_waveobs_mask(spectrum['waveobs'], segments)
        self.comparing_mask = comparing_mask == 1.0 # Wavelengths to be compared for the least square algorithm

        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 0.0   # Do not limit residuals between -1.0 and 1.0 (np.tanh(residuals/1.0))
        _t0 = default_timer()

        # Write abundances and linelist to avoid writing the same info in each iteration
        for i in range(n_stellar_components):
            alpha_in_free_params = not parinfo[3+self.parinfo_bases[i]]['fixed']
            if self.enhance_abundances or alpha_in_free_params or self.abundances is None:
                self.abundances_files.append(None)
            else:
                self.abundances_files.append(write_solar_abundances(self.abundances))

        if len(self.linelist_free_loggf) == 0:
            # Only write linelist (for optimization purposes) if there is no free loggf
            if self.code == 'synthe':
                self.linelist_file, self.molecules_files = write_atomic_linelist(self.linelist, code="synthe", tmp_dir=tmp_dir)
            elif self.code == 'sme' or self.code in ('moog', 'moog-scat'):
                # moog requires two files for the linelist
                # sme does not require files
                self.linelist_file = None
            elif self.code == 'turbospectrum' or self.code == 'spectrum':
                self.linelist_file = write_atomic_linelist(self.linelist, code=self.code, tmp_dir=tmp_dir)

        if self.code == "spectrum":
            self.isotope_file = write_isotope_data(self.isotopes, tmp_dir=tmp_dir)

        for i in range(n_stellar_components):
            # If teff, logg, MH and alpha are fixed
            if self.code not in ('sme', 'grid') and parinfo[0+self.parinfo_bases[i]]['fixed'] and parinfo[1+self.parinfo_bases[i]]['fixed'] and parinfo[2+self.parinfo_bases[i]]['fixed'] and parinfo[3+self.parinfo_bases[i]]['fixed']:
                atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, {'teff':parinfo[0+self.parinfo_bases[i]]['value'], 'logg':parinfo[1+self.parinfo_bases[i]]['value'], 'MH':parinfo[2+self.parinfo_bases[i]]['value'], 'alpha':parinfo[3+self.parinfo_bases[i]]['value']})
                self.atmosphere_layers_files.append(write_atmosphere(atmosphere_layers, parinfo[0+self.parinfo_bases[i]]['value'], parinfo[1+self.parinfo_bases[i]]['value'], parinfo[2+self.parinfo_bases[i]]['value'], code=self.code, atmosphere_filename=None, tmp_dir=tmp_dir))
            else:
                self.atmosphere_layers_files.append(None)

        super(SynthModel, self).fitData(spectrum['waveobs'][self.comparing_mask], spectrum['flux'][self.comparing_mask], weights=weights[self.comparing_mask], parinfo=parinfo, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet)

        residuals = self.combined_final_fluxes[self.comparing_mask] - spectrum['flux'][self.comparing_mask]
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))

        #### Unweighted (no errors considered):
        self.chisq = np.sum((residuals)**2)
        self.reduced_chisq = self.chisq / self.m.dof

        #### Weighted (errors considered):
        self.wchisq = np.sum((weights[self.comparing_mask] * residuals)**2)
        self.reduced_wchisq = self.wchisq / self.m.dof

        self.cache = {}

        if len(self.linelist_free_loggf) == 0 and self.code not in ('sme', 'moog', 'moog-scat', 'grid'):
            os.remove(self.linelist_file)
        if self.code == 'spectrum':
            os.remove(self.isotope_file)
        if self.code == 'synthe' and self.molecules_files is not None:
            for molecules_file in self.molecules_files:
                os.remove(molecules_file)
        for i in range(n_stellar_components):
            if self.atmosphere_layers_files[i] is not None and os.path.exists(self.atmosphere_layers_files[i]):
                os.remove(self.atmosphere_layers_files[i])
        for i in range(n_stellar_components):
            if self.abundances_files[i] is not None and os.path.exists(self.abundances_files[i]):
                os.remove(self.abundances_files[i])
        self.abundances_files = []
        self.linelist_file = None
        self.isotope_file = None
        self.molecules_files = None
        self.atmosphere_layers_files = []

        _t1 = default_timer()
        sec = timedelta(seconds=int(_t1 - _t0))
        self.calculation_time = datetime(1,1,1) + sec

    def teff(self, component=0): return self._parinfo[0+self.parinfo_bases[component]]['value']
    def logg(self, component=0): return self._parinfo[1+self.parinfo_bases[component]]['value']
    def MH(self, component=0): return self._parinfo[2+self.parinfo_bases[component]]['value']
    def alpha(self, component=0):
        if self.enhance_abundances:
            alpha_enhancement = determine_abundance_enchancements(self.MH(component=component), scale=self.scale)
        else:
            alpha_enhancement = self._parinfo[3+self.parinfo_bases[component]]['value']
        return alpha_enhancement
    def vmic(self, component=0): return self._parinfo[4+self.parinfo_bases[component]]['value']
    def vmac(self, component=0): return self._parinfo[5+self.parinfo_bases[component]]['value']
    def vsini(self, component=0): return self._parinfo[6+self.parinfo_bases[component]]['value']
    def limb_darkening_coeff(self, component=0): return self._parinfo[7+self.parinfo_bases[component]]['value']
    def R(self, component=0): return self._parinfo[8+self.parinfo_bases[component]]['value']
    def lf(self, component=0):
        if self.n_stellar_components <= 2:
            # For single stars, there is nothing to do
            # For binary stars, companion's lf is tied to primary's lf via parinfo['tied']
            return self._parinfo[9+self.parinfo_bases[component]]['value']
        else:
            # Normalization is required
            total_lf = sum( [self._parinfo[9+self.parinfo_bases[c]]['value'] for c in range(self.n_stellar_components)])
            return self._parinfo[9+self.parinfo_bases[component]]['value'] / total_lf # ensure the final sum of light fraction equals to 1
    def vrad(self, component=0):
        vrad = []
        base = PARINFO_BASE + self.parinfo_bases[component]
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            top = base+len(self.segments)
            for i in range(base, top):
                vrad.append(self._parinfo[i]['value'])
        else:
            vrad = np.zeros(len(self.segments))
        return vrad
    def free_loggf(self, component=0):
        base = PARINFO_BASE + self.parinfo_bases[component]
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            base += len(self.segments)
        loggf = []
        for i in range(base, base+len(self.linelist_free_loggf)):
            loggf.append(self._parinfo[i]['value'])
        return loggf
    def free_abundances(self, component=0):
        base = PARINFO_BASE + self.parinfo_bases[component]
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            base += len(self.segments)
        base += len(self.linelist_free_loggf)
        # Find the top for this component to infer how many free abundances there are
        if component >= len(self.parinfo_bases)-1:
            top = len(self._parinfo)
        else:
            top = self.parinfo_bases[component+1]
        n_free_abundances = top - base
        fixed_abundances = np.recarray((n_free_abundances, ), dtype=[('code', int),('Abund', float), ('element', '|U30')])
        for i in range(n_free_abundances):
            fixed_abundances['code'][i] = int(self._parinfo[base+i]['parname'])
            fixed_abundances['Abund'][i] = self._parinfo[base+i]['value']
            fixed_abundances['element'][i] = ""
        return fixed_abundances

    def eteff(self, component=0): return self.m.perror[0+self.parinfo_bases[component]]
    def elogg(self, component=0): return self.m.perror[1+self.parinfo_bases[component]]
    def eMH(self, component=0): return self.m.perror[2+self.parinfo_bases[component]]
    def ealpha(self, component=0): return self.m.perror[3+self.parinfo_bases[component]]
    def evmic(self, component=0): return self.m.perror[4+self.parinfo_bases[component]]
    def evmac(self, component=0): return self.m.perror[5+self.parinfo_bases[component]]
    def evsini(self, component=0): return self.m.perror[6+self.parinfo_bases[component]]
    def elimb_darkening_coeff(self, component=0): return self.m.perror[7+self.parinfo_bases[component]]
    def eR(self, component=0): return self.m.perror[8+self.parinfo_bases[component]]
    def elf(self, component=0):
        if self.n_stellar_components <= 2:
            # For single stars, there is nothing to do
            # For binary stars, companion's lf is tied to primary's lf via parinfo['tied']
            return self.m.perror[9+self.parinfo_bases[component]]
        else:
            ## Because lf (light fraction) is normalized in df(), error needs to be propagated accordingly:
            # 1. Get raw fitted values and covariance matrix for the LF parameters
            lf_indices = [9 + self.parinfo_bases[c] for c in range(self.n_stellar_components)]
            raw_lfs = np.array([self._parinfo[i]['value'] for i in lf_indices])
            # Extract the relevant sub-matrix from the full covariance matrix
            full_covar = self.m.covar
            lf_covar = full_covar[np.ix_(lf_indices, lf_indices)]
            # 2. Calculate the total sum of raw light fractions
            S = np.sum(raw_lfs)
            # 3. Construct the Jacobian matrix (J_ik = dl_i / dL_k)
            J = np.zeros((self.n_stellar_components, self.n_stellar_components))
            for i in range(self.n_stellar_components):
                for k in range(self.n_stellar_components):
                    if i == k:
                        J[i, k] = (S - raw_lfs[i]) / (S**2)
                    else:
                        J[i, k] = -raw_lfs[i] / (S**2)

            # 4. Propagate the error
            #    Cov_new = J * Cov_old * J.T
            covar_new = J @ lf_covar @ J.T

            # 5. The error is the square root of the diagonal element
            variance_new = covar_new[component, component]

            return np.sqrt(variance_new)
    def evrad(self, component=0):
        base = PARINFO_BASE + self.parinfo_bases[component]
        evrad = []
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            top = base+len(self.segments)
            for i in range(base, top):
                evrad.append(self.m.perror[i])
        else:
            evrad = np.zeros(len(self.segments))
        return evrad
    def efree_loggf(self, component=0):
        base = PARINFO_BASE + self.parinfo_bases[component]
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            base += len(self.segments)
        eloggf = []
        for i in range(base, base+len(self.linelist_free_loggf)):
            eloggf.append(self.m.perror[i])
        return eloggf
    def efree_abundances(self, component=0):
        base = PARINFO_BASE + self.parinfo_bases[component]
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            base += len(self.segments)
        base += len(self.linelist_free_loggf)
        eabundances = []
        for i in range(len(self._parinfo)-base):
            eabundances.append(self.m.perror[base+i])
        return eabundances

    def generate_linelist_free_loggf(self, component=0):
        linelist_free_loggf = self.linelist_free_loggf.copy()
        new_loggf = self.free_loggf(component=component)
        for i in range(len(linelist_free_loggf)):
            linelist_free_loggf['loggf'][i] = new_loggf[i]
        return linelist_free_loggf

    def transformed_free_loggf(self, component=0):
        free_loggf = {}
        free_loggf['linelist'] = self.generate_linelist_free_loggf(component=component)
        free_loggf['loggf'] = self.free_loggf(component=component)
        free_loggf['eloggf'] = self.efree_loggf(component=component)
        return free_loggf

    def transformed_free_abundances(self, component=0):
        free_abundances = self.free_abundances(component=component)
        efree_abundances = self.efree_abundances(component=component)

        transformed_abund = np.recarray((len(free_abundances), ), dtype=[('code', int),('Abund', float), ('element', '|U5'), ('[X/H]', float), ('A(X)', float), ('[X/Fe]', float), ('eAbund', float), ('e[X/H]', float), ('e[X/Fe]', float), ('eA(X)', float)])

        #function abundances, sme, solar_relative
        #abund = sme.abund
        #solar_abund, solar, labels

        #solar_relative = abund-solar

        #abund[2:*] += sme.feh ; rescale by metallicity
        #abund[1:*] = 10^abund[1:*] ; transform to linear
        #return, alog10(abund / abund[0]) + 12
        #end
        sun_log_Nh_over_Ntotal = self.abundances['Abund'][self.abundances['code'] == 1]
        for i in range(len(free_abundances)):
            sun_log_Nx_over_Ntotal = self.abundances['Abund'][self.abundances['code'] == free_abundances['code'][i]]
            x_absolute = free_abundances['Abund'][i] + 12. - sun_log_Nh_over_Ntotal # absolute, A(X)
            #x_over_fe = free_abundances['Abund'][i] - sun_log_Nx_over_Ntotal
            #x_over_h = x_over_fe + self.MH()
            x_over_h = free_abundances['Abund'][i] - sun_log_Nx_over_Ntotal
            x_over_fe = x_over_h - self.MH()

            element = self.elements[str(free_abundances['code'][i])]

            transformed_abund['code'][i] = free_abundances['code'][i]
            transformed_abund['Abund'][i] = free_abundances['Abund'][i]
            transformed_abund['element'][i] = element
            transformed_abund['[X/H]'][i] = x_over_h
            transformed_abund['[X/Fe]'][i] = x_over_fe
            transformed_abund['A(X)'][i] = x_absolute
            transformed_abund['eAbund'][i] = efree_abundances[i]
            transformed_abund['e[X/H]'][i] = efree_abundances[i]
            transformed_abund['e[X/Fe]'][i] = efree_abundances[i]
            transformed_abund['eA(X)'][i] = efree_abundances[i]
        return transformed_abund

    def print_solution(self, component=0):
        if self.use_errors:
            #error_scale_factor = self.reduced_wchisq
            #error_scale_factor = self.wchisq
            error_scale_factor = 1.
            # https://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html
            # https://www.gnu.org/software/gsl/manual/gsl-ref_38.html
            #error_scale_factor = np.max((1., self.wchisq/np.sqrt(self.m.dof)))
        else:
            error_scale_factor = 1.
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","alpha","vmic","vmac","vsini","limb","R")
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i" % (self.teff(component=component), self.logg(component=component), self.MH(component=component), self.alpha(component=component), self.vmic(component=component), self.vmac(component=component), self.vsini(component=component), self.limb_darkening_coeff(component=component), int(self.R(component=component)))
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i" % (self.eteff(component=component)*error_scale_factor, self.elogg(component=component)*error_scale_factor, self.eMH(component=component)*error_scale_factor, self.ealpha(component=component)*error_scale_factor, self.evmic(component=component)*error_scale_factor, self.evmac(component=component)*error_scale_factor, self.evsini(component=component)*error_scale_factor, self.elimb_darkening_coeff(component=component)*error_scale_factor, int(self.eR(component=component)*error_scale_factor))

        # Append free individual abundances
        abundances_header = ""
        abundances_solution = ""
        abundances_errors = ""


        if self.code != "grid":
            transformed_abund = self.transformed_free_abundances(component=component)
            for i in range(len(transformed_abund)):
                element = transformed_abund['element'][i]
                x_absolute_name = "A(" + element + ")"
                x_over_h_name = "[" + element + "/H]"
                x_over_fe_name = "[" + element + "/Fe]"
                x = transformed_abund['Abund'][i]
                x_absolute = transformed_abund['A(X)'][i]
                x_over_h = transformed_abund['[X/H]'][i]
                x_over_fe = transformed_abund['[X/Fe]'][i]
                ex = transformed_abund['eAbund'][i]*error_scale_factor
                ex_absolute = transformed_abund['eA(X)'][i]*error_scale_factor
                ex_over_h = transformed_abund['e[X/H]'][i]*error_scale_factor
                ex_over_fe = transformed_abund['e[X/Fe]'][i]*error_scale_factor
                abundances_header += "%8s\t%8s\t%8s\t%8s" % (element, x_absolute_name, x_over_h_name, x_over_fe_name)
                abundances_solution += "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (x, x_absolute, x_over_h, x_over_fe)
                abundances_errors += "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (ex, ex_absolute, ex_over_h, ex_over_fe)
        else:
            transformed_abund = []

        base = PARINFO_BASE + self.parinfo_bases[component]
        component_with_vrad_parname = len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']
        if component_with_vrad_parname:
            top = base+len(self.segments)
            component_with_free_vrad = any([not p['fixed']for p in self._parinfo[base:top]])
            if component_with_free_vrad:
                vrad_header = "          %8s\t%8s\t%8s\t%8s" % ("wave_base","wave_top","vrad","error")
                vrad_stats = ""
                for i, (vrad, evrad, segment) in enumerate(zip(self.vrad(component=component), self.evrad(component=component), self.segments)):
                    vrad_stats += "Segment   %8.2f\t%8.2f\t%8.2f\t%8.2f\n" % (segment['wave_base'], segment['wave_top'], vrad, evrad)
                print("")
                print(vrad_header)
                print(vrad_stats)
                print("")


        print("           ", header)
        print("Solution:  ", solution)
        print("Errors:    ", errors)
        print("")
        if len(transformed_abund) > 0:
            print("           ", abundances_header)
            print("Abundances:", abundances_solution)
            print("Ab. errors:", abundances_errors)
            print("")

        if self.code != "grid":
            transformed_free_loggf = self.transformed_free_loggf(component=component)
            transformed_linelist_free_loggf = transformed_free_loggf['linelist']
            loggf = transformed_free_loggf['loggf']
            eloggf = transformed_free_loggf['eloggf']
            if len(transformed_linelist_free_loggf) > 0:
                loggf_header = "          %8s\t%8s\t%8s\t%8s" % ("wave_base","wave_top","log(gf)","error")
                loggf_stats = ""
                for i in range(len(transformed_linelist_free_loggf)):
                    loggf_stats += "log(gf)     %8.4f\t%8s\t%8.3f\t%8.3f\n" % (transformed_linelist_free_loggf['wave_nm'][i], transformed_linelist_free_loggf['element'][i], loggf[i], eloggf[i])
                print("")
                print(loggf_header)
                print(loggf_stats)
                print("")

        print("Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second))
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","wchisq","rwchisq","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (self.m.dof, self.m.niter, self.m.nfev, self.wchisq, self.reduced_wchisq, self.chisq, self.reduced_chisq, self.rms)
        if self.code != "grid" and model_atmosphere_is_closest_copy(self.modeled_layers_pack, {'teff':self.teff(component=component), 'logg':self.logg(component=component), 'MH':self.MH(component=component), 'alpha':self.alpha(component=component), 'vmic': self.vmic(component=component)}):
            print("")
            print("WARNING: Model atmosphere used for the final solution was not interpolated, it is a copy of the closest model.")
        if self.code == "grid" and not valid_interpolated_spectrum_target(self.grid, {'teff':self.teff(component=component), 'logg':self.logg(component=component), 'MH':self.MH(component=component), 'alpha':self.alpha(component=component), 'vmic': self.vmic(component=component)}):
            print("")
            print("WARNING: Spectrum used for the final solution was not interpolated, it is a copy of the closest model.")
        print("")
        print("         ", header)
        print("Stats:   ", stats)
        print("Return code:", self.m.status)




def __create_param_structure(initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_lf, initial_vrad, free_params, free_abundances, linelist_free_loggf, teff_range, logg_range, MH_range, alpha_range, vmic_range, vmic_from_empirical_relation, vmac_from_empirical_relation):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = PARINFO_BASE
    free_params = [param.lower() for param in free_params]
    if "vrad" in free_params or np.any(initial_vrad != 0):
        parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0, 'mpmaxstep':0} for i in np.arange(base+len(initial_vrad)+len(free_abundances)+len(linelist_free_loggf))]
    else:
        parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0, 'mpmaxstep':0} for i in np.arange(base+len(free_abundances)+len(linelist_free_loggf))]
    ##
    min_teff = np.min(teff_range)
    max_teff = np.max(teff_range)
    #
    parinfo[0]['parname'] = "teff"
    parinfo[0]['value'] = initial_teff
    parinfo[0]['fixed'] = not parinfo[0]['parname'].lower() in free_params
    parinfo[0]['step'] = Constants.SYNTH_STEP_TEFF # For auto-derivatives
    parinfo[0]['limited'] = [True, True]
    parinfo[0]['limits'] = [min_teff, max_teff]
    if parinfo[0]['value'] > parinfo[0]['limits'][1] or parinfo[0]['value'] < parinfo[0]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[0]['parname'], parinfo[0]['value'], parinfo[0]['limits'][0], parinfo[0]['limits'][1]))
    ##
    min_logg = np.min(logg_range)
    max_logg = np.max(logg_range)
    #
    parinfo[1]['parname'] = "logg"
    parinfo[1]['value'] = initial_logg
    parinfo[1]['fixed'] = not parinfo[1]['parname'].lower() in free_params
    parinfo[1]['step'] = Constants.SYNTH_STEP_LOGG # For auto-derivatives
    #if not parinfo[1]['fixed']:
        #parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [min_logg, max_logg]
    if parinfo[1]['value'] > parinfo[1]['limits'][1] or parinfo[1]['value'] < parinfo[1]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[1]['parname'], parinfo[1]['value'], parinfo[1]['limits'][0], parinfo[1]['limits'][1]))
    ##
    min_MH = np.min(MH_range)
    max_MH = np.max(MH_range)
    #
    parinfo[2]['parname'] = "MH"
    parinfo[2]['value'] = initial_MH
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = Constants.SYNTH_STEP_MH # For auto-derivatives
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [min_MH, max_MH]
    if parinfo[2]['value'] > parinfo[2]['limits'][1] or parinfo[2]['value'] < parinfo[2]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[2]['parname'], parinfo[2]['value'], parinfo[2]['limits'][0], parinfo[2]['limits'][1]))
    #
    min_alpha = np.min(alpha_range)
    max_alpha = np.max(alpha_range)
    #
    parinfo[3]['parname'] = "alpha"
    parinfo[3]['value'] = initial_alpha
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = Constants.SYNTH_STEP_ALPHA # For auto-derivatives
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [min_alpha, max_alpha]
    if parinfo[3]['value'] > parinfo[3]['limits'][1] or parinfo[3]['value'] < parinfo[3]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[3]['parname'], parinfo[3]['value'], parinfo[3]['limits'][0], parinfo[3]['limits'][1]))
    #
    min_vmic = np.min(vmic_range)
    max_vmic = np.max(vmic_range)
    #
    parinfo[4]['parname'] = "Vmic"
    parinfo[4]['value'] = initial_vmic
    parinfo[4]['fixed'] = not parinfo[4]['parname'].lower() in free_params
    if vmic_from_empirical_relation:
        parinfo[4]['tied'] = 'estimate_vmic(p[0], p[1], p[2])'
    parinfo[4]['step'] = Constants.SYNTH_STEP_VMIC # For auto-derivatives
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [min_vmic, max_vmic]
    if parinfo[4]['value'] > parinfo[4]['limits'][1] or parinfo[4]['value'] < parinfo[4]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[4]['parname'], parinfo[4]['value'], parinfo[4]['limits'][0], parinfo[4]['limits'][1]))
    #
    parinfo[5]['parname'] = "Vmac"
    parinfo[5]['value'] = initial_vmac
    parinfo[5]['fixed'] = not parinfo[5]['parname'].lower() in free_params
    if vmac_from_empirical_relation:
        parinfo[5]['tied'] = 'estimate_vmac(p[0], p[1], p[2])'
    parinfo[5]['step'] = Constants.SYNTH_STEP_VMAC # For auto-derivatives
    parinfo[5]['limited'] = [True, True]
    parinfo[5]['limits'] = [0.0, 50.0]
    if parinfo[5]['value'] > parinfo[5]['limits'][1] or parinfo[5]['value'] < parinfo[5]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[5]['parname'], parinfo[5]['value'], parinfo[5]['limits'][0], parinfo[5]['limits'][1]))
    #
    parinfo[6]['parname'] = "Vsini"
    parinfo[6]['value'] = initial_vsini
    parinfo[6]['fixed'] = not parinfo[6]['parname'].lower() in free_params
    parinfo[6]['step'] = Constants.SYNTH_STEP_VSINI # For auto-derivatives
    parinfo[6]['limited'] = [True, True]
    parinfo[6]['limits'] = [0.0, 300.0]
    if parinfo[6]['value'] > parinfo[6]['limits'][1] or parinfo[6]['value'] < parinfo[6]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[6]['parname'], parinfo[6]['value'], parinfo[6]['limits'][0], parinfo[6]['limits'][1]))
    #
    parinfo[7]['parname'] = "limb_darkening_coeff"
    parinfo[7]['value'] = initial_limb_darkening_coeff
    parinfo[7]['fixed'] = not parinfo[7]['parname'].lower() in free_params
    parinfo[7]['step'] = Constants.SYNTH_STEP_LIMB_DARKENING_COEFF # For auto-derivatives
    parinfo[7]['limited'] = [True, True]
    parinfo[7]['limits'] = [0.0, 1.0]
    if parinfo[7]['value'] > parinfo[7]['limits'][1] or parinfo[7]['value'] < parinfo[7]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[7]['parname'], parinfo[7]['value'], parinfo[7]['limits'][0], parinfo[7]['limits'][1]))
    #
    parinfo[8]['parname'] = "R"
    parinfo[8]['value'] = initial_R
    parinfo[8]['fixed'] = not parinfo[8]['parname'].lower() in free_params
    parinfo[8]['step'] = Constants.SYNTH_STEP_R # For auto-derivatives
    parinfo[8]['limited'] = [True, True]
    parinfo[8]['limits'] = [100.0, 900000.0]
    if not parinfo[8]['fixed']:
        parinfo[8]['mpmaxstep'] = float(initial_R)/4. # Maximum change to be made in the parameter
    if parinfo[8]['value'] > parinfo[8]['limits'][1] or parinfo[8]['value'] < parinfo[8]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[8]['parname'], parinfo[8]['value'], parinfo[8]['limits'][0], parinfo[8]['limits'][1]))
    #
    parinfo[9]['parname'] = "lf" # Light Fraction (useful for binaries, ...)
    parinfo[9]['value'] = initial_lf
    parinfo[9]['fixed'] = not parinfo[9]['parname'].lower() in free_params
    parinfo[9]['mpprint'] = not parinfo[9]['fixed'] # Do not print if it is fixed
    parinfo[9]['step'] = Constants.SYNTH_STEP_LF # For auto-derivatives
    parinfo[9]['limited'] = [True, True]
    parinfo[9]['limits'] = [1e-6, 1.0]
    if not parinfo[9]['fixed']:
        parinfo[9]['mpmaxstep'] = 0.25 # Maximum change to be made in the parameter
    if parinfo[9]['value'] > parinfo[9]['limits'][1] or parinfo[9]['value'] < parinfo[9]['limits'][0]:
        raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[9]['parname'], parinfo[9]['value'], parinfo[9]['limits'][0], parinfo[9]['limits'][1]))
    # VRAD
    if "vrad" in free_params or np.any(initial_vrad != 0):
        base = PARINFO_BASE
        for i in range(len(initial_vrad)):
            parinfo[base+i]['parname'] = "vrad%03i" % (i)
            parinfo[base+i]['value'] = initial_vrad[i]
            parinfo[base+i]['fixed'] = not "vrad" in free_params
            parinfo[base+i]['mpprint'] = not parinfo[base+i]['fixed'] # Do not print if it is fixed
            parinfo[base+i]['step'] = Constants.SYNTH_STEP_VRAD # For auto-derivatives
            parinfo[base+i]['limited'] = [True, True]
            vrad_threshold = 5.0
            parinfo[base+i]['limits'] = [initial_vrad[i]-vrad_threshold, initial_vrad[i]+vrad_threshold]
            if parinfo[base+i]['value'] > parinfo[base+i]['limits'][1] or parinfo[base+i]['value'] < parinfo[base+i]['limits'][0]:
                raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[base+i]['parname'], parinfo[base+i]['value'], parinfo[base+i]['limits'][0], parinfo[base+i]['limits'][1]))
    # log(gf)
    if "vrad" in free_params or np.any(initial_vrad != 0):
        base = PARINFO_BASE + len(initial_vrad)
    else:
        base = PARINFO_BASE
    for i in range(len(linelist_free_loggf)):
        parinfo[base+i]['parname'] = str(linelist_free_loggf['wave_nm'][i])
        parinfo[base+i]['value'] = linelist_free_loggf['loggf'][i]
        parinfo[base+i]['fixed'] = not "loggf" in free_params
        parinfo[base+i]['step'] = Constants.SYNTH_STEP_LOGGF # For auto-derivatives
        parinfo[base+i]['limited'] = [True, True]
        parinfo[base+i]['limits'] = [-10., 10.]
        if parinfo[base+i]['value'] > parinfo[base+i]['limits'][1] or parinfo[base+i]['value'] < parinfo[base+i]['limits'][0]:
            raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[base+i]['parname'], parinfo[base+i]['value'], parinfo[base+i]['limits'][0], parinfo[base+i]['limits'][1]))
    if "vrad" in free_params or np.any(initial_vrad != 0):
        base = PARINFO_BASE + len(initial_vrad) + len(linelist_free_loggf)
    else:
        base = PARINFO_BASE + len(linelist_free_loggf)
    # ABUNDANCES
    for i in range(len(free_abundances)):
        parinfo[base+i]['parname'] = str(free_abundances['code'][i])
        parinfo[base+i]['value'] = free_abundances['Abund'][i]
        parinfo[base+i]['fixed'] = not parinfo[base+i]['parname'].lower() in free_params
        parinfo[base+i]['step'] = Constants.SYNTH_STEP_ABUNDANCES # For auto-derivatives
        parinfo[base+i]['limited'] = [True, True]
        parinfo[base+i]['limits'] = [-30., 0.]
        if parinfo[base+i]['value'] > parinfo[base+i]['limits'][1] or parinfo[base+i]['value'] < parinfo[base+i]['limits'][0]:
            raise Exception("Initial {} '{}' is out of range: '{}' - '{}'".format(parinfo[base+i]['parname'], parinfo[base+i]['value'], parinfo[base+i]['limits'][0], parinfo[base+i]['limits'][1]))

    return parinfo

def model_spectrum(spectrum, continuum_model, modeled_layers_pack, linelist, isotopes, abundances, free_abundances, linelist_free_loggf, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=None, linemasks=None, enhance_abundances=False, scale=None, precomputed_grid_dir=None, use_errors=True, max_iterations=20, verbose=1, code="spectrum", grid=None, use_molecules=False, vmic_from_empirical_relation=False, vmac_from_empirical_relation=False, initial_lf=1, normalize_func=None, tmp_dir=None, timeout=1800):

    """
    It matches synthetic spectrum to observed spectrum by applying a least
    square algorithm.

    - free_params is an array that can contain any combination of the following
      strings: ["teff", "logg", "MH", "alpha", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    - free_abundances can be set to 'None'
    - If segments are specified, the synthetic spectrum will be only generated for
      those regions.
    - If linemasks are specified, only those regions will be used for comparison.
    - It does not compare negative or zero fluxes
    - If enhance_abundances is True, alpha elements and CNO abundances will be scaled
      depending on the metallicity. If scale is None, by default the standard
      MARCS composition will be used (recommended).

    * timeout is for single synthesis execution and not for the whole minimization.
    """

    if verbose or verbose == 1:
        verbose = True
        quiet = False
    else:
        verbose = False
        quiet = True

    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'moog-scat', 'synthe', 'sme', 'grid']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    # Convert single values to a list with a single element
    if not isinstance(initial_teff, (list, tuple)): initial_teff = [initial_teff]
    if not isinstance(initial_logg, (list, tuple)): initial_logg = [initial_logg]
    if not isinstance(initial_MH, (list, tuple)): initial_MH = [initial_MH]
    if not isinstance(initial_alpha, (list, tuple)): initial_alpha = [initial_alpha]
    if not isinstance(initial_vmic, (list, tuple)): initial_vmic = [initial_vmic]
    if not isinstance(initial_vmac, (list, tuple)): initial_vmac = [initial_vmac]
    if not isinstance(initial_vsini, (list, tuple)): initial_vsini = [initial_vsini]
    if not isinstance(initial_limb_darkening_coeff, (list, tuple)): initial_limb_darkening_coeff = [initial_limb_darkening_coeff]
    if not isinstance(initial_R, (list, tuple)): initial_R = [initial_R]
    if not isinstance(initial_vrad, (list, tuple)): initial_vrad = [initial_vrad]
    if not isinstance(initial_lf, (list, tuple)): initial_lf = [initial_lf]
    if not isinstance(segments, (list, tuple)): segments = [segments]
    if not isinstance(linemasks, (list, tuple)): linemasks = [linemasks]

    #--------------------------------------------------------------------------------
    # All input variables need to have the same number of elements (i.e., for a spectroscopic binary, we need two initial_teff, two initial_logg...)
    names = [
        'initial_teff', 'initial_logg', 'initial_MH', 'initial_alpha',
        'initial_vmic', 'initial_vmac', 'initial_vsini',
        'initial_limb_darkening_coeff', 'initial_R', 'initial_vrad',
        'initial_lf', 'segments', 'linemasks'
    ]

    # Fetch the local variables by name
    lengths = {name: len(locals()[name]) for name in names}

    if len(set(lengths.values())) != 1:
        raise ValueError(f"Inconsistent lengths among inputs: {lengths}")
    #--------------------------------------------------------------------------------

    n_stellar_components = len(initial_teff)
    if n_stellar_components == 1:
        # For single stars, these assumptions need to be respected:
        assert initial_vrad[0] == 0
        assert initial_lf[0] == 1
    assert sum(initial_lf) == 1, "light fractions should sum up to one"

    #--------------------------------------------------------------------------------
    # Pre-process the segments list to replace all `None` values
    # 1. Define a default segment that covers the entire spectrum range
    full_spectrum_segment = np.recarray((1,), dtype=[('wave_base', float), ('wave_top', float)])
    full_spectrum_segment['wave_base'][0] = np.min(spectrum['waveobs'])
    full_spectrum_segment['wave_top'][0] = np.max(spectrum['waveobs'])
    # 2. If a segment is None, use the default full_spectrum_segment. Otherwise, use the existing segment.
    segments = [s if s is not None else full_spectrum_segment for s in segments]
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # Pre-process the linemasks list to replace all `None` values AND discard linemasks not within a segment
    for i, component_segments in enumerate(segments):
        if linemasks[i] is None:
            # If there are no linemasks, then compare all fluxes in the segments for this component
            component_linemasks = np.recarray((len(component_segments),),  dtype=[('wave_peak', float), ('wave_base', float), ('wave_top', float)])
            component_linemasks['wave_base'] = component_segments['wave_base']
            component_linemasks['wave_top'] = component_segments['wave_top']
            component_linemasks['wave_peak'] = (component_linemasks['wave_base'] + component_linemasks['wave_top']) / 2
            linemasks[i] = component_linemasks
        # Remove linemasks that are outside the segments
        component_linemasks = linemasks[i]
        linemasks[i] = _filter_linemasks_not_in_segments(component_linemasks, component_segments)
    #--------------------------------------------------------------------------------

    global_segments = None
    global_linemasks = None
    for i, component_segments in enumerate(segments):
        if global_segments is None:
            global_segments = component_segments
        else:
            global_segments = np.hstack((global_segments, component_segments))
        if global_linemasks is None:
            global_linemasks = component_linemasks
        else:
            global_linemasks = np.hstack((global_linemasks, component_linemasks))
    #
    global_segments = merge_overlapping_regions(global_segments)
    global_linemasks = merge_overlapping_regions(global_linemasks)

    # Wavelengths to be computed: component_segments
    wfilter = create_wavelength_filter(spectrum, regions=global_segments)
    filtered_spectrum = create_spectrum_structure(spectrum['waveobs'][wfilter], spectrum['flux'][wfilter], spectrum['err'][wfilter])
    #
    filtered_spectrum = normalize_spectrum(spectrum, continuum_model)

    if free_abundances is None:
        # No free abundances
        free_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float), ('element', '|U30')])
    else:
        if code == "grid":
            raise Exception("There cannot be free abundances when using a spectral grid for interpolation")
        # Add free abundances as free params
        for element in free_abundances:
            free_params.append(str(element['code']))

    atomic_dtype = _get_atomic_linelist_definition()
    if linelist_free_loggf is None:
        linelist_free_loggf = np.recarray((0, ), dtype=atomic_dtype)
    else:
        if code == "grid":
            raise Exception("There cannot be free log(gf) when using a spectral grid for interpolation")
        # Make sure both linelists have the same structure so that it will be possible to merge them
        linelist = linelist[[x[0] for x in atomic_dtype]]
        linelist_free_loggf = linelist_free_loggf[[x[0] for x in atomic_dtype]]
        # Add free loggf as free params
        free_params.append("loggf")

    if code == "grid" and precomputed_grid_dir is None and grid is None:
        raise Exception("Pre-computed grid should be specified when using 'grid' code")

    if code == "grid":
        if grid is None:
            grid = load_spectral_grid(precomputed_grid_dir)
        existing_points, grid_free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname = grid
        if "teff" not in grid_free_parameters and "teff" in free_params:
            raise Exception("teff cannot be free when using this grid")
        if "logg" not in grid_free_parameters and "logg" in free_params:
            raise Exception("logg cannot be free when using this grid")
        if "MH" not in grid_free_parameters and "MH" in free_params:
            raise Exception("MH cannot be free when using this grid")
        if "alpha" not in grid_free_parameters and ("alpha" in free_params or enhance_abundances):
            raise Exception("alpha cannot be free or be automatically enhance when using this grid")
        if "vmic" not in grid_free_parameters and ("vmic" in free_params or vmic_from_empirical_relation):
            raise Exception("vmic cannot be free or follow an empirical relation when using this grid")
    else:
        ranges = modeled_layers_pack[7]
    teff_range = ranges['teff']
    logg_range = ranges['logg']
    MH_range = ranges['MH']
    if code == "grid":
        alpha_range = ranges.get('alpha', (0.,)) # Default (0.,) if 'alpha' is not a free parameter for grid interpolation
        vmic_range = ranges.get('vmic', (0.,)) # Default (0.,) if 'vmic' is not a free parameter for grid interpolation
    else:
        alpha_range = ranges.get('alpha', (-1.5, 1.5)) # Default (0.,) if 'alpha' is not a free parameter for atmosphere interpolation
        vmic_range = ranges.get('vmic', (0.0, 50.)) # Default (0.,) if 'vmic' is not a free parameter for atmosphere interpolation


    if "alpha" in free_params and enhance_abundances:
        enhance_abundances = False
        logging.warning("'enhance_abundances' changed to False because alpha is a free parameter")

    if "vmic" in free_params and vmic_from_empirical_relation:
        vmic_from_empirical_relation = False
        logging.warning("'vmic_from_empirical_relation' changed to False because vmic is a free parameter")

    if "vmac" in free_params and vmac_from_empirical_relation:
        vmac_from_empirical_relation = False
        logging.warning("'vmac_from_empirical_relation' changed to False because vmac is a free parameter")


    #--------------------------------------------------------------------------------
    # Limit linelist for this component (unless linelist is not defined, thus using a pre-computed grid)
    if linelist is not None:
        linelist = _filter_linelist(linelist, global_segments)
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # Compare fluxes inside line masks that belong to a segment
    comparing_mask = _create_comparing_mask(filtered_spectrum['waveobs'], global_linemasks, global_segments)

    ## Fluxes
    negative_zero_flux = filtered_spectrum['flux'] <= 0.0
    bad_fluxes = np.logical_and(comparing_mask == 1, negative_zero_flux)
    num_bad_fluxes = len(np.where(bad_fluxes)[0])
    # Do not compare negative or zero fluxes
    if num_bad_fluxes > 0:
        logging.warning(f"{num_bad_fluxes} fluxes have been discarded because they are negative or zero")
        comparing_mask[negative_zero_flux] = 0.0

    ## Errors
    if use_errors and np.all(filtered_spectrum['err'][comparing_mask == 1] <= 0):
        logging.warning(f"Use of errors has been desactivated because all of them are set to zero")
        use_errors = False

    negative_zero_err = filtered_spectrum['err'] <= 0.0
    bad_errors = np.logical_and(comparing_mask == 1, negative_zero_err)
    num_bad_errors = len(np.where(bad_errors)[0])
    ## Do not compare negative or zero errors
    if use_errors and num_bad_errors > 0:
        logging.warning(f"{num_bad_errors} fluxes have been discarded because their ERRORS are negative or zero")
        comparing_mask[negative_zero_err] = 0.0

    if np.all(comparing_mask == 0):
        logging.error(f"No fluxes left to be compared!")
        raise Exception(f"No fluxes left to be compared!")

    accept_weights = np.logical_and(comparing_mask == 1., np.logical_not(negative_zero_err))
    #--------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    weights = np.ones(len(filtered_spectrum['waveobs']))
    weights[accept_weights] = 1. / filtered_spectrum['err'][accept_weights]
    weights[~accept_weights] = 0.
    #weights /= np.sum(weights) # Normalize
    #weights /= np.max(weights) # Normalize
    #weights *= len(np.where(accept_weights)[0]) # Scale to number of points
    #weights *= 10000
    weights = np.sqrt(weights)  # When squaring the flux errors, we get more reasonable parameter's errors (empirically validated)
    #--------------------------------------------------------------------------------

    # Convert initial_vrad from one per component to one per global segment
    initial_vrad = [np.ones(len(global_segments))*i_vrad for i_vrad in initial_vrad]

    parinfos = []
    parinfo_bases = []
    for i_teff, i_logg, i_MH, i_alpha, i_vmic, i_vmac, i_vsini, i_limb_darkening_coeff, i_R, i_lf, i_vrad in zip(initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_lf, initial_vrad):
        parinfo = __create_param_structure(i_teff, i_logg, i_MH, i_alpha, i_vmic, i_vmac, i_vsini, i_limb_darkening_coeff, i_R, i_lf, i_vrad, free_params, free_abundances, linelist_free_loggf, teff_range, logg_range, MH_range, alpha_range, vmic_range, vmic_from_empirical_relation, vmac_from_empirical_relation)
        parinfo_bases.append(len(parinfos))
        parinfos += parinfo # concatenate

    if n_stellar_components == 2:
        lf_indices = [9 + parinfo_bases[c] for c in range(n_stellar_components)]
        parinfos[lf_indices[1]]['tied'] = f'1 - p[{lf_indices[0]}]'
        parinfos[lf_indices[0]]['limits'] = [1e-6, 1.0-1e-6]
        parinfos[lf_indices[1]]['limits'] = [1e-6, 1.0-1e-6]

    synth_model = SynthModel(modeled_layers_pack, linelist, isotopes, linelist_free_loggf, abundances, enhance_abundances=enhance_abundances, scale=scale, precomputed_grid_dir=precomputed_grid_dir, grid=grid, normalize_func=normalize_func)

    #segments = None
    synth_model.fitData(filtered_spectrum, global_segments, comparing_mask, weights, parinfo=parinfos, parinfo_bases=parinfo_bases, use_errors=use_errors, max_iterations=max_iterations, quiet=quiet, code=code, use_molecules=use_molecules, vmic_from_empirical_relation=vmic_from_empirical_relation, vmac_from_empirical_relation=vmac_from_empirical_relation, n_stellar_components=n_stellar_components, tmp_dir=tmp_dir, timeout=timeout)

    if verbose:
        print("\n")

    if linemasks is not None:
        stats_linemasks = _get_stats_per_linemask(filtered_spectrum['waveobs'], filtered_spectrum['flux'], synth_model.combined_final_fluxes, weights, free_params, global_linemasks, verbose=verbose)
    else:
        stats_linemasks = None

    if verbose:
        print("\n")
        for component in range(n_stellar_components):
            if n_stellar_components > 1:
                print(f"[Component {component+1}/{n_stellar_components}] Light Fraction: {synth_model.lf(component=component):.2f} +/- {synth_model.elf(component=component):.2f}")
            synth_model.print_solution(component=component)

    if len(synth_model.nlte_available) > 0:
        if len(synth_model.nlte_ignored) > 0:
            print(f"NLTE for elements: {synth_model.nlte_available} | Not in region: {synth_model.nlte_ignored}")
        else:
            print(f"NLTE available for elements: {synth_model.nlte_available}")
    else:
        print("Only LTE")

    # Collect information to be returned
    all_params = []
    all_errors = []
    all_free_abundances = []
    all_free_loggf = []
    for i in range(n_stellar_components):

        params = {}
        params['teff'] = synth_model.teff()
        params['logg'] = synth_model.logg()
        params['MH'] = synth_model.MH()
        params['alpha'] = synth_model.alpha()
        params['vmic'] = synth_model.vmic()
        params['vmac'] = synth_model.vmac()
        params['vsini'] = synth_model.vsini()
        params['limb_darkening_coeff'] = synth_model.limb_darkening_coeff()
        params['R'] = synth_model.R()
        params['lf'] = synth_model.lf()
        for i, vrad in enumerate(synth_model.vrad()):
            params['vrad%04i' % (i)] = vrad

        errors = {}
        errors['teff'] = synth_model.eteff()
        errors['logg'] = synth_model.elogg()
        errors['MH'] = synth_model.eMH()
        errors['alpha'] = synth_model.ealpha()
        errors['vmic'] = synth_model.evmic()
        errors['vmac'] = synth_model.evmac()
        errors['vsini'] = synth_model.evsini()
        errors['limb_darkening_coeff'] = synth_model.elimb_darkening_coeff()
        errors['R'] = synth_model.eR()
        errors['lf'] = synth_model.elf()
        for i, evrad in enumerate(synth_model.evrad()):
            errors['vrad%04i' % (i)] = evrad

        if code != "grid":
            # Free abundances (original, transformed [X/H] [X/Fe] and errors)
            free_abundances = synth_model.transformed_free_abundances()
        else:
            free_abundances = []

        ### Scale errors using the reduced weigthed chisq (tanh)
        #   * It requires that spectrum errors are well estimated (weights are derived from them)
        # - Better fits and better SNR produce:
        #      - rwchisq(tanh) closer to 1
        #      - smaller scale factor (below 1, so the original error is divided by 1/x)
        #      - smaller errors
        if use_errors:
            #error_scale_factor = synth_model.reduced_wchisq
            #error_scale_factor = synth_model.wchisq
            error_scale_factor = 1.
            # https://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html
            # https://www.gnu.org/software/gsl/manual/gsl-ref_38.html
            #error_scale_factor = np.max((1., synth_model.wchisq/np.sqrt(synth_model.m.dof)))
        else:
            error_scale_factor = 1.
        errors['teff'] *= error_scale_factor
        errors['logg'] *= error_scale_factor
        errors['MH'] *= error_scale_factor
        errors['vmic'] *= error_scale_factor
        errors['vmac'] *= error_scale_factor
        errors['vsini'] *= error_scale_factor
        errors['limb_darkening_coeff'] *= error_scale_factor
        errors['R'] *= error_scale_factor
        for i, evrad in enumerate(synth_model.evrad()):
            errors['vrad%04i' % (i)] *= error_scale_factor
        for i in range(len(free_abundances)):
            free_abundances['eAbund'][i] *= error_scale_factor
            free_abundances['e[X/H]'][i] *= error_scale_factor
            free_abundances['e[X/Fe]'][i] *= error_scale_factor
            free_abundances['eA(X)'][i] *= error_scale_factor

        if code != "grid":
            free_loggf = synth_model.transformed_free_loggf()
        else:
            free_loggf = []

        all_params.append(params)
        all_errors.append(errors)
        all_free_abundances.append(free_abundances)
        all_free_loggf.append(free_loggf)

    if n_stellar_components == 1:
        # For backward compatibility, if the fit was for a single star, just return results not in an array
        all_params = all_params[0]
        all_errors = all_errors[0]
        all_free_abundances = all_free_abundances[0]
        all_free_loggf = all_free_loggf[0]

    status = {}
    status['days'] = synth_model.calculation_time.day-1
    status['hours'] = synth_model.calculation_time.hour
    status['minutes'] = synth_model.calculation_time.minute
    status['seconds'] = synth_model.calculation_time.second
    status['dof'] = synth_model.m.dof
    status['error'] = synth_model.m.errmsg
    status['rms'] = synth_model.rms

    # Unweighted
    status['chisq'] = synth_model.chisq
    status['rchisq'] = synth_model.reduced_chisq

    # Weighted
    status['wchisq'] = synth_model.wchisq
    status['rwchisq'] = synth_model.reduced_wchisq

    status['niter'] = synth_model.m.niter
    status['nsynthesis'] = synth_model.m.nfev
    status['status'] = synth_model.m.status

    synth_spectrum = create_spectrum_structure(filtered_spectrum['waveobs'], synth_model.combined_final_fluxes)


    return spectrum, synth_spectrum, all_params, all_errors, all_free_abundances, all_free_loggf, status, stats_linemasks


