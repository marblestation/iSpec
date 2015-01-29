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
import time
from datetime import datetime, timedelta
import numpy as np
from mpfitmodels import MPFitModel
#from continuum import fit_continuum
from abundances import write_solar_abundances, write_fixed_abundances, determine_abundances, create_free_abundances_structure
from abundances import determine_abundance_enchancements, enhance_solar_abundances
from segments import create_segments_around_lines
from atmospheres import write_atmosphere, interpolate_atmosphere_layers, valid_atmosphere_target
from lines import write_atomic_linelist, write_isotope_data
from common import mkdir_p
from spectrum import create_spectrum_structure, convolve_spectrum, correct_velocity, resample_spectrum, read_spectrum, normalize_spectrum, create_wavelength_filter, read_spectrum, write_spectrum
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import JoinableQueue
from Queue import Empty

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
import multiprocessing
from multiprocessing import Pool

import tempfile
import log
import logging


class Constants:
    ###################################
    # CONSTANTS
    ###################################
    SYNTH_STEP_TEFF = 100.
    SYNTH_STEP_LOGG = 0.1
    SYNTH_STEP_MH = 0.05
    SYNTH_STEP_VMIC = 0.5
    SYNTH_STEP_VMAC = 2.0
    SYNTH_STEP_VSINI = 2.0
    SYNTH_STEP_LIMB_DARKENING_COEFF = 0.20
    SYNTH_STEP_R = 100
    ###################################
    EW_STEP_TEFF = 500.
    EW_STEP_LOGG = 0.5
    EW_STEP_MH = 0.05
    EW_STEP_VMIC = 0.5


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, isotope_file=None, regions=None, waveobs_mask=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs.
    In case regions is specified (recarray with 'wave_base' and 'wave_top'),
    only those regions will be computed.

    No macroturbulence, rotation (vsini), limb darkening coefficient or resolution is considered
    in this process. That's why it is named as "fundamental" spectrum.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be set to 'None'.
    """
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel = microturbulence_vel, macroturbulence = 0.0, vsini = 0.0, limb_darkening_coeff = 0.20, R=0, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions)


def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.20, R=500000, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, isotope_file=None, regions=None, waveobs_mask=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs.
    In case regions is specified (recarray with 'wave_base' and 'wave_top'),
    only those regions will be computed.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be set to 'None'.
    """
    if fixed_abundances is None:
        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float), ('element', '|S30')])

    ## All synthetic spectra contain a -0.28 shift (unknown reason) so we correct it
    ## by moving the waveobs in the other sense:
    ## ** The shift has been validated comparing with radial velocities from HARPS
    ##    for a wide range of samples in the parameter space (benchmark stars)
    #rv_shift = 0.28
    #spectrum = create_spectrum_structure(waveobs)
    #spectrum = correct_velocity(spectrum, rv_shift)
    #waveobs = spectrum['waveobs']

    if waveobs_mask is None:
        if regions is None:
            waveobs_mask = np.ones(len(waveobs)) # Compute fluxes for all the wavelengths
            # Limit linelist
            wave_base = np.min(waveobs)
            wave_top = np.max(waveobs)
            lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
            linelist = linelist[lfilter]
        else:
            waveobs_mask = __create_waveobs_mask(waveobs, regions)
            # Limit linelist
            linelist = __filter_linelist(linelist, regions)

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_fixed_abund_file = False
    remove_tmp_linelist_file = False
    remove_tmp_isotope_file = False
    if atmosphere_layers_file is None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
        remove_tmp_atm_file = True
    if abundances_file is None:
        abundances_file = write_solar_abundances(abundances)
        remove_tmp_abund_file = True
    if fixed_abundances_file is None:
        fixed_abundances_file = write_fixed_abundances(fixed_abundances)
        remove_tmp_fixed_abund_file = True
    if linelist_file is None:
        linelist_file = write_atomic_linelist(linelist)
        remove_tmp_linelist_file = True
    if isotope_file is None:
        isotope_file = write_isotope_data(isotopes)
        remove_tmp_isotope_file = True
    nlayers = len(atmosphere_layers)

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    #process_communication_queue = Queue()
    process_communication_queue = JoinableQueue()

    p = Process(target=__generate_spectrum, args=(process_communication_queue, waveobs, waveobs_mask, atmosphere_layers_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file), kwargs={'microturbulence_vel': microturbulence_vel, 'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'nlayers': nlayers, 'verbose': verbose})
    p.start()
    fluxes = np.zeros(len(waveobs))
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            data = process_communication_queue.get(timeout=1)
            if type(data) == np.ndarray:
                # Results received!
                fluxes = data
                break
            elif gui_queue is not None:
                # GUI update
                # It allows communications between process in order to update the GUI progress bar
                gui_queue.put(data)
                gui_queue.join()
            process_communication_queue.task_done()
        except Empty:
            # No results, continue waiting
            pass
        num_seconds += 1
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the synthetic spectrum generation.")
        p.terminate()
    elif np.all(fluxes == 0):
        logging.error("The synthetic spectrum generation has failed for these astrophysical parameters.")
        p.terminate()
    else:
        p.join()

    if remove_tmp_atm_file:
        os.remove(atmosphere_layers_file)
    if remove_tmp_abund_file:
        os.remove(abundances_file)
    if remove_tmp_fixed_abund_file:
        os.remove(fixed_abundances_file)
    if remove_tmp_linelist_file:
        os.remove(linelist_file)
    if remove_tmp_isotope_file:
        os.remove(isotope_file)

    return fluxes


def calculate_theoretical_ew_and_depth(atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, microturbulence_vel = 2.0, atmosphere_layers_file=None, abundances_file=None, linelist_file=None, isotope_file=None, verbose=0, gui_queue=None, timeout=1800):
    """
    """

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_linelist_file = False
    remove_tmp_isotope_file = False
    if atmosphere_layers_file is None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
        remove_tmp_atm_file = True
    if abundances_file is None:
        abundances_file = write_solar_abundances(abundances)
        remove_tmp_abund_file = True
    if linelist_file is None:
        linelist_file = write_atomic_linelist(linelist)
        remove_tmp_linelist_file = True
    if isotope_file is None:
        isotope_file = write_isotope_data(isotopes)
        remove_tmp_isotope_file = True
    nlayers = len(atmosphere_layers)
    start = np.min(linelist['wave (A)']) - 0.1
    end = np.max(linelist['wave (A)']) + 0.1
    num_lines = len(linelist)

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    #process_communication_queue = Queue()
    process_communication_queue = JoinableQueue()

    p = Process(target=__calculate_ew_and_depth, args=(process_communication_queue, atmosphere_layers_file, linelist_file, isotope_file, abundances_file, num_lines), kwargs={'microturbulence_vel': microturbulence_vel, 'nlayers': nlayers, 'start': start, 'end':end, 'verbose': verbose})
    p.start()
    output_wave = np.zeros(len(linelist))
    output_code = np.zeros(len(linelist))
    output_ew = np.zeros(len(linelist))
    output_depth = np.zeros(len(linelist))
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            data = process_communication_queue.get(timeout=1)
            if type(data) is tuple:
                # Results received!
                output_wave, output_code, output_ew, output_depth = data
                break
            elif gui_queue is not None:
                # GUI update
                # It allows communications between process in order to update the GUI progress bar
                gui_queue.put(data)
                gui_queue.join()
            process_communication_queue.task_done()
        except Empty:
            # No results, continue waiting
            pass
        num_seconds += 1
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the synthetic spectrum generation.")
        p.terminate()
    else:
        p.join()

    if remove_tmp_atm_file:
        os.remove(atmosphere_layers_file)
    if remove_tmp_abund_file:
        os.remove(abundances_file)
    if remove_tmp_linelist_file:
        os.remove(linelist_file)
    if remove_tmp_isotope_file:
        os.remove(isotope_file)

    resulting_linelist = linelist.copy()
    resulting_linelist["valid_theoretical_ew_depth"] = np.abs(output_wave - linelist['wave (A)']) < 1e-5
    resulting_linelist["theoretical_ew"] = np.round(output_ew, 2)
    resulting_linelist["theoretical_depth"] = np.round(output_depth, 2)
    return resulting_linelist

def __enqueue_progress(process_communication_queue, v):
    process_communication_queue.put(("self.update_progress(%i)" % v))
    process_communication_queue.join()

def __generate_spectrum(process_communication_queue, waveobs, waveobs_mask, atmosphere_model_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.20, R=500000, nlayers=56, verbose=0):
    """
    Generate synthetic spectrum and apply macroturbulence, rotation (visini), limb darkening coeff and resolution except
    if all those parameters are set to zero, in that case the fundamental synthetic spectrum is returned.
    """
    import synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    ## The convolution (R), rotation broadening (vsini) and macroturbulence broadening (vmac),
    ## do not seem to work as expected in the SPECTRUM code, thus we set them to zero and
    ## we use a python implementation
    #fluxes = synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, 0, nlayers, verbose, update_progress_func)
    fluxes = synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file, microturbulence_vel, 0, 0, limb_darkening_coeff, 0, nlayers, verbose, update_progress_func)
    # Avoid zero fluxes, set a minimum value so that when it is convolved it
    # changes. This way we reduce the impact of the following problem:
    # SPECTRUM + MARCS makes some strong lines to have zero fluxes (i.e. 854.21nm)
    zeros = np.where(fluxes <= 1.0e-10)[0]
    fluxes[zeros] = 1.0e-10
    if R > 0:
        # Use iSpec convolution routine instead of SPECTRUM one, since iSpec is more reliable
        spectrum = create_spectrum_structure(waveobs, fluxes)
        fluxes = convolve_spectrum(spectrum, R, from_resolution=None, frame=None)['flux']
    if macroturbulence > 0:
        fluxes = __broadening_macroturbulent(waveobs, fluxes, macroturbulence, return_kernel=False)
    if vsini > 0:
        fluxes = __broadening_rotational(waveobs, fluxes, vsini)

    process_communication_queue.put(fluxes)

def __calculate_ew_and_depth(process_communication_queue, atmosphere_model_file, linelist_file, isotope_file, abundances_file, num_lines, microturbulence_vel = 2.0, nlayers=56, start=3000, end=11000, verbose=0):
    """
    start and end in Amstrom
    """
    import synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    output_wave, output_code, output_ew, output_depth = synthesizer.calculate_ew_and_depth(atmosphere_model_file, linelist_file, isotope_file, abundances_file, num_lines, microturbulence_vel, nlayers, start, end, verbose, update_progress_func)
    process_communication_queue.put((output_wave, output_code, output_ew, output_depth))


def apply_post_fundamental_effects(waveobs, fluxes, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.20, R=500000, verbose=0, gui_queue=None, timeout=1800):
    """
    Apply macroturbulence, rotation (vsini), limb darkening coefficient and/or resolution
    """

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    process_communication_queue = Queue()

    p = Process(target=__apply_post_fundamental_effects, args=(process_communication_queue, waveobs, fluxes), kwargs={'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'verbose': verbose})
    p.start()
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            data = process_communication_queue.get(timeout=1)
            if type(data) == np.ndarray:
                # Results received!
                fluxes = data
                break
            elif gui_queue is not None:
                # GUI update
                # It allows communications between process in order to update the GUI progress bar
                gui_queue.put(data)
        except Empty:
            # No results, continue waiting
            pass
        num_seconds += 1
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the application of post fundamental effects.")
        p.terminate()
    elif np.all(fluxes == 0):
        logging.error("The application of post fundamental effects has failed.")
        p.terminate()
    else:
        p.join()

    return fluxes


def __apply_post_fundamental_effects(process_communication_queue, waveobs, fluxes, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.20, R=500000, verbose=0):
    """
    Apply macroturbulence, rotation (visini), limb darkening coeff and resolution to already generated fundamental synthetic spectrum.
    """
    import synthesizer
    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    ## The convolution (R), rotation broadening (vsini) and macroturbulence broadening (vmac),
    ## do not seem to work as expected in the SPECTRUM code, thus we set them to zero and
    ## we use a python implementation
    #fluxes = synthesizer.apply_post_fundamental_effects(waveobs*10., fluxes, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, 0, verbose, update_progress_func)
    fluxes = synthesizer.apply_post_fundamental_effects(waveobs*10., fluxes, microturbulence_vel, 0., 0., limb_darkening_coeff, 0, verbose, update_progress_func)
    # Avoid zero fluxes, set a minimum value so that when it is convolved it
    # changes. This way we reduce the impact of the following problem:
    # SPECTRUM + MARCS makes some strong lines to have zero fluxes (i.e. 854.21nm)
    zeros = np.where(fluxes <= 1.0e-10)[0]
    fluxes[zeros] = 1.0e-10
    if R > 0:
        # Use iSpec convolution routine instead of SPECTRUM one, since iSpec is more reliable
        spectrum = create_spectrum_structure(waveobs, fluxes)
        fluxes = convolve_spectrum(spectrum, R, from_resolution=None, frame=None)['flux']
    if macroturbulence > 0:
        fluxes = __broadening_macroturbulent(waveobs, fluxes, macroturbulence, return_kernel=False)
    if vsini > 0:
        fluxes = __broadening_rotational(waveobs, fluxes, vsini)
    process_communication_queue.put(fluxes)


class SynthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, isotopes, abundances, enhance_abundances=True, scale=None, teff=5000, logg=3.0, MH=0.0, vmic=2.0, vmac=0.0, vsini=2.0, limb_darkening_coeff=0.0, R=0, precomputed_grid_dir=None):
        self.precomputed_grid_dir = precomputed_grid_dir
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
        self.abundances = abundances
        self.enhance_abundances = enhance_abundances
        self.scale = scale
        #
        self.calculation_time = 0
        self.waveobs = None
        self.waveobs_mask = None
        self.cache = {}
        p = [teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, R ]
        #
        self.abundances_file = None
        self.linelist_file = None
        self.isotope_file = None
        super(SynthModel, self).__init__(p)

    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]


        # Consider new abundances as fixed
        fixed_abundances = self.free_abundances()

        abundances_key = " ".join(map(str, fixed_abundances['Abund']))
        complete_key = "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %i " % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()))
        complete_key += abundances_key

        key = "%.2f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic())
        key += abundances_key

        ##### [start] Check precomputed (solar abundance)
        precomputed_file = str(self.precomputed_grid_dir) + "/unconvolved_steps/{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}.fits".format(int(self.teff()), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff())
        if self.precomputed_grid_dir is not None and abundances_key == "" and os.path.exists(precomputed_file):
            if not self.quiet:
                print "Pre-computed:", complete_key
            precomputed = read_spectrum(precomputed_file)
            convolved_precomputed = convolve_spectrum(precomputed, self.R())

            convolved_precomputed = resample_spectrum(convolved_precomputed, self.waveobs, method="bessel", zero_edges=True)
            convolved_precomputed['flux'][self.waveobs_mask == 0] = 1.
            self.last_final_fluxes = convolved_precomputed['flux'].copy()
        else:
            if self.cache.has_key(key):
                if not self.quiet:
                    print "Cache:", complete_key
                self.last_fluxes = self.cache[key].copy()
            else:
                if not self.quiet:
                    print "Generating:", complete_key


                # Enhance alpha elements + CNO abundances following MARCS standard composition
                if self.enhance_abundances:
                    alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = determine_abundance_enchancements(self.MH(), scale=self.scale)
                    abundances = enhance_solar_abundances(self.abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)
                else:
                    abundances = self.abundances

                # Atmosphere
                atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH())
                # Fundamental synthetic fluxes
                self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), self.linelist, self.isotopes, abundances, fixed_abundances, microturbulence_vel=self.vmic(), abundances_file=self.abundances_file, linelist_file=self.linelist_file, isotope_file=self.isotope_file, waveobs_mask=self.waveobs_mask, verbose=0)

                # Optimization to avoid too small changes in parameters or repetition
                self.cache[key] = self.last_fluxes.copy()

            self.last_final_fluxes = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, macroturbulence=self.vmac(), vsini=self.vsini(), limb_darkening_coeff=self.limb_darkening_coeff(), R=self.R(), verbose=0)

        return self.last_final_fluxes[self.comparing_mask]

    def fitData(self, waveobs, waveobs_mask, comparing_mask, fluxes, weights=None, parinfo=None, use_errors=False, max_iterations=20, quiet=True):
        self.use_errors = use_errors
        base = 8
        if len(parinfo) < base:
            raise Exception("Wrong number of parameters!")

        if sys.platform == "win32":
            # On Windows, the best timer is time.clock()
            default_timer = time.clock
        else:
            # On most other platforms the best timer is time.time()
            default_timer = time.time
        self.waveobs = waveobs
        self.waveobs_mask = waveobs_mask # Synthesis for wavelengths with mask different from 0.0
        self.comparing_mask = comparing_mask == 1.0 # Wavelengths to be compared for the least square algorithm
        if weights is None:
            weights = np.ones(len(waveobs))
        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 0.0   # Do not limit residuals between -1.0 and 1.0 (np.tanh(residuals/1.0))
        _t0 = default_timer()

        # Write abundances and linelist to avoid writing the same info in each iteration
        if self.enhance_abundances:
            self.abundances_file = None
        else:
            self.abundances_file = write_solar_abundances(self.abundances)
        self.linelist_file = write_atomic_linelist(self.linelist)
        self.isotope_file = write_isotope_data(self.isotopes)

        if self.use_errors:
            super(SynthModel, self).fitData(waveobs[self.comparing_mask], fluxes[self.comparing_mask], weights=weights[self.comparing_mask], parinfo=parinfo, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet)
        else:
            # Do not consider errors for minimization (all weights set to one)
            ones = np.ones(len(fluxes))
            super(SynthModel, self).fitData(waveobs[self.comparing_mask], fluxes[self.comparing_mask], weights=ones[self.comparing_mask], parinfo=parinfo, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet)

        residuals = self.last_final_fluxes[self.comparing_mask] - fluxes[self.comparing_mask]
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))

        #### Unweighted (no errors considered):
        self.chisq = np.sum((residuals)**2)
        self.reduced_chisq = self.chisq / self.m.dof

        #### Weighted (errors considered):
        self.wchisq = np.sum((weights[self.comparing_mask] * residuals)**2)
        self.reduced_wchisq = self.wchisq / self.m.dof

        self.cache = {}

        if self.abundances_file is not None:
            os.remove(self.abundances_file)
        os.remove(self.linelist_file)
        os.remove(self.isotope_file)
        self.abundances_file = None
        self.linelist_file = None
        self.isotope_file = None

        _t1 = default_timer()
        sec = timedelta(seconds=int(_t1 - _t0))
        self.calculation_time = datetime(1,1,1) + sec

    def teff(self): return self._parinfo[0]['value']
    def logg(self): return self._parinfo[1]['value']
    def MH(self): return self._parinfo[2]['value']
    def vmic(self): return self._parinfo[3]['value']
    def vmac(self): return self._parinfo[4]['value']
    def vsini(self): return self._parinfo[5]['value']
    def limb_darkening_coeff(self): return self._parinfo[6]['value']
    def R(self): return self._parinfo[7]['value']
    def free_abundances(self):
        base = 8
        fixed_abundances = np.recarray((len(self._parinfo)-base, ), dtype=[('code', int),('Abund', float), ('element', '|S30')])
        for i in xrange(len(self._parinfo)-base):
            fixed_abundances['code'] = int(self._parinfo[base+i]['parname'])
            fixed_abundances['Abund'] = self._parinfo[base+i]['value']
            fixed_abundances['element'] = ""
        return fixed_abundances

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def eMH(self): return self.m.perror[2]
    def evmic(self): return self.m.perror[3]
    def evmac(self): return self.m.perror[4]
    def evsini(self): return self.m.perror[5]
    def elimb_darkening_coeff(self): return self.m.perror[6]
    def eR(self): return self.m.perror[7]
    def efree_abundances(self):
        base = 8
        eabundances = []
        for i in xrange(len(self._parinfo)-base):
            eabundances.append(self.m.perror[base+i])
        return eabundances

    def transformed_free_abundances(self):
        free_abundances = self.free_abundances()
        efree_abundances = self.efree_abundances()

        transformed_abund = np.recarray((len(free_abundances), ), dtype=[('code', int),('Abund', float), ('element', '|S5'), ('[X/H]', float), ('A(X)', float), ('[X/Fe]', float), ('eAbund', float), ('e[X/H]', float), ('e[X/Fe]', float), ('eA(X)', float)])

        #function abundances, sme, solar_relative
        #abund = sme.abund
        #solar_abund, solar, labels

        #solar_relative = abund-solar

        #abund[2:*] += sme.feh ; rescale by metallicity
        #abund[1:*] = 10^abund[1:*] ; transform to linear
        #return, alog10(abund / abund[0]) + 12
        #end
        sun_log_Nh_over_Ntotal = self.abundances['Abund'][self.abundances['code'] == 1]
        for i in xrange(len(free_abundances)):
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

    def print_solution(self):
        if self.use_errors:
            #error_scale_factor = self.reduced_wchisq
            #error_scale_factor = self.wchisq
            error_scale_factor = 1.
            # https://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html
            # https://www.gnu.org/software/gsl/manual/gsl-ref_38.html
            #error_scale_factor = np.max((1., self.wchisq/np.sqrt(self.m.dof)))
        else:
            error_scale_factor = 1.
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","vmic","vmac","vsini","limb","R")
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i" % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()))
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i" % (self.eteff()*error_scale_factor, self.elogg()*error_scale_factor, self.eMH()*error_scale_factor, self.evmic()*error_scale_factor, self.evmac()*error_scale_factor, self.evsini()*error_scale_factor, self.elimb_darkening_coeff()*error_scale_factor, int(self.eR()*error_scale_factor))

        # Append free individual abundances
        abundances_header = ""
        abundances_solution = ""
        abundances_errors = ""

        transformed_abund = self.transformed_free_abundances()
        for i in xrange(len(transformed_abund)):
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

        print "           ", header
        print "Solution:  ", solution
        print "Errors:    ", errors
        print ""
        if len(transformed_abund) > 0:
            print "           ", abundances_header
            print "Abundances:", abundances_solution
            print "Ab. errors:", abundances_errors
            print ""

        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","wchisq","rwchisq","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (self.m.dof, self.m.niter, self.m.nfev, self.wchisq, self.reduced_wchisq, self.chisq, self.reduced_chisq, self.rms)
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status




def __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, free_abundances, teff_range, logg_range, MH_range):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 8
    free_params = [param.lower() for param in free_params]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base+len(free_abundances))]
    ##
    # Establish limits one step further away from the real limit
    min_teff = np.min(teff_range)
    max_teff = np.max(teff_range)
    if initial_teff == min_teff:
        initial_teff += Constants.SYNTH_STEP_TEFF
    elif initial_teff == min_teff:
        initial_teff -= Constants.SYNTH_STEP_TEFF
    min_teff += Constants.SYNTH_STEP_TEFF
    max_teff -= Constants.SYNTH_STEP_TEFF
    #
    parinfo[0]['parname'] = "teff"
    parinfo[0]['value'] = initial_teff
    parinfo[0]['fixed'] = not parinfo[0]['parname'].lower() in free_params
    parinfo[0]['step'] = Constants.SYNTH_STEP_TEFF # For auto-derivatives
    parinfo[0]['limited'] = [True, True]
    parinfo[0]['limits'] = [min_teff, max_teff]
    ##
    # Establish limits one step further away from the real limit
    min_logg = np.min(logg_range)
    max_logg = np.max(logg_range)
    if initial_logg == min_logg:
        initial_logg += Constants.SYNTH_STEP_LOGG
    elif initial_logg == min_logg:
        initial_logg -= Constants.SYNTH_STEP_LOGG
    min_logg += Constants.SYNTH_STEP_LOGG
    max_logg -= Constants.SYNTH_STEP_LOGG
    #
    parinfo[1]['parname'] = "logg"
    parinfo[1]['value'] = initial_logg
    parinfo[1]['fixed'] = not parinfo[1]['parname'].lower() in free_params
    parinfo[1]['step'] = Constants.SYNTH_STEP_LOGG # For auto-derivatives
    #parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [min_logg, max_logg]
    ##
    # Establish limits one step further away from the real limit
    min_MH = np.min(MH_range)
    max_MH = np.max(MH_range)
    if initial_MH == min_MH:
        initial_MH += Constants.SYNTH_STEP_MH
    elif initial_MH == min_MH:
        initial_MH -= Constants.SYNTH_STEP_MH
    min_MH += Constants.SYNTH_STEP_MH
    max_MH -= Constants.SYNTH_STEP_MH
    #
    parinfo[2]['parname'] = "MH"
    parinfo[2]['value'] = initial_MH
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = Constants.SYNTH_STEP_MH # For auto-derivatives
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [min_MH, max_MH]
    #
    parinfo[3]['parname'] = "Vmic"
    parinfo[3]['value'] = initial_vmic
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = Constants.SYNTH_STEP_VMIC # For auto-derivatives
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [0.0, 50.0]
    #
    parinfo[4]['parname'] = "Vmac"
    parinfo[4]['value'] = initial_vmac
    parinfo[4]['fixed'] = not parinfo[4]['parname'].lower() in free_params
    parinfo[4]['step'] = Constants.SYNTH_STEP_VMAC # For auto-derivatives
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [0.0, 50.0]
    #
    parinfo[5]['parname'] = "Vsini"
    parinfo[5]['value'] = initial_vsini
    parinfo[5]['fixed'] = not parinfo[5]['parname'].lower() in free_params
    parinfo[5]['step'] = Constants.SYNTH_STEP_VSINI # For auto-derivatives
    parinfo[5]['limited'] = [True, True]
    parinfo[5]['limits'] = [0.0, 50.0]
    #
    parinfo[6]['parname'] = "limb_darkening_coeff"
    parinfo[6]['value'] = initial_limb_darkening_coeff
    parinfo[6]['fixed'] = not parinfo[6]['parname'].lower() in free_params
    parinfo[6]['step'] = Constants.SYNTH_STEP_LIMB_DARKENING_COEFF # For auto-derivatives
    parinfo[6]['limited'] = [True, True]
    parinfo[6]['limits'] = [0.0, 1.0]
    #
    parinfo[7]['parname'] = "R"
    parinfo[7]['value'] = initial_R
    parinfo[7]['fixed'] = not parinfo[7]['parname'].lower() in free_params
    parinfo[7]['step'] = Constants.SYNTH_STEP_R # For auto-derivatives
    parinfo[7]['limited'] = [True, True]
    parinfo[7]['limits'] = [500.0, 300000.0]
    #
    base = 8
    for i in xrange(len(free_abundances)):
        parinfo[base+i]['parname'] = str(free_abundances['code'][i])
        parinfo[base+i]['value'] = free_abundances['Abund'][i]
        parinfo[base+i]['fixed'] = not parinfo[base+i]['parname'].lower() in free_params
        parinfo[base+i]['step'] = 0.05 # For auto-derivatives
        parinfo[base+i]['limited'] = [True, True]
        parinfo[base+i]['limits'] = [-30., 0.]

    return parinfo

def __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, teff_range, logg_range, MH_range, free_params, adjust_model_metalicity=False):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 4
    free_params = [param.lower() for param in free_params]
    #parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base)]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(base)]
    #
    parinfo[0]['parname'] = "teff"
    parinfo[0]['value'] = initial_teff
    parinfo[0]['fixed'] = not parinfo[0]['parname'].lower() in free_params
    parinfo[0]['step'] = Constants.EW_STEP_TEFF # For auto-derivatives
    #parinfo[0]['mpside'] = 2
    #parinfo[0]['mpmaxstep'] = parinfo[0]['step'] * 1.5
    parinfo[0]['limited'] = [True, True]
    parinfo[0]['limits'] = [np.min(teff_range), np.max(teff_range)]
    #
    parinfo[1]['parname'] = "logg"
    parinfo[1]['value'] = initial_logg
    parinfo[1]['fixed'] = not parinfo[1]['parname'].lower() in free_params
    parinfo[1]['step'] = Constants.EW_STEP_LOGG # For auto-derivatives
    #parinfo[1]['mpside'] = 2
    #parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    #parinfo[1]['mpmaxstep'] = parinfo[1]['step'] * 1.5
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [np.min(logg_range), np.max(logg_range)]
    #
    parinfo[2]['parname'] = "Vmic"
    parinfo[2]['value'] = initial_vmic
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = Constants.EW_STEP_VMIC # For auto-derivatives
    #parinfo[2]['mpside'] = 2
    #parinfo[2]['mpmaxstep'] = parinfo[2]['step'] * 2.0
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [0., 50.0]
    #
    parinfo[3]['parname'] = "MH"
    parinfo[3]['value'] = initial_MH
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = Constants.EW_STEP_MH # For auto-derivatives
    #parinfo[3]['mpside'] = 2
    #if not parinfo[3]['fixed']:
        #parinfo[3]['mpmaxstep'] = parinfo[3]['step'] * 1.5
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [np.min(MH_range), np.max(MH_range)]

    return parinfo


def __filter_linelist(linelist, segments):
    # Build wavelength points from regions
    lfilter = None
    for region in segments:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if lfilter is None:
            lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
        else:
            lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.))

    if lfilter is not None:
        return linelist[lfilter]
    else:
        return linelist

def __create_waveobs_mask(waveobs, segments):
    # Build wavelength points from regions
    wfilter = None
    for region in segments:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if wfilter is None:
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        else:
            wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
    waveobs_mask = np.zeros(len(waveobs))
    waveobs_mask[wfilter] = 1.0 # Compute fluxes only for selected segments

    return waveobs_mask

def __filter_linemasks_not_in_segments(linemasks, segments):
    if segments is None:
        return linemasks
    else:
        lfilter = linemasks['wave_base'] == -1
        for i, region in enumerate(linemasks):
            wave_base = region['wave_base']
            wave_top = region['wave_top']

            # Consider only lines that are inside segments
            in_segment1 = np.logical_and(segments['wave_base'] <= wave_base, segments['wave_top'] >= wave_base)
            in_segment2 = np.logical_and(segments['wave_base'] <= wave_top, segments['wave_top'] >= wave_top)
            in_segment = np.logical_and(in_segment1, in_segment2)
            if np.all(in_segment == False):
                lfilter[i] = False
            else:
                lfilter[i] = True

        return linemasks[lfilter]

def __create_comparing_mask(waveobs, linemasks, segments):
    # Build wavelength points from regions
    wfilter = None
    for region in linemasks:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        # Consider only lines that are inside segments
        if segments is not None:
            in_segment1 = np.logical_and(segments['wave_base'] <= wave_base, segments['wave_top'] >= wave_base)
            in_segment2 = np.logical_and(segments['wave_base'] <= wave_top, segments['wave_top'] >= wave_top)
            in_segment = np.logical_and(in_segment1, in_segment2)
            if np.all(in_segment == False):
                continue

        if wfilter is None:
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        else:
            wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
    waveobs_linemask = np.zeros(len(waveobs))
    waveobs_linemask[wfilter] = 1.0 # Consider fluxes only for selected line masks

    return waveobs_linemask

def __get_stats_per_linemask(waveobs, fluxes, synthetic_fluxes, weights, free_params, linemasks, verbose=False):

    results = np.recarray((len(linemasks), ), dtype=[('wave_peak', float),('wave_base', float),('wave_top', float),('chisq', float),('rchisq', float),('wchisq', float),('rwchisq', float),('rms', float)])
    results['wave_peak'] = linemasks['wave_peak']
    results['wave_base'] = linemasks['wave_base']
    results['wave_top'] = linemasks['wave_top']

    i = 0
    for region in linemasks:
        wave_peak = region['wave_peak']
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        # Do not compare negative or zero fluxes
        wfilter = np.logical_and(wfilter, fluxes > 0.0)

        # Degrees of freedom
        dof = len(waveobs[wfilter]) - len(free_params)
        if dof > 0:
            residuals = synthetic_fluxes[wfilter] - fluxes[wfilter]
            rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))

            # Unweighted
            chisq = np.sum((residuals)**2)
            reduced_chisq = chisq / dof

            # Weighted
            wchisq = np.sum((weights[wfilter] * residuals)**2)
            reduced_wchisq = wchisq / dof
        else:
            rms = -9999
            chisq = -9999
            reduced_chisq = -9999
            wchisq = -9999
            reduced_wchisq = -9999

        results['rms'][i] = rms
        # Unweighted
        results['chisq'][i] = chisq
        results['rchisq'][i] = reduced_chisq
        # Weighted
        results['wchisq'][i] = wchisq
        results['rwchisq'][i] = reduced_wchisq
        if verbose:
            header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("wave_peak","wave_base","wave_top","wchisq","rwchisq","chisq","rchisq","rms")
            stats = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (wave_peak, wave_base, wave_top, wchisq, reduced_wchisq, chisq, reduced_chisq, rms)
            if i == 0:
                print "         ", header
            print "Line     ", stats
        i += 1

    return results

def model_spectrum(spectrum, continuum_model, modeled_layers_pack, linelist, isotopes, abundances, free_abundances, initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, segments=None, linemasks=None, enhance_abundances=True, scale=None, precomputed_grid_dir=None, use_errors=True, max_iterations=20, verbose=1):
    """
    It matches synthetic spectrum to observed spectrum by applying a least
    square algorithm.

    - free_params is an array that can contain any combination of the following
      strings: ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    - free_abundances can be set to 'None'
    - If segments are specified, the synthetic spectrum will be only generated for
      those regions.
    - If linemasks are specified, only those regions will be used for comparison.
    - It does not compare negative or zero fluxes
    - If enhance_abundances is True, alpha elements and CNO abundances will be scaled
      depending on the metallicity. If scale is None, by default the standard
      MARCS composition will be used (recommended).
    """

    if verbose or verbose == 1:
        verbose = True
        quiet = False
    else:
        verbose = False
        quiet = True

    # Duplicate
    if segments is not None:
        # Wavelengths to be computed: segments
        wfilter = create_wavelength_filter(spectrum, regions=segments)
        spectrum = create_spectrum_structure(spectrum['waveobs'][wfilter], spectrum['flux'][wfilter], spectrum['err'][wfilter])
    else:
        # Wavelengths to be computed: all
        spectrum = spectrum.copy()

    # Normalization
    spectrum = normalize_spectrum(spectrum, continuum_model)

    waveobs = spectrum['waveobs']
    flux = spectrum['flux']
    err = spectrum['err']


    if free_abundances is None:
        # No fixed abundances
        free_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float), ('element', '|S30')])
    else:
        # Add free abundances as free params
        for element in free_abundances:
            free_params.append(str(element['code']))

    if segments is None:
        waveobs_mask = np.ones(len(waveobs)) # Compute fluxes for all the wavelengths
        # Limit linelist
        wave_base = np.min(waveobs)
        wave_top = np.max(waveobs)
        lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
        linelist = linelist[lfilter]
    else:
        waveobs_mask = __create_waveobs_mask(waveobs, segments)
        # Limit linelist
        linelist = __filter_linelist(linelist, segments)

    if linemasks is None:
        comparing_mask = np.ones(len(waveobs)) # Compare all fluxes
    else:
        linemasks = __filter_linemasks_not_in_segments(linemasks, segments)
        # Compare fluxes inside line masks that belong to a segment
        comparing_mask = __create_comparing_mask(waveobs, linemasks, segments)


    ## Fluxes
    negative_zero_flux = flux <= 0.0
    bad_fluxes = np.logical_and(comparing_mask == 1, negative_zero_flux)
    num_bad_fluxes = len(np.where(bad_fluxes)[0])
    # Do not compare negative or zero fluxes
    if num_bad_fluxes > 0:
        logging.warn("%i fluxes have been discarded because they are negative or zero" % num_bad_fluxes)
        comparing_mask[negative_zero_flux] = 0.0

    ## Errors
    if use_errors and np.all(err[comparing_mask == 1] <= 0):
        logging.warn("Use of errors has been desactivated because all of them are set to zero.")
        use_errors = False

    negative_zero_err = err <= 0.0
    bad_errors = np.logical_and(comparing_mask == 1, negative_zero_err)
    num_bad_errors = len(np.where(bad_errors)[0])
    ## Do not compare negative or zero errors
    if use_errors and num_bad_errors > 0:
        logging.warn("%i fluxes have been discarded because their ERRORS are negative or zero" % num_bad_errors)
        comparing_mask[negative_zero_err] = 0.0


    if np.all(comparing_mask == 0):
        logging.error("No fluxes left to be compared!")
        raise Exception("No fluxes left to be compared!")

    accept_weights = np.logical_and(comparing_mask == 1., np.logical_not(negative_zero_err))
    weights = np.ones(len(waveobs))
    weights[accept_weights] = 1. / err[accept_weights]
    weights[~accept_weights] = 0.
    #weights /= np.sum(weights) # Normalize
    #weights /= np.max(weights) # Normalize
    #weights *= len(np.where(accept_weights)[0]) # Scale to number of points
    #weights *= 10000
    weights = np.sqrt(weights)  # When squaring the flux errors, we get more reasonable parameter's errors (empirically validated)


    teff_range = modeled_layers_pack[3]
    logg_range = modeled_layers_pack[4]
    MH_range = modeled_layers_pack[5]

    parinfo = __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, free_abundances, teff_range, logg_range, MH_range)

    synth_model = SynthModel(modeled_layers_pack, linelist, isotopes, abundances, enhance_abundances=enhance_abundances, scale=scale, precomputed_grid_dir=precomputed_grid_dir)

    synth_model.fitData(waveobs, waveobs_mask, comparing_mask, flux, weights=weights, parinfo=parinfo, use_errors=use_errors, max_iterations=max_iterations, quiet=quiet)

    if verbose:
        print "\n"

    stats_linemasks = __get_stats_per_linemask(waveobs, flux, synth_model.last_final_fluxes, weights, free_params, linemasks, verbose=verbose)

    if verbose:
        print "\n"
        synth_model.print_solution()

    # Collect information to be returned
    params = {}
    params['teff'] = synth_model.teff()
    params['logg'] = synth_model.logg()
    params['MH'] = synth_model.MH()
    params['vmic'] = synth_model.vmic()
    params['vmac'] = synth_model.vmac()
    params['vsini'] = synth_model.vsini()
    params['limb_darkening_coeff'] = synth_model.limb_darkening_coeff()
    params['R'] = synth_model.R()

    errors = {}
    errors['teff'] = synth_model.eteff()
    errors['logg'] = synth_model.elogg()
    errors['MH'] = synth_model.eMH()
    errors['vmic'] = synth_model.evmic()
    errors['vmac'] = synth_model.evmac()
    errors['vsini'] = synth_model.evsini()
    errors['limb_darkening_coeff'] = synth_model.elimb_darkening_coeff()
    errors['R'] = synth_model.eR()

    # Free abundances (original, transformed [X/H] [X/Fe] and errors)
    free_abundances = synth_model.transformed_free_abundances()

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
    for i in xrange(len(free_abundances)):
        free_abundances['eAbund'][i] *= error_scale_factor
        free_abundances['e[X/H]'][i] *= error_scale_factor
        free_abundances['e[X/Fe]'][i] *= error_scale_factor
        free_abundances['eA(X)'][i] *= error_scale_factor

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

    synth_spectrum = create_spectrum_structure(waveobs, synth_model.last_final_fluxes)

    return spectrum, synth_spectrum, params, errors, free_abundances, status, stats_linemasks



class EquivalentWidthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, abundances, teff=5000, logg=3.0, MH=0.0, vmic=2.0, adjust_model_metalicity=False, enhance_abundances=True, scale=None):
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
        self.abundances = abundances
        self.enhance_abundances = enhance_abundances
        self.scale = scale
        self.adjust_model_metalicity = adjust_model_metalicity
        self.lines_for_teff = None
        self.lines_for_vmic = None
        #
        self.calculation_time = 0
        self.cache = {}
        self.m1 = None
        self.c1 = None
        self.m2 = None
        self.c2 = None
        self.fe1 = None
        self.fe2 = None
        self.fe1_std = None
        self.fe2_std = None
        self.fe1_filter = None
        self.fe2_filter = None
        p = [teff, logg, vmic]
        self._MH = MH
        self._eMH = 0.0

        self.min_MH = np.min(modeled_layers_pack[5])
        self.max_MH = np.max(modeled_layers_pack[5])
        #
        super(EquivalentWidthModel, self).__init__(p)



    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]

        key = "%.2f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic())
        if self.cache.has_key(key):
            hit_cache = True
            if not self.quiet:
                print "Cache:", key
            self.last_final_values = self.cache[key]
            spec_abund, absolute_abund, x_over_h, x_over_fe = self.cache[key]
        else:
            hit_cache = False
            if not self.quiet:
                print "Generating:", key
            # Optimization to avoid too small changes in parameters or repetition
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH())
            if self.fe1_filter is None or self.fe2_filter is None:
                ignore = np.ones(len(self.linemasks)) # Do not ignore any line since it's the first execution and it has not been done any selection
            else:
                ignore = np.zeros(len(self.linemasks))
                ignore[np.where(np.logical_or(self.fe1_filter, self.fe2_filter))[0]] = 1.0 # Do not ignore selected fe1/2 lines

            # Enhance alpha elements + CNO abundances following MARCS standard composition
            if self.enhance_abundances:
                alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = determine_abundance_enchancements(self.MH(), scale=self.scale)
                abundances = enhance_solar_abundances(self.abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)
            else:
                abundances = self.abundances

            spec_abund, absolute_abund, x_over_h, x_over_fe = determine_abundances(atmosphere_layers, \
                    self.teff(), self.logg(), self.MH(), self.linemasks, abundances, microturbulence_vel = self.vmic(), \
                    ignore=ignore, verbose=0)


            if 'EW_absolute_abund_median' in self.linemasks.dtype.names:
                # Instead of the literature solar abundance, use the solar abundance determined by iSpec (differencial analysis)
                differential_x_over_h = absolute_abund - self.linemasks['EW_absolute_abund_median']
                differential_feh = np.median(differential_x_over_h[self.linemasks['element'] == "Fe 1"])
                differential_x_over_fe = differential_x_over_h - differential_feh
                x_over_h = differential_x_over_h
                x_over_fe = differential_x_over_fe

            self.cache[key] = (spec_abund, absolute_abund, x_over_h, x_over_fe)

        # First iteration
        if self.fe1_filter is None or self.fe2_filter is None:
            self.select_good_lines(x_over_h, strict=True)
            #self.select_good_lines(x_over_h, strict=False) # Don't identify and filter outliers

        values_to_evaluate = []
        fitted_lines_params = []
        selected_x_over_h = []

        import statsmodels.api as sm
        ### Temperature
        ## y = mx + c
        x = self.linemasks['lower state (eV)'][self.fe1_filter]
        y = x_over_h[self.fe1_filter]
        x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
        linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
        self.m1 = linear_model.params[0]
        self.c1 = linear_model.params[1]
        self.fe1 = np.median(x_over_h[self.fe1_filter])
        self.fe1_std = np.std(x_over_h[self.fe1_filter])
        ##self.fe1 = np.median(linear_model.fittedvalues)
        ##self.fe1_std = np.std(linear_model.fittedvalues)
        #print "Fe 1", np.median(linear_model.fittedvalues), np.median(x_over_h[self.fe1_filter])
        #print "    ", np.std(linear_model.fittedvalues), np.std(x_over_h[self.fe1_filter])
        #import matplotlib.pyplot as plt
        #plt.scatter(x, y)
        #plt.plot(x, m1*x + c1)
        #plt.show()

        ### Vmic
        ## y = mx + c
        x = self.linemasks['ewr'][self.fe1_filter]
        y = x_over_h[self.fe1_filter]
        x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
        linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
        self.m2 = linear_model.params[0]
        self.c2 = linear_model.params[1]
        #import matplotlib.pyplot as plt
        #plt.scatter(x, y)
        #plt.plot(x, m2*x + c2)
        #plt.show()

        ### Fe2
        ## y = mx + c
        x = self.linemasks['ewr'][self.fe2_filter]
        if len(x) > 1:
            y = x_over_h[self.fe2_filter]
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
            self.fe2 = np.median(x_over_h[self.fe2_filter])
            self.fe2_std = np.std(x_over_h[self.fe2_filter])
            ##self.fe2 = np.median(linear_model.fittedvalues)
            ##self.fe2_std = np.std(linear_model.fittedvalues)
            #print "Fe 2", np.median(linear_model.fittedvalues), np.median(x_over_h[self.fe2_filter])
            #print "    ", np.std(linear_model.fittedvalues), np.std(x_over_h[self.fe2_filter])
            #import matplotlib.pyplot as plt
            #plt.scatter(x, y)
            #plt.plot(x, m2*x + c2)
            #plt.show()
        else:
            self.fe2 = np.median(x_over_h[self.fe2_filter])
            self.fe2_std = np.std(x_over_h[self.fe2_filter])

        ## Gravity
        abundance_diff = self.fe1 - self.fe2
        abundance_diff2 = self.MH() - self.fe1

        # Rounded to 3 and 2 decimals (using string convertion works better than np.round)
        values_to_evaluate.append(float("%.6f" % self.m1))
        values_to_evaluate.append(float("%.6f" % self.m2))
        values_to_evaluate.append(float("%.6f" % abundance_diff))
        if self.adjust_model_metalicity:
            values_to_evaluate.append(float("%.2f" % abundance_diff2))
        residuals = np.asarray(values_to_evaluate) - self.y

        if not hit_cache:
            print " # Element:                   Fe 1 / Fe 2\n",
            print "   Teff/Vmic slopes:            %.6f %.6f" % (self.m1, self.m2)
            print "   Abundances diff:             %.6f" % abundance_diff
            print "   Abundances diff with model:  %.6f" % abundance_diff2
            print "   Abundances stdev:            %.6f %.6f" % (np.std(x_over_h[self.fe1_filter]), np.std(x_over_h[self.fe2_filter]))
            print "   Abundances median:           %.6f %.6f" % (np.median(self.fe1), np.median(self.fe2))
            print " - Chisq:                       %.10g" % np.sum((self.weights*residuals)**2)

        fitted_lines_params.append(self.m1)
        fitted_lines_params.append(self.c1)
        fitted_lines_params.append(self.m2)
        fitted_lines_params.append(self.c2)
        selected_x_over_h.append(self.fe1_filter.copy())
        selected_x_over_h.append(self.fe2_filter.copy())
        self.last_final_values = (np.asarray(values_to_evaluate), x_over_h, selected_x_over_h, fitted_lines_params)

        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        return values_to_evaluate.copy()

    # Default procedure to be called every iteration.  It simply prints
    # the parameter values.
    import scipy
    blas_enorm32, = scipy.linalg.blas.get_blas_funcs(['nrm2'],np.array([0],dtype=np.float32))
    blas_enorm64, = scipy.linalg.blas.get_blas_funcs(['nrm2'],np.array([0],dtype=np.float64))
    def defiter(self, fcn, x, iter, fnorm=None, functkw=None,
                       quiet=0, iterstop=None, parinfo=None,
                       format=None, pformat='%.10g', dof=1):

        if quiet:
            return
        if fnorm is None:
            [status, fvec] = fcn(x, fjac=None, **functkw)
            # If the returned fvec has more than four bits I assume that we have
            # double precision
            # It is important that the machar is determined by the precision of
            # the returned value, not by the precision of the input array
            if np.array([fvec]).dtype.itemsize>4:
                self.blas_enorm = mpfit.blas_enorm64
            else:
                self.blas_enorm = mpfit.blas_enorm32
            fnorm = self.enorm(fvec)**2

        # Determine which parameters to print
        nprint = len(x)
        print "*Iter ", ('%6i' % iter),"   CHI-SQUARE = ",('%.10g' % fnorm)," DOF = ", ('%i' % dof)
        for i in range(nprint):
            if (parinfo is not None) and (parinfo[i].has_key('parname')):
                p = '   ' + parinfo[i]['parname'] + ' = '
            else:
                p = '   P' + str(i) + ' = '
            if (parinfo is not None) and (parinfo[i].has_key('mpprint')):
                iprint = parinfo[i]['mpprint']
            else:
                iprint = 1
            if iprint:
                print p + (pformat % x[i]) + '  '

        ##### Metallicity
        #self._MH = self.fe1
        #self._MH = np.min((self._MH, self.max_MH))
        #self._MH = np.max((self._MH, self.min_MH))
        #self._eMH = np.std(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])


        ##### Lines
        #values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        #self.select_good_lines(x_over_h, strict=True) # Modifies self.lines_for_teff and self.lines_for_vmic

        return 0


    def select_good_lines(self, x_over_h, strict=True):
        """
            Modifies self.fe1_filter and self.fe2_filter
        """
        # Out of range
        bad = np.logical_or(x_over_h > 1.0, x_over_h < -5)
        #### Line selection
        fe1_filter = self.linemasks['element'] == "Fe 1"
        fe2_filter = self.linemasks['element'] == "Fe 2"

        if strict and len(np.where(~bad)[0]) > 1:
            # Outliers
            import statsmodels.api as sm
            x = self.linemasks['lower state (eV)']
            y = x_over_h
            # RLM (Robust least squares)
            # Huber's T norm with the (default) median absolute deviation scaling
            # - http://en.wikipedia.org/wiki/Huber_loss_function
            # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
            linear_model = huber_t.fit()
            reject_filter1 = linear_model.weights < self.outliers_weight_limit
            reject_filter1 = np.logical_or(reject_filter1, bad)
            #import matplotlib.pyplot as plt
            #plt.scatter(self.linemasks['lower state (eV)'], x_over_h)
            #plt.scatter(self.linemasks['lower state (eV)'][reject_filter1], x_over_h[reject_filter1], color="red")
            #m1 = linear_model.params[0]
            #c1 = linear_model.params[1]
            #plt.plot(x, m1*x + c1, color="green")
            #plt.xlabel("Lower excitation energies (eV)")
            #plt.ylabel("[Fe 1/H] (dex)")
            #plt.grid()
            #plt.savefig("excitation_equilibrium.png")
            #plt.savefig("excitation_equilibrium.eps")
            #plt.savefig("excitation_equilibrium.pdf")
            #plt.show()

            # Outliers
            import statsmodels.api as sm
            x = self.linemasks['ewr']
            y = x_over_h
            # RLM (Robust least squares)
            # Huber's T norm with the (default) median absolute deviation scaling
            # - http://en.wikipedia.org/wiki/Huber_loss_function
            # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
            linear_model = huber_t.fit()
            reject_filter2 = linear_model.weights < self.outliers_weight_limit
            reject_filter2 = np.logical_or(reject_filter2, bad)
            #import matplotlib.pyplot as plt
            #plt.scatter(self.linemasks['ewr'], x_over_h)
            #plt.scatter(self.linemasks['ewr'][reject_filter2], x_over_h[reject_filter2], color="red")
            #plt.show()

            reject_filter = np.logical_or(reject_filter1, reject_filter2)

            # Discard bad lines and outliers
            clean_fe1_filter = np.logical_and(~reject_filter, fe1_filter)
            # Ensure that there are at least some lines
            if len(np.where(clean_fe1_filter)[0]) <= 1:
                # Discard only bad lines
                clean_fe1_filter = np.logical_and(~bad, fe1_filter)
                if len(np.where(clean_fe1_filter)[0]) <= 1:
                    clean_fe1_filter = fe1_filter
            if len(np.where(clean_fe1_filter)[0]) <= 1:
                raise Exception("Not enought lines for Fe 1 (%i lines)" % len(np.where(clean_fe1_filter)[0]))
            else:
                self.fe1_filter = clean_fe1_filter

            # Discard bad lines and outliers
            clean_fe2_filter = np.logical_and(~reject_filter, fe2_filter)
            # Ensure that there are at least some lines
            if len(np.where(clean_fe2_filter)[0]) <= 1:
                # Discard only bad lines
                clean_fe2_filter = np.logical_and(~bad, fe2_filter)
                if len(np.where(clean_fe2_filter)[0]) <= 1:
                    clean_fe2_filter = fe2_filter

            ## Discard ONLY bad lines for Fe 2
            #clean_fe2_filter = np.logical_and(~bad, fe2_filter)
            ## Ensure that there are at least some lines
            #if len(np.where(clean_fe2_filter)[0]) >= 1:
                #clean_fe2_filter = fe2_filter

            if len(np.where(clean_fe2_filter)[0]) <= 1:
                raise Exception("Not enought lines for Fe 1 (%i lines)" % len(np.where(clean_fe2_filter)[0]))
            else:
                self.fe2_filter = clean_fe2_filter
        else:
            ##### ACCEPT all
            if len(np.where(~bad & fe1_filter)[0]) > 0:
                self.fe1_filter = np.logical_and(fe1_filter, np.logical_not(bad))
            else:
                self.fe1_filter = fe1_filter
            if len(np.where(~bad & fe2_filter)[0]) > 0:
                self.fe2_filter = np.logical_and(fe2_filter, np.logical_not(bad))
            else:
                self.fe2_filter = fe2_filter
        print " > Selected Fe 1 lines for teff:", len(np.where(self.fe1_filter)[0]), "of", len(np.where(fe1_filter)[0])
        print " > Selected Fe 2 lines for vmic:", len(np.where(self.fe2_filter)[0]), "of", len(np.where(fe2_filter)[0])



    def fitData(self, linemasks, outliers_weight_limit=0.90, parinfo=None, max_iterations=20, quiet=True):
        base = 3
        if len(parinfo) < base:
            raise Exception("Wrong number of parameters!")

        self.outliers_weight_limit = outliers_weight_limit


        if sys.platform == "win32":
            # On Windows, the best timer is time.clock()
            default_timer = time.clock
        else:
            # On most other platforms the best timer is time.time()
            default_timer = time.time
        self.linemasks = linemasks
        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 0.0   # Not active: Residuals are limited between -1.0 and 1.0 (np.tanh(residuals/1.0))
        #chisq_limit = 4.0e-4 # 0.0004 = np.sum(np.asarray([0.01, 0.01, 0.01, 0.01])**2))
        chisq_limit = 3 # = np.sum((np.asarray([0.01, 0.01, 0.01])*100)**2) # weight 100)

        _t0 = default_timer()

        #index = np.asarray([0, 1, 2])
        #index = np.arange(len(linemasks))
        if self.adjust_model_metalicity:
            index = np.arange(4) # 4 values: zero slopes and zero difference between element1 and element2, difference with model
        else:
            index = np.arange(3) # 3 values: zero slopes and zero difference between element1 and element2
        target_values = np.zeros(len(index))
        weights = np.ones(len(index)) * 100
        #weights = np.asarray([3,1,2])
        #weights = np.asarray([100,100,1])
        #weights = np.asarray([1000,1,100])
        #weights = np.asarray([100,1,1000])
        #weights = np.asarray([1,1000,1])
        super(EquivalentWidthModel, self).fitData(index, target_values, weights=weights, parinfo=parinfo, chisq_limit=chisq_limit, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet, iterfunct=self.defiter)

        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        residuals = values_to_evaluate - target_values
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
        # Unweighted
        self.chisq = np.sum((residuals)**2)
        self.reduced_chisq = self.chisq / self.m.dof
        # Weighted
        self.wchisq = np.sum((weights * residuals)**2)
        self.reduced_wchisq = self.wchisq / self.m.dof

        #self.cache = {}

        _t1 = default_timer()
        sec = timedelta(seconds=int(_t1 - _t0))
        self.calculation_time = datetime(1,1,1) + sec

    def teff(self): return self._parinfo[0]['value']
    def logg(self): return self._parinfo[1]['value']
    def vmic(self): return self._parinfo[2]['value']

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def evmic(self): return self.m.perror[2]

    def MH(self): return self._parinfo[3]['value']
    def eMH(self): return self.m.perror[3]

    #def MH(self): return self._MH
    #def eMH(self): return self._eMH

    def print_solution(self):
        # Calculate MH
        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        #MH = np.median(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        #MH = np.std(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        MH = self.fe1
        eMH = self.fe1_std

        header = "%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","vmic")
        #solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.teff(), self.logg(), self.MH(), self.vmic())
        #errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.eteff(), self.elogg(), self.eMH(), self.evmic())
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.teff(), self.logg(), MH, self.vmic())
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.eteff(), self.elogg(), eMH, self.evmic())

        print "           ", header
        print "Solution:  ", solution
        print "Errors:    ", errors
        print ""

        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","wchisq","rwchisq","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f" % (self.m.dof, self.m.niter, self.m.nfev, self.wchisq, self.reduced_wchisq, self.chisq, self.reduced_chisq, self.rms)
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status


def model_spectrum_from_ew(linemasks, modeled_layers_pack, linelist, abundances, initial_teff, initial_logg, initial_MH, initial_vmic, free_params=["teff", "logg", "vmic"], adjust_model_metalicity=False, enhance_abundances=True, scale=None, max_iterations=20, outliers_weight_limit=0.90):
    """
    The parameter 'outliers_weight_limit' limits the outlier detection done in the first iteration,
    if it is set to 0. then no outliers are filtered. The recommended value is 0.90.
    - If enhance_abundances is True, alpha elements and CNO abundances will be scaled
      depending on the metallicity. If scale is None, by default the standard
      MARCS composition will be used (recommended).
    """
    teff_range = modeled_layers_pack[3]
    logg_range = modeled_layers_pack[4]
    MH_range = modeled_layers_pack[5]

    # Do not allow users to set free MH in free_params to avoid confusions
    # because metallicity is always free in this method, what we make by including MH in free_params
    # is turning on the adjustment in the metallicity models
    if "MH" in free_params or "mh" in free_params:
        raise Exception("Metallicity cannot be a free parameter!")

    if adjust_model_metalicity:
        free_params.append("MH")

    parinfo = __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, teff_range, logg_range, MH_range, free_params, adjust_model_metalicity=adjust_model_metalicity)


    EW_model = EquivalentWidthModel(modeled_layers_pack, linelist, abundances, MH=initial_MH, adjust_model_metalicity=adjust_model_metalicity, \
                                        enhance_abundances=enhance_abundances, scale=scale)

    lfilter = linemasks['element'] == "Fe 1"
    lfilter = np.logical_or(lfilter, linemasks['element'] == "Fe 2")
    linemasks = linemasks[lfilter]
    EW_model.fitData(linemasks, parinfo=parinfo, max_iterations=max_iterations, quiet=False, outliers_weight_limit=outliers_weight_limit)
    print "\n"
    EW_model.print_solution()

    # Calculate MH
    values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = EW_model.last_final_values
    MH = EW_model.fe1
    eMH = EW_model.fe1_std

    # Collect information to be returned
    params = {}
    params['teff'] = EW_model.teff()
    params['logg'] = EW_model.logg()
    params['MH'] = MH
    params['vmic'] = EW_model.vmic()

    errors = {}
    errors['teff'] = EW_model.eteff()
    errors['logg'] = EW_model.elogg()
    errors['MH'] = eMH
    errors['vmic'] = EW_model.evmic()

    status = {}
    values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = EW_model.last_final_values
    # Save parameters (only for Fe, if there are more elements they will not be saved)
    status['slope_excitation_potential'] = values_to_evaluate[0]
    status['slope_ewr'] = values_to_evaluate[1]
    status['abundance_diff'] = values_to_evaluate[2]
    status['fe1_lines'] = len(np.where(selected_x_over_h[0])[0])
    status['fe2_lines'] = len(np.where(selected_x_over_h[1])[0])
    status['model_MH'] = EW_model.MH()

    status['days'] = EW_model.calculation_time.day-1
    status['hours'] = EW_model.calculation_time.hour
    status['minutes'] = EW_model.calculation_time.minute
    status['seconds'] = EW_model.calculation_time.second
    status['dof'] = EW_model.m.dof
    status['error'] = EW_model.m.errmsg
    status['rms'] = EW_model.rms
    status['chisq'] = EW_model.chisq
    status['rchisq'] = EW_model.reduced_chisq
    status['niter'] = EW_model.m.niter
    status['nsynthesis'] = EW_model.m.nfev
    status['status'] = EW_model.m.status

    return params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params




def __generate_synthetic_fits(filename_out, wavelengths, segments, teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances):
    multiprocessing.current_process().daemon=False

    if valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        print "[started]", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution
        # Prepare atmosphere model
        atmosphere_layers = interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)
        # Synthesis
        synth_spectrum = create_spectrum_structure(wavelengths)
        synth_spectrum['flux'] = generate_spectrum(synth_spectrum['waveobs'], \
                atmosphere_layers, teff, logg, MH, atomic_linelist, isotopes, solar_abundances, \
                fixed_abundances=None, microturbulence_vel = vmic, \
                macroturbulence=vmac, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
                R=resolution, regions=segments, verbose=0)
        # FITS
        write_spectrum(synth_spectrum, filename_out)
        print "[finished]", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution
    else:
        raise Exception("Not valid: %i %.2f %.2f" % (teff, logg, MH))



def precompute_synthetic_grid(output_dirname, ranges, wavelengths, to_resolution, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, enhance_abundances=True, scale=None, segments=None, number_of_processes=1):
    """
    Pre-compute a synthetic grid with some reference ranges (Teff, log(g) and
    MH combinations) and all the steps that iSpec will perform in the
    astrophysical parameter determination process.

    All the non-convolved spectra will be saved in a subdir and a complete
    grid file with the reference points already convolved will be saved in a
    FITS file for fast comparison.

    The output directory can be used by the routines 'model_spectrum' and
    'estimate_initial_ap'.
    """
    reference_list_filename = output_dirname + "/reference.txt"
    reference_grid_filename = output_dirname + "/reference_grid_%i.fits" % to_resolution
    fits_dir = output_dirname + "unconvolved_steps/"
    mkdir_p(fits_dir)

    # Parallelization pool
    if number_of_processes == 1:
        pool = None
    else:
        pool = Pool(number_of_processes)

    # Create grid binary file
    elapsed = 0 # seconds

    num_ref_spec = len(ranges)
    num_spec = num_ref_spec * 8 # Reference + 7 variations in Teff, logg, MH, vmic, vmac, vsini, limb darkening coeff

    i = 0
    for teff, logg, MH in ranges:
        vmic = estimate_vmic(teff, logg, MH)
        vmac = estimate_vmac(teff, logg, MH)
        vsini = 2.0
        limb_darkening_coeff = 0.20
        resolution = 0
        # For each reference point, calculate also the variations that iSpec will perform in the first iteration
        steps =   ((teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff+Constants.SYNTH_STEP_TEFF, logg, MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff, logg+Constants.SYNTH_STEP_LOGG, MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff, logg, MH+Constants.SYNTH_STEP_MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff, logg, MH, vmic+Constants.SYNTH_STEP_VMIC, vmac, vsini, limb_darkening_coeff),
                    (teff, logg, MH, vmic, vmac+Constants.SYNTH_STEP_VMAC, vsini, limb_darkening_coeff),
                    (teff, logg, MH, vmic, vmac, vsini+Constants.SYNTH_STEP_VSINI, limb_darkening_coeff),
                    (teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff+Constants.SYNTH_STEP_LIMB_DARKENING_COEFF))

        for j, (teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff) in enumerate(steps):
            filename_out = fits_dir + "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}".format(int(teff), logg, MH, vmic, vmac, vsini, limb_darkening_coeff) + ".fits"

            if os.path.exists(filename_out):
                print "Skipping", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, "already computed"
                continue

            # Enhance alpha elements + CNO abundances following MARCS standard composition
            if enhance_abundances:
                alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = determine_abundance_enchancements(MH, scale=scale)
                abundances = enhance_solar_abundances(solar_abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)
            else:
                abundances = solar_abundances


            if pool is None:
                if sys.platform == "win32":
                    # On Windows, the best timer is time.clock()
                    default_timer = time.clock
                else:
                    # On most other platforms the best timer is time.time()
                    default_timer = time.time
                tcheck = default_timer()
                # Validate parameters
                __generate_synthetic_fits(filename_out, wavelengths, segments, teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution, modeled_layers_pack, atomic_linelist, isotopes, abundances)
                elapsed = default_timer() - tcheck

                print "-----------------------------------------------------"
                print "Remaining time:"
                print "\t", (num_spec-i)*elapsed, "seconds"
                print "\t", (num_spec-i)*(elapsed/60), "minutes"
                print "\t", (num_spec-i)*(elapsed/(60*60)), "hours"
                print "\t", (num_spec-i)*(elapsed/(60*60*24)), "days"
                print "-----------------------------------------------------"
            else:
                pool.apply_async(__generate_synthetic_fits, [filename_out, wavelengths, segments, teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution, modeled_layers_pack, atomic_linelist, isotopes, abundances])
            i += 1

    if pool is not None:
        pool.close()
        pool.join()


    reference_grid = None
    reference_list = Table()
    reference_list.add_column(Column(name='filename', dtype='|S50'))
    reference_list.add_column(Column(name='teff', dtype=int))
    reference_list.add_column(Column(name='logg', dtype=float))
    reference_list.add_column(Column(name='MH', dtype=float))
    reference_list.add_column(Column(name='vmic', dtype=float))
    reference_list.add_column(Column(name='vmac', dtype=float))
    reference_list.add_column(Column(name='vsini', dtype=float))
    reference_list.add_column(Column(name='limb_darkening_coeff', dtype=float))
    for teff, logg, MH in ranges:
        vmic = estimate_vmic(teff, logg, MH)
        vmac = estimate_vmac(teff, logg, MH)
        vsini = 2.0
        limb_darkening_coeff = 0.20
        resolution = 0
        reference_filename_out = "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}".format(int(teff), logg, MH, vmic, vmac, vsini, limb_darkening_coeff) + ".fits"
        reference_list.add_row((reference_filename_out, int(teff), logg, MH, vmic, vmac, vsini, limb_darkening_coeff))

        # Spectra in the grid is convolved to the specified resolution for fast comparison
        print "Quick grid:", reference_filename_out
        spectrum = read_spectrum(fits_dir + reference_filename_out)
        convolved_spectrum = convolve_spectrum(spectrum, to_resolution)

        if reference_grid is None:
            reference_grid = spectrum['flux']
        else:
            reference_grid = np.vstack((reference_grid, spectrum['flux']))

    ascii.write(reference_list, reference_list_filename, delimiter='\t')
    # Generate FITS file with grid for fast comparison
    primary_hdu = fits.PrimaryHDU(reference_grid)
    wavelengths_hdu = fits.ImageHDU(wavelengths, name="WAVELENGTHS")
    params_bintable_hdu = fits.BinTableHDU(reference_list._data, name="PARAMS")
    fits_format = fits.HDUList([primary_hdu, wavelengths_hdu, params_bintable_hdu])
    fits_format.writeto(reference_grid_filename, clobber=True)


def estimate_vmic(teff, logg, feh):
    """
    Estimate Microturbulence velocity (Vmic) by using an empirical relation
    considering the effective temperature, surface gravity and metallicity.

    The relation was constructed based on the UVES Gaia ESO Survey iDR1 data,
    results for the benchmark stars (Jofre et al. 2013),
    and globular cluster data from external literature sources.

    Source: http://great.ast.cam.ac.uk/GESwiki/GesWg/GesWg11/Microturbulence
    """
    t0 = 5500
    g0 = 4.0

    if logg >= 3.5:
        if teff >= 5000:
            # main sequence and subgiants (RGB)
            vmic = 1.05 + 2.51e-4*(teff-t0) + 1.5e-7*(teff-t0)**2 - 0.14*(logg-g0) - 0.05e-1*(logg-g0)**2 + 0.05*feh + 0.01*feh**2
        else:
            # main sequence
            vmic = 1.05 + 2.51e-4*(5000-t0) + 1.5e-7*(5000-t0)**2 - 0.14*(logg-g0) - 0.05e-1*(logg-g0)**2 + 0.05*feh + 0.01*feh**2
    else:
        # giants (RGB/AGB)
        vmic = 1.25 + 4.01e-4*(teff-t0) + 3.1e-7*(teff-t0)**2 - 0.14*(logg-g0) - 0.05*(logg-g0)**2 + 0.05*feh + 0.01*feh**2
    vmic = float("%.2f" % vmic)
    return vmic



def estimate_vmac(teff, logg, feh):
    """
    Estimate Microturbulence velocity (Vmic) by using an empirical relation
    considering the effective temperature, surface gravity and metallicity.

    The relation was constructed by Maria Bergemann for the Gaia ESO Survey.
    """
    t0 = 5500
    g0 = 4.0

    if logg >= 3.5:
        if teff >= 5000:
            # main sequence and subgiants (RGB)
            vmac = 3*(1.15 + 7e-4*(teff-t0) + 1.2e-6*(teff-t0)**2 - 0.13*(logg-g0) + 0.13*(logg-g0)**2 - 0.37*feh - 0.07*feh**2)
        else:
            # main sequence
            vmac = 3*(1.15 + 2e-4*(teff-t0) + 3.95e-7*(teff-t0)**2 - 0.13*(logg-g0) + 0.13*(logg-g0)**2)
    else:
        # giants (RGB/AGB)
        vmac = 3*(1.15 + 2.2e-5*(teff-t0) - 0.5e-7*(teff-t0)**2 - 0.1*(logg-g0) + 0.04*(logg-g0)**2 - 0.37*feh - 0.07*feh**2)

    vmac = float("%.2f" % vmac)
    return vmac


def estimate_initial_ap(spectrum, precomputed_dir, resolution, linemasks):
    """
    Estimate the initial atmospheric parameters by using a pre-computed grid
    at a given resolution. The comparison will be based on the linemasks.
    """
    initial_teff = 5000.
    initial_logg = 2.5
    initial_MH = 0.0
    initial_vmic = 1.0
    initial_vmac = 0.0
    initial_vsini = 2.0
    initial_limb_darkening_coeff = 0.20

    reference_grid_filename = precomputed_dir + "/reference_grid_%i.fits" % resolution
    if not os.path.exists(reference_grid_filename):
        logging.warn("Pre-computed grid does not exists for R = %i" % resolution)
    else:
        try:
            grid = fits.open(reference_grid_filename)
            grid_waveobs = np.asarray(grid['WAVELENGTHS'].data, dtype=float)
            resampled_spectrum = resample_spectrum(spectrum, grid_waveobs, method="bessel")

            fsegment = create_wavelength_filter(resampled_spectrum, regions=linemasks)
            fsegment = np.logical_and(fsegment, resampled_spectrum['flux'] > 0.0)
            isegment = np.where(fsegment == True)[0]

            # http://en.wikipedia.org/wiki/Goodness_of_fit#Example
            residuals = grid['PRIMARY'].data[:,isegment] - resampled_spectrum['flux'][isegment]
            chisq = np.sum((residuals)**2, axis=1)
            min_j = np.argmin(chisq)
            filename, initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff = grid['PARAMS'].data[min_j]

        except Exception, e:
            print "Initial parameters could not be estimated"
            print type(e), e.message
            pass
        finally:
            grid.close()

    return initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff


################################################################################
# Code from:
#   http://www.phoebe-project.org/2.0/
################################################################################

from scipy.signal import fftconvolve
from scipy.integrate import quad

#cc     = 299792458.        # speed of light              m/s


def __broadening_rotational(wave, flux, vrot, epsilon=0.6, return_kernel=False):
    r"""
    Apply rotational broadening to a spectrum assuming a linear limb darkening
    law.

    The adopted limb darkening law is the linear one, parameterize by the linear
    limb darkening parameter :envvar:`epsilon`. The default value is
    :math:`\varepsilon = 0.6`.

    The rotational kernel is defined in velocity space :math:`v` and is given by

    .. math::

        y = 1 - \left(\frac{v}{v_\mathrm{rot}}\right)^2 \\

        K_\mathrm{rot} = 2 (1-\varepsilon)\sqrt{y} + \frac{\pi}{2} \left(\frac{\varepsilon y}{\pi v_\mathrm{rot}(1-\varepsilon/3)}\right)

    **Construct a simple Gaussian line profile and convolve with vsini=65 km/s**


    >>> sigma = 0.5
    >>> wave = np.linspace(3995, 4005, 1001)
    >>> flux = 1.0 - 0.5*np.exp( - (wave-4000)**2/(2*sigma**2))

    Convolve it with a rotational velocity of :math:`v_\mathrm{rot}=65 \mathrm{km}\,\mathrm{s}^{-1}`:

    >>> vrot = 65.
    >>> flux_broad = tools.broadening_rotational(wave, flux, vrot)
    >>> flux_broad, (wave_kernel, kernel) = tools.broadening_rotational(wave, flux, vrot, return_kernel=True)

    ::

        plt.figure()
        plt.plot(wave, flux, 'k-')
        plt.plot(wave, flux_broad, 'r-', lw=2)
        plt.xlabel("Wavelength [$\AA$]")
        plt.ylabel("Normalised flux")

        plt.figure()
        plt.plot(wave_kernel, kernel, 'r-', lw=2)
        plt.xlabel("Wavelength [$\AA$]")
        plt.ylabel("Normalised flux")


    .. +----------------------------------------------------------------------+----------------------------------------------------------------------+
    .. | .. image:: ../../images/api/spectra/tools/broaden_rotational01.png   | .. image:: ../../images/api/spectra/tools/broaden_rotational02.png   |
    .. |   :scale: 50 %                                                       |   :scale: 50 %                                                       |
    .. +----------------------------------------------------------------------+----------------------------------------------------------------------+


    :parameter wave: Wavelength of the spectrum
    :type wave: array
    :parameter flux: Flux of the spectrum
    :type flux: array
    :parameter vrot: Rotational broadening
    :type vrot: float
    :parameter epsilon: linear limbdarkening parameter
    :type epsilon: float
    :parameter return_kernel: return kernel
    :type return_kernel: bool
    :return: broadened flux [, (wavelength, kernel)]
    :rtype: array [,(array, array)]
    """

    # if there is no rotational velocity, don't bother
    if vrot == 0:
        return flux

    # convert wavelength array into velocity space, this is easier. We also
    # need to make it equidistant
    velo = np.log(wave)
    delta_velo = np.diff(velo).min()
    range_velo = velo.ptp()
    n_velo = int(range_velo/delta_velo) + 1
    velo_ = np.linspace(velo[0], velo[-1], n_velo)
    flux_ = np.interp(velo_, velo, flux)
    dvelo = velo_[1]-velo_[0]
    vrot = vrot / (299792458.*1e-3)
    n_kernel = int(2*vrot/dvelo) + 1

    # The kernel might be of too low resolution, or the the wavelength range
    # might be too narrow. In both cases, raise an appropriate error
    if n_kernel == 0:
        raise ValueError(("Spectrum resolution too low for "
                          "rotational broadening"))
    elif n_kernel > n_velo:
        raise ValueError(("Spectrum range too narrow for "
                          "rotational broadening"))

    # Construct the domain of the kernel
    velo_k = np.arange(n_kernel)*dvelo
    velo_k -= velo_k[-1]/2.

    # transform the velocity array, construct and normalise the broadening
    # kernel
    y = 1 - (velo_k/vrot)**2
    kernel = (2*(1-epsilon)*np.sqrt(y) + \
                       np.pi*epsilon/2.*y) / (np.pi*vrot*(1-epsilon/3.0))
    kernel /= kernel.sum()

    # Convolve the flux with the kernel
    flux_conv = fftconvolve(1-flux_, kernel, mode='same')

    # And interpolate the results back on to the original wavelength array,
    # taking care of even vs. odd-length kernels
    if n_kernel % 2 == 1:
        offset = 0.0
    else:
        offset = dvelo / 2.0

    flux = np.interp(velo+offset, velo_, 1-flux_conv, left=1, right=1)

    # Return the results
    if return_kernel:
        lambda0 = (wave[-1]+wave[0]) / 2.0
        return flux, (velo_k*lambda0, kernel)
    else:
        return flux



def __vmacro_kernel(dlam, Ar, At, Zr, Zt):
    r"""
    Macroturbulent velocity kernel.

    See :py:func:`broadening_macroturbulent` for more information.
    """
    dlam[dlam == 0] = 1e-8
    if Zr != Zt:
        return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) * quad(lambda u: np.exp(-1/u**2),0,Zr/idlam)[0] + \
                          2*At*idlam/(np.sqrt(np.pi)*Zt**2) * quad(lambda u: np.exp(-1/u**2),0,Zt/idlam)[0])
                             for idlam in dlam])
    else:
        return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) + 2*At*idlam/(np.sqrt(np.pi)*Zt**2))\
                           * quad(lambda u: np.exp(-1/u**2),0,Zr/idlam)[0]\
                             for idlam in dlam])


def __broadening_macroturbulent(wave, flux, vmacro_rad, vmacro_tan=None,
                              return_kernel=False):
    r"""
    Apply macroturbulent broadening.

    The macroturbulent kernel is defined as in [Gray2005]_:

    .. math::

        K_\mathrm{macro}(\Delta\lambda) = \frac{2A_R\Delta\lambda}{\sqrt{\pi}\zeta_R^2}\int_0^{\zeta_R/\Delta\lambda}e^{-1/u^2}du

         & + \frac{2A_T\Delta\lambda}{\sqrt{\pi}\zeta_T^2}\int_0^{\zeta_T/\Delta\lambda}e^{-1/u^2}du

    If :envvar:`vmacro_tan` is :envvar:`None`, then the value will be put equal
    to the radial component :envvar:`vmacro_rad`.

    **Example usage**: Construct a simple Gaussian line profile and convolve with vmacro=65km/s

    Construct a simple Gaussian line profile:

    >>> sigma = 0.5
    >>> wave = np.linspace(3995, 4005, 1001)
    >>> flux = 1.0 - 0.5*np.exp( - (wave-4000)**2/(2*sigma**2))

    Convolve it with a macroturbulent velocity of :math:`v_\mathrm{macro}=65 \mathrm{km}\,\mathrm{s}^{-1}`:

    >>> vmac = 65.
    >>> flux_broad = tools.broadening_macroturbulent(wave, flux, vmac)
    >>> flux_broad, (wave_kernel, kernel) = tools.broadening_macroturbulent(wave, flux, vmac, return_kernel=True)

    ::

        plt.figure()
        plt.plot(wave, flux, 'k-')
        plt.plot(wave, flux_broad, 'r-', lw=2)
        plt.xlabel("Wavelength [$\AA$]")
        plt.ylabel("Normalised flux")

        plt.figure()
        plt.plot(wave_kernel, kernel, 'r-', lw=2)
        plt.xlabel("Wavelength [$\AA$]")
        plt.ylabel("Normalised flux")


    .. +--------------------------------------------------------------------------+-------------------------------------------------------------------------+
    .. | .. image:: ../../images/api/spectra/tools/broaden_macroturbulent01.png   | .. image:: ../../images/api/spectra/tools/broaden_macroturbulent02.png  |
    .. |   :scale: 50 %                                                           |   :scale: 50 %                                                          |
    .. +--------------------------------------------------------------------------+-------------------------------------------------------------------------+

    :parameter wave: Wavelength of the spectrum
    :type wave: array
    :parameter flux: Flux of the spectrum
    :type flux: array
    :parameter vmacro_rad: macroturbulent broadening, radial component
    :type vmacro_rad: float
    :parameter vmacro_tan: macroturbulent broadening, tangential component
    :type vmacro_tan: float
    :parameter return_kernel: return kernel
    :type return_kernel: bool
    :return: broadened flux [, (wavelength, kernel)]
    :rtype: array [,(array, array)]
    """
    if vmacro_tan is None:
        vmacro_tan = vmacro_rad

    if vmacro_rad == vmacro_tan == 0:
        return flux

    # Define central wavelength
    lambda0 = (wave[0] + wave[-1]) / 2.0

    vmac_rad = vmacro_rad/(299792458.*1e-3)*lambda0
    vmac_tan = vmacro_tan/(299792458.*1e-3)*lambda0

    # Make sure the wavelength range is equidistant before applying the
    # convolution
    delta_wave = np.diff(wave).min()
    range_wave = wave.ptp()
    n_wave = int(range_wave/delta_wave)+1
    wave_ = np.linspace(wave[0], wave[-1], n_wave)
    flux_ = np.interp(wave_, wave, flux)
    dwave = wave_[1]-wave_[0]
    n_kernel = int(5*max(vmac_rad, vmac_tan)/dwave)
    if n_kernel % 2 == 0:
        n_kernel += 1

    # The kernel might be of too low resolution, or the the wavelength range
    # might be too narrow. In both cases, raise an appropriate error
    if n_kernel == 0:
        raise ValueError(("Spectrum resolution too low for "
                          "macroturbulent broadening"))
    elif n_kernel > n_wave:
        raise ValueError(("Spectrum range too narrow for "
                          "macroturbulent broadening"))

    # Construct the broadening kernel
    wave_k = np.arange(n_kernel)*dwave
    wave_k -= wave_k[-1]/2.
    kernel = __vmacro_kernel(wave_k, 1.0, 1.0, vmac_rad, vmac_tan)
    kernel /= sum(kernel)

    flux_conv = fftconvolve(1-flux_, kernel, mode='same')

    # And interpolate the results back on to the original wavelength array,
    # taking care of even vs. odd-length kernels
    if n_kernel % 2 == 1:
        offset = 0.0
    else:
        offset = dwave / 2.0
    flux = np.interp(wave+offset, wave_, 1-flux_conv)

    # Return the results.
    if return_kernel:
        return flux, (wave_k, kernel)
    else:
        return flux





def __rotational_broadening(wave_spec,flux_spec,vrot,vmac=0.,fwhm=0.25,epsilon=0.6,
                         chard=None,stepr=0,stepi=0,alam0=None,alam1=None,
                         irel=0,cont=None,method='fortran'):
    """
    Apply rotational broadening to a spectrum assuming a linear limb darkening
    law.

    Limb darkening law is linear, default value is epsilon=0.6

    Possibility to normalize as well by giving continuum in 'cont' parameter.

    **Parameters for rotational convolution**

    C{VROT}: v sin i (in km/s):

        -  if ``VROT=0`` - rotational convolution is

                 - either not calculated,
                 - or, if simultaneously FWHM is rather large
                   (:math:`v_\mathrm{rot}\lambda/c < \mathrm{FWHM}/20.`),
                   :math:`v_\mathrm{rot}` is set to  :math:`\mathrm{FWHM}/20\cdot c/\lambda`;

        -  if ``VROT >0`` but the previous condition b) applies, the
           value of VROT is changed as  in the previous case

        -  if ``VROT<0`` - the value of abs(VROT) is used regardless of
           how small compared to FWHM it is

    C{CHARD}: characteristic scale of the variations of unconvolved stellar
    spectrum (basically, characteristic distance between two neighbouring
    wavelength points) - in A:

        - if =0 - program sets up default (0.01 A)

    C{STEPR}: wavelength step for evaluation rotational convolution;

        - if =0, the program sets up default (the wavelength
          interval corresponding to the rotational velocity
          devided by 3.)

        - if <0, convolved spectrum calculated on the original
          (detailed) SYNSPEC wavelength mesh


    **Parameters for instrumental convolution**

    C{FWHM}: WARNING: this is not the full width at half maximum for Gaussian
    instrumental profile, but the sigma (FWHM = 2.3548 sigma).

    C{STEPI}: wavelength step for evaluating instrumental convolution

          - if STEPI=0, the program sets up default (FWHM/10.)

          - if STEPI<0, convolved spectrum calculated with the previous
            wavelength mesh:
            either the original (SYNSPEC) one if vrot=0,
            or the one used in rotational convolution (vrot > 0)

    **Parameters for macroturbulent convolution**

    C{vmac}: macroturbulent velocity.

    **Wavelength interval and normalization of spectra**

    C{ALAM0}: initial wavelength
    C{ALAM1}: final wavelength
    C{IREL}: for =1 relative spectrum, =0 absolute spectrum

    @return: wavelength,flux
    @rtype: array, array
    """
    logger.info("Rot.broad with vrot={:.3f}km/s (epsilon={:.2f}), sigma={:.2f}AA, vmacro={:.3f}km/s".format(vrot,epsilon,fwhm,vmac))
    #-- first a wavelength Gaussian convolution:
    if fwhm>0:
        fwhm /= 2.3548
        #-- make sure it's equidistant
        wave_ = np.linspace(wave_spec[0],wave_spec[-1],len(wave_spec))
        flux_ = np.interp(wave_,wave_spec,flux_spec)
        dwave = wave_[1]-wave_[0]
        n = int(2*4*fwhm/dwave)
        if n==0:
            logger.info("Resolution too large, cannot broaden with instrumental profile")
        else:
            wave_k = np.arange(n)*dwave
            wave_k-= wave_k[-1]/2.
            kernel = np.exp(- (wave_k)**2/(2*fwhm**2))
            kernel /= sum(kernel)
            flux_conv = fftconvolve(1-flux_,kernel,mode='same')
            #-- this little tweak is necessary to keep the profiles at the right
            #   location
            if n%2==1:
                flux_spec = np.interp(wave_spec,wave_,1-flux_conv)
            else:
                flux_spec = np.interp(wave_spec+dwave/2,wave_,1-flux_conv)
    #-- macroturbulent profile
    if vmac>0:
        vmac = vmac/(299792458.*1e-3)*(wave_spec[0]+wave_spec[-1])/2.0
        #-- make sure it's equidistant
        wave_ = np.linspace(wave_spec[0],wave_spec[-1],len(wave_spec))
        flux_ = np.interp(wave_,wave_spec,flux_spec)
        dwave = wave_[1]-wave_[0]
        n = int(6*vmac/dwave/5)
        if n==0:
            logger.error("Resolution too large, cannot broaden with instrumental profile")
        else:
            wave_k = np.arange(n)*dwave
            wave_k-= wave_k[-1]/2.
            kernel = __vmacro_kernel(wave_k,1.,1.,vmac,vmac)
            kernel /= sum(kernel)
            flux_conv = fftconvolve(1-flux_,kernel,mode='same')
            if n%2==1:
                flux_spec = np.interp(wave_spec,wave_,1-flux_conv)
            else:
                flux_spec = np.interp(wave_spec+dwave/2,wave_,1-flux_conv)
    if vrot>0:
        #-- convert wavelength array into velocity space, this is easier
        #   we also need to make it equidistant!
        wave_ = np.log(wave_spec)
        velo_ = np.linspace(wave_[0],wave_[-1],len(wave_))
        flux_ = np.interp(velo_,wave_,flux_spec)
        dvelo = velo_[1]-velo_[0]
        vrot = vrot/(299792458.*1e-3)
        #-- compute the convolution kernel and normalise it
        n = int(2*vrot/dvelo)
        velo_k = np.arange(n)*dvelo
        velo_k -= velo_k[-1]/2.
        y = 1 - (velo_k/vrot)**2 # transformation of velocity
        G = (2*(1-epsilon)*np.sqrt(y)+np.pi*epsilon/2.*y)/(np.pi*vrot*(1-epsilon/3.0))  # the kernel
        G /= G.sum()
        #-- convolve the flux with the kernel
        flux_conv = fftconvolve(1-flux_,G,mode='same')
        if n%2==1:
            velo_ = np.arange(len(flux_conv))*dvelo+velo_[0]
        else:
            velo_ = np.arange(len(flux_conv))*dvelo+velo_[0]-dvelo/2.
        wave_conv = np.exp(velo_)
        return wave_conv,1-flux_conv


    return wave_spec,flux_spec
################################################################################

